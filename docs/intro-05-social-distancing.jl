#===
# Agent-based model

COVID-19 social distancing model

Source: [Agents.jl model zoo](https://juliadynamics.github.io/AgentsExampleZoo.jl/dev/examples/social_distancing/)
===#

using Agents
using Random

# Let us first create a simple model where balls move around in a continuous space. We need to create agents that comply with `ContinuousSpace`, i.e. they have a pos and vel fields, both of which are tuples of float numbers.

mutable struct Ball <: AbstractAgent
    id::Int                 ## Mandatory Agent identifier
    pos::NTuple{2,Float64}  ## Position, required for agents in the ContinuousSpace
    vel::NTuple{2,Float64}  ## Moving speeds
    mass::Float64           ## Can move or not
end

#---

function ball_model(; speed = 0.002, seed = 42)
    space2d = ContinuousSpace((1, 1); spacing=0.02)
    model = ABM(Ball, space2d, properties = Dict(:dt => 1.0), rng = MersenneTwister(seed))

    for i in 1:500
        pos = Tuple(rand(model.rng, 2))
        vel = sincos(2π * rand(model.rng)) .* speed
        mass = 1.0
        add_agent!(pos, model, vel, mass)
    end
    return model
end

#---

using InteractiveDynamics
using CairoMakie

agent_step!(agent::Ball, model) = move_agent!(agent, model, model.dt)

model = ball_model()

Agents.abmvideo(
    "socialdist1.mp4",
    model,
    agent_step!;
    title = "Ball Model",
    frames = 50,
    spf = 2,
    framerate = 25,
)

#--

using Base64

function display_mp4(filename)
    display("text/html", string("""<video autoplay controls><source src="data:video/x-m4v;base64,""",
        Base64.base64encode(open(read, filename)),"""" type="video/mp4"></video>"""))
end

display_mp4("socialdist1.mp4")

# As you can see the agents move in a straight line in a periodic space without interactions. Let's change that.

#===
## Billiard-like interaction

Using the continuous space API:

- `interacting_pairs()`
- `elastic_collision!()`

And we redefine the stepping function:
===#

function model_step!(model)
    for (a1, a2) in interacting_pairs(model, 0.012, :nearest)
        elastic_collision!(a1, a2, :mass)
    end
end

model2 = ball_model()

Agents.abmvideo(
    "socialdist2.mp4",
    model2,
    agent_step!,
    model_step!;
    title = "Billiard-like",
    frames = 50,
    spf = 2,
    framerate = 25,
)

display_mp4("socialdist2.mp4")

#===
## Immovable agents

For the following social distancing example, it will become crucial that some agents don't move, and can't be moved (i.e. they stay "isolated"). This is very easy to do with the elastic_collision! function, we only have to make some agents have infinite mass.
===#

model3 = ball_model()

for i in 1:400
    agent = model3[i]
    agent.mass = Inf
    agent.vel = (0.0, 0.0)
end

Agents.abmvideo(
    "socialdist3.mp4",
    model3,
    agent_step!,
    model_step!;
    title = "Billiard-like with stationary agents",
    frames = 50,
    spf = 2,
    framerate = 25,
)

display_mp4("socialdist3.mp4")

#===
## Adding Virus spread (SIR model)

The agents can be infected with a disease and transfer the disease to other agents around them.
===#

mutable struct Person <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}
    vel::NTuple{2,Float64}
    mass::Float64
    days_infected::Int  # number of days since is infected
    status::Symbol  # :S, :I or :R
    β::Float64
end

#--

const steps_per_day = 24 ## One tick per hour

#---

function init_sir(;
    infection_period = 30 * steps_per_day,
    detection_time = 14 * steps_per_day,
    reinfection_probability = 0.05,
    isolated = 0.0, # in percentage
    interaction_radius = 0.012,
    dt = 1.0,
    speed = 0.002,
    death_rate = 0.044,
    N = 1000,
    initial_infected = 5,
    seed = 42,
    βmin = 0.4,
    βmax = 0.8,
)

    properties = (;
        infection_period,
        reinfection_probability,
        detection_time,
        death_rate,
        interaction_radius,
        dt,
    )
    space = ContinuousSpace((1, 1), spacing=0.02)
    model = ABM(Person, space, properties = Dict(pairs(properties)), rng = MersenneTwister(seed))

    ## Add initial individual agents
    for ind in 1:N
        pos = Tuple(rand(model.rng, 2))
        status = ind ≤ N - initial_infected ? :S : :I
        isisolated = ind ≤ isolated * N
        mass = isisolated ? Inf : 1.0
        vel = isisolated ? (0.0, 0.0) : sincos(2π * rand(model.rng)) .* speed

        β = (βmax - βmin) * rand(model.rng) + βmin
        add_agent!(pos, model, vel, mass, 0, status, β)
    end

    return model
end

# To visualize this model, we will use black color for the susceptible, red for the infected infected and green for the recovered.

sir_model = init_sir()

sir_colors(a) = a.status == :S ? "#2b2b33" : a.status == :I ? "#bf2642" : "#338c54"

fig, abmstepper = Agents.abmplot(sir_model; ac = sir_colors)
fig

# Modify the `model_step!` function to simulate disease transmission.

function transmit!(a1::Person, a2::Person, model)

    rp = model.reinfection_probability

    ## for transmission, only 1 can have the disease (otherwise nothing happens)
    if count(a.status == :I for a in (a1, a2)) ≠ 1
        return nothing
    end

    infected, healthy = a1.status == :I ? (a1, a2) : (a2, a1)

    ## Lucky and not infected
    if rand(model.rng) > infected.β
        return nothing
    end

    ## Risk of reinfection
    if healthy.status == :R && rand(model.rng) > rp
        return nothing
    end

    ## You got virus
    healthy.status = :I

    return nothing
end

function sir_model_step!(model)
    r = model.interaction_radius
    for (a1, a2) in interacting_pairs(model, r, :all)
        transmit!(a1, a2, model)
        elastic_collision!(a1, a2, :mass)
    end
    return nothing
end

# Agent-specific functions
function update!(agent::Person)
    if agent.status == :I
        agent.days_infected += 1
    end
    return nothing
end

function recover_or_die!(agent::Person, model)
    if agent.days_infected ≥ model.infection_period
        if rand(model.rng) ≤ model.death_rate
            kill_agent!(agent, model)
        else
            agent.status = :R
            agent.days_infected = 0
        end
    end
    return nothing
end

function sir_agent_step!(agent::Person, model)
    move_agent!(agent, model, model.dt)
    update!(agent)
    recover_or_die!(agent, model)
end

# Running with default parameters.
sir_model = init_sir()

Agents.abmvideo(
    "socialdist4.mp4",
    sir_model,
    sir_agent_step!,
    sir_model_step!;
    title = "SIR model",
    frames = 200,
    ac = sir_colors,
    as = 10,
    spf = 2,
    framerate = 20,
)

display_mp4("socialdist4.mp4")

# ## Analyzing exponential spread

infected(x) = count(i == :I for i in x)
recovered(x) = count(i == :R for i in x)

# Aggregated data for number of infected and recovered indivisuals
adata = [(:status, infected), (:status, recovered)]

# Try different parameters
r1, r2 = 0.02, 0.05
β1, β2 = 0.5, 0.1
sir_model1 = init_sir(reinfection_probability = r1, βmax = β1)
sir_model2 = init_sir(reinfection_probability = r2, βmax = β1)
sir_model3 = init_sir(reinfection_probability = r1, βmax = β2)

data1, _ = run!(sir_model1, sir_agent_step!, sir_model_step!, 3000; adata)
data2, _ = run!(sir_model2, sir_agent_step!, sir_model_step!, 3000; adata)
data3, _ = run!(sir_model3, sir_agent_step!, sir_model_step!, 3000; adata)

data1[(end-10):end, :]

#---

using CairoMakie

figure = Figure()
ax = figure[1, 1] = Axis(figure; ylabel = "Infected", xlabel="Steps")
l1 = lines!(ax, data1[:, dataname((:status, infected))], color = :orange)
l2 = lines!(ax, data2[:, dataname((:status, infected))], color = :blue)
l3 = lines!(ax, data3[:, dataname((:status, infected))], color = :green)
figure[1, 2] = Legend(figure, [l1, l2, l3], ["r=$r1, beta=$β1", "r=$r2, beta=$β1", "r=$r1, beta=$β2"])

figure

#====
## Social distancing

The best way to model social distancing is to make some agents simply not move (which feels like it approximates reality better).
====#

sir_model = init_sir(isolated = 0.85)

Agents.abmvideo(
    "socialdist5.mp4",
    sir_model,
    sir_agent_step!,
    sir_model_step!;
    title = "Social Distancing",
    frames = 200,
    spf = 2,
    ac = sir_colors,
    framerate = 20,
)

display_mp4("socialdist5.mp4")

#---

r4 = 0.02
sir_model4 = init_sir(reinfection_probability = r4, βmax = β1, isolated = 0.85)

data4, _ = run!(sir_model4, sir_agent_step!, sir_model_step!, 3000; adata)

figure = Figure()
ax = figure[1, 1] = Axis(figure; ylabel = "Infected", xlabel="Steps")
l1 = lines!(ax, data1[:, dataname((:status, infected))], color = :orange)
l2 = lines!(ax, data2[:, dataname((:status, infected))], color = :blue)
l3 = lines!(ax, data3[:, dataname((:status, infected))], color = :green)
l4 = lines!(ax, data4[:, dataname((:status, infected))], color = :red)
figure[1, 2] = Legend(
    figure,
    [l1, l2, l3, l4],
    ["r=$r1, beta=$β1", "r=$r2, beta=$β1", "r=$r1, beta=$β2", "r=$r4, social distancing"],
)

figure
