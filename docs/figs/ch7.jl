# # Chapter 7
# ## Fig. 7.31
# Hasty synthetic oscillator model
using OrdinaryDiffEq
using Catalyst
using Plots

# Model
v0_731(x, y, alpha, sigma) = (1 + x^2 + alpha * sigma * x^4) / ((1 + x^2 + sigma * x^4) * (1 + y^4))

@time "Build system" rn731 = @reaction_network begin
    v0_731(x1, y1, alpha, sigma), 0 --> x1
    ay * v0_731(x1, y1, alpha, sigma), 0 --> y1
    v0_731(x2, y2, alpha, sigma), 0 --> x2
    ay * v0_731(x2, y2, alpha, sigma), 0 --> y2
    gammax, x1 --> 0
    gammay, y1 --> 0
    gammax, x2 --> 0
    gammay, y2 --> 0
    (D, D), x1 <--> x2
    (D, D), y1 <--> y2
end
