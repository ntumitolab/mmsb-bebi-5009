# # Fig 3.13
# Generalized mass action (GMA) vs. Michaelis-Menten rate laws
using Plots
plot([t -> 2t / (1+t), t -> t^0.4], 0, 4, labels=["MM" "GMA"], title="Fig. 3.13", xlabel="Substrate concentration (S / Km)", ylabel="Reaction rate (AU)")
