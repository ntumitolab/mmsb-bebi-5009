# # Fig 6.07
# Goldbeter Koshland switch
using Plots
f607(w, K1, K2) = w * (1 - w + K1)/((1 - w) * (w + K2))
yy = 0:0.001:0.999
xx1 = f607.(yy, 1, 1)
xx2 = f607.(yy, 0.1, 0.1)
xx3 = f607.(yy, 0.01, 0.01)

plot(xx1, yy, label="K1=K2=1", xlabel="Stimulus", ylabel="Response", title="Fig 6.07", xlims=(0, 3), ylims=(0, 1))
plot!(xx2, yy, label="K1=K2=0.1")
plot!(xx3, yy, label="K1=K2=0.01")
