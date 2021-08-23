using Distributions, DataFrames, CSV, GLMakie, Dates, Random

data = DataFrame(CSV.File("Input/2020_v1.csv"))

# to do, plot histogram and density, fit normal, show average and std

# in space
test = Float64.(Vector(data[1, 2:65])) 
# in time
#test = Float64.(Vector(data[:, 2]))

# fit normal to data
d = fit(Normal, test)
mean(d)
std(d)

# plot that normal distribution
x = rand(d, 10000) # this creates a vector with a mean and std from d... but its not smooth

x = quantile.(d, 0.01:0.01:0.99) # ok this works



fig = Figure(resolution = (1500, 1500))
ax1 = Axis(fig, xlabel = "value")
hist!(ax1, test, normalization=:pdf, color = (:red, 0.5), label = "hist & pdf")
density!(ax1, test, color =(:grey, 0.5), label = "density!", strokewidth=1)
density!(ax1, x, color = (:green, 0.5), label = "normal", strokewidth = 1)
axislegend(ax1, position = :rt)
fig[1,1] = ax1
fig

# to do, recreate interactive figure for 2020, using csv dataset
# add daily normal distribution of soil moisture, soil temp and modeled/measured soil resp
# add menu and many things...

# Idea for manuscript: 
# Model Rsoil with DAMM using ANN estimated distributions of Tsoil and SWC
# --> show how it perform better than using a single value of Tsoil and SWC
# --> suggest to use it for ESM/LSM subgrids

