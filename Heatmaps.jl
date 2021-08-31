using DataFrames, CSV, GLMakie, UnicodeFun, Statistics

data = DataFrame(CSV.File("Input/2020_v1.csv"))

# Average and median SWC for each location
meansSWC = []
[push!(meansSWC, mean(skipmissing(data[:, i]))) for i = 2:65]


x = [0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7] 
y = [0,1,2,2,1,0,0,1,2,2,1,0,0,1,2,2,1,0,0,1,2,2,1,0,3,4,4,3,4,3,3,4,4,3,3,4,4,3,4,3,5,6,7,7,6,5,5,6,7,7,6,5,5,6,7,7,6,5,5,6,7,7,6,5] 

x1 = x.*12.5
y1 = y.*12.5

# dry, summer, Aug 29
z = Vector(data[11600, 2:65])
#wet, winter, Jan 1
#z = Vector(data[1, 2:65])

FS = 50 

fig = Figure(resolution = (1000, 1000))
ax = Axis(fig, xlabel = "x (m)", ylabel = "y (m)",
	  xlabelsize = FS, ylabelsize = FS,
	  yticklabelsize = FS, xticklabelsize = FS,
	  #title = "Jan 1, 2020", titlesize = FS)
	 title = "Aug 29, 2020", titlesize = FS)

hmap = heatmap!(x1, y1, z, colormap = cgrad(:Spectral; categorical = true), interpolate = true)

cbar = Colorbar(fig, hmap, label = to_latex("\\theta (m^3 m^{-3})"),
	       ticklabelsize = FS, labelsize = FS,
	       width = 50)

fig[1, 1] = ax
fig[1, 2] = cbar

xlims!(ax, (0, 87.5))
ylims!(ax, (0, 87.5))

hmap.colorrange = (0.1, 0.45);
cbar.limits = (0.1, 0.45);
cbar.ticks = 0.1:0.1:0.4

fig

#save(joinpath("Output", "plot.png"), fig)
#save(joinpath("Output", "plot2.png"), fig)
# save(joinpath("Output", "plot3.png"), fig)
save(joinpath("Output", "plot4.png"), fig)






# Elevation map

data2 = DataFrame(CSV.File("Input/FieldCoordinates_c.csv"))

x2 = data2.ChamberX
y2 = data2.ChamberY
z2 = data2.OrthoHeight

x3 = x2.*12.5
y3 = y2.*12.5

FS = 50 

fig = Figure(resolution = (1000, 1000))
ax = Axis(fig, xlabel = "x (m)", ylabel = "y (m)",
	  xlabelsize = FS, ylabelsize = FS,
	  yticklabelsize = FS, xticklabelsize = FS)
	  
#hmap = heatmap!(x3, y3, z2, colormap = cgrad(:Greys_9; categorical = true), interpolate = true)

hmap = contourf!(x3, y3, z2, colormap = :grays, levels = 10)

cbar = Colorbar(fig, hmap, label = "Elevation (m)",
	       ticklabelsize = FS, labelsize = FS,
	       width = 50)

fig[1, 1] = ax
fig[1, 2] = cbar

xlims!(ax, (0, 87.5))
ylims!(ax, (0, 87.5))

fig

save(joinpath("Output", "Elevation.png"), fig)


