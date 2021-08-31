using DataFrames, CSV, LsqFit, Statistics, StatsMakie

data = DataFrame(CSV.File("Input/2020_v1.csv"))

# using DAMMmodel
include("DAMM_scaled_porosity.jl")

# to do: loop to DAMM get parameters for all of 64 locations
x = [0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7] 
y = [0,1,2,2,1,0,0,1,2,2,1,0,0,1,2,2,1,0,0,1,2,2,1,0,3,4,4,3,4,3,3,4,4,3,3,4,4,3,4,3,5,6,7,7,6,5,5,6,7,7,6,5,5,6,7,7,6,5,5,6,7,7,6,5] 

lb = [0.0, 0.0, 0.0, 0.0, 0.0] # params can't be negative
ub = [Inf, Inf, Inf, Inf, Inf] 
p_scaled = [16.8, 65.5, 0.59, 3.15, 0.0026]

SWC = []
Tsoil = []
Rsoil = []
poro_vals = []
fit = []
Param_fit = []
Modeled_data = []
poro_val = 0.4 # to avoid weird bug
for i = 1:64
	name_RS = string("RSM_Exp_Flux_", x[i], y[i])
	name_SWC = string("SWC_", x[i], y[i])
	name_Tsoil = string("Tsoil_", x[i], y[i])

	d = dropmissing(data, name_RS)
	d = dropmissing(d, name_SWC)
	push!(SWC, d[!, name_SWC])
	push!(Tsoil, d[!, name_Tsoil])
	push!(Rsoil, d[!, name_RS])
	Ind_var = hcat(Tsoil[i], SWC[i])
	push!(poro_vals, maximum(skipmissing(data[!, name_SWC])))
	poro_val = poro_vals[i]
	push!(fit, curve_fit(DAMM, Ind_var, Rsoil[i], p_scaled, lower=lb, upper=ub))
	push!(Param_fit, coef(fit[i]))
	push!(Modeled_data, DAMM(Ind_var, Param_fit[i]))
end

# means
# Average and median SWC for each location
meansSWC = []
[push!(meansSWC, mean(skipmissing(data[:, i]))) for i = 2:65]
#[push!(meansSWC, mean(data[:, i])) for i = 2:65]

# Elevation on locations, joinleft?
data2 = DataFrame(CSV.File("Input/FieldCoordinates_c.csv"))
data2.loc = @. string(data2.ChamberX, data2.ChamberY)

# param dataframe
df_p = DataFrame("x" => x,
		 "y" => y,
		 "a" => [Param_fit[i][1] for i = 1:64],
		"Ea" => [Param_fit[i][2] for i = 1:64],
		"kMsx" => [Param_fit[i][3] for i = 1:64],
		"kMo2" => [Param_fit[i][4] for i = 1:64],
		"sx" => [Param_fit[i][5] for i = 1:64],
		"SWCmean" => meansSWC)
df_p.loc = @. string(x, y)

test = innerjoin(data2, df_p, on = :loc)

CSV.write(joinpath("Output", "params.csv"), df_p)
CSV.write(joinpath("Output", "test.csv"), test)

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Elevation (m)", ylabel = "Mean SWC (m3 m-3)")
plot!(test.OrthoHeight, Float64.(test.SWCmean))
fig

# how to plot linear regression in Makie?
# plot(linear, Float64.(test.SWCmean), test.OrthoHeight)

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "mean SWC", ylabel = "kMO2")
plot!(Float64.(test.SWCmean), test.kMo2)
ylims!(-0.1, 2)
fig

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "mean SWC", ylabel = "kMsx")
plot!(Float64.(test.SWCmean), test.kMsx)
ylims!(-0.1, 2)
fig

# We can estimate errors on the fit parameters,
# to get standard error of each parameter:
sigma = stderror(fit[1])
# to get margin of error and confidence interval of each parameter at 5% significance level:
margin_of_error = margin_error(fit[1], 0.05)
confidence_inter = confidence_interval(fit[1], 0.05)
# parameter covariance matrix evaluated at the best fit point
covar = estimate_covar(fit[1])
# Calculate RMSE, sqrt of sum of square of residuals
# correlation (r2) of param from covar --> see with Julie
# v1 / sqrt(diag(v1)) %*% t(sqrt(diag(v1)))
# s_ij/si*sj





# saving all DAMM Matrix before plotting observable
using SparseArrays

DAMM_Matrix = []
y_ax = []

L = 25 # resolution
x_ax = collect(range(1, length=L, stop=L))
for i = 1:64
	poro_val = poro_vals[i]

	xD = collect(range(1, length=L, stop=1))
	[append!(xD, collect(range(i, length=L, stop=i))) for i = 2:25]
	xD = reduce(vcat, xD)
	yD = collect(range(0, length=L, stop=poro_vals[i]))
	yD = repeat(yD, outer=L)
	x_range = hcat(xD, yD)

	xD = Int.(x_range[:, 1])
	push!(y_ax, collect(range(0, length=L, stop=poro_vals[i])))
	yD = collect(range(1, length=L, stop=L))
	yD = repeat(yD, outer=L)
	yD = Int.(yD)
	
	push!(DAMM_Matrix, Matrix(sparse(xD, yD, DAMM(x_range, Param_fit[i]))))
end




using GLMakie, UnicodeFun

fig = Figure()

locations = [string(x[i], y[i]) for i in 1:64]
vals = collect(1:64)

menu = Menu(fig, options = zip(locations, vals))

loc = Node{Any}(vals[1])

data3D = @lift(Vec3f0.(Tsoil[$loc], SWC[$loc], Rsoil[$loc]))
model3D = @lift(Vec3f0.(Tsoil[$loc], SWC[$loc], Modeled_data[$loc]))
y_ax_N = @lift(y_ax[$loc])
surface3D = @lift(DAMM_Matrix[$loc])

ax3D = Axis3(fig[1,1])
ax3D.xlabel = to_latex("T_{soil} (°C)");
ax3D.ylabel = to_latex("\\theta (m^3 m^{-3})");
ax3D.zlabel = to_latex("R_{soil} (\\mumol m^{-2} s^{-1})");

p3D = scatter!(ax3D, data3D, markersize = 5000, color = :black)
p3Dm = scatter!(ax3D, model3D, markersize = 5000, color = :red)
surface!(ax3D, x_ax, y_ax_N, surface3D, colormap = Reverse(:Spectral), transparency = true, alpha = 0.2, shading = false)
wireframe!(ax3D, x_ax, y_ax_N, surface3D, overdraw = true, transparency = true, color = (:black, 0.1));

on(menu.selection) do s
	loc[] = s
	autolimits!(ax3D)
end

fig

# plot individual figures


for i = 1:64
	fig = Figure()

	data3D = Vec3f0.(Tsoil[i], SWC[i], Rsoil[i])
	model3D = Vec3f0.(Tsoil[i], SWC[i], Modeled_data[i])
	y_ax_N = y_ax[i]
	surface3D = DAMM_Matrix[i]

	ax3D = Axis3(fig[1,1])
	ax3D.xlabel = to_latex("T_{soil} (°C)");
	ax3D.ylabel = to_latex("\\theta (m^3 m^{-3})");
	ax3D.zlabel = to_latex("R_{soil} (\\mumol m^{-2} s^{-1})");

	p3D = scatter!(ax3D, data3D, markersize = 5000, color = :black)
	p3Dm = scatter!(ax3D, model3D, markersize = 5000, color = :red)
	surface!(ax3D, x_ax, y_ax_N, surface3D, colormap = Reverse(:Spectral), transparency = true, alpha = 0.2, shading = false)
	wireframe!(ax3D, x_ax, y_ax_N, surface3D, overdraw = true, transparency = true, color = (:black, 0.1));
	name = string("location_", locations[i], ".png")

	ylims!(0, 0.45)
	zlims!(0,25)

	save(joinpath("Output", "Individual plots DAMM", name), fig)
end


