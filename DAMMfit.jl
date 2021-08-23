using DataFrames, CSV, LsqFit

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

using GLMakie, UnicodeFun

i = 63
poro_val = poro_vals[i]

fig = Figure()
ax3D = Axis3(fig[1,1])
p3D = scatter!(ax3D, Tsoil[i], SWC[i], Rsoil[i], markersize = 5000, color = :black)
p3Dm = scatter!(ax3D, Tsoil[i], SWC[i], Modeled_data[i], markersize = 5000, color = :red)


L = 25 # resolution
x = collect(range(1, length=L, stop=1))
[append!(x, collect(range(i, length=L, stop=i))) for i = 2:25]
x = reduce(vcat,x)
y = collect(range(0, length=L, stop=poro_vals[i]))
y = repeat(y, outer=L)
x_range = hcat(x, y)

x = Int.(x_range[:, 1])
y_ax = collect(range(0, length=L, stop=poro_vals[i]))
y = collect(range(1, length=L, stop=L))
y = repeat(y, outer=L)
y = Int.(y)
x_ax = collect(range(1, length=L, stop=L))

using SparseArrays
surface!(ax3D, x_ax, y_ax, Matrix(sparse(x, y, DAMM(x_range, Param_fit[i]))), colormap = Reverse(:Spectral), transparency = true, alpha = 0.2, shading = false)

wireframe!(ax3D, x_ax, y_ax, Matrix(sparse(x, y, DAMM(x_range, Param_fit[i]))), overdraw = true, transparency = true, color = (:black, 0.1));


ax3D.xlabel = to_latex("T_{soil} (°C)");
ax3D.ylabel = to_latex("\\theta (m^3 m^{-3})");
ax3D.zlabel = to_latex("R_{soil} (\\mumol m^{-2} s^{-1})");

fig

# test

