using DataFrames, CSV

data = DataFrame(CSV.File("Input/2020_v1.csv"))

# using DAMMmodel
include("DAMM_scaled.jl")

# get dataframe of RSM_Exp_Flux_00

d = dropmissing(data, :RSM_Exp_Flux_73)

SWC = d.SWC_73
Tsoil = d.Tsoil_73
Rsoil = Dep_var = d.RSM_Exp_Flux_73

Ind_var = hcat(Tsoil, SWC)

# porosity, the 5th param, has to be bigger than max SWC
porosity = maximum(SWC) + 0.01
lb = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # params can't be negative
ub = [Inf, Inf, Inf, Inf, porosity, Inf]
# p0 = [1e8, 59, 3.46e-8, 2.0e-3, porosity, 0.0125] 
p_scaled = [16.8, 65.5, 0.59, 3.15, porosity, 0.0026]

output1 = DAMM(Ind_var, p_scaled)

# include("AndiMD_fitf.jl") # function to scale model param

# Fit, with scaling, model DAMM to data
#res1 = curve_fit(DAMM, Ind_var, Dep_var, p_scaled, 
#		 p->[1e7*p[1], p[2], 1e-8*p[3], 1e-3*p[4], p[5], 1e-2*p[6]],
#		 p->[1e-7*p[1], p[2], 1e8*p[3], 1e3*p[4], p[5], 1e2*p[6]],
#		 inplace=false, g_tol=1E-16, x_tol=1e-16, lower=lb, upper=ub)

using LsqFit

fit = curve_fit(DAMM, Ind_var, Dep_var, p_scaled, lower=lb, upper=ub)
Param_fit = coef(fit) 

Modeled_data = DAMM(Ind_var, Param_fit) 

using GLMakie, UnicodeFun

fig = Figure()
ax3D = Axis3(fig[1,1])
p3D = scatter!(ax3D, Tsoil, SWC, Rsoil, markersize = 5000, color = :black)
p3Dm = scatter!(ax3D, Tsoil, SWC, Modeled_data, markersize = 5000, color = :red)


L = 25 # resolution
x = collect(range(1, length=L, stop=1))
[append!(x, collect(range(i, length=L, stop=i))) for i = 2:25]
x = reduce(vcat,x)
y = collect(range(0, length=L, stop=porosity))
y = repeat(y, outer=L)
x_range = hcat(x, y)

x = Int.(x_range[:, 1])
y_ax = collect(range(0, length=L, stop=porosity))
y = collect(range(1, length=L, stop=L))
y = repeat(y, outer=L)
y = Int.(y)
x_ax = collect(range(1, length=L, stop=L))

using SparseArrays
surface!(ax3D, x_ax, y_ax, Matrix(sparse(x, y, DAMM(x_range, Param_fit))), colormap = Reverse(:Spectral), transparency = true, alpha = 0.2, shading = false)

wireframe!(ax3D, x_ax, y_ax, Matrix(sparse(x, y, DAMM(x_range, Param_fit))), overdraw = true, transparency = true, color = (:black, 0.1));


ax3D.xlabel = to_latex("T_{soil} (°C)");
ax3D.ylabel = to_latex("\\theta (m^3 m^{-3})");
ax3D.zlabel = to_latex("R_{soil} (\\mumol m^{-2} s^{-1})");

fig


