using DataFrames, CSV

data = DataFrame(CSV.File("Input/2020_v1.csv"))

using LsqFit, DAMMmodel

# get dataframe of RSM_Exp_Flux_00

d = dropmissing(data, :RSM_Exp_Flux_76)

SWC = d.SWC_76
Tsoil = d.Tsoil_76
Rsoil = Dep_var = d.RSM_Exp_Flux_76

Ind_var = hcat(Tsoil, SWC)

# porosity, the 5th param, has to be bigger than max SWC
porosity = maximum(SWC) + 0.001
lb = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # params can't be negative
ub = [Inf, Inf, Inf, Inf, Inf, Inf]
p = Param_ini = [1e8, 59, 3.46e-8, 2.0e-3, porosity, 0.0125] 
output1 = DAMM(Ind_var, p)

# Fit DAMM parameters to data,
fit = curve_fit(DAMM, Ind_var, Dep_var, Param_ini, lower=lb, upper=ub)
Param_fit = coef(fit) 


Modeled_data = DAMM(Ind_var, Param_fit) 

using GLMakie, UnicodeFun

fig = Figure()
ax3D = Axis3(fig[1,1])
p3D = scatter!(ax3D, Tsoil, SWC, Rsoil, markersize = 5000, color = :black)
p3Dm = scatter!(ax3D, Tsoil, SWC, Modeled_data, markersize = 5000, color = :red)


L = 22 # resolution
x = collect(range(1, length=L, stop=1))
[append!(x, collect(range(i, length=L, stop=i))) for i = 2:22]
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


