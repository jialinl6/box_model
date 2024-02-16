using EarthSciData, EarthSciMLBase, GasChem, ModelingToolkit, OrdinaryDiffEq, Dates, Unitful, DifferentialEquations

start = Dates.datetime2unix(Dates.DateTime(2016, 5, 1))
@parameters t
composed_ode = SuperFast(t) + FastJX(t)
tspan = (start, start+3600*24*3)
sys = structural_simplify(get_mtk(composed_ode))

vars = states(sys)  
var_dict = Dict(string(var) => var for var in vars)


sol = solve(ODEProblem(sys, [], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
using Plots 

x_t = unix2datetime.(sol[t])

pols = ["O3", "OH", "NO", "NO2", "CH4", "CH3O2", "CO","CH3OOH", "CH3O", "DMS", "SO2", "ISOP"]
var_names_p = ["superfast₊$(v)(t)" for v in pols]
pp = []
for (i, v) in enumerate(var_names_p)
    name = pols[i]
    push!(pp, Plots.plot(x_t,sol[var_dict[v]],label = "$name", size = (1000, 600), xrotation=45))
end
Plots.plot(pp..., layout=(3, 4))