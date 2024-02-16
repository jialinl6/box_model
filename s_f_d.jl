using EarthSciData, EarthSciMLBase, GasChem, ModelingToolkit, OrdinaryDiffEq, Dates, Unitful, DifferentialEquations
using DepositionMTK

ModelingToolkit.check_units(eqs...) = nothing

Base.:(+)(w::Wetdeposition, b::SuperFast) = operator_compose(b, w)
Base.:(+)(b::SuperFast, w::Wetdeposition) = w + b

Base.:(+)(d::DrydepositionG, b::SuperFast) = operator_compose(b, d)
Base.:(+)(b::SuperFast, d::DrydepositionG) = d + b

start = Dates.datetime2unix(Dates.DateTime(2016, 5, 1))
@parameters t
composed_ode = SuperFast(t) + FastJX(t) + DrydepositionG(t)
composed_ode_w = SuperFast(t) + FastJX(t) + DrydepositionG(t) + Wetdeposition(t)
tspan = (start, start+3600*24*3)
sys = structural_simplify(get_mtk(composed_ode))
sys_w = structural_simplify(get_mtk(composed_ode_w))

vars = states(sys)  # or variables(sys) depending on your system type and needs
var_dict = Dict(string(var) => var for var in vars)


sol = solve(ODEProblem(sys, [], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
sol_w = solve(ODEProblem(sys_w, [], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
using Plots 
x_t = unix2datetime.(sol[t])

pols = ["O3", "OH", "NO", "NO2", "CH4", "CH3O2", "CO","CH3OOH", "CH3O", "DMS", "SO2", "ISOP"]
var_names_p = ["superfast₊$(v)(t)" for v in pols]
pp = []
for (i, v) in enumerate(var_names_p)
    name = pols[i]
    p = Plots.plot(x_t,sol[var_dict[v]],label = "$name", size = (1000, 600), xrotation=45)
    Plots.plot!(x_t,sol_w[var_dict[v]],label = "$name+w", size = (1000, 600))
    push!(pp, p)
end
Plots.plot(pp..., layout=(3, 4))
