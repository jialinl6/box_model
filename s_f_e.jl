using EarthSciData, EarthSciMLBase, GasChem, ModelingToolkit, OrdinaryDiffEq, Dates, Unitful, DifferentialEquations

start = Dates.datetime2unix(Dates.DateTime(2016, 5, 1))

function date(t)
    return Dates.unix2datetime(t)
end

# Add unit "ppb" to Unitful 
module MyUnits
using Unitful
@unit ppb "ppb" Number 1 / 1000000000 false
end
Unitful.register(MyUnits)

struct Emission <: EarthSciMLODESystem
    sys::ODESystem
    function Emission(t)
        P = 101325 # Pa
        Tep = 298.15 # K
        R = 8.314 # Pa*m3/(mol*K)
        A = 1 # m2 asume unit area for the box model

        default_time = DateTime(2016, 5, 1)
		lon = -100.0
        lat = 30.0
        lev = 1.0
        @parameters Δz = 60 [unit = u"m"]
        @parameters uu = 1 [unit = u"ppb/s"]
        @parameters t [unit = u"s"]
        Δz2 = 60.0
        emis = NEI2016MonthlyEmis{Float64}("mrggrid_withbeis_withrwc", t, lon, lat, lev, Δz)
        fs = emis.fileset

        D = Differential(t)
        @variables NO(t) [unit = u"ppb"]
		@variables NO2(t) [unit = u"ppb"]
        @variables CH2O(t) [unit = u"ppb"]
		@variables CH4(t) [unit = u"ppb"]
		@variables CO(t) [unit = u"ppb"]
		@variables SO2(t) [unit = u"ppb"]
		@variables ISOP(t) [unit = u"ppb"]

        emis_vars = Dict{Num,Tuple{String,Float64}}(
        NO => ("NO", 30.01),
        NO2 => ("NO2", 46.0055),
        CH2O => ("FORM", 30.0260),
        CH4 => ("CH4", 16.0425),
        CO => ("CO", 28.0101),
        SO2 => ("SO2", 64.0638),
        ISOP => ("ISOP", 68.12),
        )

        eqs = [
            D(NO) ~  uu * EarthSciData.interp!(EarthSciData.DataSetInterpolator{Float64}(fs, "NO", default_time), date(t), lon, lat, lev)* 2000 / 2240 * 1e15 * R * Tep / (Δz2 * emis_vars[NO][2] * P * 3600 * 24)
            D(NO2) ~ uu * EarthSciData.interp!(EarthSciData.DataSetInterpolator{Float64}(fs, "NO2", default_time), date(t), lon, lat, lev)* 2000 / 2240 * 1e15 * R * Tep / (Δz2 * emis_vars[NO2][2] * P * 3600 * 24)
			D(CH2O) ~ uu * EarthSciData.interp!(EarthSciData.DataSetInterpolator{Float64}(fs, "FORM", default_time), date(t), lon, lat, lev)* 2000 / 2240 * 1e15 * R * Tep / (Δz2 * emis_vars[CH2O][2] * P * 3600 * 24)
			D(CH4) ~ uu * EarthSciData.interp!(EarthSciData.DataSetInterpolator{Float64}(fs, "CH4", default_time), date(t), lon, lat, lev)* 2000 / 2240 * 1e15 * R * Tep / (Δz2 * emis_vars[CH4][2] * P * 3600 * 24)
			D(CO) ~ uu * EarthSciData.interp!(EarthSciData.DataSetInterpolator{Float64}(fs, "CO", default_time), date(t), lon, lat, lev)* 2000 / 2240 * 1e15 * R * Tep / (Δz2 * emis_vars[CO][2] * P * 3600 * 24)
			D(SO2) ~ uu * EarthSciData.interp!(EarthSciData.DataSetInterpolator{Float64}(fs, "SO2", default_time), date(t), lon, lat, lev)* 2000 / 2240 * 1e15 * R * Tep / (Δz2 * emis_vars[SO2][2] * P * 3600 * 24)
			D(ISOP) ~ uu * EarthSciData.interp!(EarthSciData.DataSetInterpolator{Float64}(fs, "ISOP", default_time), date(t), lon, lat, lev)* 2000 / 2240 * 1e15 * R * Tep / (Δz2 * emis_vars[ISOP][2] * P * 3600 * 24)
        ]

        new(ODESystem(eqs, t, [NO, NO2, CH2O, CH4, CO, SO2, ISOP], [Δz, uu]; name=:Emission))
    end
end 

#visualize
# @parameters t 
# e = Emission(t)
# tspan = (start, start+3600*24*3)
# sol = solve(ODEProblem(structural_simplify(get_mtk(e)), [0.0,0,0,0,0,0,0], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
# using Plots
# plot(sol,xlabel="Time (second)", ylabel="concentration (ppb)",legend=:outerright)

ModelingToolkit.check_units(eqs...) = nothing

Base.:(+)(e::Emission, b::SuperFast) = operator_compose(b, e)
Base.:(+)(b::SuperFast, e::Emission) = e + b

@parameters t
composed_ode = SuperFast(t) + FastJX(t) + Emission(t)
composed_ode_d = SuperFast(t) + FastJX(t) + Emission(t) + DrydepositionG(t) + Wetdeposition(t)

tspan = (start, start+3600*24*3)
sys = structural_simplify(get_mtk(composed_ode))
sys_d = structural_simplify(get_mtk(composed_ode_d))

vars = states(sys)  # or variables(sys) depending on your system type and needs
var_dict = Dict(string(var) => var for var in vars)


sol = solve(ODEProblem(sys, [], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
sol_d = solve(ODEProblem(sys_d, [], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)

#visualize
using Plots 
x_t = unix2datetime.(sol[t])

pols = ["O3", "OH", "NO", "NO2", "CH4", "CH3O2", "CO","CH3OOH", "CH3O", "DMS", "SO2", "ISOP"]
var_names_p = ["superfast₊$(v)(t)" for v in pols]
pp = []
for (i, v) in enumerate(var_names_p)
    name = pols[i]
    p = Plots.plot(x_t,sol[var_dict[v]],label = "$name", size = (1000, 600), xrotation=45)
    Plots.plot!(x_t,sol_d[var_dict[v]],label = "$name+d", size = (1000, 600))
    push!(pp, p)
end
Plots.plot(pp..., layout=(3, 4))

# pols = ["O3", "OH", "NO", "NO2", "CH4", "CH3O2", "CO","CH3OOH", "CH3O", "DMS", "SO2", "ISOP"]
# var_names_p = ["superfast₊$(v)(t)" for v in pols]
# pp = []
# for (i, v) in enumerate(var_names_p)
#     name = pols[i]
#     p = Plots.plot(x_t,sol_d[var_dict[v]],label = "$name+d", size = (1000, 600), xrotation=45)
#     push!(pp, p)
# end
# Plots.plot(pp..., layout=(3, 4))