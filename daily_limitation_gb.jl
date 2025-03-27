# Coupled photosynthesis, stomatal and energy balance model
# Dynamic model with ODE solver

# Author: Silvere Vialet-Chabrand
# Date: 27/3/2025

using Statistics, Plots, Measures, DifferentialEquations, SciMLSensitivity, ComponentArrays, DataFrames, CSV
CA = ComponentArray

# Vialet-Chabrand 2019
function LatentHeat(Tk) 
    return 1.91846E6 * (Tk / (Tk - 33.91))^2 #latent heat of vaporization (Henderson-Sellers, 1984)
end

function SatVap(Tc)
    return 613.65 * exp(17.502 * Tc / (240.97 + Tc)) #Bucks equation
end

# Sharkey 2007
function Ta(T, c, dHa) #Harley et al. 1992
    return exp(c - dHa / (8.314e-3 * (273.15 + T))) #Dependence of reaction rate on temperature
end

function Tad(T, c, dHa, dHd, dS)
    return exp(c - dHa/(8.314e-3 * (273.15 + T)))/(1.0+exp((dS*(T + 273.15) - dHd)/(8.314e-3 * (T + 273.15))))
end


# FvCB 1980
function NRH(PAR, Abs, Jmax, Curv)
  Qabs = Abs*PAR
  return (Qabs+Jmax-sqrt((Qabs+Jmax)^2-4*Curv*Qabs*Jmax))/(2*Curv)
end

function FvCB(PAR, Ci, T, Pa, Vcmax25, Jmax25, gm25, Abs, Curv, Rd25)
    Pi = Ci * Pa * 1e-6 # Atmospheric pressure

    # Sharkey 2007
    Kc = Ta(T, 35.9774, 80.99) #Michaelis constant rubisco, Bernacchi 2002
    Ko = Ta(T, 12.3772, 23.72) #Inhibition constant
    O = 21.0 #Oxygen
    Km = Kc * (1.0 + O / Ko) #Diffusion resistance
    GS = Ta(T, 11.187, 24.46) #GS

    Vcmax = Vcmax25 * Ta(T, 26.355, 65.33)  # Bernacchi 2001
    Jmax = Jmax25 * Ta(T, 17.71, 43.9) # Bernacchi 2003
    gm = gm25 * Tad(T, 20.01, 49.6, 437.4, 1.4) #Bernacchi 2002 
    Rd = Rd25 * Ta(T, 18.715, 46.39) #Bernacchi 2001

    J = NRH(PAR, 0.5 * Abs, Jmax, Curv)

    # Ethier 2004
    a = -1.0 / gm
    b = (Vcmax - Rd) / gm + Pi + Km
    c = Rd * (Pi + Km) - Vcmax * (Pi - GS)
    d = b^2 - 4.0 * a * c
    Ac = (-b + sqrt(max(0.0, d))) / (2.0 * a)

    b = (J / 4.0 - Rd) / gm + Pi + 2.0 * GS
    c = Rd * (Pi + 2.0 * GS) - J / 4.0 * (Pi - GS)
    d = b^2 - 4.0 * a * c
    Aj = (-b + sqrt(max(0.0, d))) / (2.0 * a)

    # CO2 compensation point
    G = (Rd * Km + Vcmax * GS) / (Vcmax / Rd)
    return min(Ac, Aj), G, Rd
end

# Leuning 1995 stomatal behavior at low co2
function Leuning(A, Ds, g0, g1, Cs, G, Ds0)
  return g0 + g1 * A / ((Cs - G) * (1.0 + Ds / Ds0))
end

# Differential equations
function GE!(dy, y, p, t, env)
    PAR = env.PAR(t)
    Ca = env.Ca(t)
    Tair = env.Tair(t)
    RH = env.RH(t)
    Trefl = env.Trefl(t)
    Pa = env.Pa(t)
        
    ###
    # Leaf temperature (Vialet-Chabrand 2019)
    ## 
    rho = Pa / (287.058 * (Tair + 273.15))
    ea = SatVap(Tair) * RH
    SH = 0.622 * ea / (Pa - ea)
    Cs = 1005.0 + 1820.0 * SH

    gbm = p.gb * 8.314 * (Tair + 273.15) / Pa
    gsm = y[3] * 8.314 * (Tair + 273.15) / Pa

    dy[1] = 2.0 * 0.97 * 5.6703e-8 * ((Trefl + 273.15)^4 - (y[1] + 273.15)^4)
    dy[1] += p.Abssw * p.kSun * PAR
    dy[1] -= rho * Cs * 0.92 * gbm * (y[1] - Tair)

    lambda = LatentHeat(y[1] + 273.15)
    vpd = SatVap(y[1]) - ea
    gtm = 1.0 / (1.0 / gsm + 1.0 / gbm)
    dy[1] -= lambda * 0.622 * rho / Pa * gtm * vpd

    dy[1] /= p.k

    ###
    # Photosynthesis
    ##
    gtc = 1.0 / (1.6 / y[3] + 1.37 / p.gb)
    Ci = Ca - y[2] / gtc
    A, G, Rd = FvCB(PAR, Ci, y[1], Pa, p.Vcmax25, p.Jmax25, p.gm25, p.Abs, p.Curv, p.Rd25)
    dy[2] = (A - y[2]) / (y[2] < A ? p.tauAi : p.tauAd)

    ###
    # Stomatal conductance
    ##
    Csurface = Ca - 1.37 * A / p.gb
    # Ds = vpd / 1e3
    # Aphalo (1993), Eq 8
    gt = 1.0 / (1.0 / p.gb + 1.0 / y[3])
    Ds = vpd * (1.0 - gt / p.gb) * 1e-3 #Correction factor for ambient to leaf surface and in kPa
    gs = Leuning(A + Rd, Ds, p.g0, p.g1, Csurface, G, p.Ds0)
    dy[3] = (gs - y[3]) / (y[3] < gs ? p.tauGi : p.tauGd)

    return nothing
end

function pred(time, p, env)
    ode = (dy, y, p, t) -> GE!(dy, y, p, t, env)
    y0 = [env.Tair(0.0), -p.Rd25, p.g0]
    prob = ODEProblem{true, SciMLBase.FullSpecialize}(ode, y0, (0.0, time[end]), p)

    sol = solve(remake(prob, p=p), RadauIIA5(),
                saveat=time, maxiters=1e6, abstol=1e-9, reltol=1e-10,
                sensealg=InterpolatingAdjoint(checkpointing=true, autojacvec=ZygoteVJP()))
    
    return Array(sol)
end

function plotSim(time, p, env)
    tmp = pred(time, p, env)

    p1 = plot(time ./ 3600.0, view(tmp, 1, :), label="Tl obs.")
    p2 = plot(time ./ 3600.0, view(tmp, 2, :), label="A obs.")
    p3 = plot(time ./ 3600.0, view(tmp, 3, :), label="gs obs.")

    pg = plot(p1, p2, p3, layout=(3,1), size=(840, 760))
    pg |> display
end

function plotSimR(time, p, env)
    p.gb = 0.20 #note: input is double sided, divide by 2 for one sided (here amphistomatous leaf)
    tmp1 = pred(time, p, env)
    p.gb = 0.3
    tmp2 = pred(time, p, env)
    p.gb = 0.6
    tmp3 = pred(time, p, env)
    p.gb = 2
    tmp4 = pred(time, p, env)

    p1 = plot(time ./ 3600.0, env.PAR.(time), xlabel="Time (hours)", ylabel="PAR")

    p2 = plot(time ./ 3600.0, view(tmp1, 1, :), label="0.20", xlabel="Time (hours)", ylabel="Temperature", ylim=(0, 40))
    plot!(time ./ 3600.0, view(tmp2, 1, :), label="0.3")
    plot!(time ./ 3600.0, view(tmp3, 1, :), label="0.6")
    plot!(time ./ 3600.0, view(tmp4, 1, :), label="2")

    p3 = plot(time ./ 3600.0, view(tmp1, 2, :), label="0.20", xlabel="Time (hours)", ylabel="A", ylim=(0, 40))
    plot!(time ./ 3600.0, view(tmp2, 1, :), label="0.3")
    plot!(time ./ 3600.0, view(tmp3, 1, :), label="0.6")
    plot!(time ./ 3600.0, view(tmp4, 1, :), label="2")


    p4 = plot(time ./ 3600.0, view(tmp1, 3, :), label="0.20", xlabel="Time (hours)", ylabel="Gs")
    plot!(time ./ 3600.0, view(tmp2, 1, :), label="0.3")
    plot!(time ./ 3600.0, view(tmp3, 1, :), label="0.6")
    plot!(time ./ 3600.0, view(tmp4, 1, :), label="2")


    pg = plot(p1, p2, p3, p4, layout=(2,2), size=(740, 740))
    pg |> display
    
    #Make csv file for further processing   
    df = DataFrame(time = time)

    df[!, "PAR_leaf"] = env.PAR(time)
    df[!, "temp"]= env.Tair.(time)
    df[!, "RH"] = env.RH.(time)
    df[!, "Ca"] = env.Ca.(time)

    df[!, "Temp_0.20"] = view(tmp1, 1, :)
    df[!, "Temp_0.3"] = view(tmp2, 1, :)
    df[!, "Temp_0.6"] = view(tmp3, 1, :)
    df[!, "Temp_2"] = view(tmp4, 1, :)

    df[!, "A_0.20"] = view(tmp1, 2, :)
    df[!, "A_0.3"] = view(tmp2, 2, :)
    df[!, "A_0.6"] = view(tmp3, 2, :)
    df[!, "A_2"] = view(tmp4, 2, :)

    df[!, "Gs_0.20"] = view(tmp1, 3, :)
    df[!, "Gs_0.3"] = view(tmp2, 3, :)
    df[!, "Gs_0.6"] = view(tmp3, 3, :)
    df[!, "Gs_2"] = view(tmp4, 3, :)

 CSV.write("C:\\Julia_files\\Top\\simulation_results_tomato.csv", df)
end

# Kaiser 2020, Zhang 2022, McAusland 2016
function init()
    CA(gb = 0.1, # Boundary layer conductance (mol m-2 s-1)
       k = 800.0, # A amount of energy per unit area required to change the temperature of the material by 1 °K (J m–2 K–1) 
       Vcmax25 = 141, # Maximum rate of carboxylation at 25 degrees C (umol m-2 s-1)
       Jmax25 = 186, # Maximum rate of electron transport at 25 degrees C (umol m-2 s-1)
       gm25 = 3.2, # Mesophyll conductance at 25 degrees C (mol m-2 s-1 Pa-1)
       Abs = 0.84, # Leaf absorbance
       Abssw = 0.5, #Leaf absorbance short wave radiation
       Curv = 0.7, # Curvature of the J/PAR response
       Rd25 = 1.5, # Respiration at 25 degrees C (umol m-2 s-1)
       tauAi = 300.0, # Time constant of gs increase (s)
       tauAd = 1.0, # Time constant of A decrease (s)
       kSun = 1/(500/120*0.47), # Light convertion factor from umol m-2 s-1 to W m-2 #book greenhouse horticulture Stanghellini 2019
       g0 = 0.05, # Nocturnal gs (mol m-2 s-1)
       g1 = 11.0, # Slope of the relationship between A/gs
       Ds0 = 1.5, # Sensitivity factor for VPD response (KPa)
       tauGi = 900.0, # Time constant of gs increase (s)
       tauGd = 450.0) # Time constant of gs decrease (s)
end

function envir(eod)
    (PAR = x -> 1500.0 * sin.(pi / eod * x), # Light intensity (umol m-2 s-1)
     Ca = x -> 425.0, # Air CO2 concentation (umol mol-1)
     Tair = x -> 32.0, # Air temperature (degree C)
     RH = x -> 0.5, # Air relative humidity
     Trefl = x -> 33.0, # Surrounding temperature (degree C)
     Pa = x -> 101300.0) # Atmospheric pressure (Pa)
end

function Run()
    # Load datasets
    time = 0.0:10.0:(16 * 3600.0)
    p = init()
    env = envir(time[end])

    #plotSim(time, p, env)
    plotSimR(time, p, env)
    savefig("Sim_425ppm_tomato.png")
end

Run()
