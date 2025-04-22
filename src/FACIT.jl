```
    The FACIT impurity transport model calculates impurity transport coefficients. Main references include: 
    - Maget et al 2020 Plasma Phys. Control. Fusion 62 105001
    - Fajardo et al 2022 Plasma Phys. Control. Fusion 64 055017
    - Maget et al 2022 Plasma Phys. Control. Fusion 64 069501
    - Fajardo et al 2023 Plasma Phys. Control. Fusion 65 035021 

    This implementation is largely based on the Python version that can be found at: 
    https://github.com/fsciortino/Aurora/blob/master/aurora/facit.py

```
using LinearAlgebra
using IMAS
using Trapz

# -----------------
# Physical Constants
# -----------------
const mp = IMAS.mks.m_p
const me = IMAS.mks.m_e
const qe = IMAS.mks.e
const eps0 = IMAS.mks.ϵ_0
const eps_pi_fac = 3 * eps0^2 * (2pi/qe)^1.5 / qe

Base.@kwdef mutable struct FACITinput
    rho::Union{Vector{Float64},Missing} = missing
    Zimp::Union{Vector{Float64},Missing} = missing
    Aimp::Union{Float64,Missing} = missing
    Zi::Union{Float64,Missing} = missing
    Ai::Union{Float64,Missing} = missing
    Ti::Union{Vector{Float64},Missing} = missing
    Ni::Union{Vector{Float64},Missing} = missing
    Nimp::Union{Vector{Float64},Missing} = missing
    Machi::Union{Vector{Float64},Missing} = missing
    Zeff::Union{Vector{Float64},Missing} = missing
    gradTi::Union{Vector{Float64},Missing} = missing
    gradNi::Union{Vector{Float64},Missing} = missing
    gradNimp::Union{Vector{Float64},Missing} = missing
    invaspct::Union{Float64,Missing} = missing
    B0::Union{Float64,Missing} = missing
    R0::Union{Float64,Missing} = missing
    qmag::Union{Vector{Float64},Missing} = missing
    rotation_model::Union{Int,Missing} = missing
    Te_Ti::Union{Vector{Float64},Missing} = missing
    RV::Union{Matrix{Float64},Missing} = missing
    FV::Union{Vector{Float64},Missing} = missing
    fsaout::Union{Bool,Missing} = missing
    full_geom::Union{Bool,Missing} = missing
    fH::Union{Float64,Missing} = missing
    bC::Union{Float64,Missing} = missing
    sigH::Union{Float64,Missing} = missing
    TperpTpar_axis::Union{Float64,Missing} = missing
end

Base.@kwdef mutable struct FACIToutput
    Dz::Union{Vector{Float64},Missing} = missing
    Kz::Union{Vector{Float64},Missing} = missing
    Hz::Union{Vector{Float64},Missing} = missing
    Vconv::Union{Vector{Float64},Missing} = missing
    Vrz::Union{Vector{Float64},Missing} = missing
    Flux_z::Union{Vector{Float64},Missing} = missing
end

function FACITinput(
    rho::Vector{Float64}, Zimp::Union{Vector{Float64}, Float64}, Aimp::Float64,
    Zi::Float64, Ai::Float64,
    Ti::Vector{Float64}, Ni::Vector{Float64}, Nimp::Vector{Float64},
    Machi::Vector{Float64}, Zeff::Vector{Float64},
    gradTi::Vector{Float64}, gradNi::Vector{Float64}, gradNimp::Vector{Float64},
    invaspct::Float64, B0::Float64, R0::Float64, qmag::Vector{Float64};
    rotation_model::Int = 0,
    Te_Ti::Union{Float64, Vector{Float64}} = 1.0, RV::Union{Matrix{Float64},Missing}, FV::Union{Vector{Float64},Missing}, fsaout::Bool = true, full_geom::Bool = false, fH::Float64 = 0.0,
    bC::Float64 = 0.0, sigH::Float64 = 1.0, TperpTpar_axis::Float64 = 1.0
)
    nr = length(rho)
    if length(Zimp) == 1
        Zimp = fill(Zimp, nr)
    end
    if typeof(Te_Ti) == Float64
        Te_Ti = fill(Te_Ti, nr)
    end

    # Replace Nimp with Ni if mean(Nimp) < 1
    if sum(Nimp) / length(Nimp) < 1
        Nimp = copy(Ni)
    end

    eps = max.(0.005, rho .* invaspct)  # minimum epsilon
    nth = 20
    theta = range(0, 2π, length = nth)

    if ismissing(RV) || ismissing(ZV)
        # make circular geometry
        RV = R0 .* (1 .+ eps .* cos.(theta)')
        JV = B0 ./ (qmag .* RV)
    end

    Dz = zeros(nr)
    Kz = zeros(nr)
    Hz = zeros(nr)
    Vconv = zeros(nr)
    Vrz = zeros(nr)
    Flux_z = zeros(nr)

    return FACITinput(rho, Zimp, Aimp, Zi, Ai, Ti, Ni, Nimp, Machi, Zeff,
                      gradTi, gradNi, gradNimp, invaspct, B0, R0, qmag, rotation_model, Te_Ti, RV, FV,
                      fsaout, full_geom, fH, bC, sigH, TperpTpar_axis)
end

function compute_transport(input::FACITinput)
    rho = max.(input.rho, 1e-6)   # avoid rho = 0
    nr = length(input.rho)
    eps = max.(0.005, rho .* input.invaspct)

    amin = input.invaspct*input.R0
    Machi2 = input.Machi .^ 2

    grad_ln_Ni = input.gradNi ./ input.Ni
    grad_ln_Ti = input.gradTi ./ input.Ti
    grad_ln_Nimp = input.gradNimp ./ input.Nimp

    Tauii, Tauimpi, Tauiimp, Tauimpimp = collision_times(input.Zi, input.Zimp, input.Ni, input.Nimp, input.Ti, input.Ai, input.Aimp)
    ft = ftrap(eps)
    Mzstar = @. sqrt(input.Aimp / input.Ai - (input.Zimp / input.Zi) * input.Zeff / (input.Zeff + 1/input.Te_Ti)) * input.Machi
    f1, f2, f3, fG, fU, yy, fv, fdps, fhbp = facs(input.Zimp, input.Aimp, ft, Mzstar, input.rotation_model)
    K11i, K12i, K11z, K12z = KVISC(input.Nimp, input.Ni, input.Ti, input.Ai, input.Aimp, input.Zi, input.Zimp, Tauii, Tauimpimp, Tauiimp, Tauimpi, eps, ft, input.R0, input.qmag, yy)

    mi = input.Ai * mp
    g = (input.qmag .* input.R0) ./ (sqrt.(2 .* qe .* input.Ti ./ mi) .* Tauii)
    nuistar = g ./ eps.^1.5
    alpha = (input.Zimp.^2 .* input.Nimp) ./ (input.Zi^2 .* input.Ni)
    ki = ki_f(nuistar, ft, input.Zeff, input.Machi)
    
    eps2 = eps.^2
    C0z = C2(alpha, g, input.Ai, input.Aimp, f1, f2)
    nuz = 1.0 ./ (sqrt.(1 .+ input.Aimp ./ input.Ai) .* Tauimpi)

    mimp = input.Aimp * mp
    wcz = input.Zimp .* qe .* input.B0 ./ mimp
    rhoLz2 = (2 .* qe .* input.Ti ./ mimp) ./ wcz.^2

    fdps = (0.711 .+ 2.08e-3 .* input.Zimp.^1.26) ./ (1 .+ 1.06e-11 .* input.Zimp.^5.78)
    eps2 = eps.^2

    if input.rotation_model == 0 
        CgeoG = @. 2*eps2*(0.96*(1 - 0.54*ft^4.5)) #CG0
        CgeoU = 0.
        
        CclG = @. 1 + 2*eps2
        
        # deltan = zeros(nr)       # horizontal asymmetry of impurity density
        # Deltan = zeros(nr)       # vertical asymmetry of impurity density
        # nn     = ones((nr, nth)) # poloidal distribution of the impurity density
        
        e0imp = 1.0

    elseif input.rotation_model == 1
        UU = @. -(input.Zimp/input.Zi)*(C0z + ki)*(grad_ln_Ti*amin)
        GG = (grad_ln_Nimp .* amin) .- (input.Zimp ./ input.Zi) .* (grad_ln_Ni .* amin) .+ (1 .+ (input.Zimp ./ input.Zi) .* (C0z .- 1.0)) .* (grad_ln_Ti .* amin)
             
        AsymPhi, AsymN = polasym_input(input.rho, eps, input.Zeff, input.Zi, input.Te_Ti, Machi2, input.fH, input.bC, input.TperpTpar_axis, input.sigH) # need to translate this function
        if input.full_geom 
        
            b2 = BV^2/fluxavg(BV^2, JV)[:,None] # need to translate fluxavg 
    
            deltan, Deltan, nn, b2sNNavg, NV = asymmetry_iterative(regulopt, nr, theta, GG, UU, Ai, Aimp, Zi, Zimp, Te_Ti, Machi2, R0, nuz, BV, RV, JV, FV, dpsidx, AsymPhi, AsymN, b2, nat_asym)
                                                          

            b2snnavg = fluxavg(b2/nn, JV)
            nnsb2avg = fluxavg(nn/b2, JV)
                        
            CgeoG = nnsb2avg - 1/b2snnavg
            CgeoU = fluxavg(nn/NV, JV) - b2sNNavg/b2snnavg
            CclG  = nnsb2avg
            
        else
            
            deltaM = 2*(input.Aimp./input.Ai).*Machi2.*eps
            
            dminphia = input.Zimp.*(input.Te_Ti).*AsymPhi[1]
            dmajphia = input.Zimp.*(input.Te_Ti).*AsymPhi[2]
            
            dNH =  AsymN[1] 
            dNV =  AsymN[2]
            
            theta = ones(20) # fix this
            nat_asym = true # expose this

            deltan, Deltan, nn = asymmetry_analytical(input.rho, theta, GG, UU, eps, input.invaspct, input.qmag, nuz/wcz, deltaM, input.Ai, input.Aimp, input.Zi, input.Zimp, dNH, dNV, dminphia, dmajphia, nat_asym)
            
            dD2 = 0.5*(deltan^2 + Deltan^2)
        
            CgeoG = @. 2.0*eps*deltan + dD2 + 2.0*eps2
            CgeoU = @. -(eps*(dNH - deltan) - dD2 + 0.5*(deltan*dNH + Deltan*dNV))
            CclG = @. 1.0 + eps*deltan + 2*eps2
            
        e0imp = 1.0

        end
    
        
        
    elseif input.rotation_model == 2
        
        CG0 = @. 2 * eps2 * (0.96 * (1 - 0.54 * ft^4.5))
        CgeoG = fG .* CG0
        CgeoU = fU .* CG0
        CclG = 1 .+ 2 .* eps2
        
        if input.fsaout
            # e0imp = fluxavg(exp(Mzstar[...,None]^2*((RV^2 - (RV[:,0]^2)[:,None])/R0^2)), JV)
            JV = input.B0 ./ (input.qmag .* input.RV)
            e0imp = fluxavg(exp.(reshape(Mzstar, size(Mzstar)..., 1).^2).*((input.RV.^2 .- input.RV[:,1].^2) ./ input.R0^2), JV)
            @show e0imp

            # e0imp = ones(nr) # fix this 
        else
            e0imp = ones(nr)
        end
        
        
        deltan = 1/e0imp .- 1
        Deltan = zeros(nr)
        # nn     = 1 + deltan[...,None]*np.cos(theta) # poloidal distribution of the impurity density
    end

    Dz_PS = @. fdps * input.qmag^2 * rhoLz2 * nuz * (2 * eps2 * (0.96 * (1 - 0.54 * ft^4.5))) / (2 * eps2)
    Kz_PS = @. (input.Zimp / input.Zi) * Dz_PS
    Hz_PS = @. (-(1 + (input.Zimp / input.Zi) * (C0z - 1)) * Dz_PS)

    Dz_BP = @. (1.5 * (qe * input.Ti) / (input.Zimp ^ 2 * qe ^ 2 * input.FV ^ 2 * input.Nimp)) * (1 / (1 / K11i + 1 / K11z)) / e0imp
    Kz_BP = @. (input.Zimp / input.Zi) * Dz_BP
    Hz_BP = @. fhbp * ((input.Zimp / input.Zi) * (K12i / K11i - fv) - (K12z / K11z - fv)) * Dz_BP

    Dz_CL = @. (CclG * 2 * eps2 / CgeoG) * Dz_PS / (2 * input.qmag^2) # CL diffusion coefficient [m^2/s]
    Kz_CL = @. (input.Zimp / input.Zi) * Dz_CL # CL coefficient of the main ion density gradient [m^2/s]
    Hz_CL = @. -(1.0 + (input.Zimp/input.Zi)*(C0z - 1.0)) * Dz_CL # CL coefficient of the main ion temperature gradient [m^2/s]

    output = FACIToutput()

    output.Dz = @. Dz_PS + Dz_BP + Dz_CL
    output.Kz = @. Kz_PS + Kz_BP + Kz_CL
    output.Hz = @. Hz_PS + Hz_BP + Hz_CL

    output.Vrz = -output.Dz .* grad_ln_Nimp .+ output.Kz .* grad_ln_Ni .+ output.Hz .* grad_ln_Ti
    output.Vconv = output.Kz .* grad_ln_Ni .+ output.Hz .* grad_ln_Ti
    output.Flux_z = output.Vrz .* input.Nimp

    return output
end

########################
# Supporting functions #
########################

function ftrap(eps)
    return 1 .- (1 .- eps).^1.5 ./ (sqrt.(1 .+ eps) .* (1 .+ 1.46 .* sqrt.(eps)))
end

function LniiNRL(Zi, Ni, Ti)
    return 23 .- log.(Zi.^3 .* sqrt.(2 .* (Ni ./ 1e6))) .+ 1.5 .* log.(Ti)
end

function collision_times(Zi, Zimp, Ni, Nimp, Ti, Ai, Aimp)
    mi = Ai * mp
    mimp = Aimp * mp
    Lnii = LniiNRL(Zi, Ni, Ti)
    Lnimpi = 23 .- log.(Zi .* Zimp .* sqrt.((Ni ./ 1e6) .* Zi.^2 .+ (Nimp ./ 1e6) .* Zimp.^2)) .+ 1.5 .* log.(Ti)
    Lnimpimp = 23 .- log.(Zimp.^3 .* sqrt.(2 .* (Nimp ./ 1e6))) .+ 1.5 .* log.(Ti)

    Tauii = (eps_pi_fac * sqrt.(mi) .* Ti.^1.5) ./ (Zi.^4 .* Ni .* Lnii)
    Tauimpi = sqrt(Aimp / Ai) .* ((Zi.^2 .* Lnii) ./ (Zimp.^2 .* Lnimpi)) .* Tauii
    Tauiimp = ((Zi.^2 .* Ni .* Lnii) ./ (Zimp.^2 .* Nimp .* Lnimpi)) .* Tauii
    Tauimpimp = sqrt(Aimp / Ai) .* ((Zi.^4 .* Ni .* Lnii) ./ (Zimp.^4 .* Nimp .* Lnimpimp)) .* Tauii

    return Tauii, Tauimpi, Tauiimp, Tauimpimp
end

function C2(alpha, g, Ai, Aimp, f1, f2)
    return 1.5 ./ (1 .+ f1 .* (Ai / Aimp)) .- (0.29 .+ 0.68 .* alpha) ./ (0.59 .+ alpha .+ (1.34 .+ f2) .* g.^(-2))
end

function ki_f(nui_star, ft, Zeff, Machi)
    c01 = 0.53 .* (1 .+ 0.646 .* Machi.^1.5)
    c02 = 1.158 .* (1 .- 0.968 .* Machi.^1.56)
    c03 = -0.98 .* (1 .- 1.228 .* Machi.^1.7)

    l1k = 5.7 .* (1 .- ft).^6.7 .+ 0.38
    l2k = (-1.52 .+ 38.4 .* (1 .- ft).^3.02 .* ft.^2.07) .+ (-1 .+ 2.6 .* ft) .* Machi.^(2.5 .* (1 .- 0.6 .* ft))
    l3k = 0.25 .+ 1.2 .* (1 .- ft).^3.65
    l5k = 0.1 .* (1 .- ft).^1.46 .* ft.^4.33 .+ 0.051 .* (1 .- 0.82 .* ft) .* Machi.^2.5
    l4k = 0.8 .+ (1.25 .* ft .+ 0.585) .* Machi
    l6k = (-0.05 .+ 1.95 .* ft.^2.5) ./ (1 .+ 2.55 .* ft.^17) .- ((0.217 .+ 14.57 .* ft.^6.3) ./ (1 .+ 5.62 .* ft.^5.72)) .* Machi.^1.5

    ki0 = -(c01 .+ 0.055 .* (Zeff .- 1)) .* (1 .- ft) ./ ((0.53 .+ 0.17 .* (Zeff .- 1)) .* (1 .- (c02 .- 0.065 .* (Zeff .- 1)) .* ft .- c03 .* ft.^2))

    kiv = ((ki0 .+ l1k .* Zeff .* sqrt.(ft .* nui_star) .+ l2k .* nui_star.^0.25) ./ (1 .+ l3k .* sqrt.(nui_star)) .-
          l4k .* l5k .* nui_star.^2 .* ft.^6 .+ l6k .* nui_star.^0.25) ./ (1 .+ l5k .* nui_star.^2 .* ft.^6)

    return kiv
end

function facs(Z, A, ft, Mzstar, rotation_model)
    if typeof(Z) <: AbstractArray && typeof(Mzstar) <: AbstractArray
        if size(Z) != size(Mzstar)
            dims = max.(size(Z), size(Mzstar))
            Z = ones(eltype(Z), dims) .* Z
            Mzstar = ones(eltype(Mzstar), dims) .* Mzstar
        end
    elseif typeof(Z) <: AbstractArray
        Mzstar = ones(eltype(Z), size(Z)) .* Mzstar
    elseif typeof(Mzstar) <: AbstractArray
        Z = ones(eltype(Mzstar), size(Mzstar)) .* Z
    end
    
    f1 = @. (1.74*(1-0.028*A) + 10.25/(1 + A/3.0)^2.8) - 0.423*(1-0.136*A)*ft^(5/4)
    f2 = @. (88.28389935 + 10.50852772*Z)/(1 + 0.2157175*Z^2.57338463)
    f3 = @. (-4.45719179e+06 + 2.72813878e+06*Z)/(1+5.26716920e+06*Z^8.33610108e-01)
    
    if rotation_model == 2
        f2 = @. f2 * exp(-10*Mzstar^2)
        f3 = @. f3 * (1 + (1 + 1.86e6*ft^11.07*(1-ft)^7.36)*Mzstar^4)*exp(-0.8*Mzstar^2)
    end
    
    fg1 = @. -1.4*ft^7 + 2.23*(1-0.31*ft)
    fg2 = @. 2.8*(1-0.63*ft)
    fg3 = @. 3.5*(1 - ft)/fg2
    fg4 = @. 4*ft
    fg5 = @. 0.38*ft^4
    fg6 = @. 3.95*(1 + 0.424*ft*(1 - 0.65*ft))
    
    fG = @. (1 + fg1*Mzstar^fg2)^(fg3)*(1 + 0.2*Mzstar^fg4)/(1 + fg5*Mzstar^fg6)
    
    c1f = @. 2.72*(1-0.91*ft)
    c2f = @. 2.98*(1-0.471*ft)
    c3f = @. 0.95*(1-ft)^4
    c4f = @. 4*ft
    c5f = @. 0.1314*ft^2.84 + 3.178*ft^11.4
    c6f = @. -9.38*(ft-0.5)^2 + 4.64
    
    fU = @. (c1f*Mzstar^c2f)*(1 + c3f*Mzstar^c4f)/(1 + c5f*Mzstar^c6f)
        
    y11zb  = @. (8*(1-ft)^20 - 0.66*ft^4.9 + 0.94)*((1.085e4 + 9.3e3*14/Z^(5/3))/(1 + 14/Z^(5/3)))*ft^4.6*(1-ft)/
             (1 + 9.44e3*(1 - 2e-3*Z)*ft^4.16)
    y11zp  = 1.0
    y11zps = 1.0
    
    y12zb  = @. (1.8e3 + 7.54*Z^1.8)*ft^(4.05*(1 + 0.0039*Z))*(1-ft)^(0.7*(1 + 0.015*Z))/(1 + 1276*(1 + 0.053*Z)*ft^3.6)
    y12zp  = 1.0
    y12zps = 1.0
    
    y11ib  = @. ((5.91e-5*Z + 0.812 + 0.806/Z^0.44) + (0.013*Z + 0.098 - 7.03/Z^1.04)*ft + 
             (-0.047*Z + 0.79 + 13.8/Z^1.26)*ft^2 + (0.04*Z - 0.575 - 9.64/Z^1.31)*ft^3)*
             (571.6*(1 + 0.84*Z + 7.8e-7*Z^5.15)*ft^(3.43*(1 + 0.012*Z))*
             (1-ft)^(1.2*(1 + 1.6e-8*Z^5.1))/(1 + 500*ft^(8/3)) + 1e-3)
    y11ip  = 1.0
    y11ips = @. 1/((1+(99/44^6)*Z^6))
    
    y12ib  = @. (571.6*(1 + 0.84*Z + 7.8e-7*Z^5.15)*ft^(3.43*(1 + 0.012*Z))*
             (1-ft)^(1.2*(1 + 1.6e-8*Z^5.1))/(1 + 500*ft^(8/3)) + 1e-3)
    y12ip  = 1.0
    y12ips = 1.0
    
    if rotation_model == 2
        c11f = @. 1 + 14.86*(1 - ft)^16.45 + 15.27*ft^7.4
        c12f = @. 0.77*(1 + 4.11*ft)
        c13f = @. 0.01*(1 + 359*ft^2.5 + 1078*ft^12)
        l1f = @. (1 + c11f*Mzstar^c12f)*exp(-c13f*Mzstar^2)
        
        c21f = @. 6.39*(1 -ft)^15.8 + 0.1
        c22f = @. 0.943*(1 + 3.5*ft)
        l2f  = @. (1 + c21f*Mzstar^0.5)/(1 + c22f*Mzstar^(10/3))
        
        l3f = @. 1/(1 + 2*ft*Mzstar)
        
        c1f = @. 6.13 + 28.18*ft^2.13 + 336.25*(1-ft)^11.65
        c2f = @. 0.5 + 9.55*ft^1.14*(1-ft)^1.42
        c3f = @. (0.0087 + 4.49*ft^3.48)/(1 + 0.873*ft^3.48)
        c4f = @. 3.6*(1 - 0.36*ft)
        
        l4f = @. (1 + c1f*Mzstar^c2f)/(1 + c3f*Mzstar^(c4f))
        
        c1f = @. (1-ft)^8
        c2f = @. 113.5*ft^8.46
        c3f = @. 11*(1-ft)
        
        l5f = @. (1 + c1f*c2f*Mzstar^(c3f))/(1 + c2f*Mzstar^(c3f))
        
        l6f = @. (1 + 0.035*10*Mzstar^4)/(1 + 10*Mzstar^4)
        
        l7f = @. exp(-10*Mzstar^2)
        
        y11zb  = @. y11zb * l1f
        y11zp  = @. y11zp * l1f/l2f
        y11zps = @. y11zps * l1f/(l2f*l3f)
        
        y12zb  = @. y12zb * l4f
        y12zp  = @. y12zp * l4f/l5f
        y12zps = @. y12zps * l4f/(l5f*l6f)
        
        y11ib  = @. y11ib * l1f
        y11ip  = @. y11ip * l1f
        y11ips = @. y11ips * l1f
        
        y12ib  = @. y12ib * l7f
        y12ip  = @. y12ip * l7f
        y12ips = @. y12ips * l7f
    end
    
    yy = [[y11zb, y11zp, y11zps], [y12zb, y12zp, y12zps], 
          [y11ib, y11ip, y11ips], [y12ib, y12ip, y12ips]]
    
    if rotation_model == 2
        fv = @. 1.5*exp(-10*Mzstar^2)
    else
        fv = 1.5
    end
    
    fdps = @. ((0.711 + 2.08e-3*Z^1.26)/(1 + 1.06e-11*Z^5.78))
    
    if rotation_model == 2
        fhbp = @. (0.135 + 2.647e-3*Z^1.464 + 3.478e-10*Z^5.347)*
               ((1 + (3/(1 + 1e-7*Z^6))*(1/(1 + 1.2e5*ft^12))*Mzstar)/
               (1 + (3/(1 + 1e-7*Z^6))*(1.208 - 4.46*ft + 4.394*ft^2)*Mzstar^2))
    else
        fhbp = @. (1.01579172e+00 + -1.78923911e-03*Z)/(1 + 6.60170647e-13*Z^6.66398825e+00)
    end
    
    return f1, f2, f3, fG, fU, yy, fv, fdps, fhbp
end

function KVISC(nimp, ni, ti, Ai, Aimp, Zi, Zimp, Tauii, Tauimpimp, Tauiimp, Tauimpi, eps, ft, R0, qmag, yy)
    wii = @. (sqrt(2*(qe*ti)/(Ai*mp)))/(qmag*R0)  # main ion transit frequency [1/s]
    wimpimp = @. sqrt(Ai/Aimp)*wii              # impurity transit frequency [1/s]

    fac_a_P = @. nimp*(qe*ti)*sqrt(π)/(3*wimpimp)
    fac_i_P = @. ni*(qe*ti)*sqrt(π)/(3*wii)

    K11aP = @. fac_a_P*2
    K12aP = @. fac_a_P*2*3

    K11iP = @. fac_i_P*2
    K12iP = @. fac_i_P*2*3

    r00 = 1/sqrt(2)
    r01 = 1.5/sqrt(2)
    r11 = 3.75/sqrt(2)

    xai = sqrt(Aimp/Ai)
    xia = 1/xai

    x2ai = xai^2
    x2ia = xia^2

    xfac_ai = (1+x2ai)^0.5
    xfac_ia = (1+x2ia)^0.5

    qaa00 = qii00 = 8/2^1.5
    qaa01 = qii01 = 15/2^2.5
    qaa11 = qii11 = 132.5/2^3.5

    qai00 = (3+5*x2ai)/xfac_ai^3
    qia00 = (3+5*x2ia)/xfac_ia^3

    qai01 = 1.5*(3+7*x2ai)/xfac_ai^5
    qia01 = 1.5*(3+7*x2ia)/xfac_ia^5

    qai11 = (35*x2ai^3 + 38.5*x2ai^2 + 46.25*x2ai + 12.75)/xfac_ai^7
    qia11 = (35*x2ia^3 + 38.5*x2ia^2 + 46.25*x2ia + 12.75)/xfac_ia^7

    fac_qai_PS = @. (ni*Zi^2/(nimp*Zimp^2))
    fac_qia_PS = @. 1 ./ fac_qai_PS

    qa00 = @. fac_qai_PS*qai00 + qaa00 - r00
    qi00 = @. fac_qia_PS*qia00 + qii00 - r00

    qa01 = @. fac_qai_PS*qai01 + qaa01 - r01
    qi01 = @. fac_qia_PS*qia01 + qii01 - r01

    qa11 = @. fac_qai_PS*qai11 + qaa11 - r11
    qi11 = @. fac_qia_PS*qia11 + qii11 - r11

    Qa = @. 0.4*(qa00*qa11-qa01*qa01)
    Qi = @. 0.4*(qi00*qi11-qi01*qi01)

    la11 = @. qa11/Qa
    la12 = @. 3.5*(qa11+qa01)/Qa

    li11 = @. qi11/Qi
    li12 = @. 3.5*(qi11+qi01)/Qi

    fac_imp_PS = @. nimp*(qe*ti)*Tauimpimp
    fac_ion_PS = @. ni*(qe*ti)*Tauii

    K11aPS = @. fac_imp_PS*la11
    K12aPS = @. fac_imp_PS*la12

    K11iPS = @. fac_ion_PS*li11
    K12iPS = @. fac_ion_PS*li12

    fac_B = @. (ft/(1-ft))*(2*R0^2*qmag^2/(3*eps^2))

    nuDai_int  = @. (xfac_ai + x2ai*log(xai/(1+xfac_ai)))/Tauimpi
    nuD2ai_int = @. 1/(xfac_ai*Tauimpi)

    nuDia_int  = @. (xfac_ia + x2ia*log(xia/(1+xfac_ia)))/Tauiimp
    nuD2ia_int = @. 1/(xfac_ia*Tauiimp)

    nuDaa_int  = @. (sqrt(2) + log(1/(1+sqrt(2))))/Tauimpimp
    nuD2aa_int = @. 1/(sqrt(2)*Tauimpimp)

    nuDii_int  = @. (sqrt(2) + log(1/(1+sqrt(2))))/Tauii
    nuD2ii_int = @. 1/(sqrt(2)*Tauii)

    K11aB = @. fac_B*nimp*(Aimp*mp)*(nuDai_int + nuDaa_int)
    K12aB = @. fac_B*nimp*(Aimp*mp)*(nuD2ai_int + nuD2aa_int)

    K11iB = @. fac_B*ni*(Ai*mp)*(nuDia_int + nuDii_int)
    K12iB = @. fac_B*ni*(Ai*mp)*(nuD2ia_int + nuD2ii_int)

    y11zb  = yy[1][1]
    y11zp  = yy[1][2]
    y11zps = yy[1][3]

    y12zb  = yy[2][1]
    y12zp  = yy[2][2]
    y12zps = yy[2][3]

    y11ib  = yy[3][1]
    y11ip  = yy[3][2]
    y11ips = yy[3][3]

    y12ib  = yy[4][1]
    y12ip  = yy[4][2]
    y12ips = yy[4][3]

    K11a = @. y11zb*K11aB/((1 + y11zb*K11aB/(y11zp*K11aP))*(1 + y11zp*K11aP/(y11zps*K11aPS)))
    K12a = @. y12zb*K12aB/((1 + y12zb*K12aB/(y12zp*K12aP))*(1 + y12zp*K12aP/(y12zps*K12aPS)))

    K11i = @. y11ib*K11iB/((1 + y11ib*K11iB/(y11ip*K11iP))*(1 + y11ip*K11iP/(y11ips*K11iPS)))
    K12i = @. y12ib*K12iB/((1 + y12ib*K12iB/(y12ip*K12iP))*(1 + y12ip*K12iP/(y12ips*K12iPS)))

    return K11i, K12i, K11a, K12a
end

function polasym_input(rho, eps, Zeff, Zi, Te_Ti, Machi2, fH, bC, TperpTpar_axis, sigH)
    AsymPhi = Vector{Float64}[] 
    AsymN = Vector{Float64}[]

    TperpTpar = @. (TperpTpar_axis - 1.)*exp(-(rho/sigH)^2) + 1.
    
    AsymPhi1 = eps ./ (1. .+ Zeff .* (Te_Ti)) .* (fH * (TperpTpar .- 1.) * bC ./ (bC .+ TperpTpar * (1. - bC)) .+ 2 * Machi2)
    push!(AsymPhi, AsymPhi1)

    AsymPhi2 = zeros(length(AsymPhi1))
    push!(AsymPhi, AsymPhi2)
                    
    AsymN1 = -Zi.*(Te_Ti).*AsymPhi[1] + 2*eps.*Machi2
    push!(AsymN, AsymN1)

    AsymN2 = -Zi.*(Te_Ti).*AsymPhi[2]
    push!(AsymN, AsymN2)
    
    return AsymPhi, AsymN#, TperpTpar
end

function asymmetry_analytical(rho, theta, GG, UU, eps, invaspct, qmag, nuswcz, deltaM, Ai, Aimp, Zi, Zimp, dNH, dNV, dminphia, dmajphia, nat_asym)  
    UG = 1 .+ UU ./ GG
    @show size(UG)
    
    if nat_asym
        Ae = nuswcz .* qmag.^2 ./ invaspct
    else
        Ae = 0
    end
    
    AGe = Ae*GG
    HH = 1.0
    CD0 = -eps/UG
    QQ = CD0 * (dNV ./ (eps)) .* (UG .- 1.0)
    FF = CD0 * (1 .- 0.5 * dNH .* (UG .- 1.0) / (eps) )
    KK = 1.
    
    CD = FF .- 0.5 .* (dminphia .- deltaM)
    CDV = -0.5 .* (dmajphia .+ QQ)
    RD = @. sqrt((FF + 0.5*(dminphia-deltaM))^2 + 0.25*(dmajphia - QQ)^2)
    DD = @. RD^2 + AGe^2*(RD/CD0)^2
    
    num  = @. ((AGe/CD0)^2 - 1)*(FF/(CD0) + 0.5*(dminphia-deltaM)/CD0) + (AGe/CD0)*(0.5*dNV*(UG-1.0)/(eps) - 0.5*dmajphia/CD0)
    cosa = @. RD*CD0*num/DD
    
    num  = @. 2*AGe*(FF/CD0 + 0.5*(dminphia-deltaM)/CD0)+((AGe/CD0)^2-1)*(0.5*dmajphia - 0.5*dNV*CD0*(UG-1.0)/((eps)))
    sina = @. RD*num/DD
    
    deltan = @. CD + RD*cosa
    Deltan = @. CDV + RD*sqrt(KK/HH)*sina
    nn = 1 + (deltan * transpose(cos.(theta))) .+ (Deltan * transpose(sin.(theta)))
    
    return deltan, Deltan, nn
end


