using SpecialFunctions

Base.@kwdef mutable struct parameter_matrices
	Z::Union{Vector{Float64}, Missing} = missing
	mass::Union{Vector{Float64}, Missing} = missing
	dens::Union{Matrix{Float64}, Missing} = missing
	temp::Union{Matrix{Float64}, Missing} = missing
	nu::Union{Matrix{Float64}, Missing} = missing
	dlnndr::Union{Matrix{Float64}, Missing} = missing
	dlntdr::Union{Matrix{Float64}, Missing} = missing
	vth::Union{Matrix{Float64}, Missing} = missing
end

function get_ion_electron_parameters(dd::IMAS.dd)
	parameter_matrices = NEO.parameter_matrices()

	eqt = dd.equilibrium.time_slice[]
	cp1d = dd.core_profiles.profiles_1d[]

	rmin = IMAS.r_min_core_profiles(cp1d, eqt)
	a = rmin[end]

	num_ions = length(cp1d.ion)

	e = IMAS.gacode_units.e # statcoul
	k = IMAS.gacode_units.k # erg/eV
	mp = IMAS.constants.m_p * 1e3 # g
	md = 3.34358e-27 * 1e3

	loglam = 24.0 .- log.(sqrt.(cp1d.electrons.density ./ 1e6) ./ (cp1d.electrons.temperature))

	n_norm = cp1d.electrons.density ./ 1e6
	t_norm = cp1d.electrons.temperature

	m_norm = 2.0
	nu_norm = sqrt.(k .* cp1d.electrons.temperature ./ md) ./ a

	Z = Vector{Float64}(undef, num_ions + 1)
	mass = Vector{Float64}(undef, num_ions + 1)

	dens = zeros(Float64, length(cp1d.ion[1].density), num_ions + 1)
	temp = zeros(Float64, length(cp1d.ion[1].temperature), num_ions + 1)
	vth = zeros(Float64, length(cp1d.ion[1].temperature), num_ions + 1)
	nu = zeros(Float64, length(cp1d.ion[1].temperature), num_ions + 1)
	dlnndr = zeros(Float64, length(cp1d.ion[1].density), num_ions + 1)
	dlntdr = zeros(Float64, length(cp1d.ion[1].temperature), num_ions + 1)

	for i in 1:num_ions
		Z[i] = cp1d.ion[i].element[1].z_n
		mass[i] = cp1d.ion[i].element[1].a ./ m_norm

		dens[:, i] = cp1d.ion[i].density ./ 1e6 ./ n_norm
		temp[:, i] = cp1d.ion[i].temperature ./ t_norm

		nu[:, i] = (@. sqrt(2) * pi * (cp1d.ion[i].density ./ 1e6) * Z[i]^4.0 * e^4.0 * loglam / sqrt(cp1d.ion[i].element[1].a * mp) / (k * cp1d.ion[i].temperature)^1.5) ./ nu_norm

		dlnndr[:, i] = -IMAS.calc_z(rmin / a, cp1d.ion[i].density)
		dlntdr[:, i] = -IMAS.calc_z(rmin / a, cp1d.ion[i].temperature)

		vth[:, i] = sqrt.((cp1d.ion[i].temperature ./ t_norm) ./ (cp1d.ion[i].element[1].a ./ m_norm))
	end

	# tacking on electron parameters at the end 
	Z[end] = -1.0
	mass[end] = 0.00054858 / m_norm # 0.00054858 is the mass of an electron in AMU
	dens[:, end] = cp1d.electrons.density ./ 1e6 ./ n_norm
	temp[:, end] = cp1d.electrons.temperature ./ t_norm

	nu[:, end] = (@. sqrt(2) * pi * (cp1d.electrons.density ./ 1e6) * (Z[end] * e)^4.0 * loglam / sqrt(0.00054858 * mp) / (k .* cp1d.electrons.temperature)^1.5) ./ nu_norm

	dlnndr[:, end] = -IMAS.calc_z(rmin / a, cp1d.electrons.density)
	dlntdr[:, end] = -IMAS.calc_z(rmin / a, cp1d.electrons.temperature)

	vth[:, end] = sqrt.((cp1d.electrons.temperature ./ t_norm) ./ (0.00054858 ./ m_norm))

	parameter_matrices.Z = Z
	parameter_matrices.mass = mass
	parameter_matrices.dens = dens
	parameter_matrices.temp = temp
	parameter_matrices.nu = nu
	parameter_matrices.dlnndr = dlnndr
	parameter_matrices.dlntdr = dlntdr
	parameter_matrices.vth = vth

	return parameter_matrices
end

function gauss_legendre(x1::Int, x2::Int, n::Int)
	eps = 1e-15

	xm = 0.5 * (x2 + x1)
	xl = 0.5 * (x2 - x1)

	x = zeros(Float64, n)
	w = zeros(Float64, n)

	# Exception for n=1 is required:
	if n == 1
		x[1] = xm
		w[1] = 2.0 * xl
		return x, w
	end

	# Roots are symmetric. We only need to find half of them.
	m = (n + 1) / 2

	# Initialize to fail first do test
	z1 = -1.0
	pp = 0.0

	# Loop over roots.
	for i in 1:m
		i = Int(i)
		z = cos(pi * (i - 0.25) / (n + 0.5))

		while abs(z - z1) > eps
			p1 = 1.0
			p2 = 0.0

			for j in 1:n
				p3 = p2
				p2 = p1
				p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j
			end

			# p1 is the Legendre polynomial. Now compute its derivative, pp.
			pp = n * (z * p1 - p2) / (z * z - 1.0)
			z1 = z
			z = z1 - p1 / pp
		end

		x[i]     = xm - xl * z
		x[n+1-i] = xm + xl * z
		w[i]     = 2.0 * xl / ((1.0 - z * z) * pp * pp)
		w[n+1-i] = w[i]
	end

	return x, w
end

function gauss_integ(xmin::Float64, xmax::Float64, func, order::Int, n_subdiv::Int)

	x0, w0 = gauss_legendre(0, 1, order)

	dx = (xmax - xmin) / n_subdiv

	n_node = n_subdiv * order

	x = zeros(Float64, n_subdiv * order)
	w = zeros(Float64, n_subdiv * order)

	for i in 1:n_subdiv
		for j in 1:order
			p = (i - 1) * order + j
			x[p] = xmin + ((i - 1) + x0[j]) * dx
			w[p] = w0[j] * dx
		end
	end

	answer = 0.0
	for p in 1:n_node
		answer = answer + w[p] * func(x[p])
	end

	return answer
end

function get_coll_freqs(ir_loc::Int, is_loc::Int, js_loc::Int, ene::Float64)
	Z = parameter_matricesHS.Z
	dens = parameter_matricesHS.dens
	vth = parameter_matricesHS.vth

	fac = (1.0 * Z[js_loc])^2 / (1.0 * Z[is_loc])^2 * (dens[:, js_loc][ir_loc] / dens[:, is_loc][ir_loc])

	xa = sqrt(ene)
	xb = xa * (vth[:, is_loc][ir_loc] / vth[:, js_loc][ir_loc])

	if xb < 1e-4
		nu_d =
			fac * (1.0 / sqrt(pi)) *
			(
				4.0 / 3.0 * (vth[:, is_loc][ir_loc] / vth[:, js_loc][ir_loc]) - 4.0 / 15.0 * (vth[:, is_loc][ir_loc] / vth[:, js_loc][ir_loc])^3 * ene + 2.0 / 35.0 * (vth[:, is_loc][ir_loc] / vth[:, js_loc][ir_loc])^5 * ene^2 -
				2.0 / 189.0 * (vth[:, is_loc][ir_loc] / vth[:, js_loc][ir_loc])^7 * ene^3
			)
	else
		Hd_coll = exp(-xb * xb) / (xb * sqrt(pi)) + (1 - (1 / (2 * xb * xb))) * erf(xb)
		Xd_coll = 1 / xa
		nu_d = fac * Hd_coll * Xd_coll
	end

	return nu_d
end

function myHSenefunc(x::Float64)
	cp1d = ddHS.core_profiles.profiles_1d[]
	eq = ddHS.equilibrium
	eqt = eq.time_slice[]
	eq1d = eqt.profiles_1d

	nu = parameter_matricesHS.nu
	vth = parameter_matricesHS.vth

	m_to_cm = IMAS.gacode_units.m_to_cm

	n_species = length(cp1d.ion) + 1

	emin = 0.0
	emax = 16.0

	xa = 2.0 / (1.0 - sqrt(emin / emax))
	xb = -(1.0 + sqrt(emin / emax)) / (1.0 - sqrt(emin / emax))
	ene = emax * ((x - xb) / xa)^2
	de = 2.0 * sqrt(emax) / xa
	val = de * exp(-ene)

	rmin = IMAS.r_min_core_profiles(cp1d, eqt)
	a = rmin[end]
	rmaj = IMAS.interp1d(eq1d.rho_tor_norm, m_to_cm * 0.5 * (eq1d.r_outboard .+ eq1d.r_inboard)).(cp1d.grid.rho_tor_norm) ./ a
	q = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.q).(cp1d.grid.rho_tor_norm)

	eps = rmin[ir_global] / (rmaj[ir_global] .* a)

	nu_d_tot = 0.0
	for js in 1:n_species
		nu_d = get_coll_freqs(ir_global, is_globalHS, js, ene)
		nu_d_tot += nu_d * nu[:, is_globalHS][ir_global]
	end

	ft_star = (3.0 * pi / 16.0) * eps^2 * vth[:, is_globalHS][ir_global] * sqrt(2.0) * ene^1.5 / (rmaj[ir_global] * abs(q[ir_global]) * nu_d_tot)
	ftrap = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.trapped_fraction).(cp1d.grid.rho_tor_norm)
	ft_fac = 1.0 / (1.0 + ftrap[ir_global] / ft_star)


	if ietypeHS == 1
		myHSenefunc = val * nu_d_tot * ene * ft_fac
	elseif ietypeHS == 2
		myHSenefunc = val * nu_d_tot * ene * ene * ft_fac
	elseif ietypeHS == 3
		myHSenefunc = val * nu_d_tot * ene * ene * ene * ft_fac
	end

	return myHSenefunc

end

function compute_HS(ir::Int, dd::IMAS.dd)
    parameter_matrices = NEO.get_ion_electron_parameters(dd)

	global ddHS = dd
	global parameter_matricesHS = parameter_matrices

	Z = parameter_matricesHS.Z
	mass = parameter_matricesHS.mass
	dens = parameter_matricesHS.dens
	temp = parameter_matricesHS.temp
	dlnndr = parameter_matricesHS.dlnndr
	dlntdr = parameter_matricesHS.dlntdr

	eqt = ddHS.equilibrium.time_slice[]
	eq1d = eqt.profiles_1d
	cp1d = ddHS.core_profiles.profiles_1d[]

	n_species = length(cp1d.ion) + 1
	m_to_cm = IMAS.gacode_units.m_to_cm

	ftrap = IMAS.interp1d(eqt.profiles_1d.rho_tor_norm, eqt.profiles_1d.trapped_fraction).(cp1d.grid.rho_tor_norm)

	rmaj = IMAS.interp1d(eq1d.rho_tor_norm, m_to_cm * 0.5 * (eq1d.r_outboard .+ eq1d.r_inboard)).(cp1d.grid.rho_tor_norm)
	rmin = IMAS.r_min_core_profiles(cp1d, eqt)
	a = rmin[end]

	rho = IMAS.rho_s(cp1d, eqt) ./ a
	q = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.q).(cp1d.grid.rho_tor_norm)

	rmin_eqt = 0.5 * (eqt.profiles_1d.r_outboard - eqt.profiles_1d.r_inboard)
	bunit_eqt = IMAS.gradient(2 * pi * rmin_eqt, eqt.profiles_1d.phi) ./ rmin_eqt

	r, z, PSI_interpolant = IMAS.ψ_interpolant(eqt.profiles_2d[1])

	Bmag2_avgs = zeros(Real, length(eqt.profiles_1d.psi))

	for (k, psi_level0) in reverse!(collect(enumerate(eqt.profiles_1d.psi)))
		r, z, PSI_interpolant = IMAS.ψ_interpolant(eqt.profiles_2d[1])
		PSI = eqt.profiles_2d[1].psi
		pr, pz, psi_level = IMAS.flux_surface(r, z, PSI, eqt.profiles_1d.psi, eqt.global_quantities.magnetic_axis.r, eqt.global_quantities.magnetic_axis.z, psi_level0, true)

		Br, Bz = IMAS.Br_Bz(PSI_interpolant, pr, pz)
		Bp2 = Br .^ 2.0 .+ Bz .^ 2.0
		Bp_abs = sqrt.(Bp2)

		dl = vcat(0.0, sqrt.(diff(pr) .^ 2 + diff(pz) .^ 2))
		ll = cumsum(dl)
		fluxexpansion = 1.0 ./ Bp_abs
		int_fluxexpansion_dl = IMAS.integrate(ll, fluxexpansion)

		Bt = eqt.profiles_1d.f[k] ./ pr
		Btot = sqrt.(Bp2 .+ Bt .^ 2)

		Bunit = bunit_eqt[k]

		Bmag2_avgs[k] = (IMAS.flxAvg((Btot ./ Bunit) .^ 2, ll, fluxexpansion, int_fluxexpansion_dl))

	end

	Bmag2_avg = IMAS.interp1d(eq1d.rho_tor_norm, Bmag2_avgs).(cp1d.grid.rho_tor_norm)

	Nx = 100 # Can be lowered to speed up calculation time
	integ_order = 8
	omega_fac = 1.0 / Bmag2_avg[ir]

	f_cp = IMAS.interp1d(eq1d.rho_tor_norm, eq1d.f).(cp1d.grid.rho_tor_norm)
	bunit_cp = IMAS.interp1d(eq1d.rho_tor_norm, IMAS.bunit(eqt)).(cp1d.grid.rho_tor_norm)
	f = f_cp .* m_to_cm ./ bunit_cp
	HS_I_div_psip = -f[ir] * q[ir] / rmin[ir]

	nux0 = zeros(Float64, n_species)
	nux2 = zeros(Float64, n_species)
	nux4 = zeros(Float64, n_species)

	for is_global in 1:n_species
		for ietype in 1:3
			global ietypeHS = ietype
			global ir_global = ir
			global is_globalHS = is_global

			eii_val = gauss_integ(-1.0, 1.0, NEO.myHSenefunc, integ_order, Nx)

			if ietype == 1
				nux0[is_global] = eii_val * 4.0 / (3.0 * sqrt(pi))
			elseif ietype == 2
				nux2[is_global] = eii_val * 4.0 / (3.0 * sqrt(pi))
			elseif ietype == 3
				nux4[is_global] = eii_val * 4.0 / (3.0 * sqrt(pi))
			end
		end
	end

	sum_nm = 0.0
	for is_global in 1:n_species
		sum_nm += mass[is_global] * dens[:, is_global][ir] * nux0[is_global]
	end

	pflux_multi = zeros(Float64, n_species)
	eflux_multi = zeros(Float64, n_species)
	for is_global in 1:n_species
		A1 = -dlnndr[:, is_global][ir] + (1.5 * dlntdr[:, is_global][ir])
		A2 = -dlntdr[:, is_global][ir]

		pflux_multi[is_global] = 0.0
		eflux_multi[is_global] = 0.0

		L_a = nux0[is_global] * omega_fac * HS_I_div_psip^2 * rho[ir]^2 * ftrap[ir] * dens[:, is_global][ir] * temp[:, is_global][ir] * mass[is_global] / (Z[is_global] * Z[is_global] * 1.0)

		for js in 1:n_species

			L_b = nux0[js] * omega_fac * HS_I_div_psip^2 * rho[ir]^2 * ftrap[ir] * dens[:, js][ir] * temp[:, js][ir] * mass[js] / (Z[js] * Z[js] * 1.0)

			if is_global == js
				L11 = -L_a * (sum_nm - mass[is_global] * dens[:, is_global][ir] * nux0[is_global]) / sum_nm
				L12 = L11 * (nux2[is_global] / nux0[is_global])
				L22 = (nux2[is_global] / nux0[is_global]) * (L12 + L_a * (nux2[is_global] / nux0[is_global] - nux4[is_global] / nux2[is_global]))
				L21 = L12

			else
				L11 = L_a * (Z[is_global] * temp[:, js][ir]) / (Z[js] * temp[:, is_global][ir]) * mass[js] * dens[:, js][ir] * nux0[js] / sum_nm
				L12 = (nux2[js] / nux0[js]) * L11
				L21 = (nux2[is_global] / nux0[is_global]) * Z[js] / (1.0 * Z[is_global]) * (mass[is_global] * dens[:, is_global][ir] * nux0[is_global] / sum_nm) * L_b
				L22 = (nux2[is_global] / nux0[is_global]) * L12

			end

			pflux_multi[is_global] += L11 * A1 + L12 * A2
			eflux_multi[is_global] += L21 * A1 + L22 * A2

		end
	end

	return pflux_multi, eflux_multi

end

function HS_to_GB(HS_solution::Tuple{Vector{Float64}, Vector{Float64}}, dd::IMAS.dd, rho::Int)
	pflux_multi, eflux_multi = HS_solution

	eqt = dd.equilibrium.time_slice[]
	cp1d = dd.core_profiles.profiles_1d[]

	rmin = IMAS.r_min_core_profiles(cp1d, eqt)
	a = rmin[end]

	neo_rho_star = (IMAS.rho_s(cp1d, eqt)./a)[rho]

	temp_e = 1.0 # electron temperature is 1 since all NEO temps are normalized against electron temp
	dens_e = 1.0 # electron density is 1 since all NEO densities are normalized against electron density

	Gamma_neo_GB = dens_e * temp_e^1.5 * neo_rho_star .^ 2
	Q_neo_GB = dens_e * temp_e^2.5 * neo_rho_star .^ 2

	pflux_norm = pflux_multi ./ Gamma_neo_GB
	eflux_norm = eflux_multi ./ Q_neo_GB

	particle_flux_e = pflux_norm[end]
	energy_flux_e = eflux_norm[end]

	energy_flux_i = sum(eflux_norm) - energy_flux_e

	HS_fluxes = IMAS.flux_solution(particle_flux_e, 0.0, energy_flux_e, energy_flux_i)
	return HS_fluxes
end
