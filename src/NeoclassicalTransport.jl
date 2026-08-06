module NeoclassicalTransport

using IMAS
import GACODE

include("input_neo.jl")

include("models.jl")

include("neo_nn.jl")

include("hirshman_sigmar.jl")

include("chang_hinton.jl")

include("FACIT.jl")

"""
    build_serial_neo(; force=false) -> path

Build the serial (no-MPI) NEO binary from `utilities/serial_neo/` against the
user's own gacode tree and return its path. Requires a sourced GACODE
environment (`GACODE_ROOT`/`GACODE_PLATFORM` set, neo built). No-op if the
binary already exists (`force=true` rebuilds). When the package install is
read-only, builds into a per-version scratch directory instead.
"""
function build_serial_neo(; force::Bool=false)
    for var in ("GACODE_ROOT", "GACODE_PLATFORM")
        haskey(ENV, var) || error("build_serial_neo: $var is not set — source your GACODE environment first")
    end
    neo_lib = joinpath(ENV["GACODE_ROOT"], "neo", "src", "neo_lib.a")
    isfile(neo_lib) || error("build_serial_neo: $neo_lib not found — build gacode's neo first")

    srcdir = joinpath(_pkg_root(), "utilities", "serial_neo")
    builddir = srcdir
    if !_dir_writable(srcdir)
        builddir = joinpath(_model_cache_dir(), "serial_neo")
        mkpath(builddir)
        for f in ("neo_serial.f90", "Makefile")
            cp(joinpath(srcdir, f), joinpath(builddir, f); force=true)
        end
    end

    binary = joinpath(builddir, "neo_serial")
    if force || !isfile(binary)
        run(pipeline(`make -C $builddir`; stdout=devnull))
    end
    isfile(binary) || error("build_serial_neo: build did not produce $binary")
    return binary
end

"""
    use_serial_neo!(; force=false) -> path

[`build_serial_neo`](@ref) + set `ENV["NEO_EXECUTABLE"]`, so subsequent
[`run_neo`](@ref) calls invoke the serial binary directly — the one-call setup
for running full-NEO comparisons on login nodes where the `neo -e` wrapper
cannot launch.
"""
function use_serial_neo!(; force::Bool=false)
    binary = build_serial_neo(; force)
    ENV["NEO_EXECUTABLE"] = binary
    return binary
end

"""
    run_neo(input_neo::InputNEO; neo_executable=get(ENV, "NEO_EXECUTABLE", ""))

Saves input.neo file to a temporary directory, runs NEO on that directory and parses output.

By default NEO is launched through the `neo -e` GACODE wrapper, which on some systems
dispatches via `srun`/`mpirun` and therefore cannot run on a login node. Passing
`neo_executable` (or setting the `NEO_EXECUTABLE` environment variable) bypasses the
wrapper and invokes the given binary directly in the run directory — e.g. a serial
no-MPI NEO build. This path requires `GACODE_ROOT` to be set (for the input parser
`neo_parse.py`).
"""
function run_neo(input_neo::InputNEO; neo_executable::AbstractString=get(ENV, "NEO_EXECUTABLE", ""))
    folder = mktempdir()
    save_inputneo(input_neo, joinpath(folder, "input.neo"))

    if isempty(neo_executable)
        command = """
          neo -e &> command.log
          """
    else
        isfile(neo_executable) || error("run_neo: neo_executable '$neo_executable' not found")
        haskey(ENV, "GACODE_ROOT") || error("run_neo: GACODE_ROOT must be set to use neo_executable (needed for neo_parse.py)")
        # direct invocation: parse the namelist, pre-create out.neo.run (NEO
        # appends to it), then run the given binary in place
        command = """
          python "\$GACODE_ROOT/neo/bin/neo_parse.py" &> command.log
          > out.neo.run
          "$neo_executable" >> command.log 2>&1
          """
    end
    open(joinpath(folder, "command.sh"), "w") do io
        return write(io, command)
    end

    run(Cmd(`bash command.sh`; dir=folder))

    ### parse outputs ###
    # (Fortran list-directed output may use repeat tokens like `4*0.` — e.g. from
    # Cray-compiled NEO — which must be expanded, not skipped)
    tmp_fluxes = Float64[]
    open(joinpath(folder, "out.neo.transport_flux"), "r") do io
        for line in eachline(io)
            if !startswith(line, "#")
                for word in split(line)
                    m = match(r"^(\d+)\*(.+)$", word)
                    if m !== nothing
                        val = tryparse(Float64, m.captures[2])
                        val !== nothing && append!(tmp_fluxes, fill(val, parse(Int, m.captures[1])))
                    else
                        val = tryparse(Float64, word)
                        val !== nothing && push!(tmp_fluxes, val)
                    end
                end
            end
        end
    end

    loc_first_tgyro = (4 * input_neo.N_SPECIES * 2) + 1
    tgyro_fluxes = tmp_fluxes[loc_first_tgyro:end]

    # figure out indexes
    e_index = [input_neo.N_SPECIES]
    i_index = collect(1:input_neo.N_SPECIES-1)
    particle_index(index) = 2 .+ ((index .- 1) .* 4)
    energy_index(index) = 3 .+ ((index .- 1) .* 4)
    momentum_index(index) = 4 .+ ((index .- 1) .* 4)

    # sort fluxes
    electrons_energy_flux = only(tgyro_fluxes[energy_index(e_index)])
    electrons_particle_flux = only(tgyro_fluxes[particle_index(e_index)])
    ion_particle_flux = tgyro_fluxes[particle_index(i_index)]
    ion_total_energy_flux = sum(tgyro_fluxes[energy_index(i_index)])
    ion_total_momentum_flux = sum(tgyro_fluxes[momentum_index(i_index)])

    # assign fluxes to FluxSolution structure
    T = Float64
    sol = GACODE.FluxSolution{T}(electrons_energy_flux, ion_total_energy_flux, electrons_particle_flux, ion_particle_flux, ion_total_momentum_flux)

    rm(folder; force=true, recursive=true)

    return sol
end

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__; all=false, imported=false) if name != Symbol(@__MODULE__)]

end
