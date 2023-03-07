module TightBinding

include("./TBfromk2r.jl")
include("./wannier_module.jl")
include("./lattice/kpoints.jl")
include("./lattice/lattice.jl")
include("./lattice/wannier.jl")
include("./band/band.jl")
include("./kspace/k2real.jl")
#using Plots
#using .Plotfuncs
using .Wannier
using .TBfromk2r
using LinearAlgebra
using Requires
using ProgressBars
export set_Lattice,
    add_atoms!,
    add_hoppings!,
    add_diagonals!,
    hamiltonian_k,
    dispersion,
    get_position,
    calc_band,
    get_position_kspace,
    hamiltonian_k_1d,
    set_Klines,
    show_neighbors,
    add_Kpoints!,
    set_onsite!,
    set_μ!,
    write_hr,
    make_supercell,
    read_wannier
#,calc_band_plot,plotfuncs,
#plot_lattice_2d,calc_band_plot_finite,plot_DOS

export surfaceHamiltonian, σx, σy, σz, σ0
export get_DOS



function __init__()
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
        include("./Plotfuncs.jl")
        export Plotfuncs, plot_lattice_2d, calc_band_plot, plot_DOS, calc_band_plot_finite #From Plotfuncs.jl
        export plot_fermisurface_2D
    end
end





include("./Testfuncs.jl")


end # module
