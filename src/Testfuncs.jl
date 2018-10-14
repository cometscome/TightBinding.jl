function test_1D()
    la1 = set_Lattice(1,[[1.0]])
    add_atoms!(la1,[0,0])
    t = 1.0
    add_hoppings!(la1,-t,1,1,[1])
    ham1 = hamiltonian_k(la1)
    kmin = [-π]
    kmax = [π]
    nk = 20
    vec_k,energies = calc_band(kmin,kmax,nk,la1,ham1)
    println(energies)

    pls = plot(vec_k[1,:],energies[1,:],marker=:circle,label=["1D"])
    savefig("1Denergy.png")

end

function test_2Dsquare()
    println("Test1: ")
    println("2D Square lattice with nearest neighbor hoppings")
    la2 = set_Lattice(2,[[1,0],[0,1]])
    add_atoms!(la2,[0,0])

    show_neighbors(la2)

    t = 1.0
    add_hoppings!(la2,-t,1,1,[1,0])
    add_hoppings!(la2,-t,1,1,[0,1])
    ham2 = hamiltonian_k(la2)

    kmin = [-π,-π]
    kmax = [0.0,0.0]
    nk = 20
    vec_k,energies = calc_band(kmin,kmax,nk,la2,ham2)
    println("Energies on the line from (-π,π) to (0,0)")
    println(energies)
    return energies


    klines = set_Klines()
    kmin = [0,0]
    kmax = [π,0]
    add_Kpoints!(klines,kmin,kmax,"(0,0)","(pi,0)")

    kmin = [π,0]
    kmax = [π,π]
    add_Kpoints!(klines,kmin,kmax,"(pi,0)","(pi,pi)")

    kmin = [π,π]
    kmax = [0,0]
    add_Kpoints!(klines,kmin,kmax,"(pi,pi)","(pi,0)")

    pls = calc_band_plot(klines,la2)
    return pls

    #pls = plot(vec_k[:],energies[1,:],marker=:circle,label=["2D"])
    #savefig("2Denergy.png")

end

function test_pnictides()
    la = set_Lattice(2,[[1,0],[0,1]])
    add_atoms!(la,[0,0])
    add_atoms!(la,[0,0])
    t1 = -1.0
    t2 = 1.3
    t3 = -0.85
    t4 = t3
    μ = 1.45

    #dxz
    add_hoppings!(la,-t1,1,1,[1,0])
    add_hoppings!(la,-t2,1,1,[0,1])
    add_hoppings!(la,-t3,1,1,[1,1])
    add_hoppings!(la,-t3,1,1,[1,-1])

    #dyz
    add_hoppings!(la,-t2,2,2,[1,0])
    add_hoppings!(la,-t1,2,2,[0,1])
    add_hoppings!(la,-t3,2,2,[1,1])
    add_hoppings!(la,-t3,2,2,[1,-1])

    #between dxz and dyz
    add_hoppings!(la,-t4,1,2,[1,1])
    add_hoppings!(la,-t4,1,2,[-1,-1])
    add_hoppings!(la,t4,1,2,[1,-1])
    add_hoppings!(la,t4,1,2,[-1,1])

    set_μ!(la,μ) #set the chemical potential
#    add_diagonals!(la,[-μ,-μ])




    klines = set_Klines()
    kmin = [0,0]
    kmax = [π,0]
    add_Kpoints!(klines,kmin,kmax,"(0,0)","(pi,0)")

    kmin = [π,0]
    kmax = [π,π]
    add_Kpoints!(klines,kmin,kmax,"(pi,0)","(pi,pi)")

    kmin = [π,π]
    kmax = [0,0]
    add_Kpoints!(klines,kmin,kmax,"(pi,pi)","(0,0)")

    pls = calc_band_plot(klines,la)
    savefig("Fe.png")
    return pls

    ham = hamiltonian_k(la)
    kx = 0.1
    ky = 0.2
    hamk = ham([kx,ky])
    println(hamk)
    return

end

function test_2DGraphene()
    a1 = [sqrt(3)/2,1/2]
    a2= [0,1]
    la2 = set_Lattice(2,[a1,a2])
    add_atoms!(la2,[1/3,1/3])
    add_atoms!(la2,[2/3,2/3])

    show_neighbors(la2)


    t = 1.0
    add_hoppings!(la2,-t,1,2,[1/3,1/3])
    add_hoppings!(la2,-t,1,2,[-2/3,1/3])
    add_hoppings!(la2,-t,1,2,[1/3,-2/3])


    klines = set_Klines()
    kmin = [-π]
    kmax = [π]
    add_Kpoints!(klines,kmin,kmax,"-pi","pi")

    pls = calc_band_plot_finite(klines,la2,1,periodic=false)
    return pls


    pls = plot_lattice_2d(la2)
    #return pls

    klines = set_Klines()
    kmin = [0,0]
    kmax = [2π/sqrt(3),0]
    add_Kpoints!(klines,kmin,kmax,"G","K")

    kmin = [2π/sqrt(3),0]
    kmax = [2π/sqrt(3),2π/3]
    add_Kpoints!(klines,kmin,kmax,"K","M")

    kmin = [2π/sqrt(3),2π/3]
    kmax = [0,0]
    add_Kpoints!(klines,kmin,kmax,"M","G")

    #pls = plot_DOS(la2,100)
    #return pls

    pls2 = calc_band_plot(klines,la2)
    #println(energies)
    return pls2


end
