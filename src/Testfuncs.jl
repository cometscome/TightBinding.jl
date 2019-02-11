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
    return energies

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


    pls = plot_fermisurface_2D(la)
    #return pls

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
    #savefig("Fe.png")
    #return pls

    ham = hamiltonian_k(la)
    kx = 0.1
    ky = 0.2
    hamk = ham([kx,ky])
    #println(hamk)
    return hamk

end

function test_pnictides_5orbitals()
    la = set_Lattice(2,[[1,0],[0,1]])
    add_atoms!(la,[0,0])
    add_atoms!(la,[0,0])
    add_atoms!(la,[0,0])
    add_atoms!(la,[0,0])
    add_atoms!(la,[0,0])

    tmat = [
    -0.7    0 -0.4  0.2 -0.1
    -0.8    0    0    0    0
     0.8 -1.5    0    0 -0.3
       0  1.7    0    0 -0.1
    -3.0    0    0 -0.2    0
    -2.1  1.5    0    0    0
     1.3    0  0.2 -0.2    0
     1.7    0    0  0.2    0
    -2.5  1.4    0    0    0
    -2.1  3.3    0 -0.3  0.7
     1.7  0.2    0  0.2    0
     2.5    0    0  0.3    0
     1.6  1.2 -0.3 -0.3 -0.3
       0    0    0 -0.1    0
     3.1 -0.7 -0.2    0    0
    ]
    tmat = 0.1.*tmat
    imap = zeros(Int64,5,5)
    count = 0
    for μ=1:5
        for ν=μ:5
            count += 1
            imap[μ,ν] = count
        end
    end
    Is = [1,-1,-1,1,1,1,1,-1,-1,1,-1,-1,1,1,1]
    σds = [1,-1,1,1,-1,1,-1,-1,1,1,1,-1,1,-1,1]
    tmat_σy = tmat[:,:]
    tmat_σy[imap[1,2],:] = -tmat[imap[1,3],:] #(1,1)
    tmat_σy[imap[1,3],:] = -tmat[imap[1,2],:] #(1,2)
    tmat_σy[imap[1,4],:] = -tmat[imap[1,4],:] #(1,3)
    tmat_σy[imap[2,2],:] = tmat[imap[3,3],:] #(2,2)
    tmat_σy[imap[2,4],:] = tmat[imap[3,4],:] #()
    tmat_σy[imap[2,5],:] = -tmat[imap[3,5],:] #()
    tmat_σy[imap[3,3],:] = tmat[imap[2,2],:] #()
    tmat_σy[imap[3,4],:] = tmat[imap[2,4],:] #()
    tmat_σy[imap[3,5],:] = -tmat[imap[2,5],:] #()
    tmat_σy[imap[4,5],:] = -tmat[imap[4,5],:] #()

    hoppingmatrix = zeros(Float64,5,5,5,5)
    hops = [-2,-1,0,1,2]
    hopelements = [[1,0],[1,1],[2,0],[2,1],[2,2]]

    #println(tmat)

#    count = 0
    for μ = 1:5
        for ν=μ:5
            for ii=1:5
                ihop = hopelements[ii][1]
                jhop = hopelements[ii][2]
                #[a,b],[a,-b],[-a,-b],[-a,b],[b,a],[b,-a],[-b,a],[-b,-a]

                #[a,b]
                i = ihop +3
                j = jhop +3
                hoppingmatrix[μ,ν,i,j]=tmat[imap[μ,ν],ii]
                #[a,-b] = σy*[a,b] [1,1] -> [1,-1]
                if jhop != 0
                    i = ihop +3
                    j = -jhop +3
                    hoppingmatrix[μ,ν,i,j]=tmat_σy[imap[μ,ν],ii]
                end

                if μ != ν
                    #[-a,-b] = I*[a,b] [1,1] -> [-1,-1],[1,0]->[-1,0]
                    i = -ihop +3
                    j = -jhop +3
                    hoppingmatrix[μ,ν,i,j]=Is[imap[μ,ν]]*tmat[imap[μ,ν],ii]
                    #[-a,b] = I*[a,-b] = I*σy*[a,b]  #[2,0]->[-2,0]
                    if jhop != 0
                        i = -ihop +3
                        j = jhop +3
                        hoppingmatrix[μ,ν,i,j]=Is[imap[μ,ν]]*tmat_σy[imap[μ,ν],ii]
                    end
                end
                #[b,a],[b,-a],[-b,a],[-b,-a]
                if jhop != ihop
                    #[b,a] = σd*[a,b]
                    i = jhop +3
                    j = ihop +3
                    hoppingmatrix[μ,ν,i,j]=σds[imap[μ,ν]]*tmat[imap[μ,ν],ii]
                    #[-b,a] = σd*σy*[a,b]
                    if jhop != 0
                        i = -jhop +3
                        j = ihop +3
                        hoppingmatrix[μ,ν,i,j]=σds[imap[μ,ν]]*tmat_σy[imap[μ,ν],ii]
                    end

                    if μ != ν
                        #[-b,-a] = σd*[-a,-b] = σd*I*[a,b]
                        i = -jhop +3
                        j = -ihop +3
                        hoppingmatrix[μ,ν,i,j]=σds[imap[μ,ν]]*Is[imap[μ,ν]]*tmat[imap[μ,ν],ii]
                        #[b,-a] = σd*[-a,b] = σd*I*[a,-b] = σd*I*σy*[a,b]  #[2,0]->[-2,0]
                        if jhop != 0
                            i = jhop +3
                            j = -ihop +3
                            hoppingmatrix[μ,ν,i,j]=σds[imap[μ,ν]]*Is[imap[μ,ν]]*tmat_σy[imap[μ,ν],ii]
                        end
                    end
                end
            end


        end
    end

    for μ=1:5
        for ν=μ:5
            for i = 1:5
                ih = hops[i]
                for j = 1:5
                    jh = hops[j]
                    if hoppingmatrix[μ,ν,i,j] != 0.0
                        #println(ih,"\t",jh,"\t",hoppingmatrix[μ,ν,i,j] )
                        add_hoppings!(la,hoppingmatrix[μ,ν,i,j],μ,ν,[ih,jh])
                    end
                end
            end
        end
    end

    onsite = [10.75,10.96,10.96,11.12,10.62]
    set_onsite!(la,onsite)

    set_μ!(la,10.96) #set the chemical potential
#    add_diagonals!(la,[-μ,-μ])
#    ham = hamiltonian_k(la)
    #println(ham([0,0]))


    #pls = plot_fermisurface_2D(la)
    #return pls

    nk = 100
    klines = set_Klines()
    kmin = [0,0]
    kmax = [π,0]
    add_Kpoints!(klines,kmin,kmax,"(0,0)","(pi,0)",nk=nk)

    kmin = [π,0]
    kmax = [π,π]
    add_Kpoints!(klines,kmin,kmax,"(pi,0)","(pi,pi)",nk=nk)

    kmin = [π,π]
    kmax = [0,0]
    add_Kpoints!(klines,kmin,kmax,"(pi,pi)","(0,0)",nk=nk)

    pls = calc_band_plot(klines,la)
    #savefig("Fe5band.png")
    #return pls

    ham = hamiltonian_k(la)
    kx = 0.1
    ky = 0.2
    hamk = ham([kx,ky])
    #println(hamk)
    return hamk

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
    #return pls


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

    ham = hamiltonian_k(la2)
    kx = 0.1
    ky = 0.2
    hamk = ham([kx,ky])
    println(hamk)
    return hamk


end


function test_surface()
    Ax = 1
    Ay = 1
    m2x = 1
    m2y = m2x
    m0 = -2*m2x
    m(k) = m0 + 2m2x*(1-cos(k[1]))+2m2y*(1-cos(k[2]))
    Hk(k) = Ax*sin(k[1]).*σx +  Ay*sin(k[2]).*σy + m(k).*σz
    norb = 2 #The size of the matrix
    hamiltonian = surfaceHamiltonian(Hk,norb,numhop=3,L=32,kpara="kx",BC="OBC")

    nkx = 100
    kxs = range(-π,stop=π ,length=nkx)
    mat_e = zeros(Float64,nkx,32*2)
    for i=1:nkx
        kx = kxs[i]
        mat_h = hamiltonian(kx)
    #println(mat_h)
    
        e,v = eigen(Matrix(mat_h))
        #println(e)
        mat_e[i,:] = real.(e[:])
    end
    
    plot(kxs,mat_e,labels="")
    return mat_e
    
    #
    
    return
end