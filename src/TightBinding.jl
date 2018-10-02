module TightBinding
    module plotfuncs
        using Plots
        using ..TightBinding
        export plot_lattice_2d,calc_band_plot


        function plot_lattice_2d(lattice)
            dim = lattice.dim
            lw = 0.5
            ls = :dash
            colors = ["red","blue","orange","brown","yellow"]
            println("Plot the 2D lattice structure")
            if dim != 2
                println("Error! dim should be 2 to plot the lattice. dim is $dim")
                return
            end
            pls = plot(size=(500,500), legend=false)
            r0 = [0,0]
            r1 = lattice.vectors[1]
            r2 = lattice.vectors[2]

            na = 2
            nb = 2

            for ia = 0:na
                for ib=0:nb
                    rorigin = lattice.vectors[1]*ia +  lattice.vectors[2]*ib
                    plot!([r0[1]+rorigin[1], r1[1]+rorigin[1]],
                    [r0[2]+rorigin[2], r1[2]+rorigin[2]], color="red", lw=lw, ls=ls)
                    plot!([r0[1]+rorigin[1], r2[1]+rorigin[1]],
                    [r0[2]+rorigin[2], r2[2]+rorigin[2]], color="red", lw=lw, ls=ls)
                    plot!([r1[1]+rorigin[1], r1[1]+r2[1]+rorigin[1]],
                    [r1[2]+rorigin[2], r1[2]+r2[2]+rorigin[2]], color="red", lw=lw, ls=ls)
                    plot!([r2[1]+rorigin[1], r1[1]+r2[1]+rorigin[1]],
                    [r2[2]+rorigin[2], r1[2]+r2[2]+rorigin[2]], color="red", lw=lw, ls=ls)

                    for i=1:lattice.numatoms
                        x,y = get_position(lattice,lattice.positions[i])
                        plot!([x+rorigin[1]],[y+rorigin[2]],marker=:circle,color=colors[i])
                    end

                    for i=1:lattice.numhopps
                        ijpositions = lattice.hoppings[i].ijpositions
                        x,y = get_position(lattice,ijpositions)
                        ijatoms = lattice.hoppings[i].ijatoms
                        x0,y0 = get_position(lattice,lattice.positions[ijatoms[1]])
                        plot!([x0+rorigin[1], x0+x+rorigin[1]],
                        [y0+rorigin[2], y0+y+rorigin[2]], color="black", lw=1)
                    end
                end
            end
            plot!(aspect_ratio=:equal)
            return pls
        end

        function calc_band_plot(klines,lattice)
            ham = hamiltonian_k(lattice)
            numlines = klines.numlines
            dim = lattice.dim
            klength = 0.0
            vec_k = []
            numatoms = lattice.numatoms
            energies = zeros(Float64,numatoms,0)
            #plot(x,xticks=([4,5],["G","F"]))
            xticks_values = []
            xticks_labels = []

            for i=1:numlines
                push!(xticks_values,klength)
                push!(xticks_labels,klines.kpoints[i].name_start)
                kmin = klines.kpoints[i].kmin
                kmax = klines.kpoints[i].kmax
                nk = klines.kpoints[i].nk


                kmin_real = kmin[:]
                kmax_real = kmax[:]

                kdistance = sqrt(sum((kmin_real .- kmax_real).^2))
                vec_k_temp,energies_i = calc_band(kmin_real,kmax_real,nk,lattice,ham)
                vec_k_i = range(klength,length=nk,stop=klength+kdistance)
                vec_k = vcat(vec_k,vec_k_i)
                energies = hcat(energies,energies_i)
                klength += kdistance
                push!(xticks_values,klength)
                push!(xticks_labels,klines.kpoints[i].name_end)
            end

            pls = plot(vec_k[:],[energies[i,:] for i=1:numatoms],
                            xticks = (xticks_values,xticks_labels))

            return pls
        end

        function plot_DOS(lattice,nk,de)
            dim = lattice.dim
            ham = hamiltonian_k(lattice)
            n = lattice.numatoms

            if dim == 1
                energies = zeros(Float64,nk)
                dk = norm(lattice.rvectors[1])/(nk-1)
                for ik = 1:nk
                    k = (ik-1)*dk
                end
            elseif dim == 2
                energies = zeros(Float64,nk,nk)
                dk1 = norm(lattice.rvectors[1])/(nk-1)
                dk2 = norm(lattice.rvectors[2])/(nk-1)
                for ik1 = 1:nk
                    k1 = (ik1-1)*dk1
                    for ik2 = 1:nk
                        k2 = (ik2-1)*dk2
                        kx,ky = get_position_kspace(lattice,[k1,k2])
                        energy = dispersion(dim,n,ham,[kx,ky])
                        energies[ik1,ik2] = energy

                    end
                end
            elseif dim == 3
                energies = zeros(Float64,nk,nk,nk)
                dk1 = norm(lattice.rvectors[1])/(nk-1)
                dk2 = norm(lattice.rvectors[2])/(nk-1)
                dk3 = norm(lattice.rvectors[3])/(nk-1)

                for ik1 = 1:nk
                    k1 = (ik1-1)*dk1
                    for ik2 = 1:nk
                        k2 = (ik2-1)*dk2
                        for ik3=1:nk
                            k3= (ik3-1)*dk3
                            kx,ky,kz = get_position_kspace(lattice,[k1,k2,k3])
                            energy = dispersion(dim,n,ham,[kx,ky,kz])
                            energies[ik1,ik2,ik3] = energy
                        end

                    end
                end
            end
            pls = histogram(energies)

            return pls
        end

    end

    using Plots
    using .plotfuncs
    using LinearAlgebra
    export set_Lattice,add_atoms!,add_hoppings!,add_diagonals!,hamiltonian_k,
    dispersion,get_position,calc_band,get_position_kspace

    struct Hopping
        amplitude
        ijatoms::Array{Int64,1}
        ijpositions::Array{Float64,1}
    end

    struct Kpoints
        kmin::Array{<:Number,1}
        kmax::Array{<:Number,1}
        nk::Int8
        name_start
        name_end
    end

    mutable struct Klines
        numlines::Int64
        kpoints::Array{Kpoints,1}
    end

    function set_Klines()
        kpoint = Array{Kpoints,1}[]
        klines = Klines(0,kpoint)
        return klines
    end

    function add_Kpoints!(klines,kmin,kmax,name_start,name_end;nk=20)
        kpoint = Kpoints(kmin,kmax,nk,name_start,name_end)
        klines.numlines += 1
        push!(klines.kpoints,kpoint)
    end

    mutable struct Lattice
        dim::Int8
        vectors::Array{Array{Float64,1},1}
        numatoms::Int8
        positions::Array{Array{Float64,1},1}
        numhopps::Int8
        hoppings::Array{Hopping,1}
        diagonals::Array{Float64,1}
        rvectors::Array{Array{Float64,1},1} #reciplocal lattice vectors
    end

    """
        set_Lattice(dim::Integer,vectors::Array{Array{Float64,1},1})
    Initialize lattice.
    We have to call this before making Hamiltonian.
    dim: Dimension of the system.
    vector:: Primitive vectors.

    Example:

    1D system

    ```julia
    la1 = set_Lattice(1,[[1.0]])
    ```

    2D system

    ```julia
    a1 = [sqrt(3)/2,1/2]
    a2 = [0,1]
    la2 = set_Lattice(2,[a1,a2])
    ```

    3D system

    ```julia
    a1 = [1,0,0]
    a2 = [0,1,0]
    a3 = [0,0,1]
    la2 = set_Lattice(3,[a1,a2,a3])
    ```
    """
    function set_Lattice(dim,vectors)
        positions = Array{Float64,1}[]
        hoppings = Array{Hopping,1}[]
        numatoms = 0
        numhopps=0
        diagonals = Float64[]

        #making primitive vectors
        pvector_1 = zeros(Float64,3)
        pvector_2 = zeros(Float64,3)
        pvector_3 = zeros(Float64,3)
        if dim ==1
            pvector_1[1] = vectors[1][1]
            pvector_2[2] = 1.0 # 0 1 0
            pvector_3[3] = 1.0 # 0 0 1
        elseif dim == 2
            pvector_1[1:2] = vectors[1][1:2]
            pvector_2[1:2] = vectors[2][1:2]
            pvector_3[3] = 1.0 # 0 0 1
        elseif dim == 3
            pvector_1[1:3] = vectors[1][1:3]
            pvector_2[1:3] = vectors[2][1:3]
            pvector_3[1:3] = vectors[3][1:3]
        end
        #making reciplocal lattice vectors
        area = dot(pvector_1,cross(pvector_2,pvector_3))
        rvector_1 = 2π*cross(pvector_2,pvector_3)/area
        rvector_2 = 2π*cross(pvector_3,pvector_1)/area
        rvector_3 = 2π*cross(pvector_1,pvector_2)/area

        if dim ==1
            rvectors = [[rvector_1[1]]]
        elseif dim == 2
            rvectors = [rvector_1[1:2],rvector_2[1:2]]
        elseif dim == 3
            rvectors = [rvector_1[1:3],rvector_2[1:3],rvector_3[1:3]]
        end


        lattice = Lattice(dim,vectors,numatoms,positions,numhopps,hoppings,diagonals,rvectors)
        return lattice
    end

    function add_atoms!(lattice,position::Array{<:Number,1})
        lattice.numatoms += 1
        push!(lattice.positions,position)
    end

    function add_atoms!(lattice,positions::Array{<:Array})
        for i=1:length(positions)
            lattice.numatoms += 1
            push!(lattice.positions,positions[i])
        end
    end




    function add_hoppings!(lattice,amplitude,iatom,jatom,hopping)
        hop = Hopping(amplitude,[iatom,jatom],hopping)
        lattice.numhopps += 1
        push!(lattice.hoppings,hop)
    end


    function add_diagonals!(lattice,diagonal::Number)
        push!(lattice.diagonals,diagonal)
    end

    function add_diagonals!(lattice,diagonals::Array{<:Number,1})
        for i=1:length(diagonals)
            push!(lattice.diagonals,diagonals[i])
        end
    end

    function get_position(lattice,vec)
        position = zeros(Float64,lattice.dim)
        for i=1:lattice.dim
            position += vec[i]*lattice.vectors[i]
        end
        return position
    end

    function get_position_kspace(lattice,vec)
        position = zeros(Float64,lattice.dim)
        for i=1:lattice.dim
            position += vec[i]*lattice.bvectors[i]
        end
        return position
    end





    function hamiltonian_k(lattice)
        numatoms = lattice.numatoms
        diagonals = zeros(Float64,numatoms)
        vectors = lattice.vectors
        for i=1:length(lattice.diagonals)
            diagonals[i] = lattice.diagonals[i]
        end

        function calc_ham(k)
            realpart_ham = zeros(Float64,numatoms,numatoms)
            imagpart_ham = zeros(Float64,numatoms,numatoms)
            for i=1:numatoms
                realpart_ham[i,i] = diagonals[i]
            end
            if lattice.numhopps != 0
                for δ=1:lattice.numhopps
                    hopping = lattice.hoppings[δ] #Type:
                    i=hopping.ijatoms[1]
                    j=hopping.ijatoms[2]
                    hop = hopping.ijpositions
                    vec_a = zeros(Float64,lattice.dim)
                    for i=1:length(lattice.vectors)
                        vec_a += hop[i].*lattice.vectors[i]
                    end
                    ak = 0.0
                    for idim = 1:lattice.dim
                        ak += k[idim]*vec_a[idim]
                    end

                    ampcos = hopping.amplitude*cos(ak)
                    realpart_ham[i,j] +=  ampcos
                    realpart_ham[j,i] +=  ampcos
                    ampsin = hopping.amplitude*sin(ak)
                    imagpart_ham[i,j] +=  ampsin
                    imagpart_ham[j,i] +=  -ampsin

                end
            end
            if sum(abs.(imagpart_ham)) == 0.0
                ham = realpart_ham[:,:]
            else
                ham = realpart_ham[:,:]+im*imagpart_ham[:,:]
            end

            ham

        end
        k -> calc_ham(k)
    end

    function dispersion(dim,n,ham,k)
        if n==1
            return [ham(k)[1]]
        else
            energy = eigen(ham(k)).values
            return energy
        end
    end

    function calc_band(kmin_real,kmax_real,nk,lattice,ham)
        n = lattice.numatoms
        dim = lattice.dim
        energies = zeros(Float64,n,nk)
        vec_k = zeros(Float64,dim,nk)
        for idim=1:dim
            k = range(kmin_real[idim],length=nk,stop=kmax_real[idim])
            vec_k[idim,:] = k[:]
        end
        for ik = 1:nk
            energies[:,ik] = dispersion(dim,n,ham,vec_k[:,ik])
        end

        return vec_k,energies
    end

    function show_neighbers(lattice;maxn = 1,cutoff=2.0)
        dim = lattice.dim
        println("Possible hoppings")
        if dim==1
            for i1 = -maxn:maxn
                for i=1:lattice.numatoms
                    r0 = [i1]
                    for i=1:lattice.numatoms
                        ri = lattice.positions[i]
                        for j=i:lattice.numatoms
                            rj = lattice.positions[j]+r0
                            r = rationalize.(rj-ri)
                            r_pos = get_position(lattice,r)
                            if norm(r_pos) < cutoff
                                println("($i,$j), x:$(r[1])")
                            end
                        end
                    end
                end
            end

        elseif dim == 2
            for i1 = -maxn:maxn
                for i2 = -maxn:maxn
                    r0 = [i1,i2]
                    for i=1:lattice.numatoms
                        ri = lattice.positions[i]
                        for j=i:lattice.numatoms
                            rj = lattice.positions[j]+r0
                            r = rationalize.(rj-ri)
                            r_pos = get_position(lattice,r)
                            if norm(r_pos) < cutoff
                                println("($i,$j), x:$(r[1]), y:$(r[2])")
                            end
                        end
                    end
                end
            end
        elseif dim == 3
            for i1 = -maxn:maxn
                for i2 = -maxn:maxn
                    for i3 = -maxn:maxn
                        r0 = [i1,i2,i3]
                        for i=1:lattice.numatoms
                            ri = lattice.positions[i]
                            for j=i:lattice.numatoms
                                rj = lattice.positions[j]+r0
                                r = rationalize.(rj-ri)
                                r_pos = get_position(lattice,r)
                                if norm(r_pos) < cutoff
                                    println("($i,$j), x:$(r[1]), y:$(r[2]),z:$(r[3])")
                                end
                            end
                        end

                    end
                end
            end
        end
        return

    end



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
        la2 = set_Lattice(2,[[1,0],[0,1]])
        add_atoms!(la2,[0,0])

        show_neighbers(la2)

        t = 1.0
        add_hoppings!(la2,-t,1,1,[1,0])
        add_hoppings!(la2,-t,1,1,[0,1])
        ham2 = hamiltonian_k(la2)

        kmin = [-π,-π]
        kmax = [0.0,0.0]
        nk = 20
        vec_k,energies = calc_band(kmin,kmax,nk,la2,ham2)
        println(energies)


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

    function test_2DGraphene()
        a1 = [sqrt(3)/2,1/2]
        a2= [0,1]
        la2 = set_Lattice(2,[a1,a2])
        add_atoms!(la2,[1/3,1/3])
        add_atoms!(la2,[2/3,2/3])

        show_neighbers(la2)

        t = 1.0
        add_hoppings!(la2,-t,1,2,[1/3,1/3])
        add_hoppings!(la2,-t,1,2,[-2/3,1/3])
        add_hoppings!(la2,-t,1,2,[1/3,-2/3])

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

        #pls = plot_DOS(la2,10,1)
        #return pls

        pls2 = calc_band_plot(klines,la2)
        #println(energies)
        return pls2


    end


end # module
