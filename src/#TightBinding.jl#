module TightBinding
    include("./Plotfuncs.jl")
    using Plots
    using .Plotfuncs
    using LinearAlgebra
    export set_Lattice,add_atoms!,add_hoppings!,add_diagonals!,hamiltonian_k,
    dispersion,get_position,calc_band,get_position_kspace,hamiltonian_k_1d,
    set_Klines,show_neighbors,add_Kpoints!,set_onsite!,set_μ!
    #,calc_band_plot,plotfuncs,
    #plot_lattice_2d,calc_band_plot_finite,plot_DOS
    export Plotfuncs,plot_lattice_2d,calc_band_plot,plot_DOS,calc_band_plot_finite #From Plotfuncs.jl

    struct Hopping
        amplitude
        ijatoms::Array{Int64,1}
        ijpositions::Array{Float64,1}
    end

    struct Kpoints
        kmin::Array{<:Number,1}
        kmax::Array{<:Number,1}
        nk::Int64
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
        numatoms::Int64
        positions::Array{Array{Float64,1},1}
        numhopps::Int64
        hoppings::Array{Hopping,1}
        diagonals::Array{Float64,1}
        rvectors::Array{Array{Float64,1},1} #reciplocal lattice vectors
        μ::Float64
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
        μ = 0.0

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


        lattice = Lattice(dim,vectors,numatoms,positions,numhopps,hoppings,diagonals,rvectors,μ)
        return lattice
    end

    function add_atoms!(lattice,position::Array{<:Number,1})
        lattice.numatoms += 1
        push!(lattice.positions,position)
        add_diagonals!(lattice,0.0)
    end

    function add_atoms!(lattice,positions::Array{<:Array})
        for i=1:length(positions)
            lattice.numatoms += 1
            push!(lattice.positions,positions[i])
            add_diagonals!(lattice,0.0)
        end
    end

    function set_μ!(lattice,μ)
        lattice.μ = μ
    end

    function set_onsite!(lattice,onsite::Array{<:Number,1})
        lattice.diagonals = onsite[:]
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
            position += vec[i]*lattice.rvectors[i]
        end
        return position
    end


    function get_hopindex(lattice)
        hopindices = []
        for ihop=1:lattice.numhopps
            hopping = lattice.hoppings[ihop]
            i = hopping.ijatoms[1]
            ri = lattice.positions[i]
            rhop = ri+hopping.ijpositions
            iabhop = round.(Int,rhop).-1
            #j = hopping.ijatoms[2]
            push!(hopindices,iabhop)
        end
        return hopindices
    end

    function hamiltonian_k_1d(lattice,direction;nsites = 20,periodic=true)
        numatoms = lattice.numatoms
        diagonals = zeros(Float64,numatoms)
        vectors = lattice.vectors
        hopindices = get_hopindex(lattice)
        numhopps = lattice.numhopps
        μ = lattice.μ

        for i=1:length(lattice.diagonals)
            diagonals[i] = lattice.diagonals[i]-μ
        end

        function calc_ham(k)
            N = numatoms*nsites
            realpart_ham = zeros(Float64,N,N)
            imagpart_ham = zeros(Float64,N,N)
            for isite=1:nsites
                for i=1:numatoms
                    ii = (isite-1)*numatoms+i
                    realpart_ham[ii,ii] = diagonals[i]
                end
                for ihop=1:numhopps
                    δ = hopindices[ihop][direction]
                    hopping = lattice.hoppings[ihop]
                    i = hopping.ijatoms[1]
                    ii = (isite-1)*numatoms+i
                    jsite = isite + δ
                    if periodic
                        while jsite > nsites || jsite < 1
                            jsite += ifelse(jsite < 1,nsites,0)
                            jsite += ifelse(jsite > nsites,-nsites,0)
                        end
                    end
                    if 1 <= jsite <= nsites
                        j= hopping.ijatoms[2]
                        jj = (jsite-1)*numatoms+j

                        ak = 0.0
                        ik = 0
                        hop = hopping.ijpositions
                        for idim = 1:lattice.dim
                            if idim != direction
                                ik += 1
                                ak += k[ik]*hop[idim]
                            end
                        end

                        ampcos = hopping.amplitude*cos(ak)
                        realpart_ham[ii,jj] +=  ampcos
                        realpart_ham[jj,ii] +=  ampcos
                        ampsin = hopping.amplitude*sin(ak)
                        imagpart_ham[ii,jj] +=  ampsin
                        imagpart_ham[jj,ii] +=  -ampsin

                    end
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



    function hamiltonian_k(lattice)
        numatoms = lattice.numatoms
        diagonals = zeros(Float64,numatoms)
        vectors = lattice.vectors
        μ = lattice.μ
        for i=1:length(lattice.diagonals)
            diagonals[i] = lattice.diagonals[i]-μ
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

    function show_neighbors(lattice;maxn = 1,cutoff=2.0)
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



    include("./Testfuncs.jl")


end # module
