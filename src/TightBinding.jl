module TightBinding
    using LinearAlgebra
    using Plots
    export set_Lattice,add_atoms!,add_hoppings!,add_diagonals!,hamiltonian_k,
    dispersion

    struct Hopping
        amplitude
        ijatoms::Array{Int64,1}
        ijpositions::Array{Int64,1}
    end

    mutable struct Lattice
        dim::Int8
        vectors::Array{Array{Float64,1},1}
        numatoms::Int8
        positions::Array{Array{Float64,1},1}
        numhopps::Int8
        hoppings::Array{Hopping,1}
        diagonals::Array{Float64,1}
    end

    function set_Lattice(dim,vectors)
        positions = Array{Float64,1}[]
        hoppings = Array{Int64,1}[]
        numatoms = 0
        numhopps=0
        diagonals = Float64[]
        lattice = Lattice(dim,vectors,numatoms,positions,numhopps,hoppings,diagonals)
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

    function hamiltonian_k(lattice)
        numatoms = lattice.numatoms
        diagonals = zeros(Float64,numatoms)
        vectors = lattice.vectors
        for i=1:length(lattice.diagonals)
            diagonals[i] = lattice.diagonals[i]
        end

        function calc_ham(k)
#        calc_ham(k) = begin
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
                        vec_a += hop.*lattice.vectors[i]
                    end
                    ak = 0.0
                    for idim = 1:lattice.dim
                        ak += k[idim]*vec_a[idim]
                    end
                    realpart_ham[i,j] +=  hopping.amplitude*cos(ak)
                    realpart_ham[j,i] +=  realpart_ham[i,j]
                    imagpart_ham[i,j] +=  hopping.amplitude*sin(ak)
                    imagpart_ham[j,i] +=  -imagpart_ham[i,j]

                end
            end
            if sum(imagpart_ham) == 0.0
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

    function calc_band(kmin,kmax,nk,lattice,ham)
        n = lattice.numatoms
        dim = lattice.dim
        energies = zeros(Float64,n,nk)
        vec_k = zeros(Float64,dim,nk)
        for idim=1:dim
            k = range(kmin[idim],length=nk,stop=kmax[idim])
            vec_k[idim,:] = k[:]
        end
        for ik = 1:nk
            energies[:,ik] = dispersion(dim,n,ham,vec_k[ik])
        end

        return vec_k,energies
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
        #k = [0.0]
        #dispersion(la1.dim,la1.numatoms,ham1,k)
        pls = plot(vec_k[1,:],energies[1,:],marker=:circle,label=["1D"])
        savefig("1Denergy.png")

    end

greet() = print("Hello World!")

end # module
