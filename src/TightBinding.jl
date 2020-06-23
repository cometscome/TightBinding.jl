module TightBinding
    
    include("./TBfromk2r.jl")
    #using Plots
    #using .Plotfuncs
    using .TBfromk2r
    using LinearAlgebra
    using Requires
    export set_Lattice,add_atoms!,add_hoppings!,add_diagonals!,hamiltonian_k,
    dispersion,get_position,calc_band,get_position_kspace,hamiltonian_k_1d,
    set_Klines,show_neighbors,add_Kpoints!,set_onsite!,set_μ!,
    write_hr,make_supercell
    #,calc_band_plot,plotfuncs,
    #plot_lattice_2d,calc_band_plot_finite,plot_DOS
    
    export surfaceHamiltonian,σx,σy,σz,σ0
    export get_DOS

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

    function __init__()
        @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin 
            include("./Plotfuncs.jl")
            export Plotfuncs,plot_lattice_2d,calc_band_plot,plot_DOS,calc_band_plot_finite #From Plotfuncs.jl
            export plot_fermisurface_2D
        end
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

    function Base.display(lattice::Lattice)
        println("Dimension: ", lattice.dim)
        println("Unit vectors ")
        for id=1:lattice.dim
            println(lattice.vectors[id])
        end
        println("Reciplocal lattice vectors ")
        for id=1:lattice.dim
            println(lattice.rvectors[id])
        end
        println("Chemical potential: ", lattice.μ)
        println("num. of atoms ",lattice.numatoms)
        for i=1:lattice.numatoms
            println("$(i)-th atom ",lattice.positions[i], " onsite energy: ", lattice.diagonals[i])
        end
        
        println("num. of hoppings ",lattice.numhopps)
        if lattice.numhopps != 0
            println("iband jband hopping amplitude")
            for i=1:lattice.numhopps
                hopping = lattice.hoppings[i]
                iband = hopping.ijatoms[1]
                jband = hopping.ijatoms[2]
                println("$iband $jband $(hopping.ijpositions) $(hopping.amplitude)")
            end
        end
        
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

    function find_orbital(lattice::Lattice,position)
        position_sp  =position[:]
        for id=1:lattice.dim
            while position_sp[id] >= 1
                position_sp[id] -= 1
            end
            while position_sp[id] < 0
                position_sp[id] += 1
            end
        end

        for iatom = 1:lattice.numatoms
            if norm(lattice.positions[iatom] - position_sp) < 1e-15
                return iatom
            end
        end
        println(position_sp)
        error("Orbital index was not found")
    end

    function make_supercell(lattice::Lattice,supercell)
        @assert lattice.dim == length(supercell)
        numatoms = lattice.numatoms
        
        vectors = deepcopy(lattice.vectors)
        for id = 1:lattice.dim
            vectors[id][:] = supercell[id]*vectors[id][:]
        end
        
        lattice_super = set_Lattice(lattice.dim,vectors)
        diagonals = zeros(Float64,numatoms*prod(supercell))

        set_μ!(lattice_super,lattice.μ)
        
        
        if lattice.dim == 1
            icount = 0
            for i1 = 1:supercell[1]
                for iatom =1:numatoms
                    position = lattice.positions[iatom][1] + (i1-1) #*lattice.vectors[1][1]
                    position /= supercell[1]
                    add_atoms!(lattice_super,[position])
#                    icount += 1
                    icount = (i1-1)*numatoms + iatom
                    diagonals[icount] = lattice.diagonals[iatom]
                end
            end

            set_onsite!(lattice_super,diagonals)

            #=
            for iatom=1:lattice_super.numatoms
                jatom =find_orbital(lattice_super,lattice_super.positions[iatom])
                println("$jatom $iatom")
            end
            =#

            for δ=1:lattice.numhopps
                hopping = lattice.hoppings[δ] #Type:
                iband=hopping.ijatoms[1]
                jband=hopping.ijatoms[2]
                hopposition = hopping.ijpositions
                hopposition_sp = hopposition[:]/supercell[1]
                v = hopping.amplitude
                #println("$iband $jband $hopposition $hopposition_sp")
                for i1 = 1:supercell[1]
                    iiband = (i1-1)*numatoms + iband
                    atomposition_i = lattice_super.positions[iiband][:]
                    iatom = find_orbital(lattice_super,atomposition_i)
                    atomposition_j = atomposition_i[:] + hopposition_sp[:]
                    jatom =find_orbital(lattice_super,atomposition_j)
                    #println("$iatom $jatom $atomposition_i $atomposition_j")
                    add_hoppings!(lattice_super,v,iatom,jatom,hopposition_sp)
                end
            end

                

        elseif lattice.dim == 2

            for i1 = 1:supercell[1]
                for i2 = 1:supercell[2]
                    for iatom =1:numatoms
                        position_x = (lattice.positions[iatom][1] + (i1-1))/supercell[1]
                        position_y = (lattice.positions[iatom][2] + (i2-1))/supercell[2]

                        add_atoms!(lattice_super,[position_x,position_y])
    #                    icount += 1
                        icount = ((i2-1)*supercell[1] + (i1-1))*numatoms + iatom
                        diagonals[icount] = lattice.diagonals[iatom]
                    end
                end
            end
            set_onsite!(lattice_super,diagonals)

            for δ=1:lattice.numhopps
                hopping = lattice.hoppings[δ] #Type:
                iband=hopping.ijatoms[1]
                jband=hopping.ijatoms[2]
                hopposition = deepcopy(hopping.ijpositions)
                #println(hopposition)
                hopposition_sp = hopposition[:]
                hopposition_sp[1] /= supercell[1]
                hopposition_sp[2] /= supercell[2]

                v = hopping.amplitude
                #println("$iband $jband $hopposition $hopposition_sp")
                for i1 = 1:supercell[1]
                    for i2 = 1:supercell[2]
                        iiband = ((i2-1)*supercell[1]+(i1-1))*numatoms + iband
                        atomposition_i = lattice_super.positions[iiband][:]
                        iatom = find_orbital(lattice_super,atomposition_i)
                        atomposition_j = atomposition_i[:] + hopposition_sp[:]
                        jatom =find_orbital(lattice_super,atomposition_j)
                        #println("$iatom $jatom $atomposition_i $atomposition_j")
                        add_hoppings!(lattice_super,v,iatom,jatom,hopposition_sp)
                        positionindex = get_positionindex(lattice_super,atomposition_j)
                        #println("positionindex $positionindex")
                    end
                end
            end



        elseif lattice.dim == 3
            for i1 = 1:supercell[1]
                for i2 = 1:supercell[2]
                    for i3 = 1:supercell[3]
                        for iatom =1:numatoms
                            position_x = (lattice.positions[iatom][1] + (i1-1))/supercell[1]
                            position_y = (lattice.positions[iatom][2] + (i2-1))/supercell[2]
                            position_z = (lattice.positions[iatom][3] + (i3-1))/supercell[3]

                            add_atoms!(lattice_super,[position_x,position_y,position_z])
        #                    icount += 1
                            icount = (((i3-1)*supercell[2]+i2-1)*supercell[1] + (i1-1))*numatoms + iatom
                            diagonals[icount] = lattice.diagonals[iatom]
                        end
                    end
                end
            end
            set_onsite!(lattice_super,diagonals)

            for δ=1:lattice.numhopps
                hopping = lattice.hoppings[δ] #Type:
                iband=hopping.ijatoms[1]
                jband=hopping.ijatoms[2]
                hopposition = deepcopy(hopping.ijpositions)
                #println(hopposition)
                hopposition_sp = hopposition[:]
                hopposition_sp[1] /= supercell[1]
                hopposition_sp[2] /= supercell[2]
                hopposition_sp[3] /= supercell[3]

                v = hopping.amplitude
                #println("$iband $jband $hopposition $hopposition_sp")
                for i1 = 1:supercell[1]
                    for i2 = 1:supercell[2]
                        for i3 = 1:supercell[3]
                            iiband = (((i3-1)*supercell[2]+i2-1)*supercell[1]+(i1-1))*numatoms + iband
                            atomposition_i = lattice_super.positions[iiband][:]
                            iatom = find_orbital(lattice_super,atomposition_i)
                            atomposition_j = atomposition_i[:] + hopposition_sp[:]
                            jatom =find_orbital(lattice_super,atomposition_j)
                            #println("$iatom $jatom $atomposition_i $atomposition_j")
                            add_hoppings!(lattice_super,v,iatom,jatom,hopposition_sp)
                        end
                    end
                end
            end
        end

        return lattice_super
    end

    function get_positionindex(lattice::Lattice,position)
        position_sp  =position[:]
        positionindex = zeros(Int64,lattice.dim)
        for id=1:lattice.dim
            while position_sp[id] >= 1
                positionindex[id] += 1
                position_sp[id] -= 1
            end
            while position_sp[id] < 0
                positionindex[id] -= 1
                position_sp[id] += 1
            end
        end
        return positionindex

    end

    function write_hr(lattice::Lattice;filename="wannier90_hr.dat",maxhop = 10)
        numatoms = lattice.numatoms
        
        diagonals = zeros(Float64,numatoms)
        vectors = lattice.vectors
        μ = lattice.μ
        for i=1:length(lattice.diagonals)
            diagonals[i] = lattice.diagonals[i]-μ
        end
        hopdata = Array{Int64,1}[]
        valdata = Array{Float64,1}[]
        for i=1:numatoms
            v = diagonals[i]
            iband = i
            jband = i
            hop = [0,0,0,iband,jband]
            val = [real(v),-imag(v)] 
            push!(hopdata,hop)
            push!(valdata,val)
        end
        if lattice.numhopps != 0
            for δ=1:lattice.numhopps
                hopping = lattice.hoppings[δ] #Type:
                iband=hopping.ijatoms[1]
                jband=hopping.ijatoms[2]
                hopposition = hopping.ijpositions[:]

                atomposition_i = lattice.positions[iband][:]
                atomposition_j = atomposition_i[:] + hopposition[:]
                
                positionindex = get_positionindex(lattice,atomposition_j)


                v = hopping.amplitude
                if lattice.dim == 1
                    hop = [positionindex[1],0,0,iband,jband]
                elseif lattice.dim == 2
                    hop = [positionindex[1],positionindex[2],0,iband,jband]
                elseif lattice.dim == 3
                    hop = [positionindex[1],positionindex[2],positionindex[3],iband,jband]
                end
                val = [real(v),-imag(v)] 
                push!(hopdata,hop)
                push!(valdata,val)

                
                if lattice.dim == 1
                    hop = [-positionindex[1],0,0,jband,iband]
                elseif lattice.dim == 2
                    hop = [-positionindex[1],-positionindex[2],0,jband,iband]
                elseif lattice.dim == 3
                    hop = [-positionindex[1],-positionindex[2],-positionindex[3],jband,iband]
                end
                val = [real(v),imag(v)] 
                push!(hopdata,hop)
                push!(valdata,val)
                
            end
        end
        count = length(hopdata)
        if lattice.dim == 1
            hopmatrix = zeros(ComplexF64,2*maxhop+1,1,1,numatoms,numatoms)
        elseif lattice.dim == 2
            hopmatrix = zeros(ComplexF64,2*maxhop+1,2*maxhop+1,1,numatoms,numatoms)
        elseif lattice.dim == 3
            hopmatrix = zeros(ComplexF64,2*maxhop+1,2*maxhop+1,2*maxhop+1,numatoms,numatoms)
        end
        
        for i=1:count
            hop = hopdata[i]
            val = valdata[i]
            ix = hop[1] + maxhop+1
            iy = hop[2]+ maxhop+1
            iz = hop[3]+ maxhop+1
            μ = hop[4]
            ν = hop[5]
            if lattice.dim == 1
                hopmatrix[ix,1,1,μ,ν] = val[1] + im*val[2]
            elseif lattice.dim == 2
                hopmatrix[ix,iy,1,μ,ν] = val[1] + im*val[2]
            elseif lattice.dim == 3
                hopmatrix[ix,iy,iz,μ,ν] = val[1] + im*val[2]
            end
        end

        hopcount = 0
        if lattice.dim == 1
            for ix=1:2*maxhop+1
                if sum(abs.(hopmatrix[ix,1,1,:,:])) != 0
                    hopcount += 1
                end
            end
        elseif lattice.dim == 2
            for ix=1:2*maxhop+1
                for iy=1:2*maxhop+1
                    if sum(abs.(hopmatrix[ix,iy,1,:,:])) != 0
                        hopcount += 1
                    end
                end
            end
        elseif lattice.dim == 3
            for ix=1:2*maxhop+1
                for iy=1:2*maxhop+1
                    for iz = 1:2*maxhop+1
                        if sum(abs.(hopmatrix[ix,iy,iz,:,:])) != 0
                            hopcount += 1
                        end
                    end
                end
            end
        end

        count = hopcount*numatoms*numatoms

        fp = open(filename,"w")
        println(fp,"generated by TightBinding.jl")
        println(fp,numatoms)
        println(fp,count)
        num = div(count,15)
        nmod = count % 15
        for j=1:num
            for i=1:15
                print(fp,1,"\t")
            end
            println(fp,"\t")
        end
        if nmod != 0
            for i=1:nmod
                print(fp,1,"\t")
            end
            println(fp,"\t")
        end

        if lattice.dim == 1
            for i1=1:2*maxhop+1
                if sum(abs.(hopmatrix[i1,1,1,:,:])) != 0
                    ix = i1 - maxhop-1
                    iy = 0
                    iz = 0
                    for μ=1:numatoms
                        for ν=1:numatoms
                            val = [real(hopmatrix[i1,1,1,μ,ν]),imag(hopmatrix[i1,1,1,μ,ν])]
                            println(fp,"$ix\t$iy\t$iz\t$μ\t$ν\t$(val[1])\t$(val[2])")
                        end
                    end
                end
            end
        elseif lattice.dim == 2
            for i1=1:2*maxhop+1
                for i2=1:2*maxhop+1
                    if sum(abs.(hopmatrix[i1,i2,1,:,:])) != 0
                        ix = i1 - maxhop-1
                        iy = i2 - maxhop-1
                        iz = 0
                        for μ=1:numatoms
                            for ν=1:numatoms
                                val = [real(hopmatrix[i1,i2,1,μ,ν]),imag(hopmatrix[i1,i2,1,μ,ν])]
                                println(fp,"$ix\t$iy\t$iz\t$μ\t$ν\t$(val[1])\t$(val[2])")
                            end
                        end
                    end
                end
            end
        elseif lattice.dim == 3
            for i1=1:2*maxhop+1
                for i2=1:2*maxhop+1
                    for i3 = 1:2*maxhop+1
                        if sum(abs.(hopmatrix[i1,i2,i3,:,:])) != 0
                            ix = i1 - maxhop-1
                            iy = i2 - maxhop-1
                            iz = i3 - maxhop-1
                            

                            for μ=1:numatoms
                                for ν=1:numatoms
                                    val = [real(hopmatrix[i1,i2,i3,μ,ν]),imag(hopmatrix[i1,i2,i3,μ,ν])]
                                    println(fp,"$ix\t$iy\t$iz\t$μ\t$ν\t$(val[1])\t$(val[2])")
                                end
                            end
                        end
                    end
                end
            end
        end
        return 


        for i=1:count
            hop = hopdata[i]
            val = valdata[i]
            ix = hop[1]
            iy = hop[2]
            iz = hop[3]
            μ = hop[4]
            ν = hop[5]

            println(fp,"$ix\t$iy\t$iz\t$μ\t$ν\t$(val[1])\t$(val[2])")
        end
        close(fp)

        #println(hopdata)
        #println(valdata)
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


    function get_DOS(lattice,nk;nbins=100)
        dim = lattice.dim
        ham = hamiltonian_k(lattice)
        n = lattice.numatoms


        if dim == 1
            energies = zeros(Float64,n,nk)
            dk = 1/(nk-1)
            for ik = 1:nk
                k = (ik-1)*dk
            end
        elseif dim == 2
            energies = zeros(Float64,n,nk,nk)
            dk1 = 1/(nk-1)
            dk2 = 1/(nk-1)
            for ik1 = 1:nk
                k1 = (ik1-1)*dk1
                for ik2 = 1:nk
                    k2 = (ik2-1)*dk2
                    kx,ky = get_position_kspace(lattice,[k1,k2])
                    energy = dispersion(dim,n,ham,[kx,ky])
                    energies[:,ik1,ik2] = energy[:]

                end
            end
        elseif dim == 3
            energies = zeros(Float64,nk,nk,nk)
            dk1 = 1/(nk-1)
            dk2 = 1/(nk-1)
            dk3 = 1/(nk-1)

            for ik1 = 1:nk
                k1 = (ik1-1)*dk1
                for ik2 = 1:nk
                    k2 = (ik2-1)*dk2
                    for ik3=1:nk
                        k3= (ik3-1)*dk3
                        kx,ky,kz = get_position_kspace(lattice,[k1,k2,k3])
                        energy = dispersion(dim,n,ham,[kx,ky,kz])
                        energies[:,ik1,ik2,ik3] = energy[:]
                    end

                end
            end
        end
        energies = vec(energies)
#            println(energies)
        #pls = histogram(energies,bins=nbins,normalize=true)
        hist = fit(Histogram,energies,nbins=nbins)

        return hist
    end        



    include("./Testfuncs.jl")


end # module
