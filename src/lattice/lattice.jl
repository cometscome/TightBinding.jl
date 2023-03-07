struct Hopping
    amplitude::ComplexF64
    ijatoms::Array{Int64,1}
    ijpositions::Array{Float64,1}
end


mutable struct Lattice
    dim::Int8
    vectors::Array{Array{Float64,1},1}
    numatoms::Int64
    positions::Array{Array{Float64,1},1}
    numhopps::Int64
    #hoppings::Array{Hopping,1}
    hoppings::Set{Hopping}
    diagonals::Array{Float64,1}
    rvectors::Array{Array{Float64,1},1} #reciplocal lattice vectors
    μ::Float64
end

function Base.display(lattice::Lattice)
    println("Dimension: ", lattice.dim)
    println("Unit vectors ")
    for id = 1:lattice.dim
        println(lattice.vectors[id])
    end
    println("Reciplocal lattice vectors ")
    for id = 1:lattice.dim
        println(lattice.rvectors[id])
    end
    println("Chemical potential: ", lattice.μ)
    println("num. of atoms ", lattice.numatoms)
    for i = 1:lattice.numatoms
        println(
            "$(i)-th atom ",
            lattice.positions[i],
            " onsite energy: ",
            lattice.diagonals[i],
        )
    end

    println("num. of hoppings ", lattice.numhopps)

    if lattice.numhopps != 0
        println("iband jband hopping amplitude")
        for i = 1:lattice.numhopps
            hopping = lattice.hoppings[i]
            iband = hopping.ijatoms[1]
            jband = hopping.ijatoms[2]
            println("$iband $jband $(hopping.ijpositions) $(hopping.amplitude)")
        end
    end


    #=

    #println("iband jband hopping amplitude")
    for iatom=1:lattice.numatoms
        for jatom=1:lattice.numatoms
            for i=1:lattice.numhopps
                hopping = lattice.hoppings[i]
                iband = hopping.ijatoms[1]
                jband = hopping.ijatoms[2]
                if iatom == iband && jatom == jband

                    if jatom < iatom
                        println("$jband $iband $(-hopping.ijpositions) $(hopping.amplitude)")
                    else
                        println("$iband $jband $(hopping.ijpositions) $(hopping.amplitude)")
                    end
                end
            end
        end
    end
    #println("num. of hoppings ",lattice.numhopps)
    =#

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
function set_Lattice(dim, vectors)
    positions = Array{Float64,1}[]
    hoppings = Set{Hopping}()
    #hoppings = Array{Hopping,1}[]
    numatoms = 0
    numhopps = 0
    diagonals = Float64[]
    μ = 0.0

    #making primitive vectors
    pvector_1 = zeros(Float64, 3)
    pvector_2 = zeros(Float64, 3)
    pvector_3 = zeros(Float64, 3)
    if dim == 1
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
    area = dot(pvector_1, cross(pvector_2, pvector_3))
    rvector_1 = 2π * cross(pvector_2, pvector_3) / area
    rvector_2 = 2π * cross(pvector_3, pvector_1) / area
    rvector_3 = 2π * cross(pvector_1, pvector_2) / area

    if dim == 1
        rvectors = [[rvector_1[1]]]
    elseif dim == 2
        rvectors = [rvector_1[1:2], rvector_2[1:2]]
    elseif dim == 3
        rvectors = [rvector_1[1:3], rvector_2[1:3], rvector_3[1:3]]
    end


    lattice = Lattice(
        dim,
        vectors,
        numatoms,
        positions,
        numhopps,
        hoppings,
        diagonals,
        rvectors,
        μ,
    )
    return lattice
end

function Lattice(dim, vectors)
    return set_Lattice(dim, vectors)
end



function add_atoms!(lattice, position::Array{<:Number,1})
    lattice.numatoms += 1
    push!(lattice.positions, position)
    add_diagonals!(lattice, 0.0)
end

function add_atoms!(lattice, positions::Array{<:Array})
    for i = 1:length(positions)
        lattice.numatoms += 1
        push!(lattice.positions, positions[i])
        add_diagonals!(lattice, 0.0)
    end
end

function set_μ!(lattice, μ)
    lattice.μ = μ
end

function set_onsite!(lattice, onsite::Array{<:Number,1})
    lattice.diagonals = onsite[:]
end

function add_hoppings!(lattice, amplitude, iatom, jatom, hopping)
    hop = Hopping(amplitude, [iatom, jatom], hopping)
    #lattice.numhopps += 1
    push!(lattice.hoppings, hop)
    lattice.numhopps = length(lattice.hoppings)
end


function add_diagonals!(lattice, diagonal::Number)
    push!(lattice.diagonals, diagonal)
end

function add_diagonals!(lattice, diagonals::Array{<:Number,1})
    for i = 1:length(diagonals)
        push!(lattice.diagonals, diagonals[i])
    end
end

function get_position(lattice, vec)
    position = zeros(Float64, lattice.dim)
    for i = 1:lattice.dim
        position += vec[i] * lattice.vectors[i]
    end
    return position
end


function get_position_kspace(lattice, vec)
    position = zeros(Float64, lattice.dim)
    for i = 1:lattice.dim
        position += vec[i] * lattice.rvectors[i]
    end
    return position
end


function get_hopindex(lattice)
    hopindices = []
    for hopping in lattice.hoppings
    #for ihop = 1:lattice.numhopps
        #hopping = lattice.hoppings[ihop]
        i = hopping.ijatoms[1]
        ri = lattice.positions[i]
        rhop = ri + hopping.ijpositions


        positionindex = get_positionindex(lattice, rhop)
        iabhop = positionindex[:]

        #            iabhop = round.(Int,rhop).-1
        #j = hopping.ijatoms[2]
        push!(hopindices, iabhop)
    end
    return hopindices
end


function hamiltonian_k_1d(lattice, direction; nsites = 20, periodic = true)
    numatoms = lattice.numatoms
    diagonals = zeros(Float64, numatoms)
    vectors = lattice.vectors
    hopindices = get_hopindex(lattice)
    numhopps = lattice.numhopps
    μ = lattice.μ

    for i = 1:length(lattice.diagonals)
        diagonals[i] = lattice.diagonals[i] - μ
    end

    function calc_ham(k)
        N = numatoms * nsites
        realpart_ham = zeros(Float64, N, N)
        imagpart_ham = zeros(Float64, N, N)
        for isite = 1:nsites
            for i = 1:numatoms
                ii = (isite - 1) * numatoms + i
                realpart_ham[ii, ii] = diagonals[i]
            end
            ihop = 0
            for hopping in lattice.hoppings
                ihop += 1
            #for ihop = 1:numhopps
                δ = hopindices[ihop][direction]
                #hopping = lattice.hoppings[ihop]
                i = hopping.ijatoms[1]
                ii = (isite - 1) * numatoms + i
                jsite = isite + δ
                if periodic
                    while jsite > nsites || jsite < 1
                        jsite += ifelse(jsite < 1, nsites, 0)
                        jsite += ifelse(jsite > nsites, -nsites, 0)
                    end
                end
                if 1 <= jsite <= nsites
                    j = hopping.ijatoms[2]
                    jj = (jsite - 1) * numatoms + j

                    ak = 0.0
                    ik = 0
                    hop = hopping.ijpositions
                    for idim = 1:lattice.dim
                        if idim != direction
                            ik += 1
                            ak += k[ik] * hop[idim]
                        end
                    end

                    ampcos = hopping.amplitude * cos(ak)
                    realpart_ham[ii, jj] += ampcos
                    realpart_ham[jj, ii] += ampcos
                    ampsin = hopping.amplitude * sin(ak)
                    imagpart_ham[ii, jj] += ampsin
                    imagpart_ham[jj, ii] += -ampsin

                end
            end
        end
        if sum(abs.(imagpart_ham)) == 0.0
            ham = realpart_ham[:, :]
        else
            ham = realpart_ham[:, :] + im * imagpart_ham[:, :]
        end
        ham
    end
    k -> calc_ham(k)
end


function find_orbital(lattice::Lattice, supercell, position, iorbital)
    position_sp = position[:]
    for id = 1:lattice.dim
        while position_sp[id] >= 1
            position_sp[id] -= 1
        end
        while position_sp[id] < 0
            position_sp[id] += 1
        end
    end

    originalnumatoms = div(lattice.numatoms, prod(supercell))

    position_ori = position_sp[:]
    positionindex = zeros(Int64, lattice.dim)

    for id = 1:lattice.dim
        position_ori[id] *= supercell[id]
        while position_ori[id] >= 1
            position_ori[id] -= 1
            positionindex[id] += 1
        end
        while position_ori[id] < 0
            position_ori[id] += 1
            positionindex[id] += 1
        end
    end
    #println("positionindex $positionindex")

    if lattice.dim == 1
        iatom = (positionindex[1]) * originalnumatoms + iorbital
    elseif lattice.dim == 2
        iatom =
            ((positionindex[2]) * supercell[1] + positionindex[1]) * originalnumatoms +
            iorbital
    elseif lattice.dim == 3
        iatom =
            (
                (positionindex[3] * supercell[2] + positionindex[2]) * supercell[1] +
                positionindex[1]
            ) * originalnumatoms + iorbital
    end

    #println("iorbital $iorbital iatom $iatom")
    return iatom


    for iatom = 1:lattice.numatoms
        if norm(lattice.positions[iatom] - position_sp) < 1e-15
            push!(iatomcandidate, iatom)
            #                return iatom
        end
    end
    if length(iatomcandidate) != 0
        return iatomcandidate[iorbital]
    end

    println(position_sp)
    error("Orbital index was not found")
end

function make_supercell(lattice::Lattice, supercell)
    @assert lattice.dim == length(supercell)
    numatoms = lattice.numatoms

    vectors = deepcopy(lattice.vectors)
    for id = 1:lattice.dim
        vectors[id][:] = supercell[id] * vectors[id][:]
    end

    lattice_super = set_Lattice(lattice.dim, vectors)
    diagonals = zeros(Float64, numatoms * prod(supercell))

    set_μ!(lattice_super, lattice.μ)


    if lattice.dim == 1
        icount = 0
        for i1 = 1:supercell[1]
            for iatom = 1:numatoms
                position = lattice.positions[iatom][1] + (i1 - 1) #*lattice.vectors[1][1]
                position /= supercell[1]
                add_atoms!(lattice_super, [position])
                #                    icount += 1
                icount = (i1 - 1) * numatoms + iatom
                diagonals[icount] = lattice.diagonals[iatom]
            end
        end

        set_onsite!(lattice_super, diagonals)

        #=
        for iatom=1:lattice_super.numatoms
            jatom =find_orbital(lattice_super,lattice_super.positions[iatom])
            println("$jatom $iatom")
        end
        =#
        δ = 0
        for hopping in lattice.hoppings
            δ += 1
        #for δ = 1:lattice.numhopps
            #hopping = lattice.hoppings[δ] #Type:
            iband = hopping.ijatoms[1]
            jband = hopping.ijatoms[2]
            hopposition = hopping.ijpositions
            hopposition_sp = hopposition[:] / supercell[1]
            v = hopping.amplitude
            #println("$iband $jband $hopposition $hopposition_sp")
            for i1 = 1:supercell[1]
                iiband = (i1 - 1) * numatoms + iband
                atomposition_i = lattice_super.positions[iiband][:]
                iatom = find_orbital(lattice_super, supercell, atomposition_i, iband)
                #                    iatom = find_orbital(lattice_super,atomposition_i)
                atomposition_j = atomposition_i[:] + hopposition_sp[:]
                #                    jatom =find_orbital(lattice_super,atomposition_j)
                jatom = find_orbital(lattice_super, supercell, atomposition_j, jband)
                #println("$iatom $jatom $atomposition_i $atomposition_j")
                add_hoppings!(lattice_super, v, iatom, jatom, hopposition_sp)
            end
        end



    elseif lattice.dim == 2

        for i2 = 1:supercell[2]
            for i1 = 1:supercell[1]
                for iatom = 1:numatoms
                    position_x = (lattice.positions[iatom][1] + (i1 - 1)) / supercell[1]
                    position_y = (lattice.positions[iatom][2] + (i2 - 1)) / supercell[2]

                    add_atoms!(lattice_super, [position_x, position_y])
                    #                    icount += 1
                    icount = ((i2 - 1) * supercell[1] + (i1 - 1)) * numatoms + iatom
                    diagonals[icount] = lattice.diagonals[iatom]
                end
            end
        end
        set_onsite!(lattice_super, diagonals)

        δ = 0
        for hopping in lattice.hoppings
            δ += 1
        #for δ = 1:lattice.numhopps
            #hopping = lattice.hoppings[δ] #Type:
            iband = hopping.ijatoms[1]
            jband = hopping.ijatoms[2]
            hopposition = deepcopy(hopping.ijpositions)
            #println(hopposition)
            hopposition_sp = hopposition[:]
            hopposition_sp[1] /= supercell[1]
            hopposition_sp[2] /= supercell[2]

            v = hopping.amplitude
            #println("$iband $jband $hopposition $hopposition_sp")
            for i2 = 1:supercell[2]
                for i1 = 1:supercell[1]
                    iiband = ((i2 - 1) * supercell[1] + (i1 - 1)) * numatoms + iband
                    atomposition_i = lattice_super.positions[iiband][:]
                    iatom = find_orbital(lattice_super, supercell, atomposition_i, iband)
                    #println("iband iatom $iband $iatom")
                    #                        iatom = find_orbital(lattice_super,atomposition_i)

                    atomposition_j = atomposition_i[:] + hopposition_sp[:]
                    #                        jatom =find_orbital(lattice_super,atomposition_j)
                    jatom = find_orbital(lattice_super, supercell, atomposition_j, jband)
                    #println("$iatom $jatom $atomposition_i $atomposition_j")
                    add_hoppings!(lattice_super, v, iatom, jatom, hopposition_sp)
                    positionindex = get_positionindex(lattice_super, atomposition_j)
                    #println("positionindex $positionindex")
                end
            end
        end



    elseif lattice.dim == 3
        for i3 = 1:supercell[3]
            for i2 = 1:supercell[2]
                for i1 = 1:supercell[1]
                    for iatom = 1:numatoms
                        position_x = (lattice.positions[iatom][1] + (i1 - 1)) / supercell[1]
                        position_y = (lattice.positions[iatom][2] + (i2 - 1)) / supercell[2]
                        position_z = (lattice.positions[iatom][3] + (i3 - 1)) / supercell[3]

                        add_atoms!(lattice_super, [position_x, position_y, position_z])
                        #                    icount += 1
                        icount =
                            (((i3 - 1) * supercell[2] + i2 - 1) * supercell[1] + (i1 - 1)) *
                            numatoms + iatom
                        diagonals[icount] = lattice.diagonals[iatom]
                    end
                end
            end
        end
        set_onsite!(lattice_super, diagonals)

        #for δ = 1:lattice.numhopps
        δ = 0
        for hopping in lattice.hoppings
            δ += 1
            #hopping = lattice.hoppings[δ] #Type:
            iband = hopping.ijatoms[1]
            jband = hopping.ijatoms[2]
            hopposition = deepcopy(hopping.ijpositions)
            #println(hopposition)
            hopposition_sp = hopposition[:]
            hopposition_sp[1] /= supercell[1]
            hopposition_sp[2] /= supercell[2]
            hopposition_sp[3] /= supercell[3]

            v = hopping.amplitude
            #println("$iband $jband $hopposition $hopposition_sp")
            for i3 = 1:supercell[3]
                for i2 = 1:supercell[2]
                    for i1 = 1:supercell[1]
                        iiband =
                            (((i3 - 1) * supercell[2] + i2 - 1) * supercell[1] + (i1 - 1)) *
                            numatoms + iband
                        atomposition_i = lattice_super.positions[iiband][:]
                        iatom =
                            find_orbital(lattice_super, supercell, atomposition_i, iband)
                        #                            iatom = find_orbital(lattice_super,atomposition_i)
                        atomposition_j = atomposition_i[:] + hopposition_sp[:]
                        jatom =
                            find_orbital(lattice_super, supercell, atomposition_j, jband)
                        #                            jatom =find_orbital(lattice_super,atomposition_j)
                        #println("$iatom $jatom $atomposition_i $atomposition_j")
                        add_hoppings!(lattice_super, v, iatom, jatom, hopposition_sp)
                    end
                end
            end
        end
    end

    return lattice_super
end

function get_positionindex(lattice::Lattice, position)


    positionindex = zeros(Int64, lattice.dim)
    for id = 1:lattice.dim
        positionindex[id] = round(Int, position[id], RoundDown)
        #=
        if position[id] < 0
            positionindex[id] = -round(Int,-position[id],RoundDown)-1
        else
            positionindex[id] = round(Int,position[id],RoundDown)
        end
        =#
    end

    return positionindex

    position_sp = position[:]
    positionindex = zeros(Int64, lattice.dim)
    for id = 1:lattice.dim
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


function cut_hg(lattice::Lattice)
    numatoms = lattice.numatoms
    la = set_Lattice(lattice.dim, lattice.vectors)
    hopindices = get_hopindex(lattice)

    for iatom = 1:lattice.numatoms
        add_atoms!(la, lattice.positions[iatom][:])
    end
    set_onsite!(la, lattice.diagonals)
    set_μ!(la, lattice.μ)

    if lattice.numhopps != 0
        δ = 0
        
        for hopping in lattice.hoppings
            δ += 1
        #for δ = 1:lattice.numhopps
            #hopping = lattice.hoppings[δ] #Type:
            iband = hopping.ijatoms[1]
            jband = hopping.ijatoms[2]
            hoppositions = hopping.ijpositions
            v = hopping.amplitude

            double = false
            if δ > 1
                δ2 = 0
                for hopping2 in lattice.hoppings
                    δ2  += 1
                #for δ2 = 1:δ-1
                    #hopping2 = lattice.hoppings[δ2] #Type:
                    iband2 = hopping2.ijatoms[1]
                    jband2 = hopping2.ijatoms[2]
                    hoppositions2 = -hopping2.ijpositions
                    if iband == jband2 &&
                       jband == iband2 &&
                       norm(hoppositions2 - hoppositions) == 0 && δ2 != δ
                        double = true
                        break
                    end
                end
            end
            if double == false
                add_hoppings!(la, v, iband, jband, hoppositions)
            end
            println(δ,"/",lattice.numhopps)
        end
    end

    return la

end


function hamiltonian_k(lattice)
    numatoms = lattice.numatoms
    diagonals = zeros(Float64, numatoms)
    vectors = lattice.vectors
    μ = lattice.μ
    for i = 1:length(lattice.diagonals)
        diagonals[i] = lattice.diagonals[i] - μ
    end

    function calc_ham(k)
        realpart_ham = zeros(Float64, numatoms, numatoms)
        imagpart_ham = zeros(Float64, numatoms, numatoms)
        for i = 1:numatoms
            realpart_ham[i, i] = diagonals[i]
        end
        if lattice.numhopps != 0
            δ = 0
            for hopping in lattice.hoppings
                δ += 1
            #for δ = 1:lattice.numhopps
                #hopping = lattice.hoppings[δ] #Type:
                i = hopping.ijatoms[1]
                j = hopping.ijatoms[2]
                hop = hopping.ijpositions
                vec_a = zeros(Float64, lattice.dim)
                for i = 1:length(lattice.vectors)
                    vec_a += hop[i] .* lattice.vectors[i]
                end
                ak = 0.0
                for idim = 1:lattice.dim
                    ak += k[idim] * vec_a[idim]
                end
                amp = hopping.amplitude*exp(im*ak)

                ampcos = real(amp)#hopping.amplitude * cos(ak)
                #println("$ampcos $i $j")
                realpart_ham[i, j] += ampcos
                realpart_ham[j, i] += ampcos
                ampsin = imag(amp)
                #ampsin = hopping.amplitude * sin(ak)
                imagpart_ham[i, j] += ampsin
                imagpart_ham[j, i] += -ampsin

            end
        end
        if sum(abs.(imagpart_ham)) == 0.0
            ham = realpart_ham[:, :]
            for i=1:numatoms
                for j=i:numatoms
                        ham[i,j] = realpart_ham[j,i]

                end
            end
        else
            ham = realpart_ham[:, :] + im * imagpart_ham[:, :]
            for i=1:numatoms
                for j=i:numatoms
                    ham[i,j] = conj(ham[j,i])
                    if i==j
                        ham[i,j] = real(ham[i,j])
                    end
                end
            end
        end
        #println(ham)

        ham

    end
    k -> calc_ham(k)
end

function show_neighbors(lattice; maxn = 1, cutoff = 2.0)
    dim = lattice.dim
    println("Possible hoppings")
    if dim == 1
        for i1 = -maxn:maxn
            for i = 1:lattice.numatoms
                r0 = [i1]
                for i = 1:lattice.numatoms
                    ri = lattice.positions[i]
                    for j = i:lattice.numatoms
                        rj = lattice.positions[j] + r0
                        r = rationalize.(rj - ri)
                        r_pos = get_position(lattice, r)
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
                r0 = [i1, i2]
                for i = 1:lattice.numatoms
                    ri = lattice.positions[i]
                    for j = i:lattice.numatoms
                        rj = lattice.positions[j] + r0
                        r = rationalize.(rj - ri)
                        r_pos = get_position(lattice, r)
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
                    r0 = [i1, i2, i3]
                    for i = 1:lattice.numatoms
                        ri = lattice.positions[i]
                        for j = i:lattice.numatoms
                            rj = lattice.positions[j] + r0
                            r = rationalize.(rj - ri)
                            r_pos = get_position(lattice, r)
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
