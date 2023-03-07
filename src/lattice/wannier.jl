
function read_wannier(lattice::Lattice, filename)
    la_hr = set_Lattice(lattice.dim, lattice.vectors)

    for iatom = 1:lattice.numatoms
        add_atoms!(la_hr, lattice.positions[iatom][:])
    end

    set_μ!(la_hr, lattice.μ)
    diagonals = zeros(Float64, lattice.numatoms)


    data = readlines(pwd() * "/" * filename)
    i = 2
    u = data[i]
    #println(u)
    num_wann = parse(Int64, u)
    @assert num_wann == lattice.numatoms
    i = 3
    nrpts = parse(Int64, data[i])
    factors = zeros(Int64, nrpts)

    num = div(nrpts, 15)
    nmod = nrpts % 15
    count = 0
    for j = 1:num
        i += 1
        u = split(data[i])
        for k = 1:15
            count += 1
            factors[count] = parse(Int64, u[k])
        end
    end
    if nmod != 0
        i += 1
        u = split(data[i])
        for k = 1:nmod
            count += 1
            factors[count] = parse(Int64, u[k])
        end
    end
    

    println("reading a wannier format")
    #println(length(factors),"\t", nrpts)
    for k in ProgressBar(1:nrpts) #1:nrpts
        for μ = 1:num_wann
            for ν = 1:num_wann
                i += 1
                u = split(data[i])
                hops = parse.(Int64, u[1:3])
                iband = parse(Int64, u[4])
                jband = parse(Int64, u[5])
                @assert ν == iband
                @assert μ == jband
                if sum(abs.(hops)) == 0 && iband == jband
                    diagonals[iband] = parse(Float64, u[6])
                else
                    tij = parse(Float64, u[6]) - im * parse(Float64, u[7])
                    #tij = parse(Float64, u[6]) + im * parse(Float64, u[7])


                    if abs(tij) != 0
                        #atomposition_i = lattice.positions[iband][:]
                        #atomposition_j = lattice.positions[jband][:] 

                        Ri = lattice.positions[iband][:]
                        Rj = lattice.positions[jband][:]
                        for id = 1:lattice.dim
                            Rj[id] += hops[id]
                        end
                        dr = Rj - Ri


                        #for id=1:lattice.dim
                        #    atomposition_j[id] += hops[id]#*lattice.vectors[id][:]
                        #end
                        #hopping = atomposition_j- atomposition_i

                        add_hoppings!(la_hr, tij/factors[k], iband, jband, dr)
                        #add_hoppings!(la_hr, tij, iband, jband, dr)
                    end
                end

            end
        end
        #println(k,"/",nrpts)
        #println(la_hr.numhopps)
    end


    set_onsite!(la_hr, diagonals .+ lattice.μ)

    return la_hr#cut_hg(la_hr)


end

function write_hr(lattice::Lattice; filename = "wannier90_hr.dat", maxhop = 10)
    numatoms = lattice.numatoms
    hopindices = get_hopindex(lattice)


    diagonals = zeros(Float64, numatoms)
    vectors = lattice.vectors
    μ = lattice.μ
    for i = 1:length(lattice.diagonals)
        diagonals[i] = lattice.diagonals[i] - μ
    end
    hopdata = Array{Int64,1}[]
    valdata = Array{Float64,1}[]
    for i = 1:numatoms
        v = diagonals[i]
        iband = i
        jband = i
        hop = [0, 0, 0, iband, jband]
        val = [real(v), -imag(v)]
        push!(hopdata, hop)
        push!(valdata, val)
    end
    if lattice.numhopps != 0
        δ = 0
        for hopping in lattice.hoppings
            δ += 1
        #for δ = 1:lattice.numhopps
            #hopping = lattice.hoppings[δ] #Type:
            #positionindex = hopindices[δ][:]
            iband = hopping.ijatoms[1]
            jband = hopping.ijatoms[2]


            ri = lattice.positions[iband][:]
            rhop = ri + hopping.ijpositions
            positionindex = get_positionindex(lattice, rhop)

            #=
            hopping = lattice.hoppings[δ] #Type:
            iband=hopping.ijatoms[1]
            jband=hopping.ijatoms[2]
            hopposition = hopping.ijpositions[:]

            atomposition_i = lattice.positions[iband][:]
            atomposition_j = atomposition_i[:] + hopposition[:]

            positionindex = get_positionindex(lattice,atomposition_j)
            =#
            #=
            Ri = lattice.positions[iband][:]
            Rj = lattice.positions[jband][:]
            for id=1:lattice.dim
                Rj[id] += positionindex[id]
            end
            dr = Rj-Ri
            if norm(dr - hopping.ijpositions) != 0
                println("positionindex: $positionindex")
                println("Ri $Ri Rj $Rj")
                println("dr $dr ijpositions",hopping.ijpositions)
            end
            =#

            v = hopping.amplitude
            if lattice.dim == 1
                hop = [positionindex[1], 0, 0, iband, jband]
            elseif lattice.dim == 2
                hop = [positionindex[1], positionindex[2], 0, iband, jband]
            elseif lattice.dim == 3
                hop = [positionindex[1], positionindex[2], positionindex[3], iband, jband]
            end
            val = [real(v), imag(v)]
            push!(hopdata, hop)
            push!(valdata, val)

            ri = lattice.positions[jband][:]
            rhop = ri - hopping.ijpositions


            positionindex = get_positionindex(lattice, rhop)

            #=

            Ri = lattice.positions[jband][:]
            Rj = lattice.positions[iband][:]
            for id=1:lattice.dim
                Rj[id] += positionindex[id]
            end
            dr = Rj-Ri
            if norm(dr + hopping.ijpositions) != 0
                println("rhop $rhop")
                println("positionindex: $positionindex")
                println("Ri $Ri Rj $Rj")
                println("oriRj ", lattice.positions[iband][:])
                println("dr $dr ijpositions",-hopping.ijpositions)
            end
            =#


            if lattice.dim == 1
                hop = [positionindex[1], 0, 0, jband, iband]
            elseif lattice.dim == 2
                hop = [positionindex[1], positionindex[2], 0, jband, iband]
            elseif lattice.dim == 3
                hop = [positionindex[1], positionindex[2], positionindex[3], jband, iband]
            end
            val = [real(v), -imag(v)]
            push!(hopdata, hop)
            push!(valdata, val)

        end
    end
    count = length(hopdata)
    if lattice.dim == 1
        hopmatrix = zeros(ComplexF64, 2 * maxhop + 1, 1, 1, numatoms, numatoms)
    elseif lattice.dim == 2
        hopmatrix = zeros(ComplexF64, 2 * maxhop + 1, 2 * maxhop + 1, 1, numatoms, numatoms)
    elseif lattice.dim == 3
        hopmatrix = zeros(
            ComplexF64,
            2 * maxhop + 1,
            2 * maxhop + 1,
            2 * maxhop + 1,
            numatoms,
            numatoms,
        )
    end

    for i = 1:count
        hop = hopdata[i]
        val = valdata[i]
        ix = hop[1] + maxhop + 1
        iy = hop[2] + maxhop + 1
        iz = hop[3] + maxhop + 1
        μ = hop[4]
        ν = hop[5]
        if lattice.dim == 1
            hopmatrix[ix, 1, 1, μ, ν] = val[1] + im * val[2]
        elseif lattice.dim == 2
            hopmatrix[ix, iy, 1, μ, ν] = val[1] + im * val[2]
        elseif lattice.dim == 3
            hopmatrix[ix, iy, iz, μ, ν] = val[1] + im * val[2]
        end
    end

    hopcount = 0
    if lattice.dim == 1
        for ix = 1:2*maxhop+1
            if sum(abs.(hopmatrix[ix, 1, 1, :, :])) != 0
                hopcount += 1
            end
        end
    elseif lattice.dim == 2
        for ix = 1:2*maxhop+1
            for iy = 1:2*maxhop+1
                if sum(abs.(hopmatrix[ix, iy, 1, :, :])) != 0
                    hopcount += 1
                end
            end
        end
    elseif lattice.dim == 3
        for ix = 1:2*maxhop+1
            for iy = 1:2*maxhop+1
                for iz = 1:2*maxhop+1
                    if sum(abs.(hopmatrix[ix, iy, iz, :, :])) != 0
                        hopcount += 1
                    end
                end
            end
        end
    end

    count = hopcount #*numatoms*numatoms

    fp = open(filename, "w")
    println(fp, "generated by TightBinding.jl")
    println(fp, numatoms)
    println(fp, count)
    num = div(count, 15)
    nmod = count % 15
    for j = 1:num
        for i = 1:15
            print(fp, 1, "\t")
        end
        println(fp, "\t")
    end
    if nmod != 0
        for i = 1:nmod
            print(fp, 1, "\t")
        end
        println(fp, "\t")
    end

    if lattice.dim == 1
        for i1 = 1:2*maxhop+1
            if sum(abs.(hopmatrix[i1, 1, 1, :, :])) != 0
                ix = i1 - maxhop - 1
                iy = 0
                iz = 0
                for ν = 1:numatoms
                    for μ = 1:numatoms
                        val = [
                            real(hopmatrix[i1, 1, 1, μ, ν]),
                            -imag(hopmatrix[i1, 1, 1, μ, ν]),
                        ]
                        println(fp, "$ix\t$iy\t$iz\t$μ\t$ν\t$(val[1])\t$(val[2])")
                    end
                end
            end
        end
    elseif lattice.dim == 2
        for i1 = 1:2*maxhop+1
            for i2 = 1:2*maxhop+1
                if sum(abs.(hopmatrix[i1, i2, 1, :, :])) != 0
                    ix = i1 - maxhop - 1
                    iy = i2 - maxhop - 1
                    iz = 0

                    for ν = 1:numatoms
                        for μ = 1:numatoms
                            val = [
                                real(hopmatrix[i1, i2, 1, μ, ν]),
                                imag(hopmatrix[i1, i2, 1, μ, ν]),
                            ]
                            println(fp, "$ix\t$iy\t$iz\t$μ\t$ν\t$(val[1])\t$(val[2])")
                        end
                    end
                end
            end
        end
    elseif lattice.dim == 3
        for i1 = 1:2*maxhop+1
            for i2 = 1:2*maxhop+1
                for i3 = 1:2*maxhop+1
                    if sum(abs.(hopmatrix[i1, i2, i3, :, :])) != 0
                        ix = i1 - maxhop - 1
                        iy = i2 - maxhop - 1
                        iz = i3 - maxhop - 1


                        for ν = 1:numatoms
                            for μ = 1:numatoms
                                val = [
                                    real(hopmatrix[i1, i2, i3, μ, ν]),
                                    imag(hopmatrix[i1, i2, i3, μ, ν]),
                                ]
                                println(fp, "$ix\t$iy\t$iz\t$μ\t$ν\t$(val[1])\t$(val[2])")
                            end
                        end
                    end
                end
            end
        end
    end
    close(fp)
    return


    for i = 1:count
        hop = hopdata[i]
        val = valdata[i]
        ix = hop[1]
        iy = hop[2]
        iz = hop[3]
        μ = hop[4]
        ν = hop[5]

        println(fp, "$ix\t$iy\t$iz\t$μ\t$ν\t$(val[1])\t$(val[2])")
    end
    close(fp)

    #println(hopdata)
    #println(valdata)
end
