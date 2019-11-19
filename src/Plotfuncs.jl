#module Plotfuncs
    using .Plots
    using LinearAlgebra
    using StatsBase
    #using ..TightBinding
    export plot_lattice_2d,calc_band_plot,plot_DOS,calc_band_plot_finite
    export plot_fermisurface_2D
    export get_DOS


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

    """
        plot_fermisurface_2D(lattice::Lattice;Eshift = 0.0,nk = 20)

    Show the Fermi surface for 2D system. \\
    This plots the contour with E=0.0+Eshift \\
    Eshift: the energy shift from the chemical potential lattice.μ \\
    nk: the number of meshes in 2D momentum space. The total mesh is nk x nk.
    """
    function plot_fermisurface_2D(lattice;Eshift = 0.0,nk = 50)
        ham = hamiltonian_k(lattice) #Construct the Hamiltonian
        dim = lattice.dim
        n = lattice.numatoms

        if dim != 2
            println("This function only supports 2D case. lattice.dim should be 2. But dim = $dim .")
            return
        end

        k1s = range(-0.51, stop = 0.51, length = nk)
        k2s = range(-0.51, stop = 0.51, length = nk)
        kx_vec = zeros(Float64,nk)
        ky_vec = zeros(Float64,nk)
        energies = zeros(Float64,n,nk,nk)
        colors = ["red","blue","orange","brown","yellow"]
        kxmin,kymin = get_position_kspace(lattice,[-0.5,-0.5])
        kxmax,kymax = get_position_kspace(lattice,[0.5,0.5])

        for i1=1:nk
            for i2=1:nk
                kx,ky = get_position_kspace(lattice,[k1s[i1],k2s[i2]])
                kx_vec[i1] = kx
                ky_vec[i2] = ky
                energy = dispersion(dim,n,ham,[kx,ky])
                for j=1:n
                    energies[j,i1,i2] = energy[j]
                end

            end
        end
#        println("energies ",energies[1,:,:])
        #println(kx_vec)
        i = 1
        pls = contour(kx_vec,ky_vec,energies[i,:,:],levels=[Eshift],
            xlabel="kx",ylabel="ky",aspect_ratio=1,
            xlims=(kxmin,kxmax),ylims=(kymin,kymax)
            )
        for i=2:n
            pls = contour!(kx_vec,ky_vec,energies[i,:,:],levels=[Eshift],
                xlabel="kx",ylabel="ky",aspect_ratio=1,
                xlims=(kxmin,kxmax),ylims=(kymin,kymax)
                )
        end
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
            energies_i = sort(energies_i,dims=1)
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

    function calc_band_plot_finite(klines,lattice,direction;periodic = true,nsites = 20)
        ham = hamiltonian_k_1d(lattice,direction,periodic=periodic,nsites=nsites)
        numlines = klines.numlines
        dim = lattice.dim-1
        klength = 0.0
        vec_k = Float64[]
        numatoms = lattice.numatoms
        N = numatoms*nsites
        energies = zeros(Float64,0,N)
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

            vec_k_temp = zeros(Float64,dim,nk)

            for idim=1:dim
                k = range(kmin_real[idim],length=nk,stop=kmax_real[idim])
                vec_k_temp[idim,:] = k[:]
            end
            energies_i = zeros(Float64,nk,N)
            for ik=1:nk
                energies_i[ik,:] = eigen(ham(vec_k_temp[:,ik])).values[:]
            end

            vec_k_i = range(klength,length=nk,stop=klength+kdistance)
            vec_k = vcat(vec_k,vec_k_i)
            energies = vcat(energies,energies_i)
            klength += kdistance
            push!(xticks_values,klength)
            push!(xticks_labels,klines.kpoints[i].name_end)
        end

        pls = plot(vec_k[:],energies,
                        xticks = (xticks_values,xticks_labels),legend=false)

        return pls
    end


    function plot_DOS(lattice,nk;nbins=100)
        dim = lattice.dim
        hist = get_DOS(lattice,nk;nbins=nbins)
        pls = bar(hist.edges[1][1:end],hist.weights/nk^dim)
        return pls

        #=
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
        pls = histogram(energies,bins=nbins,normalize=true)

        return pls
        =#
    end

#end
