function write_energydispersion(lattice,fileheader;nka=20,nkb=20,nkc=20)
    ham = hamiltonian_k(lattice)
    kas = range(0,1,length=nka)
    kbs = range(0,1,length=nkb)
    kcs = range(0,1,length=nkc)
    b1 = lattice.rvectors[1]
    b2 = lattice.rvectors[2]
    b3 = lattice.rvectors[3]
    fps = []
    for i=1:lattice.numatoms
        fp = open(fileheader*lpad("$i",3,"0")*".txt","w")
        push!(fps,fp)
    end

    for ka in kas
        for kb in kbs
            for kc in kcs
                println("$ka $kb $kc")
                k = b1*ka + b2*kb + b3*kc
                energy = eigen(ham(k)).values
                for i=1:lattice.numatoms
                    println(fps[i],"$(k[1]) $(k[2]) $(k[3]) $(energy[i])")
                end
            end
        end
    end
    for i=1:lattice.numatoms
        close(fps[i])
    end

end


function write_energydispersion_abc(lattice,fileheader;nka=20,nkb=20,nkc=20)
    ham = hamiltonian_k(lattice)
    kas = range(0,1,length=nka)
    kbs = range(0,1,length=nkb)
    kcs = range(0,1,length=nkc)
    b1 = lattice.rvectors[1]
    b2 = lattice.rvectors[2]
    b3 = lattice.rvectors[3]
    fps = []
    for i=1:lattice.numatoms
        fp = open(fileheader*lpad("$i",3,"0")*".txt","w")
        push!(fps,fp)
    end

    for ka in kas
        for kb in kbs
            for kc in kcs
                println("$ka $kb $kc")
                k = b1*ka + b2*kb + b3*kc
                energy = eigen(ham(k)).values
                for i=1:lattice.numatoms
                    println(fps[i],"$(ka) $(kb) $(kc) $(energy[i])")
                end
            end
        end
    end
    for i=1:lattice.numatoms
        close(fps[i])
    end

end

function dispersion(dim, n, ham, k)
    if n == 1
        return real([ham(k)[1]])
    else
        #println("dis")
        #println(ham(k))
        #println("dd")
        energy = eigen(ham(k)).values
        return energy
    end
end

function calc_band(kmin_real, kmax_real, nk, lattice, ham)
    #println("cd")
    n = lattice.numatoms
    dim = lattice.dim
    energies = zeros(Float64, n, nk)
    vec_k = zeros(Float64, dim, nk)
    for idim = 1:dim
        k = range(kmin_real[idim], length = nk, stop = kmax_real[idim])
        vec_k[idim, :] = k[:]
    end
    for ik = ProgressBar(1:nk)
        #println("ik ",ik)
        #println("ik = $ik ",vec_k[:, ik])
        ene =  dispersion(dim, n, ham, vec_k[:, ik])
        #println(ene)
        energies[:, ik] = ene#dispersion(dim, n, ham, vec_k[:, ik])
    end

    return vec_k, energies
end

function get_DOS(lattice, nk; nbins = 100)
    dim = lattice.dim
    ham = hamiltonian_k(lattice)
    n = lattice.numatoms


    if dim == 1
        energies = zeros(Float64, n, nk)
        dk = 1 / (nk - 1)
        for ik = 1:nk
            k = (ik - 1) * dk
        end
    elseif dim == 2
        energies = zeros(Float64, n, nk, nk)
        dk1 = 1 / (nk - 1)
        dk2 = 1 / (nk - 1)
        for ik1 = 1:nk
            k1 = (ik1 - 1) * dk1
            for ik2 = 1:nk
                k2 = (ik2 - 1) * dk2
                kx, ky = get_position_kspace(lattice, [k1, k2])
                energy = dispersion(dim, n, ham, [kx, ky])
                energies[:, ik1, ik2] = energy[:]

            end
        end
    elseif dim == 3
        energies = zeros(Float64, nk, nk, nk)
        dk1 = 1 / (nk - 1)
        dk2 = 1 / (nk - 1)
        dk3 = 1 / (nk - 1)

        for ik1 = 1:nk
            k1 = (ik1 - 1) * dk1
            for ik2 = 1:nk
                k2 = (ik2 - 1) * dk2
                for ik3 = 1:nk
                    k3 = (ik3 - 1) * dk3
                    kx, ky, kz = get_position_kspace(lattice, [k1, k2, k3])
                    energy = dispersion(dim, n, ham, [kx, ky, kz])
                    energies[:, ik1, ik2, ik3] = energy[:]
                end

            end
        end
    end
    energies = vec(energies)
    #            println(energies)
    #pls = histogram(energies,bins=nbins,normalize=true)
    hist = fit(Histogram, energies, nbins = nbins)

    return hist
end