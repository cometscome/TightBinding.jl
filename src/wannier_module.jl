module Wannier
    using ProgressBars

    function make_hamiltonian(dr,num,mat_hop,kx,ky,kz)
        mat_ham = zeros(ComplexF64,num,num)
        for jband=1:num
            for iband=1:num
                for iz=-dr:dr
                    for iy=-dr:dr
                        for ix=-dr:dr
                            mat_ham[iband,jband] +=  mat_hop[ix+dr+1,iy+dr+1,iz+dr+1,iband,jband]*exp(im*ix*kx)*exp(im*iy*ky)*exp(im*iz*kz)
                        end
                    end
                end
            end
        end

        for i=1:num
            for j=1:num
                if i> j
                    mat_ham[i,j] = conj(mat_ham[j,i])
                end
                if i==j
                    mat_ham[i,j] = real(mat_ham[i,j])
                end
            end
        end
        return mat_ham
    end


    function wannier_load(dr,num,filename;eps = 0.0)
        mat_hop = zeros(ComplexF64,2dr+1,2dr+1,2dr+1,num,num)

        data = readlines(pwd() * "/" * filename)
        i = 2
        u = data[i]
        #println(u)
        num_wann = parse(Int64, u)
        i = 3
        nrpts = parse(Int64, data[i])
        println("num_wann = ",num_wann)
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
                    ix,iy,iz = hops[1],hops[2],hops[3]
                    iband = parse(Int64, u[4])
                    jband = parse(Int64, u[5])
                    @assert ν == iband
                    @assert μ == jband
                    tr = parse(Float64, u[6])/factors[k]
                    ti =  parse(Float64, u[7])/factors[k]
                    if abs(ix) > dr || abs(iy) > dr || abs(iz) > dr
                    else
                        mat_hop[ix+dr+1,iy+dr+1,iz+dr+1,iband,jband] = (tr+im*ti)
                        if sqrt(tr^2+ti^2) < eps 
                            mat_hop[ix+dr+1,iy+dr+1,iz+dr+1,iband,jband] = 0
                        end
                    end

                end
            end
            #println(k,"/",nrpts)
            #println(la_hr.numhopps)
        end
    
        return mat_hop


    end
end