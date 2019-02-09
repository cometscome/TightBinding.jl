module TBfromk2r
    using Dierckx
    using SparseArrays
    const eps = 1e-7

    const σ0 = [
    1 0
    0 1
    ]

    const σx = [
    0 1
    1 0
    ]

    const σy = [
    0 -im
    im 0
    ]

    const σz = [
    1 0
    0 -1
    ]

    export surfaceHamiltonian,σx,σy,σz,σ0

    function surfaceHamiltonian(Hk,norb;numhop=3,L=32,kpara="kx",BC="OBC")
        function calc_H(k)            
            if kpara == "kx"
                fx(kperp) = Hk([k,kperp])[:,:]
                mat_H = make_Hamiltonian(fx,[numhop],L,norb,BC)
            elseif kpara == "ky"
                fy(kperp) = Hk([kperp,k])[:,:]
                mat_H = make_Hamiltonian(fy,[numhop],L,norb,BC)
            else
                println("Error! not supported.")
                exit()
            end

            
            return mat_H
        end
        k ->  calc_H(k)
    end

    function make_Hamiltonian(Hk,numhop,L,norb,BC)
        Vol = prod(L)
        dim = length(L)
        mat_H = spzeros(ComplexF64,Vol*norb,Vol*norb)
        #println(Vol*norb)
        

        mat_t = make_hop(Hk,numhop,norb)  
        
        for i = 1:Vol
            for jorb = 1:norb
                for iorb = 1:norb
                    ii = (iorb-1)*Vol + i
                    #println(ii,"\t",mat_t[iorb,jorb])
                    isite = calc_i2site(i,L)

                    for hop in mat_t[iorb,jorb]
                        jsite = isite .+ hop[1]
                        if BC == "PBC"
                            for il=1:dim
                                jsite[il] += ifelse(jsite[il] < 1,L[il],0) 
                                jsite[il] += ifelse(jsite[il] > L[il],-L[il],0) 
                            end
                        end
                        flag = false

                        if dim ==1                            
                            if 1 <= jsite[1] <= L[1]
                                flag = true
                            end
                        elseif dim == 2
                            if 1 <= jsite[1] <= L[1] && 1 <= jsite[2] <= L[2]
                                flag = true
                            end
                        end
                        
                        if flag
                            j = calc_r2j(jsite,L)                        
                            
                            jj = (jorb-1)*Vol + j
                        
                            mat_H[ii,jj] = hop[2]
                        end
                    end
                end
            end
        end
        return mat_H

    end

    function calc_i2site(i,L)
        dim = length(L)
        if dim == 1
            isite = [i]
        elseif dim ==2
            ix = (i-1) % L[1] + 1
            iy = div(i-ix,L[1])+ 1
            isite = [ix,iy]
        end 
        return  isite          
    end

    function calc_r2j(hopping,L)
        dim = length(L)
        if dim == 1
            j = hopping[1]
        elseif dim == 2
            j = (hopping[2]-1)*L[1]+hopping[1]
        end
        return j
    end

    function make_hop(Hk,n::Array{Int,1},norb)                
        mat_t  = Array{Any}(undef, norb, norb)

        for i=1:norb
            for j=1:norb
                f(k) = Hk(k)[i,j]
                hop_t = calc_coeff(f,n)
                mat_t[i,j] = hop_t
                #println("hoppings matrix at ($i,$j) element")
                #for hop in hop_t
                #    println(hop[1],"\t",hop[2])
                #end
            end
        end

        #println(mat_t)
        return mat_t

    end

    function calc_coeff(fk,n::Array{Int,1})
        dim = length(n)
        
        #println(dim)
        if dim == 1
            ck = calc_coeff(fk,n[1],1)
        elseif dim == 2
            ck = calc_coeff(fk,n[1],n[2],1,1)
        else
            println("not supported yet: dimension is $dim")
            return
        end
        return ck
    end

    function calc_coeff(fk,nkx::Int,ax::Number;L=256)
        numhop = 2nkx + 1        
        ck = zeros(ComplexF64,numhop)
        
        hopping = range(-nkx, stop=nkx)
        kmax_x = 2π/ax
        kx = range(0.0,stop=kmax_x ,length=L)

        for i=1:numhop
            f(k) = fk(k)*exp(im*hopping[i]*k)
            mat_f = make_f(f,L,kmax_x)./(2π)
            spline_re = Spline1D(kx, real.(mat_f))
            spline_im = Spline1D(kx, imag.(mat_f))
            integ_re = integrate(spline_re, 0, kmax_x)
            integ_im = integrate(spline_im, 0, kmax_x)

            if abs(integ_re) > eps
                ck[i] += integ_re
            end
            if abs(integ_im) > eps
                ck[i] += im*integ_im
            end

        end

        hop_t = []
        for ix=1:numhop
        
            if abs(ck[ix]) > eps
                push!(hop_t,
                    (((hopping[ix]),ck[ix]))
                )
            end
        
        end


        return hop_t #ck ./(2π)
    end

    function make_f(f,L,kmax_x,kmax_y)
        mat_f = zeros(ComplexF64,L,L)
        kx = range(0.0,stop=kmax_x,length=L)
        ky = range(0.0,stop=kmax_y,length=L)

        for ikx =1:L
            for iky=1:L
                mat_f[ikx,iky] = f([kx[ikx],ky[iky]])#f(kx[ikx],ky[iky])
            end
        end
            
        return mat_f
    end

    function make_f(f,L,kmax_x)
        mat_f = zeros(ComplexF64,L)
        kx = range(0.0,stop=kmax_x,length=L)

        for ikx =1:L
            mat_f[ikx] = f(kx[ikx])#f(kx[ikx],ky[iky])
        end
            
        return mat_f
    end

    function calc_coeff(fk,nkx::Int,nky::Int,ax::Number,ay::Number;L=256)
        numhop_x = 2nkx + 1        
        numhop_y = 2nky + 1  
        ck = zeros(ComplexF64,numhop_x,numhop_y)
        hopping_x = range(-nkx, stop=nkx)
        hopping_y = range(-nky, stop=nky)
        kmax_x = 2π/ax
        kmax_y = 2π/ay
        kx = range(0.0,stop=kmax_x ,length=L)
        ky = range(0.0,stop=kmax_y,length=L)
        hop_t = []
        
        
        for ix=1:numhop_x
            for iy=1:numhop_y
                f(k) = fk(k)*exp(im*hopping_x[ix]*k[1])*exp(im*hopping_y[iy]*k[2])
                mat_f = make_f(f,L,kmax_x,kmax_y)./(2π)^2
                spline_re = Spline2D(kx, ky, real.(mat_f))
                spline_im = Spline2D(kx, ky, imag.(mat_f))
                integ_re = integrate(spline_re, 0, kmax_x, 0, kmax_y)
                integ_im = integrate(spline_im, 0, kmax_x, 0, kmax_y)
                if abs(integ_re) > eps
                    ck[ix,iy] += integ_re
                end
                if abs(integ_im) > eps
                    ck[ix,iy] += im*integ_im
                end
                #ck[ix,iy] = integ_re + im*integ_im                
                #println("$(hopping_x[ix]), $(hopping_y[iy]), $(ck[ix,iy])")
            end
        end

        hop_t = []
        for ix=1:numhop_x
            for iy=1:numhop_y
                if abs(ck[ix,iy]) > eps
                    push!(hop_t,
                        (((hopping_x[ix],hopping_y[iy]),ck[ix,iy]))
                    )
                end
            end
        end
        
        return hop_t

    end


end