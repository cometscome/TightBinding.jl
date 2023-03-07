function Realspace_Hamiltonian(Hk, norb; numhop = 3, L = [8,8], BC = ["OBC","OBC"])
    function calc_H(k)
        if kpara == "kx"
            fx(kperp) = Hk([k, kperp])[:, :]
            mat_H = make_Hamiltonian(fx, [numhop], L, norb, BC)
        elseif kpara == "ky"
            fy(kperp) = Hk([kperp, k])[:, :]
            mat_H = make_Hamiltonian(fy, [numhop], L, norb, BC)
        else
            println("Error! not supported.")
            exit()
        end


        return mat_H
    end
    k -> calc_H(k)
end

function calc_coeff(fk, numhop_x,numhop_y,ax,ay,L = [16,16];eps = 1e-7)
    hopping_x = range(-numhop_x,stop = numhop_x)
    hopping_y = range(-numhop_y,stop = numhop_y)
    kmax_x = 2π / ax
    kmax_y = 2π / ay
    kxs = range(0.0, stop = kmax_x, length = L[1])
    kys = range(0.0, stop = kmax_x, length = L[2])
    ck = zeros(ComplexF64,2hopping_x+1,2hopping_y+1)
    mat_f = zeros(ComplexF64,L[1],L[2])
    for dx = -numhop_x:numhop_x
        for dy = -numhop_y:numhop_y
            f(k) = fk(k) * exp(im * hopping[dx+numhop_x+1] * k[1])* exp(im * hopping[dy+numhop_y+1] * k[2])
            for ikx=1:L[1]
                kx = kxs[ikx]
                for iky=1:L[2]
                    ky = kys[iky]
                    mat_f[ikx,iky] = f([kx,ky])./ (2π)
                end
            end
            spl_re = Spline2D(kxs, kys, real.(mat_f))
            spl_im = Spline2D(kxs, kys, imag.(mat_f))
            integ_re = integrate(spl_re, 0, kmax_x,0,kmax_y)
            integ_im = integrate(spl_im, 0, kmax_x,0,kmax_y)
            if abs(integ_re) > eps
                ck[dx+numhop_x+1,dy+numhop_y+1] += integ_re
            end
            if abs(integ_im) > eps
                ck[dx+numhop_x+1,dy+numhop_y+1] += im * integ_im
            end
        end
    end

    hop_t = []
    for dx = -numhop_x:numhop_x
        for dy = -numhop_y:numhop_y
            if abs(ck[dx+numhop_x+1,dy+numhop_y+1]) > eps
                push!(hop_t, (((hopping[dx+numhop_x+1],hopping[dy+numhop_y+1]), ck[dx+numhop_x+1,dy+numhop_y+1])))
            end
        end
    end
    return hop_t
end