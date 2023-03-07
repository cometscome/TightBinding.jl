struct Kpoints
    kmin::Array{<:Number,1}
    kmax::Array{<:Number,1}
    nk::Int64
    name_start::Any
    name_end::Any
end

mutable struct Klines
    numlines::Int64
    kpoints::Array{Kpoints,1}
end


function set_Klines()
    kpoint = Array{Kpoints,1}[]
    klines = Klines(0, kpoint)
    return klines
end

function add_Kpoints!(klines, kmin, kmax, name_start, name_end; nk = 20)
    kpoint = Kpoints(kmin, kmax, nk, name_start, name_end)
    klines.numlines += 1
    push!(klines.kpoints, kpoint)
end

