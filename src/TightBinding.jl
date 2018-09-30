module TightBinding
    export set_Lattice,add_atoms!,add_hoppings!

    mutable struct Lattice
        dim::Int8
        vectors::Array{Array{Float64,1},1}
        numatoms::Int8
        positions::Array{Array{Float64,1},1}
        hoppings::Array{Array{Int64,1},1}
    end

    function set_Lattice(dim,vectors)
        positions = Array{Float64,1}[]
        hoppings = Array{Int64,1}[]
        numatoms = 0
        lattice = Lattice(dim,vectors,numatoms,positions,hoppings)
        return lattice
    end

    function add_atoms!(lattice,position::Array{<:Number,1})
        push!(lattice.positions,position)
    end

    function add_atoms!(lattice,positions::Array{<:Array})
        for i=1:length(positions)
            push!(lattice.positions,positions[i])
        end
    end

    function add_hoppings!(lattice,hopping::Array{<:Number,1})
        push!(lattice.hoppings,hopping)
    end

    function add_hoppings!(lattice,hoppings::Array{<:Array})
        for i=1:length(hoppings)
            push!(lattice.hoppings,hoppings[i])
        end
    end

greet() = print("Hello World!")

end # module
