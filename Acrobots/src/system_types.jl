macro make_type(typename, parent, fields...)
    ex = :(immutable $(esc(typename)){$(esc(:T))} <: $(parent){$(length(fields)), $(esc(:T))} end)
    push!(ex.args, Expr(:block))
    for field in fields
        push!(ex.args[3].args, :($(field)::$(esc(:T))))
    end
    constructor_expr = quote function $(esc(typename))(a::NTuple{$(length(fields)), $(esc(:T))})
            $(esc(:new)){$(esc(:T))}()
            end
        end
    constructor_expr = constructor_expr.args[2] # remove the `begin` wrapper
    for i = 1:length(fields)
        push!(constructor_expr.args[2].args[2].args, :(a[$(i)]))
    end
    push!(ex.args[3].args, constructor_expr)
    return ex
end

abstract Position{N, T} <: FixedVectorNoTuple{N, T}
abstract Velocity{N, T} <: FixedVectorNoTuple{N, T}
abstract Input{N, T} <: FixedVectorNoTuple{N, T}
abstract Output{N, T} <: FixedVectorNoTuple{N, T}
abstract State{N, T} <: FixedVectorNoTuple{N, T}

abstract DynamicalSystem{StateType, InputType, OutputType}

function isapprox{FSA <: FixedSizeArrays.FixedArray, A <: Union{Array, FixedSizeArrays.FixedArray}}(a::FSA, b::A, atol::Real)
    for i=1:length(a)
        !isapprox(a[i], b[i]; atol=atol) && return false
    end
    true
end
