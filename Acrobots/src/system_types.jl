macro make_type(typename, parent, fields...)
    ex = :(immutable $(esc(typename)){T} <: $(parent){$(length(fields)), T} end)
    push!(ex.args, Expr(:block))
    for field in fields
        push!(ex.args[3].args, :($(field)::T))
    end
    constructor_expr = quote function $(esc(typename))(a::NTuple{$(length(fields)), T})
            new{T}()
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
