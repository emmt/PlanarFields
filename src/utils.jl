
"""
    relative_precision(x::Number) -> r
    relative_precision(T::Type{<:Number}) -> r

yield the relative precision of the number `x`. Argument may also be a type. The result is
a real number which only depend on the bare numeric type of `x` or `T`. If the bare
numeric type is a floating-point type `F`, the result is `eps(F)`. If the bare numeric
type is an integer or rational type `R`, the result is `zero(R)`. Otherwise,
`eps(Float64)` is returned.

"""
relative_precision(x::Number) = relative_precision(typeof(x))
relative_precision(::Type{T}) where {T<:Number} = relative_precision(real_type(T))
relative_precision(::Type{T}) where {T<:Integer} = zero(T)
relative_precision(::Type{T}) where {T<:Rational} = zero(T)
relative_precision(::Type{T}) where {T<:AbstractFloat} = eps(T)
relative_precision(::Type{T}) where {T<:Real} = eps(Float64)
relative_precision(::Type{AbstractFloat}) = eps(Float64)

"""
    to_same_concrete_type(x1, x2, ...) -> xp1, xp2, ...

converts instances `x1`, `x2`, ... to the same concrete type. This method may be used
instead of `promote(x1,x2,...)` which does not warrant that the converted values have the
same type nor that is is concrete.

If all arguments are types, then:

    to_same_concrete_type(T1::Type, T2::Type, ...) -> T::Type

yields `T = promote_type(T1, T2, ...)` throwing an exception if `T` is not a concrete
type.

""" to_same_concrete_type

to_same_concrete_type() = throw(ArgumentError("no instance(s)/type(s) specified"))

# Implement `to_same_concrete_type` for 1, 2, and more instances.

function to_same_concrete_type(x::T) where {T}
    isconcretetype(T) || throw_not_concrete_type(T)
    return x
end

function to_same_concrete_type(x1::T, x2::T) where {T}
    isconcretetype(T) || throw_not_concrete_type(T)
    return x1, x2
end
function to_same_concrete_type(x1::T1, x2::T2) where {T1,T2}
    T = to_same_concrete_type(T1, T2)
    return as(T, x1), as(T, x2)
end
@inline function to_same_concrete_type(xs::T...) where {T}
    isconcretetype(T) || throw_not_concrete_type(T)
    return xs
end
@inline function to_same_concrete_type(xs...)
    T = to_same_concrete_type(map(typeof, xs...)...)
    return map(as(T), xs...)
end

# Implement `to_same_concrete_type` for 1, 2, and more types.

function to_same_concrete_type(::Type{T}) where {T}
    isconcretetype(T) || throw_not_concrete_type(T)
    return T
end

function to_same_concrete_type(::Type{T}, ::Type{T}) where {T}
    isconcretetype(T) || throw_not_concrete_type(T)
    return T
end

function to_same_concrete_type(::Type{T1}, ::Type{T1}) where {T1,T2}
    T = promote_type(T1, T2)
    isconcretetype(T) || throw_no_common_concrete_type(T1, T2)
    return T
end

@inline function to_same_concrete_type(::Type{T}...) where {T}
    isconcretetype(T) || throw_not_concrete_type(T)
    return T
end

@inline function to_same_concrete_type(Ts::Type...)
    T = promote_type(Ts...)
    isconcretetype(T) || throw_no_common_concrete_type(Ts...)
    return T
end

@noinline throw_not_concrete_type(T::Type) =
    throw(ArgumentError("type `$T` is not a concrete type"))

@noinline throw_no_common_concrete_type(T1::Type, T2::Type) =
    throw(ArgumentError("types `$T1` and `$T2` cannot be converted to a common concrete type"))

@noinline throw_no_common_concrete_type(Ts::Type...) =
    throw(ArgumentError(*("types `", join(Ts, "`, `", "`, and `"),
                          "` cannot be converted to a common concrete type")))

"""
    PlanarFields.nearly_equal(x, y) -> bool

yields whether numbers `x` and `y` are nearly equal relatively to the worst precision of
the two.

"""
function nearly_equal(x::T, y::T) where {T<:Number}
    if x == y
        return true
    else
        R = bare_type(T)
        if R <: Union{Integer,Rational}
            return false
        else
            prec = relative_precision(R)
            return abs(x - y) ≤ prec*max(abs(x), abs(y))
        end
    end
end

function nearly_equal(x::Tx, y::Ty) where {Tx<:Number,Ty<:Number}
    if x == y
        return true
    else
        Rx, Ry = bare_type(Tx), bare_type(Ty)
        if (Rx <: Union{Integer,Rational}) & (Ry <: Union{Integer,Rational})
            return false
        else
            # use the worst relative precision of the two
            prec = max(relative_precision(Rx), relative_precision(Ry))
            return abs(x - y) ≤ prec*max(abs(x), abs(y))
        end
    end
end
