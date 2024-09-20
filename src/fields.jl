"""
    A = PlanarField{T}(grid)
    A = PlanarField{T,S}(grid)

Allocate an object `A` for storing the values of a field sampled on the 2-dimensional
planar `grid` of equally spaced nodes. `T` is the type of the field values, `S` is the
type of node coordinates (the same as for `grid` by default). The values stored by the
returned array are undefined.

Another posibility is to specify any arguments and/or keywords that are accepted by the
[`Grid`](@ref) constructor. For example:

    A = PlanarField{T[,S]}(X, Y)
    A = PlanarField{T[,S]}(step, I, J)
    A = PlanarField{T[,S]}(I, J, step)
    A = PlanarField{T[,S]}(I, J; step)

where `X` and `Y` are the grid abscissae and ordinates, or where `I`, `J`, and `step` are
the index ranges along the array dimensions and `step` is the grid node spacing. The type
`T` of the field values must be specified. If the type `S` of the node coordinates is not
specified, it is inferred for that of `X` and `Y` or from `step`.

A planar field object `A` is an abstract 2-dimensional array whose values can be accessed
by `A[I]`. The position of a node at Cartesian index `I` is given by `step*Point(I)`.

"""
PlanarField{T}(grid::Grid) where {T} =
    PlanarField(OffsetArray(Array{T}(undef, size(grid)), axes(grid)), step(grid))

# Build a planar grid from arguments and keywords.
PlanarField{T}(args...; kwds...) where {T} = PlanarField{T}(Grid(args...; kwds...))
PlanarField{T,S}(args...; kwds...) where {T,S<:Number} = PlanarField{T}(Grid{S}(args...; kwds...))

"""
    A = PlanarField{T=eltype(vals),S=typeof(step)}(vals, step)
    A = PlanarField{T=eltype(vals),S=typeof(step)}(vals; step)

Build an object `A` representing a field sampled on a 2-dimensional planar grid of equally
spaced nodes. `vals` is the 2-dimensional array storing the sampled field values and
`step` is the node spacing. Both `vals` and `step` may have units.

If the specified type `T` of the field values is the same as `eltype(vals)` (the default),
the returned object `A` and `vals` share their values; otherwise, the values of `A` are
stored in a different array than `vals`. The method `parent(A)` yields the array
backing the storage of the values in `A`.

"""
PlanarField(vals::AbstractMatrix; step::Number) = PlanarField(vals, step)

PlanarField{T}(vals::AbstractMatrix; step::Number) where {T} = PlanarField{T}(vals, step)
PlanarField{T}(vals::AbstractMatrix{T}, step::Number) where {T} = PlanarField(vals, step)
PlanarField{T}(vals::AbstractMatrix, step::Number) where {T} =
    PlanarField(copyto!(similar(vals, T), vals), step)

PlanarField{T,S}(vals::AbstractMatrix; step::Number) where {T,S<:Number} =
    PlanarField{T,S}(vals, step)
PlanarField{T,S}(vals::AbstractMatrix, step::Number) where {T,S<:Number} =
    PlanarField{T}(vals, as(T, step))

# Copy/convert constructors.
PlanarField(A::PlanarField) = A
PlanarField{T}(A::PlanarField{T}) where {T} = A
PlanarField{T}(A::PlanarField) where {T<:Number} = PlanarField{T}(parent(A), step(A))
PlanarField{T,S}(A::PlanarField{T,S}) where {T,S<:Number} = A
PlanarField{T,S}(A::PlanarField) where {T,S<:Number} = PlanarField{T,S}(parent(A), step(A))
Base.convert(::Type{T}, A::T) where {T<:PlanarField} = A
Base.convert(::Type{T}, A::PlanarField) where {T<:PlanarField} = T(A)
Base.copy(A::PlanarField) = PlanarField(copy(parent(A)), step(A))

# The base similar method takes an additional step keyword for a planar field.
Base.similar(A::PlanarField; kwds...) = similar(A, eltype(A), axes(A); kwds...)
Base.similar(A::PlanarField, ::Type{T}; kwds...) where {T} = similar(A, T, axes(A); kwds...)
Base.similar(A::PlanarField, inds::NTuple{2,Union{Integer,ArrayAxisLike}}; kwds...) =
    similar(A, eltype(A), inds; kwds...)
function Base.similar(A::PlanarField, ::Type{T}, inds::NTuple{2,Union{Integer,ArrayAxisLike}};
                      step::Number = step(A)) where {T}
    return PlanarField{T}(Grid(step, to_axes(dims)))
end

# Retrieve the grid of the nodes.
Grid(A::PlanarField) = Grid(step(A), axes(A))
Grid{S}(A::PlanarField) where {S<:Number} = Grid{S}(step(A), axes(A))

# Accessors.
Base.parent(A::PlanarField) = getfield(A, :vals)
Base.step(A::PlanarField) = getfield(A, :step)

# Traits.
TwoDimensional.coord_type(A::PlanarField) = coord_type(typeof(A))
TwoDimensional.coord_type(::Type{<:PlanarField{T,S}}) where {T,S} = S

# Abstract array API for PlanarField objects.
Base.length(A::PlanarField) = prod(size(A))
Base.size(A::PlanarField) = map(length, axes(A))
Base.axes(A::PlanarField) = axes(parent(A))
Base.IndexStyle(::Type{PlanarField{T,S,true,A}}) where {T,S,A} = IndexLinear()
Base.IndexStyle(::Type{PlanarField{T,S,false,A}}) where {T,S,A} = IndexCartesian()

@inline function Base.getindex(A::PlanarField{T,S,true}, i::Int) where {T,S}
    @boundscheck checkbounds(A, i)
    return @inbounds getindex(parent(A), i)
end

@inline function Base.getindex(A::PlanarField{T,S,false}, I::Vararg{Int,2}) where {T,S}
    @boundscheck checkbounds(A, I...)
    return @inbounds getindex(parent(A), I...)
end

@inline function Base.setindex!(A::PlanarField{T,S,true}, x, i::Int) where {T,S}
    @boundscheck checkbounds(A, i)
    @inbounds setindex!(parent(A), x, i)
    return A
end

@inline function Base.setindex!(A::PlanarField{T,S,false}, x, I::Vararg{Int,2}) where {T,S}
    @boundscheck checkbounds(A, I...)
    @inbounds setindex!(mask(A), x, I...)
    return A
end
