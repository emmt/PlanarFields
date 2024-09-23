using PlanarFields
using TwoDimensional
using Test

@testset "PlanarFields" begin
    @testset "GridAxis" begin
        stp, rng = 0.1, -2:3
        A = @inferred GridAxis(rng, stp)
        @test A isa AbstractRange{typeof(stp)}
        @test eltype(A) === typeof(stp)
        @test length(A) == length(rng)
        @test eachindex(A) === rng
        @test axes(A) === (rng,)
        @test size(A) === (length(rng),)
        @test step(A) === stp
        @test IndexStyle(A) === IndexLinear()
        @test IndexStyle(typeof(A)) === IndexLinear()
        @test firstindex(A) === first(rng)
        @test lastindex(A) === last(rng)
        @test first(A) === step(A)*firstindex(A)
        @test last(A) === step(A)*lastindex(A)
        @test collect(A) == step(A)*collect(eachindex(A))
        @test A === @inferred GridAxis(rng; step=stp)
        @test A === @inferred GridAxis(stp, rng)
        @test A === @inferred GridAxis(stp, Int16(first(rng)):Int16(last(rng)))
        @test coord_type(A) === typeof(step(A))
        @test coord_type(typeof(A)) === typeof(step(A))

        @test A === @inferred GridAxis(A)
        @test A === @inferred GridAxis{eltype(A)}(A)
        @test A === @inferred copy(A)
        @test A === @inferred convert(GridAxis, A)
        @test A === @inferred convert(GridAxis{eltype(A)}, A)

        B = @inferred GridAxis{Float32}(A)
        @test B === @inferred convert(GridAxis{Float32}, A)
        @test B === @inferred GridAxis{Float32}(rng; step=stp)
        @test B === @inferred GridAxis{Float32}(rng, stp)
        @test B === @inferred GridAxis{Float32}(stp, rng)
        @test B isa AbstractRange{Float32}
        @test eltype(B) === Float32
        @test length(B) == length(A)
        @test eachindex(B) === eachindex(A)
        @test step(B) === eltype(B)(step(A))
        @test coord_type(B) === eltype(B)
        @test coord_type(typeof(B)) === eltype(B)
        @test collect(B) == step(B)*collect(eachindex(B))

        # Element type conversion.
        @test A === @inferred map(eltype(A), A)
        @test A === eltype(A).(A)
        @test B === @inferred map(eltype(B), A)
        @test B === eltype(B).(A)

        # Element-wise multiplication.
        C = π*B
        @test C === GridAxis(eachindex(B), π*step(B))
        @test C === B*π
        @test C === π .* B
        @test C === B .* π

        # Element-wise division.
        C = B/π
        @test C === GridAxis(eachindex(B), step(B)/π)
        @test C === π\B
        @test C === B ./ π
        @test C === π .\ B

        # Element-wise addition.
        off = 5
        @test_throws Exception B .+ π
        C = B .+ off*step(B)
        @test C === GridAxis(step(B), eachindex(B) .+ off)
        @test C === off*step(B) .+ B

        # Unary plus and minus.
        @test A === @inferred +A
        C = @inferred -A
        @test typeof(C) === typeof(A)
        @test step(C) === -step(A)
        @test eachindex(C) === eachindex(A)

        # Element-wise subtraction.
        off = 4
        @test_throws Exception B .- π
        C = B .- off*step(B)
        @test C === GridAxis(step(B), eachindex(B) .- off)
        @test C === -off*step(B) .+ B
        C = off*step(B) .- B
        @test C === GridAxis(-step(B), eachindex(B) .- off)
        @test C === -B .+ off*step(B)
    end

    @testset "Grid" begin
        stp, xrng, yrng = 0.13, -1:4, -2:5
        A = @inferred Grid(stp, xrng, yrng)
        @test A.X === GridAxis(stp, xrng)
        @test A.Y === GridAxis(stp, yrng)
        @test A isa AbstractMatrix
        @test eltype(A) === Point{typeof(stp)}
        @test axes(A) === (xrng, yrng)
        @test size(A) === (length(xrng), length(yrng))
        @test length(A) === length(xrng)*length(yrng)
        @test IndexStyle(A) === IndexCartesian()
        @test IndexStyle(typeof(A)) === IndexCartesian()
        @test first(CartesianIndices(A)) === CartesianIndex(first(xrng), first(yrng))
        @test last(CartesianIndices(A)) === CartesianIndex(last(xrng), last(yrng))
        @test first(A) === Point(step(A)*first(xrng), step(A)*first(yrng))
        @test last(A) === Point(step(A)*last(xrng), step(A)*last(yrng))
        @test coord_type(A) === typeof(step(A))
        @test coord_type(typeof(A)) === typeof(step(A))
    end
end
