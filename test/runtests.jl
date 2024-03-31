using Test
using MultiComponentReactiveMixtureProject
using Pluto
using ExampleJuggler

ExampleJuggler.verbose!(true)


notebooks = [
    "Residual_MoleFrac1D.jl",
    "Uphill_Diff1D.jl",
    "PTReactorDemo.jl"
    ]

@testset "pluto notebooks" begin
    @testplutonotebooks(joinpath(@__DIR__, "..", "notebooks"), notebooks)
end


modules = ["PTR3D.jl"]

@testset "module examples" begin
    @testmodules(joinpath(@__DIR__, "..", "scripts"), modules)
end