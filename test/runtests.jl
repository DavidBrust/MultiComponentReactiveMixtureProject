using Test
using MultiComponentReactiveMixtureProject
using ExampleJuggler

ExampleJuggler.verbose!(true)


notebooks = [
    "Residual_MoleFrac1D.jl",
    "Uphill_Diff1D.jl",
    "PTReactorDemo.jl"
    ]

# @testset "pluto notebooks" begin
#     @testplutonotebooks(joinpath(@__DIR__, "..", "notebooks"), notebooks)
# end

# scripts = ["TestSim.jl"]

# @testset "scripts + notebooks" begin
#     @testscripts(joinpath(@__DIR__, "..", "scripts"), scripts)
# end

modules = ["Test3D.jl"]

@testset "module examples" begin
    @testmodules(joinpath(@__DIR__, "..", "scripts"), modules)
end