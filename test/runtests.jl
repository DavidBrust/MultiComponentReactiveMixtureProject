using Test
using FixedBed
using ExampleJuggler

ExampleJuggler.verbose!(true)

example_dir = joinpath(@__DIR__, "..", "notebooks/Demo")

# modules = ["ExampleModule.jl"]
notebooks = ["Residual_MoleFrac1D.jl","Uphill_Diff1D.jl","PCReactorDemo.jl"]
# scripts = ["testscript.jl", "PlutoTemplate.jl", "ExamplePluto.jl"]

@testset "pluto notebooks" begin
    @testplutonotebooks(example_dir, notebooks)
end

# @testset "module examples" begin
#     @testmodules(example_dir, modules, a=2)
# end

# @testset "scripts + notebooks" begin
#     @testscripts(example_dir, scripts)
# end