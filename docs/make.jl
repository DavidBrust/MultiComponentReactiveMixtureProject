using Documenter, ExampleJuggler
using MultiComponentReactiveMixtureProject

function make_all()
    ExampleJuggler.verbose!(true)

    cleanexamples()
    notebookdir = joinpath(@__DIR__, "..", "notebooks")

    notebooks = [
        "Project documentation" => "documentation.jl",
    ]

    notebook_examples = @docplutonotebooks(notebookdir, notebooks, iframe=false)

    makedocs(; sitename = "MultiComponentReactiveMixtureProject.jl",
            modules = [MultiComponentReactiveMixtureProject],
            checkdocs = :all,
            clean = false,
            doctest = false,
            warnonly = true,
            authors = "D. Brust",
            #  repo = "https://github.com/j-fu/TwoPhaseGeothermalProject.jl",
            repo = "https://github.com/DavidBrust/PhotoCatalyticReactor.git",
            pages = [
                 "Home" => "index.md",
                 "Notebooks" => notebook_examples,
             ])

    cleanexamples()

    if !isinteractive()
        # deploydocs(; repo = "github.com/j-fu/VoronoiFVM.jl.git")
        deploydocs(; repo = "github.com/DavidBrust/PhotoCatalyticReactor.git")
    end
end

if isinteractive()
    make_all(;)
else
    make_all(;)
end