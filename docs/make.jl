push!(LOAD_PATH,"../src/")
using EquationDemoLib
# https://medium.com/coffee-in-a-klein-bottle/creating-and-deploying-your-julia-package-documentation-1d09ddc90474
using Documenter
makedocs(
        sitename = "EquationDemoLib.jl",
        modules  = [EquationDemoLib],
        remotes = nothing,
        pages=[
                "Home" => "index.md"
        ])