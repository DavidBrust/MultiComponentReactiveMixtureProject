MultiComponentReactiveMixtureProject
====================================

[![DOI](https://zenodo.org/badge/643598052.svg)](https://zenodo.org/doi/10.5281/zenodo.10901335) 


Project package for the implementation of a model for heat and multi-component, reactive gas phase transport in porous media in ``Julia``. The model is presented in [https://doi.org/10.1016/j.cej.2025.162027](https://doi.org/10.1016/j.cej.2025.162027).

The model is applied to the simulation of a [photo-thermal catalytic reactor](https://doi.org/10.1016/j.jece.2024.113372). 
Concentrated solar irradiation is used to supply heat input into the reactor to drive endothermal chemical reactions.
The reactor utilizes a catalyst coated porous material which is directly irradiated and through which the reactive gas mixture passes.
A 3D (square prismatic) domain of porous material is modelled.

The simulation utilizes the finite volume method as implemented in the
[VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl) package.

In order to use it, check out the package with ``git``:
```
git clone https://github.com/DavidBrust/MultiComponentReactiveMixtureProject
```

After checking out or updating, run from the project root directory:
```
$ julia --project -e "using Pkg; Pkg.instantiate()"
```

Computational notebooks based on the ``Pluto.jl`` package are used to illustrate this package's functionality, located in the ``notebooks`` subfolder.

#### Notebooks

- `LoschmidtTernary.jl`: simulation of diffusion in ternary gas mixture (verification of isothermal model)
- `Thermodiffusion.jl`: demonstration of thermodiffusion (reproduces figures shown in the article)
- `Uphill_Diff1D.jl`: illustration of transport against a concentration gradient (convection driven) in isothermal conditions
- `Residual_MoleFrac1D.jl`: investigation of grid refinement on residual molar fractions
- `PTReactorDemo.jl`: application to photo-thermal reactor in 2D and 3D geometry and non-isothermal conditions

#### Examples shown in the article
Running the following files allows reproducing the figures shown in the article:
- `scripts/PTR3D.jl`: script for 3D Photo-thermal reactor simulation with data export for plot (Figure 5, Section 4.1)
- `Thermodiffusion.jl`: Pluto notebook demonstrating thermodiffusion, export plot (Figure 8, Section 4.2)

To run the script `PTR3D.jl`, run from the Julia REPL in the project root:
```
include("scripts/PTR3D.jl")
PTR3D.run()
```

The script will run the simulation and will export the solution in VTK format to a new folder for plotting with [ParaView](https://www.paraview.org/). Use the provided ParaView state file `data/PTR3D.pvsm` to reproduce Figure 5 in Section 4.1.

This project has received funding from the European Union's Horizon 2020 research and innovation
programme under grant agreement No 862453.

<img src="assets/flag_yellow_eps.svg" alt="Flag of the European Union" width="150"/>