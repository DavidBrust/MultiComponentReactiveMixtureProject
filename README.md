MultiComponentReactiveMixtureProject
====================================

Project package for the implementation of a model for heat and multi-component, reactive gas phase transport in porous media.

The model is applied to the simulation of a Photo-Thermal catalytic reactor. 
Concentrated solar irradiation is used to supply heat input into the reactor to drive endothermal chemical reactions.
The reactor utilises a catalyst coated porous material which is directly irradiated and through which the reactive gas mixture passses.
A 3D (square prismatic) domain of porous material is modelled.
Model equations and correlations for effective properties of heat transfer through the porous material are taken from 
[VDI Heat Atlas](https://link.springer.com/referencework/10.1007/978-3-540-77877-6)
(chapter D6.3).
The simulation utilizes the finite volume method as implemented in the
[VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl) package.

After checking out or updating, run from the project root directory:
```
$ julia --project -e "using Pkg; Pkg.instantiate()"
```

Computational notebooks based on the ``Pluto.jl`` package are used to illustrate this package's functionality,
where notebooks located in the notebooks subfolder.

#### Notebooks

- `documentation.jl`: project documentation
- `Demo/Uphill_Diff1D.jl`: example of transport uphill a concentration gradient in iso-thermal conditions
- `Demo/Residual_MoleFrac1D.jl`: investigation of grid refinement on residual molar fractions
- `Demo/PTReactorDemo.jl.jl`: example application of photo-thermal reactor in 2D and 3D geometry and non-isothermal conditions
