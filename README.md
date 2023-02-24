# FixedBed
Project package for implementation of a catalytic reactor model.

A 2D (axisymmetric) or 3D (prismatic) domain of porous material is modelled.
Steady state is considered as well as heat and mass transfer through the porous material.
Model equations and correlations for effective properties of heat transfer through the porous
material are taken from 
[VDI Heat Atlas](https://link.springer.com/referencework/10.1007/978-3-540-77877-6)
(chapter D6.3).
The mass transport through the porous material is described with the Dusty Gas Model (DGM).
A simplified form, neglecting thermal diffusion and external forces, as presented in 
[Veldsink 1995](https://doi.org/10.1016/0923-0467(94)02929-6)
is used.


The simulation will utilize the finite volume method as implemented in the
[VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl) package.

The code development takes place in Pluto notebooks, located in the notebooks 
subfolder.