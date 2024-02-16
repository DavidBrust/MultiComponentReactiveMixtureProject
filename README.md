# MultiComponenReactiveMixtureProject
Project package for the implementation of a photo-thermal catalytic reactor model.

A 3D (square prismatic) domain of porous material is modelled.
Steady state is considered, modelling of time dependence is considered as possible future option. The physical processes considered in the model include heat tranport including irradiation exchange and has phase species mass transport through the porous material including chemical reaction.
Model equations and correlations for effective properties of heat transfer through the porous material are taken from 
[VDI Heat Atlas](https://link.springer.com/referencework/10.1007/978-3-540-77877-6)
(chapter D6.3).
The mass transport through the porous material is described with the Dusty Gas Model (DGM).
A simplified form, neglecting thermal diffusion and external forces, as presented in 
[Veldsink 1995](https://doi.org/10.1016/0923-0467(94)02929-6) is used.

The simulation utilizes the finite volume method as implemented in the
[VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl) package.

The code development takes place in Pluto notebooks, located in the notebooks 
subfolder.

The most comprehensive Pluto notebook [PorousCatalystHot3DTopFlowIrrExchange_NonAlloc.jl](https://github.com/DavidBrust/MultiComponenReactiveMixtureProject/blob/main/notebooks/PorousCatalystHot3DTopFlowIrrExchange_NonAlloc.jl) considers irradiation exchange together with plotting and post-processing.