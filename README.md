# FixedBed
Project package for implementation of a fixed bed (reactor) model.

In particular, the 2D, quasi-homogeneous steady state models for heat and mass transfer through fixed beds or porous materials are to be implemented. The model equations and correlations for effective properties used within the models is presented in:
[VDI Heat Atlas](https://link.springer.com/referencework/10.1007/978-3-540-77877-6)
(chapters M7, D6.3)


The simulation will utilize the finite volume method as implemented in the
[VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl) package.

The code development takes place in Pluto notebooks, located in the notebooks 
subfolder.