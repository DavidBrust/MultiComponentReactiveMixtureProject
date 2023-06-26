module TestKin

using LessUnitful, GridVisualize, GLMakie, ExtendableGrids

include("../notebooks/PorousCatalystHot3DTopFlowIrrExchange_NonAlloc.jl")

function plotChemCS(sol, data)
	(;gni,ip)=data
	sol_xz, grid_xz = CutPlane(data,sol)
	scalarplot(grid_xz, sol_xz[gni[:CO]] ./sol_xz[ip], show=true, Plotter=GLMakie)
end



function runSim(;data=ModelData())
    function control(sys;kwargs...)
        SolverControl(
#			    direct_umfpack(),
            gmres_umfpack(),
    #			gmres_eqnblock_umfpack(),
    #			gmres_eqnblock_iluzero(),
            sys;
            verbose="na",
            log=true,
            reltol=1.0e-8,
            reltol_linear=1.0e-5,
            kwargs...
        )
    end

    sol_,grid,sys,data_embed=main(;data=data,control,
    assembly=:edgewise);
	if sol_ isa VoronoiFVM.TransientSolution
		sol = copy(sol_(sol_.t[end]))
	else
		sol = copy(sol_)
	end

    sol,grid,sys,data_embed
end



end