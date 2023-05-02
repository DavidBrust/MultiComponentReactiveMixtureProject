module MWE

# after loading module MWE, run:
# sol,grid,sys,geom=MWE.main();

# plot Temperature field
# MWE.plot_T(sol,geom)

# calculate heat flux over symmetry boundaries, which should be 0
# MWE.FluxesSymmetryBC(sol,sys)


using VoronoiFVM
using ExtendableGrids, GridVisualize
using LessUnitful
using PyPlot

const iT=1 # index of temperature equation

const Γ_side_front = 1 # symmetry bc
const Γ_side_right = 2 # wall bc
const Γ_side_back = 3 # wall bc
const Γ_side_left = 4 # symmetry bc
const Γ_bottom = 5 # inflow bc
const Γ_top = 6 # top surface

const Tamb=298.15*ufac"K" # ambient temperature


function prism_sq(geom; nref=0,)
    (;wi,h) = geom
    hw=wi/2.0/10.0*2.0^(-nref)
    W=collect(0:hw:(wi/2.0))
    hh=h/10.0*2.0^(-nref)
    H=collect(0:hh:h)

    simplexgrid(W,W,H)
end
const grid_fun = prism_sq




function top(f,u,bnode)

    if bnode.region==Γ_top
        Glamp=100.0*ufac"kW/m^2" # irradiation flux density
        abs=0.3 # absorptivity for visible light
        emis=0.3 # emissivity
        flux_irrad = -emis*ph"σ"*(u[iT]^4-Tamb^4) + abs*Glamp
        f[iT] = -flux_irrad
    end
end

function side(f,u,bnode)
    if bnode.region==Γ_side_back || bnode.region==Γ_side_right # outer sides
        α_w=20.0*ufac"W/(m^2*K)" # wall heat transfer coefficient
        f[iT] = α_w*(u[iT]-Tamb)
    end
    # Γ_side_front and Γ_side_left are symmetry boundary conditions and
    # should have homogeneous neumann boundary conditions by default

end

function bottom(f,u,bnode)
    if bnode.region==Γ_bottom # bottom boundary
        emis=0.3
        σ=ph"σ"
        flux_irrad = -emis*σ*(u[iT]^4-Tamb^4)
        f[iT] = -flux_irrad
    end
end



function bcond(f,u,bnode)
    top(f,u,bnode)
    bottom(f,u,bnode)
    side(f,u,bnode)
end



function flux(f,u,edge)
    λ=1.0*ufac"W/(m*K)"
    f[iT]= λ*(u[iT,1]-u[iT,2])
end


Base.@kwdef mutable struct GeometryData

    # prism / 3D
    h::Float64=0.5*ufac"cm"
    wi::Float64=12.0*ufac"cm" # width/side lenght
    le::Float64=wi # prism width/side lenght

end;


function main(;geom=GeometryData())

    grid=grid_fun(geom)

    sys=VoronoiFVM.System(  grid;
                            flux=flux,
                            bcondition=bcond
                            )
    enable_species!(sys; species=[iT])
    inival=unknowns(sys)
    inival[iT,:] .= Tamb

    sol=solve(sys;inival=inival)
    sol,grid,sys,geom
end;

#### Post-processing
function plot_T(sol,geom)
    function TopPlane(sol,geom)

        grid=grid_fun(geom)
        (;wi,h)=geom
        w=wi/2

        bid = maximum(grid[BFaceRegions])+1
        bfacemask!(grid, [0,0,h],[w,w,h],bid)

        # keep x-y coordinates of parent grid
        function trans32(a,b)
            a[1]=b[1]
            a[2]=b[2]
        end
        grid_2D  = subgrid(grid, [bid], boundary=true, transform=trans32)

        sol_p = []
        for i=1:iT
            sol_i = view(sol[i, :], grid_2D)
            push!(sol_p, collect(sol_i))
        end
        sol_p, grid_2D
    end

    function CutPlane(sol,geom)


        grid=grid_fun(geom)
        (;wi,h)=geom
        w=wi/2

        bid = maximum(grid[BFaceRegions])+1
        bfacemask!(grid, [0,0.024,0],[w,0.024,h],bid)

        # keep x-z coordinates of parent grid
        function trans32(a,b)
            a[1]=b[1]
            a[2]=b[3]
        end
        grid_2D  = subgrid(grid, [bid], boundary=true, transform=trans32)

        sol_p = []
        for i=1:(iT)
            sol_i = view(sol[i, :], grid_2D)
            push!(sol_p, collect(sol_i))
        end
        sol_p, grid_2D
    end

    sol_xy, grid_xy = TopPlane(sol,geom)
    sol_xz, grid_xz = CutPlane(sol,geom)
    vis=GridVisualizer(layout=(1,2), resolution=(700, 300),Plotter=PyPlot)
    scalarplot!(vis[1,1], grid_xy, sol_xy[iT] .-273.15, colormap= :inferno, show=true)
    scalarplot!(vis[1,2], grid_xz, sol_xz[iT] .-273.15, aspect=4, colormap= :inferno, show=true)
end

function FluxesSymmetryBC(sol,sys)
    tff=TestFunctionFactory(sys)
    Γ_where_T_equal_1=[Γ_side_front,Γ_side_left] # front & left should be symmetry bcs.
    Γ_where_T_equal_0=[Γ_side_right,Γ_side_back,Γ_bottom,Γ_top]
    Tf=testfunction(tff,Γ_where_T_equal_0,Γ_where_T_equal_1)
    integrate(sys,Tf,sol)
end

function FluxesOuterSides(sol,sys)
    tff=TestFunctionFactory(sys)
    Γ_where_T_equal_1=[Γ_side_right,Γ_side_back] # front & left should be symmetry bcs.
    Γ_where_T_equal_0=[Γ_side_front,Γ_side_left,Γ_bottom,Γ_top]
    Tf=testfunction(tff,Γ_where_T_equal_0,Γ_where_T_equal_1)
    integrate(sys,Tf,sol)
end

end
