module Example1_HeatFluxCalc

using FixedBed
using VoronoiFVM
using ExtendableGrids, GridVisualize
using LessUnitful
using Plots

function prism_sq(data;nref=0, w=data.wi, h=data.h)
	
	hw=w/2.0/10.0*2.0^(-nref)
	#hl=l/2.0/10.0*2.0^(-nref)
	hh=h/10.0*2.0^(-nref)
	W=collect(0:hw:(w/2.0))
    #L=collect(0:hl:(l/2.0))
    H=collect(0:hh:h)
	
	simplexgrid(W,W,H)	
end


function main3D(data;nref=0)
	iT=data.iT

	# function return 3D velocity vector: flow upward in z-direction
    function fup(x,y,z)
        return 0,0,-data.u0
    end    
	
	function flux(f,u,edge,data)
		(;Fluid,u0,p)=data
		Tbar=0.5*(u[iT,1]+u[iT,2])
		ρf=density_idealgas(Fluid, Tbar, p)
		cf=heatcap_gas(Fluid, Tbar)
		λf=thermcond_gas(Fluid, Tbar)
		λbed=kbed(data)*λf
		
		vh=project(edge,(0,0,u0))
		conv=vh*ρf*cf/λbed
		Bp,Bm = fbernoulli_pm(conv)
		#f[iT]= λbed*(Bm*u[iT,1]-Bp*u[iT,2])
		f[iT]= λbed*(Bm*(u[iT,1]-data.Tamb)-Bp*(u[iT,2]-data.Tamb))		
	end

    function flux_rerad_top(f,u,bnode,data)
        if bnode.region==6 # top boundary
            f[iT] = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)
        end
    end

    function flux_lamp_abs(f,u,bnode,data)
        if bnode.region==6 # top boundary
            f[iT] = -data.Abs_lamp*data.G_lamp
        end
    end

    function flux_convec_top(f,u,bnode,data)
        if bnode.region==6 # top boundary
            ρf=density_idealgas(data.Fluid, u[iT], data.p)
            cf=heatcap_gas(data.Fluid, u[iT])
            f[iT] = data.u0*ρf*cf*(u[iT]-data.Tamb)
        end
    end

	function top(f,u,bnode,data)
		if bnode.region==6 # top boundary
			flux_rerad = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)
			#flux_rerad = 0
			
			ρf=density_idealgas(data.Fluid, u[iT], data.p)
			cf=heatcap_gas(data.Fluid, u[iT])
            flux_convec = data.u0*ρf*cf*(u[iT]-data.Tamb) # flow velocity is normal to top boundary
			
			f[iT] = -data.Abs_lamp*data.G_lamp + flux_rerad + flux_convec
		end
	end

	function bottom(f,u,bnode,data)
		if bnode.region==5 # bottom boundary
            f[iT] = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)
		end
	end

    function sidewalls(f,u,bnode,data)
        # boundary conditions at side walls
        #if bnode.region==2
            boundary_robin!(f,u,bnode;species=iT,region=2, factor=data.α_w, value=data.Tamb*data.α_w)
        #elseif bnode.region==3
            boundary_robin!(f,u,bnode;species=iT,region=3, factor=data.α_w, value=data.Tamb*data.α_w)
        #end
    end

	function bcondition(f,u,bnode,data)
		bottom(f,u,bnode,data)
		boundary_robin!(f,u,bnode;species=iT,region=2, factor=data.α_w, value=data.Tamb*data.α_w)
		boundary_robin!(f,u,bnode;species=iT,region=3, factor=data.α_w, value=data.Tamb*data.α_w)
		top(f,u,bnode,data)
	end
	

	
	#grid=prism(;nref=nref, w=data.wi, l=data.le, h=data.h)
	grid=prism_sq(data;nref=nref,)
	
	sys=VoronoiFVM.System(grid;
                          data=data,
                          flux=flux,
                          breaction=bcondition,
                          species=[iT],
                          )
	inival=unknowns(sys)
	inival[iT,:] .= map( (x,y,z)->(data.Tamb+500*z/data.h),grid)
	#inival[iT,:] .= data.Tamb
	sol=solve(inival,sys)
	sys,sol,data,nref



    ############################################################################
    ## analysis    
    # q_full=integrate(sys,sys.physics.breaction,sol; boundary=true)

    # q_lamp_abs=integrate(sys,flux_lamp_abs,sol; boundary=true)[6]
    # q_rerad_top=integrate(sys,flux_rerad_top,sol; boundary=true)[6]
    # q_convec_top=integrate(sys,flux_convec_top,sol; boundary=true)[6]
    # q_sidewalls=integrate(sys,sidewalls,sol; boundary=true)[1,[2,3]]

    # q_lamp_abs,q_rerad_top,q_convec_top,q_sidewalls,q_full
end


function HeatFluxes(C=[1,10])
    data=ModelData()
    iT=data.iT

    lC=length(C)
    #q_full=zeros(Float64, lC)
    q_in=zeros(Float64, lC)
    q_top_abs=zeros(Float64, lC)
    q_top_refl=zeros(Float64, lC)
    q_top_rerad=zeros(Float64, lC)
    q_top_convec=zeros(Float64, lC)
    q_sides_conv=zeros(Float64, lC)
    q_bot_rerad=zeros(Float64, lC)

    function flux_rerad_top(f,u,bnode,data)
        if bnode.region==6 # top boundary
            f[iT] = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)
        end
    end

    function flux_lamp_abs(f,u,bnode,data)
        if bnode.region==6 # top boundary
            f[iT] = -data.Abs_lamp*data.G_lamp
        end
    end

    function flux_convec_top(f,u,bnode,data)
        if bnode.region==6 # top boundary
            ρf=density_idealgas(data.Fluid, u[iT], data.p)
            cf=heatcap_gas(data.Fluid, u[iT])
            f[iT] = data.u0*ρf*cf*(u[iT]-data.Tamb)
        end
    end

    function sidewalls(f,u,bnode,data)
        # boundary conditions at side walls
            boundary_robin!(f,u,bnode;species=iT,region=2, factor=data.α_w, value=data.Tamb*data.α_w)
            boundary_robin!(f,u,bnode;species=iT,region=3, factor=data.α_w, value=data.Tamb*data.α_w)
    end

    function bottom(f,u,bnode,data)
		if bnode.region==5 # bottom boundary
            f[iT] = data.Eps_ir*ph"σ"*(u[iT]^4 - data.Tamb^4)
		end
	end

    for (i,C) in enumerate(C)
        data.G_lamp=C*1.0*ufac"kW/m^2"
        sys,sol,data,nref = main3D(data,nref=0)
        # factor 4: due to symmetry, sim domain covers only 1/4 of the frit
        q_top_abs[i] = -4*integrate(sys,flux_lamp_abs,sol; boundary=true)[6]
        q_in[i] = q_top_abs[i]/data.Abs_lamp
        q_top_refl[i] = (1.0-data.Abs_lamp)*q_in[i]
        q_top_rerad[i] = 4*integrate(sys,flux_rerad_top,sol; boundary=true)[6]
        q_top_convec[i] = 4*integrate(sys,flux_convec_top,sol; boundary=true)[6]
        q_sides_conv[i] = 4*sum(integrate(sys,sidewalls,sol; boundary=true)[1,[2,3]])
        q_bot_rerad[i] = 4*integrate(sys,bottom,sol; boundary=true)[5]

    end
    q_in,q_top_abs,q_top_refl,q_top_rerad,q_top_convec,q_sides_conv,q_bot_rerad
    
end

function PlotLosses(C=[1,10,25,50,75,100])
    q_in,q_top_abs,q_top_refl,q_top_rerad,q_top_convec,q_sides_conv,q_bot_rerad = HeatFluxes(C)
    p=plot(xguide="Irradiation / kW m⁻²", yguide="Loss contributions / -", legend=:outertopright, size=(400,250))
	areaplot!(C, [q_top_refl q_top_rerad q_top_convec q_sides_conv q_bot_rerad] ./ (q_in), fillalpha = [0.2 0.2], seriescolor = [1 2 3 4 5], label=["Refl Top" "Rerad Top" "Conv Top" "Conv Sides" "Rerad Bot"])
    savefig("img/out/loss_contributions.svg")

end

end
