using LinearAlgebra

# G0_bot, IR/vis : surface radiosities of glass underside
function flux_window_underside(f,u,bnode,data)
    (;iT,top_radiation_boundaries)=data
    # if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit
    if bnode.region in top_radiation_boundaries        

        G1_bot_vis, G1_bot_IR = PTR_radiosity_window(f,u,bnode,data)

        f[iT] = G1_bot_vis + G1_bot_IR
    end
end


# # G2_IR/vis : surface radiosities of catalyst layer 
function flux_catalyst_layer(f,u,bnode,data)
    (;iT,uc_cat,top_radiation_boundaries)=data
    if bnode.region in top_radiation_boundaries
        
        # catalyst layer properties (2)
        rho2_vis=uc_cat.rho_vis
        rho2_IR=uc_cat.rho_IR
        eps2=uc_cat.eps

        G1_bot_vis, G1_bot_IR = PTR_radiosity_window(f,u,bnode,data)

        G2_vis = rho2_vis*G1_bot_vis

        G2_IR = eps2*ph"σ"*u[iT]^4 + rho2_IR*G1_bot_IR

        f[iT] = G2_vis + G2_IR
    end
end

function flux_convection_top(f,u,bnode,data)
    (;iT,iTw,X0,uc_h,Nu,top_radiation_boundaries)=data
    if bnode.region in top_radiation_boundaries
   
        Tm=0.5*(u[iT] +  u[iTw]) # mean temperature
        # thermal conductivity at Tm and inlet composition X0 (H2/CO2 = 1/1)
        @inline _,λf=dynvisc_thermcond_mix(data, Tm, X0)
        dh=2*uc_h
        kconv=Nu*λf/dh*ufac"W/(m^2*K)"
        # positive flux in positive z-coordinate
        hflux_conv = kconv*(Tm-u[iTw])	
        f[iT] = hflux_conv
    end
end

function flux_enth_top(f,u,bnode,data)
    (;iT, solve_T_equation, inflow_boundaries, T_gas_in, X0, mfluxin, mmix0 ) =data
	if solve_T_equation && bnode.region in inflow_boundaries
        
        hin = zero(eltype(u))
        @inbounds for i=1:ng
            hin += X0[i]*enthalpy_gas_thermal(data.Fluids[i],T_gas_in)
        end

        r_hf_enth = mfluxin/mmix0 * hin
        f[iT] = r_hf_enth
    end
end

function flux_enth_bottom(f,u,bnode,data)
    (;iT, solve_T_equation, outflow_boundaries, m, mfluxin, ) =data
	if solve_T_equation && bnode.region in outflow_boundaries
        
        mmix = zero(eltype(u))
        hout = zero(eltype(u))
        @inbounds for i=1:ng
            mmix += u[i]*m[i]
            hout += u[i]*enthalpy_gas_thermal(data.Fluids[i],u[iT])
        end

        # total mass is constant, outflow cross section = inflow cross section
        r_hf_enth = mfluxin/mmix * hout
        f[iT] = r_hf_enth
    end
end

# Post-processing (pp) wrapping function to handle time-dependent ramping inside b.c. functions
function PTR_top_post(f,u,bnode,data)
    bnode.time = Inf
    PTR_top(f,u,bnode,data)
end

function PTR_side_post(f,u,bnode,data)
    bnode.time = Inf
    PTR_side(f,u,bnode,data)
end

function PTR_bottom_post(f,u,bnode,data)
    bnode.time = Inf
    PTR_bottom(f,u,bnode,data)
end

function DMS_boutflow_post(f,u,edge,data)
    edge.time = Inf
    DMS_boutflow(f,u,edge,data)
end


function flux_radiation_frit_bottom(f,u,bnode,data)
    (;iT,iTp,lc_frit,lc_plate,bottom_radiation_boundaries) = data
    #if bnode.region==Γ_bottom
    if bnode.region in bottom_radiation_boundaries		
		
		# irradiation exchange between porous frit (3) and Al bottom plate (4)
		# porous frit properties (3)
		eps3 = lc_frit.eps
        rho3_IR = lc_frit.rho_IR
		# Al bottom plate properties (4)
		eps4 = lc_plate.eps
        rho4_IR = lc_plate.rho_IR
	
        T3 = u[iT] # frit temperature
        T4 = u[iTp] # plate temperature
        
	    G34_IR = (eps3*ph"σ"*T3^4 + rho3_IR*eps4*ph"σ"*T4^4)/(1-rho3_IR*rho4_IR)

        f[iT] = G34_IR
    end
end

function flux_radiation_plate_bottom(f,u,bnode,data)
    (;iT,iTp,lc_frit,lc_plate,bottom_radiation_boundaries) = data
    # if bnode.region==Γ_bottom
    if bnode.region in bottom_radiation_boundaries		
		
		# irradiation exchange between porous frit (3) and Al bottom plate (4)
		# porous frit properties (3)
		eps3 = lc_frit.eps
        rho3_IR = lc_frit.rho_IR
		# Al bottom plate properties (4)
		eps4 = lc_plate.eps
        rho4_IR = lc_plate.rho_IR		

        T3 = u[iT] # frit temperature
        T4 = u[iTp] # plate temperature

	    G34_IR = (eps3*ph"σ"*T3^4 + rho3_IR*eps4*ph"σ"*T4^4)/(1-rho3_IR*rho4_IR)

        G43_IR = rho4_IR*G34_IR + eps4*ph"σ"*T4^4

        f[iT] = G43_IR
    end
end

function flux_convection_bottom(f,u,bnode,data)
    (;iT,iTp,lc_h,Nu,bottom_radiation_boundaries)=data

    # if bnode.region==Γ_bottom        
    if bnode.region in bottom_radiation_boundaries
        ng=ngas(data)

        X=zeros(eltype(u), ng)
        # mole_frac!(bnode,data,X,u)
        MoleFrac!(X,u,data)

        # local plate temperature
        Tplate = u[iTp]
        Tm=0.5*(u[iT] + Tplate)
        @inline _,λf=dynvisc_thermcond_mix(data, Tm, X)

        dh=2*lc_h
        kconv=Nu*λf/dh*ufac"W/(m^2*K)"
        
        # positive flux in negative z coord. -> pointing towards bottom plate
        hflux_conv = kconv*(Tm-Tplate)	
        f[iT] = hflux_conv
   end
end

function flux_side(f,u,bnode,data)
	(;iT,k_nat_conv,delta_gap,Tamb,side_boundaries)=data
	ng=ngas(data)

	if bnode.region in side_boundaries
	    X=MVector{ng,eltype(u)}(undef)
        @inline MoleFrac!(X,u,data)
        @inline _,λf=dynvisc_thermcond_mix(data, u[iT], X)

        f[iT] = (u[iT]-Tamb)/(delta_gap/λf+1/k_nat_conv) 
    end
end


function HeatFluxes_EB_I(t,solt,sys,data)

    sol = solt(t)
    HeatFluxes_EB_I(sol,sys,data)

end

# stationary energy balance: dE_dt = 0
function HeatFluxes_EB_I(sol,sys,data)

    (;iT, inflow_boundaries, top_radiation_boundaries, outflow_boundaries, bottom_radiation_boundaries, side_boundaries) = data

    top_boundaries = unique(reduce(vcat, [inflow_boundaries,top_radiation_boundaries]))

    bottom_boundaries = unique(reduce(vcat, [bottom_radiation_boundaries,outflow_boundaries]))
    
    # Top boundary
    Q_irrad_12 = integrate(sys,flux_window_underside,sol; boundary=true)[iT,top_radiation_boundaries]
    Q_irrad_12 = sum(Q_irrad_12)

    Q_irrad_21 = integrate(sys,flux_catalyst_layer,sol; boundary=true)[iT,top_radiation_boundaries]
    Q_irrad_21 = sum(Q_irrad_21)

    H_enth_in_2 = integrate(sys,flux_enth_top,sol; boundary=true)[iT,top_boundaries]
    H_enth_in_2 = sum(H_enth_in_2)
    
    Q_conve_21 = integrate(sys,flux_convection_top,sol; boundary=true)[iT,top_radiation_boundaries]
    Q_conve_21 = sum(Q_conve_21)

    # Bottom boundary
    Q_irrad_34 = integrate(sys,flux_radiation_frit_bottom,sol; boundary=true)[iT,bottom_radiation_boundaries]
    Q_irrad_34 = sum(Q_irrad_34)
    
    Q_irrad_43 = integrate(sys,flux_radiation_plate_bottom,sol; boundary=true)[iT,bottom_radiation_boundaries]
    Q_irrad_43 = sum(Q_irrad_43)

    H_enth_out_3 = integrate(sys,flux_enth_bottom,sol; boundary=true)[iT,bottom_boundaries]
	H_enth_out_3 = sum(H_enth_out_3)

    Q_conve_34 = integrate(sys,flux_convection_bottom,sol; boundary=true)[iT,bottom_radiation_boundaries]
    Q_conve_34 = sum(Q_conve_34)

    # Side boundaries
    Q_sides = integrate(sys,flux_side,sol; boundary=true)[iT,side_boundaries]
    Q_sides = sum(Q_sides)

    H_chemical = sum(integrate(sys,sys.physics.reaction,sol), dims=2)[iT]
    # H_thermal = H_enth_in_2 - H_enth_out_3

    (
        Q_irrad_12=Q_irrad_12,
        Q_irrad_21=-Q_irrad_21,        
        Q_conve_21=-Q_conve_21,
        Q_irrad_34=-Q_irrad_34,
        Q_irrad_43=Q_irrad_43,        
        Q_conve_34=Q_conve_34,
        Q_sides=-Q_sides,
        H_enth_in_2=H_enth_in_2,
        H_enth_out_3=-H_enth_out_3,
        H_chemical=-H_chemical
    )
end

function HeatFluxes_EB_I_inner(t,solt,sys,data)
    sol = solt(t)
    HeatFluxes_EB_I_inner(sol,sys,data)
end

function HeatFluxes_EB_I_inner(sol,sys,data)

    (;iT, inflow_boundaries,top_radiation_boundaries, side_boundaries,bottom_radiation_boundaries, outflow_boundaries) = data

	top_boundaries = unique(reduce(vcat, [inflow_boundaries,top_radiation_boundaries]))

	bottom_boundaries = unique(reduce(vcat, [bottom_radiation_boundaries,outflow_boundaries]))
	
	hf_top = integrate(sys,MultiComponentReactiveMixtureProject.PTR_top_post,sol; boundary=true)[iT,top_boundaries]
	hf_top = -sum(hf_top)

	hf_top_enth = integrate(sys,MultiComponentReactiveMixtureProject.flux_enth_top,sol; boundary=true)[iT,top_boundaries]
	hf_top_enth = sum(hf_top_enth)
	
	hf_side = integrate(sys,MultiComponentReactiveMixtureProject.PTR_side_post,sol; boundary=true)[iT,side_boundaries]
	hf_side = sum(hf_side)

	hf_bottom_rad = integrate(sys,MultiComponentReactiveMixtureProject.PTR_bottom_post,sol; boundary=true)[iT,bottom_boundaries]
	hf_bottom_rad = sum(hf_bottom_rad)

	hf_bottom_enth = integrate(sys,MultiComponentReactiveMixtureProject.flux_enth_bottom,sol; boundary=true)[iT,bottom_boundaries]
	hf_bottom_enth = sum(hf_bottom_enth)

	hf_reaction = sum(integrate(sys,sys.physics.reaction,sol), dims=2)[iT]

	# # [hf_top, hf_top_enth, hf_side, hf_bottom_enth, hf_bottom_rad], hf_top-hf_side-hf_bottom_enth-hf_bottom_rad

	# Heatflows = Dict(
	# 	"Top radiation" => hf_top-hf_top_enth,
	# 	"Top enthalpy inflow" => hf_top_enth,
	# 	"Side convection" => -hf_side,
	# 	"Bottom radiation" => -hf_bottom_rad,
	# 	"Bottom enthalpy outflow" => -hf_bottom_enth,
	# 	"Reaction enthalpy" => -hf_reaction,
	# )

    (
        Q_rad_top_in = hf_top-hf_top_enth,
        H_thermal_top_in = hf_top_enth,
        Q_conv_sides_out = -hf_side,
        Q_rad_bottom_out = -hf_bottom_rad,
        H_thermal_bottom_out = -hf_bottom_enth,
        H_chemical_out = -hf_reaction,
    )
end

function BoundaryFlows_Integrals(sol, sys, data) 
	(;outflow_boundaries,inflow_boundaries)=data
	
	tfact=TestFunctionFactory(sys)	
	tf_out=testfunction(tfact,inflow_boundaries,outflow_boundaries)
	tf_in=testfunction(tfact,outflow_boundaries,inflow_boundaries)
		
    inflow_rate = integrate(sys, tf_in, sol)
    outflow_rate = integrate(sys, tf_out, sol)

    reaction_rate = sum(integrate(sys, sys.physics.reaction, sol), dims=2)
		
    return (
        inflow_rate,
        outflow_rate,
        reaction_rate
    )
end

function BoundaryFlows_Integrals(solt::TransientSolution, sys, data)
	(;outflow_boundaries,inflow_boundaries)=data
	
	tfact=TestFunctionFactory(sys)	
	tf_out=testfunction(tfact,inflow_boundaries,outflow_boundaries)
	tf_in=testfunction(tfact,outflow_boundaries,inflow_boundaries)
		
	inflow_rate=Vector{Float64}[]
	outflow_rate=Vector{Float64}[]
	reaction_rate=Vector{Float64}[]
	stored_amount=Vector{Float64}[]

	for i=2:length(solt)
		
		ifr=integrate(sys,tf_in,solt[i],solt[i-1],solt.t[i]-solt.t[i-1])
		ofr=integrate(sys,tf_out,solt[i],solt[i-1],solt.t[i]-solt.t[i-1])

		push!(inflow_rate,ifr)
		push!(outflow_rate,ofr)

		rr = vec(sum(integrate(sys,sys.physics.reaction,solt[i]), dims=2))
		amount = vec(sum(integrate(sys,sys.physics.storage,solt[i]), dims=2))

		push!(reaction_rate, rr)
		push!(stored_amount, amount)

   	end

	# integrals
	I_in=zeros(num_species(sys))
	I_out=zeros(num_species(sys))
	I_reac=zeros(num_species(sys))
	for i=1:length(solt)-1
		I_in .+= inflow_rate[i]*(solt.t[i+1]-solt.t[i])
		I_out .+= outflow_rate[i]*(solt.t[i+1]-solt.t[i])
		I_reac .+= reaction_rate[i]*(solt.t[i+1]-solt.t[i])
	end

	return (
		inflow_rate=inflow_rate,
		outflow_rate=outflow_rate,
		reaction_rate=reaction_rate, 
		stored_amount=stored_amount,
		I_in=I_in,
		I_out=I_out,
		I_reac=I_reac
	)

end

@doc raw"""
Helper function to calculate flux integrals over in- and outflow boundaries in 
    the Darcy-Maxwell-Stefan (DMS) model.
"""
function BoundaryFluxes(sol,sys,data)	
	(;inflow_boundaries,outflow_boundaries)=data
	tfact=TestFunctionFactory(sys)

	tf_in=testfunction(tfact,outflow_boundaries,inflow_boundaries)
	tf_out=testfunction(tfact,inflow_boundaries,outflow_boundaries)
	
	(;in=integrate(sys,tf_in,sol),out=integrate(sys,tf_out,sol) )
end

function BoundaryFluxes(t,solt,sys,data)	
	(;inflow_boundaries,outflow_boundaries)=data
	tfact=TestFunctionFactory(sys)

	tf_in=testfunction(tfact,outflow_boundaries,inflow_boundaries)
	tf_out=testfunction(tfact,inflow_boundaries,outflow_boundaries)

    ind = findfirst(x -> x.==t, solt.t)

    ifr=integrate(sys, tf_in, solt[ind], solt[ind-1], solt.t[ind]-solt.t[ind-1])
    ofr=integrate(sys, tf_out, solt[ind], solt[ind-1], solt.t[ind]-solt.t[ind-1])
	
	(;in=ifr, out=ofr)
end

@doc raw"""
Helper function to print a summary based on calculated flux integrals over in- 
    and outflow boundaries in the Darcy-Maxwell-Stefan (DMS) model.
"""
function Print_summary(sol,grid,sys,data)
    (;ip,m,mfluxin,mmix0,X0,ng,inflow_boundaries) = data
    in_,out_=BoundaryFluxes(sol,sys,data)

    nout(i) = out_[i]/m[i]
    local Ain = 0.0
	for boundary in inflow_boundaries
		Ain += bareas(boundary,sys,grid)
	end
    nin(i) = mfluxin/mmix0*Ain*X0[i]

    RI=sum(integrate(sys,sys.physics.reaction,sol),dims=2) # reaction integral

    println("Total mass inflows and outflows:")
    @printf "IN: %2.6e \t OUT: %2.6e \t REACT: %2.6e kg/hr \nSUM: %2.6e kg/hr\n\n" in_[ip]/ufac"kg/hr" out_[ip]/ufac"kg/hr" RI[ip]/ufac"kg/hr" (in_[ip]+out_[ip]+RI[ip])/ufac"kg/hr"

    println("Molar species inflows, outflows and reaction integrals:")
    for i = 1:ng
        @printf "%i\tIN: %2.6e \t OUT: %2.6e \t REACT: %2.6e mol/hr \n\tSUM: %2.6e mol/hr\n" i nin(i)/ufac"mol/hr" nout(i)/ufac"mol/hr" -RI[i]/m[i]/ufac"mol/hr" (nin(i)+nout(i)-RI[i]/m[i])/ufac"mol/hr"
    end
end

function Print_summary_ext_(in,out,HeatFluxes,data)
    (;gn,gni,m,nflowin,X0) = data
	ng=ngas(data)

	nout(i) = -out[i]/m[i]
	nin(i) = nflowin*X0[i]
	nout_dry = 0.0
	
	@printf "Molar species in- & outflows: (Norm conditions T=293.15 K, p=%2.5f bara)\n" 1.0ufac"atm"/ufac"bar"
	for i = 1:ng
		@printf "%s\tIN: %2.2f\t OUT: %2.2f mol/hr (%2.2f NL/min) \n" gn[i] nin(i)/ufac"mol/hr" nout(i)/ufac"mol/hr" nout(i)*ph"R"*293.15/(1.0ufac"atm")/ufac"l/minute"
		if i != gni[:H2O] 
			nout_dry += nout(i)
		end
	end


	println("\nDry Product Molar Fractions:")
	for i=1:ng
		if i != gni[:H2O] && i != gni[:N2] 
		@printf "%3s: %2.1f%%\n" gn[i] nout(i)/nout_dry*100
		end
	end
	
	println("\nConversion:")
	for i in(gni[:CO2],gni[:H2])
		@printf "X%4s: %2.2f\n" gn[i] (nin(i)-nout(i))/nin(i)
	end
	println("\nYield & Selectivity (CO2 based):")
	for i in(gni[:CO],gni[:CH4])
		@printf "Y%4s: %2.2f \tS%4s: %2.2f\n" gn[i] nout(i)/nin(gni[:CO2]) gn[i] nout(i)/(nin(gni[:CO2])-nout(gni[:CO2]))
	end


	# fluxes = HeatFluxes_EB_I(t,solt,grid,sys,data)

    if !isnothing(HeatFluxes)
        HeatFluxes_outer, HeatFluxes_inner = HeatFluxes
        
        # energy balance with outer radiation fluxes
        println("\nEnergy Balancing [W]:\nOuter radiation balance: (based on surface radiosities)\n")
        ns = keys(HeatFluxes_outer)
        vs = values(HeatFluxes_outer)
        for i=1:length(HeatFluxes_outer)
            @printf "%15s: %6.1f \n" String(ns[i]) vs[i]
        end
        @printf "%*s\n" 22 "="
        @printf "%15s: %6.1f \n" "Sum" sum(HeatFluxes_outer)

        # energy balance with inner radiation fluxes
        println("\nEnergy Balancing [W]:\nInner radiation balance: (based on net absorbed/emitted radiation fluxes)\n")
        ns = keys(HeatFluxes_inner)
        vs = values(HeatFluxes_inner)
        for i=1:length(HeatFluxes_inner)
            @printf "%21s: %6.1f \n" String(ns[i]) vs[i]
        end
        @printf "%*s\n" 28 "="
        @printf "%21s: %6.1f \n" "Sum" sum(HeatFluxes_inner)
    end
end

function Print_summary_ext(sol,sys,data)
    (;solve_T_equation) = data

	in_,out_ = BoundaryFluxes(sol,sys,data)
    if solve_T_equation
        HeatFluxes_outer = HeatFluxes_EB_I(sol,sys,data)
        HeatFluxes_inner = HeatFluxes_EB_I_inner(sol,sys,data)
        HeatFluxes = (HeatFluxes_outer, HeatFluxes_inner)
    else
        HeatFluxes = nothing
    end

	Print_summary_ext_(in_,out_,HeatFluxes,data)
end

@doc raw"""
Helper function to print an extended summary based on calculated flux integrals over in- 
    and outflow boundaries in the Darcy-Maxwell-Stefan (DMS) model.
"""
function Print_summary_ext(t,solt,sys,data)

    in_,out_ = BoundaryFluxes(t,solt,sys,data)	
    HeatFluxes = HeatFluxes_EB_I(t,solt,sys,data)

	Print_summary_ext_(in_,out_,HeatFluxes,data)
end

@doc raw"""
Helper function to export to VTK format for visualization 3D solutions of the
    photo thermal catalytic reactor (PTR) model. Exported are the pressure, 
    species molar fractions and temperature fields in the 3D domain.
"""
function WriteSolution3D(sol,grid,data;desc="")
    (;dim,ip,iT,ibf,gn,nflowin,solve_T_equation) = data
    ng=ngas(data)
    _t = now()
    tm = "$(hour(_t))_$(minute(_t))_$(second(_t))"
    desc = isempty(desc) ? desc : "_"*desc
    path = "../data/out/$(Date(_t))/$(tm)$(desc)"
    try
        mkpath(path)
    catch e
        println("Directory " * path * " already exists.")
    end
    #mkdir(string(Date(_t)))
    
    VoronoiFVM.writeVTK("$(path)/$(tm)_3D_ptot_$(data.nom_flux/ufac"kW/m^2")suns_$(nflowin/ufac"mol/hr").vtu", grid; point_data = sol[ip,:])
	if solve_T_equation
        VoronoiFVM.writeVTK("$(path)/$(tm)_3D_T_$(data.nom_flux/ufac"kW/m^2")suns_$(nflowin/ufac"mol/hr").vtu", grid; point_data = sol[iT,:] .-273.15)
    end
    for i=1:ng
        VoronoiFVM.writeVTK("$(path)/$(tm)_3D_x$(gn[i])_$(data.nom_flux/ufac"kW/m^2")suns_$(nflowin/ufac"mol/hr").vtu", grid; point_data = sol[i,:])
    end
	if dim == 3
		VoronoiFVM.writeVTK("$(path)/$(tm)_3D_irradiation_flux_$(data.nom_flux/ufac"kW/m^2")suns_$(nflowin/ufac"mol/hr").vtu", grid; point_data = sol[ibf,:])
	end
end

function WriteSolution3D(solt::TransientSolution,grid,data;desc="")
    (;dim,ip,iT,ibf,gn,nflowin,solve_T_equation) = data
    ng=ngas(data)
    _t = now()
    tm = "$(hour(_t))_$(minute(_t))_$(second(_t))"
    desc = isempty(desc) ? desc : "_"*desc
    path = "../data/out/$(Date(_t))/$(tm)$(desc)"
    try
        mkpath(path)
    catch e
        println("Directory " * path * " already exists.")
    end

    
    for (i,t) in enumerate(solt.t)
                
        # pressure
        VoronoiFVM.writeVTK("$(path)/3D_ptot_$(data.nom_flux/ufac"kW/m^2")suns_$(nflowin/ufac"mol/hr")_$(i).vtu", grid; point_data = solt(t)[ip,:])

        # temperature
        if solve_T_equation            
            VoronoiFVM.writeVTK("$(path)/3D_T_$(data.nom_flux/ufac"kW/m^2")suns_$(nflowin/ufac"mol/hr")_$(i).vtu", grid; point_data = solt(t)[iT,:] .-273.15)
        end

        # species molar fractions
        for j=1:ng
            VoronoiFVM.writeVTK("$(path)/3D_x$(gn[j])_$(data.nom_flux/ufac"kW/m^2")suns_$(nflowin/ufac"mol/hr")_$(i).vtu", grid; point_data = solt(t)[j,:])           
        end

        # irradiation boundary flux
        if dim == 3
            VoronoiFVM.writeVTK("$(path)/3D_irradiation_flux_$(data.nom_flux/ufac"kW/m^2")suns_$(nflowin/ufac"mol/hr")_$(i).vtu", grid; point_data = solt(t)[ibf,:])
        end

    end
    
end

function getIndices(grid)
    d1 = 0.025*sqrt(2)/2
    d2 = 0.05*sqrt(2)/2
    
    # coordinate origin in lower right corner
    #               ^ x
    #               |
    #               |
    #       <-------o
    #       y

    # uc
    T_03 = [0.08,0.08,0.005] # y,x,z
    T_04 = T_03 .+ [d2, -d2, 0]
    T_05 = T_03 .+ [-d1, d1, 0]
    T_06 = T_03 .+ [d1, -d1, 0]
    T_07 = T_03 .+ [-d2, -d2, 0]
    # lc
    T_12 = [0.08,0.08,0.0] # y,x,z
    T_13 = T_12 .+ [-d1, d1, 0]
    T_14 = T_12 .+ [-d2, -d2, 0]	
    
    
    function index(c, coords)
        cols = eachcol(coords)
        findmin(norm.([c - col for col in cols]))[2]
    end

    [index(c,grid[Coordinates]) for c in [T_03,T_04,T_05,T_06,T_07,T_12,T_13,T_14]]
    # [index(c,grid[Coordinates]) for c in [T_03]] # only get center temperature
end

function probe_Temps(solt,grid,data)

    sol = solt(solt.t[end])
    (;iT,dim) = data

    if dim == 2
        return sol[iT,91] - 273.15
    elseif dim == 3
        # return sol[iT,3035] .- 273.15
        return sol[iT,getIndices(grid)] .- 273.15
    end
    
end

function WriteTemperatures(solt,grid,data;desc="")

    (;nom_flux, nflowin) = data
    T_03,T_04,T_05,T_06,T_07,T_12,T_13,T_14 = probe_Temps(solt,grid,data)
              
    df = DataFrame()
    df[!, :nom_flux] = [nom_flux/ufac"kW/m^2"]
    df[!, :nflowin]  = [nflowin/ufac"mol/hr"]
    df[!, :T_03] .= T_03
    df[!, :T_04] .= T_04
    df[!, :T_05] .= T_05
    df[!, :T_06] .= T_06
    df[!, :T_07] .= T_07
    df[!, :T_03_Uc] .= 33.5 # average value of standard uncertainty in calculation of cat. surface T
    
    # df[!, :T_12] .= T_12
    # df[!, :T_13] .= T_13
    # df[!, :T_14] .= T_14

    _t = now()
    tm = "$(hour(_t))_$(minute(_t))_$(second(_t))"
    desc = isempty(desc) ? desc : "_"*desc
    path = "../data/out/$(Date(_t))/$(tm)$(desc)"
    try
        mkpath(path)
    catch e
        println("Directory " * path * " already exists.")
    end

    # CSV.write("data/out/2024-01-26/Tc_Uc.csv", df)
    CSV.write("$(path)/Sim_T_probe_$(data.nom_flux/ufac"kW/m^2")suns_$(nflowin/ufac"mol/hr").csv", df)
    # "$(path)/Sim_T_probe_$(data.nom_flux/ufac"kW/m^2")suns_$(nflowin/ufac"mol/hr").csv"
end