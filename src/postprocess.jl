# G0_bot, IR/vis : surface radiosities of glass underside
function flux_window_underside(f,u,bnode,data)
    (;iT,irradiated_boundaries)=data
    # if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit
    if bnode.region in irradiated_boundaries        

        G1_bot_vis, G1_bot_IR = radiosity_window(f,u,bnode,data)

        f[iT] = G1_bot_vis + G1_bot_IR
    end
end


# # G2_IR/vis : surface radiosities of catalyst layer 
function flux_catalyst_layer(f,u,bnode,data)
    (;iT,uc_cat,irradiated_boundaries)=data
    if bnode.region in irradiated_boundaries
        
        # catalyst layer properties (2)
        rho2_vis=uc_cat.rho_vis
        rho2_IR=uc_cat.rho_IR
        eps2=uc_cat.eps

        G1_bot_vis, G1_bot_IR = radiosity_window(f,u,bnode,data)

        G2_vis = rho2_vis*G1_bot_vis

        G2_IR = eps2*ph"σ"*u[iT]^4 + rho2_IR*G1_bot_IR

        f[iT] = G2_vis + G2_IR
    end
end

function flux_convection_top(f,u,bnode,data)
    (;iT,iTw,X0,uc_h,Nu,irradiated_boundaries)=data
    if bnode.region in irradiated_boundaries
   
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

function flux_radiation_frit_bottom(f,u,bnode,data)
    (;iT,iTp,lc_frit,lc_plate,outlet_boundaries) = data
    #if bnode.region==Γ_bottom
    if bnode.region in outlet_boundaries
		
		
		# irradiation exchange between porous frit (1) and Al bottom plate (2)
		# porous frit properties (1)
		eps1=lc_frit.eps;  rho1_IR=lc_frit.rho_IR; 	
		# Al bottom plate properties (2)
		eps2=lc_plate.eps; rho2_IR=lc_plate.rho_IR;
	
        # local plate temperature
        Tplate = u[iTp]
        G1_IR = (eps1*ph"σ"*u[iT]^4 + rho1_IR*eps2*ph"σ"*Tplate^4)/(1-rho1_IR*rho2_IR)

        f[iT] = G1_IR
    end
end

function flux_radiation_plate_bottom(f,u,bnode,data)
    (;iT,iTp,lc_frit,lc_plate,outlet_boundaries) = data
    # if bnode.region==Γ_bottom
    if bnode.region in outlet_boundaries		
		
		# irradiation exchange between porous frit (1) and Al bottom plate (2)
		# porous frit properties (1)
		eps1=lc_frit.eps;  rho1_IR=lc_frit.rho_IR; 	
		# Al bottom plate properties (2)
		eps2=lc_plate.eps; rho2_IR=lc_plate.rho_IR;
		
        # local plate temperature
        Tplate = u[iTp]
        G2_IR = (eps2*ph"σ"*Tplate^4 + rho2_IR*eps1*ph"σ"*u[iT]^4)/(1-rho1_IR*rho2_IR)

        f[iT] = G2_IR
    end
end

function flux_convection_bottom(f,u,bnode,data)
    (;iT,iTp,lc_h,Nu,outlet_boundaries)=data

    # if bnode.region==Γ_bottom        
    if bnode.region in outlet_boundaries
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

function TotalFlows(t,solt,grid,sys,data)
    (;iT,Fluids,inlet_boundaries,outlet_boundaries)=data

    tfact=TestFunctionFactory(sys)	
	tf_out=testfunction(tfact,inlet_boundaries,outlet_boundaries)
	tf_in=testfunction(tfact,outlet_boundaries,inlet_boundaries)

    inflow_rate=Float64[]
	outflow_rate=Float64[]
    
    # for i=2:length(solt)
	# 	ifr=integrate(sys,tf_in,solt[i],solt[i-1],solt.t[i]-solt.t[i-1])[iT]
	# 	ofr=integrate(sys,tf_out,solt[i],solt[i-1],solt.t[i]-solt.t[i-1])[iT]
	# 	push!(inflow_rate,ifr)
	# 	push!(outflow_rate,ofr)        
   	# end

    # ind = findfirst(x->x.== t, solt.t) -1
    ind = findfirst(x->x.== t, solt.t)


    ifr=integrate(sys,tf_in,solt[ind],solt[ind-1],solt.t[ind]-solt.t[ind-1])[iT]
    ofr=integrate(sys,tf_out,solt[ind],solt[ind-1],solt.t[ind]-solt.t[ind-1])[iT]
   
    ifr, ofr
    # inflow_rate[ind], outflow_rate[ind]
end


function HeatFluxes_EB_I(t,solt,grid,sys,data)

    (;iT,inlet_boundaries,irradiated_boundaries,outlet_boundaries,side_boundaries,dt_hf_irrad)=data

    inflow, outflow = TotalFlows(t,solt,grid,sys,data)
    dE_dt = inflow + outflow
    # tend = solt.t[end]
    # sol = solt(tend)
    sol = solt(t)


    # calc_hf = tend >= dt_hf_irrad[1]
    calc_hf = t >= dt_hf_irrad[1]

    Q_conve_10=integrate(sys,flux_convection_top,sol; boundary=true)[iT,irradiated_boundaries]
    Q_conve_10 = calc_hf ? sum(Q_conve_10) : 0.0
    
    Q_irrad_01=integrate(sys,flux_window_underside,sol; boundary=true)[iT,irradiated_boundaries]
    Q_irrad_01 = calc_hf ? sum(Q_irrad_01) : 0.0

    Q_irrad_10=integrate(sys,flux_catalyst_layer,sol; boundary=true)[iT,irradiated_boundaries]
    Q_irrad_10 = calc_hf ? sum(Q_irrad_10) : 0.0

    Q_irrad_34=integrate(sys,flux_radiation_frit_bottom,sol; boundary=true)[iT,outlet_boundaries]
    Q_irrad_34 = calc_hf ? sum(Q_irrad_34) : 0.0
    
    Q_irrad_43=integrate(sys,flux_radiation_plate_bottom,sol; boundary=true)[iT,outlet_boundaries]
    Q_irrad_43 = calc_hf ? sum(Q_irrad_43) : 0.0

    Q_conve_34=integrate(sys,flux_convection_bottom,sol; boundary=true)[iT,outlet_boundaries]
    Q_conve_34 = calc_hf ? sum(Q_conve_34) : 0.0

    # Qsides=integrate(sys,MultiComponentReactiveMixtureProject.PCR_side,sol; boundary=true)[iT,side_boundaries]
    Q_sides=integrate(sys,flux_side,sol; boundary=true)[iT,side_boundaries]
    Q_sides = sum(Q_sides)

    H_reaction = sum(integrate(sys,sys.physics.reaction,sol), dims=2)[iT]

    H_thermal = dE_dt - (Q_irrad_01 - Q_irrad_10) - (Q_irrad_43 - Q_irrad_34) + Q_conve_10 + Q_conve_34 + Q_sides

    (
        Q_irrad_01=Q_irrad_01,
        Q_irrad_10=-Q_irrad_10,
        H_thermal=H_thermal,
        H_reaction=-H_reaction,
        Q_conve_10=-Q_conve_10,
        Q_conve_34=-Q_conve_34,
        Q_irrad_34=-Q_irrad_34,
        Q_irrad_43=Q_irrad_43,
        Q_sides=-Q_sides
    )
    

end

function BoundaryFlows_Integrals(solt, sys, data)
	(;outlet_boundaries,inlet_boundaries)=data
	
	tfact=TestFunctionFactory(sys)	
	tf_out=testfunction(tfact,inlet_boundaries,outlet_boundaries)
	tf_in=testfunction(tfact,outlet_boundaries,inlet_boundaries)
		
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