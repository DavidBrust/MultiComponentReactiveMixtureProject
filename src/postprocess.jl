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

function TotalFlows(t,solt,grid,sys,data)
    (;iT,Fluids,inlet_boundaries,outlet_boundaries,mfluxin,mmix0,nflowin, T_gas_in,X0)=data

    tfact=TestFunctionFactory(sys)	
	tf_out=testfunction(tfact,inlet_boundaries,outlet_boundaries)
	tf_in=testfunction(tfact,outlet_boundaries,inlet_boundaries)

    inflow_rate=Float64[]
	outflow_rate=Float64[]
    
    for i=2:length(solt)
		ifr=integrate(sys,tf_in,solt[i],solt[i-1],solt.t[i]-solt.t[i-1])[iT]
		ofr=integrate(sys,tf_out,solt[i],solt[i-1],solt.t[i]-solt.t[i-1])[iT]
		push!(inflow_rate,ifr)
		push!(outflow_rate,ofr)        
   	end

    ind = findfirst(x->x.== t, solt.t) -1 
    inflow_rate[ind], outflow_rate[ind]


end


function HeatFluxes_EB_I(t,solt,grid,sys,data)

    (;iT,inlet_boundaries,irradiated_boundaries,outlet_boundaries,side_boundaries,dt_hf_irrad)=data

    inf_t, outf_t = TotalFlows(t,solt,grid,sys,data)
    dE_dt = inf_t + outf_t
    # tend = solt.t[end]
    # sol = solt(tend)
    sol = solt(t)


    # calc_hf = tend >= dt_hf_irrad[1]
    calc_hf = t >= dt_hf_irrad[1]

    Qconv_10=integrate(sys,flux_convection_top,sol; boundary=true)[iT,irradiated_boundaries]
    Qconv_10 = calc_hf ? sum(Qconv_10) : 0.0
    
    QG_01=integrate(sys,flux_window_underside,sol; boundary=true)[iT,irradiated_boundaries]
    QG_01 = calc_hf ? sum(QG_01) : 0.0

    QG_10=integrate(sys,flux_catalyst_layer,sol; boundary=true)[iT,irradiated_boundaries]
    QG_10 = calc_hf ? sum(QG_10) : 0.0

    QG_34=integrate(sys,flux_radiation_frit_bottom,sol; boundary=true)[iT,outlet_boundaries]
    QG_34 = calc_hf ? sum(QG_34) : 0.0
    
    QG_43=integrate(sys,flux_radiation_plate_bottom,sol; boundary=true)[iT,outlet_boundaries]
    QG_43 = calc_hf ? sum(QG_43) : 0.0

    Qconv_34=integrate(sys,flux_convection_bottom,sol; boundary=true)[iT,outlet_boundaries]
    Qconv_34 = calc_hf ? sum(Qconv_34) : 0.0

    Qsides=integrate(sys,FixedBed.PCR_side,sol; boundary=true)[iT,side_boundaries]
    Qsides = calc_hf ? sum(Qsides) : 0.0

    dH_dot = dE_dt - (QG_01 - QG_10) - (QG_43 - QG_34) + Qconv_10 + Qconv_34 + Qsides

    (
        QG_01=QG_01,
        QG_10=-QG_10,
        dH_dot=dH_dot,
        Qconv_10=-Qconv_10,
        Qconv_34=-Qconv_34,
        QG_34=-QG_34,
        QG_43=QG_43,
        Qsides=-Qsides
    )
    

end