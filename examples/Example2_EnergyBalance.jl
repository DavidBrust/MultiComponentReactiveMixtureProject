module Example2_EnergyBalance

using DataFrames, CSV


include("../notebooks/PorousCatalystHot3DTopFlowIrrExchange.jl")

function RunSim()
    data=ModelData(
        isreactive=0
        )
    sol_embed,grid,sys,data_embed = main(data=data);
    sol = sol_embed(1.0)
    sol,grid,sys,data_embed
end;

function flux_in_aperture(f,u,bnode,data)
    if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit
        iT=data.iT
        f[iT] = data.Glamp
    end
end

function flux_radiation_top(f,u,bnode,data,T2,T3)
    if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit

        (;ng,iT,Glamp,Tglass,uc_window,uc_cat,uc_frit,vf_uc_window_cat,vf_uc_window_frit)=data
        σ=ph"σ"

        # irradiation exchange between quartz window (1), cat surface (2) & frit surface (3)
        # window properties (1)
        tau1_vis=uc_window.tau_vis
        rho1_vis=uc_window.rho_vis
        tau1_IR=uc_window.tau_IR
        rho1_IR=uc_window.rho_IR
        eps1=uc_window.eps
        # catalyst layer properties (2)
        rho2_vis=uc_cat.rho_vis
        rho2_IR=uc_cat.rho_IR
        eps2=uc_cat.eps
        # uncoated frit properties (3)
        rho3_vis=uc_frit.rho_vis
        rho3_IR=uc_frit.rho_IR
        eps3=uc_frit.eps
        #view factors
        ϕ12=vf_uc_window_cat
        ϕ13=vf_uc_window_frit


        # surface brigthness of quartz window inwards / towards catalyst facing surface (1) in vis & IR	
        G1_bot_vis = tau1_vis*Glamp/(1-rho1_vis*(ϕ12*rho2_vis+ϕ13*rho3_vis))

		#G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*u[iT]^4+ϕ13*eps3*σ*u[iT]^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))
        G1_bot_IR = (eps1*σ*Tglass^4 + rho1_IR*(ϕ12*eps2*σ*T2^4+ϕ13*eps3*σ*T3^4))/(1-rho1_IR*(ϕ12*rho2_IR+ϕ13*rho3_IR))

                    

        G2_vis = rho2_vis*G1_bot_vis
        #G2_IR = eps2*σ*u[iT]^4 + rho2_IR*G1_bot_IR
        G2_IR = eps2*σ*T2^4 + rho2_IR*G1_bot_IR

        G3_vis = rho3_vis*G1_bot_vis
        #G3_IR = eps3*σ*u[iT]^4 + rho3_IR*G1_bot_IR
        G3_IR = eps3*σ*T3^4 + rho3_IR*G1_bot_IR

        # surface brigthness of quartz window outward facing surface (1) in vis & IR	
        G1_top_vis = rho1_vis*Glamp + tau1_vis*(ϕ12*G2_vis + ϕ13*G3_vis)
        G1_top_IR = tau1_IR*(ϕ12*G2_IR + ϕ13*G3_IR) + eps1*σ*Tglass^4

    
        iT=data.iT
        f[iT] = G1_top_vis + G1_top_IR
    end
end


# conductive heat flux through top chamber
function flux_conduction_top(f,u,bnode,data)
    if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit
        (;iT,Tglass,X0,uc_h,vf_uc_window_cat,vf_uc_window_frit)=data

		# mean temperature
        Tm=0.5*(u[iT] + Tglass)

        # thermal conductivity at Tm and inlet composition X0 (H2/CO2 = 1/1)
        _,λf=dynvisc_thermcond_mix(data, Tm, X0)

        # positive flux in positive z-coordinate
        q_cond = -λf*(Tglass-u[iT])/uc_h

        f[iT] = q_cond

    end
end

# conductive heat flux through top chamber
function flux_conduction_top(f,u,bnode,data,T2,T3)
    if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit
        (;iT,Tglass,X0,uc_h,vf_uc_window_cat,vf_uc_window_frit)=data

        #view factors
	    ϕ12=vf_uc_window_cat
	    ϕ13=vf_uc_window_frit

		# mean temperature
        Tm2=0.5*(T2 + Tglass)
		Tm3=0.5*(T3 + Tglass)

        # thermal conductivity at Tm and inlet composition X0 (H2/CO2 = 1/1)
        _,λf2=dynvisc_thermcond_mix(data, Tm2, X0)
		_,λf3=dynvisc_thermcond_mix(data, Tm3, X0)

        # positive flux in positive z-coordinate
        q_cond = -λf2*(Tglass-T2)/uc_h*ϕ12 -λf3*(Tglass-T3)/uc_h*ϕ13

        f[iT] = q_cond

    end
end


function flux_enthalpy_bottom(f,u,bnode,data)
    if bnode.region==Γ_bottom
		(;ng,iT,ip,Fluids,ubot,Tamb) = data
		
        X=zeros(eltype(u), ng)
        mole_frac!(bnode,data,X,u)

        hmix=enthalpy_mix(Fluids, u[iT], X)
		flux_enth_bot=ubot*u[ip]/(ph"R"*u[iT])*hmix

		#flux_enth_bot=ubot*u[ip]/(ph"R"*u[iT])*cf*(u[iT]-Tamb)
        f[iT]=flux_enth_bot
    end
end

function flux_enthalpy_top(f,u,bnode,data)
    if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit
        
		(;iT,Fluids,u0,X0,pn,Tn,Tamb) = data
		
        # for i=1:ng
        #     # sign convention: outward pointing surface normal as positive
        #     MolarFlow = -u0*X0[i]*pn/(ph"R"*Tn)
        #     f[iT] += MolarFlow * enthalpy_gas(Fluids[i], Tamb)
        # end

        hmix=enthalpy_mix(Fluids, Tamb, X0)
		flux_enth_top=u0*pn/(ph"R"*Tn)*hmix

        f[iT]=flux_enth_top

    end
end


function flux_enthalpy_reaction(f,u,node,data)
    if node.region == 2 && data.isreactive # catalyst layer
		(;ng,iT,Fluids,kinpar,mcats) = data
		(;nuij) = kinpar
		
		pi = u[1:ng]./ufac"bar"
		# negative sign: sign convention of VoronoiFVM: source term < 0 
		# ri returns reaction rate in mol/(h gcat)
		RR = -mcats*ri(kinpar,u[iT],pi)*ufac"mol/hr"*ufac"1/g"
		## Xu & Froment 1989 kinetics
		# R1: CH4 + H2O = CO + 3 H2
		# R2: CO + H2O = CO2 + H2
		# R3: CH4 + 2 H2O = CO2 + 4 H2

		for i=1:ng
			f[i] = sum(nuij[i,:] .* RR)
		end
		
		
		ΔHiT = zeros(size(nuij,2))
		for i=1:ng
			ΔHiT .-= enthalpy_gas(Fluids[i],u[iT])*nuij[i,:]
		end

		# temperature eq. / heat source
		f[iT] = sum(RR .* ΔHiT)
	end
end

# formulation for net radiative transport in bottom chamber, energy balance around
# underside of frit
function flux_radiation_bottom(f,u,bnode,data)
    if bnode.region==Γ_bottom
		(;iT,lc_frit,lc_plate,Tplate) = data
		
		# irradiation exchange between porous frit (1) and Al bottom plate (2)
		# porous frit properties (1)
		eps1=lc_frit.eps; alpha1_IR=lc_frit.alpha_IR; rho1_IR=lc_frit.rho_IR; 	
		# Al bottom plate properties (2)
		eps2=lc_plate.eps; rho2_IR=lc_plate.rho_IR;
		
		#(;eps1,alpha1_IR,rho1_IR,rho2_IR,eps2,rho2_IR) = lower_chamber
		σ=ph"σ"
		flux_irrad = -eps1*σ*u[iT]^4 + alpha1_IR/(1-rho1_IR*rho2_IR)*(eps2*σ*Tplate^4+rho2_IR*eps1*σ*u[iT]^4)
        f[iT] = -flux_irrad
    end
end

function flux_conduction_bottom(f,u,bnode,data)
    if bnode.region==Γ_bottom
        (;ng,iT,Tplate,lc_h)=data
        X=zeros(eltype(u), ng)    	

        # height of bottom chamber/ distance over which heat conduction occurs
        Tsurf = u[iT]

        # mean temperature
        Tm=0.5*(u[iT] + Tplate)
        mole_frac!(bnode,data,X,u)
        
        # thermal conductivity at Tm and outlet composition
        _,λf=dynvisc_thermcond_mix(data, Tm, X)

        f[iT] = λf*(Tsurf-Tplate)/lc_h

    end
end

function sidewalls(f,u,bnode,data)
    (;iT,α_w,Tamb) = data
    boundary_robin!(f,u,bnode;species=iT,region=Γ_side_right, factor=α_w, value=Tamb*α_w)
    boundary_robin!(f,u,bnode;species=iT,region=Γ_side_back, factor=α_w, value=Tamb*α_w)
end

function areas(sol,sys,grid,data)
	iT = data.iT
	function area(f,u,bnode,data)
		# repurpose temperature index to hold area information
		f[iT] = one(eltype(u))
	end

	integrate(sys,area,sol; boundary=true)[iT,:]
end

function T_avg(sol,sys,grid,data)

	iT = data.iT	
	function T_avg_(f,u,bnode,data)			
		f[iT] = u[iT]		
	end

	areas_=areas(sol,sys,grid,data)
	
	T_int=integrate(sys,T_avg_,sol; boundary=true)[iT,:]
	T_int./areas_	
end


function HeatFluxes()

    data=ModelData()
    
    sol_embed,grid,sys,data_embed = main(data=data);
    sol = sol_embed(1.0) # select final solution at the end of parameter embedding

    (;iT,lc_plate,Tglass,Tplate,Tamb,α_nat_conv)=data_embed
    σ=ph"σ"

    Tfrit_avg,Tcat_avg = T_avg(sol,sys,grid,data_embed)[[Γ_top_frit,Γ_top_cat]]
    
    q_in = 4*integrate(sys,flux_in_aperture,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    q_in = sum(q_in)

    q_rad_top = 4*integrate(
        sys,
        (f,u,bnode,data) -> flux_radiation_top(f,u,bnode,data_embed,Tcat_avg,Tfrit_avg),sol; boundary=true
        )[iT,[Γ_top_cat,Γ_top_frit]]
    q_rad_top = sum(q_rad_top)

    q_cond_top = 4*integrate(sys,flux_conduction_top,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    q_cond_top = sum(q_cond_top)

    q_sidewalls = 4*integrate(sys,sidewalls,sol; boundary=true)[iT,[Γ_side_right,Γ_side_back]]        
    q_sidewalls = sum(q_sidewalls)

    q_rad_bottom = 4*integrate(sys,flux_radiation_bottom,sol; boundary=true)[iT,Γ_bottom]
    q_cond_bottom = 4*integrate(sys,flux_conduction_bottom,sol; boundary=true)[iT,Γ_bottom]

    # alternative way to calculate heat transfer over system boundaries via
    # outer energy balance
    # bottom boundary
    areas_=areas(sol,sys,grid,data_embed)
    q_rad_bottom_2 = 4*areas_[Γ_bottom]*lc_plate.eps*σ*Tplate^4
    q_conv_bottom = 4*areas_[Γ_bottom]*α_nat_conv*(Tplate-Tamb)
    # top boundary
    q_conv_top = 4*sum(areas_[[Γ_top_cat,Γ_top_frit]])*α_nat_conv*(Tglass-Tamb)

    q_enthalpy_top_in = 4*integrate(sys,flux_enthalpy_top,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    q_enthalpy_top_in = sum(q_enthalpy_top_in)

    
    q_enthalpy_bottom_out = 4*integrate(sys,flux_enthalpy_bottom,sol; boundary=true)[iT,Γ_bottom]
    
    q_delta_enthalpy = q_enthalpy_bottom_out - q_enthalpy_top_in
    
    q_delta_enthalpy_reaction = 4*integrate(sys,flux_enthalpy_reaction,sol)[iT,2]



    (
        q_in=q_in,
        q_rad_top=q_rad_top,
        q_cond_top=q_cond_top,
        q_sidewalls=q_sidewalls,
        q_rad_bottom=q_rad_bottom,
        q_cond_bottom=q_cond_bottom,
        q_delta_enthalpy=q_delta_enthalpy,
        q_delta_enthalpy_reaction=q_delta_enthalpy_reaction,

        q_rad_bottom_2=q_rad_bottom_2,
        q_conv_bottom=q_conv_bottom,

        q_conv_top=q_conv_top
    )
end

function loss_plot(losses)
    (;q_in,q_rad_top,q_cond_top,q_sidewalls,q_rad_bottom,q_cond_bottom,q_delta_enthalpy,q_delta_enthalpy_reaction,q_rad_bottom_2,q_conv_bottom,q_conv_top) = losses

    q_delta_enthalpy_sensible = q_delta_enthalpy-q_delta_enthalpy_reaction

    lbs = ["Top Rad", "Top Cond", "Conv Sides", "Bottom Rad", "Bottom Cond", "HSensible","HReaction"]
    #lbs = ["Top Rad", "Top Cond", "Top Conv", "Conv Sides", "Bottom Rad", "Bottom Cond", "Bottom Rad 2", "Bottom Conv", "HSensible", "HReaction"]
    x = [string(x) for x in 1:length(lbs)]
    lbsx = [x*": "* l for (l,x) in zip(lbs,x)]
    total_loss = sum([q_rad_top,q_cond_top,q_sidewalls,q_rad_bottom,q_cond_bottom,q_delta_enthalpy])
    y = [q_rad_top,q_cond_top,q_sidewalls,q_rad_bottom,q_cond_bottom,q_delta_enthalpy_sensible,q_delta_enthalpy_reaction] / q_in * 100.0
    #y = [q_rad_top,q_cond_top,q_conv_top,q_sidewalls,q_rad_bottom,q_cond_bottom,q_rad_bottom_2,q_conv_bottom,q_delta_enthalpy_sensible,q_delta_enthalpy_reaction] / q_in * 100.0
    
    p=Plots.plot(size=(400,300), yguide="Heat flow contribution / %")
    Plots.plot!(p, permutedims(x), permutedims(y); st=:bar, label=permutedims(lbsx))    
    Plots.annotate!(p, x, y, round.(y, sigdigits=2), :bottom)
    Plots.annotate!(p, 2.0, 30.0, "Σ flows = "*string(round(total_loss/q_in*100.0,sigdigits=3))*" %", :left)
    

end



end
