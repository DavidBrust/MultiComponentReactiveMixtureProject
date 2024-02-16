module EnergyBalance

using DataFrames, CSV, Revise, LessUnitful


#include("../notebooks/PTReactorBase.jl")
# include("../notebooks/Demo/ThermalDemo.jl")



# function flux_in_profile(f,u,bnode,data)
#     if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit
#         (;iT,FluxIntp)=data
        
#         # @views x,y,z = bnode.coord[:,bnode.index]
#         # different ordering of axes in VoronoiFVM than in interpolated data
#         @views y,x,_ = bnode.coord[:,bnode.index] 
#         f[iT] = FluxIntp(x,y)
#     end
# end


function radiosity_window(f,u,bnode,data)
    (;iTw,G_lamp,Tamb,uc_window,irradiated_boundaries)=data

    eps1=uc_window.eps

    Tglass = u[iTw] # local tempererature of quartz window
    G1_bot_IR = eps1*ph"σ"*(Tglass^4-Tamb^4)+ph"σ"*Tamb^4
	G1_bot_vis = 0.0
    if bnode.region in irradiated_boundaries		
		G1_bot_vis += G_lamp # flux profile measured behind quarz in plane of cat layer
    end
    return G1_bot_vis,G1_bot_IR
end

# G0_bot, IR/vis : surface radiosities of glass underside
function flux_window_underside(f,u,bnode,data)
    (;iT,irradiated_boundaries)=data
    # if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit
    if bnode.region in irradiated_boundaries        

        G1_bot_vis, G1_bot_IR = radiosity_window(f,u,bnode,data)

        f[iT] = G1_bot_vis + G1_bot_IR
    end
end

# function irrad_top(f,u,bnode,data)
#     (;iT,irradiated_boundaries,uc_cat)=data

#     if bnode.region in irradiated_boundaries
#         alpha2_vis=uc_cat.alpha_vis
#         alpha2_IR=uc_cat.alpha_IR
#         eps2=uc_cat.eps       

#         G1_vis, G1_IR = radiosity_window(f,u,bnode,data)

#         hflux_irrad = (-eps2*ph"σ"*u[iT]^4 + alpha2_vis*G1_vis + alpha2_IR*G1_IR) 
#         f[iT] += hflux_irrad
#     end
# end

# function irrad_bottom(f,u,bnode,data)
#     (;iT,iTp,outlet_boundaries,lc_frit,lc_plate)=data

#     if bnode.region in outlet_boundaries
#         rho1_IR=lc_frit.rho_IR
# 		alpha1_IR=lc_frit.alpha_IR
# 		eps1=lc_frit.eps
# 		rho2_IR=lc_plate.rho_IR
# 		eps2=lc_plate.eps      

#         Tplate = u[iTp]
        
#         hflux_irrad = -eps1*ph"σ"*u[iT]^4 + alpha1_IR/(1-rho1_IR*rho2_IR)*(eps2*ph"σ"*Tplate^4+rho2_IR*eps1*ph"σ"*u[iT]^4)
#         f[iT] += hflux_irrad
#     end
# end


# G0_bot, IR only
# function flux_window_underside_IR(f,u,bnode,data)
#     if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit
#         (;iT)=data

#         _, G1_bot_IR = radiosity_window(f,u,bnode,data)

#         f[iT] = G1_bot_IR
#     end
# end

# # G0_bot, vis only
# function flux_window_underside_vis(f,u,bnode,data)
#     if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit
#         (;iT)=data

#         G1_bot_vis, _ = radiosity_window(f,u,bnode,data)

#         f[iT] = G1_bot_vis
#     end
# end

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

# # G2 IR only
# function flux_catalyst_layer_IR(f,u,bnode,data)
#     if bnode.region==Γ_top_cat
#         (;iT,uc_cat)=data
#         # catalyst layer properties (2)
#         rho2_vis=uc_cat.rho_vis
#         rho2_IR=uc_cat.rho_IR
#         eps2=uc_cat.eps

#         G1_bot_vis, G1_bot_IR = radiosity_window(f,u,bnode,data)
        
#         G2_IR = eps2*ph"σ"*u[iT]^4 + rho2_IR*G1_bot_IR

#         f[iT] = G2_IR
#     end
# end

# # G2 vis only
# function flux_catalyst_layer_vis(f,u,bnode,data)
#     if bnode.region==Γ_top_cat
#         (;iT,uc_cat)=data
#         # catalyst layer properties (2)
#         rho2_vis=uc_cat.rho_vis
#         rho2_IR=uc_cat.rho_IR
#         eps2=uc_cat.eps

#         G1_bot_vis, G1_bot_IR = radiosity_window(f,u,bnode,data)
        
#         G2_vis = rho2_vis*G1_bot_vis

#         f[iT] = G2_vis
#     end
# end

# # G3_IR/vis : surface radiosities of frit 
# function flux_frit(f,u,bnode,data)
#     if bnode.region==Γ_top_frit        
#         (;iT,uc_frit)=data
#         # uncoated frit properties (3)
#         rho3_vis=uc_frit.rho_vis
#         rho3_IR=uc_frit.rho_IR
#         eps3=uc_frit.eps

#         G1_bot_vis, G1_bot_IR = radiosity_window(f,u,bnode,data)

#         G3_vis = rho3_vis*G1_bot_vis
        
#         G3_IR = eps3*ph"σ"*u[iT]^4 + rho3_IR*G1_bot_IR
        
#         f[iT] = G3_vis + G3_IR
#     end
# end

# # G3 IR only
# function flux_frit_IR(f,u,bnode,data)
#     if bnode.region==Γ_top_frit        
#         (;iT,uc_frit)=data
#         # uncoated frit properties (3)
#         rho3_vis=uc_frit.rho_vis
#         rho3_IR=uc_frit.rho_IR
#         eps3=uc_frit.eps

#         G1_bot_vis, G1_bot_IR = radiosity_window(f,u,bnode,data)
       
#         G3_IR = eps3*ph"σ"*u[iT]^4 + rho3_IR*G1_bot_IR
        
#         f[iT] =  G3_IR
#     end
# end

# # G3 vis only
# function flux_frit_vis(f,u,bnode,data)
#     if bnode.region==Γ_top_frit        
#         (;iT,uc_frit)=data
#         # uncoated frit properties (3)
#         rho3_vis=uc_frit.rho_vis
#         rho3_IR=uc_frit.rho_IR
#         eps3=uc_frit.eps

#         G1_bot_vis, G1_bot_IR = radiosity_window(f,u,bnode,data)

#         G3_vis = rho3_vis*G1_bot_vis
                
#         f[iT] = G3_vis
#     end
# end

# # tau0*G1_IR/vis : transmitted fraction through window of surface radiosityies of catalyst layer 
# function flux_transmission_top_catalyst_layer(f,u,bnode,data)
#     if bnode.region==Γ_top_cat
#         (;iT,uc_cat,uc_window)=data
#         # window properties (1)
#         tau1_vis=uc_window.tau_vis
#         tau1_IR=uc_window.tau_IR
#         # catalyst layer properties (2)
#         rho2_vis=uc_cat.rho_vis
#         rho2_IR=uc_cat.rho_IR
#         eps2=uc_cat.eps

#         G1_bot_vis, G1_bot_IR = radiosity_window(f,u,bnode,data)

#         G2_vis = rho2_vis*G1_bot_vis

#         G2_IR = eps2*ph"σ"*u[iT]^4 + rho2_IR*G1_bot_IR

#         f[iT] = tau1_vis*G2_vis + tau1_IR*G2_IR
#     end
# end

# # tau0*G2_IR/vis : transmitted fraction through window of surface radiosityies of frit
# function flux_transmission_top_frit(f,u,bnode,data)
#     if bnode.region==Γ_top_frit
#         (;iT,uc_frit,uc_window)=data

#         # window properties (1)
#         tau1_vis=uc_window.tau_vis
#         tau1_IR=uc_window.tau_IR
#         # uncoated frit properties (3)
#         rho3_vis=uc_frit.rho_vis
#         rho3_IR=uc_frit.rho_IR
#         eps3=uc_frit.eps

#         G1_bot_vis, G1_bot_IR = radiosity_window(f,u,bnode,data)
        
#         G3_vis = rho3_vis*G1_bot_vis

#         G3_IR = eps3*ph"σ"*u[iT]^4 + rho3_IR*G1_bot_IR

#         f[iT] = tau1_vis*G3_vis + tau1_IR*G3_IR
#     end
# end

# # convective heat flux through top chamber
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

# # convective heat flux exiting upper window surface
# function flux_convection_window_top(f,u,bnode,data)
#     if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit
#         (;iT,iTw,Tamb,k_nat_conv)=data
#         # local quartz window temperature
#         Tglass = u[iTw]

#         f[iT] = k_nat_conv*(Tglass-Tamb)
#     end
# end

# # emission heat flux exiting upper window surface
# function flux_radiation_window_top(f,u,bnode,data)
#     if bnode.region==Γ_top_cat || bnode.region==Γ_top_frit
#         (;iT,uc_window,iTw)=data
#         # local quartz window temperature
#         Tglass = u[iTw]
        
#         f[iT] = uc_window.eps*ph"σ"*Tglass^4
#     end
# end

# # surface radiosity of underside of frit, facing towards the bottom plate
# # solely in IR, consists of emission and reflection
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

# # conductive heat flux through bottom chamber
# function flux_conduction_bottom(f,u,bnode,data)
#     if bnode.region==Γ_bottom
#         (;iT,Tplate,lc_h)=data
#         ng=ngas(data)

#         X=zeros(eltype(u), ng)
# 		mole_frac!(bnode,data,X,u)

# 		# mean temperature
#         Tm=0.5*(u[iT] + Tplate)
#         # thermal conductivity at Tm and outlet composition X
#         _,λf=dynvisc_thermcond_mix(data, Tm, X)

#         # positive flux in negative z coord. -> pointing towards bottom plate
#         q_cond = -λf*(Tplate-u[iT])/lc_h 

#         f[iT] = q_cond

#     end
# end

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

# # convective heat flux exiting upper window surface
# function flux_convection_outer_plate(f,u,bnode,data)
#     if bnode.region==Γ_bottom
#         (;iT,iTp,Tamb,k_nat_conv)=data
#         # local Al bottom plate temperature
#         Tplate = u[iTp]

#         f[iT] = k_nat_conv*(Tplate-Tamb)
#     end
# end

# # emission heat flux exiting upper window surface
# function flux_emission_outer_plate(f,u,bnode,data)
#     if bnode.region==Γ_bottom
#         (;iT,lc_plate,iTp)=data
#         # local Al bottom plate temperature
#         Tplate = u[iTp]
        
#         f[iT] = lc_plate.eps*ph"σ"*Tplate^4
#     end
# end

# function flux_abs_top_catalyst_layer(f,u,bnode,data)
#     if bnode.region==Γ_top_cat
#         (;iT,uc_cat,uc_window)=data

#         # window properties (1)
#         alpha1_vis=uc_window.alpha_vis
#         alpha1_IR=uc_window.alpha_IR
#         # catalyst layer properties (2)
#         rho2_vis=uc_cat.rho_vis
#         rho2_IR=uc_cat.rho_IR
#         eps2=uc_cat.eps

#         G1_bot_vis, G1_bot_IR = radiosity_window(f,u,bnode,data)

#         G2_vis = rho2_vis*G1_bot_vis		
        
#         G2_IR = eps2*ph"σ"*u[iT]^4 + rho2_IR*G1_bot_IR

#         f[iT] = alpha1_vis*G2_vis + alpha1_IR*G2_IR
#     end
# end

# function flux_abs_top_frit(f,u,bnode,data)
#     if bnode.region==Γ_top_frit

#         (;iT,uc_frit,uc_window)=data

#         # window properties (1)
#         alpha1_vis=uc_window.alpha_vis
#         alpha1_IR=uc_window.alpha_IR
#         # uncoated frit properties (3)
#         rho3_vis=uc_frit.rho_vis
#         rho3_IR=uc_frit.rho_IR
#         eps3=uc_frit.eps

#         G1_bot_vis, G1_bot_IR = radiosity_window(f,u,bnode,data)
        
#         G3_vis = rho3_vis*G1_bot_vis

#         G3_IR = eps3*ph"σ"*u[iT]^4 + rho3_IR*G1_bot_IR

#         f[iT] = alpha1_vis*G3_vis + alpha1_IR*G3_IR
#     end
# end

# function flux_abs_bottom(f,u,bnode,data)
#     if bnode.region==Γ_bottom

#         (;iT,iTp,lc_frit,lc_plate) = data
		
# 		# irradiation exchange between porous frit (1) and Al bottom plate (2)
# 		# porous frit properties (1)
#         eps1=lc_frit.eps
#         rho1_IR=lc_frit.rho_IR
#         # Al bottom plate properties (2)
#         # plate properties
# 		alpha2_IR=lc_plate.alpha_IR
# 		eps2=lc_plate.eps
# 		rho2_IR=lc_plate.rho_IR	
        
#         # local Al bottom plate temperature
#         Tplate = u[iTp]
#         G1_IR = (eps1*ph"σ"*u[iT]^4 + rho1_IR*eps2*ph"σ"*Tplate^4)/(1-rho1_IR*rho2_IR)

#         f[iT] = alpha2_IR*G1_IR
#     end
# end

function flux_enthalpy_top(f,u,bnode,data)
    (;iT,Fluids,X0,mfluxin,mmix0,T_gas_in,inlet_boundaries) = data

    if bnode.region in inlet_boundaries
        flux_enth_top = mfluxin/mmix0 * enthalpy_mix(Fluids, T_gas_in, X0)
        f[iT]=flux_enth_top
    end
end

function flux_mass_top(f,u,bnode,data)
    if bnode.region==Γ_top_in        
		(;ip,mfluxin) = data
        f[ip]=mfluxin
    end
end

# function flux_enthalpy_bottom(f,u,bnode,data)
#     if bnode.region==Γ_bottom
# 		(;iT,ip,Fluids,ubot,Tamb) = data
#         ng=ngas(data)
		
#         X=zeros(eltype(u), ng)
#         mole_frac!(bnode,data,X,u)

#         hmix=enthalpy_mix(Fluids, u[iT], X)
# 		flux_enth_bot=ubot*u[ip]/(ph"R"*u[iT])*hmix

# 		#flux_enth_bot=ubot*u[ip]/(ph"R"*u[iT])*cf*(u[iT]-Tamb)
#         f[iT]=flux_enth_bot
#     end
# end

# function flux_enthalpy_reaction(f,u,node,data)
    
#     if node.region == 2 && data.isreactive # catalyst layer
# 		(;iT,Fluids,kinpar,lcats) = data
# 		(;nuij,rni) = kinpar
#         ng=ngas(data)
		
# 		pi = u[1:ng]./ufac"bar"
# 		# negative sign: sign convention of VoronoiFVM: source term < 0 
# 		# unit conversion factor needs to be applied to Xu&Froment 1989 (mol/(hr*g))
# 		# RR = -lcats*ri(kinpar,u[iT],pi)
#         RR = -lcats*ri(data,u[iT],pi)
# 		## Xu & Froment 1989 kinetics
# 		# R1: CH4 + H2O = CO + 3 H2
# 		# R2: CO + H2O = CO2 + H2
# 		# R3: CH4 + 2 H2O = CO2 + 4 H2

# 		## Vazquez 2017 / S3P kinetics
# 		# R1: CO + H2O = CO2 + H2
# 		# R2: CH4 + 2 H2O = CO2 + 4 H2
# 		# R3: CH4 + H2O = CO + 3 H2
		
#         DHi = 0.0
#         RR_ = 0.0
#         if kinpar == XuFroment1989
#             DHi = -kinpar.ΔHi[:R2] # reaction is written in the reverse direction
# 			# unit conversion factor, needs to be applied to Xu&Froment 1989
#             RR_ = RR[rni[:R2]]*ufac"mol/(hr*g)" 
#         elseif kinpar == S3P
#             DHi = -kinpar.ΔHi[:R1] # reaction is written in the reverse direction
#             RR_ = RR[rni[:R1]]
#         end

#         f[iT] = RR_*DHi
# 	end
# end

# function areas(sol,sys,grid,data)
# 	iT = data.iT
# 	function area(f,u,bnode,data)
# 		# repurpose temperature index to hold area information
# 		f[iT] = one(eltype(u))
# 	end

# 	integrate(sys,area,sol; boundary=true)[iT,:]
# end

# function T_avg(sol,sys,grid,data)

# 	iT = data.iT	
# 	function T_avg_(f,u,bnode,data)			
# 		f[iT] = u[iT]		
# 	end

# 	areas_=areas(sol,sys,grid,data)
	
# 	T_int=integrate(sys,T_avg_,sol; boundary=true)[iT,:]
# 	T_int./areas_	
# end

# function HeatFluxes_Test(solt,grid,sys,data)
#     (;iT,irradiated_boundaries,outlet_boundaries) = data
#     inf_t, outf_t = TotalFlows(solt,grid,sys,data)

#     dE_dt = inf_t + outf_t
#     sol = solt(solt.t[end])
#     Qabs_inner_top = integrate(sys,irrad_top,sol; boundary=true)[iT,irradiated_boundaries]
#     Qabs_inner_top = sum(Qabs_inner_top)

#     QG_01=integrate(sys,flux_window_underside,sol; boundary=true)[iT,irradiated_boundaries]
#     QG_01=sum(QG_01)

#     QG_10=integrate(sys,flux_catalyst_layer,sol; boundary=true)[iT,irradiated_boundaries]
#     QG_10=sum(QG_10)

#     Qabs_inner_bottom = integrate(sys,irrad_bottom,sol; boundary=true)[iT,outlet_boundaries]
#     Qabs_inner_bottom = sum(Qabs_inner_bottom)

#     QG_34=integrate(sys,flux_radiation_frit_bottom,sol; boundary=true)[iT,outlet_boundaries]
#     QG_34=sum(QG_34)
#     # #QG_43=integrate(sys,flux_radiation_plate_bottom,sol; boundary=true)[iT,Γ_bottom]
#     QG_43=integrate(sys,flux_radiation_plate_bottom,sol; boundary=true)[iT,outlet_boundaries]
#     QG_43=sum(QG_43)

#     dH_dot = dE_dt - Qabs_inner_top
#     #dH_dot, Qabs_irr

#     # Qabs_inner_top, QG_01 - QG_10, QG_01, QG_10
#     Qabs_inner_bottom, QG_43 - QG_34, QG_43, QG_34
# end


function TotalFlows(solt,grid,sys,data)
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

    inflow_rate[end], outflow_rate[end]


end


function HeatFluxes_EB_I(solt,grid,sys,data)

    (;iT,inlet_boundaries,irradiated_boundaries,outlet_boundaries,side_boundaries)=data

    inf_t, outf_t = TotalFlows(solt,grid,sys,data)
    dE_dt = inf_t + outf_t
    sol = solt(solt.t[end])


    Qconv_10=integrate(sys,flux_convection_top,sol; boundary=true)[iT,irradiated_boundaries]
    Qconv_10=sum(Qconv_10)
    
    QG_01=integrate(sys,flux_window_underside,sol; boundary=true)[iT,irradiated_boundaries]
    QG_01=sum(QG_01)

    QG_10=integrate(sys,flux_catalyst_layer,sol; boundary=true)[iT,irradiated_boundaries]
    QG_10=sum(QG_10)

        
    QG_34=integrate(sys,flux_radiation_frit_bottom,sol; boundary=true)[iT,outlet_boundaries]
    QG_34=sum(QG_34)
    
    QG_43=integrate(sys,flux_radiation_plate_bottom,sol; boundary=true)[iT,outlet_boundaries]
    QG_43=sum(QG_43)

    Qconv_34=integrate(sys,flux_convection_bottom,sol; boundary=true)[iT,outlet_boundaries]
    Qconv_34=sum(Qconv_34)

    # # Qsides=integrate(sys,side,sol; boundary=true)[iT,[Γ_side_right,Γ_side_back,Γ_side_front,Γ_side_left]]
    Qsides=integrate(sys,MultiComponentReactiveMixtureProject.PTR_side,sol; boundary=true)[iT,side_boundaries]
    Qsides=sum(Qsides)

    dH_dot = dE_dt - (QG_01 - QG_10) - (QG_43 - QG_34) + Qconv_10 + Qconv_34 + Qsides

    # EB_bot= -Hout -Qcond_34 -QG_34 +QG_43    

    # EB = EB_top +EB_bot -Qsides    
    (
        QG_01=QG_01,
        # # QG_02=QG_02,
        QG_10=-QG_10,
        dH_dot=dH_dot,
        # # QG_20=-QG_20,
        Qconv_10=-Qconv_10,
        Qconv_34=-Qconv_34,
        QG_34=-QG_34,
        QG_43=QG_43,
        Qsides=-Qsides
    )
    

end

function HeatFluxes_EB_II(sol,grid,sys,data)


    # (;iT,Tglass,Tamb,k_nat_conv,uc_window)=data

    # Atop = sum(areas(sol,sys,grid,data)[[Γ_top_cat,Γ_top_frit]])
    # Qconv0 = Atop*k_nat_conv*(Tglass-Tamb)
    # Qemit0 = Atop*uc_window.eps*ph"σ"*Tglass^4

    (;iT,uc_window)=data
    Qconv0=integrate(sys,flux_convection_window_top,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    Qconv0=sum(Qconv0)
    Qemit0=integrate(sys,flux_radiation_window_top,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    Qemit0=sum(Qemit0)

    Qin=integrate(sys,flux_in_profile,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    Qin=sum(Qin)
    Qrefl0 = uc_window.rho_vis*Qin

    Qtrans_10=integrate(sys,flux_transmission_top_catalyst_layer,sol; boundary=true)[iT,Γ_top_cat]
    Qtrans_20=integrate(sys,flux_transmission_top_frit,sol; boundary=true)[iT,Γ_top_frit]

    # Qcond_10=integrate(sys,flux_conduction_top,sol; boundary=true)[iT,Γ_top_cat]
    Qconv_10=integrate(sys,flux_convection_top,sol; boundary=true)[iT,Γ_top_cat]
    # Qcond_20=integrate(sys,flux_conduction_top,sol; boundary=true)[iT,Γ_top_frit]
    Qconv_20=integrate(sys,flux_convection_top,sol; boundary=true)[iT,Γ_top_frit]

    QG_01=integrate(sys,flux_window_underside,sol; boundary=true)[iT,Γ_top_cat]
    QG_02=integrate(sys,flux_window_underside,sol; boundary=true)[iT,Γ_top_frit]

    QG_10=integrate(sys,flux_catalyst_layer,sol; boundary=true)[iT,Γ_top_cat]
    QG_20=integrate(sys,flux_frit,sol; boundary=true)[iT,Γ_top_frit]


    # EB_bot = -QG_01 +Qcond_10 +QG_10 -QG_02 +Qcond_20 +QG_20
    # EB_bot = -QG_01 +Qconv_10 +QG_10 -QG_02 +Qconv_20 +QG_20
    # EB_top = -Qconv0 -Qemit0 +Qin -Qrefl0 -Qtrans_10 -Qtrans_20
    # EB = EB_bot + EB_top

    (
        QG_01=-QG_01,
        Qconv_10=Qconv_10,
        QG_10=QG_10,
        QG_02=-QG_02,
        Qconv_20=Qconv_20,
        QG_20=QG_20,
        Qconv0=-Qconv0,
        Qemit0=-Qemit0,
        Qin=Qin,
        Qrefl0=-Qrefl0,
        Qtrans_10=-Qtrans_10,
        Qtrans_20=-Qtrans_20
    )


end


function HeatFluxes_EB_III(sol,grid,sys,data)



    # (;iT,Tamb,k_nat_conv,Tplate,lc_plate)=data
    # Abot = areas(sol,sys,grid,data)[Γ_bottom]
    # Qconv4 = Abot*k_nat_conv*(Tplate-Tamb)
	# Qemit4 = Abot*lc_plate.eps*ph"σ"*Tplate^4
    
    (;iT)=data

    Qconv4 = integrate(sys,flux_convection_outer_plate,sol; boundary=true)[iT,Γ_bottom]
	Qemit4 = integrate(sys,flux_emission_outer_plate,sol; boundary=true)[iT,Γ_bottom]

    # Qcond_34=integrate(sys,flux_conduction_bottom,sol; boundary=true)[iT,Γ_bottom]
    Qconv_34=integrate(sys,flux_convection_bottom,sol; boundary=true)[iT,Γ_bottom]

    QG_34=integrate(sys,flux_radiation_frit_bottom,sol; boundary=true)[iT,Γ_bottom]
    QG_43=integrate(sys,flux_radiation_plate_bottom,sol; boundary=true)[iT,Γ_bottom]


    (
        Qconv_34=Qconv_34,
        QG_43=-QG_43,
        QG_34=QG_34,
        Qconv4=-Qconv4,
        Qemit4=-Qemit4
    )
    
    
end

# outer energy balance: also consider transmitted and reflected irradiation, that
# will not change the inner energy of the system
function HeatFluxes_EB_IV_outer(sol,grid,sys,data)
    # (;iT,k_nat_conv,Tglass,Tamb,uc_window)=data

    # balance of upper side of window
    # Atop = sum(areas(sol,sys,grid,data)[[Γ_top_cat,Γ_top_frit]])
    # Qconv0 = Atop*k_nat_conv*(Tglass-Tamb)
    # Qemit0 = Atop*uc_window.eps*ph"σ"*Tglass^4
    (;iT,uc_window)=data

    Qconv0=integrate(sys,flux_convection_window_top,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    Qconv0=sum(Qconv0)
    Qemit0=integrate(sys,flux_radiation_window_top,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    Qemit0=sum(Qemit0)
    Qin=integrate(sys,flux_in_profile,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    Qin=sum(Qin)
    Qrefl0 = uc_window.rho_vis*Qin

    Qtrans_10=integrate(sys,flux_transmission_top_catalyst_layer,sol; boundary=true)[iT,Γ_top_cat]
    Qtrans_20=integrate(sys,flux_transmission_top_frit,sol; boundary=true)[iT,Γ_top_frit]
    # EB_top = -Qconv0 -Qemit0 +Qin -Qrefl0 -Qtrans_10 -Qtrans_20

    # balance of lower side of window        
    # Qcond_10=integrate(sys,flux_conduction_top,sol; boundary=true)[iT,Γ_top_cat]
    Qconv_10=integrate(sys,flux_convection_top,sol; boundary=true)[iT,Γ_top_cat]
    # Qcond_20=integrate(sys,flux_conduction_top,sol; boundary=true)[iT,Γ_top_frit]
    Qconv_20=integrate(sys,flux_convection_top,sol; boundary=true)[iT,Γ_top_frit]

    QG_10=integrate(sys,flux_catalyst_layer,sol; boundary=true)[iT,Γ_top_cat]
    QG_20=integrate(sys,flux_frit,sol; boundary=true)[iT,Γ_top_frit]

    QG_0=integrate(sys,flux_window_underside,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    QG_0=sum(QG_0)

    QG_0_IR=integrate(sys,flux_window_underside_IR,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    QG_0_IR=sum(QG_0_IR)

    QG_0_vis=integrate(sys,flux_window_underside_vis,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    QG_0_vis=sum(QG_0_vis)

    
    (
        Qconv0=-Qconv0,
        Qemit0=-Qemit0,
        Qin=Qin,
        Qrefl0=-Qrefl0,
        Qtrans_10=-Qtrans_10,
        Qtrans_20=-Qtrans_20,
        Qconv_10=Qconv_10,
        Qconv_20=Qconv_20,
        QG_10=QG_10,
        QG_20=QG_20,
        QG_0=-QG_0
    )


end


function HeatFluxes_EB_IV_inner(sol,grid,sys,data)
    (;iT)=data

        
    Qabs_10 = integrate(sys,flux_abs_top_catalyst_layer,sol; boundary=true)[iT,Γ_top_cat]
    Qabs_20 = integrate(sys,flux_abs_top_frit,sol; boundary=true)[iT,Γ_top_frit]

    # Qcond_10=integrate(sys,flux_conduction_top,sol; boundary=true)[iT,Γ_top_cat]
    Qconv_10=integrate(sys,flux_convection_top,sol; boundary=true)[iT,Γ_top_cat]
    # Qcond_20=integrate(sys,flux_conduction_top,sol; boundary=true)[iT,Γ_top_frit]
    Qconv_20=integrate(sys,flux_convection_top,sol; boundary=true)[iT,Γ_top_frit]

    # Atop = sum(areas(sol,sys,grid,data)[[Γ_top_cat,Γ_top_frit]])

    # Qconv0 = Atop*k_nat_conv*(Tglass-Tamb)
    # Qemit0 = Atop*uc_window.eps*ph"σ"*Tglass^4		

    Qconv0=integrate(sys,flux_convection_window_top,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    Qconv0=sum(Qconv0)
    Qemit0=integrate(sys,flux_radiation_window_top,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    Qemit0=sum(Qemit0)

    (
        Qconv0=-Qconv0,
        Qconv_10=Qconv_10,
        Qconv_20=Qconv_20,
        Qemit0=-2*Qemit0,
        Qabs_10=Qabs_10,
        Qabs_20=Qabs_20
    )

end


function HeatFluxes_EB_V(sol,grid,sys,data)

    # (;iT,Tamb,k_nat_conv,Tplate,lc_plate)=data
    # Abot = areas(sol,sys,grid,data)[Γ_bottom]
    # Qconv4 = Abot*k_nat_conv*(Tplate-Tamb)
	# Qemit4 = Abot*lc_plate.eps*ph"σ"*Tplate^4
    
    (;iT)=data

    Qconv4 = integrate(sys,flux_convection_outer_plate,sol; boundary=true)[iT,Γ_bottom]
	Qemit4 = integrate(sys,flux_emission_outer_plate,sol; boundary=true)[iT,Γ_bottom]

    # Qcond_34=integrate(sys,flux_conduction_bottom,sol; boundary=true)[iT,Γ_bottom]
    Qconv_34=integrate(sys,flux_convection_bottom,sol; boundary=true)[iT,Γ_bottom]

    Qabs_34 = integrate(sys,flux_abs_bottom,sol; boundary=true)[iT,Γ_bottom]
    
    (
        Qconv_34=Qconv_34,
        Qabs_34=Qabs_34,
        Qconv4=-Qconv4,
        Qemit4=-2*Qemit4,
    )

end

function HeatFluxes_EB_VI(sol,grid,sys,data)

    # (;iT,k_nat_conv,Tglass,Tplate,Tamb,uc_window,lc_plate)=data    
    (;iT,uc_window)=data    


    #### top boundary ####
    # Atop = sum(areas(sol,sys,grid,data)[[Γ_top_cat,Γ_top_frit]])

    # Qconv0 = Atop*k_nat_conv*(Tglass-Tamb)
    # Qemit0 = Atop*uc_window.eps*ph"σ"*Tglass^4
    Qconv0=integrate(sys,flux_convection_window_top,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    Qconv0=sum(Qconv0)
    Qemit0=integrate(sys,flux_radiation_window_top,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    Qemit0=sum(Qemit0)

    Qin=integrate(sys,flux_in_profile,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    Qin=sum(Qin)
    Qrefl0 = uc_window.rho_vis*Qin
    Qtrans_10=integrate(sys,flux_transmission_top_catalyst_layer,sol; boundary=true)[iT,Γ_top_cat]
    Qtrans_20=integrate(sys,flux_transmission_top_frit,sol; boundary=true)[iT,Γ_top_frit]

    #### side boundary (including enthalpy fluxes) ####
    Hin=integrate(sys,flux_enthalpy_top,sol; boundary=true)[iT,[Γ_top_cat,Γ_top_frit]]
    Hin=sum(Hin)
    Hout=integrate(sys,flux_enthalpy_bottom,sol; boundary=true)[iT,Γ_bottom]

    Qsides=integrate(sys,side,sol; boundary=true)[iT,[Γ_side_right,Γ_side_back,Γ_side_front,Γ_side_left]]
    Qsides=sum(Qsides)

    #### bottom boundarty ####
    # Abot = areas(sol,sys,grid,data)[Γ_bottom]

    # Qconv4 = Abot*k_nat_conv*(Tplate-Tamb)
	# Qemit4 = Abot*lc_plate.eps*ph"σ"*Tplate^4		
    Qconv4 = integrate(sys,flux_convection_outer_plate,sol; boundary=true)[iT,Γ_bottom]
	Qemit4 = integrate(sys,flux_emission_outer_plate,sol; boundary=true)[iT,Γ_bottom]

    DH = Hout -Hin    
    DH_reaction = integrate(sys,flux_enthalpy_reaction,sol)[iT,2]

    # calc global energy balance VI
    # EB_top = -Qconv0 -Qemit0 +Qin -Qrefl0 -Qtrans_10 -Qtrans_20
    # EB_sides = Hin -Hout -Qsides
    # EB_bot = -Qconv4 -Qemit4
    
    # EB = EB_top +EB_sides +EB_bot

    (
        Qin=Qin,
        Qconv0=Qconv0,
        Qemit0=Qemit0,
        Qrefl0=Qrefl0,
        Qtrans_10=Qtrans_10,
        Qtrans_20=Qtrans_20,
        Qsides=Qsides,
        DH=DH,
        DH_reaction=DH_reaction,
        Qconv4=Qconv4,
        Qemit4=Qemit4
    )
end



function HeatFluxes_EB_VI_plot(heatflows)
    (;Qin,Qconv0,Qemit0,Qrefl0,Qtrans_10,Qtrans_20,Qsides,DH,DH_reaction,Qconv4,Qemit4) = heatflows

    DH_sens = DH-DH_reaction

    labels = ["Top Trans", "Top Refl", "Top Conv", "Top Emit",  "Sides Conv", "Bottom Conv", "Bottom Emit", "ΔHSensible","ΔHReaction"]
    
    x = [string(x) for x in 1:length(labels)]

    labelsx = [x*": "* l for (l,x) in zip(labels,x)]

    flows=[(Qtrans_10+Qtrans_20),Qrefl0,Qconv0,Qemit0,Qsides,Qconv4,Qemit4,DH_sens,DH_reaction]

    y = flows / Qin * 100.0
        
    p=Plots.plot(size=(400,300), yguide="Heat flow contribution / %", grid=false, 
    top_margin=5mm)
    Plots.plot!(p, permutedims(x), permutedims(y); st=:bar, label=permutedims(labelsx))    
    Plots.annotate!(p, x, y, round.(y, digits=1), :bottom)
    Plots.annotate!(p, 12.0, 22.5,
        "Σ flows = "*string(round(sum(flows)/Qin*100.0,digits=2))*" %", :right
        )
    # Plots.savefig(p,"img/out/HeatFlowDistr.svg")
    Plots.savefig(p,"img/out/HeatFlows_baseCase230607.pdf")

end


end
