!    !==================== test variables
!     real fwsoil_yr,omega_yr,topfws_yr, difference,diff_yr,diff_d
            
!       integer, parameter :: iiterms=7            ! 9 for Duke forest FACE
!       integer, parameter :: ilines=150000         ! the maxmum records of Duke Face, 1998~2007
!       real, parameter:: times_storage_use=720.   ! 720 hours, 30 days
!       integer  lines,idays,MCMC,do_co2_da
!       integer,dimension(ilines):: year_seq,doy_seq,hour_seq
!       real forcing_data(iiterms,ilines),input_data(iiterms,ilines)
! !   *** before ..int
! !      real Simu_dailyflux(12,10000)
!       real Simu_dailyflux14(14,10000)
!       real Simu_dailywater(31,10000)
! !      real obs_spruce(12,1000)
! !      integer pheno,phenoset
! !      site specific parameters
! !   *** 
! !   *** after ..int
!       real Simu_dailyflux(12,80000),Simu_soilwater(10,40000),Simu_soiltemp(11,40000)
!       real Simu_watertable(1,40000),Simu_dailysoilt(11,80000),Simu_dailywatertable(1,80000)
!       real Simu_dailyice(10,80000)
!       real Simu_snowdepth(1,80000)
!       real water_table(ilines),snow_in(ilines)      
!       real obs_spruce(12,1000),obs_soilwater(5,40000),obs_soilt(11,40000)
!       integer pheno,phenoset,day_mod,num
! !   ***      
      
!       real lat,longi,rdepth,LAIMAX,LAIMIN
!       real wsmax,wsmin,co2ca,CO2treat
!       real tau_L,tau_W,tau_R
!       real tau_F,tau_C,tau_Micr,tau_Slow,tau_Pass
!       real TauC(8)
! !      the variables that should be initialized in the begining
!       real Q_soil
!       real QC(8) !  leaf,wood,root,fine lit.,coarse lit.,Micr,Slow,Pass
!       real Pool1,Pool2,Pool3,Pool4,Pool5,Pool6,Pool7,Pool8
!       real out1_yr,out2_yr,out3_yr,out4_yr,out5_yr,out6_yr,out7_yr,out8_yr
!       real OutC(8)
!       real Rh_pools(5)
! !      for soil conditions
!       real WILTPT,FILDCP,infilt
!       real Rsoilabs
!       real fwsoil,topfws,omega
! !      for plant growth and allocation
!       real NSC,NSCmin,NSCmax,add               ! none structural carbon pool
!       real Growth,Groot,Gshoot,GRmax           ! growth rate of plant,root,shoot,and max of root
!       real St,Sw,Ss,Sn,Srs,Sps,fnsc,Weight     ! scaling factors for growth
! !      variables for canopy model


!       real evap,transp,ET,G

! !      real evap,transp,ET
! !!   *** ..int
      
!       real wind,eairp,esat,rnet
!       real Pa_air 
!       real gpp,gpp_ra,NPP,NEE,NEP,gpp_d,NPP_d
!       real evap_d,transp_d
!       real,dimension(3):: tauL,rhoL,rhoS,reffbm,reffdf,extkbm,extkdm
!       real,dimension(2):: Radabv
!       real Qcan(3,2)
! !      parameters for photosynthesis model
!       real stom_n,a1,Ds0,Vcmx0,Vcmax0,extkU,xfang,alpha
!       real pi,emleaf,emsoil
!       real Rconst,sigma,cpair,Patm,Trefk,H2OLv0,airMa,H2OMw,chi,Dheat
!       real wleaf,gsw0,eJmx0,theta,conKc0,conKo0,Ekc,Eko,o2ci
!       real Eavm,Edvm,Eajm,Edjm,Entrpy,gam0,gam1,gam2
! !     for nitrogen sub-model
!       real CNmin,CNmax,NSNmax,NSNmin
!       real NSN
! !      QNleaf,QNwood,QNroot,QNfine,QNcoarse,QNmicr,QNslow,QNpass
! !      CN_leaf,CN_wood,CN_root,CN_fine,CN_coarse,CN_micr,CN_slowC,CN_pass
    !   real QN(8),CN0(8),CN(8),OutN(8),QNplant,QNminer
!       real N_leaf,N_wood,N_root,N_deficit
!       real N_LF,N_WF,N_RF
!       real N_uptake,N_leach,N_vol,N_fixation,N_deposit,N_fert
!       real N_up_d,N_fix_d,N_dep_d,N_leach_d,N_vol_d
!       real N_up_yr,N_fix_yr,N_dep_yr,N_leach_yr,N_vol_yr
!       real N_miner,alphaN
!       real SNvcmax,SNgrowth,SNRauto,SNrs
! !   ***   .. int add pars for soil thermal
!       real diff_s,diff_snow,albedo_snow,resht,thd_snow_depth,shcap_snow,condu_snow,depth_ex
!       real infilt_rate
!       real b_bound,fa,fsub,rho_snow,decay_m,condu_b
! !   ***      
! !      additional arrays to allow output of info for each layer
!       real,dimension(5):: RnStL,QcanL,RcanL,AcanL,EcanL,HcanL
!       real,dimension(5):: GbwcL,GswcL,hG,hIL
!       real,dimension(5):: Gaussx,Gaussw,Gaussw_cum 
! !      for phenology
!       real LAI,bmroot,bmstem,bmleaf,bmplant,totlivbiom,ht
!       real SLA,SLAx,L_fall,L_add,litter,seeds
!       real GDDonset,GDD5,accumulation,storage,stor_use,store
!       real RaL,RaS,RaR  !allocation to respiration
!       real alpha_L,alpha_W,alpha_R ! allocation ratio to Leaf, stem, and Root
!       real Q10,Rl0,Rs0,Rr0         ! parameters for auto respiration
!       real Rgrowth,Rnitrogen,Rmain,Rauto !respirations
!       real RmLeaf,RmStem,RmRoot          ! maintanence respiration
!       real RgLeaf,RgStem,RgRoot          ! growth respiration
!       real RaLeaf,RaStem,RaRoot
!       real Rsoil,Rhetero,Rtotal
!       real Ra_Nfix,Rh_Nfix
!       real gpp_yr,NPP_yr,NEE_yr,RaL_yr,RaR_yr,RaS_yr,Rh_yr
!       real Rh4_yr,Rh5_yr,Rh6_yr,Rh7_yr,Rh8_yr,Ra_yr
!       real R_Ntr_yr
!       real NPPL_yr,NPPR_yr,NPPS_yr,NPP_L,NPP_R,NPP_W
!       real Rootmax,Stemmax,SapS,SapR,StemSap,RootSap
!       REAL ws,wdepth
! !      climate variables for every day
! !   *** before ..int
! !      real Ta,Tair,Ts,Tsoil,Ttreat
! !   *** after ..int
!       real Ta,Tair,Ts,Tsoil,Ttreat,water_table_depth,snow_depth      
! !   ***
!       real doy,hour,Dair,Rh,radsol
!       real PAR
! !      output daily means of driving variables
!       real CO2air_d_avg,SWdown_d_avg,Psurf_d_avg
!       real Rain_d_avg,Tair_d_avg,Wind_d_avg
! !   *** ..int
!       real Simu_dailyCH4(16,80000)
! !   ***
! !      output from canopy model
!       real evap_yr,transp_yr
!       real,dimension(10):: thksl,wupl,evapl,wcl,FRLEN   ! wsc is the output from soil water module
!       real wsc(10)
!       real runoff,runoff_d,runoff_yr,rain,rain_d,rain_yr
!       real ws1,ws2,dws,net_dws
!       real Esoil,Hcrop,ecstot,Anet,DEPH2O,Acanop
!       real Hcanop,Hcanop_d
!       real Raplant,Glmax,Gsmax,Rh_d
!       real GLmx,Gsmx,GRmx
! !      output for ORNL model comparison
!       real CO2h,PARh,ATh,STh,VPDh,SWh
!       real RECOh
!       real ETh,Th,Eh,INTh,ROh,DRAINh,LEh,SHh
!       real LWH,Rgrowth_d,abvLitter,blLitter
! !     daily output
!       real PAR_d,AT_d,ST_d,VPD_d
!       real SW_d,NEP_d,NEE_d,RECO_d
!       real Ra_d,RLEAV_d,RWOOD_d,RROOT_d,RHET_d,RSOIL_d,ET_d,T_d
!       real E_d,INT_d,RO_d,DRAIN_d,LE_d,SH_d,CL_d,CW_d,CFR_d,TNC_d
!       real CSOIL_d,GL_d,GW_d,GR_d,LFALL_d,LMA_d,NCAN_d,NWOOD_d
!       real GL_yr,GR_yr,GW_yr
!       real NFR_d,NSOIL_d,NUP_d,NMIN_d,NVOL_d,NLEACH_d
!       real N_LG_d,N_WG_d,N_RG_d
!       real N_LF_d,N_WF_d,N_RF_d
!       real WFALL_D,RFALL_D
!       real Simu_lit
      
! !   *** added for ..int  
!       ! for soil temp
!       real sftmp,Tsnow,Twater,Tice,ice_tw,water_tw 
!       real,dimension(10):: Tsoill,ice,liq_water
!       real,dimension(11):: testout
!       real soilt_d_simu(11),soilt_d_obs(7),watertable_d_obs,ice_d_simu(10)
!       integer obs_counter(7) 
!       real zwt_d,snow_depth_e,snow_dsim,melt,dcount,dcount_soil
!       character(len=80) outfile

!       integer dlayer      
! !   *** added for ..int      
      
! !      NEE observation
!       real NEE_annual,Cumol2gram
!       real NEE_annual_array(30)
!       integer year_array(30),year_obs
! !     for loops
!       integer jrain,W_flag(7)
!       integer onset !flag of phenological stage
!       integer year,yr,days,i,j,m,n,yrs_eq,hoy,iyr,daily
!       integer k1
!       integer lines_NEE,yr_NEE
!       integer istat1,istat2,istat3,istat4
!       integer dtimes,yr_length
!       integer num_scen,isite
!       integer idoy,ihour,ileaf,first_year
!       integer dylim,yrlim
!       real zwt,phi

! !   *** ..int
!       !*****for methane subroutine      MS
!       integer, parameter :: nlayers=10
!       real CH4(nlayers),CH4_V(nlayers),CH4V_d(nlayers)                                   !MS
! !      real CH4(nlayers),CH4_V(nlayers+1),CH4V_d(nlayers)                                   !MS
!       real ProCH4(nlayers),Pro_sum,Pro_sum_d,Pro_sum_yr
!       real OxiCH4(nlayers),Oxi_sum,Oxi_sum_d,Oxi_sum_yr
!       real simuCH4,simuCH4_d,simuCH4_yr
!       real Fdifu(nlayers+1),Fdifu1_d,Fdifu1_yr
!       real Ebu_sum,Ebu_sum_d,Ebu_sum_yr
!       real Pla_sum,Pla_sum_d,Pla_sum_yr
!       real S_omega
!       real r_me,Q10pro,kCH4,Omax,CH4_thre,Tveg,Tpro_me,Toxi
!       logical do_snow,do_soilphy
!       real Ebu_sum_sat, Ebu_sum_unsat