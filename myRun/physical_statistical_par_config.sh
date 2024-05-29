#!/bin/tcsh

set data_offset_histos 		= 7
set data_offset_bin 		= 10
set data_offset_cSFRD 		= 5
set data_offset_sfr2z 		= 98
set data_offset_ngal 		= 10
set data_offset_CUTE 		= 6
set data_offset_analyseTargetSelection 		= 10

set nbins_mstar			= 60
set nbins_sfr 			= 60
set nbins_ssfr 			= 60
set nbins_cSFRD 		= 18
set nbins_sfr2z 		= 39
set nbins_sfr2mstar 		= 80
set nbins_ssfr2mstar 		= 80
set nbins_mhalo 		= 54
set nbins_mstar2mhalo	 	= 60
set nbins_mstar2mhalofunc	= 25
set nbins_ngal2mhalo 		= 12
set nbins_mcold2mstar 		= 80
set nbins_mbh2mstarsph 		= 80
set nbins_twoPCF		= 20
set nbins_zgas2mstar		= 25
set nbins_zgas2mcold		= 35


#set observational data should be included in the plots
set obs_mstar 		    = 'Moustakas13_SDSS_GALEX_z010'
#set obs_mstar 		    = 'SAGE_load_add2,SAGE_load_add,Moustakas13_SDSS_GALEX_z010'
set obs_mstar 		    = 'SMF_CMASS_spall_por_mer_050_DR12v4'
#set obs_mstar 		    = 'False'
#,SMF_wis_055_DR12v4,SMF_gra_055_DR12v4'
#set obs_mstar 		    = 'DEEP2_FireFly'
#set obs_mstar 		    = 'ZFOURGE_150_z_250_SF'
#set obs_mstar		    = 'Rodriguez15SMF_CMASS'
#set obs_mstar		    = 'Baldry_08_SMF'
#set obs_mstar 		    = 'Maraston13_BOSS_CMASS_model'
#set obs_mstar		    = 'False'

#set obs_sfr 		    = 'CARNage_sfr_Gruppioni15'
set obs_sfr 		    = 'Gruppioni15_z_010'
set obs_sfr 		    = 'False'
set obs_ssfr 		    = 'False'
set obs_cSFRD		    = 'Behroozi13_sfr2z,Maudau14_sfr2z,Driver18_sfr2z,Gruppioni15_cSFRD'
#set obs_cSFRD		    = 'False'
set obs_sfr2z		    = 'False'

set obs_ssfr2mstar 	    = 'Elbaz11_dots,Gal_ssfr2mstar,Gal2_ssfr2mstar'
#set obs_ssfr2mstar 	    = 'zevol1,zevol2,zevol3,zevol4,zevol5,Elbaz11'
#set obs_ssfr2mstar 	    = 'False'
set obs_mhalo 		    = 'False'
set obs_mstar2mhalo 		    = 'False'
#set obs_mstar2mhalo 	    = 'Behroozi13c'
#,Rodriguez15BMstar2Mhalo,Rodriguez16_SMHM_obs_plus_dy'
#set obs_mstar2mhalo 	    = 'Behroozi10'
set obs_mstar2mhalofunc	    = 'False'
set obs_mstar2mhalovsSFR    = 'False'
set obs_mhalo2ngal 	    = 'False'
set obs_zgas2mstar	    = 'Tremonti_04_OH_mstar'
set obs_zgas2mcold	    = 'False'
set obs_oh2mstar	    = 'Tremonti_04_OH_mstar'
set obs_oh2mstar	    = 'SDSS_ELG_z01'
set obs_zgas2mstar	    = 'SDSS_ELG_z01'
set obs_mstar2mcold 	    = 'CARNage_cold_gas_fraction_Peeples14'
set obs_mstar2mcold 	    = 'Boselli14_Mcold2Mstar'
set obs_mbh2mstarsph	    = 'CARNage_mbh_mbulge_Kormendy13,CARNage_mbh_mbulge_McConnell13'
set obs_twoPCF	    	= 'False'
set obs_mags	    	= 'Blanton03_LF_u'
set obs_HOD		= 'load_the_shit'
#set obs_mstar2rhalf	= 'False'


if ($HOME == /home/doris) then
	#mcold/mstar vs. mstar
	set obs_analyse			= 'myCatalog,Boselli14_Mcold2Mstar,CARNage_cold_gas_fraction_Peeples11'
	#ssfr vs. mstar
	set obs_analyse			= 'myCatalog,Elbaz11_function,Gal2_ssfr2mstar,Gal_ssfr2mstar'
else
	set obs_analyse			= 'CMASS_spall_por_mer_056_DR12v4,SMF_CMASS_spall_por_mer_050_DR12v4'
endif

set obs_analyse		= 'CMASS_Galacticus_down3_low_catalog,CMASS_Galacticus_down3_passive_catalog,CMASS_Galacticus_down3_red_catalog,CMASS_Galacticus_down3_lowZ_catalog,CMASS_Galacticus_down3_highZ_catalog'
#set obs_analyse		= 'myCatalog2,myCatalog6,myCatalog7,Salinas21_LSBGs_Fig6_jstar_mstar,Salinas21_LSBGs_Fig6_fit,Salinas21_LSBGs_Fig6_fit_dy_lower,Salinas21_LSBGs_Fig6_fit_dy_higher'
#set obs_analyse		= 'myCatalog2,Salinas21_LSBGs_Fig6_jstar_mstar,Salinas21_LSBGs_Fig6_fit,Salinas21_LSBGs_Fig6_fit_dy_lower,Salinas21_LSBGs_Fig6_fit_dy_higher'
#set obs_analyse		= 'LG_catalog,red_blue_cut,Montero16_RS'
#set obs_analyse		= 'myCatalog,angM0,angM,angM2'
#set obs_analyse			= 'CMASS_Galacticus_color_catalog,CMASS_SAGv4_color_catalog,SAG_run2_catalog,CMASS_dmesa_cut'
#,CMASS_dmesa_cut,red_blue_cut,Montero16_RS'
#set obs_analyse			= 'CMASS_Galacticus_color_catalog,CMASS_spall_por_mer_050_DR12v4'
#set obs_analyse		= 'CMASS_Galacticus_color_catalog,CMASS_Galacticus_down3_catalog'
#,CMASS_Galacticus_down3_crossed_catalog'
#,CMASS_Galacticus_mass_catalog'
#set obs_analyse		= 'CMASS_Galacticus_color_catalog,CMASS_Galacticus_down3_catalog,CMASS_Galacticus_mass_catalog,CMASS_spall_por_mer_050_DR12v4'
#,CMASS_cut'
#,CMASS_spall_por_PS_050_DR12v4'
#,CMASS_dmesa_cut,red_blue_cut,Montero16_RS'
#set obs_analyse		= 'CMASS_Galacticus_density_catalog,CMASS_Galacticus_down3_catalog,CMASS_Galacticus_down_catalog,CMASS_spall_por_mer_050_DR12v4,CMASS_cut'
#,CMASS_dmesa_cut,red_blue_cut,Montero16_RS'
#set obs_analyse		= 'CMASS_SAG_color_catalog,CMASS_SAG_density_catalog,CMASS_SAG_mass_catalog,SAG_catalog'
#,CMASS_spall_por_mer_056_DR12v4,CMASS_dmesa_cut,red_blue_cut,Montero16_RS'
#set obs_plotXY			= 'Bower06_xir2_CUT2_mstar,DeLucia07_xir2_CUT2_mstar'
#set obs_plotXY			= 'SAGE_z,OBS_SAGE_z,Montero09_LF_z'
#set obs_plotXY			= 'Montero16_LF_i,CMASS_DR12_LF_RS_050'
#,Bernadi16_LF_CMASS'
set obs_plotXY			= 'SMF_CMASS_spall_PS_056_DR12v4,SMF_CMASS_spall_por_SF_056_DR12v4,SMF_CMASS_spall_por_mer_050_DR12v4'

set obs_plotXY			= 'wp_CMASS_DR12_N_cents'
#set obs_plotXY			= 'wp_CMASS_DR12_all_1125_1136,wp_CMASS_DR12_all_1136_1148,wp_CMASS_DR12_all_1148'
#,wp_CMASS_DR12_N_all_1145,wp_CMASS_DR12_N_all_1155'

#set obs_plotXY			= 'wp_CMASS_DR12_N_all_1125,wp_CMASS_DR12_N_all_1135,wp_CMASS_DR12_N_all_1145,wp_CMASS_DR12_N_all_1155,wp_CMASS_DR12_N_all_1165'
#set obs_plotXY			= 'SMF_wis_055_DR12v4,Guo18_SMF_CMASS'
#,Rodriguez16_SMF_z055_plus_dy'
#set obs_plotXY			= 'Behroozi13_sfr2z'
#set obs_plotXY			= 'LF_CMASS_spall_por_mer_050_DR12v4_Mi'
#,LF_CMASS_spall_por_PS_050_DR12v4_g,MLF_CMASS_spall_por_PS_050_DR12v4_g'
#set obs_mstar			= 'Maraston13_BOSS_CMASS_model'
#,SMF_por_merged_055_DR12v4'
#set obs_plotXY			= 'SMF_CMASS_spall_por_mer_050_DR12v4,SMF_CMASS_spall_por_PS_050_DR12v4,SMF_CMASS_spall_por_SF_050_DR12v4,SMF_CMASS_spall_gra_050_DR12v4,SMF_CMASS_spall_wis_050_DR12v4'
#set obs_plotXY			= 'Reid14_HOD_fit_all,Reid14_HOD_fit_cents,Reid14_HOD_fit_sats,BigMDPL_LC_HOD_all,BigMDPL_LC_HOD_cents,BigMDPL_LC_HOD_sats,HOD_all_dummy,HOD_cents_dummy,HOD_sats_dummy'
set obs_plotXY			= 'Reid14_HOD_fit_all,Reid14_HOD_fit_cents,Reid14_HOD_fit_sats,MD_LC_HOD_all,MD_LC_HOD_cents,MD_LC_HOD_sats,HOD_all_dummy,HOD_cents_dummy,HOD_sats_dummy'
set obs_plotXY			= 'MD_LC_HOD_all,MD_LC_HOD_cents,MD_LC_HOD_sats,HOD_all_dummy,HOD_cents_dummy,HOD_sats_dummy'
#set obs_plotXY			= 'mod_SMF_CMASS_spall_por_mer_050_DR12v4'
#set obs_plotXY			= 'SFRF_CMASS_spall_por_mer_050_DR12v4_cross'
#set obs_plotXY			= 'Zehavi11_22_21_wp'
#,Zehavi11_20_19_wp,Zehavi11_21_20_wp,Zehavi11_22_21_wp'
set obs_plotXY			= 'False'
#set obs_plotXY			= 'CUT1_dummy,CUT2_dummy,CUT3_dummy'
#set obs_plotXY			= 'format_MTs_Rockstar'
#set obs_plotXY			= 'load_the_shit'
#set obs_plotXY			= 'loadDiverse'
#set obs_plotXY			= 'Elbaz11_function'


#set physical boundary parameter used for calculations
set h_min_SMF_mstar			= 5e8		
set h_max_SMF_mstar			= 1e13
set h_min_SFRF_sfr 			= 1e-5		
set h_max_SFRF_sfr 			= 10000		
set h_min_sSFRF_ssfr			= 1e-15	
set h_max_sSFRF_ssfr			= 1e-7		
set h_min_sfr2z_mstar			= min		
set h_max_sfr2z_mstar			= max
set h_min_mstar2sfr			= min		
set h_max_mstar2sfr 			= max				
set h_min_mstar2ssfr 			= min		
set h_max_mstar2ssfr 			= max

set h_min_oh2mstar_mstar		= 1e9		
set h_max_oh2mstar_mstar		= 3.16e11
		
set h_min_mhalo_mhalo			= 8.9125094e9		
#set h_max_mhalo_mhalo			= 1.7782794e15
set h_max_mhalo_mhalo			= 2.2387211e15		
set h_min_mstar2mhalo_mhalo		= min		
set h_max_mstar2mhalo_mhalo		= max
set h_min_mstar2mhalo_sfr		= min		
set h_max_mstar2mhalo_sfr		= max		
set h_min_mh2mstfu_mhalo 		= min	
set h_max_mh2mstfu_mhalo 		= max	
set h_min_mh2mstfu_mstar2mhalo 		= min	
set h_max_mh2mstfu_mstar2mhalo 		= max
set h_min_ngal2mhalo_mstar 		= min		
set h_max_ngal2mhalo_mstar 		= max

set h_min_mhalo_no_mhalo		= 1e11		
set h_max_mhalo_no_mhalo		= 1e15
set h_min_mstar2mhalo_no_mhalo		= 8.9125094e9		
set h_max_mstar2mhalo_no_mhalo		= 2.2387211e15
set h_min_mstar2mhalo_no_mstar		= min		
set h_max_mstar2mhalo_no_mstar		= max		
set h_min_mh2mstfu_no_mhalo 		= min	
set h_max_mh2mstfu_no_mhalo 		= max	
set h_min_mh2mstfu_no_mstar2mhalo 	= 0.00001	
set h_max_mh2mstfu_no_mstar2mhalo 	= 1e7			
set h_min_ngal2mhalo_no_mstar 		= min		
set h_max_ngal2mhalo_no_mstar 		= max
set h_min_mstar2mhalovssfr_sfr		= 1e-3		
set h_max_mstar2mhalovssfr_sfr		= 1000
set h_min_mstar2mhalovssfr_ssfr	= 1e-15		
set h_max_mstar2mhalovssfr_ssfr	= 1e-6
set h_min_mstar2mhalovssfr_mhalo	= min		
set h_max_mstar2mhalovssfr_mhalo	= max		
		
set h_min_zgas2mstar_mstar 		= 1e9		
set h_max_zgas2mstar_mstar 		= 3.16e11		
set h_min_mcold2mstar_mstar		= min		
set h_max_mcold2mstar_mstar		= max		
set h_min_mbh2mstarsph_mstar_spheroid	= min		
set h_max_mbh2mstarsph_mstar_spheroid	= max

set h_min_twoPCF_galgalsep		= 5		
set h_max_twoPCF_galgalsep		= 200

set h_min_HOD_mhalo			= min		
set h_max_HOD_mhalo			= max

set h_min_mstar2rhalf_mstar		= min		
set h_max_mstar2rhalf_mstar		= max
set h_min_mstar2rhalf_rhalf		= min		
set h_max_mstar2rhalf_rhalf		= max

########################################################

set mc_min_SMF_mstar			= 5e8	
set mc_max_SMF_mstar			= 1e13
#OII
#set mc_min_SMF_mstar			= 5e8	
#set mc_max_SMF_mstar			= 1e13


set mc_min_SFRF_mstar			= min		
set mc_max_SFRF_mstar			= max

set mc_min_SFRF_sfr 			= 1e-5	
set mc_max_SFRF_sfr 			= 10000
	
set mc_min_sSFRF_mstar			= min
set mc_max_sSFRF_mstar			= max			
set mc_min_sSFRF_ssfr			= 1e-15	
set mc_max_sSFRF_ssfr			= 1e-7
		
set mc_min_sfr2z_mstar			= 7.4e8		
set mc_max_sfr2z_mstar			= max
set mc_min_mstar2ssfr_mstar 	= 1e8
set mc_max_mstar2ssfr_mstar 	= max		
set mc_min_mstar2ssfr_ssfr 		= min		
set mc_max_mstar2ssfr_ssfr 		= max
set mc_min_mstar2sfr_mstar 		= 1e3
set mc_max_mstar2sfr_mstar 		= max		
set mc_min_mstar2sfr_sfr 		= min		
set mc_max_mstar2sfr_sfr 		= max

set mc_min_oh2mstar_mstar 		= 1e9
set mc_max_oh2mstar_mstar 		= 3.16e11		
set mc_min_oh2mstar_sfr 		= min	
set mc_max_oh2mstar_sfr 		= max


set mc_min_mhalo_mhalo			= 8.9125094e9		
#set mc_max_mhalo_mhalo			= 1.7782794e15
set mc_max_mhalo_mhalo			= 2.2387211e15	
set mc_min_mstar2mhalo_mhalo		= 1e6		
set mc_max_mstar2mhalo_mhalo		= max		
set mc_min_mh2mstfu_mhalo 		= min	
set mc_max_mh2mstfu_mhalo 		= max	
set mc_min_mh2mstfu_mstar2mhalo		= 0.00001	
set mc_max_mh2mstfu_mstar2mhalo		= 1e7
set mc_min_ngal2mhalo_mstar 		= min		
set mc_max_ngal2mhalo_mstar 		= max

set mc_min_mhalo_no_mhalo		= 1e11		
set mc_max_mhalo_no_mhalo		= 1e15
set mc_min_mstar2mhalo_no_mhalo		= 8.9125094e9		
set mc_max_mstar2mhalo_no_mhalo		= 2.2387211e15
set mc_min_mstar2mhalo_no_mstar		= min		
set mc_max_mstar2mhalo_no_mstar		= max
set mc_min_mstar2mhalovssfr_sfr		= 1e-3		
set mc_max_mstar2mhalovssfr_sfr		= 1000
set mc_min_mstar2mhalovssfr_ssfr	= 1e-15		
set mc_max_mstar2mhalovssfr_ssfr	= 1e-6	
set mc_min_mstar2mhalovssfr_mhalo	= min		
set mc_max_mstar2mhalovssfr_mhalo	= 1e11
	
set mc_min_mh2mstfu_no_mhalo 		= min	
set mc_max_mh2mstfu_no_mhalo 		= max	
set mc_min_mh2mstfu_no_mstar2mhalo 	= 0.00001	
set mc_max_mh2mstfu_no_mstar2mhalo 	= 1e7			
set mc_min_ngal2mhalo_no_mstar 		= min		
set mc_max_ngal2mhalo_no_mstar 		= max
		
set mc_min_zgas2mstar_mstar 		= 1e9		
set mc_max_zgas2mstar_mstar 		= 3.16e11		
set mc_min_mcold2mstar_mstar		= 100		
set mc_max_mcold2mstar_mstar 		= max
set mc_min_mcold2mstar_mcold		= 100		
set mc_max_mcold2mstar_mcold 		= max		
set mc_min_mbh2mstarsph_mstar_spheroid	= 100		
set mc_max_mbh2mstarsph_mstar_spheroid	= max	
set mc_min_mbh2mstarsph_mbh		= 100		
set mc_max_mbh2mstarsph_mbh		= max

set mc_min_twoPCF_boxsize		= 100		
set mc_max_twoPCF_boxsize		= 900
set mc_min_twoPCF_mstar			= 1e9		
set mc_max_twoPCF_mstar 		= max
set mc_twoPCF_ngal	 		= 50000

set mc_min_r_i		 		= min
set mc_max_r_i		 		= 2
set mc_min_dmesa	 		= 0.55
set mc_max_dmesa	 		= max
set mc_min_g_r		 		= min
set mc_max_g_r		 		= max
set mc_min_i_lt_dmesa 			= min
set mc_max_i_lt_dmesa 			= max


set mc_min_HOD_mhalo			= min		
set mc_max_HOD_mhalo			= max

set mc_min_mstar2rhalf_mstar		= min		
set mc_max_mstar2rhalf_mstar		= max
set mc_min_mstar2rhalf_rhalf		= 0.01		
set mc_max_mstar2rhalf_rhalf		= max

set plot_array = ()
foreach item ($5)
	set plot_array = ($plot_array $item)
end

set name_array = ()
foreach item ($4)
	set name_array = ($name_array $item)
end

set col_array = ()
foreach item ($7)
	set col_array = ($col_array $item)
end

set analyse_tarsel_array = ()
foreach item ($9)
	set analyse_tarsel_array = ($analyse_tarsel_array $item)
end
#echo $name_array
#echo $plot_array

set i=1

if ("$plot_array[$i]" == 1) then
	#echo "$name_array[$i]" 'plot_'"$name_array[$i]"'_config.txt '0 0 $obs_plotXY 		>> $1
	#LF
	#echo "$name_array[$i]" 'plot_'"$name_array[$i]"'_config.txt '$nbins_mstar $data_offset_histos $obs_mstar >> $1
	#2PCF
	echo "$name_array[$i]" 'plot_'"$name_array[$i]"'_config.txt 30' $data_offset_CUTE $obs_plotXY >> $1

	echo 'histo_min_mass_'"$name_array[$i]"': '$h_min_SMF_mstar >> $3
	echo 'histo_max_mass_'"$name_array[$i]"': '$h_max_SMF_mstar >> $3
	echo 'mycond_min_mass_'"$name_array[$i]"': '$mc_min_SMF_mstar >> $3
	echo 'mycond_max_mass_'"$name_array[$i]"': '$mc_max_SMF_mstar >> $3
endif
@ i++

if ("$plot_array[$i]" == 1) then

	set k=1
	foreach item ($6)

		if ("$col_array[$k]" != 99) then

			if ($item =~ mstar) then
				set obs = $obs_mstar
			else if ($item =~ mA* || $item =~ MA*) then
				set obs = $obs_mags
			else
				set obs = 'False'
			endif

			echo "$name_array[$i]"'_'$item 'plot_'"$item"'_config.txt '25 $data_offset_histos $obs >> $1
		endif
		@ k++	
	end
	echo 'mycond_min_'"$name_array[$i]"'_r_i: '$mc_min_r_i 			>> $3
	echo 'mycond_max_'"$name_array[$i]"'_r_i: '$mc_max_r_i 			>> $3
	echo 'mycond_min_mc_min_'"$name_array[$i]"'_d_mesa: '$mc_min_dmesa 	>> $3
	echo 'mycond_max_'"$name_array[$i]"'_d_mesa: '$mc_max_dmesa 		>> $3
	echo 'mycond_min_'"$name_array[$i]"'_g_r: '$mc_min_g_r 			>> $3
	echo 'mycond_max_'"$name_array[$i]"'_g_r: '$mc_max_g_r			>> $3
	echo 'mycond_min_'"$name_array[$i]"'_i_lt_dmesa: '$mc_min_i_lt_dmesa 	>> $3
	echo 'mycond_max_'"$name_array[$i]"'_i_lt_dmesa: '$mc_max_i_lt_dmesa 	>> $3
endif
@ i++

if ("$plot_array[$i]" == 1) then
	echo "$name_array[$i]" 'plot_'"$name_array[$i]"'_config.txt '0 0 'False' 		>> $1
	echo '-: -' >> $3
endif
@ i++

if ("$plot_array[$i]" == 1) then
	echo "$name_array[$i]" 'plot_'"$name_array[$i]"'_config.txt '0 0 $obs_analyse 		>> $1
	echo '-: -' >> $3
endif
@ i++

#if ("$plot_array[$i]" == 1) then
	#echo "$name_array[$i]" 'plot_'"$name_array[$i]"'_config.txt '0 $data_offset_analyseTargetSelection $obs_analyse		>> $1
	#echo '-: -' >> $3
#endif
#@ i++

if ("$plot_array[$i]" == 1) then

	set k=1
	foreach item ($8)

		if ("$analyse_tarsel_array[$k]" == 1) then

			if ($item =~ mstar) then
				set obs = $obs_mstar
			else if ($item =~ mA* || $item =~ MA*) then
				set obs = $obs_mags
			else
				set obs = 'False'
			endif
			#echo "$name_array[$i]"'_'$item 'plot_'"$item"'_config.txt '0 $data_offset_analyseTargetSelection $obs_analyse
			echo "$name_array[$i]"'_'$item 'plot_'"$item"'_config.txt '0 $data_offset_analyseTargetSelection $obs_analyse >> $1
		endif
		@ k++	
	end
	echo 'mycond_min_'"$name_array[$i]"'_r_i: '$mc_min_r_i 			>> $3
	echo 'mycond_max_'"$name_array[$i]"'_r_i: '$mc_max_r_i 			>> $3
	echo 'mycond_min_mc_min_'"$name_array[$i]"'_d_mesa: '$mc_min_dmesa 	>> $3
	echo 'mycond_max_'"$name_array[$i]"'_d_mesa: '$mc_max_dmesa 		>> $3
	echo 'mycond_min_'"$name_array[$i]"'_g_r: '$mc_min_g_r 			>> $3
	echo 'mycond_max_'"$name_array[$i]"'_g_r: '$mc_max_g_r			>> $3
	echo 'mycond_min_'"$name_array[$i]"'_i_lt_dmesa: '$mc_min_i_lt_dmesa 	>> $3
	echo 'mycond_max_'"$name_array[$i]"'_i_lt_dmesa: '$mc_max_i_lt_dmesa 	>> $3
endif
@ i++

if ("$plot_array[$i]" == 1) then
	echo "$name_array[$i]" 'plot_'"$name_array[$i]"'_config.txt '$nbins_mstar $data_offset_histos $obs_mstar >> $1


	echo 'histo_min_SMF_mcold: '$h_min_SMF_mstar >> $3
	echo 'histo_max_SMF_mcold: '$h_max_SMF_mstar >> $3
	echo 'histo_min_SMF_mcold_disk: '$h_min_SMF_mstar >> $3
	echo 'histo_max_SMF_mcold_disk: '$h_max_SMF_mstar >> $3

	echo 'histo_min_SMF_mstar: '$h_min_SMF_mstar >> $3
	echo 'histo_max_SMF_mstar: '$h_max_SMF_mstar >> $3

	echo 'mycond_min_SMF_sfr: '$mc_min_SMF_mstar >> $3
	echo 'mycond_max_SMF_sfr: '$mc_max_SMF_mstar >> $3
	echo 'mycond_min_SMF_mcold: '$mc_min_SMF_mstar >> $3
	echo 'mycond_max_SMF_mcold: '$mc_max_SMF_mstar >> $3
	echo 'mycond_min_SMF_mcold_disk: '$mc_min_SMF_mstar >> $3
	echo 'mycond_max_SMF_mcold_disk: '$mc_max_SMF_mstar >> $3

	echo 'mycond_min_SMF_mstar: '$mc_min_SMF_mstar >> $3
	echo 'mycond_max_SMF_mstar: '$mc_max_SMF_mstar >> $3

endif
@ i++

if ("$plot_array[$i]" == 1) then
      	echo "$name_array[$i]"' plot_'"$name_array[$i]"'_config.txt '$nbins_sfr $data_offset_histos $obs_sfr>> $1

	echo 'histo_min_SFRF_sfr: '$h_min_SFRF_sfr >> $3
	echo 'histo_max_SFRF_sfr: '$h_max_SFRF_sfr >> $3
	echo 'mycond_min_SFRF_mstar: '$mc_min_SFRF_mstar >> $3
	echo 'mycond_max_SFRF_mstar: '$mc_max_SFRF_mstar >> $3
	echo 'mycond_min_SFRF_sfr: '$mc_min_SFRF_sfr >> $3
	echo 'mycond_max_SFRF_sfr: '$mc_max_SFRF_sfr >> $3
endif
@ i++

if ("$plot_array[$i]" == 1) then
      	echo "$name_array[$i]"' plot_'"$name_array[$i]"'_config.txt '$nbins_ssfr $data_offset_histos $obs_ssfr	>> $1

	echo 'histo_min_sSFRF_ssfr: '$h_min_sSFRF_ssfr >> $3
	echo 'histo_max_sSFRF_ssfr: '$h_max_sSFRF_ssfr >> $3
	echo 'mycond_min_sSFRF_mstar: '$mc_min_sSFRF_mstar >> $3
	echo 'mycond_max_sSFRF_mstar: '$mc_max_sSFRF_mstar >> $3
	echo 'mycond_min_sSFRF_ssfr: '$mc_min_sSFRF_ssfr >> $3
	echo 'mycond_max_sSFRF_ssfr: '$mc_max_sSFRF_ssfr >> $3
endif
@ i++

if ("$plot_array[$i]" == 1) then
      	echo "$name_array[$i]"' plot_'"$name_array[$i]"'_config.txt '$nbins_cSFRD $data_offset_cSFRD $obs_cSFRD>> $1

	echo 'histo_min_'"$name_array[$i]"'_mstar: '$h_min_sfr2z_mstar >> $3
	echo 'histo_max_'"$name_array[$i]"'_mstar: '$h_max_sfr2z_mstar >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mstar: '$mc_min_sfr2z_mstar >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mstar: '$mc_max_sfr2z_mstar >> $3
endif
@ i++

if ("$plot_array[$i]" == 1) then
      	echo "$name_array[$i]"' plot_'"$name_array[$i]"'_config.txt '$nbins_sfr2z $data_offset_sfr2z $obs_sfr2z>> $1

	echo 'histo_min_'"$name_array[$i]"'_mstar: '$h_min_sfr2z_mstar >> $3
	echo 'histo_max_'"$name_array[$i]"'_mstar: '$h_max_sfr2z_mstar >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mstar: '$mc_min_sfr2z_mstar >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mstar: '$mc_max_sfr2z_mstar >> $3
endif
@ i++

if ("$plot_array[$i]" == 1) then
      	echo "$name_array[$i]"' plot_'"$name_array[$i]"'_config.txt '$nbins_sfr2mstar $data_offset_bin $obs_sfr2mstar    >> $1

	echo 'histo_min_'"$name_array[$i]"': '$h_min_mstar2sfr >> $3
	echo 'histo_max_'"$name_array[$i]"': '$h_max_mstar2sfr >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mstar: '$mc_min_mstar2sfr_mstar >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mstar: '$mc_max_mstar2sfr_mstar >> $3
	echo 'mycond_min_'"$name_array[$i]"'_sfr: '$mc_min_mstar2sfr_sfr >> $3
	echo 'mycond_max_'"$name_array[$i]"'_sfr: '$mc_max_mstar2sfr_sfr >> $3
endif
@ i++

if ("$plot_array[$i]" == 1) then
      	echo "$name_array[$i]"' plot_'"$name_array[$i]"'_config.txt '$nbins_ssfr2mstar $data_offset_histos $obs_ssfr2mstar    >> $1

	echo 'histo_min_'"$name_array[$i]"'_mstar: '$h_min_mstar2ssfr >> $3
	echo 'histo_max_'"$name_array[$i]"'_mstar: '$h_max_mstar2ssfr >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mstar: '$mc_min_mstar2ssfr_mstar >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mstar: '$mc_max_mstar2ssfr_mstar >> $3
	echo 'mycond_min_'"$name_array[$i]"'_sfr: '$mc_min_mstar2ssfr_ssfr >> $3
	echo 'mycond_max_'"$name_array[$i]"'_sfr: '$mc_max_mstar2ssfr_ssfr >> $3
endif
@ i++

if ("$plot_array[$i]" == 1) then
      	echo "$name_array[$i]"' plot_'"$name_array[$i]"'_config.txt '$nbins_zgas2mstar $data_offset_bin $obs_oh2mstar    >> $1

	echo 'histo_min_'"$name_array[$i]"'_mstar: '$h_min_oh2mstar_mstar >> $3
	echo 'histo_max_'"$name_array[$i]"'_mstar: '$h_max_oh2mstar_mstar >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mstar: '$mc_min_oh2mstar_mstar >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mstar: '$mc_max_oh2mstar_mstar >> $3
	echo 'mycond_min_'"$name_array[$i]"'_sfr: '$mc_min_oh2mstar_sfr >> $3
	echo 'mycond_max_'"$name_array[$i]"'_sfr: '$mc_max_oh2mstar_sfr >> $3
endif
@ i++

if ("$plot_array[$i]" == 1) then
      	echo "$name_array[$i]"' plot_'"$name_array[$i]"'_config.txt '$nbins_mhalo $data_offset_histos $obs_mhalo   >> $1

	echo 'histo_min_'"$name_array[$i]"'_mhalo: '$h_min_mhalo_mhalo >> $3
	echo 'histo_max_'"$name_array[$i]"'_mhalo: '$h_max_mhalo_mhalo >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mhalo: '$mc_min_mhalo_mhalo >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mhalo: '$mc_max_mhalo_mhalo >> $3
	echo 'histo_min_'"$name_array[$i]"'_mhalo_cents: '$h_min_mhalo_mhalo >> $3
	echo 'histo_max_'"$name_array[$i]"'_mhalo_cents: '$h_max_mhalo_mhalo >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mhalo_cents: '$mc_min_mhalo_mhalo >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mhalo_cents: '$mc_max_mhalo_mhalo >> $3
	echo 'histo_min_'"$name_array[$i]"'_mhalo_200c: '$h_min_mhalo_mhalo >> $3
	echo 'histo_max_'"$name_array[$i]"'_mhalo_200c: '$h_max_mhalo_mhalo >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mhalo_200c: '$mc_min_mhalo_mhalo >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mhalo_200c: '$mc_max_mhalo_mhalo >> $3
	echo 'histo_min_'"$name_array[$i]"'_mhalo_cents_200c: '$h_min_mhalo_mhalo >> $3
	echo 'histo_max_'"$name_array[$i]"'_mhalo_cents_200c: '$h_max_mhalo_mhalo >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mhalo_cents_200c: '$mc_min_mhalo_mhalo >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mhalo_cents_200c: '$mc_max_mhalo_mhalo >> $3

endif
@ i++

if ("$plot_array[$i]" == 1 ) then
  	echo "$name_array[$i]"' plot_'"$name_array[$i]"'_config.txt '$nbins_mstar2mhalo $data_offset_bin $obs_mstar2mhalo    >> $1

	echo 'histo_min_'"$name_array[$i]"'_mhalo: '$h_min_mstar2mhalo_mhalo >> $3
	echo 'histo_max_'"$name_array[$i]"'_mhalo: '$h_max_mstar2mhalo_mhalo >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mhalo: '$mc_min_mstar2mhalo_mhalo >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mhalo: '$mc_max_mstar2mhalo_mhalo >> $3

endif
@ i++

if ("$plot_array[$i]" == 1) then
	echo "$name_array[$i]"' plot_'"$name_array[$i]"'_config.txt '$nbins_mstar2mhalofunc $data_offset_histos $obs_mstar2mhalofunc     >> $1

	echo 'histo_min_'"$name_array[$i]"'_mhalo: '$h_min_mh2mstfu_mhalo >> $3
	echo 'histo_max_'"$name_array[$i]"'_mhalo: '$h_max_mh2mstfu_mhalo >> $3
	echo 'histo_min_'"$name_array[$i]"'_mstar2mhalo: '$h_min_mh2mstfu_mstar2mhalo >> $3
	echo 'histo_max_'"$name_array[$i]"'_mstar2mhalo: '$h_max_mh2mstfu_mstar2mhalo >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mhalo: '$mc_min_mh2mstfu_mhalo >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mhalo: '$mc_max_mh2mstfu_mhalo >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mstar2mhalo: '$mc_min_mh2mstfu_mstar2mhalo >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mstar2mhalo: '$mc_max_mh2mstfu_mstar2mhalo >> $3
endif
@ i++

if ("$plot_array[$i]" == 1) then
	echo 'ngal2mhalo plot_ngal2mhalo_config.txt '$nbins_ngal2mhalo $data_offset_ngal $obs_mhalo2ngal >> $1

	echo 'histo_min_ngal2mhalo_mstar: '$h_min_ngal2mhalo_mstar >> $3
	echo 'histo_max_ngal2mhalo_mstar: '$h_max_ngal2mhalo_mstar >> $3
	echo 'mycond_min_ngal2mhalo_mstar: '$mc_min_ngal2mhalo_mstar >> $3
	echo 'mycond_max_ngal2mhalo_mstar: '$mc_max_ngal2mhalo_mstar >> $3
endif
@ i++

if ("$plot_array[$i]" == 1) then
      	echo "$name_array[$i]"' plot_'"$name_array[$i]"'_config.txt '$nbins_mhalo $data_offset_histos $obs_mhalo   >> $1

	echo 'histo_min_'"$name_array[$i]"'_mhalo: '$h_min_mhalo_no_mhalo >> $3
	echo 'histo_max_'"$name_array[$i]"'_mhalo: '$h_max_mhalo_no_mhalo >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mhalo: '$mc_min_mhalo_no_mhalo >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mhalo: '$mc_max_mhalo_no_mhalo >> $3
	echo 'histo_min_'"$name_array[$i]"'_mhalo_200c: '$h_min_mhalo_no_mhalo >> $3
	echo 'histo_max_'"$name_array[$i]"'_mhalo_200c: '$h_max_mhalo_no_mhalo >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mhalo_200c: '$mc_min_mhalo_no_mhalo >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mhalo_200c: '$mc_max_mhalo_no_mhalo >> $3

endif
@ i++

if ("$plot_array[$i]" == 1)  then
	echo "$name_array[$i]"' plot_'"$name_array[$i]"'_config.txt '$nbins_mstar2mhalo $data_offset_bin $obs_mstar2mhalo    >> $1

	echo 'histo_min_'"$name_array[$i]"'_mstar: '$h_min_mstar2mhalo_no_mstar >> $3
	echo 'histo_max_'"$name_array[$i]"'_mstar: '$h_max_mstar2mhalo_no_mstar >> $3
	echo 'histo_min_'"$name_array[$i]"'_mhalo: '$h_min_mstar2mhalo_no_mhalo >> $3
	echo 'histo_max_'"$name_array[$i]"'_mhalo: '$h_max_mstar2mhalo_no_mhalo >> $3
	echo 'histo_min_'"$name_array[$i]"'_mhalo_200c: '$h_min_mstar2mhalo_no_mhalo >> $3
	echo 'histo_max_'"$name_array[$i]"'_mhalo_200c: '$h_max_mstar2mhalo_no_mhalo >> $3
	echo 'histo_min_'"$name_array[$i]"'_mbasic: '$h_min_mstar2mhalo_no_mhalo >> $3
	echo 'histo_max_'"$name_array[$i]"'_mbasic: '$h_max_mstar2mhalo_no_mhalo >> $3
	echo 'histo_min_'"$name_array[$i]"'_mbasic_200c: '$h_min_mstar2mhalo_no_mhalo >> $3
	echo 'histo_max_'"$name_array[$i]"'_mbasic_200c: '$h_max_mstar2mhalo_no_mhalo >> $3
	echo 'histo_min_'"$name_array[$i]"'_mhalo_cents: '$h_min_mstar2mhalo_no_mhalo >> $3
	echo 'histo_max_'"$name_array[$i]"'_mhalo_cents: '$h_max_mstar2mhalo_no_mhalo >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mstar: '$mc_min_mstar2mhalo_no_mstar >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mstar: '$mc_max_mstar2mhalo_no_mstar >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mhalo: '$mc_min_mstar2mhalo_no_mhalo >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mhalo: '$mc_max_mstar2mhalo_no_mhalo >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mhalo_200c: '$mc_min_mstar2mhalo_no_mhalo >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mhalo_200c: '$mc_max_mstar2mhalo_no_mhalo >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mbasic: '$mc_min_mstar2mhalo_no_mhalo >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mbasic: '$mc_max_mstar2mhalo_no_mhalo >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mbasic_200c: '$mc_min_mstar2mhalo_no_mhalo >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mbasic_200c: '$mc_max_mstar2mhalo_no_mhalo >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mhalo_cents: '$mc_min_mstar2mhalo_no_mhalo >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mhalo_cents: '$mc_max_mstar2mhalo_no_mhalo >> $3

endif
@ i++

if ("$plot_array[$i]" == 1) then
	echo 'mstar2mhalofunc_no plot_mstar2mhalofunc_no_config.txt '$nbins_mstar2mhalofunc $data_offset_histos $obs_mstar2mhalofunc     >> $1

	echo 'histo_min_mstar2mhalofunc_no_mhalo: '$h_min_mh2mstfu_no_mhalo >> $3
	echo 'histo_max_mstar2mhalofunc_no_mhalo: '$h_max_mh2mstfu_no_mhalo >> $3
	echo 'histo_min_mstar2mhalofunc_no_mstar2mhalo: '$h_min_mh2mstfu_no_mstar2mhalo >> $3
	echo 'histo_max_mstar2mhalofunc_no_mstar2mhalo: '$h_max_mh2mstfu_no_mstar2mhalo >> $3
	echo 'mycond_min_mstar2mhalofunc_no_mhalo: '$mc_min_mh2mstfu_no_mhalo >> $3
	echo 'mycond_max_mstar2mhalofunc_no_mhalo: '$mc_max_mh2mstfu_no_mhalo >> $3
	echo 'mycond_min_mstar2mhalofunc_no_mstar2mhalo: '$mc_min_mh2mstfu_no_mstar2mhalo >> $3
	echo 'mycond_max_mstar2mhalofunc_no_mstar2mhalo: '$mc_max_mh2mstfu_no_mstar2mhalo >> $3

endif
@ i++

if ("$plot_array[$i]" == 1) then
	echo 'ngal2mhalo_no plot_ngal2mhalo_no_config.txt '$nbins_ngal2mhalo $data_offset_ngal $obs_mhalo2ngal >> $1

	echo 'histo_min_ngal2mhalo_no_mstar: '$h_min_ngal2mhalo_no_mstar >> $3
	echo 'histo_max_ngal2mhalo_no_mstar: '$h_max_ngal2mhalo_no_mstar >> $3
	echo 'mycond_min_ngal2mhalo_no_mstar: '$mc_min_ngal2mhalo_no_mstar >> $3
	echo 'mycond_max_ngal2mhalo_no_mstar: '$mc_max_ngal2mhalo_no_mstar >> $3
endif
@ i++

if ("$plot_array[$i]" == 1) then
	echo 'zgas2mstar plot_zgas2mstar_config.txt '$nbins_zgas2mstar $data_offset_bin $obs_zgas2mstar >> $1

	echo 'histo_min_zgas2mstar_mstar: '$h_min_zgas2mstar_mstar >> $3
	echo 'histo_max_zgas2mstar_mstar: '$h_max_zgas2mstar_mstar >> $3
	echo 'mycond_min_zgas2mstar_mstar: '$mc_min_zgas2mstar_mstar >> $3
	echo 'mycond_max_zgas2mstar_mstar: '$mc_max_zgas2mstar_mstar >> $3

endif
@ i++

if ("$plot_array[$i]" == 1) then
	echo 'zgas2mcold plot_zgas2mcold_config.txt '$nbins_zgas2mcold $data_offset_bin $obs_zgas2mcold >> $1

	echo 'histo_min_zgas2mcold_mcold: '$h_min_zgas2mstar_mstar >> $3
	echo 'histo_max_zgas2mcold_mcold: '$h_max_zgas2mstar_mstar >> $3
	echo 'mycond_min_zgas2mcold_mcold: '$mc_min_zgas2mstar_mstar >> $3
	echo 'mycond_max_zgas2mcold_mcold: '$mc_max_zgas2mstar_mstar >> $3

endif
@ i++

if ("$plot_array[$i]" == 1) then
	echo 'mcold2mstar plot_mcold2mstar_config.txt '$nbins_mcold2mstar $data_offset_bin $obs_mstar2mcold >> $1

	echo 'histo_min_mcold2mstar_mstar: '$h_min_mcold2mstar_mstar >> $3
	echo 'histo_max_mcold2mstar_mstar: '$h_max_mcold2mstar_mstar >> $3
	echo 'mycond_min_mcold2mstar_mstar: '$mc_min_mcold2mstar_mstar >> $3
	echo 'mycond_max_mcold2mstar_mstar: '$mc_max_mcold2mstar_mstar >> $3
	echo 'mycond_min_mcold2mstar_mcold: '$mc_min_mcold2mstar_mcold >> $3
	echo 'mycond_max_mcold2mstar_mcold: '$mc_max_mcold2mstar_mcold >> $3
endif
@ i++

if ("$plot_array[$i]" == 1) then
	echo 'mbh2mstarsph plot_mbh2mstarsph_config.txt '$nbins_mbh2mstarsph $data_offset_histos $obs_mbh2mstarsph >> $1

	echo 'histo_min_mbh2mstarsph_mstar_spheroid: '$h_min_mbh2mstarsph_mstar_spheroid >> $3
	echo 'histo_max_mbh2mstarsph_mstar_spheroid: '$h_max_mbh2mstarsph_mstar_spheroid >> $3
	echo 'mycond_min_mbh2mstarsph_mstar_spheroid: '$mc_min_mbh2mstarsph_mstar_spheroid >> $3
	echo 'mycond_max_mbh2mstarsph_mstar_spheroid: '$mc_max_mbh2mstarsph_mstar_spheroid >> $3
	echo 'mycond_min_mbh2mstarsph_mbh: '$mc_min_mbh2mstarsph_mbh >> $3
	echo 'mycond_max_mbh2mstarsph_mbh: '$mc_max_mbh2mstarsph_mbh >> $3
endif
@ i++

if ("$plot_array[$i]" == 1) then
	echo 'twoPCF plot_twoPCF_config.txt '$nbins_twoPCF $data_offset_histos $obs_twoPCF >> $1

	echo 'mycond_min_twoPCF_mstar: '$mc_min_twoPCF_mstar >> $3
	echo 'mycond_max_twoPCF_mstar: '$mc_max_twoPCF_mstar >> $3
	echo 'mycond_min_twoPCF_boxsize: '$mc_min_twoPCF_boxsize >> $3
	echo 'mycond_max_twoPCF_boxsize: '$mc_max_twoPCF_boxsize >> $3
	echo 'histo_min_twoPCF_galgalsep: '$h_min_twoPCF_galgalsep >> $3
	echo 'histo_max_twoPCF_galgalsep: '$h_max_twoPCF_galgalsep >> $3
	echo 'mycond_twoPCF_ngal: '$mc_twoPCF_ngal >> $3
endif
@ i++

if ("$plot_array[$i]" == 1) then
	echo 'HOD plot_HOD_config.txt '$nbins_sfr $data_offset_histos $obs_HOD >> $1

	echo 'mycond_min_HOD_mhalo: '$mc_min_HOD_mhalo >> $3
	echo 'mycond_max_HOD_mhalo: '$mc_max_HOD_mhalo >> $3
	echo 'histo_min_HOD_mhalo: '$h_min_HOD_mhalo >> $3
	echo 'histo_max_HOD_mhalo: '$h_max_HOD_mhalo >> $3

endif
@ i++

if ("$plot_array[$i]" == 1 ) then
  	echo "$name_array[$i]"' plot_'"$name_array[$i]"'_config.txt '$nbins_mstar2mhalo $data_offset_bin $obs_mstar2mhalovsSFR    >> $1

	echo 'histo_min_'"$name_array[$i]"'_sfr: '$h_min_mstar2mhalovssfr_sfr >> $3
	echo 'histo_max_'"$name_array[$i]"'_sfr: '$h_max_mstar2mhalovssfr_sfr >> $3
	echo 'mycond_min_'"$name_array[$i]"'_sfr: '$mc_min_mstar2mhalovssfr_sfr >> $3
	echo 'mycond_max_'"$name_array[$i]"'_sfr: '$mc_max_mstar2mhalovssfr_sfr >> $3
	echo 'histo_min_'"$name_array[$i]"'_ssfr: '$h_min_mstar2mhalovssfr_ssfr >> $3
	echo 'histo_max_'"$name_array[$i]"'_ssfr: '$h_max_mstar2mhalovssfr_ssfr >> $3
	echo 'mycond_min_'"$name_array[$i]"'_ssfr: '$mc_min_mstar2mhalovssfr_ssfr >> $3
	echo 'mycond_max_'"$name_array[$i]"'_ssfr: '$mc_max_mstar2mhalovssfr_ssfr >> $3
	echo 'histo_min_'"$name_array[$i]"'_mhalo: '$h_min_mstar2mhalovssfr_mhalo >> $3
	echo 'histo_max_'"$name_array[$i]"'_mhalo: '$h_max_mstar2mhalovssfr_mhalo >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mhalo: '$mc_min_mstar2mhalovssfr_mhalo >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mhalo: '$mc_max_mstar2mhalovssfr_mhalo >> $3
	echo 'histo_min_'"$name_array[$i]"'_mhalo_cents: '$h_min_mstar2mhalovssfr_mhalo >> $3
	echo 'histo_max_'"$name_array[$i]"'_mhalo_cents: '$h_max_mstar2mhalovssfr_mhalo >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mhalo_cents: '$mc_min_mstar2mhalovssfr_mhalo >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mhalo_cents: '$mc_max_mstar2mhalovssfr_mhalo >> $3
	echo 'histo_min_'"$name_array[$i]"'_mhalo_200c: '$h_min_mstar2mhalovssfr_mhalo >> $3
	echo 'histo_max_'"$name_array[$i]"'_mhalo_200c: '$h_max_mstar2mhalovssfr_mhalo >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mhalo_200c: '$mc_min_mstar2mhalovssfr_mhalo >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mhalo_200c: '$mc_max_mstar2mhalovssfr_mhalo >> $3
	echo 'histo_min_'"$name_array[$i]"'_mhalo_cents_200c: '$h_min_mstar2mhalovssfr_mhalo >> $3
	echo 'histo_max_'"$name_array[$i]"'_mhalo_cents_200c: '$h_max_mstar2mhalovssfr_mhalo >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mhalo_cents_200c: '$mc_min_mstar2mhalovssfr_mhalo >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mhalo_cents_200c: '$mc_max_mstar2mhalovssfr_mhalo >> $3
endif
@ i++
if ("$plot_array[$i]" == 1 ) then
  	echo "$name_array[$i]"' plot_'"$name_array[$i]"'_config.txt '25 $data_offset_bin $obs_mstar2rhalf    >> $1

	echo 'histo_min_'"$name_array[$i]"'_mstar: '$h_min_mstar2rhalf_mstar >> $3
	echo 'histo_max_'"$name_array[$i]"'_mstar: '$h_max_mstar2rhalf_mstar >> $3
	echo 'mycond_min_'"$name_array[$i]"'_mstar: '$mc_min_mstar2rhalf_mstar >> $3
	echo 'mycond_max_'"$name_array[$i]"'_mstar: '$mc_max_mstar2rhalf_mstar >> $3

	echo 'histo_min_'"$name_array[$i]"'_rhalf_mass: '$h_min_mstar2rhalf_rhalf >> $3
	echo 'histo_max_'"$name_array[$i]"'_rhalf_mass: '$h_max_mstar2rhalf_rhalf >> $3
	echo 'mycond_min_'"$name_array[$i]"'_rhalf_mass: '$mc_min_mstar2rhalf_rhalf >> $3
	echo 'mycond_max_'"$name_array[$i]"'_rhalf_mass: '$mc_max_mstar2rhalf_rhalf >> $3

endif
echo 'physical and statistical parameter config ... DONE!'
