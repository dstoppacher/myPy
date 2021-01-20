#!/bin/tcsh

######################################################################################################################
#
# HOW TO USE:
# ==========
# 
# - generate a soft-link to the models to be plotted in $PATH as specified below
#   -> those links need to direct to the folder that contains the galaxy catalogues
#
# - choose which snapshots to plot by setting $MYSNAPS
#   -> use the mapping between snapid and zred to pick the correct ones
#
# - choose which plot to actually make by setting plotconfig to the desired template
######################################################################################################################

echo ' '
echo 'DONT WORRY DORIS, YOU ARE DOING WELL TODAY ;-) '
echo ' '

#Set the name of the simulation: e.g. MDPL, Bolshoi, etc.!
#+++++++++++++++++++++++++++++++++++++++++++
set simulation_name = 'OBS'

#Set the name of the telescope: e.g. SDSS, etc.!
#+++++++++++++++++++++++++++++++++++++++++++
set telescope_name = 'SDSS'

#Set the name of the cosmology used in the simulation: e.g. Planck, etc.!
#+++++++++++++++++++++++++++++++++++++++++++
set cosmology = 'WMAP7'
#wMAP7: flat LCDM (H0=70, Omega_m=0.274, and Omega_L=0.726. (White et al. 2011))'

#Set any note to they output-file!
#+++++++++++++++++++++++++++++++++++++++++++
set my_annotation = ''

#Select with unit system should be used
#DS: the system standard units are chosen which are mass [Msun], sfr [Msun/yr], pos [comMpc], mean_age_stars [Gyr]
#MD: units according to Multidark data base mass [Msun/h], sfr [Msun/Gyr/h], pos [comMpc/h], mean_age_stars [Gyr]
#+++++++++++++++++++++++++++++++++++++++++++
set UNIT_CODE 	 = 'DS'

#Specify plot configurations: single: every catalog in a single plot, single_mix: mix of catalogs in one plot!
#+++++++++++++++++++++++++++++++++++++++++++
set plot_single = False
set plot_single_mix = False
set subplots = False


#Enter target selection code: e.g. 'CMASS', 'LRG', etc --> name of the subcat (hdf5 or ASCII) you will create with filterData method --> use the same tarsel_code if you want to filter an alrady selected catalog (subfilter)!
#+++++++++++++++++++++++++++++++++++++++++++
set add_tarsel_code = 'allmasses_spall'
set add_tarsel_code = '0.6_spall_'

#set add_tarsel_code = 'portsmouth'
#set add_tarsel_code = 'antonio'
#set add_tarsel_code = 'spall_'

set add_tarsel_code_end = '_DR12v4_compl_sample'
#set add_tarsel_code_end = '_DR12'

#set tarsel_code_default = 'allmasses-spall-complete'
#set tarsel_code_default = 'portsmouth-passive_kroup'

#set tarsel_code_default = $add_tarsel_code'spall_reduced'$add_tarsel_code_end
#set tarsel_code_default = $add_tarsel_code'por_PS'$add_tarsel_code_end
#set tarsel_code_default = $add_tarsel_code'por_SF'$add_tarsel_code_end
set tarsel_code_default = $add_tarsel_code'por_merged'$add_tarsel_code_end'_cross'
#set tarsel_code_default = $add_tarsel_code'gra_dust'$add_tarsel_code_end


#CMASS 0.55
#set tarsel_code_default = 'spall_wis_DR12v4_compl_sample'
#set tarsel_code_default = 'spall_gra_nodust_DR12v4_compl_sample'
#set tarsel_code_default = 'spall_por_merged_DR12v4_compl_sample'
#set tarsel_code_default = 'spall_por_SF_DR12v4_compl_sample'
#set tarsel_code_default = 'spall_por_PS_DR12v4_compl_sample'

#CMASS 0.43 - 0.7
#set tarsel_code_default = '0.7_wis_DR12v4_compl_sample'
#set tarsel_code_default = '0.7_gra_nodust_DR12v4_compl_sample'
#set tarsel_code_default = '0.7_por_merged_DR12v4_compl_sample'
#set tarsel_code_default = '0.7_por_SF_DR12v4_compl_sample'
#set tarsel_code_default = '0.7_por_PS_DR12v4_compl_sample'

#set tarsel_code_default = ''
#Choose different catalogues to read!
#+++++++++++++++++++++++++++++++++++++++++++

set add_tarsel_code = ''
#set add_tarsel_code = 'CMASS'
#set add_tarsel_code = 'densCut_CMASS_mstar'

set tarsel_code_Galacticus = ''$add_tarsel_code
set tarsel_code_SAG = 'v3_new_mags'$add_tarsel_code
set tarsel_code_SAGE = 'v3_1260'$add_tarsel_code

set output_filename_code_Galacticus = 'maras'
set output_filename_code_SAG = 'maras'
set output_filename_code_SAGE = 'maras'

#set tarsel_code_SAGE = 'v2_mags_350Mpc_run_1182'

#CHANGE IMF:
set change_IMF = 'Kroupa_Croton+16'
set change_IMF = 'False'
#give an additional identifier to the filename which is generated automatically --> use the same output_filename_code to read the histo-files in again!
#+++++++++++++++++++++++++++++++++++++++++++

#NOTE standard: [plot_key]_[catname+boxsize]_z_[redshift(0.2f)].[fileformat]
#set add_output_filename = 'z_-99.99_'
set add_output_filename = 'spall_'

set add_output_filename_end = '_DR12v4_compl_sample_wg_fixed_35bins_1e10_-0.2'
#set add_output_filename_end = '_DR12'

#set output_filename_code_default = $add_output_filename'por_merged'$add_output_filename_end
#set output_filename_code_default = $add_output_filename'por_PS_North'$add_output_filename_end
#set output_filename_code_default = $add_output_filename'por_SF'$add_output_filename_end
#set output_filename_code_default = $add_output_filename'gra_dust'$add_output_filename_end
#set output_filename_code_default = $add_output_filename'wis'$add_output_filename_end'
#set output_filename_code_default = 'Antonio_full'
set output_filename_code_default = 'wg_fixed_35bins_1e10_-0.2'

#name of the register you want to write/read data, if default register is not desired! Change 'default' to 'your name'
set set_register_name = ''
#++++++++++++++++++++++++++++++++++++++++++

#Set file job-id numbers which should be selected: NOTE: only Galacticus job0-job103==1Snapshot/redshift
#+++++++++++++++++++++++++++++++++++++++++++
set start_fileID 	= 0
set nr_files_snapshot 	= 1
set end_fileID 	= False

#choose a selection name which describes the current selection you want to apply
#+++++++++++++++++++++++++++++++++++++++++++
set selection_name = ''
#'CMASS 2PCF subcatalogue extraction for CUTE'

#Enter 'SAM' for reading a SAM catalouge or anything else to read a random HDF5 file
#+++++++++++++++++++++++++++++++++++++++++++

#NOTE: Load subcat=True to load the extracted subcatalogues of the SAM!
if ($HOME == /home/claudia || $HOME == /home/doriss) then
	set HOME_MODE = True
	set HDF5_read_code = 'HDF5'
	set load_subcat = True
	set PLOT_CUM = 'False'
else
	set HOME_MODE = False
	set HDF5_read_code = 'HDF5'
	if ($HDF5_read_code == 'HDF5') then
		set load_subcat = True
	else
		set load_subcat = False
	endif
endif

#If create_subcat == True, than a subsample of the catalog will be extracted, or the list of MYSNAPS will be worked trough like '0.0 0.1 0.5' (usefull e.g. with plotXY)
set create_subcat = True
set convert_sph_to_cart_coords = False
set PLOT_CUM = 'Truee'
set set_home_mode = Truee
set obs_volume = '1.7595e9'

if ($HOME_MODE == $set_home_mode) then
	set MANUAL_REDSHIFT_INPUT_default = 0.0
	set SOFTLINK_CODE = F
	set plot_all = 'True'			
	set LOAD_FROM_FILE = 'True'
	#set output_filename_code_default = 'tarsel'
else
	set MANUAL_REDSHIFT_INPUT_default = 0.5
	set SOFTLINK_CODE = CM
	set plot_all = 'False'		
	set LOAD_FROM_FILE = 'False'
endif

#Choose the Snapshots/redshifts --> see list for the SAMS below OR enter a key which is representative for the catalog e.g. 'Gal' for 'Galacticus_1Gpc.hdf5'
#+++++++++++++++++++++++++++++++++++++++++++
set MYSNAPS_62    			= '0.000'
set MYSNAPS_125				= '107' 
set MYSNAPS_400Mpc			= '96'

#redshifts Galacticus
#set MYSNAPS_1Gpc			= '79 75 70 67 65 60 57 52 50 47 45 42 40 38 36 34 32 30 28 26 24 22 20 18 16 14 13 12 11 10 9 8 7 6 5 4 3 2 1'


#redshifts SAGE
#set MYSNAPS_1Gpc			= '0.000 0.093 0.142 0.490 0.592 0.740 0.901 1.220 2.095 3.037 4.038 5.017 6.022 7.026 8.372'

#set MYSNAPS_1Gpc			= '125 124 122 121 119 115 110 107 104 101 98 95 92 89 86 83 80 77 75 70 67 65 60 57 52 50 47 45 42 40 38 36 34 32 30 28 26 24 20'
#set MYSNAPS_1Gpc			= '0.000 0.093 0.142 0.490 0.557 0.740 0.859 1.077 1.270 1.425 1.650 1.896 2.095 2.382 2.614 3.037 3.411 3.929 4.038 4.385 4.627 4.754 4.882 5.017 5.289 5.720 5.873 6.022 6.184 6.342 6.508 6.849 7.026 7.203 7.389 7.764 7.961 8.166 8.372'
#Wisconsin
#set MYSNAPS_default			= 'bc03-DR12'
#Granada no dust
#set MYSNAPS_default			= 'nodust'
#Granada dust
#set MYSNAPS_default			= '_dust'
#Portsmouth
#set MYSNAPS_default			= 'passive_krou-DR12'
#set MYSNAPS_default			= 'starforming_krou-DR12'
#set MYSNAPS_default			= 'merged'
#spall
set MYSNAPS_default			= 'reduced'
#default
set MYSNAPS_default			= 'dr12v4'
set MYSNAPS_default			= 'merged'
#set MYSNAPS_default			= 'full'
#for sfr2mstar plot: 0.0 0.6 0.75 0.9 1.2


#Galacticus SMPL 400 Mpc/h: 96 z=0.0, 92 z=0.1, 85 z=0.14, 47 z=0.55'
#Galacticus MDPL2 1 Gpc: 79 z=0, 75 z=0.09, 73 z=0.14, 61 z=0.49, 59 z=0.56, 55=0.7, 52=.82, 50 =0.9, 29 z=2.03, 24 z=2.38, min 1 z=8.0
#SAG MDPL2 1 Gpc: 	125 z=0,  122 z=0.07, 119 z=0.14, 107=0.49, 104 z=0.59, 101=0.7, min 75 z=2.03
#SAG cal01 nifty 125: 	107 z=0,  98 z=0.15, min 4 z=7.73

if ($HDF5_read_code == 'SAMHDF5' || $HDF5_read_code == 'BINARY_SAGE' || $HDF5_read_code == 'FITS2HDF5') then

	if ($SOFTLINK_CODE == A) then 
		set MYSNAPS_1Gpc	= '75'
		set MYSNAPS_1Gpc	= '79 75 70 67 65 60 57 52 50 47 45 42 40'
# 38 36 34 32 30 28 26 24 22 20 18 16 14 13 12 11 10 9 8 7 6 5 4 3 2 1'

		#Galacticus 400 Mpc/h snapshot list for sfr2z
		#'96 94 92 90 87 85 84 80 75 70 65 50 47 45 42 40 38 36 30 28 26 24 22 20 18 16 14 12 10 8 6 5 4 3 2 1'
	else if ($SOFTLINK_CODE == B || $SOFTLINK_CODE == I) then
		#set MYSNAPS_1Gpc	= '125 92 89 86 83 80 77 75'

		#redshifts  SAGv2 
		set MYSNAPS_1Gpc	= '125'
# 122 119 116 113 110 107 104 101 98 95 92 89 86 83 80 77 75'

		#redshifts  SAGv3
		#set MYSNAPS_default	= '125 124 122 121 119 115 110 107 104 101 98 95 92 89 86 83 80 77 75 70 67 65 60 57 52 50 47 45 42 40 38 36 34 32 30 28 26 24 20'
	
		#SAG calibration01 nifty 125 Mpc/h snapshot list for sfr2z
		#107 98 92 80 70 61 47 35 27 21 17 14 12 9 6 4

	else

		set MYSNAPS_1Gpc	= '0.000'
		#set MYSNAPS_1Gpc	= '0.000 0.093 0.142 0.490 0.557 0.740 0.859 1.077 1.270 1.425 1.650 1.896 2.095 2.382 2.614 3.037 3.411 3.929 4.038 4.385 4.627 4.754 4.882 5.017 5.289 5.720 5.873 6.022 6.184 6.342 6.508 6.849 7.026 7.203 7.389 7.764 7.961 8.166 8.372'

	endif
else
	if ($SOFTLINK_CODE == A) then 
		set MYSNAPS_1Gpc	= '59'
	else if ($SOFTLINK_CODE == B) then 
		set MYSNAPS_1Gpc	= '105'
	else if ($SOFTLINK_CODE == I) then 
		set MYSNAPS_1Gpc	= '95'
	else
		set MYSNAPS_1Gpc	= '0.54'
	endif	
endif

#Granada DR12 BOSS nodust
#set fits_map_names	= 'DEC;RA;Z;LOGMASS;LOGMASS_MEDIAN;LOGMASS_ERR;AGE;SSFR;SSFR_MIN;SSFR_MAX;MAGSCALED;MAGSCALED1;MAGSCALED2;MAGSCALED3;MAGSCALED4;CMODELABSMAG;CMODELABSMAG1;CMODELABSMAG2;CMODELABSMAG3;CMODELABSMAG4;KE;KE1;KE2;KE3;KE4;BOSS_TARGET1;SPECPRIMARY;WARNING'
#set fits_map_names_mapping = 'DEC;RA;Z;mstar;mstar_5000;mstar_err;age;ssfr;ssfr_min;ssfr_max;mAB_total_u;mAB_total_g;mAB_total_r;mAB_total_i;mAB_total_z;cMAB_total_u;cMAB_total_g;cMAB_total_r;cMAB_total_i;cMAB_total_z;kcorr_u;kcorr_g;kcorr_r;kcorr_i;kcorr_z;boss_target1;specprimary;warning'
#Granada DR12 BOSS dust
#set fits_map_names_mapping = 'DEC;RA;Z;mstar;mstar_5000;mstar_err;age;ssfr;ssfr_min;ssfr_max;mAB_dA_total_u;mAB_dA_total_g;mAB_dA_total_r;mAB_dA_total_i;mAB_dA_total_z;cMAB_dA_total_u;cMAB_dA_total_g;cMAB_dA_total_r;cMAB_dA_total_i;cMAB_dA_total_z;kcorr_u;kcorr_g;kcorr_r;kcorr_i;kcorr_z;boss_target1;specprimary;warning'

#Portsmouth 
#set fits_map_names	= 'PLATE;FIBERID;DEC;RA;Z;Z_ERR;LOGMASS;MINLOGMASS;MAXLOGMASS;AGE;SFR;MINSFR;MAXSFR;MAGSCALED;MAGSCALED1;MAGSCALED2;MAGSCALED3;MAGSCALED4;MAGSCALED_ERR;MAGSCALED_ERR1;MAGSCALED_ERR2;MAGSCALED_ERR3;MAGSCALED_ERR4;ABSMAGK;BOSS_TARGET1;SPECPRIMARY'
#set fits_map_names_mapping = 'haloid;hostid;DEC;RA;Z;Z_err;mstar;mstar_min;mstar_max;age;sfr;sfr_min;sfr_max;mAB_dA_total_u;mAB_dA_total_g;mAB_dA_total_r;mAB_dA_total_i;mAB_dA_total_z;mAB_dA_total_err_u;mAB_dA_total_err_g;mAB_dA_total_err_r;mAB_dA_total_err_i;mAB_dA_total_err_z;Mag_dA_total_K;boss_target1;specprimary'

#set fits_map_names	= 'DEC_1;RA_1;Z_1;LOGMASS;WEIGHT_FKP;WEIGHT_SEEING;WEIGHT_CP;WEIGHT_STAR;WEIGHT_SYSTOT;WEIGHT_NOZ;MAGSCALED;MAGSCALED1;MAGSCALED2;MAGSCALED3;MAGSCALED4'
#set fits_map_names_mapping = 'DEC;RA;Z;mstar;weight_fkp;weight_seeing;weight_cp;weight_star;weight_systot;weight_noz;mAB_dA_total_u;mAB_dA_total_g;mAB_dA_total_r;mAB_dA_total_i;mAB_dA_total_z'

#Wisconsin DR12 BOSS
##set fits_map_names	= 'DEC;RA;Z;MSTELLAR_MEDIAN;MSTELLAR_P97P5;MSTELLAR_P84;MSTELLAR_P16;MSTELLAR_P2P5;MSTELLAR_ERR'
#set fits_map_names_mapping = 'DEC;RA;Z;mstar;mstar_9750;mstar_8400;mstar_0160;mstar_0025;mstar_err'

#DR13v4 spall from Sergio
#set fits_map_names	= 'DEC;RA;Z;LOGMASS;WEIGHT_FKP;WEIGHT_SEEING;WEIGHT_CP;WEIGHT_STAR;WEIGHT_SYSTOT;WEIGHT_NOZ;AGE;U;G;R;I'
#set fits_map_names_mapping = 'DEC;RA;Z;mstar;weight_fkp;weight_seeing;weight_cp;weight_star;weight_systot;weight_noz;age;mAB_dA_total_u;mAB_dA_total_g;mAB_dA_total_r;mAB_dA_total_i'
set fits_map_names	= 'PLATE;FIBERID;DEC;RA;Z;LOGMASS_PORT_MER;WEIGHT_FKP;WEIGHT_SEEING;WEIGHT_CP;WEIGHT_STAR;WEIGHT_SYSTOT;WEIGHT_NOZ;AGE;U;G;R;I'
set fits_map_names_mapping = 'haloid;hostid;DEC;RA;Z;mstar;weight_fkp;weight_seeing;weight_cp;weight_star;weight_systot;weight_noz;age;mAB_dA_total_u;mAB_dA_total_g;mAB_dA_total_r;mAB_dA_total_i'

#SPALL from Antonio
set fits_map_names	= 'DEC;RA;Z;U;G;R;I'
#;PLATE;FIBER_ID;MJD;THETA;THETA_ERR;VDISP;VDISP_ERR;VDISP_SHU_MEAN;VDISP_SHU_ERR'
set fits_map_names_mapping = 'DEC;RA;Z;mAB_dA_total_u;mAB_dA_total_g;mAB_dA_total_r;mAB_dA_total_i'


#SDSS DR7 MPA-JHU
#set fits_map_names	= '0;1;2;3'
#set fits_map_names_mapping = 'DEC;RA;Z;mstar;'

#Observational catalog details
#CMASS DR12
set skycoverage = 9376
set zmin	= 0.5
set zmax	= 0.6
set col_name_weights	= 'weight_tot'

#Give infos about the sample of galaxies selected: 'full': no cuts what so ever, or standard cuts. 'else': a subsample will be selected --> choose number of galaxies which should be selected!
#+++++++++++++++++++++++++++++++++++++++++++
set sample_info	= 'full'

#Select ngalaxies for random catalogue:
set ngalaxies_random = 500000

#Select 'skip_reading_data' if for example with 2PCF the hdf5 should not be read or is not exsiting!
#+++++++++++++++++++++++++++++++++++++++++++
set skip_reading_data = yees
set use_store_register 	= no
set filter_density = 'False'
#Contrears CUTS
#set filter_density_ngal = 46750000,11770000,530000
#set filter_density_cut_name = CUT1_Contreras+13,CUT2_Contreras+13,CUT3_Contreras+13

set filter_density_ngal = 340000
#3.4x10-4 Mpc-3 h3
set filter_density_cut_name = densCut_CMASS

set filter_density_select_highest = True
set filter_density_write_coordinates = False
set ASCII_TO_HDF5 = True
set HDF5_TO_ASCII = False

set which_kcorrect = 'approx'

#specify details for 'twoPCF'!
#+++++++++++++++++++++++++++++++++++++++++++
set twoPCF_path_to_data	 = '/store/erebos/doris/'
set twoPCF_which	 = 'BOX'
set twoPCF_calculate	 = 'wp'
set twoPCF_pimax 	 = '100'
set twoPCF_nthreads 	 = '8'

set BAO = 'Truee'

if ($BAO == 'True') then
	set output_filename_code_Galacticus = 'maras_BAO_peak3'
	set output_filename_code_SAG = 'maras_BAO_peak3'
	set output_filename_code_SAGE = 'maras_BAO_peak'
	set twoPCF_rmin 	 = '120'
	set twoPCF_rmax		 = '200'
	set twoPCF_nbins	 = '40'

else
	set twoPCF_rmin 	 = '0.01'
	set twoPCF_rmax		 = '200'
	set twoPCF_nbins	 = '60'
endif



set twoPCF_ngal_random_sample = 'False'
set filter_halomass_Sergio = 'False'

#specify details for 'selectRegion'!
#+++++++++++++++++++++++++++++++++++++++++++

set selectRegion 		 = 'False'
#enter 'SPHERE' to calculate the volume of a sphere and anything for a BOX around the coordinates
set selectRegion_volume_type	 = 'SPHERE'
set selectRegion_filename 	 = 'regions.txt'

if ($HOME == /home/claudia || $HOME == /home/doris) then
	set selectRegion_path_to_data 	 = $HOME'/anaconda/pro/data/300/'
	set selectRegion_output_filename = $selectRegion_path_to_data 
else
	set selectRegion_path_to_data 	 = '/store/erebos/doris/workshop_300/'
	set selectRegion_output_filename = '/store/erebos/doris/workshop_300/'
endif
set selectRegion_periodic_conds	 = 'False'
set selectRegion_units	 	 = 'h-1Mpc'
set selectRegion_region_name 	 = 0
set selectRegion_col_id_x_pos 	 = 1
set selectRegion_col_id_y_pos 	 = 2
set selectRegion_col_id_z_pos 	 = 3
set selectRegion_col_id_radius 	 = 4

#specify details for 'createProgenitorList'!
#+++++++++++++++++++++++++++++++++++++++++++

set createProgenitorList		 = 'False'
#enter 'SPHERE' to calculate the volume of a sphere and anything for a BOX around the coordinates

if ($HOME == /home/claudia || $HOME == /home/doris) then
	set createProgenitorList_path_to_data 	 = $HOME'/anaconda/pro/data/sussing_trees_link/'
	set createProgenitorList_output_filename = $HOME'/anaconda/pro/data/sussing_trees/'
else
	set createProgenitorList_path_to_data 	 = '/store/erebos/cnvega/'
	set createProgenitorList_output_filename = '/store/erebos/doris/progenitorList/'
endif


#Choose the method you want to compute/plot:
#+++++++++++++++++++++++++++++++++++++++++++

set plotXY			= 0
set plotOnly			= 0
set filterData			= 0
set matchHaloCat		= 0
set analyseTargetSelection	= 0

set SMF 			= 0
set SFRF 			= 0
set sSFRF			= 1
set sfr2z			= 0
set sfr2mstar			= 0
set ssfr2mstar			= 0
set oh2mstar			= 0
set HMF 			= 0
set mstar2mhalo			= 0
set mstar2mhalofunc		= 0
set ngal2mhalo 			= 0
set HMF_no			= 0
set mstar2mhalo_no		= 0
set mstar2mhalofunc_no		= 0
set mstar2mhalovsSFR		= 0
set ngal2mhalo_no		= 0
set zgas2mstar			= 0
set zgas2mcold			= 0
set mcold2mstar			= 0
set mbh2mstarsph		= 0
#CUTE files form Alexander doris@taurus:/home/aknebe/Projects/MultiDarkGalaxies/CUTE/2PCF/
set twoPCF			= 0
set HOD				= 0

#Set manually redshift if home mode or load subcat are activated:
#+++++++++++++++++++++++++++++++++++++++++++

set histo_norm_y = False

#Default settings:
if ($HOME_MODE == $set_home_mode || $load_subcat == True) then

	if ($SMF == 1 || $HMF == 1 || $HMF_no == 1) then
		set SAG_REDSHIFT_INPUT = 0.56
		set SAG_125_REDSHIFT_INPUT = 1.22
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 0.56
		set GALACTICUS_400Mpc_REDSHIFT_INPUT = 0.1
		set SAGE_REDSHIFT_INPUT = 0.56

	else if ($SFRF == 1) then
		set SAG_125_REDSHIFT_INPUT = 0.15
		set SAG_REDSHIFT_INPUT = 1.22
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 1.22
		set GALACTICUS_400Mpc_REDSHIFT_INPUT = 0.14
		set SAGE_REDSHIFT_INPUT = 0.14

	else if ($sSFRF == 1) then
		set SAG_125_REDSHIFT_INPUT = 0.0
		set SAG_REDSHIFT_INPUT = 0.0
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 0.0
		set GALACTICUS_400Mpc_REDSHIFT_INPUT = 0.0
		set SAGE_REDSHIFT_INPUT = 0.0

	else if ($sfr2z == 1) then
		set SAG_REDSHIFT_INPUT = 8.58
		set SAG_125_REDSHIFT_INPUT = 7.37
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 8.0
		set SAGE_REDSHIFT_INPUT = 8.37

	else if ($ssfr2mstar == 1) then
		set SAG_125_REDSHIFT_INPUT = 0.0
		set SAG_REDSHIFT_INPUT = 0.0
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 0.0
		set GALACTICUS_400Mpc_REDSHIFT_INPUT = 0.0
		set SAGE_REDSHIFT_INPUT = 0.0

	else if ($sfr2mstar == 1) then
		set SAG_125_REDSHIFT_INPUT = 0.0
		set SAG_REDSHIFT_INPUT = 0.59
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 0.0
		set GALACTICUS_400Mpc_REDSHIFT_INPUT = 0.0
		set SAGE_REDSHIFT_INPUT = 0.0

	else if ($oh2mstar == 1) then
		set SAG_125_REDSHIFT_INPUT = 0.0
		set SAG_REDSHIFT_INPUT = 0.09
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 0.09
		set GALACTICUS_400Mpc_REDSHIFT_INPUT = 0.0
		set SAGE_REDSHIFT_INPUT = 0.09

	else if ($zgas2mstar == 1 || $zgas2mcold == 1) then
		set SAG_125_REDSHIFT_INPUT = 0.0
		set SAG_REDSHIFT_INPUT = 0.09
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 0.09
		set GALACTICUS_400Mpc_REDSHIFT_INPUT = 0.0
		set SAGE_REDSHIFT_INPUT = 0.09

	else if ($mcold2mstar == 1) then
		set SAG_125_REDSHIFT_INPUT = 0.0
		set SAG_REDSHIFT_INPUT = 0.0
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 0.0
		set GALACTICUS_400Mpc_REDSHIFT_INPUT = 0.0
		set SAGE_REDSHIFT_INPUT = 0.0

	else if ($mbh2mstarsph == 1) then
		set SAG_125_REDSHIFT_INPUT = 0.0
		set SAG_REDSHIFT_INPUT = 0.0
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 0.0
		set GALACTICUS_400Mpc_REDSHIFT_INPUT = 0.0
		set SAGE_REDSHIFT_INPUT = 0.0

	else if ($mstar2mhalo == 1 || $mstar2mhalo_no == 1 || $mstar2mhalovsSFR == 1 ) then
		set SAG_REDSHIFT_INPUT = 0.09
		set SAG_125_REDSHIFT_INPUT = 0.51
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 0.09
		set GALACTICUS_400Mpc_REDSHIFT_INPUT = 0.55
		set SAGE_REDSHIFT_INPUT = 0.09

	else if ($filterData == 1 ) then
		set SAG_125_REDSHIFT_INPUT = 0.0
		set SAG_REDSHIFT_INPUT = 0.09
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 0.09
		set SAGE_REDSHIFT_INPUT = 0.56


	
	else if ($analyseTargetSelection == 1 ) then
		#OII-sfr paper: 0.0, 0.07, 0.59, 0.7, 0.94, 1.22
		set SAG_REDSHIFT_INPUT = 0.56
		set SAG_125_REDSHIFT_INPUT = 0.0
		#OII-sfr paper: 0.0, 0.09, 0.59, 0.74, 0.9, 1.22
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 0.0
		set SAGE_REDSHIFT_INPUT = 0.0


	else if ($plotOnly == 1 ) then
		set SAG_REDSHIFT_INPUT = 0.0
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 0.56
		set SAGE_REDSHIFT_INPUT = 0.0
	
	else if ($twoPCF == 1 ) then
		set SAG_REDSHIFT_INPUT = 0.0
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 0.09
		set SAGE_REDSHIFT_INPUT = 0.00

	else

		set SAG_125_REDSHIFT_INPUT = False
		set SAG_REDSHIFT_INPUT = 0.0
		set SAGE_REDSHIFT_INPUT = 0.56
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 0.09
		set GALACTICUS_400Mpc_REDSHIFT_INPUT = False


	endif
else

	set SAG_125_REDSHIFT_INPUT = False
	set SAG_REDSHIFT_INPUT = False
	set SAGE_REDSHIFT_INPUT = 0.56
	set GALACTICUS_1Gpc_REDSHIFT_INPUT = False
	set GALACTICUS_400Mpc_REDSHIFT_INPUT = False

endif

#Set method specific option:
#+++++++++++++++++++++++++++++++++++++++++++
#filterData
if ($filterData == 1) then
	set create_subcat = True
endif


#2PCF:
#choose '_box' for CUTE_box or default for normal CUTE!
set choose_CUTE = _box
set CF_NR_bins = 128
#unit = MPC/h
set CF_R_max = 200
set corr_estimator = monopol
set param_file_CUTE = 'myparams'$choose_CUTE'.txt'
set CF_box_size = 1475.58

#analyseTargetSelection: Specify the plot you want to make
if ($analyseTargetSelection == 1) then
	set filterSample 		= 'False'
	set calcHistoSample 		= 'False'
	if ($calcHistoSample == 'True') then
		set analyse_tarsel_histo	= 1
	else
		set analyse_tarsel_histo	= 0

	endif

	set plotSample			= 'True'
else
	set filterSample 		= 'False'
	set calcHistoSample 		= 'False'
	set analyse_tarsel_histo	= 0
	set plotSample			= 'False'
endif

set analyse_tarsel_dmesa_i	= 0
set analyse_tarsel_gminr_mstar 	= 0
set analyse_tarsel_rmini_mstar	= 0
set analyse_tarsel_gmini_mstar	= 0
set analyse_tarsel_uminr_mstar	= 0
set analyse_tarsel_gminz_mstar	= 0
set analyse_tarsel_rmini_i	= 0
set analyse_tarsel_i_mstar	= 0
set analyse_tarsel_r_mstar	= 0
set analyse_tarsel_g_mstar	= 0

#Absolute Magnitude are signed with an I before: e.g. Ii=MAB_i band
set analyse_tarsel_Ii_i		= 0
set analyse_tarsel_Ii_mstar	= 0
set analyse_tarsel_Ir_r		= 0
set analyse_tarsel_Ir_mstar	= 0
set analyse_tarsel_Ig_g		= 0
set analyse_tarsel_Ig_mstar	= 0
set analyse_tarsel_Irmini_i	= 0

set analyse_tarsel_sfr_mstar	= 0
set analyse_tarsel_mbh_mstar_spheroid	= 0
set analyse_tarsel_mcold_mstar	= 0
set analyse_tarsel_Mzgas_mstar	= 0
set analyse_tarsel_Mzgas_mcold	= 0
set analyse_tarsel_gminr_r	= 0
set analyse_tarsel_uminr_r	= 0

set analyse_tarsel_gminr_uming	= 0
set analyse_tarsel_rmini_gminr	= 0

#Choose variables to compute/plot/process:
#+++++++++++++++++++++++++++++++++++++++++++

#99='DO NOT CHOOSE COLUMN', others=col_id
#COLUMN-CODE IN HALO CATALOGUE
set halo_nhalos		= 99
set halo_haloid		= 99
set halo_hostid		= 99
set halo_desIndex	= 99
set halo_nodemass	= 99
set halo_vmax		= 99
set halo_rsca		= 99
set halo_angmom		= 99
set halo_spin		= 99

set halo_redshift	= 99

set halo_x_pos	= 99
set halo_y_pos	= 99
set halo_z_pos	= 99

set halo_x_vel	= 99
set halo_y_vel	= 99
set halo_z_vel	= 99

set halo_x_vel_disp	= 99
set halo_y_vel_disp	= 99
set halo_z_vel_disp	= 99

#Standard properties
set ngalaxies		= 99
set haloid		= 0				
set hostid		= 1				
set orphan		= 99				
set mhalo		= 99

#only for Galacticus
set mhalo_sat		= 99
set nodeIndex		= 99
set parentIndex		= 99
set satelliteNodeIndex	= 99
set satelliteIndex	= 99
set siblingIndex	= 99
set satelliteMergeTime	= 99
set isolated		= 99
set timeLastIsolated	= 99
set firstProgenitorID	= 99
set npros		= 99

#only for SAG and SAGE
set vmax		= 99

#only for SAG
set vpeak		= 99
set NFW_con		= 99

#only SAG
set rhalf_bulge		= 99
set rhalf_disk		= 99

#only SAGE and Galacticus
set rdisk		= 99

#only Galacticus
set rbulge		= 99
set rhalf_mass		= 99

#only SAG and SAGE
set mean_age_stars	= 99

#only Galacticus and SAGE		
set spinParameter	= 99

set mbh			= 99

set mstar_spheroid	= 99
set mstar_disk		= 99
set mstar		= 2

set mstar_IC		= 99

set mcold_disk		= 99

#only Galacitcus and SAG
set mcold_spheroid	= 99
set mcold		= 99

set mhot		= 99

set sfr_spheroid	= 99
set sfr_disk		= 99
set sfr			= 20

#only SAG
set sfr_spheroid_inst	= 99
set sfr_quies_inst	= 99


set Mzgas_disk		= 99

#only SAG
set Mzgas_spheroid	= 99
set Mzgas		= 99


#only SAG and SAGE
set Mzstar_spheroid	= 99
set Mzstar_disk		= 99
set Mzstar		= 99
set Mzhot_halo		= 99

#only Galacticus
set zgas_spheroid	= 99

#only Galacticus and SAGE
set zgas_disk		= 99

#only Galacticus			 
set zstar_spheroid	= 99
set zstar_disk		= 99
set zhot_halo		= 99

#only SAG
set OH_gas_disk		= 99
set OH_gas_bulge	= 99
set OH_gas_disk_bulge	= 99

set x_pos		= 99
set y_pos		= 99
set z_pos		= 99
set x_vel		= 99
set y_vel		= 99
set z_vel		= 99


	
#Observational data
set DEC			= 3
set RA			= 4
set Z			= 5
set Z_err		= 99
set age			= 99
set weight		= 99
set weight_fkp		= 99
set weight_seeing	= 99
set weight_star		= 99
set weight_cp		= 99
set weight_noz		= 99
set weight_systot	= 99
set weight_tot		= 99

set mstar_err		= 99
set mstar_9750		= 99
set mstar_8400		= 99
set mstar_5000		= 99
set mstar_0160		= 99
set mstar_0025		= 99
set mstar_max		= 99
set mstar_min		= 99

set mstar_char		= 99

set sfr_min		= 99
set sfr_max		= 99

set ssfr		= 21
set ssfr_min		= 99
set ssfr_max		= 99

set kcorr_u		= 99
set kcorr_g		= 99
set kcorr_r		= 99
set kcorr_i		= 99
set kcorr_z		= 99

set mAB_dA_total_u	= 99
set mAB_dA_total_g	= 99
set mAB_dA_total_r	= 99
set mAB_dA_total_i	= 99
set mAB_dA_total_z	= 99

set MAB_dA_total_u	= 99
set MAB_dA_total_g	= 99
set MAB_dA_total_r	= 99
set MAB_dA_total_i	= 99
set MAB_dA_total_z	= 99

set Mag_dA_total_K	= 99

set mAB_dA_total_err_u	= 99
set mAB_dA_total_err_g	= 99
set mAB_dA_total_err_r	= 99
set mAB_dA_total_err_i	= 99
set mAB_dA_total_err_z	= 99

set mAB_total_u		= 99
set mAB_total_g		= 99
set mAB_total_r		= 99
set mAB_total_i		= 99
set mAB_total_z		= 99

set MAB_total_u	= 99
set MAB_total_g	= 99
set MAB_total_r	= 99
set MAB_total_i	= 99
set MAB_total_z	= 99

set mAB_total_cut_r_i	= 99
set mAB_total_cut_dmesa	= 99
set mAB_total_cut_i_lt_dmesa	= 99

set mAB_dA_total_cut_r_i	= 99
set mAB_dA_total_cut_dmesa	= 99
set mAB_dA_total_cut_i_lt_dmesa	= 99

set cmAB_dA_total_u	= 99
set cmAB_dA_total_g	= 99
set cmAB_dA_total_r	= 99
set cmAB_dA_total_i	= 99
set cmAB_dA_total_z	= 99

set cmAB_total_u	= 99
set cmAB_total_g	= 99
set cmAB_total_r	= 99
set cmAB_total_i	= 99
set cmAB_total_z	= 99

set cMAB_total_u	= 99
set cMAB_total_g	= 99
set cMAB_total_r	= 99
set cMAB_total_i	= 99
set cMAB_total_z	= 99

set cMAB_dA_total_u	= 99
set cMAB_dA_total_g	= 99
set cMAB_dA_total_r	= 99
set cMAB_dA_total_i	= 99
set cMAB_dA_total_z	= 99

set boss_target1	= 99
set specprimary		= 99
set warning		= 99

if ($tarsel_code_default =~ *'compl'*) then

	set mAB_dA_total_u	= 6
	set mAB_dA_total_g	= 7
	set mAB_dA_total_r	= 8
	set mAB_dA_total_i	= 9

	set weight_fkp		= 10
	set weight_seeing	= 11
	set weight_star		= 12
	set weight_cp		= 13
	set weight_noz		= 14
	set weight_systot	= 15
	set weight_tot		= 16

	set mAB_dA_total_cut_r_i	= 17
	set mAB_dA_total_cut_dmesa	= 18
	set mAB_dA_total_cut_i_lt_dmesa	= 19

else if ($tarsel_code_default =~ *'wis'*) then

	#Observational data
	set mstar_err		= 9
	set mstar_9750		= 5
	set mstar_8400		= 6
	set mstar_0160		= 7
	set mstar_0025		= 8

else if ($tarsel_code_default =~ *'gra'*) then
	#Observational data
	set mstar_err		= 5
	set mstar_5000		= 6
	set ssfr		= 7
	set ssfr_min		= 8
	set ssfr_max		= 9
	set kcorr_u		= 10
	set kcorr_g		= 11
	set kcorr_r		= 12
	set kcorr_i		= 13
	set kcorr_z		= 14
	
	set Z_err		= 28
	set boss_target1	= 29
	set specprimary		= 30
	set warning		= 31

	if ($tarsel_code_default =~ *'nodust'*) then

		set mAB_total_cut_r_i	= 15
		set mAB_total_cut_dmesa	= 16
		set mAB_total_cut_i_lt_dmesa	= 17

		set mAB_total_u	= 18
		set mAB_total_g	= 19
		set mAB_total_r	= 20
		set mAB_total_i	= 21
		set mAB_total_z	= 22

		set cMAB_total_u	= 23
		set cMAB_total_g	= 24
		set cMAB_total_r	= 25
		set cMAB_total_i	= 26
		set cMAB_total_z	= 27


	else 

		set mAB_dA_total_cut_r_i	= 15
		set mAB_dA_total_cut_dmesa	= 16
		set mAB_dA_total_cut_i_lt_dmesa	= 17

		set mAB_dA_total_u	= 18
		set mAB_dA_total_g	= 19
		set mAB_dA_total_r	= 20
		set mAB_dA_total_i	= 21
		set mAB_dA_total_z	= 22

		set cMAB_dA_total_u	= 23
		set cMAB_dA_total_g	= 24
		set cMAB_dA_total_r	= 25
		set cMAB_dA_total_i	= 26
		set cMAB_dA_total_z	= 27

	endif

else if ($tarsel_code_default =~ *'por'*) then

	#Observational data
	set mstar_max		= 6
	set mstar_min		= 34
	set sfr			= 7
	set sfr_min		= 8
	set sfr_max		= 9

	set mAB_dA_total_u	= 10
	set mAB_dA_total_g	= 11
	set mAB_dA_total_r	= 12
	set mAB_dA_total_i	= 13
	set mAB_dA_total_z	= 14

	set Mag_dA_total_K	= 15

	set mAB_dA_total_err_u	= 16
	set mAB_dA_total_err_g	= 17
	set mAB_dA_total_err_r	= 18
	set mAB_dA_total_err_i	= 19
	set mAB_dA_total_err_z	= 20

	set mAB_dA_total_cut_r_i	= 21
	set mAB_dA_total_cut_dmesa	= 22
	set mAB_dA_total_cut_i_lt_dmesa	= 23

	set Z_err		= 27
	set boss_target1	= 24
	set specprimary		= 25
	set warning		= 26

	set weight_fkp		= 28
	set weight_seeing	= 29
	set weight_star		= 30
	set weight_cp		= 31
	set weight_noz		= 32
	set weight_systot	= 33


endif
#----------------------------------------

#Luminosities
set L_SDSS_spheroid_u	= 99
set L_SDSS_spheroid_g	= 99
set L_SDSS_spheroid_r	= 99
set L_SDSS_spheroid_i	= 99
set L_SDSS_spheroid_z	= 99

set L_SDSS_dA_spheroid_u	= 99
set L_SDSS_dA_spheroid_g	= 99
set L_SDSS_dA_spheroid_r	= 99
set L_SDSS_dA_spheroid_i	= 99
set L_SDSS_dA_spheroid_z	= 99

set L_SDSS_disk_u	= 99
set L_SDSS_disk_g	= 99
set L_SDSS_disk_r	= 99
set L_SDSS_disk_i	= 99
set L_SDSS_disk_z	= 99

set L_SDSS_dA_disk_u	= 99
set L_SDSS_dA_disk_g	= 99
set L_SDSS_dA_disk_r	= 99
set L_SDSS_dA_disk_i	= 99
set L_SDSS_dA_disk_z	= 99

set L_SDSS_u	= 99
set L_SDSS_g	= 99
set L_SDSS_r	= 99
set L_SDSS_i	= 99
set L_SDSS_z	= 99

set L_SDSS_dA_u	= 99
set L_SDSS_dA_g	= 99
set L_SDSS_dA_r	= 99
set L_SDSS_dA_i	= 99
set L_SDSS_dA_z	= 99

set L_SDSS_dA_total_u	= 99
set L_SDSS_dA_total_g	= 99
set L_SDSS_dA_total_r	= 99
set L_SDSS_dA_total_i	= 99
set L_SDSS_dA_total_z	= 99
#----------------------------------------

#Magnitudes Standard
set Mag_dA_total_B	= 99

#----------------------------------------

#Further Magnitudes
set mag_u		= 99
set mag_g		= 99
set mag_r		= 99
set mag_i		= 99
set mag_z		= 99

set mAB_u		= 99
set mAB_g		= 99
set mAB_r		= 99
set mAB_i		= 99
set mAB_z		= 99

set mAB_u	= 99
set mAB_g	= 99
set mAB_r	= 99
set mAB_i	= 99
set mAB_z	= 99

set mAB_dA_u		= 99
set mAB_dA_g		= 99
set mAB_dA_r		= 99
set mAB_dA_i		= 99
set mAB_dA_z		= 99

set mAs_u	= 99
set mAs_g	= 99
set mAs_r	= 99
set mAs_i	= 99
set mAs_z	= 99

set mAs_dA_u		= 99
set mAs_dA_g		= 99
set mAs_dA_r		= 99
set mAs_dA_i		= 99
set mAs_dA_z		= 99

set mAs_dA_total_u		= 99
set mAs_dA_total_g		= 99
set mAs_dA_total_r		= 99
set mAs_dA_total_i		= 99
set mAs_dA_total_z		= 99
#----------------------------------------

#Target selection cuts
set mAB_total_cut_g_r	= 99

set mAB_dA_total_cut_g_r	= 99
set mAB_dA_total_cut_g_i	= 99
set mAB_dA_total_cut_u_g	= 99
set mAB_dA_total_cut_u_i	= 99
set mAB_dA_total_cut_i_z	= 99

set mAB_dA_total_cut_i_lt_dmesa_sparse	= 99

set mAB_dA_total_cut_cmesa	= 99
set mAB_dA_total_cut_cpar	= 99
set mAB_dA_total_cut_r_lt_cpar	= 99

set mAs_dA_total_cut_r_i	= 99
set mAs_dA_total_cut_dmesa	= 99
set mAs_dA_total_cut_g_r	= 99
set mAs_dA_total_cut_i_lt_dmesa	= 99
#----------------------------------------

#Emission lines
set OII_3727_ext	= 99
set OII_3727		= 99
set OII_3729_ext	= 99
set OII_3729		= 99
#----------------------------------------

#Continuum Emission lines
set OII_cont_3727_ext	= 99
set OII_cont_3727	= 99
set OII_cont_3729_ext	= 99
set OII_cont_3729	= 99

if ($SOFTLINK_CODE == Aa) then

	set L_SDSS_dA_total_u	= 12
	set L_SDSS_dA_total_g	= 13
	set L_SDSS_dA_total_r	= 14
	set L_SDSS_dA_total_i	= 99
	set L_SDSS_dA_total_z	= 99

	set MAB_dA_total_u	= 15
	set MAB_dA_total_g	= 16
	set MAB_dA_total_r	= 17
	set MAB_dA_total_i	= 99
	set MAB_dA_total_z	= 99

	set MAB_total_u	= 99
	set MAB_total_g	= 99
	set MAB_total_r	= 99
	set MAB_total_i	= 99
	set MAB_total_z	= 99

else if ($SOFTLINK_CODE == Bb) then

	set MAB_dA_total_u	= 99
	set MAB_dA_total_g	= 99
	set MAB_dA_total_r	= 99
	set MAB_dA_total_i	= 99
	set MAB_dA_total_z	= 99

	set MAB_total_u	= 17
	set MAB_total_g	= 18
	set MAB_total_r	= 19
	set MAB_total_i	= 20
	set MAB_total_z	= 21

	set OII_3727_ext	= 99
	set OII_3727		= 99
	set OII_3729_ext	= 99
	set OII_3729		= 99

endif





#DON'T TOUCH BEYOND THIS LINE!!!!!

# |		|		|		|		|		|		|		|
# |		|		|		|		|		|		|		|
# v		v		v		v		v		v		v		v
######################################################################################################################

echo '###########################################################'

set mycomp = $HOME'/'

set MAIN_PATH				=  $mycomp'anaconda/pro/'
set PLOTPATH				=  $MAIN_PATH'myRun/plot_config/'


#########################################################################
#	GENERATE A PATH HANDLER		"MY_PATH_HANDLER_FILE"		#
#########################################################################

set MY_PATH_HANDLER_FILE		=  $MAIN_PATH'myRun/run_config/mypath_handler.txt'
rm -f $MY_PATH_HANDLER_FILE
touch $MY_PATH_HANDLER_FILE

#################################################################################
#	GENERATE FILE WITH BASIS PHYSICAL SPECIFICATIONS "MY_PHYSICS_SPEC"	#
#################################################################################

set MY_PHYSICS_SPECS			= $MAIN_PATH'myRun/run_config/my_physics_specs.txt'
rm -f $MY_PHYSICS_SPECS
touch $MY_PHYSICS_SPECS

echo 'plot_all= '$plot_all		>> $MY_PHYSICS_SPECS
echo 'single= '$plot_single		>> $MY_PHYSICS_SPECS
echo 'single_mix= '$plot_single_mix	>> $MY_PHYSICS_SPECS
echo 'subplots= '$subplots		>> $MY_PHYSICS_SPECS
echo 'load_from_file= '$LOAD_FROM_FILE	>> $MY_PHYSICS_SPECS

echo 'MY_PHYSICS_SPECS' $MY_PHYSICS_SPECS	 		>> $MY_PATH_HANDLER_FILE

set method_name_array = (plotXY plotOnly filterData matchHaloCat analyseTargetSelection SMF SFRF sSFRF sfr2z sfr2mstar ssfr2mstar oh2mstar HMF mstar2mhalo mstar2mhalofunc ngal2mhalo HMF_no mstar2mhalo_no mstar2mhalofunc_no ngal2mhalo_no zgas2mstar zgas2mcold mcold2mstar mbh2mstarsph twoPCF HOD mstar2mhalovsSFR)

set method_array = ($plotXY $plotOnly $filterData $matchHaloCat $analyseTargetSelection $SMF $SFRF $sSFRF $sfr2z $sfr2mstar $ssfr2mstar $oh2mstar $HMF $mstar2mhalo $mstar2mhalofunc $ngal2mhalo $HMF_no $mstar2mhalo_no $mstar2mhalofunc_no $ngal2mhalo_no $zgas2mstar $zgas2mcold $mcold2mstar $mbh2mstarsph $twoPCF $HOD $mstar2mhalovsSFR)

#echo 559 method_name_array $method_array

#########################################################################
#	GENERATE FILTER AND CUT VALUE FILE 	"CUT_CONDS"		#
#########################################################################

set halo_col_name_array = (halo_nhalos halo_haloid halo_hostid halo_desIndex halo_nodemass halo_vmax halo_rsca halo_angmom halo_spin halo_redshift halo_x_pos halo_y_pos halo_z_pos halo_x_vel halo_y_vel halo_z_vel halo_x_vel_disp halo_y_vel_disp halo_z_vel_disp)

set halo_col_array = ($halo_nhalos $halo_haloid $halo_hostid $halo_desIndex $halo_nodemass $halo_vmax $halo_rsca $halo_angmom $halo_spin $halo_redshift $halo_x_pos $halo_y_pos $halo_z_pos $halo_x_vel $halo_y_vel $halo_z_vel $halo_x_vel_disp $halo_y_vel_disp $halo_z_vel_disp)

#set tarsel_col_name_array = (mAB_dA_total_cut_r_i mAB_dA_total_cut_dmesa mAB_dA_total_cut_g_r)

#set tarsel_col_array = ($mAB_dA_total_cut_r_i $mAB_dA_total_cut_dmesa $mAB_dA_total_cut_g_r)


set col_name_array = (ngalaxies haloid hostid satelliteNodeIndex parentIndex orphan mhalo vmax vpeak spinParameter NFW_con zgas_spheroid zgas_disk zstar_disk zstar_spheroid zhot_halo mcold_spheroid mcold_disk mcold mbh mstar_spheroid mstar_disk mstar mhot Mzgas_spheroid Mzgas_disk Mzgas Mzstar_spheroid Mzstar_disk Mzstar Mzhot_halo sfr_spheroid sfr_disk sfr x_pos y_pos z_pos x_vel y_vel z_vel L_SDSS_spheroid_u L_SDSS_spheroid_g L_SDSS_spheroid_r L_SDSS_spheroid_i L_SDSS_spheroid_z L_SDSS_dA_spheroid_u L_SDSS_dA_spheroid_g L_SDSS_dA_spheroid_r L_SDSS_dA_spheroid_i L_SDSS_dA_spheroid_z L_SDSS_disk_u L_SDSS_disk_g L_SDSS_disk_r L_SDSS_disk_i L_SDSS_disk_z L_SDSS_dA_disk_u L_SDSS_dA_disk_g L_SDSS_dA_disk_r L_SDSS_dA_disk_i L_SDSS_dA_disk_z L_SDSS_u L_SDSS_g L_SDSS_r L_SDSS_i L_SDSS_z L_SDSS_dA_u L_SDSS_dA_g L_SDSS_dA_r L_SDSS_dA_i L_SDSS_dA_z L_SDSS_dA_total_u L_SDSS_dA_total_g L_SDSS_dA_total_r L_SDSS_dA_total_i L_SDSS_dA_total_z mag_u mag_g mag_r mag_i mag_z mAB_u mAB_g mAB_r mAB_i mAB_z mAB_total_u mAB_total_g mAB_total_r mAB_total_i mAB_total_z mAB_dA_u mAB_dA_g mAB_dA_r mAB_dA_i mAB_dA_z mAB_dA_total_u mAB_dA_total_g mAB_dA_total_r mAB_dA_total_i mAB_dA_total_z MAB_dA_total_u MAB_dA_total_g MAB_dA_total_r MAB_dA_total_i MAB_dA_total_z MAB_total_u MAB_total_g MAB_total_r MAB_total_i MAB_total_z mAs_u mAs_g mAs_r mAs_i mAs_z mAs_dA_u mAs_dA_g mAs_dA_r mAs_dA_i mAs_dA_z mAs_dA_total_u mAs_dA_total_g mAs_dA_total_r mAs_dA_total_i mAs_dA_total_z mAB_total_cut_r_i mAB_total_cut_dmesa mAB_total_cut_g_r mAB_total_cut_i_lt_dmesa 'mAB_dA_total_cut_r_i' 'mAB_dA_total_cut_g_r' 'mAB_dA_total_cut_g_i' 'mAB_dA_total_cut_u_g' 'mAB_dA_total_cut_u_i' 'mAB_dA_total_cut_i_z' mAB_dA_total_cut_dmesa mAB_dA_total_cut_cmesa mAB_dA_total_cut_cpar mAB_dA_total_cut_i_lt_dmesa mAB_dA_total_cut_i_lt_dmesa_sparse mAB_dA_total_cut_r_lt_cpar mAs_dA_total_cut_r_i mAs_dA_total_cut_dmesa mAs_dA_total_cut_g_r mAs_dA_total_cut_i_lt_dmesa RA DEC Z age mean_age_stars mhalo_sat OH_gas_disk OH_gas_bulge OH_gas_disk_bulge weight OII_3727_ext OII_3727 OII_3729_ext OII_3729 nodeIndex satelliteIndex siblingIndex satelliteMergeTime isolated timeLastIsolated firstProgenitorID npros Mag_dA_total_B rhalf_bulge rhalf_disk rbulge rdisk rhalf_mass mstar_IC OII_cont_3727_ext OII_cont_3727 OII_cont_3729_ext OII_cont_3729 sfr_spheroid_inst sfr_quies_inst mstar_err mstar_9750 mstar_8400 mstar_5000 mstar_0160 mstar_0025 mstar_max mstar_min ssfr ssfr_max ssfr_min kcorr_u kcorr_g kcorr_r kcorr_i kcorr_z sfr_min sfr_max Mag_dA_total_K mAB_dA_total_err_u mAB_dA_total_err_g mAB_dA_total_err_r mAB_dA_total_err_i mAB_dA_total_err_z cmAB_dA_total_u cmAB_dA_total_g cmAB_dA_total_r cmAB_dA_total_i cmAB_dA_total_z cmAB_total_u cmAB_total_g cmAB_total_r cmAB_total_i cmAB_total_z weight_fkp weight_cp weight_noz weight_seeing weight_star weight_systot cMAB_dA_total_u cMAB_dA_total_g cMAB_dA_total_r cMAB_dA_total_i cMAB_dA_total_z cMAB_total_u cMAB_total_g cMAB_total_r cMAB_total_i cMAB_total_z boss_target1 specprimary warning Z_err weight_tot mstar_char)

set col_array = ($ngalaxies $haloid $hostid $satelliteNodeIndex $parentIndex $orphan $mhalo $vmax $vpeak $spinParameter $NFW_con $zgas_spheroid $zgas_disk $zstar_disk $zstar_spheroid $zhot_halo $mcold_spheroid $mcold_disk $mcold $mbh $mstar_spheroid $mstar_disk $mstar $mhot $Mzgas_spheroid $Mzgas_disk $Mzgas $Mzstar_spheroid $Mzstar_disk $Mzstar $Mzhot_halo $sfr_spheroid $sfr_disk $sfr $x_pos $y_pos $z_pos $x_vel $y_vel $z_vel $L_SDSS_spheroid_u $L_SDSS_spheroid_g $L_SDSS_spheroid_r $L_SDSS_spheroid_i $L_SDSS_spheroid_z $L_SDSS_dA_spheroid_u $L_SDSS_dA_spheroid_g $L_SDSS_dA_spheroid_r $L_SDSS_dA_spheroid_i $L_SDSS_dA_spheroid_z $L_SDSS_disk_u $L_SDSS_disk_g $L_SDSS_disk_r $L_SDSS_disk_i $L_SDSS_disk_z $L_SDSS_dA_disk_u $L_SDSS_dA_disk_g $L_SDSS_dA_disk_r $L_SDSS_dA_disk_i $L_SDSS_dA_disk_z $L_SDSS_u $L_SDSS_g $L_SDSS_r $L_SDSS_i $L_SDSS_z $L_SDSS_dA_u $L_SDSS_dA_g $L_SDSS_dA_r $L_SDSS_dA_i $L_SDSS_dA_z $L_SDSS_dA_total_u $L_SDSS_dA_total_g $L_SDSS_dA_total_r $L_SDSS_dA_total_i $L_SDSS_dA_total_z $mag_u $mag_g $mag_r $mag_i $mag_z $mAB_u $mAB_g $mAB_r $mAB_i $mAB_z $mAB_total_u $mAB_total_g $mAB_total_r $mAB_total_i $mAB_total_z $mAB_dA_u $mAB_dA_g $mAB_dA_r $mAB_dA_i $mAB_dA_z $mAB_dA_total_u $mAB_dA_total_g $mAB_dA_total_r $mAB_dA_total_i $mAB_dA_total_z $MAB_dA_total_u $MAB_dA_total_g $MAB_dA_total_r $MAB_dA_total_i $MAB_dA_total_z $MAB_total_u $MAB_total_g $MAB_total_r $MAB_total_i $MAB_total_z $mAs_u $mAs_g $mAs_r $mAs_i $mAs_z $mAs_dA_u $mAs_dA_g $mAs_dA_r $mAs_dA_i $mAs_dA_z $mAs_dA_total_u $mAs_dA_total_g $mAs_dA_total_r $mAs_dA_total_i $mAs_dA_total_z $mAB_total_cut_r_i $mAB_total_cut_dmesa $mAB_total_cut_g_r $mAB_total_cut_i_lt_dmesa $mAB_dA_total_cut_r_i $mAB_dA_total_cut_g_r $mAB_dA_total_cut_g_i $mAB_dA_total_cut_u_g $mAB_dA_total_cut_u_i $mAB_dA_total_cut_i_z $mAB_dA_total_cut_dmesa $mAB_dA_total_cut_cmesa $mAB_dA_total_cut_cpar $mAB_dA_total_cut_i_lt_dmesa $mAB_dA_total_cut_i_lt_dmesa_sparse $mAB_dA_total_cut_r_lt_cpar $mAs_dA_total_cut_r_i $mAs_dA_total_cut_dmesa $mAs_dA_total_cut_g_r $mAs_dA_total_cut_i_lt_dmesa $RA $DEC $Z $age $mean_age_stars $mhalo_sat $OH_gas_disk $OH_gas_bulge $OH_gas_disk_bulge $weight $OII_3727_ext $OII_3727 $OII_3729_ext $OII_3729 $nodeIndex $satelliteIndex $siblingIndex $satelliteMergeTime $isolated $timeLastIsolated $firstProgenitorID $npros $Mag_dA_total_B $rhalf_bulge $rhalf_disk $rbulge $rdisk $rhalf_mass $mstar_IC $OII_cont_3727_ext $OII_cont_3727 $OII_cont_3729_ext $OII_cont_3729 $sfr_spheroid_inst $sfr_quies_inst $mstar_err $mstar_9750 $mstar_8400 $mstar_5000 $mstar_0160 $mstar_0025 $mstar_max $mstar_min $ssfr $ssfr_max $ssfr_min $kcorr_u $kcorr_g $kcorr_r $kcorr_i $kcorr_z $sfr_min $sfr_max $Mag_dA_total_K $mAB_dA_total_err_u $mAB_dA_total_err_g $mAB_dA_total_err_r $mAB_dA_total_err_i $mAB_dA_total_err_z $cmAB_dA_total_u $cmAB_dA_total_g $cmAB_dA_total_r $cmAB_dA_total_i $cmAB_dA_total_z $cmAB_total_u $cmAB_total_g $cmAB_total_r $cmAB_total_i $cmAB_total_z $weight_fkp $weight_cp $weight_noz $weight_seeing $weight_star $weight_systot $cMAB_dA_total_u $cMAB_dA_total_g $cMAB_dA_total_r $cMAB_dA_total_i $cMAB_dA_total_z $cMAB_total_u $cMAB_total_g $cMAB_total_r $cMAB_total_i $cMAB_total_z $boss_target1 $specprimary $warning $Z_err $weight_tot $mstar_char)

set total_col_name_array=()
set total_col_array=()

foreach item ($halo_col_name_array)
	set total_col_name_array = ($total_col_name_array $item)
end
foreach item ($halo_col_array)
	set total_col_array = ($total_col_array $item)
end

foreach item ($col_name_array)
	set total_col_name_array = ($total_col_name_array $item)
end
foreach item ($col_array)
	set total_col_array = ($total_col_array $item)
end

set analyse_tarsel_name_array = (analyse_tarsel_dmesa_i analyse_tarsel_g-r_mstar analyse_tarsel_r-i_mstar analyse_tarsel_g-i_mstar analyse_tarsel_u-r_mstar analyse_tarsel_g-z_mstar analyse_tarsel_r-i_i analyse_tarsel_i_mstar analyse_tarsel_r_mstar analyse_tarsel_g_mstar analyse_tarsel_Ii_i analyse_tarsel_Ii_mstar analyse_tarsel_Ir_r analyse_tarsel_Ir_mstar analyse_tarsel_Ig_g analyse_tarsel_Ig_mstar analyse_tarsel_histo analyse_tarsel_sfr_mstar analyse_tarsel_g-r_r analyse_tarsel_u-r_r analyse_tarsel_mbh_mstar_spheroid analyse_tarsel_mcold_mstar 'analyse_tarsel_g-r_u-g' analyse_tarsel_Mzgas_mstar analyse_tarsel_Mzgas_mcold 'analyse_tarsel_r-i_g-r')

set analyse_tarsel_array = ($analyse_tarsel_dmesa_i $analyse_tarsel_gminr_mstar $analyse_tarsel_rmini_mstar $analyse_tarsel_gmini_mstar $analyse_tarsel_uminr_mstar $analyse_tarsel_gminz_mstar $analyse_tarsel_rmini_i $analyse_tarsel_i_mstar $analyse_tarsel_r_mstar $analyse_tarsel_g_mstar $analyse_tarsel_Ii_i $analyse_tarsel_Ii_mstar $analyse_tarsel_Ir_r $analyse_tarsel_Ir_mstar $analyse_tarsel_Ig_g $analyse_tarsel_Ig_mstar $analyse_tarsel_histo $analyse_tarsel_sfr_mstar $analyse_tarsel_gminr_r $analyse_tarsel_uminr_r $analyse_tarsel_mbh_mstar_spheroid $analyse_tarsel_mcold_mstar $analyse_tarsel_gminr_uming $analyse_tarsel_Mzgas_mstar $analyse_tarsel_Mzgas_mcold $analyse_tarsel_rmini_gminr)

#echo $analyse_tarsel_name_array 
#echo $analyse_tarsel_array

#echo "$total_col_array"
#echo "$total_col_name_array"

if ($output_filename_code_default =~ CMASs*) then
	echo 'set cut values --> CMASS'
	set set_cut_values=(`$MAIN_PATH'myRun/set_cut_values_config_CMASS.sh' $MAIN_PATH "$total_col_name_array" "$total_col_array" $SOFTLINK_CODE $UNIT_CODE`)
else if ($output_filename_code_default =~ LOW*) then
	echo 'set cut values --> LOWZ'
	set set_cut_values=(`$MAIN_PATH'myRun/set_cut_values_config_LOWZ.sh' $MAIN_PATH "$total_col_name_array" "$total_col_array" $SOFTLINK_CODE $UNIT_CODE`)
else if ($output_filename_code_default == 'LRG') then
	echo 'set cut values --> LRG'
	set set_cut_values=(`$MAIN_PATH'myRun/set_cut_values_config_LRG.sh' $MAIN_PATH "$total_col_name_array" "$total_col_array" $SOFTLINK_CODE $UNIT_CODE`)
else if ($output_filename_code_default =~ *EW) then
	echo 'set cut values --> filter'
	set set_cut_values=(`$MAIN_PATH'myRun/set_cut_values_config_filter.sh' $MAIN_PATH "$total_col_name_array" "$total_col_array" $SOFTLINK_CODE $UNIT_CODE`)
else
	echo 'set cut values --> default'
	set set_cut_values=(`$MAIN_PATH'myRun/set_cut_values_config.sh' $MAIN_PATH "$total_col_name_array" "$total_col_array" $SOFTLINK_CODE $UNIT_CODE`)
endif

echo 'MYCUT_VALUES_CONFIGFILE' "$set_cut_values[1]"	>> $MY_PATH_HANDLER_FILE
echo 'MY_UNIT_CORRECT_FILE' "$set_cut_values[2]"	>> $MY_PATH_HANDLER_FILE
echo 'MY_NAME_CONV_MAP_FILE' "$set_cut_values[3]"	>> $MY_PATH_HANDLER_FILE

#########################################################################
#	GENERATE COL NAME ID AND FILE 	"MYCOLID_DATAFILE"		#
#########################################################################

set MYCOLID_DATAFILE    		=  $MAIN_PATH'myRun/run_config/col_name_id_file.txt'

#rm -f $MYCOLID_DATAFILE
#touch $MYCOLID_DATAFILE  

set i=1
foreach item ($col_name_array) 	

	echo $item "$col_array[$i]"  >> $MYCOLID_DATAFILE 
	@ i++
end

echo 'MYCOLID_DATAFILE' $MYCOLID_DATAFILE 	>> $MY_PATH_HANDLER_FILE

#########################################################################
#	GENERATE RUN SPECIFIC FILES 	"RUNFILE"			#
#########################################################################

set PATH_TO_SOFTLINKS   		=  $MAIN_PATH'myRun/mylinks/'$SOFTLINK_CODE
set SOFTLINK_DIRECTORY_CODE 		= `echo $SOFTLINK_CODE | head -1 | awk -F'_' '{print $1}'`


#########################################################################
#	PLOT CONFIGURATION FILES 	"plotname_PLOTCONFIG"		#
#########################################################################

set plot_config_array=(`$MAIN_PATH'myRun/root_config.sh' $MAIN_PATH $PLOTPATH "$method_name_array" "$method_array" '' ''`)

if ($plotOnly == 1) then
	set plot_plotOnly_config_array=(`$MAIN_PATH'myRun/root_config.sh' $MAIN_PATH $PLOTPATH "$total_col_name_array" "$total_col_array" 'plotOnly' '99'`)
endif

if ($analyseTargetSelection == 1) then
	set plot_analyseTargetSelection_config_array=(`$MAIN_PATH'myRun/root_config.sh' $MAIN_PATH $PLOTPATH "$analyse_tarsel_name_array" "$analyse_tarsel_array" 'analyseTargetSelection' '0'`)
endif

#########################################################################
#	GENERATING PLOT SPECIFIC CONFIGURATIONS		"MAP_PLOT_FILE" #
#########################################################################

set RUN_NR			= 'RUN2'
set MAP_PLOT_FILE    		=  $MAIN_PATH'myRun/run_config/'$RUN_NR'_name_binning_offset_map.txt'

rm -f $MAP_PLOT_FILE
touch $MAP_PLOT_FILE

echo 'MAP_PLOT_FILE' $MAP_PLOT_FILE	>> $MY_PATH_HANDLER_FILE


#########################################################################
#	GENERATE HISTO KEYWORDS FILE 		"HISTO_CONFIGFILE"	#
#########################################################################

set HISTO_CONFIGFILE = $MAIN_PATH'myRun/run_config/histo_configfile_'$SOFTLINK_CODE'.txt'

rm -f $HISTO_CONFIGFILE
touch $HISTO_CONFIGFILE

echo 'HISTO_CONFIGFILE' $HISTO_CONFIGFILE	>> $MY_PATH_HANDLER_FILE

$MAIN_PATH'myRun/physical_statistical_par_config.sh' $MAP_PLOT_FILE 'False' $HISTO_CONFIGFILE "$method_name_array" "$method_array" "$total_col_name_array" "$total_col_array" "$analyse_tarsel_name_array" "$analyse_tarsel_array"


#########################################################################
#	GENERATE MARKER AND COLOR MAP FILE	"MARKER_COLOR_MAP_FILE"	#
#########################################################################

set MARKER_COLOR_MAP_FILE = $MAIN_PATH'myRun/run_config/marker_color_map.txt'

rm -f $MARKER_COLOR_MAP_FILE
touch $MARKER_COLOR_MAP_FILE
echo 'MARKER_COLOR_MAP_FILE' $MARKER_COLOR_MAP_FILE	>> $MY_PATH_HANDLER_FILE

#########################################################################
#	GENERATE RUNFILE			"RUNFILE"		#
#########################################################################

set RUNFILE =  $MAIN_PATH'myRun/input_cats/input_cat.txt'

rm -f $RUNFILE
touch $RUNFILE

echo 'RUNFILE' $RUNFILE					>> $MY_PATH_HANDLER_FILE

#########################################################################
#	GENERATE CUTE IMPUT FILES 	"CUTE_INPUT_FILES"		#
#########################################################################

#set MYCUTE_PARAMFILE   		=  $mycomp'anaconda/CUTE/'$param_file_CUTE

#rm -f MYCUTE_PARAMFILE
#touch MYCUTE_PARAMFILE  

#$MAIN_PATH'myRun/set_CUTE_config.sh' $MYCUTE_PARAMFILE 
#echo 'MYCUTE_PARAMFILE' $MYCUTE_PARAMFILE	>> $MY_PATH_HANDLER_FILE


#########################################################################
#	LINESTYLE AND MARKER COMPLIATION	"MATPLOT_STYLE"		#
#########################################################################

set MATPLOT_LINESTYLE =  $MAIN_PATH'myRun/input_cats/matplot_linestyle.txt'

rm -f $MATPLOT_LINESTYLE
touch $MATPLOT_LINESTYLE

echo 'MATPLOT_LINESTYLE' $MATPLOT_LINESTYLE			>> $MY_PATH_HANDLER_FILE

set MATPLOT_COL =  $MAIN_PATH'myRun/input_cats/matplot_col.txt'

rm -f $MATPLOT_COL
touch $MATPLOT_COL

echo 'MATPLOT_COL' $MATPLOT_COL					>> $MY_PATH_HANDLER_FILE

set MATPLOT_MARKERSTYLE =  $MAIN_PATH'myRun/input_cats/matplot_markerstyle.txt'

rm -f $MATPLOT_MARKERSTYLE
touch $MATPLOT_MARKERSTYLE

echo 'MATPLOT_MARKERSTYLE' $MATPLOT_MARKERSTYLE			>> $MY_PATH_HANDLER_FILE


set MATPLOT_MARKERCOL =  $MAIN_PATH'myRun/input_cats/matplot_markercol.txt'

rm -f $MATPLOT_MARKERCOL
touch $MATPLOT_MARKERCOL

echo 'MATPLOT_MARKERCOL' $MATPLOT_MARKERCOL			>> $MY_PATH_HANDLER_FILE

set MATPLOT_COLORMAP =  $MAIN_PATH'myRun/input_cats/matplot_colormap.txt'

rm -f $MATPLOT_COLORMAP
touch $MATPLOT_COLORMAP

echo 'MATPLOT_COLORMAP' $MATPLOT_COLORMAP			>> $MY_PATH_HANDLER_FILE

set MATPLOT_COLORBAR =  $MAIN_PATH'myRun/input_cats/matplot_colorbar.txt'

rm -f $MATPLOT_COLORBAR
touch $MATPLOT_COLORBAR

echo 'MATPLOT_COLORBAR' $MATPLOT_COLORBAR			>> $MY_PATH_HANDLER_FILE

set MATPLOT_PLOTTYPE =  $MAIN_PATH'myRun/input_cats/matplot_plottype.txt'

rm -f $MATPLOT_PLOTTYPE
touch $MATPLOT_PLOTTYPE

echo 'MATPLOT_PLOTTYPE' $MATPLOT_PLOTTYPE			>> $MY_PATH_HANDLER_FILE
#########################################################################
#	GENERATE MODEL CONFIGURATION FILES	"CONFIGFILE"		#
#########################################################################



set CALI_COUNT_uc 	= 0
set CALI_COUNT 		= 0
set BOX_COUNT_62 	= 0
set BOX_COUNT_125 	= 0
set BOX_COUNT_400 	= 0

set i=0
set j=0


foreach MODEL (`ls -d $PATH_TO_SOFTLINKS/*`)



	
	#########################################################################
	# FIND / SET MODEL NAME AND PATH TO CATALOG				#
	#########################################################################
	echo 'MODEL:' $MODEL
	set NR_MODELS = `ls -l $PATH_TO_SOFTLINKS/* | egrep -c '^'`
	set NAME = `basename $MODEL`
	set PATH_TO_CATALOG = $MAIN_PATH'data/'

	set CAT_NAME_IN_PLOT$j = `echo $NAME | head -1 | awk -F'_' '{print $1}'`
	set dummy = `expr $j + 1`
	set CAT_NAME_IN_PLOT$dummy = 'dummy'

	
	if ( $NAME == SAG_1Gpc) then
		set MANUAL_REDSHIFT_INPUT = $SAG_REDSHIFT_INPUT
	else if ( $NAME =~ SAG*125 ) then
		set MANUAL_REDSHIFT_INPUT = $SAG_125_REDSHIFT_INPUT
	else if ( $NAME =~ SAGE_* ) then
		set MANUAL_REDSHIFT_INPUT = $SAGE_REDSHIFT_INPUT
		
	else if ( $NAME =~ Galacticus_1Gp* ) then
		set MANUAL_REDSHIFT_INPUT = $GALACTICUS_1Gpc_REDSHIFT_INPUT

	else if ( $NAME =~ Galacticus_40* ) then
		set MANUAL_REDSHIFT_INPUT = $GALACTICUS_400Mpc_REDSHIFT_INPUT	
	else
		set MANUAL_REDSHIFT_INPUT = $MANUAL_REDSHIFT_INPUT_default
	endif
		
	set REDSHIFT_PLOT$j = $MANUAL_REDSHIFT_INPUT
	set dummy = `expr $j + 1`
	set REDSHIFT_PLOT$dummy = 'dummy'

	if ($NAME =~ *uc*) 	then
		set CALI_CODE 		= 'uncalibrated'
		set CALI_COUNT_uc 	= 1
	else
		set CALI_CODE 		= 'calibrated'
		set CALI_COUNT 		= 1
	endif

	if ($NAME =~ *125) 	then
		set MYSNAPS 		= `echo $MYSNAPS_125`
		set BOX_SIZE 		= '125Mpc/h'
		set YLIM 		= '-5'
		set BOX_COUNT_62 	= 1
	
	else if ($NAME =~ *1Gp*) 	then
		set MYSNAPS		= `echo $MYSNAPS_1Gpc`
		set BOX_SIZE 		= '1Gpc/h'

	else if ($NAME =~ *62*) 	then
		set MYSNAPS 		= `echo $MYSNAPS_62`
		set BOX_SIZE 		= '62.5Mpc/h'
		set YLIM 		= '-5'
		set BOX_COUNT_125 	= 1

	else if ($NAME =~ *400*) 	then
		set MYSNAPS 		= `echo $MYSNAPS_400Mpc`
		set BOX_SIZE 		= '400Mpc/h'
		set YLIM 		= '-5'
		set BOX_COUNT_400Mpc 	= 1

	else
		set MYSNAPS 		= `echo $MYSNAPS_default`
		set BOX_SIZE 		= '0.0'
		set YLIM 		= '-5'
		set BOX_COUNT		= 1
	endif	
	echo $NAME 'manual z:' $MANUAL_REDSHIFT_INPUT 'MYSNAPS: ' $MYSNAPS, 'BOX_SIZE:' $BOX_SIZE

	#########################################################################
	#	GENERATE PHYSICAL UNITS FILE		"PHYSICAL_UNITS"	#
	#########################################################################

	set PHYSICAL_UNITS =  $MAIN_PATH'myRun/input_cats/'$NAME'_physical_units.txt'

	rm -f $PHYSICAL_UNITS
	touch $PHYSICAL_UNITS

	echo 'PHYSICAL_UNITS' $PHYSICAL_UNITS					>> $MY_PATH_HANDLER_FILE

	#########################################################################
	#	GENERATE SAMs SCALE FACTOR MAPPING  "SAM_SCALE_FACTOR_MAPPING"	#
	#########################################################################

	set SAM_SCALE_FACTOR_MAPPING = $MAIN_PATH'myRun/run_config/SAM_SCALE_FACTOR_MAPPING.txt'

	rm -f $SAM_SCALE_FACTOR_MAPPING
	touch $SAM_SCALE_FACTOR_MAPPING

	echo 'SAM_SCALE_FACTOR_MAPPING' $SAM_SCALE_FACTOR_MAPPING		>> $MY_PATH_HANDLER_FILE
	
	#########################################################################
	# SET MODEL SPECIFIC OUTPUT FORMATS / COLOURS / MARKERS	/ BOXSIZE	#
	#########################################################################

	set result=(`$MAIN_PATH'myRun/cat_config.sh' $NAME $MARKER_COLOR_MAP_FILE "$MYSNAPS" $SAM_SCALE_FACTOR_MAPPING $UNIT_CODE`)

	set NAME_SL="$result[1]"
	set SUBDIRNAME="$result[2]"
	set SIM_LITTLE_H="$result[3]"
	set GENERAL_LITTLE_H_CORR="$result[4]"
	set H_CORR="$result[5]"
	set MPC_CORR="$result[6]"
	set MASS_CORR="$result[7]"
	set GYR_CORR="$result[8]"
	set LUM_CORR="$result[9]"
	set NO_CORR="$result[10]"
	set CONV_TO_AB_MAG="$result[11]"
	set KMS1_CORR="$result[12]"
	set COMV_CORR="$result[13]"
	set APPLY_K_CORR_APP="$result[14]"
	set APPLY_Z_BOOST="$result[15]"
	set LITTLE_H_1_CORR="$result[16]"
	set e3_CORR="$result[17]"
	set MAB_CORR="$result[18]"
	set LITTLE_H_2_CORR="$result[19]"
	set APPLY_K_CORR_ABS="$result[20]"
	set IMF_CORR="$result[21]"
	set GYR1_CORR="$result[22]"

	#########################################################################
	# FIND / SET DATA PATH / SOFTLINK-REDIRECTION				#
	#########################################################################
	echo ' '	
	echo $NAME
	echo '+++++++++++++++++++++++'
	echo ' '
	echo 'NAME_SL:' $NAME_SL, 'SUBDIRNAME:' $SUBDIRNAME

	set test_sl = `ls -la $PATH_TO_CATALOG$NAME  | grep "\->" | tail -1 | awk -F'->' '{print $2}'`


	if ( $test_sl != '' ) then
		set SOFTLINK_TO_DATA =  `ls -la $PATH_TO_CATALOG$NAME  | grep "\->" | tail -1 | awk -F'->' '{print $2}'` 	
		set MODEL =  `echo $SOFTLINK_TO_DATA | awk -F$NAME_SL '{print $1}'`
		echo 'SL exists! SOFTLINK_TO_DATA DER DA WO WIRKLICH IST:' $SOFTLINK_TO_DATA
	else
		set SOFTLINK_TO_DATA = $PATH_TO_CATALOG$NAME'/'
		echo 'SL NOT exists! Path to CATLOG:' $SOFTLINK_TO_DATA
	endif


	#########################################################################
	# SET DETAILS FOR REDSHIFT AND FILEFORMAT MAPPING, BOX_SIZE		#
	#########################################################################

	set HDF5_CODE = 'False'
	set BINARY_CODE = 'False'
	set FITS_CODE = 'False'
	set CROSSMATCH_CODE = 'False'

	if ($NAME =~ SAGE*) then
		if ($HDF5_read_code =~ *DF5) then
			set HDF5_CODE = 'True'
		else 
			set BINARY_CODE = 'True'
		endif

	else if ($HDF5_read_code =~ HDF* || $HDF5_read_code =~ SAMHDF*) then
		set HDF5_CODE = 'True'
	else if ($HDF5_read_code == CROSSMATCH) then
		set CROSSMATCH_CODE = 'True'	
	else
		if (`ls $SOFTLINK_TO_DATA/* | grep "hdf5"` =~ *hdf5) 		set HDF5_CODE = 'True'
		if (`ls $SOFTLINK_TO_DATA/* | grep "h5"` =~ *h5) 		set HDF5_CODE = 'H5'
		if (`ls $SOFTLINK_TO_DATA/* | grep "dat"` =~ *niftydat) 	set BINARY_CODE = 'True'
		if (`ls $SOFTLINK_TO_DATA/* | grep "fits"` =~ *fits) 		set FITS_CODE = 'True'

	endif
	
	#echo 'HDF5_CODE:' $HDF5_CODE, 'BINARY_CODE:' $BINARY_CODE

	if ($HDF5_CODE == 'True') then
		echo 'hdf5 --> from softlink'
		set CUSTOM_SNAPNAME			= 'hdf5'
		set SNAPIDZ_MAPPING 			= False
		if ($NAME =~ *sample*) set MANUAL_REDSHIFT_INPUT	= 0.0
	else if ($HDF5_CODE == 'H5') then
		echo '--> H5, SNAPIDZ_MAPPING --> TRUE'
		set CUSTOM_SNAPNAME			= 'hdf5'
		set SNAPIDZ_MAPPING 			= True
		set HDF5_CODE = 'True'
	else if ($BINARY_CODE == 'True') then
		echo '--> binary'
		set CUSTOM_SNAPNAME			= 'binary'
		set SNAPIDZ_MAPPING 			= False

	else if ($CROSSMATCH_CODE == 'True') then
		echo '--> CROSSMATCH two catalogs'
		set CUSTOM_SNAPNAME			= 'crossmatch'
		set SNAPIDZ_MAPPING 			= False

	else if ($FITS_CODE == 'True') then
		echo '--> FITS'
		set CUSTOM_SNAPNAME			= 'fits'
		set SNAPIDZ_MAPPING 			= False

	else
		echo '--> ASCII'
		set CUSTOM_SNAPNAME			= 'txt'
		set SNAPIDZ_MAPPING 			= True
		if ($NAME =~ *sample* || $NAME =~ *DR*) then
			set MANUAL_REDSHIFT_INPUT	= $MANUAL_REDSHIFT_INPUT_default
			set SNAPIDZ_MAPPING 		= False
		endif
	endif






	#########################################################################
	# CREATE CONFIG-FILE							#
	#########################################################################

	set CONFIGFILE = `echo $PATH_TO_CATALOG$NAME'/'$NAME"_config.txt"`
	echo $CONFIGFILE

	rm -f $CONFIGFILE
	touch $CONFIGFILE

	echo name_of_cataloge':'$NAME 				>> $RUNFILE

	echo $SOFTLINK_TO_DATA				>> $RUNFILE
	echo $PATH_TO_CATALOG$NAME'/' 				>> $RUNFILE

	#echo 'SOFTLINK_TO_DATA:' $SOFTLINK_TO_DATA
	#echo 'PATH_TO_CATALOG$NAME:' $PATH_TO_CATALOG$NAME'/'  


	#########################################################################
	# SET CONFIG-FILE AND SNAPSHOT INPUTS					#
	#########################################################################

	if ($NAME =~ SAG* || $NAME =~ suss* || $CUSTOM_SNAPNAME == 'crossmatch') then
		set FILENAME_LIST = `echo $PATH_TO_CATALOG$NAME'/'$NAME"_file_names.txt"`
		echo $FILENAME_LIST

		rm -f $FILENAME_LIST
		touch $FILENAME_LIST
	endif

	echo $MYSNAPS
	set a=0
	foreach ISNAP ($MYSNAPS)
		echo 'MYSNAPS:' $MYSNAPS 'NAME:' $NAME 'ISNAP:' $ISNAP 'CUSTOM_SNAPNAME:' $CUSTOM_SNAPNAME

		if ($NAME =~ *1Gpc* || $NAME =~ *400* && $NAME !~ SAG_* || $NAME =~ suss*) then
			set SNAPID 		= `echo $ISNAP`
		else if ($NAME =~ *DR* || $NAME =~ *sample* || $NAME =~ CMAS* ) then
			set SNAPID 		= $ISNAP
		else
			set SNAPID 		= `echo $ISNAP | awk '{printf ("%03d\n",$1)}'`
		endif

		set c=1
		set k=1
		foreach item ($method_array)
			if ($item == 1)  then
				if ("$method_name_array[$c]" == plotOnly && $plotOnly == 1) then				
					
					foreach plotOnly ($plot_plotOnly_config_array)
						#echo 'HERE:' "$plot_plotOnly_config_array[$k]"
						echo redshift$a'= '$SNAPID >> "$plot_plotOnly_config_array[$k]"
						@ k++
					end
				
				else if ("$method_name_array[$c]" == analyseTargetSelection && $analyseTargetSelection == 1) then				
					
					foreach plotOnly ($plot_config_array)
						#echo 'HERE:' "$plot_config_array[$k]"
						echo redshift$a'= '$SNAPID >> "$plot_config_array[$k]"
						@ k++
					end
				else
					echo redshift$a'= '$SNAPID >> "$plot_config_array[$c]"
				endif
			endif
			@ c++
		end


		if ( $CUSTOM_SNAPNAME == 'hdf5' || $CUSTOM_SNAPNAME == 'binary' || $CUSTOM_SNAPNAME == 'crossmatch') then  
			set key = $NAME
		else
			set key = `ls $MODEL/*$SNAPID*`
			
		endif
		
		echo 'key:' $key, $SNAPID
		foreach CATALOG ($key)

			if ( $CUSTOM_SNAPNAME == 'hdf5'  && $PATH_TO_CATALOG$NAME'/' != $SOFTLINK_TO_DATA ) then
				set SNAPNAME = `echo $NAME`
				echo 'option 1'
			else if ( $CUSTOM_SNAPNAME == 'hdf5' && $NAME =~ suss* ) then
				set SNAPNAME = `echo $key`
				echo 'option 2'
			else if ( $CUSTOM_SNAPNAME == 'hdf5' && $NAME !~ SAG_* && $PATH_TO_CATALOG$NAME'/' == $SOFTLINK_TO_DATA ) then
				set SNAPNAME = `echo $key`
				echo 'option 3'
			else if ( $CUSTOM_SNAPNAME == 'hdf5' && $NAME =~ SAG_* && $PATH_TO_CATALOG$NAME'/' == $SOFTLINK_TO_DATA ) then
				set SNAPNAME = `basename $key`
				echo 'option 4'
			else if ( $CUSTOM_SNAPNAME == 'binary' && $NAME == SAGE_1Gpc && $CUSTOM_SNAPNAME == 'fits') then
				set SNAPNAME = `echo $key`
				echo 'option 5'
			else if ( $CUSTOM_SNAPNAME == 'crossmatch') then
				set SNAPNAME = `echo $PATH_TO_CATALOG$NAME`
				echo 'option 6 crossmatch'
			else
				set SNAPNAME = `basename $CATALOG`
				echo 'else'
			endif

			echo 'Snapname:' $SNAPNAME 'softlink:' $SOFTLINK_TO_DATA $NAME


			#Create a file list for SAG or SAGE!
			if ($NAME =~ SAG*) then

				if ($NAME =~ SAG_*) then
					set PATH_TO_SNAPSHOT = `basename $SOFTLINK_TO_DATA*$SNAPID | awk '{print $1}'`'/'

					foreach File (`ls -d $SOFTLINK_TO_DATA$PATH_TO_SNAPSHOT*`)
						#echo 'File:' $File
						echo `basename $File` >> $FILENAME_LIST
					end
				else
					if ($HDF5_read_code == 'SAMHDF5') then
						
						echo 'SAGE:' $HDF5_read_code $PATH_TO_CATALOG$NAME'/' $SNAPID				
						foreach File (`ls -d $PATH_TO_CATALOG$NAME'/'*$SNAPID*`)
							echo 'File:' $File
							echo $File >> $FILENAME_LIST
						end


					else
						set PATH_TO_SNAPSHOT = ''
						#echo 'SAGE:' $SOFTLINK_TO_DATA$PATH_TO_SNAPSHOT $ISNAP
					
						foreach File (`ls -d $SOFTLINK_TO_DATA$PATH_TO_SNAPSHOT*$SNAPID*`)
							#echo 'File:' $File
							echo `basename $File` >> $FILENAME_LIST
						end
					endif
				endif
			else if ($NAME =~ suss*) then
				#echo '$SOFTLINK_TO_DATA ~$SNAPID' $SOFTLINK_TO_DATA $NAME $SNAPID
				foreach File (`ls -d $SOFTLINK_TO_DATA*suss*`)
					echo 'File:' `basename $File`
					echo `basename $File` >> $FILENAME_LIST
				end						

			else if ($CUSTOM_SNAPNAME == crossmatch) then
				#echo '$SOFTLINK_TO_DATA ~$SNAPID' $SOFTLINK_TO_DATA $NAME $SNAPID
				foreach File (`ls -d $SNAPNAME'/crossmatch/'*`)
					echo 'File:' `basename $File`
					echo `basename $File` >> $FILENAME_LIST
				end
			endif

			#create a snapshot list from which the python-routines read the different redshifts
			if ($a == 0) then
				set SNAP_LIST 	= $SNAPID
			else
				set SNAP_LIST = `echo $SNAP_LIST,$SNAPID`

				#Keep that lines or the CATALOG_config.txt will be printed double!
				rm -f $CONFIGFILE
				touch $CONFIGFILE
			endif

			echo 'SNAP_LIST:' $SNAP_LIST 'SNAPNAME:' $SNAPNAME

			if ($SNAPNAME =~ *config* || $SNAPNAME =~ *idzred* || $SNAPNAME =~ *READ*) then
				echo $SNAPNAME': NOT A SNAPSHOT FILE'
			else
				echo  $SNAPNAME 			>> $RUNFILE
				set b = `expr $a \* $NR_MODELS`
		
				set z_i = `expr $j + $b`

				set c=1
				set k=1
				foreach item ($method_array)

					if ($item == 1) then

		 				if ("$method_name_array[$c]" == plotOnly && $plotOnly == 1) then				
							foreach plotOnly ($plot_plotOnly_config_array)
								echo 'plot_legend'$z_i'= '$NAME >> "$plot_plotOnly_config_array[$k]"
								@ k++
							end
		 			else if ("$method_name_array[$c]" == analyseTargetSelection && $analyseTargetSelection == 1) then				
							foreach analyseTargetSelection ($plot_analyseTargetSelection_config_array)
								echo 'plot_legend'$z_i'= '$NAME >> "$plot_analyseTargetSelection_config_array[$k]"
								@ k++
							end
					else
							echo 'plot_legend'$z_i'= '$NAME >> "$plot_config_array[$c]"
						endif
						
					endif
					@ c++
				end

				@ i++

				echo 'halocat_code= False'			>> $CONFIGFILE
				if ($NAME =~ SAGE_* && $CUSTOM_SNAPNAME !~ *df5) then
					echo 'data_format= BINARY_SAGE'              	>> $CONFIGFILE
					echo 'format_info= choose_all'     		>> $CONFIGFILE
					echo 'delimiter= \t'		     	 	>> $CONFIGFILE
					echo 'fileformat= '				>> $CONFIGFILE
				else if ($CUSTOM_SNAPNAME =~ *match) then
					echo 'data_format= '$HDF5_read_code  		>> $CONFIGFILE
					echo 'format_info= choose_all'     		>> $CONFIGFILE
					echo 'delimiter= \t'		     	 	>> $CONFIGFILE
					echo 'fileformat= 	'			>> $CONFIGFILE
				else if ($CUSTOM_SNAPNAME =~ *its) then
					echo 'data_format= '$HDF5_read_code  		>> $CONFIGFILE
					echo 'format_info= choose_all'     		>> $CONFIGFILE
					echo 'delimiter= \t'		     	 	>> $CONFIGFILE
					echo 'fileformat= 	'			>> $CONFIGFILE
				else if ($SNAPNAME =~ *dat || $CUSTOM_SNAPNAME =~ *dat) then
					echo 'data_format= BINARY'              	>> $CONFIGFILE
					echo 'format_info= choose_all'     		>> $CONFIGFILE
					echo 'delimiter= \t'		     	 	>> $CONFIGFILE
					echo 'fileformat= niftydat'			>> $CONFIGFILE
				else if ($SNAPNAME =~ *xt || $CUSTOM_SNAPNAME =~ *xt) then
					echo 'data_format= CATASCII' 		 	>> $CONFIGFILE
					echo 'format_info= shaped'	     	>> $CONFIGFILE
					echo 'delimiter= '		 	>> $CONFIGFILE
					echo 'fileformat= txt'			>> $CONFIGFILE
				else if ($SNAPNAME =~ *MDPL*) then
					echo 'data_format= ASCII'               	>> $CONFIGFILE
					echo 'format_info= unshaped'	     	>> $CONFIGFILE
					echo 'delimiter= 	'		 	>> $CONFIGFILE
					echo 'fileformat= txt'				>> $CONFIGFILE
				else if ($SNAPNAME =~ *tree*) then
					echo 'data_format= SAMHDF5'               	>> $CONFIGFILE
					echo 'format_info= unshaped'	     	>> $CONFIGFILE
					echo 'delimiter= 	'		 	>> $CONFIGFILE
					echo 'fileformat= h5'				>> $CONFIGFILE
					echo 'halocat_code= True'			>> $CONFIGFILE
				else if ($SNAPNAME =~ *hdf5 || $CUSTOM_SNAPNAME =~ *hdf5) then
					if ($HDF5_read_code == 'SAMHDF5') then

						echo 'data_format= SAMHDF5'            	>> $CONFIGFILE
					else 
						echo 'data_format= HDF5'               	>> $CONFIGFILE
					endif

					echo 'format_info= unshaped'	     	        >> $CONFIGFILE
					echo 'delimiter= \t'		     	 	>> $CONFIGFILE
					echo 'fileformat= hdf5'				>> $CONFIGFILE
				endif

				#needen for unshaped ASCII-file read in (old Nifty catalogs old!)
				echo 'mydata_block_offset= 2'               	 	>> $CONFIGFILE
				echo 'nr_col= 55'              	 		>> $CONFIGFILE

				#Take care of proper NR OF ROWS TO BE READ:
				if ($NAME =~ *125*) then
					echo 'nr_rows= 550000'		>> $CONFIGFILE
				else if ($NAME =~ *400*) then
					echo 'nr_rows= 200000000'   	>> $CONFIGFILE
				else if ($NAME =~ *62* || $HOME_MODE == True || $HOME =~ /home/*) then
					echo 'nr_rows= 5000000'   >> $CONFIGFILE
				else
					echo 'nr_rows= 200000000'   >> $CONFIGFILE
					#echo 'nr_rows= 4600000'   >> $CONFIGFILE
				endif					
				
				#Unit correction specifics
				echo 'correct_little_h= '$GENERAL_LITTLE_H_CORR	>> $CONFIGFILE
				echo 'unit_corr_h= '$H_CORR				>> $CONFIGFILE
				echo 'unit_corr_Mpc= '$MPC_CORR	   			>> $CONFIGFILE
				echo 'unit_corr_mass= '$MASS_CORR   			>> $CONFIGFILE
				echo 'unit_corr_Gyr= '$GYR_CORR	   			>> $CONFIGFILE
				echo 'unit_corr_Gyr-1= '$GYR1_CORR 			>> $CONFIGFILE
				echo 'unit_corr_lum= '$LUM_CORR	   			>> $CONFIGFILE
				#non-orphan correction currently not set!
				echo 'unit_corr_no_corr= '$NO_CORR			>> $CONFIGFILE
				echo 'unit_corr_id_num= -'				>> $CONFIGFILE
				echo 'unit_corr_int_num= -'				>> $CONFIGFILE
				echo 'unit_corr_kms-1= '$KMS1_CORR	   		>> $CONFIGFILE
				echo 'unit_corr_MAB_CORR= '$MAB_CORR			>> $CONFIGFILE
				echo 'unit_corr_COMV= '$COMV_CORR	   		>> $CONFIGFILE
				echo 'conv_to_AB_mag= '$CONV_TO_AB_MAG   		>> $CONFIGFILE
				echo 'apply_k_corr_app= '$APPLY_K_CORR_APP		>> $CONFIGFILE
				echo 'apply_k_corr_abs= '$APPLY_K_CORR_ABS		>> $CONFIGFILE
				echo 'apply_z_boost= '$APPLY_Z_BOOST			>> $CONFIGFILE
				echo 'unit_corr_h-1= '$LITTLE_H_1_CORR			>> $CONFIGFILE
				echo 'unit_corr_1e3= '$e3_CORR				>> $CONFIGFILE
				echo 'unit_corr_h-2= '$LITTLE_H_2_CORR			>> $CONFIGFILE

				set i=1
				foreach item ($col_name_array) 	

					echo 'col_'$item'= '"$col_array[$i]"  		>> $CONFIGFILE
					@ i++
				end

				set i=1
				foreach item ($halo_col_name_array) 	

					echo 'halo_col_'$item'= '"$halo_col_array[$i]"  >> $CONFIGFILE
					@ i++
				end
				set i=1
				#foreach item ($tarsel_col_name_array) 	

					#echo 'tarsel_col_'$item'= '"$tarsel_col_array[$i]"  >> $CONFIGFILE
					#@ i++
				#end
				#Take care of proper BOXSIZE:
				echo 'box_size= '$BOX_SIZE  >> $CONFIGFILE
				if ($NAME =~ *125*) 		echo 'box_size= 125.0'  >> $CONFIGFILE
				if ($NAME =~ *1Gpc*) 		echo 'box_size= 1000.0' >> $CONFIGFILE
				if ($NAME =~ *62*) 		echo 'box_size= 62.5'  	>> $CONFIGFILE
				if ($NAME =~ *400*) 		echo 'box_size= 400.0' 	>> $CONFIGFILE
				if ($NAME =~ *sub) 		echo 'box_size= 112.0' 	>> $CONFIGFILE
				
				if ($SNAPIDZ_MAPPING == 'True') then

					set INPUTFILENAME_PART1 = `ls $MODEL/*$SNAPID* | head -1 | awk -F$SNAPID '{print $1}'`		
					set INPUTFILENAME_PART2 = `ls $MODEL/*$SNAPID* | tail -1 | awk -F$SNAPID '{print $2}'`

					if ($INPUTFILENAME_PART1 =~ *'/') then 
						set INPUTFILENAME_PART1 = ''
					else
						set INPUTFILENAME_PART1 = `basename $INPUTFILENAME_PART1`
					endif
				else

					if ($CUSTOM_SNAPNAME == 'crossmatch') then
						set INPUTFILENAME_PART1 = $SNAPNAME
					else if ($CUSTOM_SNAPNAME == 'hdf5' && $NAME =~ SAG_* && $PATH_TO_CATALOG$NAME'/' == $SOFTLINK_TO_DATA || $NAME =~ *DR*) then
						set INPUTFILENAME_PART1 = `basename $key`
					else
						set INPUTFILENAME_PART1 = $key
					endif

					set INPUTFILENAME_PART2 = ''
				endif				

				#Take care of proper naming SNAPIDZRED.txt files:
				if (! -e $PATH_TO_CATALOG$NAME/snapidzred.txt) then

					if ($SNAPNAME =~ *hdf5) then
						echo 'HDF5: no snapidzred.txt file existing'
						set INPUTFILENAME_PART1 = $INPUTFILENAME_PART1 
						set SNAPID = $MYSNAPS
						echo $SNAPID, $INPUTFILENAME_PART1

					#else if ($SNAPNAME =~ *$CUSTOM_SNAPNAME* ) then
						#echo 'MDPL: no snapidzred.txt file existing'
						#set INPUTFILENAME_PART1 = $INPUTFILENAME_PART1 
						#set INPUTFILENAME_PART2 = `basename $INPUTFILENAME_PART2`
						#echo $SNAPID, $INPUTFILENAME_PART1

					else if (-e `ls $PATH_TO_CATALOG$NAME/*.lst`) then
						cp $PATH_TO_CATALOG$NAME/*.lst $PATH_TO_CATALOG$NAME/snapidzred.txt
						echo copied $PATH_TO_CATALOG$NAME/*.lst to $PATH_TO_CATALOG$NAME/snapidzred.txt
					
					else if (-e `ls $PATH_TO_CATALOG$NAME/suss*`) then
						cp $PATH_TO_CATALOG$NAME/suss* $PATH_TO_CATALOG$NAME/snapidzred.txt
						echo copied $PATH_TO_CATALOG$NAME/suss* to $PATH_TO_CATALOG$NAME/snapidzred.txt			
					else 
						#echo 'NO MAPPING'
						set SNAPID = `echo $MYSNAPS` 
					endif
				endif
				echo 'unit_code= '$UNIT_CODE	>> $CONFIGFILE
				echo 'plot_cum= '$PLOT_CUM	>> $CONFIGFILE
				echo 'use_store_register= '$use_store_register	>> $CONFIGFILE					
				echo 'inputfilename_part1= '$INPUTFILENAME_PART1  		>> $CONFIGFILE
				echo 'snapid= '$SNAP_LIST		 			>> $CONFIGFILE
				echo 'inputfilename_part2= '$INPUTFILENAME_PART2  		>> $CONFIGFILE
				echo 'softlink_directory_code= '$SOFTLINK_DIRECTORY_CODE  	>> $CONFIGFILE
				echo 'hubble_par= '$SIM_LITTLE_H	 		>> $CONFIGFILE
				echo 'skycoverage= '$skycoverage	 		>> $CONFIGFILE
				echo 'zmin= '$zmin	 				>> $CONFIGFILE
				echo 'zmax= '$zmax	 				>> $CONFIGFILE
				echo 'name_weights= '$col_name_weights	 				>> $CONFIGFILE

				set nr_of_cols_to_load = 0
				set c=1
				foreach item ($method_array)
					set nr_of_cols_to_load = `expr $nr_of_cols_to_load + $item`
					@ c++
				end

				#read in specifics
				echo 'nr_of_cols_to_load= '$nr_of_cols_to_load 	 	>> $CONFIGFILE
				echo 'create_subcat= '$create_subcat			>> $CONFIGFILE
				echo 'convert_sph_to_cart_coords= '$convert_sph_to_cart_coords	>> $CONFIGFILE

				echo 'use_snapidz_mapping= '$SNAPIDZ_MAPPING   	 	>> $CONFIGFILE
				echo 'manual_input_redshift= '$MANUAL_REDSHIFT_INPUT	>> $CONFIGFILE
				echo 'load_from_file= '$LOAD_FROM_FILE	   		>> $CONFIGFILE
				echo 'sample_info= '$sample_info	   	 	>> $CONFIGFILE
				echo 'ngal_random_sample= '$ngalaxies_random		>> $CONFIGFILE
				echo 'fits_map_names= '$fits_map_names			>> $CONFIGFILE
				echo 'fits_map_names_mapping= '$fits_map_names_mapping	>> $CONFIGFILE
				echo 'UNIT_CODE= '$UNIT_CODE		   	 	>> $CONFIGFILE
				echo 'cosmology= '$cosmology		   	 	>> $CONFIGFILE
				echo 'my_annotation= '$my_annotation	   	 	>> $CONFIGFILE

				echo 'filter_density= '$filter_density	>> $CONFIGFILE
				echo 'filter_density_ngal= '$filter_density_ngal	>> $CONFIGFILE
				echo 'filter_density_cut_name= '$filter_density_cut_name	>> $CONFIGFILE
				echo 'filter_density_select_highest= '$filter_density_select_highest	>> $CONFIGFILE
				echo 'filter_density_write_coordinates= '$filter_density_write_coordinates	>> $CONFIGFILE
				echo 'histo_norm_y= '$histo_norm_y	>> $CONFIGFILE

				echo 'plot_all= '$plot_all		>> $CONFIGFILE
				echo 'plot_single= '$plot_single	>> $CONFIGFILE
				echo 'single_mix= '$plot_single_mix	>> $CONFIGFILE
				echo 'subplots= '$subplots		>> $CONFIGFILE

				if ($NAME == Galacticus_1Gpc) then
					set tarsel_code = $tarsel_code_Galacticus
				else if ($NAME =~ Galacticus_1Gpc*) then
					set tarsel_code = $tarsel_code_Galacticus_CMASS
				else if ($NAME == SAG_1Gpc) then
					set tarsel_code = $tarsel_code_SAG
				else if ($NAME =~ *125) then
					set tarsel_code = $tarsel_code_SAG_125
				else if ($NAME == SAGE_1Gpc) then
					set tarsel_code = $tarsel_code_SAGE
				else
					set tarsel_code = $tarsel_code_default
				endif

				echo 'tarsel_code= '$tarsel_code			>> $CONFIGFILE
				echo 'simulation_name= '$simulation_name		>> $CONFIGFILE
				echo 'telescope_name= '$telescope_name			>> $CONFIGFILE
				echo 'load_subcat= '$load_subcat			>> $CONFIGFILE

				if ($NAME == Galacticus_1Gpc) then
					set output_filename_code = $output_filename_code_Galacticus
				else if ($NAME == Galacticus_1Gpc_CMASS) then
					set output_filename_code = $output_filename_code_Galacticus_CMASS
				else if ($NAME == SAG_1Gpc) then
					set output_filename_code = $output_filename_code_SAG
				else if ($NAME == SAGE_1Gpc) then
					set output_filename_code = $output_filename_code_SAGE
				else if ($NAME == SAGE_1Gpc_CMASS) then
					set output_filename_code = $output_filename_code_SAGE_CMASS
				else if ($NAME == SAGE_1Gpc_CMASS_centrals) then
					set output_filename_code = $output_filename_code_SAGE_CMASS_centrals
				else if ($NAME == SAGE_1Gpc_CMASS_MHALO) then
					set output_filename_code = $output_filename_code_SAGE_CMASS_MHALO
				else
					set output_filename_code = $output_filename_code_default
				endif
				echo 'change_IMF= '$change_IMF			 	>> $CONFIGFILE
				echo 'which_kcorrect= '$which_kcorrect			>> $CONFIGFILE

				if ($change_IMF != 'False') then
					if ($output_filename_code == '') then
						set output_filename_code = $output_filename_code$change_IMF
					else
						set output_filename_code = $change_IMF
					endif
				endif

				echo 'output_filename_code= '$output_filename_code 	>> $CONFIGFILE
				echo 'set_register_name= '$set_register_name	 	>> $CONFIGFILE
 				echo 'selection_name= '$selection_name			>> $CONFIGFILE
				echo 'start_fileID= '$start_fileID			>> $CONFIGFILE
				echo 'nr_files_snapshot= '$nr_files_snapshot		>> $CONFIGFILE
				echo 'end_fileID= '$end_fileID				>> $CONFIGFILE
				echo 'skip_reading_data= '$skip_reading_data		>> $CONFIGFILE	
				#twoPCF specifics
				echo 'CF_NR_bins= '$CF_NR_bins				>> $CONFIGFILE
				echo 'CF_R_max= '$CF_R_max				>> $CONFIGFILE
				echo 'choose_CUTE= '$choose_CUTE			>> $CONFIGFILE
				echo 'corr_estimator= '$corr_estimator			>> $CONFIGFILE
				echo 'CF_box_size= '$CF_box_size			>> $CONFIGFILE
				#analyseTargetSelection specifics
				echo 'filterSample= '$filterSample 			>> $CONFIGFILE
				echo 'calcHistoSample= '$calcHistoSample		>> $CONFIGFILE
				echo 'plotSample= '$plotSample				>> $CONFIGFILE

				echo 'OUTPUT_ASCII= '$HDF5_TO_ASCII				>> $CONFIGFILE
				echo 'filter_halomass_Sergio= '$filter_halomass_Sergio		>> $CONFIGFILE

				echo 'twoPCF_path_to_data= '$twoPCF_path_to_data		>> $CONFIGFILE
				echo 'twoPCF_which= '$twoPCF_which				>> $CONFIGFILE
				echo 'twoPCF_calculate= '$twoPCF_calculate			>> $CONFIGFILE
				echo 'twoPCF_pimax= '$twoPCF_pimax				>> $CONFIGFILE
				echo 'twoPCF_nthreads= '$twoPCF_nthreads			>> $CONFIGFILE
				echo 'twoPCF_rmin= '$twoPCF_rmin				>> $CONFIGFILE
				echo 'twoPCF_rmax= '$twoPCF_rmax				>> $CONFIGFILE
				echo 'twoPCF_nbins= '$twoPCF_nbins				>> $CONFIGFILE
				echo 'twoPCF_ngal_random_sample= '$twoPCF_ngal_random_sample	>> $CONFIGFILE

				echo 'selectRegion= '$selectRegion			>> $CONFIGFILE
				echo 'selectRegion_volume_type= '$selectRegion_volume_type >> $CONFIGFILE
				echo 'selectRegion_filename= '$selectRegion_filename	>> $CONFIGFILE
				echo 'selectRegion_path_to_data= '$selectRegion_path_to_data	>> $CONFIGFILE
				echo 'selectRegion_output_filename= '$selectRegion_output_filename >> $CONFIGFILE
				echo 'selectRegion_region_name= '$selectRegion_region_name	>> $CONFIGFILE
				echo 'selectRegion_units= '$selectRegion_units	>> $CONFIGFILE
				echo 'selectRegion_col_id_x_pos= '$selectRegion_col_id_x_pos	>> $CONFIGFILE
				echo 'selectRegion_col_id_y_pos= '$selectRegion_col_id_y_pos	>> $CONFIGFILE
				echo 'selectRegion_col_id_z_pos= '$selectRegion_col_id_z_pos	>> $CONFIGFILE
				echo 'selectRegion_col_id_radius= '$selectRegion_col_id_radius	>> $CONFIGFILE
				echo 'selectRegion_periodic_conds= '$selectRegion_periodic_conds >> $CONFIGFILE

				echo 'createProgenitorList= '$createProgenitorList			>> $CONFIGFILE
				echo 'createProgenitorList_path_to_data= '$createProgenitorList_path_to_data	>> $CONFIGFILE
				echo 'createProgenitorList_output_filename= '$createProgenitorList_output_filename >> $CONFIGFILE
				set c=1
				foreach item ($analyse_tarsel_array)
					if ($item == 1) then
						#echo $item "$analyse_tarsel_name_array[$c]"
						echo "$analyse_tarsel_name_array[$c]"'= '"$analyse_tarsel_array[$c]" >> $CONFIGFILE
					endif
					@ c++
				end

			endif

		end # CATALOG

	echo '##########################################'
	echo ' '

	@ a++
	end # ISNAP
    @ j++
end # Model

set CONTROL_CATS 	= `expr $BOX_COUNT_62 + $BOX_COUNT_125`
set CONTROL_CATS 	= `expr $CALI_COUNT_uc + $CALI_COUNT`
if ($CONTROL_CATS == 2) set CALI_CODE = ' '

set c=1
set i=1
set a=1

if ($tarsel_code_Galacticus == '') set tarsel_code_Galacticus = '-'
if ($tarsel_code_SAG == '') set tarsel_code_SAG = '-'
if ($tarsel_code_default == '') set tarsel_code_default = '-'

foreach item ($method_array)
	#echo "$method_name_array[$c]" $item
	if ("$method_name_array[$c]" == plotOnly && $item == 1) then
		foreach item2 ($total_col_name_array)

			if ("$total_col_array[$i]" != 99) then
				#echo 1291 $item2 "$total_col_array[$i]" "$method_name_array[$c]" $tarsel_code "$plot_plotOnly_config_array[$a]"
				
				$MAIN_PATH'myRun/plot_config.sh' "plotOnly" "$plot_plotOnly_config_array[$a]" $BOX_SIZE $CALI_CODE $item2 $CAT_NAME_IN_PLOT0 $CAT_NAME_IN_PLOT1 $REDSHIFT_PLOT0 $REDSHIFT_PLOT1 $SOFTLINK_CODE $a $tarsel_code_Galacticus $tarsel_code_SAG $MATPLOT_LINESTYLE $MATPLOT_COL $MATPLOT_MARKERSTYLE $MATPLOT_MARKERCOL
				@ a++
			endif
			@ i++
				
		end

	else if ("$method_name_array[$c]" == analyseTargetSelection && $item == 1) then
			foreach item2 ($analyse_tarsel_name_array)

				if ("$analyse_tarsel_array[$i]" == 1) then
					#echo 1421 $item2 "$analyse_tarsel_name_array[$i]"  $tarsel_code "$plot_analyseTargetSelection_config_array[$a]"
				
				$MAIN_PATH'myRun/plot_config.sh' "$analyse_tarsel_name_array[$i]" "$plot_analyseTargetSelection_config_array[$a]" $BOX_SIZE $CALI_CODE $CAT_NAME_IN_PLOT0 $CAT_NAME_IN_PLOT1 $REDSHIFT_PLOT0 $REDSHIFT_PLOT1 $SOFTLINK_CODE $c $tarsel_code_Galacticus $tarsel_code_SAG $MATPLOT_LINESTYLE $MATPLOT_COL $MATPLOT_MARKERSTYLE $MATPLOT_MARKERCOL $MATPLOT_COLORMAP $MATPLOT_COLORBAR $MATPLOT_PLOTTYPE
				@ a++
			endif
			@ i++
				
		end

	else 

		if ($item == 1)  then

			$MAIN_PATH'myRun/plot_config.sh' "$method_name_array[$c]" "$plot_config_array[$c]" $BOX_SIZE $CALI_CODE $CAT_NAME_IN_PLOT0 $CAT_NAME_IN_PLOT1 $REDSHIFT_PLOT0 $REDSHIFT_PLOT1 $SOFTLINK_CODE $c $tarsel_code_Galacticus $tarsel_code_SAG $MATPLOT_LINESTYLE $MATPLOT_COL $MATPLOT_MARKERSTYLE $MATPLOT_MARKERCOL $MATPLOT_COLORMAP $MATPLOT_COLORBAR $MATPLOT_PLOTTYPE
		endif
	
	endif
	@ c++
end


echo 'GENERATING CONFIG FILES ... DONE!'
echo ' '
