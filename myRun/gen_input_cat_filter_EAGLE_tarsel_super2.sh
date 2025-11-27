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
set simulation_name = 'MDPL2'
set simulation_name = 'Gadget3'

#Set the number of particles as x for x^3 for total numbers of particles in simulation: e.g. 1504 for EAGLE 1504^3!
#+++++++++++++++++++++++++++++++++++++++++++
set resolution = '1504'

#Set the name of the telescope: e.g. SDSS, etc.!
#+++++++++++++++++++++++++++++++++++++++++++
set telescope_name = 'SDSS'

#Set the name of the cosmology used in the simulation: e.g. Planck, etc.!
#+++++++++++++++++++++++++++++++++++++++++++
set cosmology = 'Planck'

#Set the name of the cosmology used in the simulation: e.g. Planck, etc.!
#+++++++++++++++++++++++++++++++++++++++++++
set IMF = 'Chabrier03'

#Set any note to they output-file!
#+++++++++++++++++++++++++++++++++++++++++++
set my_annotation = 'catalog from Patricia (with and without rhalf_gas) crossmatched with them from Yetli with envr info'
set my_annotation = 'EAGLE data base centrals catalog crossmatched with various from Patricia Yetli Silvio including effective radius of the gas and angular momentum of gas and stars'
set my_annotation = 'catalog from Yetli with envr'
set my_annotation = 'EAGLE data base properties crossmatched with Patricia3 - for Silvio - Yetli ev catalogs -- EAGLE-super'

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
set tarsel_code_default = 'REF100_generaltable_z0-header_crossmatched'
#set tarsel_code_default = 'REF25_generaltable_z0'
set tarsel_code_default = 'centrals_super_corr_ropt2'
#set tarsel_code_default = 'tree_roots'
#set tarsel_code_default = 'tree_roots_history_assembly'
#set tarsel_code_default = 'centrals_crossmatched_bins_regplot'
#set tarsel_code_default = ''
#HF-PaperII_dens6'

set tarsel_code_Galacticus = ''
set tarsel_code_SAG = ''
set tarsel_code_SAGE = ''

#Choose different catalogues to read!
#+++++++++++++++++++++++++++++++++++++++++++

set add_tarsel_code = ''


set output_filename_code_Galacticus = ''
set output_filename_code_SAG = ''
set output_filename_code_SAGE = ''

#NOTE standard: [plot_key]_[catname+boxsize]_z_[redshift(0.2f)].[fileformat]
#set output_filename_code_default = 'SFH_300_main_cents'
#set output_filename_code_default = 'crossmatched_mstar_gt_1e9_existing'
set output_filename_code_default = 'parents_dens'
#set output_filename_code_default = 'dens7'
#set output_filename_code_default = 'super'
set output_filename_code_default = 'kine'

#Should the header of an filtered ASCII-output file be written so that Topcat uses automatically the names to name the columns? name1(1) name2(2)
set set_header_Topcat_format = 'yes'

#CHANGE IMF:
set change_IMF = 'Kroupa_Croton+16'
set change_IMF = 'False'
#give an additional identifier to the filename which is generated automatically --> use the same output_filename_code to read the histo-files in again!
#+++++++++++++++++++++++++++++++++++++++++++

#name of the register you want to write/read data, if default register is not desired! Change 'default' to 'your name'
set set_register_name = False
#++++++++++++++++++++++++++++++++++++++++++

#Set file job-id numbers which should be selected: NOTE: only Galacticus job0-job103==1Snapshot/redshift
#+++++++++++++++++++++++++++++++++++++++++++
set start_fileID 	= 0
set nr_files_snapshot 	= 1
set end_fileID 	= False

#choose a selection name which describes the current selection you want to apply
#+++++++++++++++++++++++++++++++++++++++++++
set selection_name = False

set plot_custom_array = ('')

#Enter 'SAM' for reading a SAM catalouge or anything else to read a random HDF5 file
#+++++++++++++++++++++++++++++++++++++++++++

#NOTE: Load subcat=True to load the extracted subcatalogues of the SAM!
if ($HOME == /home/claudia || $HOME == /home/doris) then
	set HOME_MODE = True
	set HDF5_read_code = 'HDF5'
	set load_subcat = False
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
set create_subcat = False
set convert_sph_to_cart_coords = False
set PLOT_CUM = 'Truee'
set set_home_mode = Truee

set PLOT_CUSTOM_LOOP = 'False'
set plot_custom_array = ()

set MANUAL_REDSHIFT_INPUT_default = 0.0

if ($HOME_MODE == $set_home_mode) then
	set SOFTLINK_CODE = A
	set plot_all = 'True'			
	set LOAD_FROM_FILE = 'True'

else
	set SOFTLINK_CODE = E
	set plot_all = 'False'		
	set LOAD_FROM_FILE = 'False'
endif

#Choose the Snapshots/redshifts --> see list for the SAMS below OR enter a key which is representative for the catalog e.g. 'Gal' for 'Galacticus_1Gpc.hdf5'
#+++++++++++++++++++++++++++++++++++++++++++
set MYSNAPS_62    			= '0.000'
set MYSNAPS_100Mpc			= '28'
set MYSNAPS_125				= '107' 
set MYSNAPS_400Mpc			= '96'
set MYSNAPS_500Mpc			= 'Guo2013_MSI_sn45'

#redshifts/snapnr Galacticus
#set MYSNAPS_1Gpc			= '79 75 70 67 65 60 57 52 50 47 45 42 40 38 36 34 32 30 28 26 24 22 20 18 16 14 13 12 11 10 9 8 7 6 5 4 3 2 1'

#redshifts/snapnr SAG
#set MYSNAPS_1Gpc			= '125 124 122 121 119 115 110 107 104 101 98 95 92 89 86 83 80 77 75 70 67 65 60 57 52 50 47 45 42 40 38 36 34 32 30 28 26 24 20'

#redshifts/snapnr SAGE
#set MYSNAPS_1Gpc			= '0.000 0.093 0.142 0.490 0.592 0.740 0.901 1.220 2.095 3.037 4.038 5.017 6.022 7.026 8.372'
#set MYSNAPS_1Gpc			= '0.000 0.093 0.142 0.490 0.557 0.740 0.859 1.077 1.270 1.425 1.650 1.896 2.095 2.382 2.614 3.037 3.411 3.929 4.038 4.385 4.627 4.754 4.882 5.017 5.289 5.720 5.873 6.022 6.184 6.342 6.508 6.849 7.026 7.203 7.389 7.764 7.961 8.166 8.372'


set MYSNAPS_default			= '28'

#redshifts of commonly used models!
#Galacticus SMPL 400 Mpc/h: 96 z=0.0, 92 z=0.1, 85 z=0.14, 47 z=0.55'
#Galacticus MDPL2 1 Gpc/h: 79 z=0, 75 z=0.09, 73 z=0.14, 61 z=0.49, 59 z=0.56, 55=0.7, 52=.82, 50 =0.9, 49=0.94, 29 z=2.03, 24 z=2.38, min 1 z=8.0
#SAG MDPL2 1 Gpc/h: 	125 z=0,  122 z=0.07, 119 z=0.14, 107=0.49, 104 z=0.59, 101=0.7, 95=0.94, min 75 z=2.03
#SAG cal01 nifty 125Mpc/h: 	107 z=0,  98 z=0.15, min 4 z=7.73
#EAGLE 100 Mpc/h:     	28 z=0

if ($HDF5_read_code == 'SAMHDF5' || $HDF5_read_code == 'BINARY_SAGE') then

	if ($SOFTLINK_CODE == A || $SOFTLINK_CODE =~ A*) then 
		#SFR_OII
		set MYSNAPS_1Gpc	= '59'
# 75 73 58 55 49 43 29'

		#cSFRF
		#set MYSNAPS_1Gpc	= '79 75 70 67 65 60 57 52 50 47 45 42 40 38 36 34 32 30 28 26 24 22 20 18 16 14 13 12 11 10 9 8 7 6 5 4 3 2 1'

		#SFH 300
		set MYSNAPS_1Gpc	= '79 78 77 76 75 74 73 72 71 70 69 68 67 66 65 64 63 62 61 60 59 58 57 56 55 54 53 52 51 50 49 48 47 46 45 44 43 42 41 40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5'
		#all Galacticus
		#set MYSNAPS_1Gpc	= '59 58 57 56 55'
# 54 53 52 51 50 49 48 47 46 45 44 43 42 41 40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5'
		#set MYSNAPS_1Gpc	= '59 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1'

		#Galacticus 400 Mpc/h snapshot list for sfr2z
		#'96 94 92 90 87 85 84 80 75 70 65 50 47 45 42 40 38 36 30 28 26 24 22 20 18 16 14 12 10 8 6 5 4 3 2 1'

	else if ($SOFTLINK_CODE == C) then
		#Galacticus 400Mpc (SMDPL)
		set MYSNAPS_400Mpc	= '47 46 45 44 43 42 42 40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1'

	else if ($SOFTLINK_CODE == B || $SOFTLINK_CODE == I) then
		#set MYSNAPS_1Gpc	= '125 92 89 86 83 80 77 75'

		#redshifts  SAGv2 
		set MYSNAPS_1Gpc	= '125 122 119 116 113 110 107 104 101 98 95 92 89 86 83 80 77 75'

		#redshifts  SAGv3
		#set MYSNAPS_default	= '125 124 122 121 119 115 110 107 104 101 98 95 92 89 86 83 80 77 75 70 67 65 60 57 52 50 47 45 42 40 38 36 34 32 30 28 26 24 20'
	
		#SAG calibration01 nifty 125 Mpc/h snapshot list for sfr2z
		#107 98 92 80 70 61 47 35 27 21 17 14 12 9 6 4

	else

		set MYSNAPS_1Gpc	= '0.557'
		#set MYSNAPS_1Gpc	= '0.000 0.093 0.142 0.490 0.557 0.740 0.859 1.077 1.270 1.425 1.650 1.896 2.095 2.382 2.614 3.037 3.411 3.929 4.038 4.385 4.627 4.754 4.882 5.017 5.289 5.720 5.873 6.022 6.184 6.342 6.508 6.849 7.026 7.203 7.389 7.764 7.961 8.166 8.372'

	endif
else
	if ($SOFTLINK_CODE == A) then 
		set MYSNAPS_1Gpc	= '29'
	else if ($SOFTLINK_CODE == B) then 
		set MYSNAPS_1Gpc	= '121'
	else if ($SOFTLINK_CODE == I) then 
		set MYSNAPS_1Gpc	= '75'
	else
		set MYSNAPS_1Gpc	= '2.028'
	endif	
endif


#For extracting properties from observational catalogs
#+++++++++++++++++++++++++++++++++++++++++++
#CMASS catalog
#set fits_map_names	= 'DEC_1;RA_1;Z_1;LOGMASS;AGE;MAGSCALED'
#set fits_map_names_mapping = 'DEC;RA;Z;mstar;age;mAB_dA_total_u;mAB_dA_total_g;mAB_dA_total_r;mAB_dA_total_i;mAB_dA_total_z;'

#Portsmouth starforming Salpeter
set fits_map_names	= 'DEC;RA;Z;LOGMASS;AGE;SFR;MAGSCALED'
set fits_map_names_mapping = 'DEC;RA;Z;mstar;age;sfr;mAB_dA_total_u;mAB_dA_total_g;mAB_dA_total_r;mAB_dA_total_i;mAB_dA_total_z;'

#SDSS DR7 MPA-JHU
set fits_map_names	= '0;1;2;3'
set fits_map_names_mapping = 'DEC;RA;Z;mstar;'

#Observational catalog details
#CMASS DR12
set skycoverage = 9376
set zmin	= 0.5
set zmax	= 0.6

#Give infos about the sample of galaxies selected: 'full': no cuts what so ever, or standard cuts. 'else': a subsample will be selected --> choose number of galaxies which should be selected!
#+++++++++++++++++++++++++++++++++++++++++++
set sample_info	= 'full'

#Select ngalaxies for random catalogue:
set ngalaxies_random = 'False'



#Select 'skip_reading_data' if for example with 2PCF the hdf5 should not be read or is not exsiting!
#+++++++++++++++++++++++++++++++++++++++++++
set skip_reading_data = yess
set use_store_register 	= yees


#Put if 'True', put number of galaxies to select e.g. 'filter_density_ngal = 340000' for number density of CMASS 3.4x10-4 Mpc-3 h3 most
#+++++++++++++++++++++++++++++++++++++++++++
set filter_density = 'False'
set filter_density_ngal = 340000
set filter_density_cut_name = densCut_CMASS

set filter_density_select_highest = True
set filter_density_write_coordinates = False


#Any additional outputs?
#+++++++++++++++++++++++++++++++++++++++++++
set ASCII_TO_HDF5 = False
set HDF5_TO_ASCII = True

#Use k-corrections to correct luminosities and magnitudes: 'approx' is KCORRECT from Chiligarian
set which_kcorrect = 'approx'

#specify details for 'twoPCF' using CORRFUNC (Manodeep Sinha)!
#+++++++++++++++++++++++++++++++++++++++++++
#set twoPCF_path_to_data	 = '/store/erebos/doris/'
set twoPCF_path_to_data	 = '/z/doris/anaconda/pro/data/LGALAXIES/'
set twoPCF_which	 = 'BOX'
set twoPCF_calculate	 = 'xi'
set twoPCF_pimax 	 = '60'
set twoPCF_nthreads 	 = '8'
set twoPCF_rmin 	 = '0.5'
set twoPCF_rmax		 = '200'
set twoPCF_nbins	 = '25'

set twoPCF_ngal_random_sample = 'False'

#specify details for 'selectRegion' (e.g. TheThreeHundred Cluster Project)!
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

set calc_fast_histo		= True

set plotXY				= 0
set plotOnly			= 0
set filterData			= 1
set analyseTargetSelection	= 0

set SMF 				= 0
set SFRF 				= 0
set sSFRF				= 0
set cSFRD				= 0
set sfr2z				= 0
set sfr2mstar			= 0
set ssfr2mstar			= 0
set oh2mstar			= 0
set HMF 				= 0
set mstar2mhalo			= 0
set mstar2mhalofunc		= 0
set ngal2mhalo 			= 0
set HMF_no				= 0
set mstar2mhalo_no		= 0
set mstar2mhalofunc_no		= 0
set mstar2mhalovsSFR		= 0
set ngal2mhalo_no			= 0
set zgas2mstar			= 0
set zgas2mcold			= 0
set mcold2mstar			= 0
set mbh2mstarsph			= 0
set twoPCF				= 0
set HOD				= 0
set mstar2rhalf			= 0

#Set manually redshift if home mode or load subcat are activated:
#+++++++++++++++++++++++++++++++++++++++++++

set histo_norm_y = False

#Default settings:
if ($HOME_MODE == $set_home_mode || $load_subcat == True) then

	if ($SMF == 1 || $HMF == 1 || $HMF_no == 1) then
		set SAG_REDSHIFT_INPUT = 0.09
		set SAG_125_REDSHIFT_INPUT = 1.22
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 0.56
		set GALACTICUS_400Mpc_REDSHIFT_INPUT = 0.1
		set SAGE_REDSHIFT_INPUT = 0.56
		set SAG_v2_REDSHIFT_INPUT = 0.59

	else if ($SFRF == 1) then
		set SAG_125_REDSHIFT_INPUT = 0.15
		set SAG_REDSHIFT_INPUT = 0.14
		set SAG_v2_REDSHIFT_INPUT = 0.07
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 0.14
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
		set SAG_v2_REDSHIFT_INPUT = 0.59
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
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 0.56
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
		set GALACTICUS_1Gpc_REDSHIFT_INPUT = 0.56
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

#analyseTargetSelection: Specify the plot you want to make
if ($analyseTargetSelection == 1) then
	set filterSample 		= 'False'
	set calcHistoSample 		= 'False'
	if ($calcHistoSample == 'True') then
		set analyse_tarsel_histo	= 1
		set plotSample			= 'False'
	else
		set analyse_tarsel_histo	= 0
		set plotSample			= 'True'

	endif
else
	set filterSample 		= 'False'
	set calcHistoSample 		= 'False'
	set analyse_tarsel_histo	= 0
	set plotSample			= 'False'
endif

set analyse_tarsel_dmesa_i	= 0
set analyse_tarsel_dmesa_mstar	= 0
set analyse_tarsel_gminr_mstar 	= 0
set analyse_tarsel_rmini_mstar	= 0
set analyse_tarsel_gmini_mstar	= 0
set analyse_tarsel_uminr_mstar	= 0
set analyse_tarsel_gminz_mstar	= 0
set analyse_tarsel_rmini_i	= 0
set analyse_tarsel_gmini_i	= 0
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
set analyse_tarsel_ssfr_mstar	= 0
set analyse_tarsel_mbh_mstar_spheroid	= 0
set analyse_tarsel_mcold_mstar	= 0
set analyse_tarsel_Mzgas_mstar	= 0
set analyse_tarsel_Mzgas_mcold	= 0
set analyse_tarsel_gminr_r	= 0
set analyse_tarsel_uminr_r	= 0

set analyse_tarsel_gminr_uming	= 0
set analyse_tarsel_rmini_gminr	= 0

set analyse_tarsel_rdisk_rbulge	= 0
set analyse_tarsel_rhalfmass_mstar	= 0
set analyse_tarsel_rhalfmass_rbulge	= 0
set analyse_tarsel_rhalfmass_rdisk	= 0
set analyse_tarsel_rhalfdisk_mstar	= 0
set analyse_tarsel_rhalfbulge_mstar	= 0
set analyse_tarsel_rdisk_mstar	= 0
set analyse_tarsel_spinParameter_mstar	= 0
set analyse_tarsel_spinParameter_Mzstar	= 0
set analyse_tarsel_Mzstar_mstar	= 0
set analyse_tarsel_dmesa_ssfr	= 0
set analyse_tarsel_ssfr_sfr	= 0
set analyse_tarsel_spinParameter_dmesa	= 0
set analyse_tarsel_gminr_ssfr	= 0
set analyse_tarsel_NFW_mhalo200c	= 0
set analyse_tarsel_ssfr_mhalo	= 0
set analyse_tarsel_ssfr_i	= 0
set analyse_tarsel_ssfr_r	= 0
set analyse_tarsel_ssfr_mcold	= 0
set analyse_tarsel_zcold_mhalo	= 0
set analyse_tarsel_zcold_sfr	= 0
set analyse_tarsel_ssfr_zcold	= 0
set analyse_tarsel_mstar_mhalo	= 0

#Choose variables to compute/plot/process:
#+++++++++++++++++++++++++++++++++++++++++++

#99='DO NOT CHOOSE COLUMN', others=col_id
#COLUMN-CODE IN HALO CATALOGUE
set halo_nhalos		= 999
set halo_haloid		= 999
set halo_hostid		= 999
set halo_desIndex	= 999
set halo_nodemass	= 999
set halo_vmax		= 999
set halo_rsca		= 999
set halo_angM		= 999
set halo_spin		= 999

set halo_redshift	= 999

set halo_x_pos	= 999
set halo_y_pos	= 999
set halo_z_pos	= 999

set halo_x_vel	= 999
set halo_y_vel	= 999
set halo_z_vel	= 999

set halo_x_vel_disp	= 999
set halo_y_vel_disp	= 999
set halo_z_vel_disp	= 999

###############################################################################
#	IDNUMBERS
###############################################################################

set haloid		= 1		
set hostid		= 999
set galaxyID		= 3
set jsub 		= 5

set nodeIndex		= 999
set parentIndex		= 999
set satelliteNodeIndex	= 999
set satelliteIndex	= 999
set siblingIndex	= 999

set fofID		= 0
set progFofID		= 999
set mainLeafID		= 999
set topLeafID		= 999
set firstProgenitorID	= 999
set lastProgID		= 999

###############################################################################
#	HALOS
###############################################################################
	
set redshift 		= 999
	
set mhalo		= 999
#only for Galacticus
set mhalo_sat		= 999

set mhalo_30kpc		= 999
set mhalo_200c		= 4
set mhalo_cents		= 999
set mhalo_cents_200c	= 999
set mbasic		= 999
set mbasic_200c		= 999
set mhalo_fof		= 999
set mhalo_50kpc		= 999
set mhalo_70kpc		= 999

#Standard properties
set ngalaxies		= 999
set nsats		= 999
set nSubhalos		= 999
set npros		= 999

set timeLastIsolated	= 999
set satelliteMergeTime	= 999

###############################################################################
#	MERGER
###############################################################################

set lastMerger		= 84
set lastMinorM		= 999
set lastMajorM		= 999
set ratio_lastM 	= 85

###############################################################################
#	POSITIONS
###############################################################################

set x_pos		= 999
set y_pos		= 999
set z_pos		= 999

set x_pos_subhalo		= 7
set y_pos_subhalo		= 8
set z_pos_subhalo		= 9

set x_pos_cof		= 999
set y_pos_cof		= 999
set z_pos_cof		= 999


###############################################################################
#	STRUCTURE & KINEMATICS
###############################################################################

set NFW_con		= 13

set x_vel		= 25
set y_vel		= 26
set z_vel		= 27

set vmax		= 41

set vbulge		= 999
set vdisk		= 999

set vpeak		= 999
set vdisp 		= 999
set vdisp_30kpc	= 14
set vdisp_50kpc = 999
set vdisp_70kpc = 999

set vpec_norm 	= 999

#EAGLE Silvio
set vdisp_r502D_edgeOn	= 15
set vdisp_r502D_random	= 999
set vdisp_2r502D_edgeOn	= 16
set vdisp_2r502D_random	= 999

set delta_vdisp_edgeOn	= 24

set VvS_r502D_edgeOn	= 17
set VvS_r502D_random	= 999
set VvS_2r502D_edgeOn	= 18
set VvS_2r502D_random	= 999

set delta_VvS_edgeOn	= 69

set ell_r502D_edgeOn	= 19
set ell_r502D_random	= 999
set ell_2r502D_edgeOn	= 20
set ell_2r502D_random	= 999

set Sersic_n = 21

set ell_30kpc_T19 = 10
set ellDM_30kpc_T19 = 11
set vdisp_30kpc_T19 = 12
set Vrot2Vdisp_30kpc_T19 = 89
set kappaCoRot_30kpc_T19 = 90
set medOrbCirc_30kpc_T19 = 28
set triax_30kpc_T19 = 33

set v200c	= 999

###############################################################################
#	PARTICLES & PARTICLE MASSES
###############################################################################

#only for hydros
#number of particles
set np_disk			= 999
set np_gas			= 999
set np_stars			= 999

set np_disk_1comma5ropt		= 999
set np_gas_1comma5ropt		= 999
set np_stars_1comma5ropt	= 999

#EAGLE total partical masses for different components
set mPart_DM 	= 999
set mPart_stars = 999
set mPart_gas 	= 999
set mPart_bh 	= 999

#EAGLE starforming gas
set mPart_gasSF		= 999
#EAGLE non-starforming gasmbh
set mPart_gasNSF		= 999

###############################################################################
#	FLAGS
###############################################################################

set orphan		= 999
set isolated		= 999

set flag_SB 				= 999
set flag_SB_envr			= 999
set flag_SB_B 				= 999
set flag_SB_B_new			= 29
set flag_SB_r 				= 999
set flag_SFE 				= 999
set flag_SBplusSFE 			= 999
set flag_sample 			= 999
set flag_mstar_bin 			= 999
set flag_mhalo_bin 			= 999
set flag_sfe_bin 			= 999
set flag_age_bin			= 999
set flag_SHMR_bin			= 999

set flag_t50_st2ha			= 999
set flag_t50_bh2ha			= 999
set flag_t50_bh2st			= 999

set flag_t70_st2ha			= 999
set flag_t70_bh2ha			= 999
set flag_t70_bh2st			= 999

set flag_tform_halo			= 999
set flag_tform_stars			= 999
set flag_tform_bh			= 999

set flag_mhalo				= 999
set flag_majorM				= 30

set flag_dens				= 999
set flag_dens3				= 999
set flag_dens6				= 999
set flag_dens3b				= 999
set flag_dens6b				= 999
set flag_dens3b_new			= 999
set flag_dens6b_new			= 2

set flag_std_fit_t70_stars		= 999
set flag_std_fit_t70_halo		= 999
set flag_std_fit_t70_bh			= 999

set flag_std_fit_t50_stars		= 999
set flag_std_fit_t50_halo		= 999
set flag_std_fit_t50_bh			= 999

set flag_use_inter_stars		= 999
set flag_use_inter_halo			= 999
set flag_use_inter_bh			= 999

set flag_use_fit_stars			= 999
set flag_use_fit_halo			= 999
set flag_use_fit_bh			= 999

###############################################################################
#	ENVIRONMENT
###############################################################################

set envr			= 6

#Environment as in VWeb (Cui et al. 2019)
set env_512			= 999
set env_1024		= 999

#Densities EAGLE for Silvio
set sDensity_N10	= 32
set sDensity_N7		= 999
set sDensity_N5		= 34


###############################################################################
#	RADII & DIAMETER
###############################################################################

#only SAG
set rhalf_bulge		= 999
set rhalf_disk		= 999

#only SAGE and Galacticus
set rdisk		= 999

#only Galacticus
set rbulge			= 999
set rhalf_mass		= 999


set ropt			= 35

set reff_gas					= 999
set reff_gas_disk 				= 999
set reff_gas_1comma5ropt		= 36
set reff_gas_disk_1comma5ropt 	= 88

set r200c			= 37
set rhalf_DM		= 999
set rhalf_DM_2D		= 999
set rhalf_stars		= 97
set rhalf_stars_2D	= 999
set rhalf_stars_30kpc		= 38
set rhalf_stars_30kpc_2D	= 999
set rhalf_gas		= 999
set rhalf_gas_2D	= 999
set rhalf_bh		= 999
set rhalf_bh_2D		= 999
set rhalf_stars_1comma5ropt		= 39
set rhalf_gas_1comma5ropt		= 999

set Dgas_disk_1comma5ropt		= 999

set rhalf_sfr				= 999
set rVmax					= 40
set rcenter2r200c 			= 999

set rVoid					= 42
set DtoVoidCent 			= 43

###############################################################################
#	AGES
###############################################################################

#only SAG and SAGE
set mean_age_stars	= 999

#only Galacticus
set mean_age_stars_disk	= 999
set age_sfr_int_disk 	= 999

set mean_age_stars_spheroid	= 999
set age_sfr_int_spheroid = 999

set age_stars_rband_r502D	= 44
set age_stars_rband_2r502D	= 45

set delta_age_stars_rband	= 91

set z_mean_birth_stars	= 999
set age_mean_stars		= 999

set t50_stars				= 46
set t70_stars				= 999

set t50_halo				= 47
set t70_halo				= 999

set t50_bh					= 48
set t70_bh					= 999

set t50_curve_fit_stars				= 999
set t70_curve_fit_stars				= 999

set t50_curve_fit_halo				= 999
set t70_curve_fit_halo				= 999

set t50_curve_fit_bh					= 999
set t70_curve_fit_bh					= 999

set t50_cs_fit_stars				= 999
set t70_cs_fit_stars				= 999

set t50_cs_fit_halo				= 999
set t70_cs_fit_halo				= 999

set t50_cs_fit_bh					= 999
set t70_cs_fit_bh					= 999

set t50_stars_std				= 999
set t70_stars_std				= 999

set t50_halo_std				= 999
set t70_halo_std				= 999

set t50_bh_std					= 999
set t70_bh_std					= 999

set tmin_rVmax					= 999
set tmax_rVmax					= 999

set tmin_angM_stars					= 999
set tmax_angM_stars					= 999

set tmin_angM_SFgas					= 999
set tmax_angM_SFgas					= 999

set tmin_angM_NSFgas					= 999
set tmax_angM_NSFgas					= 999

set tmin_angM_norm_stars					= 999
set tmax_angM_norm_stars					= 999

set tmin_angM_norm_SFgas					= 999
set tmax_angM_norm_SFgas					= 999

set tmin_angM_norm_NSFgas					= 999
set tmax_angM_norm_NSFgas					= 999

set tmin_sfr					= 999
set tmax_sfr					= 999

set tmin_ssfr					= 999
set tmax_ssfr					= 999

set tmin_bh_acc_rate					= 999
set tmax_bh_acc_rate					= 999

set tmin_vmax					= 999
set tmax_vmax					= 999

set tmin_vdisp					= 999
set tmax_vdisp					= 999

set tmin_mgas					= 999
set tmax_mgas					= 999

set tmin_SHMR					= 999
set tmax_SHMR					= 999

set tmin_mgas_SF					= 999
set tmax_mgas_SF					= 999

set tmin_mgas_NSF					= 999
set tmax_mgas_NSF					= 999

set tmin_vel					= 999
set tmax_vel					= 999

set tmin_age					= 999
set tmax_age					= 999

###############################################################################
#	DELTAS
###############################################################################

set delta_tform_stars	= 999
set delta_tform_halo	= 999
set delta_tform_bh		= 999

set delta_tform_curve_fit_stars		= 999
set delta_tform_curve_fit_halo		= 999
set delta_tform_curve_fit_bh		= 999

set delta_tform_cs_fit_stars	= 999
set delta_tform_cs_fit_halo		= 999
set delta_tform_cs_fit_bh		= 999

set delta_tvar_bh_acc_rate		= 999
set delta_tvar_rVmax			= 999
set delta_tvar_angM_stars		= 999
set delta_tvar_angM_SFgas		= 999
set delta_tvar_angM_NSFgas		= 999
set delta_tvar_angM_norm_stars		= 999
set delta_tvar_angM_norm_SFgas		= 999
set delta_tvar_angM_norm_NSFgas		= 999
set delta_tvar_sfr	= 999
set delta_tvar_ssfr	= 999
set delta_tvar_vmax	= 999
set delta_tvar_vdisp	= 999
set delta_tvar_mgas_SF	= 999
set delta_tvar_mgas_NSF	= 999
set delta_tvar_mgas	= 999
set delta_tvar_SHMR	= 999
set delta_tvar_vel	= 999
set delta_tvar_age	= 999

set delta_t50				= 999

set delta_t50_st2ha			= 999
set delta_t50_bh2ha			= 999
set delta_t50_bh2st			= 999

set delta_t70_st2ha			= 999
set delta_t70_bh2ha			= 999
set delta_t70_bh2st			= 999

set delta_inter_fit_t50_stars		= 999
set delta_inter_fit_t50_halo		= 999
set delta_inter_fit_t50_bh			= 999

set delta_inter_fit_t70_stars		= 999
set delta_inter_fit_t70_halo		= 999
set delta_inter_fit_t70_bh			= 999

set nr_strikes_stars				= 999
set nr_strikes_halo					= 999
set nr_strikes_bh					= 999

###############################################################################
#	ENERGIES
###############################################################################

set Etot		= 999
set Ekin		= 999
set Emech		= 999
set Etherm_gas	= 999

set Etot_gasSF		= 999
set Ekin_gasSF		= 999
set Etherm_gasSF 	= 999

set Etot_gasNSF		= 999
set Ekin_gasNSF		= 999
set Etherm_gasNSF 	= 999

###############################################################################
#	ANGULAR MOMENTA
###############################################################################

set angM_disk		= 999
set angM_spheroid	= 999

set angM_stars		= 999
set angM_SFgas		= 999
set angM_NSFgas		= 999

#normalised
set angM_norm_1comma5ropt_stars	= 49
set angM_norm_1comma5ropt_gas	= 50
set angM_norm_1comma5ropt_bar	= 999

set angM_norm_stars		= 999
set angM_norm_SFgas		= 999
set angM_norm_NSFgas	= 999


###############################################################################
#	SPIN 
###############################################################################

#EAGLE total gas in particles
set spinGas_x = 999
set spinGas_y = 999
set spinGas_z = 999

#EAGLE starforming gas
set spinGasSF_x = 999
set spinGasSF_y = 999
set spinGasSF_z = 999

#EAGLE non-starforming gas
set spinGasNSF_x = 999
set spinGasNSF_y = 999
set spinGasNSF_z = 999

#EAGLE total stars in particles
set spinStars_x = 999
set spinStars_y = 999
set spinStars_z = 999

#only Galacticus and SAGE		
set spinParameter	= 999

#EAGLE stellar
set lambda_r502D_edgeOn	= 54
set lambda_r502D_random	= 999

set lambda_2r502D_edgeOn	= 55
set lambda_2r502D_random	= 999

set delta_lambda_edgeOn		= 56
set delta_lambda_random		= 999

###############################################################################
#	BLACK HOLES
###############################################################################

set mbh			= 57
set bhcount		= 999
set bheff		= 999
set bh_acc_rate		= 999

###############################################################################
#	STELLAR
###############################################################################

set mstar_spheroid	= 999
set mstar_disk		= 999
set mstar		= 999

set mstar_30kpc 	= 58
set mstar_50kpc		= 999
set mstar_70kpc		= 999

set mstar_1comma5ropt 	= 59
set mstar_half_reff_gas_disk 	= 99
set mstar_birth 	= 99

set BvT				= 999
set BvT_Lagos17b	= 61
set DvT_stars		= 999
set DvT_stars_alt	= 999
set DvT_gas			= 999
set DvT_30kpc_T19	= 96

set DvT_stars_c0comma4		= 999
set DvT_stars_c0comma5		= 999
set DvT_gas_c0comma4		= 999
set DvT_gas_c0comma5		= 999

set DvT_stars_c0comma4_1comma5ropt		= 999
set DvT_stars_c0comma5_1comma5ropt		= 62
set DvT_gas_c0comma4_1comma5ropt		= 999
set DvT_gas_c0comma5_1comma5ropt		= 63

set mstar_IC		= 999
set mstarPlusIC		= 999

set SHMR 			 = 999
set SHMR_30kpc 		 = 64
set SHMR_1comma5ropt = 65

###############################################################################
#	SURFACE DENSITIES & BRIGHTNESSES
###############################################################################

#Surface Densities
#----------------------------------------

set Sigma_gas_1comma5ropt			= 51
set Sigma_gas_reff_1comma5ropt		= 102
set Sigma_gas_reff_disk_1comma5ropt	= 60
set Sigma_stars_1comma5ropt		= 52
set Sigma_stars_30kpc			= 101
set Sigma_sfr_1comma5ropt		= 53
set Sigma_sfr_30kpc				= 98

set Sigma_HI_30kpc				= 99
set Sigma_H2_30kpc				= 100
set Sigma_HIH2_30kpc			= 83

#star formation efficiency
#----------------------------------------

set sfe_gas_reff_disk_1comma5ropt	= 106
set sfe_gas_reff_1comma5ropt		= 107
set sfe_gas_1comma5ropt				= 103
set sfe_HI_30kpc					= 104
set sfe_HIH2_30kpc					= 105

#Surface Brightnesses
#----------------------------------------

set SB_mu_eff_gas_disk_B			= 999
set SB_mu_ropt_B					= 999
set SB_mu_1comma5ropt_B				= 66
set SB_mu_corr_ropt_B				= 67
set SB_mu_eff_stars_1comma5ropt_B	= 999
set SB_mu_eff_stars_30kpc_B			= 999

set SB_mu_eff_gas_disk_r			= 999
set SB_mu_opt_r						= 999
set SB_mu_eff_stars_1comma5ropt_r	= 999
set SB_mu_eff_stars_30kpc_r			= 999

set SB_mu_eff_gas_disk_V			= 999
set SB_mu_opt_V						= 999
set SB_mu_eff_stars_1comma5ropt_V	= 999
set SB_mu_eff_stars_30kpc_V			= 999

set SB_mu_eff_gas_disk_g			= 999
set SB_mu_opt_g						= 999
set SB_mu_eff_stars_1comma5ropt_g	= 999
set SB_mu_eff_stars_30kpc_g			= 999

###############################################################################
#	MASSES OF GAS
###############################################################################

#Masses of gas
#----------------------------------------
set Mgas_disk		= 999

#only SAG and Galacticus
set Mgas_spheroid	= 999
set Mgas		= 999

#only hydros
set Mgas_HI		= 999
set Mgas_H2		= 999

set Mgas_HI_60kpc_GK11		= 999
set Mgas_H2_60kpc_GK11		= 999

set Mgas_HI_60kpc_K13		= 999
set Mgas_H2_60kpc_K13		= 999

set Mgas_HI_30kpc_GK11		= 93
set Mgas_H2_30kpc_GK11		= 94
set Mgas_HIH2_30kpc_GK11	= 95

set Mgas_HI_30kpc_K13		= 999
set Mgas_H2_30kpc_K13		= 999

set mcold_disk		= 999

#only Galacitcus and SAG
set mcold_spheroid	= 999
set mcold		= 999

set mhot		= 999

set mgas_30kpc		= 70
set mgas_50kpc		= 999
set mgas_70kpc		= 999

set mgas_1comma5ropt 	= 71
set mgas_half_reff_gas_disk	= 999

set mgas_SF_30kpc		= 999
set mgas_SF		= 999
set mgas_NSF		= 999


###############################################################################
#	STAR FORMATION
###############################################################################

set sfr_spheroid	= 999
set sfr_disk		= 999
set sfr			= 999
set ssfr		= 999

set sfr_total	= 87
set sfr_30kpc	= 72
set ssfr_30kpc 	= 73
set sfr_50kpc	= 999
set sfr_70kpc	= 999
set sfr_1comma5ropt	= 74
set ssfr_1comma5ropt = 75

set sfr_int_disk	= 999
set sfr_int_spheroid	= 999

set sfr_int_quies	= 999
set sfr_gasSF		= 999
set ssfr_gasSF		= 999


###############################################################################
#	METALLICITIES & MASSES OF METALS
###############################################################################

#only SAG and SAGE
set Mzstar_spheroid	= 999
set Mzstar_disk		= 999
set Mzstar			= 999
set Mzhot_halo		= 999

#only SAG
set OH_gas_disk			= 999
set OH_gas_bulge		= 999
set OH_gas_disk_bulge	= 999

#metallicities / abundances
#----------------------------------------
#only Galacticus
set zgas_spheroid	= 999

#only Galacticus and SAGE
set zgas_disk		= 999
set zcold			= 999

#EAGLE starforming gas
set zgasSF_H		= 999
set zgasSF_O		= 999
set zgasSF_N		= 999

set OH_gas_30kpc	= 86

set metalfrac_SF 	= 999

set zcold_gasSF		= 999
set z_grad			= 999

#EAGLE non-starforming gas
set zgasNSF_H		= 999
set zgasNSF_O		= 999
set zgasNSF_N		= 999
set metalfrac_NSF 	= 999

#only Galacticus			 
set zstar_spheroid	= 999
set zstar_disk		= 999
set zhot_halo		= 999

###############################################################################
#	OBSERVATIONS
###############################################################################

set DEC			= 999
set RA			= 999
set Z			= 999
set age			= 999
set weight_tot	= 999

###############################################################################
#	OTHER PROPERTIES
###############################################################################

#Gas depletion time mcold/sfr
set Tcons					= 999
set Tcons_30kpc				= 22
set Tcons_1comma5ropt		= 999
set Tcons_HI_30kpc			= 999
set Tcons_HIH2_30kpc		= 999
set Tcons_gasSF				= 999

#cold gas fraction mcold/mstar
set cgf					= 999
set cgf_30kpc			= 999
set cgf_HI_30kpc		= 999
set cgf_HIH2_30kpc		= 23
set cgf_1comma5ropt		= 999
set cgf_gasSF			= 999

#total baryon fraction: mcold/(mstar+mcold)
set fbar				= 999
set fbar_30kpc			= 68
set fbar_1comma5ropt	= 999
set fbar_HI_30kpc		= 999
set fbar_HIH2_30kpc		= 999
set fbar_gasSF			= 999

#atomic gas fraction according to Rosas-Guevara+22 (page 5) and Obreschkow+16
set fatom_30kpc_GK11	= 31
set fmol_30kpc_GK11		= 92

###############################################################################
#	LAST TIME CENTRAL (ltc) PROPERTIES
###############################################################################

set lastTimeCentral = 999

set mstar_ltc	= 999
set sfr_ltc		= 999
set ssfr_ltc	= 999

set VvS_r502D_edgeOn_ltc	= 999
set VvS_r502D_random_ltc	= 999
set VvS_2r502D_edgeOn_ltc	= 999
set VvS_2r502D_random_ltc	= 999

set lambda_r502D_random_ltc	= 999
set lambda_2r502D_random_ltc	= 999

set angM_r502D_random_ltc	= 999
set angM_2r502D_random_ltc	= 999

###############################################################################
#	LUMINOSITIES, FLUXES & MAGNITUDES
###############################################################################

#Luminosities
#----------------------------------------

set L_bolo				= 999

set L_SDSS_spheroid_u	= 999
set L_SDSS_spheroid_g	= 999
set L_SDSS_spheroid_r	= 999
set L_SDSS_spheroid_i	= 999
set L_SDSS_spheroid_z	= 999

set L_SDSS_dA_spheroid_u	= 999
set L_SDSS_dA_spheroid_g	= 999
set L_SDSS_dA_spheroid_r	= 999
set L_SDSS_dA_spheroid_i	= 999
set L_SDSS_dA_spheroid_z	= 999

set L_SDSS_disk_u	= 999
set L_SDSS_disk_g	= 999
set L_SDSS_disk_r	= 999
set L_SDSS_disk_i	= 999
set L_SDSS_disk_z	= 999

set L_SDSS_dA_disk_u	= 999
set L_SDSS_dA_disk_g	= 999
set L_SDSS_dA_disk_r	= 999
set L_SDSS_dA_disk_i	= 999
set L_SDSS_dA_disk_z	= 999

set L_SDSS_u	= 999
set L_SDSS_g	= 999
set L_SDSS_r	= 999
set L_SDSS_i	= 999
set L_SDSS_z	= 999

set L_SDSS_dA_u	= 999
set L_SDSS_dA_g	= 999
set L_SDSS_dA_r	= 999
set L_SDSS_dA_i	= 999
set L_SDSS_dA_z	= 999

set L_SDSS_dA_total_u	= 999
set L_SDSS_dA_total_g	= 999
set L_SDSS_dA_total_r	= 999
set L_SDSS_dA_total_i	= 999
set L_SDSS_dA_total_z	= 999

#Fluxes Standard
#----------------------------------------

set F_dA_SDSS_u				= 999
set F_dA_SDSS_g				= 999
set F_dA_SDSS_r				= 999
set F_dA_SDSS_i				= 999
set F_dA_SDSS_z				= 999

set F_dA_Johnson_U				= 999
set F_dA_Johnson_B				= 999
set F_dA_Johnson_V				= 999
set F_dA_Johnson_R				= 999


#Magnitudes Standard
#----------------------------------------
set mAB_dA_SDSS_u	= 999
set mAB_dA_SDSS_g	= 999
set mAB_dA_SDSS_r	= 999
set mAB_dA_SDSS_i	= 999
set mAB_dA_SDSS_z	= 999

set mAB_dA_Johnson_U	= 999
set mAB_dA_Johnson_B	= 999
set mAB_dA_Johnson_V	= 999
set mAB_dA_Johnson_R	= 999

set mAB_SDSS_u	= 999
set mAB_SDSS_g	= 999
set mAB_SDSS_r	= 999
set mAB_SDSS_i	= 999
set mAB_SDSS_z	= 999

set mAB_Johnson_U	= 999
set mAB_Johnson_B	= 999
set mAB_Johnson_V	= 999
set mAB_Johnson_R	= 999

#------------------------------------------

set MAB_SDSS_u		= 999
set MAB_SDSS_g		= 999
set MAB_SDSS_r		= 999
set MAB_SDSS_i		= 999
set MAB_SDSS_z		= 999

set MAB_dA_SDSS_u	= 999
set MAB_dA_SDSS_g	= 76
set MAB_dA_SDSS_r	= 77
set MAB_dA_SDSS_i	= 78
set MAB_dA_SDSS_z	= 999

set MAB_Johnson_U	= 999
set MAB_Johnson_B	= 999
set MAB_Johnson_V	= 999
set MAB_Johnson_R	= 999

set MAB_dA_Johnson_U	= 999
set MAB_dA_Johnson_B	= 79
set MAB_dA_Johnson_V	= 80
set MAB_dA_Johnson_R	= 999
#----------------------------------------

#Further Magnitudes
#----------------------------------------
set mag_u		= 999
set mag_g		= 999
set mag_r		= 999
set mag_i		= 999
set mag_z		= 999

set mAB_u		= 999
set mAB_g		= 999
set mAB_r		= 999
set mAB_i		= 999
set mAB_z		= 999

set MAB_total_u		= 999
set MAB_total_g		= 999
set MAB_total_r		= 999
set MAB_total_i		= 999
set MAB_total_z		= 999

set mAB_total_u		= 999
set mAB_total_g		= 999
set mAB_total_r		= 999
set mAB_total_i		= 999
set mAB_total_z		= 999

set mAB_u	= 999
set mAB_g	= 999
set mAB_r	= 999
set mAB_i	= 999
set mAB_z	= 999

set mAB_dA_u		= 999
set mAB_dA_g		= 999
set mAB_dA_r		= 999
set mAB_dA_i		= 999
set mAB_dA_z		= 999

set mAs_u	= 999
set mAs_g	= 999
set mAs_r	= 999
set mAs_i	= 999
set mAs_z	= 999

set mAs_dA_u		= 999
set mAs_dA_g		= 999
set mAs_dA_r		= 999
set mAs_dA_i		= 999
set mAs_dA_z		= 999

set mAs_dA_total_u		= 999
set mAs_dA_total_g		= 999
set mAs_dA_total_r		= 999
set mAs_dA_total_i		= 999
set mAs_dA_total_z		= 999
#----------------------------------------

#Target selection cuts
#----------------------------------------
set mAB_total_cut_r_i	= 999
set mAB_total_cut_g_r	= 999
set mAB_total_cut_g_i	= 999
set mAB_total_cut_dmesa	= 999
set mAB_total_cut_i_lt_dmesa	= 999

set mAB_dA_total_cut_r_i	= 81
set mAB_dA_total_cut_g_r	= 999
set mAB_dA_total_cut_g_i	= 999
set mAB_dA_total_cut_u_g	= 999
set mAB_dA_total_cut_u_i	= 999
set mAB_dA_total_cut_i_z	= 999

set mAB_dA_total_cut_B_V	= 82

set mAB_dA_total_cut_dmesa	= 999
set mAB_dA_total_cut_i_lt_dmesa	= 999
set mAB_dA_total_cut_i_lt_dmesa_sparse	= 999

set mAB_dA_total_cut_cmesa	= 999
set mAB_dA_total_cut_cpar	= 999
set mAB_dA_total_cut_r_lt_cpar	= 999

set mAs_dA_total_cut_r_i	= 999
set mAs_dA_total_cut_dmesa	= 999
set mAs_dA_total_cut_g_r	= 999
set mAs_dA_total_cut_i_lt_dmesa	= 999
#----------------------------------------

#Emission lines
#----------------------------------------
set OII_3727_ext	= 999
set OII_3727		= 999
set OII_3729_ext	= 999
set OII_3729		= 999
#----------------------------------------

#Continuum Emission lines
#----------------------------------------
set OII_cont_3727_ext	= 999
set OII_cont_3727	= 999
set OII_cont_3729_ext	= 999
set OII_cont_3729	= 999


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
echo 'plot_custom_loop= '$PLOT_CUSTOM_LOOP	>> $MY_PHYSICS_SPECS

echo 'MY_PHYSICS_SPECS' $MY_PHYSICS_SPECS	 		>> $MY_PATH_HANDLER_FILE

set method_name_array = (plotXY plotOnly filterData analyseTargetSelection SMF SFRF sSFRF cSFRD sfr2z sfr2mstar ssfr2mstar oh2mstar HMF mstar2mhalo mstar2mhalofunc ngal2mhalo HMF_no mstar2mhalo_no mstar2mhalofunc_no ngal2mhalo_no zgas2mstar zgas2mcold mcold2mstar mbh2mstarsph twoPCF HOD mstar2mhalovsSFR mstar2rhalf)

set method_array = ($plotXY $plotOnly $filterData $analyseTargetSelection $SMF $SFRF $sSFRF $cSFRD $sfr2z $sfr2mstar $ssfr2mstar $oh2mstar $HMF $mstar2mhalo $mstar2mhalofunc $ngal2mhalo $HMF_no $mstar2mhalo_no $mstar2mhalofunc_no $ngal2mhalo_no $zgas2mstar $zgas2mcold $mcold2mstar $mbh2mstarsph $twoPCF $HOD $mstar2mhalovsSFR $mstar2rhalf)

#echo 559 method_name_array $method_array

#########################################################################
#	GENERATE FILTER AND CUT VALUE FILE 	"CUT_CONDS"		#
#########################################################################

set halo_col_name_array = (halo_nhalos halo_haloid halo_hostid halo_desIndex halo_nodemass halo_vmax halo_rsca halo_angM halo_spin halo_redshift halo_x_pos halo_y_pos halo_z_pos halo_x_vel halo_y_vel halo_z_vel halo_x_vel_disp halo_y_vel_disp halo_z_vel_disp)

set halo_col_array = ($halo_nhalos $halo_haloid $halo_hostid $halo_desIndex $halo_nodemass $halo_vmax $halo_rsca $halo_angM $halo_spin $halo_redshift $halo_x_pos $halo_y_pos $halo_z_pos $halo_x_vel $halo_y_vel $halo_z_vel $halo_x_vel_disp $halo_y_vel_disp $halo_z_vel_disp)

set col_name_array = (ngalaxies haloid hostid satelliteNodeIndex parentIndex orphan mhalo vmax vpeak spinParameter NFW_con zgas_spheroid zgas_disk zstar_disk zstar_spheroid zhot_halo mcold_spheroid mcold_disk mcold mbh mstar_spheroid mstar_disk mstar mhot Mgas_spheroid Mgas_disk Mgas Mzstar_spheroid Mzstar_disk Mzstar Mzhot_halo sfr_spheroid sfr_disk sfr x_pos y_pos z_pos x_vel y_vel z_vel L_SDSS_spheroid_u L_SDSS_spheroid_g L_SDSS_spheroid_r L_SDSS_spheroid_i L_SDSS_spheroid_z L_SDSS_dA_spheroid_u L_SDSS_dA_spheroid_g L_SDSS_dA_spheroid_r L_SDSS_dA_spheroid_i L_SDSS_dA_spheroid_z L_SDSS_disk_u L_SDSS_disk_g L_SDSS_disk_r L_SDSS_disk_i L_SDSS_disk_z L_SDSS_dA_disk_u L_SDSS_dA_disk_g L_SDSS_dA_disk_r L_SDSS_dA_disk_i L_SDSS_dA_disk_z L_SDSS_u L_SDSS_g L_SDSS_r L_SDSS_i L_SDSS_z L_SDSS_dA_u L_SDSS_dA_g L_SDSS_dA_r L_SDSS_dA_i L_SDSS_dA_z L_SDSS_dA_total_u L_SDSS_dA_total_g L_SDSS_dA_total_r L_SDSS_dA_total_i L_SDSS_dA_total_z mag_u mag_g mag_r mag_i mag_z mAB_u mAB_g mAB_r mAB_i mAB_z mAB_total_u mAB_total_g mAB_total_r mAB_total_i mAB_total_z mAB_dA_u mAB_dA_g mAB_dA_r mAB_dA_i mAB_dA_z mAB_dA_SDSS_u mAB_dA_SDSS_g mAB_dA_SDSS_r mAB_dA_SDSS_i mAB_dA_SDSS_z MAB_dA_SDSS_u MAB_dA_SDSS_g MAB_dA_SDSS_r MAB_dA_SDSS_i MAB_dA_SDSS_z MAB_total_u MAB_total_g MAB_total_r MAB_total_i MAB_total_z mAs_u mAs_g mAs_r mAs_i mAs_z mAs_dA_u mAs_dA_g mAs_dA_r mAs_dA_i mAs_dA_z mAs_dA_total_u mAs_dA_total_g mAs_dA_total_r mAs_dA_total_i mAs_dA_total_z mAB_total_cut_r_i mAB_total_cut_dmesa mAB_total_cut_g_r mAB_total_cut_i_lt_dmesa 'mAB_dA_total_cut_r_i' 'mAB_dA_total_cut_g_r' 'mAB_dA_total_cut_g_i' 'mAB_dA_total_cut_u_g' 'mAB_dA_total_cut_u_i' 'mAB_dA_total_cut_i_z' mAB_dA_total_cut_dmesa mAB_dA_total_cut_cmesa mAB_dA_total_cut_cpar mAB_dA_total_cut_i_lt_dmesa mAB_dA_total_cut_i_lt_dmesa_sparse mAB_dA_total_cut_r_lt_cpar mAs_dA_total_cut_r_i mAs_dA_total_cut_dmesa mAs_dA_total_cut_g_r mAs_dA_total_cut_i_lt_dmesa RA DEC Z age mean_age_stars mhalo_sat OH_gas_disk OH_gas_bulge OH_gas_disk_bulge weight_tot OII_3727_ext OII_3727 OII_3729_ext OII_3729 nodeIndex satelliteIndex siblingIndex satelliteMergeTime isolated timeLastIsolated firstProgenitorID npros rhalf_bulge rhalf_disk rbulge rdisk rhalf_mass mstar_IC OII_cont_3727_ext OII_cont_3727 OII_cont_3729_ext OII_cont_3729 sfr_int_spheroid sfr_quies_int mhalo_200c mbasic mbasic_200c mhalo_cents_200c mhalo_cents 'mstar+IC' ssfr env_512 env_1024 zcold lastProgID mainLeafID fofID bhcount vbulge vdisk age_sfr_int_disk sfr_int_disk sfr_int_spheroid age_sfr_int_spheroid sfr_int_spheroid mean_age_stars_disk mean_age_stars_spheroid cgf angM_spheroid angM_disk Tcons fbar bheff BvT Mgas_HI Mgas_H2 envr lastMinorM lastMajorM DvT_stars DvT_stars_alt DvT_gas np_stars np_gas np_disk ropt rhalf_stars r200c sfr_gasSF bh_acc_rate mAB_total_cut_g_i SHMR rhalf_gas rhalf_stars_2D vdisp_r502D_edgeOn vdisp_r502D_random vdisp_2r502D_random VvS_r502D_edgeOn VvS_r502D_random ell_r502D_edgeOn ell_r502D_random Sersic_n Mgas_HI_60kpc_GK11 Mgas_HI_60kpc_K13 Mgas_H2_60kpc_GK11 Mgas_H2_60kpc_K13 vdisp_2r502D_edgeOn VvS_2r502D_edgeOn VvS_2r502D_random ell_2r502D_edgeOn ell_2r502D_random age_stars_rband_r502D age_stars_rband_2r502D angM_r502D_random_ltc angM_2r502D_random_ltc mstar_ltc sfr_ltc ssfr_ltc VvS_r502D_random_ltc VvS_2r502D_random_ltc VvS_r502D_edgeOn_ltc VvS_2r502D_edgeOn_ltc lastTimeCentral rcenter2r200c sDensity_N10 sDensity_N7 sDensity_N5 lastMerger ratio_lastM rhalf_sfr lambda_r502D_edgeOn lambda_r502D_random lambda_2r502D_edgeOn lambda_2r502D_random lambda_r502D_random_ltc lambda_2r502D_random_ltc nsats Mgas_HI_30kpc_GK11 Mgas_H2_30kpc_GK11 Mgas_HI_30kpc_K13 Mgas_H2_30kpc_K13 mstar_30kpc angM_norm_1.5ropt_stars angM_norm_1.5ropt_gas reff_gas_disk DvT_stars_c0.4 DvT_stars_c0.5 DvT_gas_c0.4 DvT_gas_c0.5 np_disk_1.5ropt np_stars_1.5ropt np_gas_1.5ropt jsub mstar_1.5ropt rhalf_stars_1.5ropt rhalf_gas_1.5ropt mhalo_fof galaxyID nSubhalos x_pos_subhalo y_pos_subhalo z_pos_subhalo mPart_DM mPart_stars mPart_gas mPart_bh rhalf_DM spinGas_x spinGas_y spinGas_z rhalf_gas_2D rhalf_bh rhalf_bh_2D Etot Ekin Emech Etherm_gas mstar_birth rVmax vdisp rhalf_DM_2D sfr_1.5ropt ssfr_1.5ropt sfr_30kpc BvT_Lagos17b DvT_30kpc_T19 ell_30kpc_T19 vdisp_30kpc_T19 Vrot2Vdisp_30kpc_T19 ellDM_30kpc_T19 sfr_total zgasSF_H zgasSF_O zgasSF_N zgasNSF_H zgasNSF_O zgasNSF_N reff_gas mgas_30kpc mgas_50kpc mgas_70kpc mstar_50kpc mstar_70kpc sfr_50kpc sfr_70kpc vdisp_50kpc vdisp_70kpc mhalo_50kpc mhalo_70kpc spinGasSF_x spinGasSF_y spinGasSF_z spinGasNSF_x spinGasNSF_y spinGasNSF_z mPart_gasSF mPart_gasNSF metalfrac_SF metalfrac_NSF Etot_gasSF Ekin_gasSF Etherm_gasSF Etot_gasNSF Ekin_gasNSF Etherm_gasNSF reff_gas_1.5ropt reff_gas_disk_1.5ropt DvT_stars_c0.4_1.5ropt DvT_stars_c0.5_1.5ropt DvT_gas_c0.4_1.5ropt DvT_gas_c0.5_1.5ropt mgas_1.5ropt mstar_half_reff_gas_disk mgas_half_reff_gas_disk Sigma_HI_30kpc Sigma_H2_30kpc Sigma_HIH2_30kpc sfe_HIH2_30kpc SB_mu_eff_stars_1.5ropt_r SB_mu_eff_stars_30kpc_r SB_mu_eff_gas_disk_B SB_mu_ropt_B ssfr_30kpc ssfr_gasSF zcold_gasSF Sigma_gas_1.5ropt Sigma_gas_reff_disk_1.5ropt Sigma_stars_1.5ropt Sigma_stars_30kpc Sigma_sfr_1.5ropt Sigma_sfr_30kpc Tcons_30kpc Tcons_1.5ropt Tcons_HIH2_30kpc Tcons_gasSF cgf_30kpc cgf_1.5ropt cgf_HIH2_30kpc cgf_gasSF fbar_30kpc fbar_1.5ropt fbar_HIH2_30kpc fbar_gasSF SHMR_1.5ropt rhalf_stars_30kpc rhalf_stars_30kpc_2D F_dA_Johnson_U F_dA_Johnson_B F_dA_Johnson_V F_dA_Johnson_R F_dA_SDSS_u F_dA_SDSS_g F_dA_SDSS_r F_dA_SDSS_i F_dA_SDSS_z MAB_Johnson_U MAB_Johnson_B MAB_Johnson_V MAB_Johnson_R MAB_dA_Johnson_U MAB_dA_Johnson_B MAB_dA_Johnson_V MAB_dA_Johnson_R mAB_Johnson_U mAB_Johnson_B mAB_Johnson_V mAB_Johnson_R mAB_dA_Johnson_U mAB_dA_Johnson_B mAB_dA_Johnson_V mAB_dA_Johnson_R mAB_SDSS_u mAB_SDSS_g mAB_SDSS_r mAB_SDSS_i mAB_SDSS_z MAB_SDSS_u MAB_SDSS_g MAB_SDSS_r MAB_SDSS_i MAB_SDSS_z z_mean_birth_stars age_mean_stars sfe_HI_30kpc Tcons_HI_30kpc cgf_HI_30kpc fbar_HI_30kpc SHMR_30kpc SB_mu_eff_stars_1.5ropt_B SB_mu_eff_stars_30kpc_B SB_mu_eff_gas_disk_r SB_mu_opt_r SB_mu_eff_stars_1.5ropt_V SB_mu_eff_stars_30kpc_V SB_mu_eff_gas_disk_V SB_mu_opt_V SB_mu_eff_stars_1.5ropt_g SB_mu_eff_stars_30kpc_g SB_mu_eff_gas_disk_g SB_mu_opt_g Sigma_gas_reff_1.5ropt SB_mu_1.5ropt_B Mgas_HIH2_30kpc_GK11 sfe_gas_1.5ropt mAB_dA_total_cut_B_V sfe_gas_reff_1.5ropt sfe_gas_reff_disk_1.5ropt flag_SB flag_SB_B flag_SB_r flag_SFE flag_sample flag_SB+SFE flag_SB_envr DtoVoidCent delta_age_stars_rband t50_stars t70_stars OH_gas_30kpc z_grad L_bolo t50_halo t70_halo v200c delta_tform_stars delta_tform_halo flag_mstar_bin flag_mhalo_bin flag_sfe_bin flag_age_bin delta_t50 delta_lambda_edgeOn delta_lambda_random delta_tform_bh t50_bh t70_bh delta_t50_st2ha delta_t50_bh2ha flag_t50_st2ha flag_t50_bh2ha fatom_30kpc_GK11 tmin_rVmax tmax_rVmax delta_tvar_rVmax tmin_angM_stars tmax_angM_stars delta_tvar_angM_stars tmin_angM_SFgas tmax_angM_SFgas delta_tvar_angM_SFgas tmin_angM_NSFgas tmax_angM_NSFgas delta_tvar_angM_NSFgas angM_stars angM_SFgas angM_NSFgas fmol_30kpc_GK11 flag_SHMR_bin tmin_sfr tmax_sfr tmin_ssfr tmax_ssfr tmin_bh_acc_rate tmax_bh_acc_rate tmin_vmax tmax_vmax tmin_vdisp tmax_vdisp tmin_mgas_SF tmax_mgas_SF tmin_mgas_NSF tmax_mgas_NSF delta_tvar_sfr delta_tvar_ssfr delta_tvar_vmax delta_tvar_vdisp delta_tvar_mgas_SF delta_tvar_mgas_NSF tmin_SHMR tmax_SHMR tmin_mgas tmax_mgas delta_tvar_SHMR delta_tvar_mgas mhalo_30kpc spinStars_x spinStars_y spinStars_z mgas_SF mgas_NS vpec_norm x_pos_cof y_pos_cof z_pos_cof redshift topLeafID progFofID delta_tvar_bh_acc_rate vdisp_30kpc tmin_vel tmax_vel delta_tvar_vel tmin_age tmax_age delta_tvar_age angM_norm_stars angM_norm_SFgas angM_norm_NSFgas delta_tvar_angM_norm_stars delta_tvar_angM_norm_SFgas delta_tvar_angM_norm_NSFgas tmin_angM_norm_stars tmax_angM_norm_stars tmin_angM_norm_SFgas tmax_angM_norm_SFgas tmin_angM_norm_NSFgas tmax_angM_norm_NSFgas flag_majorM rVoid t50_curve_fit_stars t70_curve_fit_stars t50_cs_fit_stars t70_cs_fit_stars t50_stars_std t70_stars_std delta_tform_curve_fit_stars delta_tform_cs_fit_stars t50_curve_fit_halo t70_curve_fit_halo t50_cs_fit_halo t70_cs_fit_halo t50_halo_std t70_halo_std delta_tform_curve_fit_halo delta_tform_cs_fit_halo t50_curve_fit_bh t70_curve_fit_bh t50_cs_fit_bh t70_cs_fit_bh t50_bh_std t70_bh_std delta_tform_curve_fit_bh delta_tform_cs_fit_bh flag_std_fit_t50_stars flag_std_fit_t50_halo flag_std_fit_t50_bh delta_inter_fit_t50_stars delta_inter_fit_t50_halo delta_inter_fit_t50_bh flag_std_fit_t70_stars flag_std_fit_t70_halo flag_std_fit_t70_bh delta_inter_fit_t70_stars delta_inter_fit_t70_halo delta_inter_fit_t70_bh flag_use_inter_stars flag_use_inter_halo flag_use_inter_bh flag_use_fit_stars flag_use_fit_halo flag_use_fit_bh nr_strikes_stars nr_strikes_halo nr_strikes_bh delta_t50_bh2st delta_t70_st2ha delta_t70_bh2ha delta_t70_bh2st flag_t50_bh2st flag_t70_st2ha flag_t70_bh2ha flag_t70_bh2st flag_tform_stars flag_mhalo flag_tform_bh flag_tform_halo flag_dens flag_dens3 flag_dens6 flag_dens3b flag_dens6b angM_norm_1.5ropt_bar delta_VvS_edgeOn delta_vdisp_edgeOn flag_SB_B_new flag_dens3b_new flag_dens6b_new SB_mu_corr_ropt_B kappaCoRot_30kpc_T19 medOrbCirc_30kpc_T19 triax_30kpc_T19)

set col_array = ($ngalaxies $haloid $hostid $satelliteNodeIndex $parentIndex $orphan $mhalo $vmax $vpeak $spinParameter $NFW_con $zgas_spheroid $zgas_disk $zstar_disk $zstar_spheroid $zhot_halo $mcold_spheroid $mcold_disk $mcold $mbh $mstar_spheroid $mstar_disk $mstar $mhot $Mgas_spheroid $Mgas_disk $Mgas $Mzstar_spheroid $Mzstar_disk $Mzstar $Mzhot_halo $sfr_spheroid $sfr_disk $sfr $x_pos $y_pos $z_pos $x_vel $y_vel $z_vel $L_SDSS_spheroid_u $L_SDSS_spheroid_g $L_SDSS_spheroid_r $L_SDSS_spheroid_i $L_SDSS_spheroid_z $L_SDSS_dA_spheroid_u $L_SDSS_dA_spheroid_g $L_SDSS_dA_spheroid_r $L_SDSS_dA_spheroid_i $L_SDSS_dA_spheroid_z $L_SDSS_disk_u $L_SDSS_disk_g $L_SDSS_disk_r $L_SDSS_disk_i $L_SDSS_disk_z $L_SDSS_dA_disk_u $L_SDSS_dA_disk_g $L_SDSS_dA_disk_r $L_SDSS_dA_disk_i $L_SDSS_dA_disk_z $L_SDSS_u $L_SDSS_g $L_SDSS_r $L_SDSS_i $L_SDSS_z $L_SDSS_dA_u $L_SDSS_dA_g $L_SDSS_dA_r $L_SDSS_dA_i $L_SDSS_dA_z $L_SDSS_dA_total_u $L_SDSS_dA_total_g $L_SDSS_dA_total_r $L_SDSS_dA_total_i $L_SDSS_dA_total_z $mag_u $mag_g $mag_r $mag_i $mag_z $mAB_u $mAB_g $mAB_r $mAB_i $mAB_z $mAB_total_u $mAB_total_g $mAB_total_r $mAB_total_i $mAB_total_z $mAB_dA_u $mAB_dA_g $mAB_dA_r $mAB_dA_i $mAB_dA_z $mAB_dA_SDSS_u $mAB_dA_SDSS_g $mAB_dA_SDSS_r $mAB_dA_SDSS_i $mAB_dA_SDSS_z $MAB_dA_SDSS_u $MAB_dA_SDSS_g $MAB_dA_SDSS_r $MAB_dA_SDSS_i $MAB_dA_SDSS_z $MAB_total_u $MAB_total_g $MAB_total_r $MAB_total_i $MAB_total_z $mAs_u $mAs_g $mAs_r $mAs_i $mAs_z $mAs_dA_u $mAs_dA_g $mAs_dA_r $mAs_dA_i $mAs_dA_z $mAs_dA_total_u $mAs_dA_total_g $mAs_dA_total_r $mAs_dA_total_i $mAs_dA_total_z $mAB_total_cut_r_i $mAB_total_cut_dmesa $mAB_total_cut_g_r $mAB_total_cut_i_lt_dmesa $mAB_dA_total_cut_r_i $mAB_dA_total_cut_g_r $mAB_dA_total_cut_g_i $mAB_dA_total_cut_u_g $mAB_dA_total_cut_u_i $mAB_dA_total_cut_i_z $mAB_dA_total_cut_dmesa $mAB_dA_total_cut_cmesa $mAB_dA_total_cut_cpar $mAB_dA_total_cut_i_lt_dmesa $mAB_dA_total_cut_i_lt_dmesa_sparse $mAB_dA_total_cut_r_lt_cpar $mAs_dA_total_cut_r_i $mAs_dA_total_cut_dmesa $mAs_dA_total_cut_g_r $mAs_dA_total_cut_i_lt_dmesa $RA $DEC $Z $age $mean_age_stars $mhalo_sat $OH_gas_disk $OH_gas_bulge $OH_gas_disk_bulge $weight_tot $OII_3727_ext $OII_3727 $OII_3729_ext $OII_3729 $nodeIndex $satelliteIndex $siblingIndex $satelliteMergeTime $isolated $timeLastIsolated $firstProgenitorID $npros $rhalf_bulge $rhalf_disk $rbulge $rdisk $rhalf_mass $mstar_IC $OII_cont_3727_ext $OII_cont_3727 $OII_cont_3729_ext $OII_cont_3729 $sfr_int_spheroid $sfr_int_quies $mhalo_200c $mbasic $mbasic_200c $mhalo_cents_200c $mhalo_cents $mstarPlusIC $ssfr $env_512 $env_1024 $zcold $lastProgID $mainLeafID $fofID $bhcount $vbulge $vdisk $age_sfr_int_disk $sfr_int_disk $sfr_int_spheroid $age_sfr_int_spheroid $sfr_int_spheroid $mean_age_stars_disk $mean_age_stars_spheroid $cgf $angM_spheroid $angM_disk $Tcons $fbar $bheff $BvT $Mgas_HI $Mgas_H2 $envr $lastMinorM $lastMajorM $DvT_stars $DvT_stars_alt $DvT_gas $np_stars $np_gas $np_disk $ropt $rhalf_stars $r200c $sfr_gasSF $bh_acc_rate $mAB_total_cut_g_i $SHMR $rhalf_gas $rhalf_stars_2D $vdisp_r502D_edgeOn $vdisp_r502D_random $vdisp_2r502D_random $VvS_r502D_edgeOn $VvS_r502D_random $ell_r502D_edgeOn $ell_r502D_random $Sersic_n $Mgas_HI_60kpc_GK11 $Mgas_HI_60kpc_K13 $Mgas_H2_60kpc_GK11 $Mgas_H2_60kpc_K13 $vdisp_2r502D_edgeOn $VvS_2r502D_edgeOn $VvS_2r502D_random $ell_2r502D_edgeOn $ell_2r502D_random $age_stars_rband_r502D $age_stars_rband_2r502D $angM_r502D_random_ltc $angM_2r502D_random_ltc $mstar_ltc $sfr_ltc $ssfr_ltc $VvS_r502D_random_ltc $VvS_2r502D_random_ltc $VvS_r502D_edgeOn_ltc $VvS_2r502D_edgeOn_ltc $lastTimeCentral $rcenter2r200c $sDensity_N10 $sDensity_N7 $sDensity_N5 $lastMerger $ratio_lastM $rhalf_sfr $lambda_r502D_edgeOn $lambda_r502D_random $lambda_2r502D_edgeOn $lambda_2r502D_random $lambda_r502D_random_ltc $lambda_2r502D_random_ltc $nsats $Mgas_HI_30kpc_GK11 $Mgas_H2_30kpc_GK11 $Mgas_HI_30kpc_K13 $Mgas_H2_30kpc_K13 $mstar_30kpc $angM_norm_1comma5ropt_stars $angM_norm_1comma5ropt_gas $reff_gas_disk $DvT_stars_c0comma4 $DvT_stars_c0comma5 $DvT_gas_c0comma4 $DvT_gas_c0comma5 $np_disk_1comma5ropt $np_stars_1comma5ropt $np_gas_1comma5ropt $jsub $mstar_1comma5ropt $rhalf_stars_1comma5ropt $rhalf_gas_1comma5ropt $mhalo_fof $galaxyID $nSubhalos $x_pos_subhalo $y_pos_subhalo $z_pos_subhalo $mPart_DM $mPart_stars $mPart_gas $mPart_bh $rhalf_DM $spinGas_x $spinGas_y $spinGas_z $rhalf_gas_2D $rhalf_bh $rhalf_bh_2D $Etot $Ekin $Emech $Etherm_gas $mstar_birth $rVmax $vdisp $rhalf_DM_2D $sfr_1comma5ropt $ssfr_1comma5ropt $sfr_30kpc $BvT_Lagos17b $DvT_30kpc_T19 $ell_30kpc_T19 $vdisp_30kpc_T19 $Vrot2Vdisp_30kpc_T19 $ellDM_30kpc_T19 $sfr_total $zgasSF_H $zgasSF_O $zgasSF_N $zgasNSF_H $zgasNSF_O $zgasNSF_N $reff_gas $mgas_30kpc $mgas_50kpc $mgas_70kpc $mstar_50kpc $mstar_70kpc $sfr_50kpc $sfr_70kpc $vdisp_50kpc $vdisp_70kpc $mhalo_50kpc $mhalo_70kpc $spinGasSF_x $spinGasSF_y $spinGasSF_z $spinGasNSF_x $spinGasNSF_y $spinGasNSF_z $mPart_gasSF $mPart_gasNSF $metalfrac_SF $metalfrac_NSF $Etot_gasSF $Ekin_gasSF $Etherm_gasSF $Etot_gasNSF $Ekin_gasNSF $Etherm_gasNSF $reff_gas_1comma5ropt $reff_gas_disk_1comma5ropt $DvT_stars_c0comma4_1comma5ropt $DvT_stars_c0comma5_1comma5ropt $DvT_gas_c0comma4_1comma5ropt $DvT_gas_c0comma5_1comma5ropt $mgas_1comma5ropt $mstar_half_reff_gas_disk $mgas_half_reff_gas_disk $Sigma_HI_30kpc $Sigma_H2_30kpc $Sigma_HIH2_30kpc $sfe_HIH2_30kpc $SB_mu_eff_stars_1comma5ropt_r $SB_mu_eff_stars_30kpc_r $SB_mu_eff_gas_disk_B $SB_mu_ropt_B $ssfr_30kpc $ssfr_gasSF $zcold_gasSF $Sigma_gas_1comma5ropt $Sigma_gas_reff_disk_1comma5ropt $Sigma_stars_1comma5ropt $Sigma_stars_30kpc $Sigma_sfr_1comma5ropt $Sigma_sfr_30kpc $Tcons_30kpc $Tcons_1comma5ropt $Tcons_HIH2_30kpc $Tcons_gasSF $cgf_30kpc $cgf_1comma5ropt $cgf_HIH2_30kpc $cgf_gasSF $fbar_30kpc $fbar_1comma5ropt $fbar_HIH2_30kpc $fbar_gasSF $SHMR_1comma5ropt $rhalf_stars_30kpc $rhalf_stars_30kpc_2D $F_dA_Johnson_U $F_dA_Johnson_B $F_dA_Johnson_V $F_dA_Johnson_R $F_dA_SDSS_u $F_dA_SDSS_g $F_dA_SDSS_r $F_dA_SDSS_i $F_dA_SDSS_z $MAB_Johnson_U $MAB_Johnson_B $MAB_Johnson_V $MAB_Johnson_R $MAB_dA_Johnson_U $MAB_dA_Johnson_B $MAB_dA_Johnson_V $MAB_dA_Johnson_R $mAB_Johnson_U $mAB_Johnson_B $mAB_Johnson_V $mAB_Johnson_R $mAB_dA_Johnson_U $mAB_dA_Johnson_B $mAB_dA_Johnson_V $mAB_dA_Johnson_R $mAB_SDSS_u $mAB_SDSS_g $mAB_SDSS_r $mAB_SDSS_i $mAB_SDSS_z $MAB_SDSS_u $MAB_SDSS_g $MAB_SDSS_r $MAB_SDSS_i $MAB_SDSS_z $z_mean_birth_stars $age_mean_stars $sfe_HI_30kpc $Tcons_HI_30kpc $cgf_HI_30kpc $fbar_HI_30kpc $SHMR_30kpc $SB_mu_eff_stars_1comma5ropt_B $SB_mu_eff_stars_30kpc_B $SB_mu_eff_gas_disk_r $SB_mu_opt_r $SB_mu_eff_stars_1comma5ropt_V $SB_mu_eff_stars_30kpc_V $SB_mu_eff_gas_disk_V $SB_mu_opt_V $SB_mu_eff_stars_1comma5ropt_g $SB_mu_eff_stars_30kpc_g $SB_mu_eff_gas_disk_g $SB_mu_opt_g $Sigma_gas_reff_1comma5ropt $SB_mu_1comma5ropt_B $Mgas_HIH2_30kpc_GK11 $sfe_gas_1comma5ropt $mAB_dA_total_cut_B_V $sfe_gas_reff_1comma5ropt $sfe_gas_reff_disk_1comma5ropt $flag_SB $flag_SB_B $flag_SB_r $flag_SFE $flag_sample $flag_SBplusSFE $flag_SB_envr $DtoVoidCent $delta_age_stars_rband $t50_stars $t70_stars $OH_gas_30kpc $z_grad $L_bolo $t50_halo $t70_halo $v200c $delta_tform_stars $delta_tform_halo $flag_mstar_bin $flag_mhalo_bin $flag_sfe_bin $flag_age_bin $delta_t50 $delta_lambda_edgeOn $delta_lambda_random $delta_tform_bh $t50_bh $t70_bh $delta_t50_st2ha $delta_t50_bh2ha $flag_t50_st2ha $flag_t50_bh2ha $fatom_30kpc_GK11 $tmin_rVmax $tmax_rVmax $delta_tvar_rVmax $tmin_angM_stars $tmax_angM_stars $delta_tvar_angM_stars $tmin_angM_SFgas $tmax_angM_SFgas $delta_tvar_angM_SFgas $tmin_angM_NSFgas $tmax_angM_NSFgas $delta_tvar_angM_NSFgas $angM_stars $angM_SFgas $angM_NSFgas $fmol_30kpc_GK11 $flag_SHMR_bin $tmin_sfr $tmax_sfr $tmin_ssfr $tmax_ssfr $tmin_bh_acc_rate $tmax_bh_acc_rate $tmin_vmax $tmax_vmax $tmin_vdisp $tmax_vdisp $tmin_mgas_SF $tmax_mgas_SF $tmin_mgas_NSF $tmax_mgas_NSF $delta_tvar_sfr $delta_tvar_ssfr $delta_tvar_vmax $delta_tvar_vdisp $delta_tvar_mgas_SF $delta_tvar_mgas_NSF $tmin_SHMR $tmax_SHMR $tmin_mgas $tmax_mgas $delta_tvar_SHMR $delta_tvar_mgas $mhalo_30kpc $spinStars_x $spinStars_y $spinStars_z $mgas_SF $mgas_NSF $vpec_norm $x_pos_cof $y_pos_cof $z_pos_cof $redshift $topLeafID $progFofID $delta_tvar_bh_acc_rate $vdisp_30kpc $tmin_vel $tmax_vel $delta_tvar_vel $tmin_age $tmax_age $delta_tvar_age $angM_norm_stars $angM_norm_SFgas $angM_norm_NSFgas $delta_tvar_angM_norm_stars $delta_tvar_angM_norm_SFgas $delta_tvar_angM_norm_NSFgas $tmin_angM_norm_stars $tmax_angM_norm_stars $tmin_angM_norm_SFgas $tmax_angM_norm_SFgas $tmin_angM_norm_NSFgas $tmax_angM_norm_NSFgas $flag_majorM $rVoid $t50_curve_fit_stars $t70_curve_fit_stars $t50_cs_fit_stars $t70_cs_fit_stars $t50_stars_std $t70_stars_std $delta_tform_curve_fit_stars $delta_tform_cs_fit_stars $t50_curve_fit_halo $t70_curve_fit_halo $t50_cs_fit_halo $t70_cs_fit_halo $t50_halo_std $t70_halo_std $delta_tform_curve_fit_halo $delta_tform_cs_fit_halo $t50_curve_fit_bh $t70_curve_fit_bh $t50_cs_fit_bh $t70_cs_fit_bh $t50_bh_std $t70_bh_std $delta_tform_curve_fit_bh $delta_tform_cs_fit_bh $flag_std_fit_t50_stars $flag_std_fit_t50_halo $flag_std_fit_t50_bh $delta_inter_fit_t50_stars $delta_inter_fit_t50_halo $delta_inter_fit_t50_bh $flag_std_fit_t70_stars $flag_std_fit_t70_halo $flag_std_fit_t70_bh $delta_inter_fit_t70_stars $delta_inter_fit_t70_halo $delta_inter_fit_t70_bh $flag_use_inter_stars $flag_use_inter_halo $flag_use_inter_bh $flag_use_fit_stars $flag_use_fit_halo $flag_use_fit_bh $nr_strikes_stars $nr_strikes_halo $nr_strikes_bh $delta_t50_bh2st $delta_t70_st2ha $delta_t70_bh2ha $delta_t70_bh2st $flag_t50_bh2st $flag_t70_st2ha $flag_t70_bh2ha $flag_t70_bh2st $flag_tform_stars $flag_mhalo $flag_tform_bh $flag_tform_halo $flag_dens $flag_dens3 $flag_dens6 $flag_dens3b $flag_dens6b $angM_norm_1comma5ropt_bar $delta_VvS_edgeOn $delta_vdisp_edgeOn $flag_SB_B_new $flag_dens3b_new $flag_dens6b_new $SB_mu_corr_ropt_B $kappaCoRot_30kpc_T19 $medOrbCirc_30kpc_T19 $triax_30kpc_T19)

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

set analyse_tarsel_name_array = (analyse_tarsel_dmesa_i analyse_tarsel_g-r_mstar analyse_tarsel_r-i_mstar analyse_tarsel_g-i_mstar analyse_tarsel_u-r_mstar analyse_tarsel_g-z_mstar analyse_tarsel_r-i_i analyse_tarsel_i_mstar analyse_tarsel_r_mstar analyse_tarsel_g_mstar analyse_tarsel_Ii_i analyse_tarsel_Ii_mstar analyse_tarsel_Ir_r analyse_tarsel_Ir_mstar analyse_tarsel_Ig_g analyse_tarsel_Ig_mstar analyse_tarsel_histo analyse_tarsel_sfr_mstar analyse_tarsel_g-r_r analyse_tarsel_u-r_r analyse_tarsel_mbh_mstar_spheroid analyse_tarsel_mcold_mstar 'analyse_tarsel_g-r_u-g' analyse_tarsel_Mzgas_mstar analyse_tarsel_Mzgas_mcold 'analyse_tarsel_r-i_g-r' 'analyse_tarsel_g-i_i' analyse_tarsel_dmesa_mstar analyse_tarsel_rdisk_rbulge analyse_tarsel_rhalfmass_mstar analyse_tarsel_rhalfmass_rbulge analyse_tarsel_rhalfmass_rdisk analyse_tarsel_rdisk_mstar analyse_tarsel_rhalfdisk_mstar analyse_tarsel_rhalfbulge_mstar analyse_tarsel_spinParameter_mstar analyse_tarsel_spinParameter_Mzstar analyse_tarsel_Mzstar_mstar analyse_tarsel_dmesa_ssfr analyse_tarsel_spinParameter_dmesa 'analyse_tarsel_g-r_ssfr' analyse_tarsel_ssfr_mstar analyse_tarsel_ssfr_sfr analyse_tarsel_ssfr_mhalo analyse_tarsel_zcold_mhalo analyse_tarsel_zcold_sfr analyse_tarsel_mstar_mhalo analyse_tarsel_ssfr_zcold analyse_tarsel_ssfr_mcold analyse_tarsel_ssfr_i analyse_tarsel_ssfr_r)

set analyse_tarsel_array = ($analyse_tarsel_dmesa_i $analyse_tarsel_gminr_mstar $analyse_tarsel_rmini_mstar $analyse_tarsel_gmini_mstar $analyse_tarsel_uminr_mstar $analyse_tarsel_gminz_mstar $analyse_tarsel_rmini_i $analyse_tarsel_i_mstar $analyse_tarsel_r_mstar $analyse_tarsel_g_mstar $analyse_tarsel_Ii_i $analyse_tarsel_Ii_mstar $analyse_tarsel_Ir_r $analyse_tarsel_Ir_mstar $analyse_tarsel_Ig_g $analyse_tarsel_Ig_mstar $analyse_tarsel_histo $analyse_tarsel_sfr_mstar $analyse_tarsel_gminr_r $analyse_tarsel_uminr_r $analyse_tarsel_mbh_mstar_spheroid $analyse_tarsel_mcold_mstar $analyse_tarsel_gminr_uming $analyse_tarsel_Mzgas_mstar $analyse_tarsel_Mzgas_mcold $analyse_tarsel_rmini_gminr $analyse_tarsel_gmini_i $analyse_tarsel_dmesa_mstar $analyse_tarsel_rdisk_rbulge $analyse_tarsel_rhalfmass_mstar $analyse_tarsel_rhalfmass_rbulge $analyse_tarsel_rhalfmass_rdisk $analyse_tarsel_rdisk_mstar $analyse_tarsel_rhalfdisk_mstar $analyse_tarsel_rhalfbulge_mstar $analyse_tarsel_spinParameter_mstar $analyse_tarsel_spinParameter_Mzstar $analyse_tarsel_Mzstar_mstar $analyse_tarsel_dmesa_ssfr $analyse_tarsel_spinParameter_dmesa $analyse_tarsel_gminr_ssfr $analyse_tarsel_ssfr_mstar $analyse_tarsel_ssfr_sfr $analyse_tarsel_ssfr_mhalo $analyse_tarsel_zcold_mhalo $analyse_tarsel_zcold_sfr $analyse_tarsel_mstar_mhalo $analyse_tarsel_ssfr_zcold $analyse_tarsel_ssfr_mcold $analyse_tarsel_ssfr_i $analyse_tarsel_ssfr_r)


if ($output_filename_code_default =~ CMAS*) then
	echo 'set cut values --> CMASS'
	set set_cut_values=(`$MAIN_PATH'myRun/set_cut_values_config_CMASS.sh' $MAIN_PATH "$total_col_name_array" "$total_col_array" $SOFTLINK_CODE $UNIT_CODE`)
else
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

set plot_custom_config_array=(`$MAIN_PATH'myRun/root_config.sh' $MAIN_PATH $PLOTPATH "$plot_custom_array" "$plot_custom_array" 'plotXY' ''`)

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

set MATPLOT_ALPHA =  $MAIN_PATH'myRun/input_cats/matplot_alpha.txt'

rm -f $MATPLOT_ALPHA
touch $MATPLOT_ALPHA

echo 'MATPLOT_ALPHA' $MATPLOT_ALPHA				>> $MY_PATH_HANDLER_FILE

set MATPLOT_PLOTLEGEND =  $MAIN_PATH'myRun/input_cats/matplot_plotlegend.txt'

rm -f $MATPLOT_PLOTLEGEND
touch $MATPLOT_PLOTLEGEND

echo 'MATPLOT_PLOTLEGEND' $MATPLOT_PLOTLEGEND				>> $MY_PATH_HANDLER_FILE

set MATPLOT_ADD_XAXIS =  $MAIN_PATH'myRun/input_cats/matplot_add_xaxis.txt'

rm -f $MATPLOT_ADD_XAXIS
touch $MATPLOT_ADD_XAXIS

echo 'MATPLOT_ADD_XAXIS' $MATPLOT_ADD_XAXIS				>> $MY_PATH_HANDLER_FILE

set MATPLOT_ADD_YAXIS =  $MAIN_PATH'myRun/input_cats/matplot_add_yaxis.txt'

rm -f $MATPLOT_ADD_YAXIS
touch $MATPLOT_ADD_YAXIS

echo 'MATPLOT_ADD_YAXIS' $MATPLOT_ADD_YAXIS				>> $MY_PATH_HANDLER_FILE

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
	else if ( $NAME =~ SAG*v2 ) then
		set MANUAL_REDSHIFT_INPUT = $SAG_v2_REDSHIFT_INPUT
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

	else if ($NAME =~ *100*) 	then
		set MYSNAPS 		= `echo $MYSNAPS_100Mpc`
		set BOX_SIZE 		= '100Mpc/h'
		set YLIM 		= '-5'
		set BOX_COUNT_100Mpc 	= 1

	else if ($NAME =~ *400*) 	then
		set MYSNAPS 		= `echo $MYSNAPS_400Mpc`
		set BOX_SIZE 		= '400Mpc/h'
		set YLIM 		= '-5'
		set BOX_COUNT_400Mpc 	= 1

	else if ($NAME =~ *500*) 	then
		set MYSNAPS 		= `echo $MYSNAPS_500Mpc`
		set BOX_SIZE 		= '500Mpc/h'
		set YLIM 		= '-5'
		set BOX_COUNT_500Mpc 	= 1

	else
		set MYSNAPS 		= `echo $MYSNAPS_default`
		set BOX_SIZE 		= 'None'
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

	if ($NAME =~ SAGE*) then
		if ($HDF5_read_code =~ *DF5) then
			set HDF5_CODE = 'True'
		else 
			set BINARY_CODE = 'True'
		endif
	
	else if ($HDF5_read_code =~ HDF* || $HDF5_read_code =~ SAMHDF*) then
		set HDF5_CODE = 'True'

	else if ($HDF5_read_code == False || $NAME =~ EAGLE*) then
		set HDF5_CODE = 'False'

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

	else if ($FITS_CODE == 'True') then
		echo '--> FITS'
		set CUSTOM_SNAPNAME			= 'fits'
		set SNAPIDZ_MAPPING 			= False
	else
		echo '--> ASCII'
		set CUSTOM_SNAPNAME			= 'txt'
		set SNAPIDZ_MAPPING 			= True
		if ($NAME =~ *sample* || $NAME =~ *DR* || $NAME =~ EAGLE*) then
			set MANUAL_REDSHIFT_INPUT	= $MANUAL_REDSHIFT_INPUT_default
			set SNAPIDZ_MAPPING 		= False
		endif
		if ($NAME =~ EAGLE*) then
			#set CUSTOM_SNAPNAME			= 'EAGLE_ASCII_full_trees'
			set CUSTOM_SNAPNAME			= 'EAGLE_ASCII_full'
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

	if ($NAME =~ SAG* || $NAME =~ suss* || $NAME =~ LGAL* || $NAME =~ EAGLE*) then
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
		else if ($NAME =~ *DR* || $NAME =~ *sample* || $NAME =~ CMAS*) then
			set SNAPID 		= $ISNAP
		else
			set SNAPID 		= `echo $ISNAP | awk '{printf ("%03d\n",$1)}'`
		endif

		set c=1
		set k=1
		set x=1
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
				else if ($PLOT_CUSTOM_LOOP == True && $plotXY == 1) then

					foreach plot_custom_item ($plot_custom_array)
						echo redshift$a'= '$SNAPID >> "$plot_custom_config_array[$x]"
						@ x++
				else
					echo redshift$a'= '$SNAPID >> "$plot_config_array[$c]"

				endif
			endif
			@ c++
		end


		if ( $CUSTOM_SNAPNAME == 'hdf5' || $CUSTOM_SNAPNAME == 'binary') then  
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
						echo 'File:' $File
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

			else if ($NAME =~ EAGL*) then
				echo 'EAGLE: $SOFTLINK_TO_DATA ~$SNAPID' $SOFTLINK_TO_DATA $NAME $SNAPID
				foreach File (`ls -d $SOFTLINK_TO_DATA/*$SNAPID*`)
					echo 'File:' $File
					echo $File >> $FILENAME_LIST
				end


			else if ($NAME =~ suss*) then
				#echo '$SOFTLINK_TO_DATA ~$SNAPID' $SOFTLINK_TO_DATA $NAME $SNAPID
				foreach File (`ls -d $SOFTLINK_TO_DATA*suss*`)
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
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_EAGLE_cross.sh', shell=True)
			echo 'SNAP_LIST:' $SNAP_LIST 'SNAPNAME:' $SNAPNAME

			if ($SNAPNAME =~ *config* || $SNAPNAME =~ *idzred* || $SNAPNAME =~ *READ*) then
				echo $SNAPNAME': NOT A SNAPSHOT FILE'
			else
				echo  $SNAPNAME 			>> $RUNFILE
				set b = `expr $a \* $NR_MODELS`
		
				set z_i = `expr $j + $b`

				set c=1
				set k=1
				set x=1
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
					else if ($PLOT_CUSTOM_LOOP == True && $plotXY == 1) then
						foreach plot_custom_item ($plot_custom_array)
							echo 'plot_legend'$z_i'= '$NAME >> "$plot_custom_config_array[$x]"
							@ x++
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
					echo 'fileformat= ""'				>> $CONFIGFILE
				else if ($NAME =~ EAGLE* && $CUSTOM_SNAPNAME == hdf5 ) then
					echo 'data_format= HDF5'	 		>> $CONFIGFILE
					echo 'format_info= choose_all'     		>> $CONFIGFILE
					echo 'delimiter= ""'		     	 	>> $CONFIGFILE
					echo 'fileformat= hdf5'		>> $CONFIGFILE
				else if ($NAME =~ EAGLE* && $CUSTOM_SNAPNAME =~ EAGLE_ASCI*) then
					echo 'data_format= '$CUSTOM_SNAPNAME  		>> $CONFIGFILE
					echo 'format_info= choose_all'     		>> $CONFIGFILE
					echo 'delimiter= ""'		     	 	>> $CONFIGFILE
					echo 'fileformat= dat'				>> $CONFIGFILE					
				else if ($CUSTOM_SNAPNAME =~ *its) then
					echo 'data_format= '$HDF5_read_code  		>> $CONFIGFILE
					echo 'format_info= choose_all'     		>> $CONFIGFILE
					echo 'delimiter= \t'		     	 	>> $CONFIGFILE
					echo 'fileformat= ""'			>> $CONFIGFILE
				else if ($SNAPNAME =~ *dat || $CUSTOM_SNAPNAME =~ *dat) then
					echo 'data_format= BINARY'              	>> $CONFIGFILE
					echo 'format_info= choose_all'     		>> $CONFIGFILE
					echo 'delimiter= \t'		     	 	>> $CONFIGFILE
					echo 'fileformat= niftydat'			>> $CONFIGFILE
				else if ($SNAPNAME =~ *xt || $CUSTOM_SNAPNAME =~ *xt) then
					echo 'data_format= CATASCII' 		 	>> $CONFIGFILE
					echo 'format_info= shaped'	     	>> $CONFIGFILE
					echo 'delimiter= ""'		 	>> $CONFIGFILE
					echo 'fileformat= txt'			>> $CONFIGFILE
				else if ($SNAPNAME =~ *MDPL*) then
					echo 'data_format= ASCII'               	>> $CONFIGFILE
					echo 'format_info= unshaped'	     	>> $CONFIGFILE
					echo 'delimiter= 	'		 	>> $CONFIGFILE
					echo 'fileformat= txt'				>> $CONFIGFILE
				else if ($SNAPNAME =~ *tree*) then
					echo 'data_format= SAMHDF5'               	>> $CONFIGFILE
					echo 'format_info= unshaped'	     	>> $CONFIGFILE
					echo 'delimiter= ""'		 	>> $CONFIGFILE
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
				if ($NAME =~ *125* || $NAME =~ *100*) then
					echo 'nr_rows= 550000'		>> $CONFIGFILE
				else if ($NAME =~ *400*) then
					echo 'nr_rows= 200000000'   	>> $CONFIGFILE
				else if ($NAME =~ *62* || $HOME_MODE == True || $HOME =~ /home/*) then
					echo 'nr_rows= 5000000'   >> $CONFIGFILE
				else
					echo 'nr_rows= 200000000'   >> $CONFIGFILE
					echo 'nr_rows= 10000'   >> $CONFIGFILE
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
				if ($NAME =~ *62*) 		echo 'box_size= 62.5'  	>> $CONFIGFILE
				if ($NAME =~ *100*) 		echo 'box_size= 100.0'  >> $CONFIGFILE
				if ($NAME =~ *sub) 		echo 'box_size= 112.0' 	>> $CONFIGFILE
				if ($NAME =~ *125*) 		echo 'box_size= 125.0'  >> $CONFIGFILE
				if ($NAME =~ *400*) 		echo 'box_size= 400.0' 	>> $CONFIGFILE
				if ($NAME =~ *500*) 		echo 'box_size= 500.0' 	>> $CONFIGFILE
				if ($NAME =~ *1Gpc*) 		echo 'box_size= 1000.0' >> $CONFIGFILE

				
				if ($SNAPIDZ_MAPPING == 'True') then

					set INPUTFILENAME_PART1 = `ls $MODEL/*$SNAPID* | head -1 | awk -F$SNAPID '{print $1}'`		
					set INPUTFILENAME_PART2 = `ls $MODEL/*$SNAPID* | tail -1 | awk -F$SNAPID '{print $2}'`

					if ($INPUTFILENAME_PART1 =~ *'/') then 
						set INPUTFILENAME_PART1 = ''
					else
						set INPUTFILENAME_PART1 = `basename $INPUTFILENAME_PART1`
					endif
				else
					if ($CUSTOM_SNAPNAME == 'hdf5' && $NAME =~ SAG_* && $PATH_TO_CATALOG$NAME'/' == $SOFTLINK_TO_DATA || $NAME =~ *DR*) then
						set INPUTFILENAME_PART1 = `basename $key`
					else
						set INPUTFILENAME_PART1 = $key
					endif

					set INPUTFILENAME_PART2 = '/'
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
				echo 'plot_cusom_loop= '$PLOT_CUSTOM_LOOP	   		>> $CONFIGFILE
				echo 'sample_info= '$sample_info	   	 	>> $CONFIGFILE
				echo 'ngal_random_sample= '$ngalaxies_random		>> $CONFIGFILE
				echo 'fits_map_names= '$fits_map_names			>> $CONFIGFILE
				echo 'fits_map_names_mapping= '$fits_map_names_mapping	>> $CONFIGFILE
				echo 'UNIT_CODE= '$UNIT_CODE		   	 	>> $CONFIGFILE
				echo 'resolution= '$resolution		   	 	>> $CONFIGFILE
				echo 'cosmology= '$cosmology		   	 	>> $CONFIGFILE
				echo 'IMF= '$IMF			   	 	>> $CONFIGFILE
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
					set tarsel_code = $tarsel_code_Galacticus
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

				echo 'set_header_Topcat_format= '$set_header_Topcat_format			 	>> $CONFIGFILE
				echo 'change_IMF= '$change_IMF			 	>> $CONFIGFILE
				echo 'which_kcorrect= '$which_kcorrect			>> $CONFIGFILE
				echo 'calc_fast_histo= '$calc_fast_histo 	 	>> $CONFIGFILE
	
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

				#analyseTargetSelection specifics
				echo 'filterSample= '$filterSample 			>> $CONFIGFILE
				echo 'calcHistoSample= '$calcHistoSample		>> $CONFIGFILE
				echo 'plotSample= '$plotSample				>> $CONFIGFILE

				echo 'OUTPUT_ASCII= '$HDF5_TO_ASCII				>> $CONFIGFILE

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
set x=1

if ($tarsel_code_Galacticus == '') set tarsel_code_Galacticus = '-'
if ($tarsel_code_SAG == '') set tarsel_code_SAG = '-'
if ($tarsel_code_default == '') set tarsel_code_default = '-'

foreach item ($method_array)
	#echo "$method_name_array[$c]" $item
	if ("$method_name_array[$c]" == plotOnly && $item == 1) then
		foreach item2 ($total_col_name_array)

			if ("$total_col_array[$i]" != 999) then
				#echo 1291 $item2 "$total_col_array[$i]" "$method_name_array[$c]" $tarsel_code "$plot_plotOnly_config_array[$a]"
				
				$MAIN_PATH'myRun/plot_config.sh' "plotOnly" "$plot_plotOnly_config_array[$a]" $BOX_SIZE $CALI_CODE $item2 $CAT_NAME_IN_PLOT0 $CAT_NAME_IN_PLOT1 $REDSHIFT_PLOT0 $REDSHIFT_PLOT1 $SOFTLINK_CODE $a $tarsel_code_Galacticus $tarsel_code_SAG $MATPLOT_LINESTYLE $MATPLOT_COL $MATPLOT_MARKERSTYLE $MATPLOT_MARKERCOL $MATPLOT_COLORMAP $MATPLOT_COLORBAR $MATPLOT_PLOTTYPE $MATPLOT_ALPHA $MATPLOT_PLOTLEGEND $MATPLOT_ADD_XAXIS $MATPLOT_ADD_YAXIS $PLOT_CUSTOM_LOOP
				@ a++
			endif
			@ i++
				
		end

	else if ("$method_name_array[$c]" == analyseTargetSelection && $item == 1) then
			foreach item2 ($analyse_tarsel_name_array)

				if ("$analyse_tarsel_array[$i]" == 1) then
					#echo 1421 $item2 "$analyse_tarsel_name_array[$i]"  $tarsel_code "$plot_analyseTargetSelection_config_array[$a]"
				
				$MAIN_PATH'myRun/plot_config.sh' "$analyse_tarsel_name_array[$i]" "$plot_analyseTargetSelection_config_array[$a]" $BOX_SIZE $CALI_CODE $CAT_NAME_IN_PLOT0 $CAT_NAME_IN_PLOT1 $REDSHIFT_PLOT0 $REDSHIFT_PLOT1 $SOFTLINK_CODE $c $tarsel_code_Galacticus $tarsel_code_SAG $MATPLOT_LINESTYLE $MATPLOT_COL $MATPLOT_MARKERSTYLE $MATPLOT_MARKERCOL $MATPLOT_COLORMAP $MATPLOT_COLORBAR $MATPLOT_PLOTTYPE $MATPLOT_ALPHA $MATPLOT_PLOTLEGEND $MATPLOT_ADD_XAXIS $MATPLOT_ADD_YAXIS $PLOT_CUSTOM_LOOP
				@ a++
			endif
			@ i++
				
		end

	

	else if ( $PLOT_CUSTOM_LOOP == True && $item == 1) then
		foreach custom_item ($plot_custom_array)
			$MAIN_PATH'myRun/plot_config.sh' 'plotXY' "$plot_custom_config_array[$x]" $BOX_SIZE $CALI_CODE $CAT_NAME_IN_PLOT0 $CAT_NAME_IN_PLOT1 $REDSHIFT_PLOT0 $REDSHIFT_PLOT1 $SOFTLINK_CODE $c $tarsel_code_Galacticus $tarsel_code_SAG $MATPLOT_LINESTYLE $MATPLOT_COL $MATPLOT_MARKERSTYLE $MATPLOT_MARKERCOL $MATPLOT_COLORMAP $MATPLOT_COLORBAR $MATPLOT_PLOTTYPE $MATPLOT_ALPHA $MATPLOT_PLOTLEGEND $MATPLOT_ADD_XAXIS $MATPLOT_ADD_YAXIS $custom_item
			@ x++

	else
	
		if ($item == 1)  then
			$MAIN_PATH'myRun/plot_config.sh' "$method_name_array[$c]" "$plot_config_array[$c]" $BOX_SIZE $CALI_CODE $CAT_NAME_IN_PLOT0 $CAT_NAME_IN_PLOT1 $REDSHIFT_PLOT0 $REDSHIFT_PLOT1 $SOFTLINK_CODE $c $tarsel_code_Galacticus $tarsel_code_SAG $MATPLOT_LINESTYLE $MATPLOT_COL $MATPLOT_MARKERSTYLE $MATPLOT_MARKERCOL $MATPLOT_COLORMAP $MATPLOT_COLORBAR $MATPLOT_PLOTTYPE $MATPLOT_ALPHA $MATPLOT_PLOTLEGEND $MATPLOT_ADD_XAXIS $MATPLOT_ADD_YAXIS $PLOT_CUSTOM_LOOP
			
		endif
	
	endif
	@ c++
end


echo 'GENERATING CONFIG FILES ... DONE!'
echo ' '
