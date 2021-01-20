#!/bin/tcsh

set MYVALUES_CONFIGFILE		=  $1'myRun/run_config/myvalues_config.txt'

rm -f $MYVALUES_CONFIGFILE 
touch $MYVALUES_CONFIGFILE

set MY_UNIT_CORRECT_FILE		=  $1'myRun/run_config/myunit_correct_file_config.txt'

rm -f $MY_UNIT_CORRECT_FILE 
touch $MY_UNIT_CORRECT_FILE

set MY_NAME_CONV_MAP_FILE		=  $1'myRun/run_config/myname_conv_map_file.txt'

rm -f $MY_NAME_CONV_MAP_FILE 
touch $MY_NAME_CONV_MAP_FILE

set array = ()
foreach item ($3)
	set array = ($array $item)
end

set i=1
foreach item ($2)
	if ("$array[$i]" != 99) then

		#set the columns how the should be named in the output-file
		#if False: DS code will be used
		#if not False, this name will be used to store in the output-files
		set output_col_name = $item
		if ($5 == 'MD') then
			if ($item == hostid) set output_col_name	= 'HostHaloID'			
			if ($item == haloid) set output_col_name	= 'MainHaloID'		
			if ($item == orphan) set output_col_name	= 'GalaxyType'
			if ($item == mhalo) set output_col_name		= 'HaloMass'
			if ($item == mhalo_200c) set output_col_name	= 'HaloMass_200c'
			if ($item == mhalo_cents) set output_col_name	= 'HaloMass_cents'
			if ($item == mhalo_cents_200c) set output_col_name	= 'HaloMass_cents_200c'
			if ($item == mbasic) set output_col_name	= 'HaloMass_basic'
			if ($item == mbasic_200c) set output_col_name	= 'HaloMass_basic_200c'
			if ($item == vmax) set output_col_name		= 'Vmax'
			if ($item == vpeak) set output_col_name		= 'Vpeak'		
			if ($item == NFW_con) set output_col_name	= 'NFWconcentration'
			if ($item == spinParameter) set output_col_name	= 'SpinParameter'

			if ($item == x_pos) set output_col_name		= 'X'
			if ($item == y_pos) set output_col_name		= 'Y'		
			if ($item == z_pos) set output_col_name		= 'Z'
			if ($item == x_vel) set output_col_name		= 'Vx'
			if ($item == y_vel) set output_col_name		= 'Vy'
			if ($item == z_vel) set output_col_name		= 'Vz'
		
			if ($item == mstar) set output_col_name		= 'Mstar'
			if ($item == mstar_spheroid) set output_col_name = 'MstarSpheroid'
			if ($item == mstar_disk) set output_col_name	= 'MstarDisk'

			if ($item == mcold) set output_col_name		= 'Mcold'
			if ($item == mcold_spheroid) set output_col_name = 'McoldSpheroid'
			if ($item == mcold_disk) set output_col_name	= 'McoldDisk'
			if ($item == mstar_IC) set output_col_name	= 'MstarIC'

			if ($item == mhot) set output_col_name		= 'Mhot'
			if ($item == mhot_outflow) set output_col_name		= 'MhotOutflow'
			if ($item == mbh) set output_col_name 		= 'Mbh'

			if ($item == mstar) set output_col_name		= 'Mstar'
			if ($item == mstar_spheroid) set output_col_name = 'MstarSpheroid'
			if ($item == mstar_disk) set output_col_name	= 'MstarDisk'

			if ($item == mcold) set output_col_name		= 'Mcold'
			if ($item == mcold_spheroid) set output_col_name = 'McoldSpheroid'
			if ($item == mcold_disk) set output_col_name	= 'McoldDisk'

			if ($item == sfr) set output_col_name		= 'SFR'
			if ($item == sfr_spheroid) set output_col_name  = 'SFRspheroid'
			if ($item == sfr_disk) set output_col_name	= 'SFRdisk'

			if ($item == L_SDSS_dA_total_u) set output_col_name		= 'LstarSDSSu'
			if ($item == L_SDSS_dA_total_g) set output_col_name		= 'LstarSDSSg'
			if ($item == L_SDSS_dA_total_r) set output_col_name		= 'LstarSDSSr'
			if ($item == L_SDSS_dA_total_i) set output_col_name		= 'LstarSDSSi'
			if ($item == L_SDSS_dA_total_z) set output_col_name		= 'LstarSDSSz'

			if ($item == L_SDSS_dA_spheroid_u) set output_col_name		= 'LbulgeSDSSu'
			if ($item == L_SDSS_dA_spheroid_g) set output_col_name		= 'LbulgeSDSSg'
			if ($item == L_SDSS_dA_spheroid_r) set output_col_name		= 'LbulgeSDSSr'
			if ($item == L_SDSS_dA_spheroid_i) set output_col_name		= 'LbulgeSDSSi'
			if ($item == L_SDSS_dA_spheroid_z) set output_col_name		= 'LbulgeSDSSz'

			if ($item == L_SDSS_dA_disk_u) set output_col_name		= 'LdiskSDSSu'
			if ($item == L_SDSS_dA_disk_g) set output_col_name		= 'LdiskSDSSg'
			if ($item == L_SDSS_dA_disk_r) set output_col_name		= 'LdiskSDSSr'
			if ($item == L_SDSS_dA_disk_i) set output_col_name		= 'LdiskSDSSi'
			if ($item == L_SDSS_dA_disk_z) set output_col_name		= 'LdiskSDSSz'

			if ($item == L_RGO_dA_total_B) set output_col_name		= 'LstarRGOB'

			if ($item == MAB_dA_total_u) set output_col_name		= 'MagStarSDSSu'
			if ($item == MAB_dA_total_g) set output_col_name		= 'MagStarSDSSg'
			if ($item == MAB_dA_total_r) set output_col_name		= 'MagStarSDSSr'
			if ($item == MAB_dA_total_i) set output_col_name		= 'MagStarSDSSi'
			if ($item == MAB_dA_total_z) set output_col_name		= 'MagStarSDSSz'


			if ($item == MAB_bulge_u) set output_col_name		= 'MagBulgeSDSSu'
			if ($item == MAB_bulge_g) set output_col_name		= 'MagBulgeSDSSg'
			if ($item == MAB_bulge_r) set output_col_name		= 'MagBulgeSDSSr'
			if ($item == MAB_bulge_i) set output_col_name		= 'MagBulgeSDSSi'
			if ($item == MAB_bulge_z) set output_col_name		= 'MagBulgeSDSSz'

			if ($item == Mag_dA_total_B) set output_col_name		= 'MagStarRGOB'

			if ($item == mAB_dA_total_u) set output_col_name		= 'mABStarSDSSu'
			if ($item == mAB_dA_total_g) set output_col_name		= 'mABStarSDSSg'
			if ($item == mAB_dA_total_r) set output_col_name		= 'mABStarSDSSr'
			if ($item == mAB_dA_total_i) set output_col_name		= 'mABStarSDSSi'
			if ($item == mAB_dA_total_z) set output_col_name		= 'mABStarSDSSz'

			if ($item == Mzgas) set output_col_name  	 = 'MZgas'
			if ($item == Mzgas_spheroid) set output_col_name = 'MZgasSpheroid'
			if ($item == Mzgas_disk) set output_col_name	 = 'MZgasDisk'

			if ($item == Mzstar) set output_col_name  	 = 'MZstar'
			if ($item == Mzstar_spheroid) set output_col_name = 'MZstarSpheroid'
			if ($item == Mzstar_disk) set output_col_name	 = 'MZstarDisk'

			if ($item == Mzhot_halo) set output_col_name	 = 'MZhotHalo'

			if ($4 == A || $4 =~ A* || $4 == C) then
			#Galacticus: wrong name in the columns for z~, those are actually masses of metals!!!
				if ($item == zstar_spheroid) set output_col_name = 'MZstarSpheroid'
				if ($item == zstar_disk) set output_col_name	 = 'MZstarDisk'
				if ($item == zgas_spheroid) set output_col_name  = 'MZgasSpheroid'
				if ($item == zgas_disk) set output_col_name	 = 'MZgasDisk'
				if ($item == zhot_halo) set output_col_name	 = 'MZhotHalo'
			else
				if ($item == zstar_spheroid) set output_col_name = 'ZstarSpheroid'
				if ($item == zstar_disk) set output_col_name	 = 'ZstarDisk'
				if ($item == zgas_spheroid) set output_col_name  = 'ZgasSpheroid'
				if ($item == zgas_disk) set output_col_name  	 = 'ZgasDisk'
				if ($item == zgas_hot) set output_col_name	 = 'ZhotHalo'
			endif

			if ($item == mean_age_stars) set output_col_name = 'MeanAgeStars'
			if ($item == Z) set output_col_name = 'redshift'

		else

			if ($4 == A || $4 =~ A* || $4 == C) then
			#Galacticus: wrong name in the columns for z~, those are actually masses of metals!!!
				if ($item == zstar_spheroid) set output_col_name = 'MZstar_spheroid'
				if ($item == zstar_disk) set output_col_name	 = 'MZstar_disk'
				if ($item == zgas_spheroid) set output_col_name  = 'MZgas_spheroid'
				if ($item == zgas_disk) set output_col_name	 = 'MZgas_disk'
				if ($item == zhot_halo) set output_col_name	 = 'MZhot_halo'
			endif	
		endif



		echo $item $output_col_name	>> $MY_NAME_CONV_MAP_FILE

#CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN
		#MIN is smaller or equal (<=)
		#default values: 
		set cut_value_min		= 'min'
		#echo 'i:' $i $item "$array[$i]"
		 		 
		#if ($item == mhalo) set cut_value_min		= 1e3
		#if ($item == hostid) set cut_value_min		= -2

		#if ($item == mstar) set cut_value_min	 	= 5e7
		#if ($item == mstar) set cut_value_min	 	= 1e7
		#if ($item == zcold) set cut_value_min	 	= 8.5
		#if ($item == sfr) set cut_value_min	 	= 1e-4
		#if ($item == 'mstar+IC') set cut_value_min	 	= 5e10
		#if ($item == MAB_dA_total_r) set cut_value_min	= -100.0
		#if ($item == MAB_total_r) set cut_value_min	= -21.0
		#if ($item =~ L_*) set cut_value_min	= 1.0
		#if ($item =~ *rphan) set cut_value_min		= 1
		#if ($item == pop) set cut_value_min		= 1.1
		#CMASS
		#if ($item == mAB_dA_total_i) set cut_value_min	= 17.5
		#if ($item =~ *_cut_dmesa) set cut_value_min	= 0.55
		#if ($item =~ *_cut_i_lt_dmes*) set cut_value_min = 1.0

		#if ($item == mAB_total_i) set cut_value_min	= 17.5


#CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN -- CUT MIN



#CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX
		#max is greater >
		#default values: 
		set cut_value_max		= 'max'
		#if ($item == Z ) set cut_value_max		= 0.6		 
		#if ($item =~ *rphan) set cut_value_max		= 0
		#if ($item == pop) set cut_value_max		= 1.1
		#if ($item == mstar) set cut_value_max		= 1e11
		#if ($item == MAB_dA_total_u) set cut_value_min	= 'inf'
		#if ($item == MAB_dA_total_r) set cut_value_max	= 100.0
		#if ($item == MAB_total_r) set cut_value_max	= -100.0

		#if ($item == zcold) set cut_value_max	= 8.5

		#if ($item =~ *pos) set cut_value_max	= 516.45
		#CMASS
		#if ($item == mAB_dA_total_i) set cut_value_max	= 19.9
		#if ($item =~ *_cut_r_i) set cut_value_max	= 2.0

		#if ($item == boss_target1) set cut_value_max	= 127
		#if ($item == warning) set cut_value_max	= 0

		#if ($item == mAB_total_i) set cut_value_max	= 19.9

#CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX -- CUT MAX




		#default values:
		set unit			= '-'
		#if ($item =~ *n) set unit			= 'idnumber'
		if ($item =~ *ngal*) set unit			= 'idnumber'
		if ($item =~ *nhalo*) set unit			= 'idnumber'
		if ($item =~ *aloi*) set unit			= 'idnumber'
		if ($item =~ *Index) set unit			= 'idnumber'
		if ($item =~ *osti*) set unit			= 'idnumber'
		if ($item =~ *ID) set unit				= 'idnumber'		 					 
		if ($item =~ *rpha*) set unit			= '0=cent,1=sat,2=orph'
		if ($item == isolated) set unit			= '0=non-iso,1=iso'
		if ($item =~ *targe*) set unit			= 'flag'
		if ($item == debugR) set unit			= 'debugQuant'


		if ($5 == 'DS') then
		#units without any little-h according to Doris datapipline mass [Msun], sfr [Msun/yr], pos [comMpc], mean_age_stars [Gyr]			 		 
			if ($item =~ mhal*) set unit			= 'Msun'
			if ($item =~ mpseud*) set unit			= 'Msun'
			if ($item =~ mbasi*) set unit			= 'Msun'
			if ($item =~ Mz*) set unit			= 'Msun'				 		
			if ($item =~ mcol*) set unit			= 'Msun'
			if ($item =~ mhot*) set unit			= 'Msun'		
			if ($item == mbh) set unit			= 'Msun'			
			if ($item =~ msta*) set unit	 		= 'Msun'
			if ($item == mbasic) set unit	 		= 'Msun'
			if ($item == delta_mhalo) set unit	= 'Msun'		 
			if ($item =~ sf*) set unit			= 'Msunyr-1'
			if ($item =~ bh_*) set unit			= 'Msunyr-1'
			if ($item =~ ssf*) set unit			= 'yr-1'
			if ($item =~ age_*) set unit			= '-'
			if ($item =~ mean_age_*) set unit		= 'Gyr'			
			if ($item =~ *pos) set unit			= 'comvMpc'
			if ($item =~ r*) set unit			= 'Mpc'
			if ($item =~ *_a) set unit			= 'Mpc'
			if ($item =~ *_a*) set unit			= 'Mpc'
			if ($item == xoff) set unit			= 'Mpc'
			if ($item == delta_rvir) set unit	= 'Mpc'
			if ($item =~ OII*) set unit			= '1e40ergs-1'
			if ($item =~ OII_cont*) set unit		= '1e40ergs-1A-1'

			if ($item =~ ang*) set unit			= 'MsunMpckms-1'
			if ($item =~ *_ang) set unit		= 'MsunMpckms-1'

			if ($4 == A || $4 =~ A* || $4 == C) then
			#Galacticus: wrong name in the columns for z~, those are actually masses of metals!!!
				if ($item =~ L_*) set unit		= '4.4659e13WHz-1'
			endif		
		else
		#MD: units according to Multidark data base mass [Msun/h], sfr [Msun/Gyr/h], pos [comMpc/h], mean_age_stars [Gyr], L_* (luminosities) [4.4659e13WHz-1]
		#CM: units according to CMASS paper data base mass [Msun/h], sfr [Msun/yr/h], pos [comMpc/h], mean_age_stars [Gyr], L_* (luminosities) [4.4659e13WHz-1]
			if ($item =~ mhal*) set unit			= 'h-1Msun'
			if ($item =~ mbasi*) set unit			= 'h-1Msun'
			if ($item =~ Mz*) set unit			= 'h-1Msun'			
			if ($item =~ mcol*) set unit			= 'h-1Msun'
			if ($item =~ mhot*) set unit			= 'h-1Msun'		
			if ($item == mbh) set unit			= 'h-1Msun'			
			if ($item =~ msta*) set unit	 		= 'h-1Msun'
			if ($item == mbasic) set unit	 		= 'h-1Msun'

			if ($5 == 'CM') then
				if ($item =~ msta*) set unit		= 'h-1Msun'
				if ($item =~ sf*) set unit		= 'h-1Msunyr-1'
				if ($item =~ ssf*) set unit		= 'Gyr-1'
			else	 
				if ($item =~ sf*) set unit		= 'h-1MsunGyr-1'
			endif
			if ($item =~ bh_acc*) set unit			= 'h-1MsunGyr-1'		
			if ($item =~ *pos) set unit			= 'h-1comvMpc'
			if ($item =~ r*) set unit			= 'h-1Mpc'
			if ($item =~ OII*) set unit			= '1e40h-2ergs-1'
			if ($item =~ OII_cont*) set unit		= '1e40h-2ergs-1A-1'

			if ($item =~ ang*) set unit			= 'h-1Msunkpckms-1'

			if ($4 == A || $4 =~ A* || $4 == C) then
			#Galacticus: wrong name in the columns for z~, those are actually masses of metals!!!
				if ($item =~ zgas*) set unit			= 'h-1Msun'
				if ($item =~ zsta*) set unit			= 'h-1Msun'
				if ($item =~ zho*) set unit			= 'h-1Msun'
				if ($item =~ L_*) set unit			= '4.4659e13WHz-1'
			endif
		endif

		if ($item =~ env_*) set unit		= '0=vo,1=sh,2=fi,3=no'
		if ($item == pop) set unit			= '-'
		if ($item =~ zga*) set unit			= 'abFrac'
		if ($item =~ zsta*) set unit		= 'abFrac'
		if ($item =~ zho*) set unit			= 'abFrac'
		if ($item == cgf) set unit			= '-'
		if ($item == fbary) set unit		= '-'
		if ($item =~ *vs*) set unit			= '-'
		if ($item =~ zcol*) set unit		= '-'
		if ($item == zstar) set unit		= '-'
		if ($item == bh_rad_eff) set unit	= 'frac'

		if ($item =~ *sca) set unit			= 'Mpc'					 
		if ($item == RA) set unit			= 'deg'
		if ($item == DEC) set unit			= 'deg'
		if ($item =~ time*) set unit		= 'Gyr'
		if ($item =~ *MergeTime) set unit	= 'Gyr'
		if ($item == Tcons) set unit		= 'Gyr'
		if ($item =~ v*) set unit			= 'kms-1'
		if ($item =~ *vel) set unit			= 'kms-1'
		if ($item =~ *cut_i_lt_dmesa) set unit		= '0=False,1=True'
		if ($item =~ *cut_r_lt_cpar) set unit		= '0=False,1=True'

		if ($4 == I) then
			if ($item == orphan) set unit			= '1=cents,0=sats'
		endif			


		#default values:
		set format			= '%0.10E'
		if ($item == pop) set format			= '%i'
		if ($item =~ env_*) set format			= '%i'
		if ($item =~ *ngal*) set format			= '%i'
		if ($item =~ npr*) set format			= '%i'
		if ($item =~ *nhalo*) set format		= '%i'
		if ($item =~ *aloi*) set format			= '%i'
		if ($item =~ *osti*) set format			= '%i'
		if ($item =~ *Index) set format			= '%i'
		if ($item =~ *rpha*) set format			= '%i'
		if ($item =~ *ID) set format			= '%i'
		if ($item == isolated) set format		= '%i'
		if ($item =~ *targe*) set format		= '%i'
		if ($item =~ *arning*) set format		= '%i'
		if ($item =~ *primar*) set format		= '%i'				 		 
		if ($item == RA) set format			= '%0.10f'
		if ($item == DEC) set format			= '%0.10f'
		if ($item == Z) set format			= '%0.10f'			
		if ($item =~ *pos) set format			= '%f'
		if ($item =~ *vel) set format			= '%0.6E'		
		if ($item =~ mA*) set format			= '%0.7f'
		if ($item =~ MA*) set format			= '%0.7f'
		if ($item =~ kcorr*) set format			= '%0.7f'
		if ($item == NFW_con) set format		= '%0.4f'
		if ($item =~ *cut_*) set format			= '%0.7f'
		if ($item =~ *cut_i_lt_dmesa) set format	= '%d'
		if ($item == bhcount) set format		= '%i'
		if ($item =~ ang*) set format			= '%0.10e'
		if ($item == vbulge) set format			= '%0.10e'
		if ($item == vdisk) set format			= '%0.10e'
		if ($item =~ zcol*) set format			= '%0.10f'
		if ($item == zstar) set format			= '%0.10f'
		if ($item == Tcons) set format			= '%0.10f'

		#default values:
		if ($4 == F || $4 == B || $4 == Ff || $4 == Bb || $4 == I || $4 == R || $4 == C256 || $4 == C512) then
			set data_type			= 'float32'
		else
			set data_type			= 'float64'
		endif 				 		 

		if ($item =~ n*) set data_type			= 'int32'
		if ($item =~ *nhal*) set data_type		= 'int32'
		if ($item =~ npr*) set data_type		= 'int32'
		if ($item =~ *ngal*) set data_type		= 'int32'
		if ($item =~ *aloi*) set data_type		= 'int64'
		if ($item =~ *osti*) set data_type		= 'int64'
		if ($item =~ *Index) set data_type		= 'int64'
		if ($item =~ *ID) set data_type			= 'int64'
		if ($item == bhcount) set data_type		= 'int8'
		if ($item =~ *rpha*) set data_type		= 'int8'
		if ($item == isolated) set data_type	= 'int8'
		if ($item =~ *targe*) set data_type		= 'int8'



		set corr_type		= 'no_corr'
		if ($item =~ n*) set corr_type			= 'int_num'
		if ($item =~ *nhal*) set corr_type		= 'int_num'
		if ($item =~ npr*) set corr_type		= 'int_num'
		if ($item =~ *ngal*) set corr_type		= 'int_num'
		if ($item =~ *targe*) set corr_type		= 'int_num'
		if ($item =~ *aloi*) set corr_type		= 'id_num'
		if ($item =~ *Index) set corr_type		= 'id_num'
		if ($item =~ *ID) set corr_type			= 'id_num'			 
		if ($item =~ *osti*) set corr_type		= 'id_num'
		if ($item =~ *rpha*) set corr_type		= 'id_num'
		if ($item == isolated) set corr_type		= 'id_num'

		if ($5 == 'DS') then
		#units without any little-h according to Doris datapipline mass [Msun], sfr [Msun yr-1], pos [comMpc], mean_age_stars [Gyr], radius [Mpc], angular momentum [Msun Mpc km s-1]
			
			if ($4 == A || $4 =~ A* || $4 == C) then
			#Galacticus 				 		 
				if ($item =~ sf*) set corr_type			= 'Gyr-1'
				if ($item =~ bh_acc*) set corr_type		= 'Gyr-1'		
				if ($item =~ *pos) set corr_type		= 'COMV'
				if ($item =~ ang*) set corr_type		= 'Mpc'		
		
			else if ($4 == B || $4 == E) then
			#SAG
				if ($item =~ mhal*) set corr_type		= 'h'
				if ($item =~ Mz*) set corr_type 		= 'h'	 		
				if ($item =~ mcol*) set corr_type		= 'h'			
				if ($item == mbh) set corr_type			= 'h'			
				if ($item =~ msta*) set corr_type 		= 'h'
				if ($item =~ mho*) set corr_type 		= 'h'
				if ($item =~ OII*) set corr_type		= 'h-2'
		
				if ($item =~ *pos) set corr_type		= 'Mpc,h'
				if ($item =~ rhalf*) set corr_type		= 'Mpc,h'
				if ($item == r200c) set corr_type		= 'Mpc,h'	
				if ($item =~ sf*) set corr_type			= 'h'

			else if ($4 == F) then
			#SAGE
				if ($item =~ mhal*) set corr_type		= 'mass,h'
				if ($item =~ Mz*) set corr_type 		= 'mass,h'	 		
				if ($item =~ mcol*) set corr_type		= 'mass,h'			
				if ($item == mbh) set corr_type			= 'mass,h'			
				if ($item =~ msta*) set corr_type 		= 'mass,h'
				if ($item =~ mho*) set corr_type 		= 'mass,h'

				if ($item =~ *pos) set corr_type		= 'h'
				if ($item =~ r*) set corr_type			= 'h'
				if ($item =~ mean_age_stars) set corr_type	= '1e3,h-1'
			
			else if ($4 == L) then
			#LGALAXIES 500 Mpc
				if ($item =~ mhal*) set corr_type		= 'mass,h'
				if ($item =~ mcol*) set corr_type		= 'mass,h' 					
				if ($item =~ msta*) set corr_type 		= 'mass,h'
				if ($item =~ *pos) set corr_type 		= 'h'
				if ($item =~ r*) set corr_type 			= 'h'

			else if ($4 == I) then
			#IllustrisTNG300-1 205 Mpc
				if ($item =~ mhal*) set corr_type		= 'h'				
				if ($item =~ *pos) set corr_type 		= 'h'

			else if ($4 == R || $4 == C256 || $4 == C512) then
			#ROCKSTAR on CHOLLA / CHOLLA HYDRO
				if ($item =~ mhal*) set corr_type		= 'h'				
				if ($item =~ mpseud*) set corr_type 	= 'h'
				if ($4 == R) then
					if ($item =~ *pos) set corr_type 		= 'h'
				else
					if ($item =~ *pos) set corr_type 		= 'h,Mpc'
				endif
				if ($item =~ r*) set corr_type 			= 'h,Mpc'
				if ($item =~ _a*) set corr_type 		= 'h,Mpc'
				if ($item =~ _a) set corr_type 			= 'h,Mpc'
				if ($item == xoff) set corr_type 		= 'h,Mpc'
				if ($item =~ *_ang) set corr_type 		= 'h-2'
			else
					
				if ($item =~ msta*) set corr_type 		= 'mass,h'

			endif

		else
		#MD: units according to Multidark data base mass [Msun/h], sfr [Msun/Gyr/h], pos [comMpc/h], mean_age_stars [Gyr], L_* (luminosities) [4.4659e13WHz-1]
		#CM: units according to CMASS paper data base mass [Msun/h], sfr [Msun/yr/h], pos [comMpc/h], mean_age_stars [Gyr], L_* (luminosities) [ergs-1Hz-1]
			if ($4 == A || $4 =~ A* || $4 == C) then
			#Galacticus 				 		 
				if ($item =~ mhal*) set corr_type		= 'h'
				if ($item =~ mbasi*) set corr_type		= 'h'
				if ($item =~ Mz*) set corr_type 		= 'h'	 		
				if ($item =~ mcol*) set corr_type		= 'h'			
				if ($item == mbh) set corr_type		= 'h'			
				if ($item =~ msta*) set corr_type 		= 'h'
				if ($item =~ mho*) set corr_type 		= 'h'
				if ($item == mbasic) set corr_type 		= 'h'

				#Galacticus: wrong name in the columns for z~, those are actually masses of metals!!!
				if ($item =~ zga*) set corr_type		= 'h'
				if ($item =~ zsta*) set corr_type		= 'h'
				if ($item =~ zho*) set corr_type		= 'h'

				if ($item =~ sf*) set corr_type			= 'h'
				if ($item =~ bh_acc*) set corr_type			= 'h'			
				if ($item =~ *pos) set corr_type		= 'COMV,h'
				if ($item =~ r* ) set corr_type			= 'h'
				if ($item =~ ang* ) set corr_type		= 'h'

		
			else if ($4 == B || $4 == E) then
			#SAG
				if ($5 == 'MD') then
					if ($item =~ sf*) set corr_type		= 'Gyr-1'
				endif	
				if ($item =~ *pos) set corr_type		= 'Mpc'
				if ($item == r200c) set corr_type		= 'Mpc'
				if ($item =~ rhalf*) set corr_type		= 'Mpc'

			else if ($4 == F) then
			#SAGE
				if ($item =~ mhal*) set corr_type		= 'mass'
				if ($item =~ Mz*) set corr_type 		= 'mass'	 		
				if ($item =~ mcol*) set corr_type		= 'mass'			
				if ($item == mbh) set corr_type			= 'mass'			
				if ($item =~ msta*) set corr_type 		= 'mass'
				if ($item =~ mho*) set corr_type 		= 'mass'
				if ($5 == 'MD') then
					if ($item =~ sf*) set corr_type		= 'Gyr-1,h'
				else
					if ($item =~ sf*) set corr_type		= 'h'
				endif	

				if ($item =~ mean_age_stars) set corr_type	= '1e3,h-1'
			endif				
		endif
		
		if ($item == halo_nodemass) set corr_type	= 'h'
		if ($item == halo_vmax) set corr_type		= 'h'
		if ($item == halo_rsca) set corr_type		= 'h'
		if ($item == halo_x_pos) set corr_type		= 'h'
		if ($item == halo_y_pos) set corr_type		= 'h'	
		if ($item == halo_z_pos) set corr_type		= 'h'

		if ($item =~ *dmesa) then
			set name_in_plot			= '$d_{\perp}$'
		else if ($item =~ *cmesa) then
			set name_in_plot			= '$c_{\perp}$'		
		else if ($item =~ *cpar) then
			set name_in_plot			= '$c_{\parallel}$'
		else if ($item =~ *r_lt_dmesa) then
			set name_in_plot			= '$r<d_{\perp}$'	
		else if ($item =~ *i_lt_cpar) then
			set name_in_plot			= '$i<c_{\parallel}$'		
		else if ($item == mAB_dA_total_u) then
			set name_in_plot			= '$AB$'
		else if ($item == mAB_total_u) then
			set name_in_plot			= 'mAB_u_without_ext'
		else if ($item == mAB_dA_total_g) then
			set name_in_plot			= '$AB'
		else if ($item == mAB_total_g) then
			set name_in_plot			= 'mAB_g_without_ext'
		else if ($item == mAB_dA_total_r) then
			set name_in_plot			= 'mAB_r_ext'
		else if ($item == mAB_total_r) then
			set name_in_plot			= 'mAB_r_without_ext'
		else if ($item == mAB_dA_total_i) then
			set name_in_plot			= '$AB$'
		else if ($item == mAB_total_i) then
			set name_in_plot			= 'mAB_i_without_ext'
		else if ($item == mAB_dA_total_z) then
			set name_in_plot			= '$z$'
		else if ($item == mAB_total_z) then
			set name_in_plot			= 'mAB_z_without_ext'
		else if ($item =~ *cut_g_r) then
			set name_in_plot			= '$g-r$'
		else if ($item =~ *cut_g_i) then
			set name_in_plot			= '$g-i$'
		else if ($item =~ *cut_g_z) then
			set name_in_plot			= '$g-z$'
		else if ($item =~ *cut_r_i) then
			set name_in_plot			= '$r-i$'
		else if ($item =~ msta*) then
			set name_in_plot			= 'mstar'
		else if ($item == MAB_dA_total_g) then
			set name_in_plot			= 'MAB_g_ext'
		else if ($item == MAB_total_g) then
			set name_in_plot			= 'MAB_g_without_ext'
		else if ($item == MAB_dA_total_r) then
			set name_in_plot			= 'MAB_r_ext'
		else if ($item == MAB_total_r) then
			set name_in_plot			= 'MAB_r_without_ext'
		else if ($item == MAB_dA_total_i) then
			set name_in_plot			= 'MAB_i_ext'
		else if ($item == MAB_total_i) then
			set name_in_plot			= 'MAB_i_without_ext'
		else
			set name_in_plot			= 'default'
		endif
	

		#if ($item == 'mhalo_sat' || $item == mbasic || $item =~ *Index || $item =~ L_SDSS_*_sph* || $item =~ L_SDSS_*_disk* ) then
		#if ( $item == 'mhalo_sat' || $item == mbasic || $item =~ *Index || $item == sfr || $item == mstar || $item == 'Mzstar' || $item == 'Mzgas') then
		if ($item =~ L_SDSS_dA_sph* || $item =~ L_SDSS_dA_disk*) then
			set exclude		= 'yes'
		else
			set exclude		= 'no'
		endif

		echo $item $cut_value_min $cut_value_max $unit $data_type $format "$array[$i]" $name_in_plot $exclude  $output_col_name >> $MYVALUES_CONFIGFILE

		echo $item $corr_type $data_type >> $MY_UNIT_CORRECT_FILE


	endif
	@ i++
end

echo $MYVALUES_CONFIGFILE $MY_UNIT_CORRECT_FILE $MY_NAME_CONV_MAP_FILE



