#!/bin/tcsh

set NAME_SL		= False
set GENERAL_LITTLE_H	= 1
set SIM_LITTLE_H 	= 0.6777
set LITTLE_H_CORR	= 1
set LITTLE_H_1_CORR	= 1
set LITTLE_H_2_CORR	= 1
set GYR_CORR		= 1
set GYR1_CORR		= 1
set e3_CORR		= 1
set MAB_CORR		= 0
set LUM_UNITS_CORR	= 1
set MPC_CORR		= 1
set MASS_CORR		= 1
set KMS1_CORR		= 1
set IMF_CORR		= 1
set CONV_TO_AB_MAG	= False
set COMV_CORR		= False
set APPLY_K_CORR_APP	= False
set APPLY_K_CORR_ABS	= False
set APPLY_Z_BOOST	= False
set COLOR_CODE		=  'k'
set MARKER_CODE		=  'no'
set MARKER_FACECOLOR_CODE		=  'no'
set LINESTYLE_CODE	=  '-'
set COLOR_MAP		= 'Paired'

if ($1 =~ DLB*) then
	set NAME_SL		=  'False'
	set MARKER_CODE		=  '.'
	set COLOR_CODE		=  '005b96'
	set SUBDIRID		=  'False'
	set LITTLE_H		= 0.6777

	echo $1 $item $SAM_redshift $SAM_scale_factor	>> $4

	echo $1 $MARKER_CODE $COLOR_CODE $LINESTYLE_CODE $COLOR_MAP $MARKER_FACECOLOR_CODE	 >> $2

else if ($1 =~ Illustris*) then
	set NAME_SL		=  'False'
	set MARKER_CODE		=  '.'
	set COLOR_CODE		=  '005b96'
	set SUBDIRID		=  'False'
	set LINESTYLE_CODE	=  '-'
	set MARKER_FACECOLOR_CODE		= 'w'
	set COLOR_MAP		= 'Paired'

	set SIM_LITTLE_H 	= 0.6774
	set LITTLE_H_CORR	= 0.6774
	set MPC_CORR		= 1e3

	foreach item ($3)

		if ($item =~ snapshot99*) then
			set SAM_scale_factor = 1.0
			set SAM_redshift = 0.0

		else if ($item =~ snapshot91*) then
			set SAM_scale_factor = 0.909091 
			set SAM_redshift = 0.1

		else if ($item =~ snapshot67*) then
			set SAM_scale_factor =  0.666667
			set SAM_redshift = 0.5

		else if ($item =~ snapshot50*) then
			set SAM_scale_factor =  0.5
			set SAM_redshift = 1.0

		else if ($item =~ snapshot33*) then
			set SAM_scale_factor =  0.333333
			set SAM_redshift = 2.0

		else if ($item =~ snapshot25*) then
			set SAM_scale_factor =   0.249377
			set SAM_redshift = 3.01
		else
			set SAM_scale_factor = 'None'
			set SAM_redshift = 'None'
		endif

		echo $1 $item $SAM_redshift $SAM_scale_factor	>> $4
	end


	echo $1 $MARKER_CODE $COLOR_CODE $LINESTYLE_CODE $COLOR_MAP $MARKER_FACECOLOR_CODE	 >> $2

else if ($1 =~ Galacticus*) then

	if ($1 =~ *run2) then
		set NAME_SL = 'abenson'
	else
		set NAME_SL = 'latest'
	endif
	
	if ($1 =~ *CMASS) then
		set COLOR_CODE		=  'w'
		set COLOR_CODE		=  '#225588'
		set MARKER_CODE		=  'o'
		set LINESTYLE_CODE	=  '--'
		set COLOR_MAP		= 'Paired'
	else if ($1 =~ *CMASS_full) then
		set COLOR_CODE		=  'w'
		set COLOR_CODE		=  'r'
		set MARKER_CODE		=  'no'
		set LINESTYLE_CODE	=  '--'
		set COLOR_MAP		= 'Paired'
	else if ($1 =~ *full) then
		#set COLOR_CODE		=  'w'
		set COLOR_CODE		=  '#225588'
		set MARKER_CODE		=  'no'
		set LINESTYLE_CODE	=  '--'
		set COLOR_MAP		= 'Paired'
	else
		set COLOR_CODE		=  '#225588'
		set COLOR_CODE		=  'w'
		set MARKER_CODE		=  'no'
		set LINESTYLE_CODE	=  '-'
		set COLOR_MAP		= 'Paired'
	endif

	#OII-sfr paper:
	set LINESTYLE_CODE	=  '-'
	set MARKER_FACECOLOR_CODE		= 'w'
	set SUBDIRID		= '103-103'

	if ($5 == 'DS') then
		set LITTLE_H_CORR	= 1
		set GYR_CORR		= 1e-9
		set GYR1_CORR		= 1e9
		set MPC_CORR		= 1e-3
	else
		set LITTLE_H_CORR	= 1.475579165
	endif

	#L = 4.4659e13 [W Hz-1] --> in ergs s-1 = L * 1e7 =  4.4659e20 
	#set LUM_UNITS_CORR	= 2.23919030e-21 this is 1/LUM_UNITS_CORR
	set LUM_UNITS_CORR	= 4.4659e20
	set CONV_TO_AB_MAG	= True
	set COMV_CORR		= True
	set APPLY_K_CORR_ABS	= False
	set APPLY_Z_BOOST	= True

	#Scale Factor: The precision of the redshift changes the scale factor heavily. Therefore enter manually the scale factor from the SAMs:
	
	#Galacticus MDPL2 1 Gpc: 79 z=0, 75 z=0.09, 73 z=0.14, 59 z=0.56, min 1 z=8.0

	foreach item ($3)

		if ($1 =~ *1Gp*) then

			if ($item == region0001) then
				set SAM_scale_factor = 1.0
				set SAM_redshift = 0.00
			#125
			else if ($item == 79) then
				set SAM_scale_factor = 1.0
				set SAM_redshift = 0.0
			#124
			else if ($item == 78) then
				set SAM_scale_factor = 0.980392157
				set SAM_redshift = 0.02
			#123
			else if ($item == 77) then
				set SAM_scale_factor = 0.952380952
				set SAM_redshift = 0.05
			#122
			else if ($item == 76) then
				set SAM_scale_factor = 0.934579439
				set SAM_redshift = 0.07
			#121
			else if ($item == 75) then
				set SAM_scale_factor = 0.9151999999999999
				set SAM_redshift = 0.09
			#120
			else if ($item == 74) then
				set SAM_scale_factor = 0.892857142857
				set SAM_redshift = 0.12
			#119
			else if ($item  == 73) then
				set SAM_scale_factor = 0.8754999999999998
				set SAM_redshift = 0.14
			#118
			else if ($item  == 72) then
				set SAM_scale_factor = 0.854700855
				set SAM_redshift = 0.17
			#117
			else if ($item == 71) then
				set SAM_scale_factor = 0.837520938
				set SAM_redshift = 0.19
			#116
			else if ($item == 70) then
				set SAM_scale_factor = 0.819672131
				set SAM_redshift = 0.22
			#115
			else if ($item == 69) then
				set SAM_scale_factor = 0.8
				set SAM_redshift = 0.25
			#114
			else if ($item == 68) then
				set SAM_scale_factor = 0.78125
				set SAM_redshift = 0.28
			#113	
			else if ($item == 67) then
				set SAM_scale_factor = 0.766871166
				set SAM_redshift = 0.30
			#112	
			else if ($item == 66) then
				set SAM_scale_factor = 0.751879699
				set SAM_redshift = 0.33
			#111	
			else if ($item == 65) then
				set SAM_scale_factor = 0.735294118
				set SAM_redshift = 0.36
			#110	
			else if ($item == 64) then
				set SAM_scale_factor = 0.71942446
				set SAM_redshift = 0.39
			#109
			else if ($item == 63) then
				set SAM_scale_factor = 0.699300699
				set SAM_redshift = 0.43
			#108
			else if ($item == 62) then
				set SAM_scale_factor = 0.684931507
				set SAM_redshift = 0.46
			#107
			else if ($item == 61) then
				set SAM_scale_factor = 0.67114094
				set SAM_redshift = 0.49
			#106
			else if ($item == 60) then
				set SAM_scale_factor = 0.657894737
				set SAM_redshift = 0.52
			#105
			else if ($item  == 59) then 
				set SAM_scale_factor = 0.6421
				set SAM_redshift = 0.56
			#104
			else if ($item  == 58) then 
				set SAM_scale_factor = 0.628140704
				set SAM_redshift = 0.59
			#103
			else if ($item  == 57) then 
				set SAM_scale_factor = 0.614250614
				set SAM_redshift = 0.63
			#102	
			else if ($item == 56) then
				set SAM_scale_factor = 0.600961538
				set SAM_redshift = 0.66
			#101	
			else if ($item == 55) then
				set SAM_scale_factor = 0.587544066
				set SAM_redshift = 0.7
			#100
			else if ($item  == 54) then 
				set SAM_scale_factor = 0.5747
				set SAM_redshift = 0.74
			#99
			else if ($item == 53) then
				set SAM_scale_factor = 0.562113547
				set SAM_redshift = 0.78
			#98
			else if ($item == 52) then
				set SAM_scale_factor = 0.549752611
				set SAM_redshift = 0.82
			#97
			else if ($item == 51) then
				set SAM_scale_factor = 0.537923615
				set SAM_redshift = 0.86
			#96
			else if ($item  == 50) then 
				set SAM_scale_factor = 0.526
				set SAM_redshift = 0.9
			#95
			else if ($item  == 49) then 
				set SAM_scale_factor = 0.514403292
				set SAM_redshift = 0.94
			#94
			else if ($item == 48) then
				set SAM_scale_factor = 0.503271263
				set SAM_redshift = 0.99
			#93
			else if ($item == 47) then
				set SAM_scale_factor = 0.492125984
				set SAM_redshift = 1.03
			#92
			else if ($item  == 46) then 
				set SAM_scale_factor = 0.481463649
				set SAM_redshift = 1.08
			#91
			else if ($item  == 45) then 
				set SAM_scale_factor = 0.470809793
				set SAM_redshift = 1.12
			#90
			else if ($item  == 44) then 
				set SAM_scale_factor = 0.460405157
				set SAM_redshift = 1.17
			#89
			else if ($item  == 43) then 
				set SAM_scale_factor = 0.4505
				set SAM_redshift = 1.22
			#88
			else if ($item  == 42) then 
				set SAM_scale_factor = 0.440528634
				set SAM_redshift = 1.27
			#87
			else if ($item  == 41) then 
				set SAM_scale_factor = 0.430848772
				set SAM_redshift = 1.32
			#86
			else if ($item  == 40) then 
				set SAM_scale_factor = 0.42158516
				set SAM_redshift = 1.372

			#85
			else if ($item  == 39) then 
				set SAM_scale_factor = 0.411522634
				set SAM_redshift = 1.43
			#84
			else if ($item  == 38) then 
				set SAM_scale_factor = 0.403225806
				set SAM_redshift = 1.48
			#83
			else if ($item  == 37) then 
				set SAM_scale_factor = 0.393700787
				set SAM_redshift = 1.54
			#82
			else if ($item  == 36) then 
				set SAM_scale_factor = 0.386100386
				set SAM_redshift = 1.59
			#81
			else if ($item  == 35) then 
				set SAM_scale_factor = 0.377358491
				set SAM_redshift = 1.65
			#80
			else if ($item  == 34) then 
				set SAM_scale_factor = 0.36900369
				set SAM_redshift = 1.71
			#79
			else if ($item  == 33) then 
				set SAM_scale_factor = 0.36101083
				set SAM_redshift = 1.77
			#78
			else if ($item  == 32) then
				set SAM_scale_factor = 0.35335689
				set SAM_redshift = 1.83
			#77
			else if ($item == 31) then
				set SAM_scale_factor = 0.344827586
				set SAM_redshift = 1.9
			#76
			else if ($item == 30) then
				set SAM_scale_factor = 0.337837838
				set SAM_redshift = 1.96
			#75
			else if ($item  == 29) then 
				set SAM_scale_factor = 0.3303
				set SAM_redshift = 2.03
			#74
			else if ($item  == 28) then 
				set SAM_scale_factor = 0.322580645
				set SAM_redshift = 2.1
			#73
			else if ($item  == 27) then 
				set SAM_scale_factor = 0.315457413
				set SAM_redshift = 2.17
			#72
			else if ($item  == 26) then 
				set SAM_scale_factor = 0.308641975
				set SAM_redshift = 2.24
			#71
			else if ($item  == 25) then 
				set SAM_scale_factor = 0.302114804
				set SAM_redshift = 2.31
			#70
			else if ($item  == 24) then 
				set SAM_scale_factor = 0.2957
				set SAM_redshift = 2.38
			#69
			else if ($item  == 23) then 
				set SAM_scale_factor = 0.289017341
				set SAM_redshift = 2.46
			#68
			else if ($item  == 22) then 
				set SAM_scale_factor = 0.282485876
				set SAM_redshift = 2.54
			#67
			else if ($item  == 21) then 
				set SAM_scale_factor = 0.27700831
				set SAM_redshift = 2.61
			#66
			else if ($item  == 20) then 
				set SAM_scale_factor = 0.27027027
				set SAM_redshift = 2.7
			#65
			else if ($item  == 19) then 
				set SAM_scale_factor = 0.264550265
				set SAM_redshift = 2.78
			#64
			else if ($item  == 18) then 
				set SAM_scale_factor = 0.259067358
				set SAM_redshift = 2.86
			#63
			else if ($item  == 17) then 
				set SAM_scale_factor = 0.253164557
				set SAM_redshift = 2.95
			#62
			else if ($item  == 16) then 
				set SAM_scale_factor = 0.247524752
				set SAM_redshift = 3.04
			#61
			else if ($item  == 15) then 
				set SAM_scale_factor = 0.242130751
				set SAM_redshift = 3.13
			#60
			else if ($item  == 14) then 
				set SAM_scale_factor = 0.236966825
				set SAM_redshift = 3.22
			#59
			else if ($item  == 13) then 
				set SAM_scale_factor = 0.232018561
				set SAM_redshift = 3.31
			#58
			else if ($item  == 12) then 
				set SAM_scale_factor = 0.22675737
				set SAM_redshift = 3.41
			#57
			else if ($item  == 11) then 
				set SAM_scale_factor = 0.22172949
				set SAM_redshift = 3.51
			#56
			else if ($item  == 10) then 
				set SAM_scale_factor = 0.21691974
				set SAM_redshift = 3.61
			#55
			else if ($item  == 9) then 
				set SAM_scale_factor = 0.211864407
				set SAM_redshift = 3.72
			#54
			else if ($item  == 8) then 
				set SAM_scale_factor = 0.20746888
				set SAM_redshift = 3.82
			#53
			else if ($item  == 7) then 
				set SAM_scale_factor = 0.202839757
				set SAM_redshift = 3.93
			#52
			else if ($item  == 6) then 
				set SAM_scale_factor = 0.198412698
				set SAM_redshift = 4.04
			#51
			else if ($item  == 5) then 
				set SAM_scale_factor = 0.194174757
				set SAM_redshift = 4.15
			#50
			else if ($item  == 4) then 
				set SAM_scale_factor = 0.149253731
				set SAM_redshift = 5.7
			#49
			else if ($item  == 3) then 
				set SAM_scale_factor = 0.142857143
				set SAM_redshift = 6.0
			#48
			else if ($item  == 2) then 
				set SAM_scale_factor = 0.125
				set SAM_redshift = 7.0
			#47
			else if ($item  == 1) then 
				set SAM_scale_factor = 0.111111111
				set SAM_redshift = 8.0

			else
				set SAM_scale_factor = 'None'
				set SAM_redshift = 'None'
			endif
		
		else if ($1 =~ *400*) then

			#116
			if ($item == 96) then
				set SAM_scale_factor = 1.0
				set SAM_redshift = 0.00

			#112
			else if ($item == 92) then
				set SAM_scale_factor = 0.9074
				set SAM_redshift = 0.1

			#106		
			else if ($item  == 85) then
				set SAM_scale_factor = 0.8754999999999998
				set SAM_redshift = 0.14

			#71
			else if ($item  == 51) then 
				set SAM_scale_factor = 0.6704
				set SAM_redshift = 0.49

			#67
			else if ($item  == 47) then 
				set SAM_scale_factor = 0.6464
				set SAM_redshift = 0.55

			#48
			else if ($item  == 28) then 
				set SAM_scale_factor = 0.5
				set SAM_redshift = 1.0

			#38
			else if ($item  == 18) then 
				set SAM_scale_factor =  0.331126
				set SAM_redshift = 2.02

			#37
			else if ($item  == 17) then 
				set SAM_scale_factor = 0.3180
				set SAM_redshift = 2.14

			#31
			else if ($item  == 11) then 
				set SAM_scale_factor = 0.248139 
				set SAM_redshift = 3.03

			else if ($item  == 6) then 
				set SAM_scale_factor = 0.202020202
				set SAM_redshift = 3.95

			else if ($item  == 4) then 
				set SAM_scale_factor = 0.149253731
				set SAM_redshift = 5.7

			endif

		else
			if ($item  == 1) then
				set SAM_scale_factor = 0.1111111111111111
				set SAM_redshift = 8.0
			else
				set SAM_scale_factor = 'None'
				set SAM_redshift = 'None'
			endif
		endif


		echo $1 $item $SAM_redshift $SAM_scale_factor	>> $4
	end


	echo $1 $MARKER_CODE $COLOR_CODE $LINESTYLE_CODE $COLOR_MAP $MARKER_FACECOLOR_CODE	 >> $2

else if ($1 =~ Galform*) then
	if ($1 =~ *gp15*) then
		set MARKER_CODE		=  '<'
		set COLOR_CODE		=  'c1073f'
		echo $1 $MARKER_CODE $COLOR_CODE			>> $2
	else if ($1 =~ *gp14*) then
		set MARKER_CODE		=  '^'
		set COLOR_CODE		=  '660066'
		echo $1 $MARKER_CODE $COLOR_CODE			>> $2
	else if ($1 =~ *kf*) then
		set MARKER_CODE		=  'v'
		set COLOR_CODE		=  '72021a'
		echo $1 $MARKER_CODE $COLOR_CODE			>> $2
	else if ($1 =~ *kb*) then
		set MARKER_CODE		=  '>'
		set COLOR_CODE		=  '828fd0'
		echo $1 $MARKER_CODE $COLOR_CODE			>> $2
	else
		set MARKER_CODE		=  '+'
		set COLOR_CODE		=  'f47835'
		echo $1 $MARKER_CODE $COLOR_CODE			>> $2
	endif
	set NAME_SL		=  'False'
	set SUBDIRID		=  'False'
	set LITTLE_H		= 0.6777


else if ($1 =~ GalICS*) then
	set NAME_SL		=  'False'
	set MARKER_CODE		=  'x'
	set COLOR_CODE		=  'ffe200'
	set SUBDIRID		=  'False'
	set LITTLE_H		= 0.6777

	echo $1 $MARKER_CODE $COLOR_CODE			>> $2


else if ($1 =~ LGALAXIES*) then
	set NAME_SL		=  'False'
	set MARKER_CODE		=  'o'
	set SUBDIRID		=  'False'
	set LITTLE_H		= 0.6777
	set SIM_LITTLE_H 	= 0.6777
	set COLOR_MAP		= 'Paired'
	set COLOR_CODE		=  'CC6677'
	set LINESTYLE_CODE	=  '-'
	set MARKER_FACECOLOR_CODE = 'w'
	set SUBDIRID		= '-'


	#Corrections:
	set LITTLE_H_CORR	= 0.6777
	set MASS_CORR		= 'power10'

	set APPLY_K_CORR_ABS	= False
	foreach item ($3)

		if ($item =~ *45) then
			set SAM_scale_factor = 0.636314
			set SAM_redshift = 0.56

		else if ($item =~ *MergerTree*) then
			set SAM_scale_factor = 1.0
			set SAM_redshift = 0.0
		else
			set SAM_scale_factor = 1.0
			set SAM_redshift = 0.0
		endif

		echo $1 $item $SAM_redshift $SAM_scale_factor	>> $4

	end

	echo $1 $MARKER_CODE $COLOR_CODE $LINESTYLE_CODE $COLOR_MAP $MARKER_FACECOLOR_CODE	 >> $2

else if ($1 =~ MICE*) then
	set NAME_SL		=  'False'
	set MARKER_CODE		=  'h'
	set COLOR_CODE		=  'f47835'
	set SUBDIRID		=  'False'
	set LITTLE_H		= 0.6777

	echo $1 $MARKER_CODE $COLOR_CODE			>> $2


else if ($1 =~ MORGANA*) then
	set NAME_SL		=  'False'
	set MARKER_CODE		=  '8'
	set COLOR_CODE		=  'eac117'
	set SUBDIRID		=  'False'
	set LITTLE_H		= 0.6777

	echo $1 $MARKER_CODE $COLOR_CODE			>> $2

else if ($1 == ROCKSTAR_50Mpc || $1 =~ Cholla* ) then
	#Rockstar run on Cholla 50 Mpc/h box
	set NAME_SL		=  'False'
	set MARKER_CODE		=  '.'
	set COLOR_CODE		=  '005b96'
	set SUBDIRID		=  'False'
	set LITTLE_H		= 0.6766
	set SIM_LITTLE_H 	= 0.6766

	if ($5 == 'DS') then
		set LITTLE_H_CORR	= 0.6766
		set LITTLE_H_2_CORR	= 0.45778756
		set MPC_CORR		= 1e3
	endif

	foreach item ($3)
		if ($item  == 299) then 
			set SAM_scale_factor = 0.166667
			set SAM_redshift = 5.0
		else if ($item  == 73) then 
			set SAM_scale_factor = 1
			set SAM_redshift = 0.0
		else if ($item  =~ 73*) then 
			set SAM_scale_factor = 1
			set SAM_redshift = 0.0
		else if ($item  == 72) then 
			set SAM_scale_factor = 0.98850
			set SAM_redshift = 0.01164
		else if ($item  =~ 72*) then 
			set SAM_scale_factor = 0.98850
			set SAM_redshift = 0.01164

		else if ($item  == 51) then 
			set SAM_scale_factor = 0.553423
			set SAM_redshift = 0.81
		else if ($item  =~ 51*) then 
			set SAM_scale_factor = 0.553423
			set SAM_redshift = 0.81
		else if ($item  == 50) then 
			set SAM_scale_factor = 0.535851
			set SAM_redshift = 0.87
		else if ($item  =~ 50*) then 
			set SAM_scale_factor = 0.535851
			set SAM_redshift = 0.87
		else if ($item  == 49) then 
			set SAM_scale_factor = 0.518436
			set SAM_redshift = 0.93
		else if ($item  =~ 49*) then 
			set SAM_scale_factor = 0.518436
			set SAM_redshift = 0.93
		else if ($item  == 48) then 
			set SAM_scale_factor = 0.501043
			set SAM_redshift = 1.0
		else if ($item  =~ 48*) then 
			set SAM_scale_factor = 0.501043
			set SAM_redshift = 1.0
		else if ($item  == 47) then 
			set SAM_scale_factor = 0.483629
			set SAM_redshift = 1.07
		else if ($item  =~ 47*) then 
			set SAM_scale_factor = 0.483629
			set SAM_redshift = 1.07
		else if ($item  == 46) then 
			set SAM_scale_factor = 0.466233
			set SAM_redshift = 1.14
		else if ($item  =~ 46*) then 
			set SAM_scale_factor = 0.466233
			set SAM_redshift = 1.14
		else if ($item  == 45) then 
			set SAM_scale_factor = 0.448888
			set SAM_redshift = 1.23
		else if ($item  =~ 45*) then 
			set SAM_scale_factor = 0.448888
			set SAM_redshift = 1.23
		else if ($item  == 44) then 
			set SAM_scale_factor = 0.431490
			set SAM_redshift = 1.32
		else if ($item  =~ 44*) then 
			set SAM_scale_factor = 0.431490
			set SAM_redshift = 1.32
		else if ($item  == 43) then 
			set SAM_scale_factor = 0.413924
			set SAM_redshift = 1.42
		else if ($item  =~ 43*) then 
			set SAM_scale_factor = 0.413924
			set SAM_redshift = 1.42
		else if ($item  == 42) then 
			set SAM_scale_factor = 0.396173
			set SAM_redshift = 1.52
		else if ($item  =~ 42*) then 
			set SAM_scale_factor = 0.396173
			set SAM_redshift = 1.52
		else if ($item  == 41) then 
			set SAM_scale_factor = 0.378935
			set SAM_redshift = 1.64
		else if ($item  =~ 41*) then 
			set SAM_scale_factor = 0.378935
			set SAM_redshift = 1.64
		else if ($item  == 40) then 
			set SAM_scale_factor = 0.362396
			set SAM_redshift = 1.76
		else if ($item  =~ 40*) then 
			set SAM_scale_factor = 0.362396
			set SAM_redshift = 1.76
		else if ($item  == 39) then 
			set SAM_scale_factor = 0.346561
			set SAM_redshift = 1.89
		else if ($item  =~ 39*) then 
			set SAM_scale_factor = 0.346561
			set SAM_redshift = 1.89
		else if ($item  == 38) then 
			set SAM_scale_factor = 0.330164
			set SAM_redshift = 2.03
		else if ($item  =~ 38*) then 
			set SAM_scale_factor = 0.330164
			set SAM_redshift = 2.03
		else if ($item  == 37) then 
			set SAM_scale_factor = 0.313620
			set SAM_redshift = 2.19
		else if ($item  =~ 37*) then 
			set SAM_scale_factor = 0.313620
			set SAM_redshift = 2.19
		else if ($item  == 36) then 
			set SAM_scale_factor = 0.297190
			set SAM_redshift = 2.36
		else if ($item  =~ 36*) then 
			set SAM_scale_factor = 0.297190
			set SAM_redshift = 2.36
		else if ($item  == 35) then 
			set SAM_scale_factor = 0.280856
			set SAM_redshift = 2.56
		else if ($item  =~ 35*) then 
			set SAM_scale_factor = 0.280856
			set SAM_redshift = 2.56
		else if ($item  == 34) then 
			set SAM_scale_factor = 0.264317
			set SAM_redshift = 2.78
		else if ($item  =~ 34*) then 
			set SAM_scale_factor = 0.264317
			set SAM_redshift = 2.78
		else if ($item  == 33) then 
			set SAM_scale_factor = 0.246770
			set SAM_redshift = 3.05
		else if ($item  =~ 33*) then 
			set SAM_scale_factor = 0.246770
			set SAM_redshift = 3.05
		else if ($item  == 32) then 
			set SAM_scale_factor = 0.228770
			set SAM_redshift = 3.37
		else if ($item  =~ 32*) then 
			set SAM_scale_factor = 0.228770
			set SAM_redshift = 3.37
		else if ($item  == 31) then 
			set SAM_scale_factor = 0.211173
			set SAM_redshift = 3.74
		else if ($item  =~ 31*) then 
			set SAM_scale_factor = 0.211173
			set SAM_redshift = 3.74
		else if ($item  == 30) then 
			set SAM_scale_factor = 0.193994
			set SAM_redshift = 4.15
		else if ($item  =~ 30*) then 
			set SAM_scale_factor = 0.193994
			set SAM_redshift = 4.15
		else if ($item  == 29) then 
			set SAM_scale_factor = 0.176213
			set SAM_redshift = 4.67
		else if ($item  =~ 29*) then 
			set SAM_scale_factor = 0.176213
			set SAM_redshift = 4.67
		else if ($item  == 28) then 
			set SAM_scale_factor = 0.159674
			set SAM_redshift = 5.26
		else if ($item  =~ 28*) then 
			set SAM_scale_factor = 0.159674
			set SAM_redshift = 5.26
		else if ($item  == 27) then 
			set SAM_scale_factor = 0.144680
			set SAM_redshift = 5.91
		else if ($item  =~ 27*) then 
			set SAM_scale_factor = 0.144680
			set SAM_redshift = 5.91
		else if ($item  == 26) then 
			set SAM_scale_factor = 0.131083
			set SAM_redshift = 6.63
		else if ($item  =~ 26*) then 
			set SAM_scale_factor = 0.131083
			set SAM_redshift = 6.63
		else if ($item  == 25) then 
			set SAM_scale_factor = 0.118752
			set SAM_redshift = 7.42
		else if ($item  =~ 25*) then 
			set SAM_scale_factor = 0.118752
			set SAM_redshift = 7.42
		else if ($item  == 24) then 
			set SAM_scale_factor = 0.107571
			set SAM_redshift = 8.27
		else if ($item  =~ 24*) then 
			set SAM_scale_factor = 0.107571
			set SAM_redshift = 8.27
		else if ($item  == 23) then 
			set SAM_scale_factor = 0.097434
			set SAM_redshift = 9.26
		else if ($item  =~ 23*) then 
			set SAM_scale_factor = 0.097434
			set SAM_redshift = 9.26
		else if ($item  == 22) then 
			set SAM_scale_factor = 0.088247
			set SAM_redshift = 10.33
		else if ($item  =~ 22*) then 
			set SAM_scale_factor = 0.088247
			set SAM_redshift = 10.33
		else if ($item  == 21) then 
			set SAM_scale_factor = 0.079918
			set SAM_redshift = 11.51
		else if ($item  =~ 21*) then 
			set SAM_scale_factor = 0.98850
			set SAM_redshift = 11.51
		else
			set SAM_scale_factor = 'None'
			set SAM_redshift = 'None'
		endif

		echo $1 $item $SAM_redshift $SAM_scale_factor	>> $4
	end

	echo $1 $MARKER_CODE $COLOR_CODE $LINESTYLE_CODE $COLOR_MAP $MARKER_FACECOLOR_CODE	 >> $2

else if ($1 =~ SantaCruz*) then
	set NAME_SL		=  'False'
	set MARKER_CODE		=  'H'
	set COLOR_CODE		=  '#00d27f'
	set SUBDIRID		=  'False'
	set LITTLE_H		= 0.6777

	echo $1 $MARKER_CODE $COLOR_CODE			>> $2


else if ($1 =~ SAG_*) then

	if ($1 =~ *125) then
		set COLOR_CODE		=  '003964'
		set COLOR_CODE		=  'r'
		set MARKER_CODE		=  'no'
		set LINESTYLE_CODE	=  '-'
		set COLOR_MAP		= 'Paired'
	else if ($1 =~ *v2) then
		set NAME_SL		= 'SAG-7.86-fullmodel'
		set COLOR_CODE		=  'w'
		set MARKER_CODE		=  'o'
		#OII-sfr paper:
		set LINESTYLE_CODE	=  '-'
		set COLOR_MAP		= 'Paired'
	else
		set NAME_SL		= 'SAG-7.128'
		set COLOR_MAP		= 'Paired'
		set COLOR_CODE		=  '005900'
		set COLOR_CODE		=  'w'
		set MARKER_CODE		=  'no'
		set LINESTYLE_CODE	=  '-'
	endif


	set MARKER_FACECOLOR_CODE	= 'w'
	set SUBDIRID		= '1-128'
	#DO NOT USE THAT ANYMORE BECAUSE THE CATALOGS ARE NOW HIGHLY DIVERSE!!!!
	set MPC_CORR		= 1e3
	set APPLY_K_CORR_APP	= False

	if ($5 == 'DS') then
		set LITTLE_H_CORR	= 0.6777
		set LITTLE_H_2_CORR	= 0.45927729
	else
		set LITTLE_H_CORR	= 1
		set GYR_CORR		= 1e-9
	endif

	#echo 564

	foreach item ($3)
		#echo 567
		if ($1 =~ *1Gp*) then
			#125
			if ($item == 125) then
				set SAM_scale_factor = 1.0
				set SAM_redshift = 0.0
				#echo 573
			#124
			else if ($item == 124) then
				set SAM_scale_factor = 0.980392157
				set SAM_redshift = 0.02
			#123
			else if ($item == 123) then
				set SAM_scale_factor = 0.952380952
				set SAM_redshift = 0.05
			#122
			else if ($item == 122) then
				set SAM_scale_factor = 0.934579439
				set SAM_redshift = 0.07
			#121
			else if ($item == 121) then
				set SAM_scale_factor = 0.9151999999999999
				set SAM_redshift = 0.09
			#120
			else if ($item == 120) then
				set SAM_scale_factor = 0.892857142857
				set SAM_redshift = 0.12
			#119
			else if ($item  == 119) then
				set SAM_scale_factor = 0.8754999999999998
				set SAM_redshift = 0.14
			#118
			else if ($item  == 118) then
				set SAM_scale_factor = 0.854700855
				set SAM_redshift = 0.17
			#117
			else if ($item == 117) then
				set SAM_scale_factor = 0.837520938
				set SAM_redshift = 0.19
			#116
			else if ($item == 116) then
				set SAM_scale_factor = 0.819672131
				set SAM_redshift = 0.22
			#115
			else if ($item == 115) then
				set SAM_scale_factor = 0.8
				set SAM_redshift = 0.25
			#114
			else if ($item == 114) then
				set SAM_scale_factor = 0.78125
				set SAM_redshift = 0.28
			#113	
			else if ($item == 113) then
				set SAM_scale_factor = 0.766871166
				set SAM_redshift = 0.30
			#112	
			else if ($item == 112) then
				set SAM_scale_factor = 0.751879699
				set SAM_redshift = 0.33
			#111	
			else if ($item == 111) then
				set SAM_scale_factor = 0.735294118
				set SAM_redshift = 0.36
			#110	
			else if ($item == 110) then
				set SAM_scale_factor = 0.71942446
				set SAM_redshift = 0.39
			#109
			else if ($item == 109) then
				set SAM_scale_factor = 0.699300699
				set SAM_redshift = 0.43
			#108
			else if ($item == 108) then
				set SAM_scale_factor = 0.684931507
				set SAM_redshift = 0.46
			#107
			else if ($item == 107) then
				set SAM_scale_factor = 0.67114094
				set SAM_redshift = 0.49
			#106
			else if ($item == 106) then
				set SAM_scale_factor = 0.657894737
				set SAM_redshift = 0.52
			#105
			else if ($item  == 105) then 
				set SAM_scale_factor = 0.6421
				set SAM_redshift = 0.56
			#104
			else if ($item  == 104) then 
				set SAM_scale_factor = 0.628140704
				set SAM_redshift = 0.59
			#103
			else if ($item  == 103) then 
				set SAM_scale_factor = 0.614250614
				set SAM_redshift = 0.63
			#102	
			else if ($item == 102) then
				set SAM_scale_factor = 0.600961538
				set SAM_redshift = 0.66
			#101	
			else if ($item == 101) then
				set SAM_scale_factor = 0.587544066
				set SAM_redshift = 0.7
			#100
			else if ($item  == 100) then 
				set SAM_scale_factor = 0.5747
				set SAM_redshift = 0.74
			#99
			else if ($item == 99) then
				set SAM_scale_factor = 0.562113547
				set SAM_redshift = 0.78
			#98
			else if ($item == 98) then
				set SAM_scale_factor = 0.549752611
				set SAM_redshift = 0.82
			#97
			else if ($item == 97) then
				set SAM_scale_factor = 0.537923615
				set SAM_redshift = 0.86
			#96
			else if ($item  == 96) then 
				set SAM_scale_factor = 0.526
				set SAM_redshift = 0.9
			#95
			else if ($item  == 95) then 
				set SAM_scale_factor = 0.514403292
				set SAM_redshift = 0.94
			#94
			else if ($item == 94) then
				set SAM_scale_factor = 0.503271263
				set SAM_redshift = 0.99
			#93
			else if ($item == 93) then
				set SAM_scale_factor = 0.492125984
				set SAM_redshift = 1.03
			#92
			else if ($item  == 92) then 
				set SAM_scale_factor = 0.481463649
				set SAM_redshift = 1.08
			#91
			else if ($item  == 91) then 
				set SAM_scale_factor = 0.470809793
				set SAM_redshift = 1.12
			#90
			else if ($item  == 90) then 
				set SAM_scale_factor = 0.460405157
				set SAM_redshift = 1.17
			#89
			else if ($item  == 89) then 
				set SAM_scale_factor = 0.4505
				set SAM_redshift = 1.22
			#88
			else if ($item  == 88) then 
				set SAM_scale_factor = 0.440528634
				set SAM_redshift = 1.27
			#87
			else if ($item  == 87) then 
				set SAM_scale_factor = 0.430848772
				set SAM_redshift = 1.32
			#86
			else if ($item  == 86) then 
				set SAM_scale_factor = 0.42158516
				set SAM_redshift = 1.37
			#85
			else if ($item  == 85) then 
				set SAM_scale_factor = 0.411522634
				set SAM_redshift = 1.43
			#84
			else if ($item  == 84) then 
				set SAM_scale_factor = 0.403225806
				set SAM_redshift = 1.48
			#83
			else if ($item  == 83) then 
				set SAM_scale_factor = 0.393700787
				set SAM_redshift = 1.54
			#82
			else if ($item  == 82) then 
				set SAM_scale_factor = 0.386100386
				set SAM_redshift = 1.59
			#81
			else if ($item  == 81) then 
				set SAM_scale_factor = 0.377358491
				set SAM_redshift = 1.65
			#80
			else if ($item  == 80) then 
				set SAM_scale_factor = 0.36900369
				set SAM_redshift = 1.71
			#79
			else if ($item  == 79) then 
				set SAM_scale_factor = 0.36101083
				set SAM_redshift = 1.77
			#78
			else if ($item  == 78) then
				set SAM_scale_factor = 0.35335689
				set SAM_redshift = 1.83
			#77
			else if ($item == 77) then
				set SAM_scale_factor = 0.344827586
				set SAM_redshift = 1.9
			#76
			else if ($item == 76) then
				set SAM_scale_factor = 0.337837838
				set SAM_redshift = 1.96
			#75
			else if ($item  == 75) then 
				set SAM_scale_factor = 0.3303
				set SAM_redshift = 2.03
			#74
			else if ($item  == 74) then 
				set SAM_scale_factor = 0.322580645
				set SAM_redshift = 2.1
			#73
			else if ($item  == 73) then 
				set SAM_scale_factor = 0.315457413
				set SAM_redshift = 2.17
			#72
			else if ($item  == 72) then 
				set SAM_scale_factor = 0.308641975
				set SAM_redshift = 2.24
			#71
			else if ($item  == 71) then 
				set SAM_scale_factor = 0.302114804
				set SAM_redshift = 2.31
			#70
			else if ($item  == 70) then 
				set SAM_scale_factor = 0.2957
				set SAM_redshift = 2.38
			#69
			else if ($item  == 69) then 
				set SAM_scale_factor = 0.289017341
				set SAM_redshift = 2.46
			#68
			else if ($item  == 68) then 
				set SAM_scale_factor = 0.282485876
				set SAM_redshift = 2.54
			#67
			else if ($item  == 67) then 
				set SAM_scale_factor = 0.27700831
				set SAM_redshift = 2.61
			#66
			else if ($item  == 66) then 
				set SAM_scale_factor = 0.27027027
				set SAM_redshift = 2.7
			#65
			else if ($item  == 65) then 
				set SAM_scale_factor = 0.264550265
				set SAM_redshift = 2.78
			#64
			else if ($item  == 64) then 
				set SAM_scale_factor = 0.259067358
				set SAM_redshift = 2.86
			#63
			else if ($item  == 63) then 
				set SAM_scale_factor = 0.253164557
				set SAM_redshift = 2.95
			#62
			else if ($item  == 62) then 
				set SAM_scale_factor = 0.247524752
				set SAM_redshift = 3.04
			#61
			else if ($item  == 61) then 
				set SAM_scale_factor = 0.242130751
				set SAM_redshift = 3.13
			#60
			else if ($item  == 60) then 
				set SAM_scale_factor = 0.236966825
				set SAM_redshift = 3.22
			#59
			else if ($item  == 59) then 
				set SAM_scale_factor = 0.232018561
				set SAM_redshift = 3.31
			#58
			else if ($item  == 58) then 
				set SAM_scale_factor = 0.22675737
				set SAM_redshift = 3.41
			#57
			else if ($item  == 57) then 
				set SAM_scale_factor = 0.22172949
				set SAM_redshift = 3.51
			#56
			else if ($item  == 56) then 
				set SAM_scale_factor = 0.21691974
				set SAM_redshift = 3.61
			#55
			else if ($item  == 55) then 
				set SAM_scale_factor = 0.211864407
				set SAM_redshift = 3.72
			#54
			else if ($item  == 54) then 
				set SAM_scale_factor = 0.20746888
				set SAM_redshift = 3.82
			#53
			else if ($item  == 53) then 
				set SAM_scale_factor = 0.202839757
				set SAM_redshift = 3.93
			#52
			else if ($item  == 52) then 
				set SAM_scale_factor = 0.198412698
				set SAM_redshift = 4.04
			#51
			else if ($item  == 51) then 
				set SAM_scale_factor = 0.194174757
				set SAM_redshift = 4.15
			#50
			else if ($item  == 50) then 
				set SAM_scale_factor = 0.189900
				set SAM_redshift = 4.27
			#49
			else if ($item  == 49) then 
				set SAM_scale_factor = 0.185700
				set SAM_redshift = 4.39
			#48
			else if ($item  == 48) then 
				set SAM_scale_factor = 0.181600
				set SAM_redshift = 4.51
			#47
			else if ($item  == 47) then 
				set SAM_scale_factor = 0.177700
				set SAM_redshift = 4.63
			#46
			else if ($item == 46) then
				set SAM_scale_factor = 0.173800
				set SAM_redshift = 4.75
			#45
			else if ($item == 45) then
				set SAM_scale_factor = 0.170000
				set SAM_redshift = 4.88
			#44
			else if ($item == 44) then
				set SAM_scale_factor = 0.166200
				set SAM_redshift = 5.02
			#43
			else if ($item  == 43) then 
				set SAM_scale_factor = 0.162600
				set SAM_redshift = 5.15
			#42
			else if ($item == 42) then
				set SAM_scale_factor = 0.159000
				set SAM_redshift = 5.29
			#41
			else if ($item == 41) then
				set SAM_scale_factor = 0.155500
				set SAM_redshift = 5.43
			#40
			else if ($item == 40) then
				set SAM_scale_factor = 0.152100
				set SAM_redshift = 5.57
			#39
			else if ($item  == 39) then 
				set SAM_scale_factor = 0.148800
				set SAM_redshift = 5.72
			#38
			else if ($item  == 38) then 
				set SAM_scale_factor = 0.145500
				set SAM_redshift = 5.87
			#37
			else if ($item  == 37) then 
				set SAM_scale_factor = 0.142400
				set SAM_redshift = 6.02
			#36
			else if ($item  == 36) then 
				set SAM_scale_factor = 0.139200
				set SAM_redshift = 6.18
			#35
			else if ($item == 35) then
				set SAM_scale_factor = 0.136200
				set SAM_redshift = 6.34
			#34
			else if ($item == 34) then
				set SAM_scale_factor = 0.133200
				set SAM_redshift = 6.51
			#33
			else if ($item == 33) then
				set SAM_scale_factor = 0.130300
				set SAM_redshift = 6.67
			#32
			else if ($item  == 32) then 
				set SAM_scale_factor = 0.127400
				set SAM_redshift = 6.85
			#31
			else if ($item == 31) then
				set SAM_scale_factor = 0.124600
				set SAM_redshift = 7.03
			#30
			else if ($item == 30) then
				set SAM_scale_factor = 0.121900
				set SAM_redshift = 7.2
			#29
			else if ($item == 29) then
				set SAM_scale_factor = 0.119200
				set SAM_redshift = 7.39
			#28
			else if ($item == 28) then
				set SAM_scale_factor = 0.116550117
				set SAM_redshift = 7.58
			#27
			else if ($item  == 27) then 
				set SAM_scale_factor = 0.114155251
				set SAM_redshift = 7.76
			#26
			else if ($item  == 26) then 
				set SAM_scale_factor = 0.111607143
				set SAM_redshift = 7.96
			#25
			else if ($item  == 25) then 
				set SAM_scale_factor = 0.109051254
				set SAM_redshift = 8.17
			#24
			else if ($item  == 24) then 
				set SAM_scale_factor = 0.106723586
				set SAM_redshift = 8.37
			#23
			else if ($item == 23) then
				set SAM_scale_factor = 0.104384134
				set SAM_redshift = 8.58
			#22
			else if ($item == 22) then
				set SAM_scale_factor = 0.102040816
				set SAM_redshift = 8.8
			#21
			else if ($item == 21) then
				set SAM_scale_factor = 0.0999001
				set SAM_redshift = 9.01
			#20
			else if ($item  == 20) then 
				set SAM_scale_factor = 0.09765625
				set SAM_redshift = 9.24
			#19
			else if ($item == 019) then
				set SAM_scale_factor = 0.095510984
				set SAM_redshift = 9.47
			#18
			else if ($item == 018) then
				set SAM_scale_factor = 0.093370682
				set SAM_redshift = 9.71
			#17
			else if ($item == 017) then
				set SAM_scale_factor = 0.091407678
				set SAM_redshift = 9.94
			#16
			else if ($item  == 016) then 
				set SAM_scale_factor = 0.089365505
				set SAM_redshift = 10.19
			#15
			else if ($item  == 015) then 
				set SAM_scale_factor = 0.087412587
				set SAM_redshift = 10.44
			#14
			else if ($item  == 014) then 
				set SAM_scale_factor = 0.085470085
				set SAM_redshift = 10.7
			#13
			else if ($item == 013) then
				set SAM_scale_factor = 0.08361204
				set SAM_redshift = 10.96
			#12
			else if ($item == 012) then
				set SAM_scale_factor = 0.081766149
				set SAM_redshift = 11.23
			#11
			else if ($item == 011) then
				set SAM_scale_factor = 0.08
				set SAM_redshift = 11.5
			#10
			else if ($item  == 010) then 
				set SAM_scale_factor = 0.078308536
				set SAM_redshift = 11.77
			#9
			else if ($item == 009) then
				set SAM_scale_factor = 0.076569678
				set SAM_redshift = 12.06
			#8
			else if ($item == 008) then
				set SAM_scale_factor = 0.074906367
				set SAM_redshift = 12.35
			#7
			else if ($item == 007) then
				set SAM_scale_factor = 0.073206442
				set SAM_redshift = 12.66
			#6
			else if ($item == 006) then
				set SAM_scale_factor = 0.071428571
				set SAM_redshift = 13.0
			#5
			else if ($item == 005) then
				set SAM_scale_factor = 0.070077085
				set SAM_redshift = 6.51
			#4
			else if ($item == 004) then
				set SAM_scale_factor = 0.068493151
				set SAM_redshift = 13.27
			#3
			else if ($item  == 003) then 
				set SAM_scale_factor = 0.066979236
				set SAM_redshift = 13.93
			#2
			else if ($item == 002) then
				set SAM_scale_factor = 0.065616798
				set SAM_redshift = 14.24
			#1
			else if ($item == 001) then
				set SAM_scale_factor = 0.064102564
				set SAM_redshift = 14.6
			endif
			

		else if ($1 =~ *125) then

			if ($item == 107) then
				set SAM_scale_factor = 1.0
				set SAM_redshift = 0.0

			else if ($item == 98) then
				set SAM_scale_factor = 0.869565217
				set SAM_redshift = 0.15

			else if ($item == 92) then
				set SAM_scale_factor = 0.793650794
				set SAM_redshift = 0.26

			else if ($item == 80) then
				set SAM_scale_factor = 0.662251656
				set SAM_redshift = 0.51

			else if ($item == 70) then
				set SAM_scale_factor = 0.568181818
				set SAM_redshift = 0.76

			else if ($item == 61) then
				set SAM_scale_factor = 0.495049505
				set SAM_redshift = 1.02
		
			else if ($item == 47) then
				set SAM_scale_factor = 0.4
				set SAM_redshift = 1.5

			else if ($item == 35) then
				set SAM_scale_factor = 0.333333333
				set SAM_redshift = 2.0

			else if ($item == 27) then
				set SAM_scale_factor = 0.285714286
				set SAM_redshift = 2.5

			else if ($item == 21) then
				set SAM_scale_factor = 0.248756219
				set SAM_redshift = 3.02

			else if ($item == 17) then
				set SAM_scale_factor = 0.220264317
				set SAM_redshift = 3.54

			else if ($item == 14) then
				set SAM_scale_factor = 0.198019802
				set SAM_redshift = 4.05
		
			else if ($item == 12) then
				set SAM_scale_factor = 0.183486239
				set SAM_redshift = 4.45

			else if ($item == 9) then
				set SAM_scale_factor = 0.162337662
				set SAM_redshift = 5.16

			else if ($item == 6) then
				set SAM_scale_factor = 0.137174211
				set SAM_redshift = 6.29
		
			else if ($item == 4) then
				set SAM_scale_factor = 0.119474313
				set SAM_redshift = 7.37

			endif

		else
			set SAM_scale_factor = 'None'
			set SAM_redshift = 'None'
		endif

		echo $1 $item $SAM_redshift $SAM_scale_factor	>> $4
	end

	echo $1 $MARKER_CODE $COLOR_CODE $LINESTYLE_CODE $COLOR_MAP $MARKER_FACECOLOR_CODE  >> $2

else if ($1 =~ SAGE*) then
	set NAME_SL		=  'Link2_MDPL_SAGE_1Gpc'

	set COLOR_CODE		=  'E50000'
	set COLOR_CODE		=  'w'	
	set MARKER_CODE		=  'no'
	set LINESTYLE_CODE	=  ':'
	set COLOR_MAP		= 'Paired'

	#OII-sfr paper:
	set LINESTYLE_CODE	=  '-'

	set MARKER_FACECOLOR_CODE	= 'w'	
	set SUBDIRID		=  'False'

	set MASS_CORR		= 1e-10

	if ($5 == 'DS') then
		set LITTLE_H_CORR	= 0.6777
	else
		set LITTLE_H_CORR	= 1.475579165
		set GYR_CORR		= 1e-9

	endif

	set LITTLE_H_1_CORR	= 0.6777
	set e3_CORR		= 1000

	foreach item ($3)

		if ($item == 0.000) then
			set SAM_scale_factor = 1.0
			set SAM_redshift = 0.00

		else if ($item == 1205.hdf5) then
			set SAM_scale_factor = 0.6421
			set SAM_redshift = 0.56

		else if ($item =~ 1235*) then
			#350Mpc/h
			set SAM_scale_factor = 0.9151999999999999
			set SAM_redshift = 0.09

		else if ($item =~ 1238*) then
			set SAM_scale_factor = 0.6421
			set SAM_redshift = 0.56

		else if ($item =~ 1244*) then
			#350Mpc/h
			set SAM_scale_factor = 0.6421
			set SAM_redshift = 0.56

		else if ($item == 1265) then
			#250 Mpc/h
			set SAM_scale_factor = 0.4505
			set SAM_redshift = 1.22

		else if ($item == 1266) then
			#250 Mpc/h
			set SAM_scale_factor = 0.5747
			set SAM_redshift = 0.74

		else if ($item == 1267) then
			#250 Mpc/h
			set SAM_scale_factor = 0.628140704
			set SAM_redshift = 0.59

		else if ($item == 1268) then
			#250 Mpc/h
			set SAM_scale_factor = 0.3303
			set SAM_redshift = 2.03

		else if ($item == 1269) then
			#1000 Mpc/h
			set SAM_scale_factor = 0.628140704
			set SAM_redshift = 0.59

		else if ($item == 1270) then
			#1000 Mpc/h
			set SAM_scale_factor = 0.67114094
			set SAM_redshift = 0.49

		else if ($item =~ 1271*) then
			#250 Mpc/h apparent mags
			set SAM_scale_factor = 0.514403292
			set SAM_redshift = 0.94

		else if ($item == 1272) then
			#250 Mpc/h
			set SAM_scale_factor = 0.526
			set SAM_redshift = 0.9

		else if ($item == 1273) then
			#250 Mpc/h
			set SAM_scale_factor = 0.9151999999999999
			set SAM_redshift = 0.09

		else if ($item =~ 1283*) then
			#250 Mpc/h absolute mags
			set SAM_scale_factor = 0.514403292
			set SAM_redshift = 0.94
		#124
		else if ($item == 0.022) then
			set SAM_scale_factor = 0.980392157
			set SAM_redshift = 0.02
		#123
		else if ($item == 0.045) then
			set SAM_scale_factor = 0.952380952
			set SAM_redshift = 0.05

		#122
		else if ($item == 0.069) then
			set SAM_scale_factor = 0.934579439
			set SAM_redshift = 0.07

		else if ($item == 1182.hdf5) then
			set SAM_scale_factor = 0.9151999999999999
			set SAM_redshift = 0.09

		else if ($item  == 75) then 
			#dummy entry to plot more smoothly
			set SAM_scale_factor = 1.0
			set SAM_redshift = 0.0
		#121
		else if ($item == 0.093) then
			set SAM_scale_factor = 0.9151999999999999
			set SAM_redshift = 0.09
		#120
		else if ($item == 0.117) then
			set SAM_scale_factor = 0.892857142857
			set SAM_redshift = 0.12
		#119
		else if ($item  == 0.142) then
			set SAM_scale_factor = 0.8754999999999998
			set SAM_redshift = 0.14
		#118
		else if ($item  == 0.168) then
			set SAM_scale_factor = 0.854700855
			set SAM_redshift = 0.17
		#117
		else if ($item == 0.194) then
			set SAM_scale_factor = 0.837520938
			set SAM_redshift = 0.19
		#116
		else if ($item == 0.221) then
			set SAM_scale_factor = 0.819672131
			set SAM_redshift = 0.22
		#115
		else if ($item == 0.248) then
			set SAM_scale_factor = 0.8
			set SAM_redshift = 0.25
		#114
		else if ($item == 0.276) then
			set SAM_scale_factor = 0.78125
			set SAM_redshift = 0.28
		#113	
		else if ($item == 0.304) then
			set SAM_scale_factor = 0.766871166
			set SAM_redshift = 0.30
		#112	
		else if ($item == 0.334) then
			set SAM_scale_factor = 0.751879699
			set SAM_redshift = 0.33
		#111	
		else if ($item == 0.364) then
			set SAM_scale_factor = 0.735294118
			set SAM_redshift = 0.36
		#110	
		else if ($item == 0.394) then
			set SAM_scale_factor = 0.71942446
			set SAM_redshift = 0.39
		#109
		else if ($item == 0.425) then
			set SAM_scale_factor = 0.699300699
			set SAM_redshift = 0.43
		#108
		else if ($item == 0.457) then
			set SAM_scale_factor = 0.684931507
			set SAM_redshift = 0.46
		#107
		else if ($item == 0.490) then
			set SAM_scale_factor = 0.67114094
			set SAM_redshift = 0.49
		#106
		else if ($item == 0.523) then
			set SAM_scale_factor = 0.657894737
			set SAM_redshift = 0.52

		else if ($item == 1247) then
			set SAM_scale_factor = 0.66666666
			set SAM_redshift = 0.5
		#105
		else if ($item  == 0.557) then 
			set SAM_scale_factor = 0.6421
			set SAM_redshift = 0.56
		#104
		else if ($item  == 0.592) then 
			set SAM_scale_factor = 0.628140704
			set SAM_redshift = 0.59
		#103
		else if ($item  == 0.628) then 
			set SAM_scale_factor = 0.614250614
			set SAM_redshift = 0.63
		#102	
		else if ($item == 0.664) then
			set SAM_scale_factor = 0.600961538
			set SAM_redshift = 0.66
		#101	
		else if ($item == 0.702) then
			set SAM_scale_factor = 0.587544066
			set SAM_redshift = 0.7
		#100
		else if ($item  == 0.740) then 
			set SAM_scale_factor = 0.5747
			set SAM_redshift = 0.74
		#99
		else if ($item == 0.779) then
			set SAM_scale_factor = 0.562113547
			set SAM_redshift = 0.78
		#98
		else if ($item == 0.819) then
			set SAM_scale_factor = 0.549752611
			set SAM_redshift = 0.82
		#97
		else if ($item == 0.859) then
			set SAM_scale_factor = 0.537923615
			set SAM_redshift = 0.86
		#96
		else if ($item  == 0.901) then 
			set SAM_scale_factor = 0.526
			set SAM_redshift = 0.9
		#95
		else if ($item  == 0.944) then 
			set SAM_scale_factor = 0.514403292
			set SAM_redshift = 0.94
		#94
		else if ($item == 0.987) then
			set SAM_scale_factor = 0.503271263
			set SAM_redshift = 0.99
		#93
		else if ($item == 1.032) then
			set SAM_scale_factor = 0.492125984
			set SAM_redshift = 1.03
		#92
		else if ($item  == 1.077) then 
			set SAM_scale_factor = 0.481463649
			set SAM_redshift = 1.08
		#91
		else if ($item  == 1.124) then 
			set SAM_scale_factor = 0.470809793
			set SAM_redshift = 1.12
		#90
		else if ($item  == 1.172) then 
			set SAM_scale_factor = 0.460405157
			set SAM_redshift = 1.17
		#89
		else if ($item  == 1.220) then 
			set SAM_scale_factor = 0.4505
			set SAM_redshift = 1.22
		#88
		else if ($item  == 1.270) then 
			set SAM_scale_factor = 0.440528634
			set SAM_redshift = 1.27
		#87
		else if ($item  == 1.321) then 
			set SAM_scale_factor = 0.430848772
			set SAM_redshift = 1.32
		#86
		else if ($item  == 1.372) then 
			set SAM_scale_factor = 0.42158516
			set SAM_redshift = 1.372	
		#85
		else if ($item  == 1.425) then 
			set SAM_scale_factor = 0.411522634
			set SAM_redshift = 1.43
		#84
		else if ($item  == 1.480) then 
			set SAM_scale_factor = 0.403225806
			set SAM_redshift = 1.48
		#83
		else if ($item  == 1.535) then 
			set SAM_scale_factor = 0.393700787
			set SAM_redshift = 1.54
		#82
		else if ($item  == 1.593) then 
			set SAM_scale_factor = 0.386100386
			set SAM_redshift = 1.59
		#81
		else if ($item  == 1.650) then 
			set SAM_scale_factor = 0.377358491
			set SAM_redshift = 1.65
		#80
		else if ($item  == 1.710) then 
			set SAM_scale_factor = 0.36900369
			set SAM_redshift = 1.71
		#79
		else if ($item  == 1.771) then 
			set SAM_scale_factor = 0.36101083
			set SAM_redshift = 1.77
		#78
		else if ($item  == 1.833) then
			set SAM_scale_factor = 0.35335689
			set SAM_redshift = 1.83
		#77
		else if ($item == 1.896) then
			set SAM_scale_factor = 0.344827586
			set SAM_redshift = 1.9
		#76
		else if ($item == 1.961) then
			set SAM_scale_factor = 0.337837838
			set SAM_redshift = 1.96
		#75
		else if ($item  == 2.028) then 
			set SAM_scale_factor = 0.3303
			set SAM_redshift = 2.03
		#74
		else if ($item  == 2.095) then 
			set SAM_scale_factor = 0.322580645
			set SAM_redshift = 2.1
		#73
		else if ($item  == 2.165) then 
			set SAM_scale_factor = 0.315457413
			set SAM_redshift = 2.17
		#72
		else if ($item  == 2.235) then 
			set SAM_scale_factor = 0.308641975
			set SAM_redshift = 2.24
		#71
		else if ($item  == 2.308) then 
			set SAM_scale_factor = 0.302114804
			set SAM_redshift = 2.31
		#70
		else if ($item  == 2.382) then 
			set SAM_scale_factor = 0.295857988
			set SAM_redshift = 2.38
		#69
		else if ($item  == 2.458) then 
			set SAM_scale_factor = 0.289017341
			set SAM_redshift = 2.46
		#68
		else if ($item  == 2.535) then 
			set SAM_scale_factor = 0.282485876
			set SAM_redshift = 2.54
		#67
		else if ($item  == 2.614) then 
			set SAM_scale_factor = 0.27700831
			set SAM_redshift = 2.61
		#66
		else if ($item  == 2.695) then 
			set SAM_scale_factor = 0.27027027
			set SAM_redshift = 2.7
		#65
		else if ($item  == 2.778) then 
			set SAM_scale_factor = 0.264550265
			set SAM_redshift = 2.78
		#64
		else if ($item  == 2.862) then 
			set SAM_scale_factor = 0.259067358
			set SAM_redshift = 2.86
		#63
		else if ($item  == 2.949) then 
			set SAM_scale_factor = 0.253164557
			set SAM_redshift = 2.95
		#62
		else if ($item  == 3.037) then 
			set SAM_scale_factor = 0.247524752
			set SAM_redshift = 3.04
		#61
		else if ($item  == 3.127) then 
			set SAM_scale_factor = 0.242130751
			set SAM_redshift = 3.13
		#60
		else if ($item  == 3.221) then 
			set SAM_scale_factor = 0.236966825
			set SAM_redshift = 3.22
		#59
		else if ($item  == 3.314) then 
			set SAM_scale_factor = 0.232018561
			set SAM_redshift = 3.31
		#58
		else if ($item  == 3.411) then 
			set SAM_scale_factor = 0.22675737
			set SAM_redshift = 3.41
		#57
		else if ($item  == 3.511) then 
			set SAM_scale_factor = 0.22172949
			set SAM_redshift = 3.51
		#56
		else if ($item  == 3.610) then 
			set SAM_scale_factor = 0.21691974
			set SAM_redshift = 3.61
		#55
		else if ($item  == 3.715) then 
			set SAM_scale_factor = 0.211864407
			set SAM_redshift = 3.72
		#54
		else if ($item  == 3.819) then 
			set SAM_scale_factor = 0.20746888
			set SAM_redshift = 3.82
		#53
		else if ($item  == 3.929) then 
			set SAM_scale_factor = 0.202839757
			set SAM_redshift = 3.93
		#52
		else if ($item  == 4.038) then 
			set SAM_scale_factor = 0.198412698
			set SAM_redshift = 4.04
		#51
		else if ($item  == 4.152) then 
			set SAM_scale_factor = 0.194174757
			set SAM_redshift = 4.15
		#50
		else if ($item  == 4.266) then 
			set SAM_scale_factor = 0.189900
			set SAM_redshift = 4.27
		#49
		else if ($item  == 4.385) then 
			set SAM_scale_factor = 0.185700
			set SAM_redshift = 4.39
		#48
		else if ($item  == 4.507) then 
			set SAM_scale_factor = 0.181600
			set SAM_redshift = 4.51
		#47
		else if ($item  == 4.627) then 
			set SAM_scale_factor = 0.177700
			set SAM_redshift = 4.63
		#46
		else if ($item == 4.754) then
			set SAM_scale_factor = 0.173800
			set SAM_redshift = 4.75
		#45
		else if ($item == 4.882) then
			set SAM_scale_factor = 0.170000
			set SAM_redshift = 4.88
		#44
		else if ($item == 5.017) then
			set SAM_scale_factor = 0.166200
			set SAM_redshift = 5.02
		#43
		else if ($item  == 5.150) then 
			set SAM_scale_factor = 0.162600
			set SAM_redshift = 5.15
		#42
		else if ($item == 5.289) then
			set SAM_scale_factor = 0.159000
			set SAM_redshift = 5.29
		#41
		else if ($item == 5.431) then
			set SAM_scale_factor = 0.155500
			set SAM_redshift = 5.43
		#40
		else if ($item == 5.575) then
			set SAM_scale_factor = 0.152100
			set SAM_redshift = 5.57
		#39
		else if ($item  == 5.720) then 
			set SAM_scale_factor = 0.148800
			set SAM_redshift = 5.72
		#38
		else if ($item  == 5.873) then 
			set SAM_scale_factor = 0.145500
			set SAM_redshift = 5.87
		#37
		else if ($item  == 6.022) then 
			set SAM_scale_factor = 0.142400
			set SAM_redshift = 6.02
		#36
		else if ($item  == 6.184) then 
			set SAM_scale_factor = 0.139200
			set SAM_redshift = 6.18
		#35
		else if ($item == 6.342) then
			set SAM_scale_factor = 0.136200
			set SAM_redshift = 6.34
		#34
		else if ($item == 6.508) then
			set SAM_scale_factor = 0.133200
			set SAM_redshift = 6.51
		#33
		else if ($item == 6.675) then
			set SAM_scale_factor = 0.130300
			set SAM_redshift = 6.67
		#32
		else if ($item  == 6.849) then 
			set SAM_scale_factor = 0.127400
			set SAM_redshift = 6.85
		#31
		else if ($item == 7.026) then
			set SAM_scale_factor = 0.124600
			set SAM_redshift = 7.03
		#30
		else if ($item == 7.203) then
			set SAM_scale_factor = 0.121900
			set SAM_redshift = 7.2
		#29
		else if ($item == 7.389) then
			set SAM_scale_factor = 0.119200
			set SAM_redshift = 7.39
		#28
		else if ($item == 7.576) then
			set SAM_scale_factor = 0.116550117
			set SAM_redshift = 7.58
		#27
		else if ($item  == 7.764) then 
			set SAM_scale_factor = 0.114155251
			set SAM_redshift = 7.76
		#26
		else if ($item  == 7.961) then 
			set SAM_scale_factor = 0.111607143
			set SAM_redshift = 7.96
		#25
		else if ($item  == 8.166) then 
			set SAM_scale_factor = 0.109051254
			set SAM_redshift = 8.17
		#24
		else if ($item  == 8.372) then 
			set SAM_scale_factor = 0.106723586
			set SAM_redshift = 8.37
		#23
		else if ($item == 8.579) then
			set SAM_scale_factor = 0.104384134
			set SAM_redshift = 8.58
		#22
		else if ($item == 8.794) then
			set SAM_scale_factor = 0.102040816
			set SAM_redshift = 8.8
		#21
		else if ($item == 9.010) then
			set SAM_scale_factor = 0.0999001
			set SAM_redshift = 9.01
		#20
		else if ($item  == 9.235) then 
			set SAM_scale_factor = 0.09765625
			set SAM_redshift = 9.24
		#19
		else if ($item == 9.471) then
			set SAM_scale_factor = 0.095510984
			set SAM_redshift = 9.47
		#18
		else if ($item == 9.707) then
			set SAM_scale_factor = 0.093370682
			set SAM_redshift = 9.71
		#17
		else if ($item == 9.941) then
			set SAM_scale_factor = 0.091407678
			set SAM_redshift = 9.94
		#16
		else if ($item  == 10.186) then 
			set SAM_scale_factor = 0.089365505
			set SAM_redshift = 10.19
		#15
		else if ($item  == 10.442) then 
			set SAM_scale_factor = 0.087412587
			set SAM_redshift = 10.44
		#14
		else if ($item  == 10.696) then 
			set SAM_scale_factor = 0.085470085
			set SAM_redshift = 10.7
		#13
		else if ($item == 10.962) then
			set SAM_scale_factor = 0.08361204
			set SAM_redshift = 10.96
		#12
		else if ($item == 11.225) then
			set SAM_scale_factor = 0.081766149
			set SAM_redshift = 11.23
		#11
		else if ($item == 11.500) then
			set SAM_scale_factor = 0.08
			set SAM_redshift = 11.5
		#10
		else if ($item  == 11.771) then 
			set SAM_scale_factor = 0.078308536
			set SAM_redshift = 11.77
		#9
		else if ($item == 12.055) then
			set SAM_scale_factor = 0.076569678
			set SAM_redshift = 12.06
		#8
		else if ($item == 12.351) then
			set SAM_scale_factor = 0.074906367
			set SAM_redshift = 12.35
		#7
		else if ($item == 12.661) then
			set SAM_scale_factor = 0.073206442
			set SAM_redshift = 12.66
		#6
		else if ($item == 12.966) then
			set SAM_scale_factor = 0.071428571
			set SAM_redshift = 13.0
		#5
		else if ($item == 13.265) then
			set SAM_scale_factor = 0.070077085
			set SAM_redshift = 6.51
		#4
		else if ($item == 13.599) then
			set SAM_scale_factor = 0.068493151
			set SAM_redshift = 13.27
		#3
		else if ($item  == 13.925) then 
			set SAM_scale_factor = 0.066979236
			set SAM_redshift = 13.93
		#2
		else if ($item == 14.244) then
			set SAM_scale_factor = 0.065616798
			set SAM_redshift = 14.24
		#1
		else if ($item == 14.601) then
			set SAM_scale_factor = 0.064102564
			set SAM_redshift = 14.6

		else
			set SAM_scale_factor = 'None'
			set SAM_redshift = 'None'
		endif

		
		echo $1 $item $SAM_redshift $SAM_scale_factor	>> $4
	end

	echo $1 $MARKER_CODE $COLOR_CODE $LINESTYLE_CODE $COLOR_MAP $MARKER_FACECOLOR_CODE	 >> $2


else if ($1 =~ ySAM*) then
	set NAME_SL		=  'False'
	set MARKER_CODE		=  'p'
	set COLOR_CODE		=  'ef4f91'
	set SUBDIRID		=  'False'
	set LITTLE_H		= 0.6777

	echo $1 $MARKER_CODE $COLOR_CODE			>> $2

else
	set NAME_SL		=  'False'
	set MARKER_CODE		=  'p'
	set COLOR_CODE		=  'e4354b'
	set COLOR_CODE		=  'k'
	set SUBDIRID		=  'False'
	set LINESTYLE_CODE	=  '-'
	set COLOR_MAP		= 'YlOrBr'
	set COLOR_MAP		= 'no'
	set MARKER_FACECOLOR_CODE = 'w'
	set LITTLE_H		= 0.7
	set SIM_LITTLE_H	= 0.7

	#Kindersicherung!
	foreach item ($3)
		#echo $1 $item -99.99 1.0	>> $4
		#echo $1 $item 0.0 1.0	>> $4
		#echo $1 $item 0.43 0.68027	>> $4
		#echo $1 $item 0.47 0.6993	>> $4
		echo $1 $item 0.5 0.6666667	>> $4
		#echo $1 $item 0.51 0.66225	>> $4
		#echo $1 $item 0.52 0.657895	>> $4
		#echo $1 $item 0.523 0.656599	>> $4
		#echo $1 $item 0.53 0.653595	>> $4
		#echo $1 $item 0.54 0.64935	>> $4
		#echo $1 $item 0.55 0.64516129	>> $4
		#echo $1 $item 0.554 0.6435	>> $4
		#echo $1 $item 0.67 0.5988	>> $4

	end

	echo $1 $MARKER_CODE $COLOR_CODE $LINESTYLE_CODE $COLOR_MAP $MARKER_FACECOLOR_CODE	>> $2

endif

echo $NAME_SL $SUBDIRID $SIM_LITTLE_H $GENERAL_LITTLE_H $LITTLE_H_CORR $MPC_CORR $MASS_CORR $GYR_CORR $LUM_UNITS_CORR '1' $CONV_TO_AB_MAG $KMS1_CORR $COMV_CORR $APPLY_K_CORR_APP $APPLY_Z_BOOST $LITTLE_H_1_CORR $e3_CORR $MAB_CORR $LITTLE_H_2_CORR $APPLY_K_CORR_ABS $IMF_CORR $GYR1_CORR

