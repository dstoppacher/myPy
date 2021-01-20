#!/bin/tcsh

#####################################################################################################
#	How to change a keyword:
#	------------------------------------------
#
#	'keyword_change1=default'
#
#	instead of default use the name of the keywords you want to change e.g. plot_legend0 
#	and write the new string of the "plot_legend0" should be assigend after a comma:
#
#	'keyword_change1=plot_legend0,my new legend'
#
#####################################################################################################

#echo '1:' $1 '2:' $2 '3:' $3 '4:' $4 '5:' $5 '6:' $6 '7:' $7 '8:' $8 '9:' $9 '10:' $10 '11:' $11 '12:' $12 '13:' $13 '14:' $14 '15:' $15 '16:' $16 '17:' $17 '18:' $18 '19:' $19 '20:' $20 '21:' $21 '22:' $22 '23:' $23

echo 'plot-config:' $24
set demo = 'False'
set reduced = 'False'
#Set DEFAULT settings for all plots
echo 'plot_type= default'					>> $2
echo 'title= ' 							>> $2
#CMASS paper l=0.12, r=0.98, b=0.98, t=0.99
#l=0.18 for contour plots
echo 'adjust_left= 0.14'					>> $2
echo 'adjust_right= 0.99'					>> $2
echo 'adjust_bottom= 0.14'					>> $2
echo 'adjust_top= 0.99'					>> $2
echo 'share_axis= no'						>> $2
echo 'share_axis_x= no'						>> $2
echo 'share_axis_y= no'						>> $2
echo 'share_axis_offset_top= 0.0'					>> $2
echo 'share_axis_offset_right= 0.0'				>> $2
echo 'share_which_axis= '						>> $2
echo 'share_which_col= '						>> $2
echo 'share_axis_sub= no'						>> $2
echo 'share_which_axis_sub= '					>> $2
echo 'share_which_col_sub= '					>> $2
echo 'share_x_range_min= 1'					>> $2
echo 'share_x_range_max= 10'					>> $2
echo 'share_y_range_min= 1'					>> $2
echo 'share_y_range_max= 10'					>> $2
echo 'log_scale_share_x= no'					>> $2
echo 'log_scale_share_y= no'					>> $2
echo 'add_axis_steps_sub= 6'					>> $2
echo 'add_range_min= min'					>> $2
echo 'add_range_max= max'					>> $2
echo 'add_axis_plot_exact_ticks= True'				>> $2
echo 'add_axis= no'						>> $2
echo 'add_axis_offset_top= 0.0'					>> $2
echo 'add_axis_offset_right= 0.0'				>> $2
echo 'add_which_axis= '						>> $2
echo 'add_which_col= '						>> $2
echo 'add_axis_steps_x= 6'					>> $2
echo 'add_axis_steps_y= 12'					>> $2
echo 'add_axis_sub= no'						>> $2
echo 'add_which_axis_sub= '					>> $2
echo 'add_which_col_sub= '					>> $2
echo 'log_scale_add_x= no'					>> $2
echo 'log_scale_add_y= no'					>> $2
echo 'add_axis_steps_sub= 6'					>> $2
echo 'add_range_min= min'					>> $2
echo 'add_range_max= max'					>> $2
echo 'add_axis_plot_exact_ticks= True'				>> $2
echo 'x_title= '	         				>> $2
echo 'y_title= '	         				>> $2
echo 'z_title= False'	         				>> $2
echo 'subtitle= False'	         				>> $2

#space from bottom to x-label, default 0.01, 0.02 (from 0-1 % of figure size)
echo 'xlabel_pad= 0.02'	         			>> $2
#space from left to y-label, default 0.0, 0.01 (from 0-1 % of figure size)
echo 'ylabel_pad= 0.01'	         			>> $2

#pos of x-label with reference to the axis, default 0.54, (0.5 is center, from 0-1 % of figure size)
echo 'xlabel_pos= 0.54'	         			>> $2
#pos of y-label with reference to the axis, default 0.53, (0.5 is center, from 0-1 % of figure size)
echo 'ylabel_pos= 0.53'	         			>> $2         

echo 'format_xticks= True'	         			>> $2
echo 'format_yticks= True'	         			>> $2

echo 'no_xticks= False'     		    			>> $2
echo 'no_last_xticks= False'       	  			>> $2
echo 'no_first_xticks= False'        	 			>> $2
echo 'no_first_last_xticks= False'         			>> $2

echo 'no_yticks= False'        					>> $2
echo 'no_last_yticks= False'         				>> $2
echo 'no_first_yticks= False'         				>> $2
echo 'no_first_last_yticks= False'         			>> $2

echo 'float_format_x= 0f'         				>> $2
#for 10^x labeling use 3f
echo 'float_format_y= 3f'         				>> $2

echo 'tight_layout= False'		 		 	>> $2

echo 'caption1= '				 		>> $2
echo 'caption2= '				 		>> $2
echo 'caption3= '				 		>> $2

echo 'hratio1= 1'				 		 >> $2
echo 'hratio2= 2'				 		 >> $2

echo 'x_range_min= min'			 			 >> $2
echo 'x_range_max= max'			 		 	 >> $2
echo 'y_range_min= min'		 			 	 >> $2
echo 'y_range_max= max'				 		 >> $2
echo 'z_range_min= min'		 			 	 >> $2
echo 'z_range_max= max'				 		 >> $2

echo 'keyword_change1= default'		 			>> $2
echo 'keyword_change2= default'		 			>> $2
echo 'keyword_change3= default'		 			>> $2
echo 'keyword_change4= default'		 			>> $2
echo 'keyword_change5= default'		 			>> $2
echo 'keyword_change6= default'		 			>> $2
echo 'keyword_change7= default'		 			>> $2
echo 'keyword_change8= default'		 			>> $2
echo 'keyword_change9= default'		 			>> $2
echo 'keyword_change10= default'	 			>> $2
echo 'keyword_change11= default'	 			>> $2
echo 'keyword_change12= default'	 			>> $2
echo 'keyword_change12= default'	 			>> $2
echo 'keyword_change13= default'	 			>> $2
echo 'keyword_change14= default'	 			>> $2

echo 'title_fontsize= 20'				 	 >> $2
#default 32
echo 'axis_label_fontsize= 32'				 	 >> $2
#default 26
echo 'text_fontsize= 28'					>> $2
#default 30				 	 
echo 'axis_ticks_fontsize= 30'				 	 >> $2
#default 28
echo 'legend_fontsize= 28'				 	 >> $2
echo 'legend_ncols= 1'						>> $2
echo 'legend_fancy= True'					>> $2
echo 'legend_handletextpad= 0.8'					>> $2

echo 'marker= no'					 	 >> $2
#Default 15
echo 'markersize= 15'				 		 >> $2
echo 'markersize_offset= 5'				 		 >> $2
echo 'marker_add_subplot= o'					 >> $2
echo 'markersize_subplot= 13'				 	 >> $2

echo 'mylw= 10'				 		 	 >> $2
echo 'mylw_offset= 4'				 		 	 >> $2
echo 'alpha_offset= 0'		 		 	 >> $2
#which line should be thinner (start with 0)
echo 'set_thin_line= 99'				 	>> $2
#How much thinner that the #mylw'
echo 'thin_line= 4'				 		 >> $2
#which line should be thinner (start with 0)
echo 'set_thin_line2= 99'				 	>> $2
#How much thinner that the #mylw'
echo 'thin_line2= 4'				 		 >> $2

echo 'mymarkeredgewidth= 2.5'				 	 >> $2
echo 'mymarkerfacecolor= w'				 	 >> $2
echo 'myframe_lw= 3'					 	 >> $2
echo 'lw_offset= -6'				 		 >> $2
echo 'size_x= 12'				 		 >> $2
echo 'size_y= 9'				 		 >> $2

echo 'nr_plots_x= 1'				 		 >> $2
echo 'nr_plots_y= 1'				 		 >> $2
echo 'nr_cats= 	1'				 		 >> $2

#colorbar/hexbin/countour details
echo 'cb_min= 1'				 		 >> $2
echo 'cb_max= 2.5'				 		 >> $2
echo 'cb_steps= 0.5'				 		 >> $2
echo 'bins= log'				 		 >> $2
echo 'min_count= 10'				 		 >> $2
echo 'gridsize= 100'				 		 >> $2
echo 'colorbar_anchor_right= 0.4'		 		 >> $2

set plot_CUT = 'False'
echo 'plot_CUT= '$plot_CUT >> $2

echo 'print_redshift= ' >> $2
echo 'print_redshift2= ' >> $2
echo 'print_redshift3= ' >> $2
echo 'quality_dpi= 100'			 		 >> $2

#set ranges for histogram panels on top and left of a contour plot
echo 'histo_panels= no'			>> $2
echo 'histo_num_density_nbins_top= 30'			>> $2
echo 'histo_num_density_nbins_left= 30'			>> $2
echo 'histo_num_density_min= 0.0'				>> $2
echo 'histo_num_density_max= 0.22'				>> $2
echo 'histo_num_density_number_major_ticks= 6'			>> $2

echo 'plot_lines_min_max_population= no'		 	 >> $2

######################################################################################################################################################
#select key for plotXY: 'SMF', 'LF', 'wp', 'wp_mstar', 'sfr2z', 'HOD', or 'rhalf'

set OII_ls_set = 'False'
set plotXY_key = 'sfr2z'
set hod_plot_key = 'False'
set use_cb_colors = 6
echo 'legend_ncols= 1'			>> $2
echo 'legend_fancy= False'		>> $2
######################################################################################################################################################

#marker style table: 	MATPLOT_MARKERSTYLE
echo '0= '	>> $15
echo '1= '	>> $15
echo '2= '	>> $15
echo '3= '	>> $15
echo '4= '	>> $15
echo '5= '	>> $15
echo '6= '	>> $15
echo '7= '	>> $15
echo '8= '	>> $15
echo '9= '	>> $15
echo '10= '	>> $15
echo '11= '	>> $15
echo '12= '	>> $15
echo '13= '	>> $15
echo '14= '	>> $15
echo '15= '	>> $15
echo '16= '	>> $15

#marker colour table: 	MATPLOT_MARKERCOL
echo '0= w'		>> $16
echo '1= k'		>> $16
echo '2= #ffe200'	>> $16
echo '3= w'		>> $16
echo '4= w'		>> $16
echo '5= w'		>> $16
echo '6= w' 		>> $16
echo '7= w' 		>> $16
echo '8= #ffe200'	>> $16
echo '9= #ffe200'	>> $16
echo '10= #ffe200'	>> $16
echo '11= #ffe200'	>> $16
echo '12= #ffe200'	>> $16
echo '13= #ffe200'	>> $16
echo '14= #ffe200'	>> $16
echo '15= #ffe200'	>> $16
echo '16= #ffe200'	>> $16

#ADD_X_AXIS
#yes= an additional x-axis of that data set will be plotted on top
echo '0= no'		>> $22
echo '1= no'		>> $22
echo '2= no'		>> $22
echo '3= no'		>> $22
echo '4= no'		>> $22
echo '5= no'		>> $22
echo '6= no'		>> $22
echo '7= no'		>> $22
echo '8= no'		>> $22
echo '9= no'		>> $22
echo '10= no'		>> $22
echo '11= no'		>> $22
echo '12= no'		>> $22
echo '13= no'		>> $22
echo '14= no'		>> $22
echo '15= no'		>> $22
echo '16= no'		>> $22

#ADD_Y_AXIS
#yes= an additional y-axis of that data set will be plotted at the right hand side
echo '0= no'		>> $23
echo '1= no'		>> $23
echo '2= no'		>> $23
echo '3= no'		>> $23
echo '4= no'		>> $23
echo '5= no'		>> $23
echo '6= no'		>> $23
echo '7= no'		>> $23
echo '8= no'		>> $23
echo '9= no'		>> $23
echo '10= no'		>> $23
echo '11= no'		>> $23
echo '12= no'		>> $23
echo '13= no'		>> $23
echo '14= no'		>> $23
echo '15= no'		>> $23
echo '16= no'		>> $23

#plot legend: 	MATPLOT_PLOTLEGEND
#no= no legend will be plotted for that entry
echo '0= yes'		>> $21
echo '1= yes'		>> $21
echo '2= yes'		>> $21
echo '3= yes'		>> $21
echo '4= yes'		>> $21
echo '5= yes'		>> $21
echo '6= yes'		>> $21
echo '7= yes'		>> $21
echo '8= yes'		>> $21
echo '9= yes'		>> $21
echo '10= yes'		>> $21
echo '11= yes'		>> $21
echo '12= yes'		>> $21
echo '13= yes'		>> $21
echo '14= yes'		>> $21
echo '15= yes'		>> $21
echo '16= yes'		>> $21

echo 'print_redshift_sep= ' >> $2
echo 'legend_sep_xpos= 0.55'	>> $2
echo 'legend_sep_ypos= 0.14'	>> $2
echo 'legend_sep_offset= 0.038'	>> $2

#marker alpha: 	MATPLOT_ALPHA
echo '0= 0.5'		>> $20
echo '1= 1'		>> $20
echo '2= 1'		>> $20
echo '3= 0.9'		>> $20
echo '4= 0.5'		>> $20
echo '5= 1'		>> $20
echo '6= 1'		>> $20
echo '7= 0.4'		>> $20
echo '8= 0.3'		>> $20
echo '9= 1'		>> $20
echo '10= 1'		>> $20
echo '11= 1'		>> $20
echo '12= 1'		>> $20
echo '13= 1'		>> $20
echo '14= 1'		>> $20
echo '15= 1'		>> $20
echo '16= 1'		>> $20

#marker colourmap: 	MATPLOT_COLORMAP
#no= no colormap will be used
echo '0= binary'		>> $17
echo '1= coolwarm'	>> $17
echo '2= no'		>> $17
echo '3= no'		>> $17
echo '4= no'		>> $17
echo '5= no'		>> $17
echo '6= no'		>> $17
echo '7= no'		>> $17
echo '8= no'		>> $17
echo '9= no'		>> $17
echo '10= no'		>> $17
echo '11= no'		>> $17
echo '12= no'		>> $17
echo '13= no'		>> $17
echo '14= no'		>> $17
echo '15= no'		>> $17
echo '16= no'		>> $17

#marker colourbar: 	MATPLOT_COLORBAR
#no= no colorbar will be plotted
#ALWAYS USE COLOUR BAR WITH HEXBINS!!!
echo '0= no'		>> $18
echo '1= no'		>> $18
echo '2= no'		>> $18
echo '3= no'		>> $18
echo '4= no'		>> $18
echo '5= no'		>> $18
echo '6= no'		>> $18
echo '7= no'		>> $18
echo '8= no'		>> $18
echo '9= no'		>> $18
echo '10= no'		>> $18
echo '11= no'		>> $18
echo '12= no'		>> $18
echo '13= no'		>> $18
echo '14= no'		>> $18
echo '15= no'		>> $18
echo '16= no'		>> $18

#subplot type: 		MATPLOT_PLOTTYPE
set key = 'histosf'
if ($24 =~ *hist*) then
	echo '0= barplot' 	>> $19
	echo '1= barplot' 	>> $19
	echo '2= barplot' 	>> $19
	echo '3= barplot' 	>> $19
	echo '4= barplot' 	>> $19
	echo '5= barplot' 	>> $19
	echo '6= barplot' 	>> $19
	echo '7= barplot' 	>> $19
	echo '8= barplot' 	>> $19
	echo '9= barplot' 	>> $19
else if ($key =~ *violi*) then
	echo '0= violin' 	>> $19
	echo '1= violin' 	>> $19
	echo '2= violin' 	>> $19
	echo '3= violin' 	>> $19
	echo '4= violin' 	>> $19
	echo '5= violin' 	>> $19
	echo '6= violin' 	>> $19
	echo '7= violin' 	>> $19
	echo '8= violin' 	>> $19
	echo '9= violin' 	>> $19
	echo '10= violin' 	>> $19
else
	echo '0= hexbins' 	>> $19
	echo '1= contour' 	>> $19
	echo '2= contour' 	>> $19
	echo '3= contour' 	>> $19
	echo '4= contour' 	>> $19
	echo '5= contour' 	>> $19
	echo '6= contour' 	>> $19
endif
echo '7= default' 	>> $19
echo '8= default' 	>> $19
echo '9= default' 	>> $19
echo '10= default' 	>> $19
echo '11= default' 	>> $19
echo '12= default' 	>> $19
echo '13= default' 	>> $19
echo '14= default' 	>> $19
echo '15= default' 	>> $19
echo '16= default' 	>> $19

#echo 'CATNAME $5' $5
if ($5 =~ 'Gal'*) then
	set catname = 'Gal'
else if ($5 == 'SAG') then
	set catname = 'SAG'
else if ($5 =~ 'SAGE'*) then
	set catname = 'SAGE'
else
	set catname = $5
endif

#echo '$5' $5 'catname:' $catname

#Choose dimension of the plot: yes:2D-plot, else: 3D-plot
echo '2D-plot= yes'					>> $2

#Choose the standard implementation of colors for the SAMs: True! Use a user defint color set: False!
echo 'use_cat_colors= False'					>> $2
echo 'subplot_use_cat_colors= False'				>> $2

echo 'subplot_is_obs= True'	 		 		 >> $2
echo 'xticks= default'			 			>> $2
echo 'xticks_minor= 0.1'				>> $2
echo 'xticks_minor_offset= 0'				>> $2

echo 'yticks= default'			 			>> $2
echo 'yticks_minor= 0.1'				>> $2
echo 'yticks_minor_offset= 0'				>> $2

#put minor ticks in the plot, if 'no' then no minor ticks will be plotted, if a float is given this will be the spacing of the thicks
echo 'minor_ticks_x_space= 0.1'			 	>> $2
echo 'minor_ticks_x_set_top= on'				>> $2
echo 'minor_ticks_x_set_bottom= on'				>> $2

echo 'minor_ticks_y_space= 0.1'			 	>> $2
echo 'minor_ticks_y_set_left= on'				>> $2
echo 'minor_ticks_y_set_right= on'				>> $2

echo 'minor_ticks_x_share_space= no'			 	>> $2
echo 'minor_ticks_y_share_space= no'			 	>> $2

if ($HOME == /home/doriss) then
	echo 'minor_ticks_x_space= 10'			 	>> $2
	echo 'minor_ticks_y_space= 10'			 	>> $2

endif

echo 'cut_part1= ' 				>> $2
echo 'cut_part2= '			>> $2
echo 'linestyle= '			 			 >> $2
echo 'linestyle_change= '		 			 >> $2
echo 'linestyle_sub= '	 			 >> $2
echo 'linestyle_sub_change= '		 			 >> $2

echo 'log_scale_x= no'						 >> $2
echo 'log_scale_y= no'						 >> $2
echo 'log_scale_x_sub= no'					 >> $2
echo 'log_scale_y_sub= no'					 >> $2

echo 'error_bars_x= no'					 	 >> $2
echo 'error_bars_y= no'					 >> $2
echo 'filled_between= True'					 >> $2

#if 'False', a shaded region if plotted, if 'True' then two lines in the same linestyle as the data is plotted without an filled
echo 'fill_no_facecolor= False' 				>> $2
echo 'errorbars= False'				 		 >> $2

echo 'error_bars_x_sub= no'					 >> $2
echo 'error_bars_y_sub= yes'					 >> $2
echo 'subplot_errorbars= False' 		 			 >> $2
echo 'filled_between_subplot= True' 		 		 >> $2

#filling intensity of filled errorbars 0=0% 1=100% of the color chosen
echo 'alpha= 0.05'					 >> $2
echo 'alpha_sub= 0.05'					 >> $2

echo 'show_plot= yes'			 			 >> $2
echo 'print_plot= yes'			 			 >> $2

#if pint_pdf=yes, a PDF is created otherwise a svg-file with 100 dpi
echo 'print_pdf= yes'			 	 >> $2
echo 'figure_number= '				 >> $2

#define an optional scalebar which width give you the scale of the axis
echo 'add_scalebar= False'			>> $2
echo 'scalebar_loc_code= 6'			>> $2
echo 'scalebar_width_x= 0.1'			>> $2
echo 'scalebar_width_y= 0.1'			>> $2
echo 'plot_grid= no'			>> $2

#Plot a horizontal line: False=no line plotted, other: at that position you entered
echo 'hline_ypos= False'					 >> $2	

#Plot a vertical line: False=no line plotted, other: at that position you entered
echo 'vline_xpos= False'					 >> $2	

#Choose contour-plot specifications (NOTE: Only for catalogs with ngal>5e6, else default values ware taken!)
echo 'nlevels= 5'			>> $2
echo 'contour_histo_nbins= 50'		>> $2
#SAG v3
#echo 'nlevels= 7'			>> $2
#echo 'contour_histo_nbins= 60'		>> $2
echo 'contour_log= yes'		>> $2
echo 'contour_log_wich_axis= x'		>> $2
echo 'plot_colorbar= yes'		>> $2
echo 'set_percent_to_zero= 100'		>> $2

#Scale color bar of contour plot with max(x)-max(x)/scale_color_bar, if scale_color_bar chlose to cero --> no scaling
echo 'scale_color_bar= 0.0'		>> $2
echo 'plot_colormap= yes'		>> $2

echo 'plot_legend= yes'		 			>> $2
echo 'legend_position= upper right'			 	 >> $2

#use location coordinates in float (better for contour plots), otherwise use 'loc-code': e.g. "upper right"
echo 'use_loc_co= yes'			>> $2
#is 'use_loc_co=yes, use the *_x and *_y to enter manually the location coordinates (0.0-1.0 values)
#lower left
echo 'use_loc_co_x= 0.03'			>> $2
echo 'use_loc_co_y= 0.03'			>> $2
#lower right
#echo 'use_loc_co_x= 0.50'			>> $2
#echo 'use_loc_co_y= 0.05'			>> $2
#upper right
echo 'use_loc_co_x= 0.47'			>> $2
echo 'use_loc_co_y= 0.60'			>> $2
#upper left
#echo 'use_loc_co_x= 0.03'			>> $2
#echo 'use_loc_co_y= 0.55'			>> $2


#True if you want to plot 1 or more small plots into a bigger plot, else False
echo 'plot_into_plot= False'		>> $2

#default x/y-position for redshift default=0.15/0.85 for one plot/figure
echo 'z_print_position_x= 0.16'		>> $2
echo 'z_print_position_y= 0.16'		>> $2
echo 'z_print_position_x2= 0.16'		>> $2
echo 'z_print_position_y2= 0.16'		>> $2
echo 'z_print_position_x3= 0.16'		>> $2
echo 'z_print_position_y3= 0.16'		>> $2

set prefix = 'LG'
#set prefix = 'SAGE'
set prefix = 'z='

echo 'keyword_change1= plot_legend0,'$prefix'-down' 	>> $2
echo 'keyword_change1= plot_legend0,'$prefix'-dens $M_*$'  	>> $2
#echo 'keyword_change1= plot_legend0,'$prefix'-dens $SFR$'  	>> $2
#echo 'keyword_change1= plot_legend0,'$prefix'-dens $sSFR$' 	>> $2
echo 'keyword_change1= plot_legend0,'$prefix'-dens $V_{max}$' 	>> $2
echo 'keyword_change1= plot_legend0,'$prefix'-dens $V_{max_{sats}}$' 	>> $2

echo 'keyword_change1= plot_legend0,'$prefix'-down' 	>> $2
echo 'keyword_change2= plot_legend1,'$prefix'-dens $M_*$'  	>> $2
echo 'keyword_change2= plot_legend1,'$prefix'-dens $SFR$ (<)'  	>> $2
echo 'keyword_change3= plot_legend2,'$prefix'-dens $sSFR$ (<)' 	>> $2
#echo 'keyword_change3= plot_legend2,'$prefix'-dens $V_{max}$' 	>> $2
#echo 'keyword_change1= plot_legend0,'$prefix'-dens $V_{max_{sats}}$' 	>> $2

if ($1 == plotXY && $plotXY_key == 'HOD') then

	echo 'title= HOD$ '$3' '$4 >> $2
	echo 'x_title= $\log_{10}$ ($M_{200c}$ $[M_{\odot}]$)'	 >> $2
	echo 'y_title= $\log_{10}$ <$N_{gal}$>'		 >> $2
	#echo 'y_title= '		 >> $2
	#echo 'no_yticks= True'        					>> $2
	#echo 'no_first_yticks= True'        					>> $2
	#echo 'x_title= '		 >> $2
	echo 'no_xticks= False'        					>> $2
	#echo 'no_first_xticks= True'       	  			>> $2

	echo 'legend_ncols= 1'			>> $2	
	set use_cb_colors = 1
	echo 'error_bars_y= yes'					 >> $2
	echo 'size_y= 8'				 		 >> $2
	echo 'size_x= 9'				 		 >> $2

	set hod_plot_key = 'CMASS'
	if ($hod_plot_key == 'CMASS') then
		#echo 'xticks= 12,13,14,15'		>> $2
		echo 'minor_ticks_x_space= 0.1'			 	>> $2
		echo 'minor_ticks_y_space= 0.1'			 	>> $2
		echo 'x_range_min= 12.5'				 	 	 >> $2
		echo 'x_range_max= 15.0'				 	 >> $2
		echo 'y_range_min= -1.7'			 			 >> $2
		echo 'y_range_max= 1.2'			 		 >> $2
		#echo 'hline_ypos= 0.0'					 >> $2
		echo 'no_last_yticks= True'       	  			>> $2

		echo 'keyword_change1= plot_legend0,Gal-dens  (15% sats)' 	>> $2
		echo 'keyword_change4= plot_legend3,Gal400-dens  (14% sats)' 	>> $2

	else
		echo 'xticks= 10,11,12,13,14,15'		>> $2
		echo 'minor_ticks_x_space= 0.1'			 	>> $2
		echo 'minor_ticks_y_space= 0.2'			 	>> $2
		echo 'x_range_min= 10'				 	 	 >> $2
		echo 'x_range_max= 15.15'				 	 >> $2
		echo 'y_range_min= -1.7'			 			 >> $2
		echo 'y_range_min= -3.0'			 			 >> $2
		echo 'y_range_max= 3.2'			 		 >> $2
		echo 'hline_ypos= 0.0'					 >> $2

		echo 'keyword_change1= plot_legend0,Gal-MD1000' 	>> $2
		echo 'keyword_change4= plot_legend3,Gal-SMD400' 	>> $2
		echo 'keyword_change7= plot_legend6,Ill-TNG300'  	>> $2
	endif

	echo 'use_loc_co_x= 0.22'		>> $2
	#echo 'use_loc_co_x= 0.32'		>> $2
	echo 'use_loc_co_y= 0.73'		>> $2
	#echo 'use_loc_co_y= 0.60'		>> $2

	if ($hod_plot_key == 'single') then
		echo '0= yes'		>> $21
		echo '1= no'		>> $21
		echo '2= no'		>> $21
		echo '3= yes'		>> $21
		echo '4= no'		>> $21
		echo '5= no'		>> $21
		echo '6= yes'		>> $21
		echo '7= no'		>> $21
		echo '8= no'		>> $21
		echo '9= sep'		>> $21
		echo '10= sep'		>> $21
		echo '11= sep'		>> $21
	else
		echo '0= yes'		>> $21
		echo '1= no'		>> $21
		echo '2= no'		>> $21
		echo '3= yes'		>> $21
		echo '4= no'		>> $21
		echo '5= no'		>> $21
		echo '6= yes'		>> $21
		echo '7= no'		>> $21
		echo '8= no'		>> $21
		echo '9= sep'		>> $21
		echo '10= sep'		>> $21
		echo '11= sep'		>> $21
	endif

	echo 'z_print_position_x= 0.25'		>> $2
	echo 'z_print_position_x= 0.75'		>> $2
	echo 'z_print_position_y= 0.85'		>> $2
	echo 'z_print_position_y= 0.9'		>> $2
	echo 'print_redshift= '$prefix'-density ' 			>> $2
	echo 'z_print_position_x2= 0.25'		>> $2
	echo 'z_print_position_y2= 0.80'		>> $2

	set halotyp = 'Rockstar'
	#set halotyp = 'FOF'
	#Gal-dens
	#echo 'print_redshift2= 15% sats, '$halotyp' halos' 			>> $2
	#Gal2-dens
	echo 'print_redshift2= 10% sats, '$halotyp' halos' 			>> $2
	#Gal2-dens ID in Gal-dens
	echo 'print_redshift2= 6.5% sats, '$halotyp' halos' 			>> $2
	#Gal2-dens ID from Gal-dens
	echo 'print_redshift2= 9% sats, '$halotyp' halos' 			>> $2
	#Gal2-dens mhalo
	echo 'print_redshift2= 2% sats, '$halotyp' halos' 			>> $2
	#Gal2-dens mstar
	echo 'print_redshift2= 8% sats, '$halotyp' halos' 			>> $2
	#Gal400-dens
	#echo 'print_redshift2= 14% sats, '$halotyp' halos' 			>> $2
	#TNG300 z=0.0
	echo 'print_redshift2= 8% sats' 			>> $2
	echo 'print_redshift2= CUT3-$sSFR$' 			>> $2
	echo 'print_redshift2= ' 			>> $2	
	echo 'legend_fontsize= 22'		 	 >> $2

	echo 'adjust_left= 0.18'					>> $2
	echo 'adjust_right= 0.99'					>> $2
	echo 'adjust_bottom= 0.14'					>> $2
	echo 'adjust_top= 0.99'					>> $2

	echo 'mylw= 6'				 		 >> $2
	echo 'lw_offset= -2'				 		 >> $2

	echo 'markersize= 14'				 	 >> $2
	echo 'mymarkeredgewidth= 3.0'				 	 >> $2

	echo 'text_fontsize= 28'					>> $2

	echo 'print_redshift_sep= ' >> $2
	echo 'legend_sep_xpos= 0.60'	>> $2
	echo 'legend_sep_ypos= 0.84'	>> $2
	echo 'legend_sep_offset= 0.035'	>> $2

	echo 'legend_sep_xpos= 0.60'	>> $2
	echo 'legend_sep_ypos= 0.16'	>> $2
	echo 'legend_sep_offset= 0.03'	>> $2
	
	echo 'plot_legend= yes'	 				 >> $2

else if ($1 == plotXY && $OII_ls_set == 'True') then
	#OH histos

	echo 'x_title= $\log_{10}$ (O/H)'	 	>> $2
	echo 'x_title= $\log_{10}$ (sfr [$M_{\odot}$ $yr^{-1}$])'	 	>> $2
	echo 'y_title= $\log_{10}$ (sfr inst [$M_{\odot}$ $yr^{-1}$])'	 	>> $2
	#echo 'y_title= $\log_{10}$ ($\Phi$  [$Mpc^{-3}$ $dex^{-1}$])'	>> $2
	echo 'keyword_change1= plot_legend0,O/H disk+bulge' 	>> $2
	echo 'keyword_change2= plot_legend1,(O/H disk+bulge)/16' 	>> $2
	echo 'keyword_change3= plot_legend2,(O/H disk+bulge)/16 + 0.25' 	>> $2
	echo 'keyword_change4= plot_legend3,apprx. gas' 	>> $2
	echo 'keyword_change5= plot_legend4,apprx. gas+star' 	>> $2

	echo 'keyword_change1= plot_legend0,sfr' 	>> $2
	echo 'keyword_change2= plot_legend1,sfr inst spheroid' 	>> $2
	echo 'keyword_change3= plot_legend2,sfr inst quiescent' 	>> $2
	echo 'keyword_change4= plot_legend3,sfr inst sph+quies' 	>> $2

	echo 'x_range_min= -8'			 		 >> $2
	echo 'x_range_max= 3.5'			 		 >> $2

	echo 'y_range_min= -8'		 			 	 >> $2
	echo 'y_range_max= 3.5'			 		 >> $2

	echo 'use_loc_co= yes'					>> $2
	echo 'use_loc_co_x= 0.1'			>> $2
	echo 'use_loc_co_y= 0.05'			>> $2
	echo 'print_redshift= SAGv3, z=1.0' 			>> $2

else if ($1 == plotXY && $plotXY_key == 'LF') then
	#CMASS-paper
	#SAGE 350h-1Mpc 1244

	set mag = 'M_RS'
	echo 'title=  '$3' '$4	 			 >> $2

	echo 'y_range_min= -7.0'		 			 	 >> $2
	echo 'y_range_max= -3.5'			 		 >> $2

	if ($mag == 'M_RS') then
		echo 'subplot_errorbars= False' 		 			 >> $2
		echo 'filled_between_subplot= True' 		 		 >> $2
		echo 'x_title= $M^{0.555}_i$'	 	>> $2
		echo 'x_range_min= -24.5'			 		 >> $2
		echo 'x_range_max= -21.3'			 		 >> $2

		echo 'use_loc_co_x= 0.45'			>> $2
		echo 'use_loc_co_y= 0.18'			>> $2
		echo 'z_print_position_x= 0.55'		>> $2
		echo 'z_print_position_y= 0.16'		>> $2


	else if ($mag == 'Mhalo') then
		echo 'x_range_min= 10.5'			 		 >> $2
		echo 'x_range_max= 15.5'			 		 >> $2
		echo 'y_range_min= -3.0'		 			 	 >> $2
		echo 'y_range_max= 3'			 		 >> $2
	
	else if ($mag == 'M') then
		echo 'x_title= $M_{AB_i}$'	 	>> $2
		echo 'x_range_min= -26'			 		 >> $2
		echo 'x_range_max= -20'			 		 >> $2

	else if ($mag == 'm') then
		echo 'x_title= $m_{AB}$'	 	>> $2
		echo 'x_range_min= 17'			 		 >> $2
		echo 'x_range_max= 24'			 		 >> $2

	endif 
	echo 'y_title= $\log_{10}$ ($\Phi$  [$Mpc^{-3}$ $mag^{-1}$])'	>> $2

	echo 'keyword_change1= plot_legend0,'$prefix'-cols' 	>> $2
	echo 'keyword_change2= plot_legend1,'$prefix'-dens'  	>> $2
	echo 'keyword_change3= plot_legend2,'$prefix'-mass' 	>> $2
			 		
	echo 'use_loc_co_x= 0.45'			>> $2
	echo 'use_loc_co_y= 0.16'			>> $2
	echo 'z_print_position_x= 0.16'		>> $2
	echo 'z_print_position_y= 0.85'		>> $2
	echo 'legend_fontsize= 20'				 	 >> $2

	echo 'error_bars_y= yes'					 >> $2
	echo 'error_bars_y_sub= yes'					 >> $2
	echo 'subplot_errorbars= False' 		 			 >> $2
	echo 'filled_between_subplot= True' 		 		 >> $2
	echo 'no_yticks= False'					>> $2

	#echo 'keyword_change1= plot_legend0,Gal-all centrals' 	>> $2
	#echo 'keyword_change2= plot_legend1,Gal-all sats' 	>> $2
	#echo 'keyword_change3= plot_legend2,Gal-cols centrals'  	>> $2
	#echo 'keyword_change4= plot_legend3,Gal-cols sats'	>> $2
	#echo 'keyword_change5= plot_legend4,Gal-dens centrals'	>> $2
	#echo 'keyword_change6= plot_legend5,Gal-dens sats'	>> $2


else if ($1 == plotXY && $plotXY_key == 'rhalf') then
	#CMASS mstar2rhalf		 			 >> $2 
	echo 'print_redshift= z=0.56' >> $2
	echo 'legend_fontsize= 20'				 	 >> $2
	echo 'x_range_min= 11'			 			 >> $2
	echo 'x_range_max= 12.3'			 		 	 >> $2
	echo 'y_range_min= 1'			 		 	 >> $2
	echo 'y_range_max= 4.2'			 		 >> $2
	set catname = 'Galacticus'
	echo 'minor_ticks_x_space= 0.05'			 	>> $2
	echo 'error_bars_y= yes'					 >> $2

	echo 'no_yticks= False'					>> $2
	echo 'use_loc_co_x= 0.40'				>> $2
	echo 'use_loc_co_y= 0.15'			>> $2
	echo 'x_title= $\log_{10}$ ($M_{_*}$ [$M_{\odot}$])'	         >> $2
	echo 'y_title= $\log_{10}$ ($r^{1/2}_{mass}$ [$kpc$])'       >> $2

	echo 'keyword_change1= plot_legend0,color    BOSS-CMASS cuts' 	>> $2
	echo 'keyword_change2= plot_legend1,density\t$n_z=0.89$x$10^{-4}$' 	>> $2
	echo 'keyword_change3= plot_legend2,full'  	>> $2
	echo 'keyword_change4= plot_legend3,full'	>> $2
	echo 'keyword_change5= plot_legend4,POR PS CMASS DR12'	>> $2
	echo 'keyword_change6= plot_legend5,POR SF CMASS DR12'	>> $2
	echo 'keyword_change7= plot_legend6,POR MER CMASS DR12'	>> $2

	echo 'z_print_position_x= 0.15'		>> $2
	echo 'z_print_position_y= 0.9'		>> $2

	echo 'add_axis= no'						 >> $2
	echo 'add_which_axis= y'					>> $2
	echo 'add_which_col= 7'						>> $2
	echo 'log_scale_add= no'						>> $2
	echo 'title_add0= $\log_{10}$ ($r_{bulge}$ [$kpc$]) color'		 	 >> $2
	echo 'title_add1= $\log_{10}$ ($r_{bulge}$ [$kpc$]) density'		 	 >> $2
	echo 'title_add2= $\log_{10}$ ($r_{bulge}$ [$kpc$]) full'		 	 >> $2
	echo 'add_range_min= min'					 >> $2
	echo 'add_range_max= max'					 >> $2
	echo 'add_range_min1= min'					 >> $2
	echo 'add_range_max1= max'					 >> $2
	echo 'add_range_min2= min'					 >> $2
	echo 'add_range_max2= max'					 >> $2

	echo 'minor_ticks_y_space= 0.1'			 	>> $2
	echo 'minor_ticks_y_set_right= on'				>> $2
	echo 'major_ticks_y_set_right= on'				>> $2

else if ($1 == plotXY && ($plotXY_key == 'SMF' || $plotXY_key == 'HMF' || $plotXY_key =~ *'SFR' || $plotXY_key == 'Mhalo' || $plotXY_key == 'Zcold' || $plotXY_key =~ *'-'* || $24 =~ *hist*)) then
	echo 'thin_line= 6'							>> $2
	echo 'set_thin_line= 99'				 	>> $2
	echo 'plot_legend= yes'		 			>> $2
	echo 'legend_fontsize= 28'				 	 >> $2
	set use_cb_colors = 9
	set sample = 'lo'
	if ($sample == 'low') then
		echo '0= '	>> $13						
		echo '1= '	>> $13
		echo '0= no'	>> $21						
		echo '1= no'	>> $21

	else
		echo '0= -'	>> $13						
		echo '1= --'	>> $13
		echo '0= yes'	>> $21						
		echo '1= yes'	>> $21
	endif


	#CMASS
	echo 'x_range_min= 10.8'			 			 >> $2
	echo 'x_range_max= 12.3'			 		 	 >> $2
	echo 'y_range_min= -7.2'			 		 	 >> $2
	echo 'y_range_max= -3.7'			 		 >> $2

	echo 'y_range_min= -7'			 		 	 >> $2
	echo 'y_range_max= -0.5'			 		 >> $2

	echo 'y_range_min= 0'			 		 	 >> $2
	echo 'y_range_max= 0.5'			 		 >> $2

	echo 'log_scale_y= no'						 >> $2
	echo 'log_scale_y_sub= no'					 >> $2
	#echo $plotXY_key
	#echo 'here echo $24: ' $24
	#LGALAXIES
	if ($plotXY_key == 'SMF' || $24 =~ *mstar*) then
		echo 'x_range_min= 9'			 			 >> $2
		echo 'x_range_max= 12'			 		 	 >> $2
	else if ($plotXY_key == 'SFR' || $24 =~ *_sfr*) then
		echo 'x_range_min= -3.5'			 		>> $2
		echo 'x_range_max= 2.5'			 		 	 >> $2
	else if ($plotXY_key == 'sSFR'|| $24 =~ *ssfr*) then
		echo 'x_range_min= -13'			 			 >> $2
		echo 'x_range_max= -7.5'			 		 	 >> $2
	else if ($plotXY_key == 'SHMF' || $24 =~ *SHMF*) then
		echo 'x_range_min= -4'			 		 	 >> $2
		echo 'x_range_max= -1'						>> $2
	else if ($plotXY_key == 'HMF' || $24 =~ *mhalo_200c*) then
		echo 'x_range_min= 9.9'			 			 >> $2
		echo 'x_range_max= 15.3'			 		 	 >> $2
	else if ($plotXY_key == 'Mhalo' || $24 =~ *mhalo_200c*) then
		echo 'x_range_min= 10.5'			 			 >> $2
		echo 'x_range_max= 14.5'			 		 	 >> $2
	else if ($plotXY_key == 'Mhalo' || $24 =~ *mhalo*) then
		echo 'x_range_min= 10.5'			 			 >> $2
		echo 'x_range_max= 14.7'			 		 	 >> $2
	else if ($plotXY_key == 'Mcold' || $24 =~ *mcold*) then
		echo 'x_range_min= 8'			 			 >> $2
		echo 'x_range_max= 11.5'			 		 	 >> $2
	else if ($plotXY_key == 'Mbh' || $24 =~ *mbh*) then
		echo 'x_range_min= 5'			 			 >> $2
		echo 'x_range_max= 9.3'			 		 	 >> $2
	else if ($plotXY_key == 'MZgas' || $24 =~ *Mzgas*) then
		echo 'x_range_min= 5'			 			 >> $2
		echo 'x_range_max= 10'			 		 	 >> $2
	else if ($plotXY_key == 'cgf' || $24 =~ *cgf*) then
		echo 'x_range_min= -5'			 			 >> $2
		echo 'x_range_max= 2'			 		 	 >> $2
	else if ($plotXY_key == 'Zcold' || $24 =~ *zcold*) then
		echo 'x_range_min= 8'			 			 >> $2
		echo 'x_range_max= 11'			 		 	 >> $2
	else if ($plotXY_key == 'g-i' || $24 =~ *g-i*) then
		echo 'x_range_min= 0'			 			 >> $2
		echo 'x_range_max= 3.7'			 		 	 >> $2
	else if ($plotXY_key == 'r-i' || $24 =~ *r-i*) then
		echo 'x_range_min= 0'			 			 >> $2
		echo 'x_range_max= 2'			 		 	 >> $2
	else if ($plotXY_key == 'rhalfmass' || $24 =~ *rhalfmass*) then
		echo 'x_range_min= -5'			 			 >> $2
		echo 'x_range_max= 1'			 		 	 >> $2
	else if ($plotXY_key == 'rbulgevsrdisk' || $24 =~ *rbulgevsrdisk*) then
		echo 'x_range_min= -1'			 			 >> $2
		echo 'x_range_max= 5'			 		 	 >> $2
	else if ($plotXY_key == 'vmax' || $24 =~ *vmax*) then
		echo 'x_range_min= 1.75'			 			 >> $2
		echo 'x_range_max= 3.25'			 		 	 >> $2
	else if ($plotXY_key == 'vdisp' || $24 =~ *vdisp*) then
		echo 'x_range_min= 1.75'			 			 >> $2
		echo 'x_range_max= 3.25'			 		 	 >> $2
	else if ($plotXY_key =~ 'mean_age'* || $24 =~ *mean_age*) then
		echo 'x_range_min= 0'			 			 >> $2
		echo 'x_range_max= 8'			 		 	 >> $2
	endif

	echo 'y_range_min= 0'			 		 	 >> $2
	echo 'y_range_max= 0.4'			 		 >> $2
	echo 'float_format_y= 2f'         				>> $2
	echo 'minor_ticks_x_space= 0.1'			 	>> $2
	echo 'minor_ticks_y_space= 0.025'			 	>> $2

	if ($plotXY_key == 'SMF'|| $24 =~ *mstar*) then
		echo 'x_title= $\log_{10}$ ($M_{Star}$ [$M_{\odot}$])'	         >> $2
		echo 'xticks= 9.5,10,10.5,11,11.5'		>> $2
	else if ($plotXY_key == 'SFR' || $24 =~ *_sfr*) then
		echo 'x_title= $\log_{10}$ ($SFR$ [$M_{\odot}$ $yr^{-1}$])'	         >> $2
		echo 'xticks= -3,-2,-1,0,1,2'		>> $2
	else if ($plotXY_key == 'sSFR'|| $24 =~ *ssfr*) then
		echo 'x_title= $\log_{10}$ ($sSFR$ [$yr^{-1}$])'	         >> $2
		echo 'xticks= -13,-12,-11,-10,-9,-8'		>> $2
	else if ($key == 'plotXY_SHMF' || $24 =~ *SHMF*) then
 		echo 'y_title= $M_{*}$/$M_{vir}$'  >> $2
		echo 'xticks= -4,-3,-2,-1'		>> $2
	else if ($plotXY_key == 'HMF' || $24 =~ *mhalo_200c*) then
		echo 'x_title= $\log_{10}$ ($M_{200c}$ [$M_{\odot}$])'	         >> $2
		echo 'xticks= 10,11,12,13,14,15'		>> $2
	else if ($plotXY_key == 'Mhalo' || $24 =~ *mhalo_200c*) then
		echo 'x_title= $\log_{10}$ ($M_{200c}$ [$M_{\odot}$])'	         >> $2
		echo 'xticks= 10.5,11,11.5,12,12.5,13,13.5,14,14.5'		>> $2
	else if ($plotXY_key == 'Mhalo' || $24 =~ *mhalo*) then
		echo 'x_title= $\log_{10}$ ($M_{Vir}$ [$M_{\odot}$])'	         >> $2
		echo 'xticks= 10.5,11,11.5,12,12.5,13,13.5,14,14.5'		>> $2
	else if ($plotXY_key == 'Mcold' || $24 =~ *mcold*) then
		echo 'x_title= $\log_{10}$ ($M_{Cold}$ [$M_{\odot}$])'	         >> $2
	else if ($plotXY_key == 'Mzgas' || $24 =~ *Mzgas*) then
		echo 'x_title= $\log_{10}$ ($M_{Z_{gas}}$ [$M_{\odot}$])'	         >> $2
	else if ($plotXY_key == 'cgf' || $24 =~ *cgf*) then
		echo 'x_title= $\log_{10}$ ($M_{Cold}/M_*$)'	         >> $2
	else if ($plotXY_key == 'Mbh' || $24 =~ *mbh*) then
		echo 'x_title= $\log_{10}$ ($M_{BH}$ [$M_{\odot}$])'	         >> $2
		echo 'xticks= 5,6,7,8,9'		>> $2
	else if ($plotXY_key == 'MZgas' || $24 =~ *Mzgas*) then
		echo 'x_title= $\log_{10}$ ($M_{Cold}$ [$M_{\odot}$])'	         >> $2
		#echo 'xticks= 11,11.5,12,12.5,13,13.5,14'		>> $2
	else if ($plotXY_key == 'Zcold' || $24 =~ *zcold*) then
		echo 'x_title= $Z_{cold}$'	         >> $2
		echo 'xticks= 8,9,10,11'		>> $2
	else if ($plotXY_key == 'g-i' || $24 =~ *g-i*) then
		echo 'x_title= $g-i$'	         >> $2
		echo 'xticks= 0,0.5,1,1.5,2,2.5,3,3.5'		>> $2
	else if ($plotXY_key == 'r-i' || $24 =~ *r-i*) then
		echo 'x_title= $r-i$'	         >> $2
		echo 'xticks= 0,0.5,1,1.5,2'		>> $2
	else if ($plotXY_key == 'vmax' || $24 =~ *vmax*) then
		echo 'x_title= $\log_{10}$ ($V_{max}$ [$kms^-1$])'	         >> $2
		echo 'xticks= 1.5,2,2.5,3'		>> $2
	else if ($plotXY_key == 'vdisp' || $24 =~ *vdisp*) then
		echo 'x_title=  $\log_{10}$ ($V_{disp}$ [$kms^-1$])'	         >> $2
		echo 'xticks= 1.5,2,2.5,3'		>> $2
	else if ($plotXY_key == 'rbulgevsrdisk' || $24 =~ *rbulgevsrdisk*) then
		echo 'x_title=  $\log_{10}$ ($r_{bulge}/r_{disk}$)'	         >> $2
		echo 'xticks= -1,0,1,2,3,4,5'		>> $2
	else if ($plotXY_key == 'rhalfmass' || $24 =~ *rhalfmass*) then
		echo 'x_title=  $\log_{10}$ ($r_{1/2}$ [$Mpc$])'	         >> $2
		echo 'xticks= -5,-4,-3,-2,-1,0,1'		>> $2
	else if ($plotXY_key == 'mean_age_stars_disk' || $24 =~ *mean_age_stars_disk*) then
		echo 'x_title= $age_{disk}$ [$Gyr$]'	         >> $2
		echo 'xticks= 0,1,2,3,4,5,6,7,8'		>> $2
	else if ($plotXY_key == 'mean_age_stars_spheroid' || $24 =~ *mean_age_stars_spheroid*) then
		echo 'x_title= $age_{bulge}$ [$Gyr$]'	         >> $2
		echo 'xticks= 0,1,2,3,4,5,6,7,8'		>> $2
	endif

	echo 'error_bars_y= no'					 >> $2
	echo 'error_bars_y_sub= yes'					 >> $2
	echo 'no_yticks= False'					>> $2
	echo 'no_last_yticks= False'					>> $2
	echo 'xticks_minor_offset= 0'				>> $2

	echo 'y_title= $\log_{10}$ ($\Phi$/Mpc$^{-3}$ dlog$_{10}$M$_{\star}$)  '       >> $2
	echo 'y_title= $\log_{10}$ ($\Phi$ [$Mpc^{-3}$ $dex^{-1}$])  '       >> $2
	echo 'y_title= f'       >> $2
	#echo 'y_title= N'		 >> $2
	#echo 'y_title= '		 >> $2
	#echo 'no_yticks= True'    	 >> $2

	#echo 'x_title= '		 >> $2
	#echo 'no_xticks= True'    	 >> $2

	echo 'z_print_position_x= 0.14'		>> $2
	echo 'z_print_position_y= 0.86'		>> $2
	echo 'use_loc_co_x= 0.35'		>> $2
	#8 cats = 0.80
	echo 'use_loc_co_y= 0.85'		>> $2

	if ($prefix == 'Gal') then
		echo 'keyword_change1= plot_legend0,'$prefix'-cols' 	>> $2
		echo 'keyword_change2= plot_legend1,'$prefix'-all'  	>> $2
		echo 'keyword_change3= plot_legend2,'$prefix'-dens'  	>> $2
		echo 'keyword_change4= plot_legend3,'$prefix'-mass' 	>> $2
	else if ($prefix == 'sLG') then
		echo 'keyword_change1= plot_legend0,'$prefix'-down' 	>> $2
		echo 'keyword_change2= plot_legend1,'$prefix'-all'  	>> $2
		echo 'keyword_change3= plot_legend2,'$prefix'-dens $M_*$'  	>> $2
		#echo 'keyword_change3= plot_legend2,'$prefix'-down w/o orph-sats'  	>> $2
		echo 'keyword_change4= plot_legend3,'$prefix'-dens $V_{max}$' 	>> $2
		#echo 'keyword_change4= plot_legend3,'$prefix'-all w/o orph-sats' 	>> $2
		#echo 'keyword_change4= plot_legend3,'$prefix'-dens $M_{Halo}$ cents'  	>> $2
	else if ($prefix == 'z=') then

		#echo 'y_range_min= 10'			 		 	 >> $2
		#echo 'y_range_max= 1e4'			 		 >> $2



		echo 'keyword_change1= plot_legend0,'$prefix'0.56' 	>> $2
		echo 'keyword_change2= plot_legend1,'$prefix'0.75'  	>> $2
		echo 'keyword_change3= plot_legend2,'$prefix'1'  	>> $2
		echo 'keyword_change4= plot_legend3,'$prefix'1.5' 	>> $2
		echo 'keyword_change5= plot_legend4,'$prefix'2' 	>> $2
		echo 'keyword_change6= plot_legend5,'$prefix'2.5' 	>> $2
		echo 'keyword_change7= plot_legend6,'$prefix'3' 	>> $2
		echo 'keyword_change8= plot_legend7,'$prefix'3.5' 	>> $2
		echo 'keyword_change9= plot_legend8,'$prefix'4' 	>> $2

		if ($3 =~ 400*) then
			if ($24 =~ *isto1) then
				echo 'keyword_change1= plot_legend0,'$prefix'0.55' 	>> $2
				echo 'keyword_change2= plot_legend1,'$prefix'0.74' 	>> $2
				echo 'keyword_change3= plot_legend2,'$prefix'0.84' 	>> $2
				echo 'keyword_change4= plot_legend3,'$prefix'0.89' 	>> $2
				echo 'keyword_change5= plot_legend4,'$prefix'1.0' 	>> $2
				echo 'keyword_change6= plot_legend5,'$prefix'1.26' 	>> $2
			else if ($24 =~ *isto2) then
				echo 'keyword_change1= plot_legend0,'$prefix'1.36' 	>> $2
				echo 'keyword_change2= plot_legend1,'$prefix'1.56' 	>> $2
				echo 'keyword_change3= plot_legend2,'$prefix'2.02' 	>> $2
				echo 'keyword_change4= plot_legend3,'$prefix'2.56' 	>> $2
				echo 'keyword_change5= plot_legend4,'$prefix'3.03' 	>> $2
				echo 'keyword_change6= plot_legend5,'$prefix'4.15' 	>> $2
			else
				echo 'keyword_change1= plot_legend0,'$prefix'0.55' 	>> $2
				echo 'keyword_change2= plot_legend1,'$prefix'0.74' 	>> $2
				echo 'keyword_change3= plot_legend2,'$prefix'1.0' 	>> $2
				echo 'keyword_change4= plot_legend3,'$prefix'2.02' 	>> $2
				echo 'keyword_change5= plot_legend4,'$prefix'3.03' 	>> $2
				echo 'keyword_change6= plot_legend5,'$prefix'4.15' 	>> $2
			endif
		else
			if ($24 =~ *isto1) then
				echo 'keyword_change1= plot_legend0,'$prefix'0.56' 	>> $2
				echo 'keyword_change2= plot_legend1,'$prefix'0.74' 	>> $2
				echo 'keyword_change3= plot_legend2,'$prefix'0.82' 	>> $2
				echo 'keyword_change4= plot_legend3,'$prefix'0.9' 	>> $2
				echo 'keyword_change5= plot_legend4,'$prefix'0.99' 	>> $2
				echo 'keyword_change6= plot_legend5,'$prefix'1.27' 	>> $2
			else if ($24 =~ *isto2) then
				echo 'keyword_change1= plot_legend0,'$prefix'1.37' 	>> $2
				echo 'keyword_change2= plot_legend1,'$prefix'1.54' 	>> $2
				echo 'keyword_change3= plot_legend2,'$prefix'2.03' 	>> $2
				echo 'keyword_change4= plot_legend3,'$prefix'2.53' 	>> $2
				echo 'keyword_change5= plot_legend4,'$prefix'3.04' 	>> $2
				echo 'keyword_change6= plot_legend5,'$prefix'4.15' 	>> $2
			else
				echo 'keyword_change1= plot_legend0,'$prefix'0.56' 	>> $2
				echo 'keyword_change2= plot_legend1,'$prefix'0.74' 	>> $2
				echo 'keyword_change3= plot_legend2,'$prefix'0.99' 	>> $2
				echo 'keyword_change4= plot_legend3,'$prefix'2.03' 	>> $2
				echo 'keyword_change5= plot_legend4,'$prefix'3.04' 	>> $2
				echo 'keyword_change6= plot_legend5,'$prefix'4.15' 	>> $2
			endif
		endif


	else
		echo 'keyword_change1= plot_legend0,'$prefix'-down' 	>> $2
		echo 'keyword_change2= plot_legend1,'$prefix'-dens $V_{max}$'  	>> $2
		#echo 'keyword_change3= plot_legend2,'$prefix'-dens $SFR$'  	>> $2
		#echo 'keyword_change4= plot_legend3,'$prefix'-dens $sSFR$' 	>> $2
		echo 'keyword_change3= plot_legend2,'$prefix'-dens $M_*$' 	>> $2
		#echo 'keyword_change6= plot_legend5,'$prefix'-dens $V_{max}$ sats' 	>> $2

		echo 'keyword_change1= plot_legend0,sfr>1e-4' 	>> $2
		echo 'keyword_change2= plot_legend1,MD-paper'  	>> $2
		echo 'keyword_change3= plot_legend2,test'  	>> $2

	endif

	#echo 'print_redshift= z=0.56'		>> $2

	#echo 'keyword_change4= plot_legend3,Gal-sample2 centrals' 	>> $2
	#echo 'keyword_change5= plot_legend4,Gal-all centrals' 	>> $2
	#echo 'keyword_change6= plot_legend5,Gal-sample2 orphans' 	>> $2
	#echo 'keyword_change7= plot_legend6,Gal-all orphans' 	>> $2
	#echo 'keyword_change8= plot_legend7,Gal-sample2 no-sats' 	>> $2
	#echo 'keyword_change9= plot_legend8,Gal-all no-sats' 	>> $2

	#echo 'keyword_change1= plot_legend0,all' 	>> $2
	#echo 'keyword_change2= plot_legend1,g-i>2.35'  	>> $2
	#echo 'keyword_change3= plot_legend2,g-i>2.15'  	>> $2
	#echo 'keyword_change4= plot_legend3,g-i>2.0' 	>> $2

	#echo 'keyword_change1= plot_legend0,Gal-cols' 	>> $2
	#echo 'keyword_change2= plot_legend1,down'  	>> $2
	#echo 'keyword_change3= plot_legend2,down2'  	>> $2
	#echo 'keyword_change4= plot_legend3,down3' 	>> $2

	#OII
	set key_OII = 'False'
	if ($key_OII == 'True') then
		echo 'x_range_min= 9'			 			 >> $2
		echo 'x_range_max= 13'			 		 	 >> $2
		echo 'y_range_min= -7.5'			 		 	 >> $2
		echo 'y_range_max= -1'			 		 >> $2
		echo 'float_format_y= 0f'         				>> $2
		echo 'adjust_right= 0.96'					>> $2
		echo 'adjust_top= 0.97'					>> $2
		echo 'adjust_left= 0.12'					>> $2
		echo 'adjust_bottom= 0.16'					>> $2
		echo 'size_x= 14'				 		 >> $2
		echo 'size_y= 7.5'				 		 >> $2
	endif

	echo 'adjust_right= 0.98'					>> $2
	echo 'adjust_top= 0.83'					>> $2
	echo 'adjust_left= 0.14'					>> $2
	echo 'adjust_bottom= 0.14'					>> $2
	echo 'size_x= 13'				 		 >> $2
	echo 'size_y= 9'				 		 >> $2

	#echo 'keyword_change1= plot_legend0,MDPL2' 	>> $2
	#echo 'keyword_change2= plot_legend1,SMDPL'  	>> $2
	#echo 'keyword_change3= plot_legend2,TNG300-1'  	>> $2
	#echo 'keyword_change4= plot_legend3,Gal1000' 	>> $2
	#echo 'keyword_change5= plot_legend4,Gal400'  	>> $2
	#echo 'keyword_change6= plot_legend5,Illustris'  	>> $2

	#echo 'float_format_y= 1f'         				>> $2
	#echo 'minor_ticks_x_space= 0.1'			 	>> $2
	#echo 'minor_ticks_y_space= 0.25'			 	>> $2
	#echo 'size_x= 12'				 		 >> $2
	#echo 'size_y= 10'				 		 >> $2
	#echo 'adjust_right= 0.99'					>> $2
	#echo 'adjust_top= 0.99'					>> $2
	#echo 'adjust_left= 0.14'					>> $2
	#echo 'adjust_bottom= 0.12'					>> $2

	#echo 'no_last_yticks= True'					>> $2

	#echo 'z_print_position_x= 0.8'		>> $2
	#echo 'z_print_position_y= 0.16'		>> $2
	#echo 'use_loc_co_x= 0.45'		>> $2
	#echo 'use_loc_co_y= 0.74'		>> $2
	#set use_cb_colors = 3

############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
else if ($1 == plotXY && $plotXY_key == 'sfr2z') then
	set no_x = noxd
	set no_y = noyd
	set share_y = 'no'
	set no_top = False
	set no_right = False
	set xIsz = 'True'
	echo 'legend_fontsize= 28'				 	 >> $2
	echo 'text_fontsize= 35'				 	 >> $2
	echo 'plot_legend= no'			>> $2
	echo 'print_redshift= True'				 	 >> $2
	echo 'no_xticks= False'        		>> $2
	#echo 'no_first_xticks= True'        		>> $2
	#echo 'no_last_yticks= True'        		>> $2
	set small_plot = 'True'
	set reduced = 'True'
	set demo = 'False'
	#choose which sample to plot: 1=M1/M2, 2a=M2 additional samples/original, 2b=M2 additional samples / reduced for paper, 300=cluster, 'more', 'pop'
	set sample = 1b
	set key = $24

	set yerr = 1
	echo 'fill_no_facecolor= False' 				>> $2

	#echo 'vline_xpos= 0.7,1.38,2.12,3.54'					 >> $2

	#echo 'print_redshift2= (b)'				>> $2
	echo 'z_print_position_x2= 0.18'		>> $2
	echo 'z_print_position_y2= 0.90'		>> $2

	if ( $small_plot == 'True') then
		echo 'size_x= 13'				 		 >> $2
		echo 'size_y= 4.5'				 		 >> $2
		echo 'adjust_left= 0.16'					>> $2
		echo 'adjust_right= 0.99'					>> $2
		echo 'adjust_bottom= 0.20'					>> $2
		echo 'adjust_top= 0.70'					>> $2

		echo 'z_print_position_x= 0.85'		>> $2
		echo 'z_print_position_y= 0.82'		>> $2

		echo 'use_loc_co_x= 0.45'		>> $2
		echo 'use_loc_co_y= 0.55'		>> $2

	else if ($small_plot == 'res') then
		echo 'size_x= 4.5'				 		 >> $2
		echo 'size_y= 3.5'				 		 >> $2
		echo 'adjust_left= 0.02'					>> $2
		echo 'adjust_right= 0.98'					>> $2
		echo 'adjust_bottom= 0.29'					>> $2
		echo 'adjust_top= 0.98'					>> $2

		echo 'z_print_position_x= 0.04'		>> $2
		echo 'z_print_position_y= 0.83'		>> $2


		echo 'x_range_min= 12.7'			 		 	 >> $2
		echo 'x_range_max= 14.1'			 			 >> $2

		echo 'log_scale_x= no'						 >> $2
		echo 'log_scale_x_sub= no'					 >> $2

		echo 'x_title= $M_{vir}$ [$M_{\odot}$]'  >> $2

		#mhalo high-Zcold, 1) filaments, 2)knots
		echo 'vline_xpos= 12.92,13.09,13.26,12.82,12.97,13.12'					 >> $2
		#echo 'z_print_position_x= 0.54'		>> $2

		#mhalo low-Zcold, 1) filaments, 2)knots
		echo 'vline_xpos= 13.69,13.85,14.01,13.42,13.55,13.70'					 >> $2

		echo 'xticks= 13,13.5,14'		>> $2
		echo 'minor_ticks_x_space= 0.1'			 	>> $2
		echo 'no_yticks= True' >> $2
		#echo 'no_xticks= True' >> $2
		#echo 'x_title= '  >> $2
		echo 'text_fontsize= 28'				 	 >> $2

	else if ($demo == 'True') then
		echo 'size_x= 12'				 		 >> $2
		echo 'size_y= 8'				 		 >> $2
		echo 'adjust_left= 0.16'					>> $2
		echo 'adjust_right= 0.99'					>> $2
		echo 'adjust_bottom= 0.13'					>> $2
		echo 'adjust_top= 0.97'					>> $2

		#echo 'z_print_position_x= 0.51'		>> $2
		#echo 'z_print_position_y= 0.68'		>> $2

		#echo 'z_print_position_x= 0.2'		>> $2
		#echo 'z_print_position_y= 0.15'		>> $2

		echo 'use_loc_co_x= 0.60'		>> $2
		echo 'use_loc_co_y= 0.63'		>> $2
	else
		echo 'size_x= 12'				 		 >> $2
		echo 'size_y= 8'				 		 >> $2
		echo 'adjust_left= 0.16'					>> $2
		echo 'adjust_right= 0.99'					>> $2
		echo 'adjust_bottom= 0.14'					>> $2
		echo 'adjust_top= 0.99'					>> $2

		echo 'z_print_position_x= 0.17'		>> $2
		echo 'z_print_position_y= 0.75'		>> $2

		#echo 'z_print_position_x= 0.2'		>> $2
		#echo 'z_print_position_y= 0.15'		>> $2

		echo 'use_loc_co_x= 0.25'		>> $2
		echo 'use_loc_co_y= 0.87'		>> $2
	endif

	echo 'ylabel_pos= 0.55'	         			>> $2    

	if ( $xIsz == 'True') then
		echo 'x_title= z'				 			 >> $2 
		echo 'log_scale_x= yes'						 >> $2
		echo 'log_scale_x_sub= yes'					 >> $2

		echo 'minor_ticks_x_space= no'			 	>> $2
		echo 'minor_ticks_y_space= no'			 	>> $2
		echo 'xticks_minor= 0'					>> $2

		if ($sample == 3) then
			echo 'x_range_min= 0.093'			 			 >> $2
			echo 'x_range_max= 0.5542'			 			 >> $2

			if ($no_x == 'nox') then
				echo 'x_title= ' >> $2
				echo 'no_xticks= True'        					>> $2
				echo 'xticks= 0.1,0.12,0.15,0.2,0.25,0.3,0.4,0.5'		>> $2
				#echo 'z_print_position_x= 0.005'		>> $2
				#echo 'z_print_position_y= 0.72'		>> $2
			else

				echo 'xticks= 0.1,0.12,0.15,0.2,0.25,0.3,0.4,0.5'		>> $2
			endif
		else if ($sample == 300) then
			echo 'log_scale_x= symlog'						 >> $2
			echo 'x_range_min= 0'			 		 	 >> $2 
			echo 'x_range_max= 1.5'			 			 >> $2

			echo 'xticks_minor= 0'					>> $2
			echo 'minor_ticks_x_space= no'			 	>> $2
			echo 'xticks= 0,0.1,0.5,0.75,1,1.5,2,3,4'		>> $2
			echo 'xticks= 0,0.1,0.25,0.5,0.75,1.0'		>> $2
			echo 'no_first_xticks= True'		 			 >> $2

			if ($no_x == 'nox') then
				echo 'x_title= ' >> $2
				echo 'no_xticks= True'        					>> $2			
			endif
		else
			echo 'x_range_min= 0.5574'			 			 >> $2
			#echo 'x_range_min= 0.5924'			 			 >> $2
			echo 'x_range_min= 0.55'			 			 >> $2
			echo 'x_range_max= 4.15'			 			 >> $2
			echo 'xticks= 0.55,0.7,0.85,1,1.2,1.5,2,3,4'		>> $2
			if ($no_x == 'nox') then
				echo 'x_title= ' >> $2
				echo 'no_xticks= True'        					>> $2			
			endif
		endif

	else if ($xIsz == 'Gyr') then
		echo 'x_title= lookback time [$Gyr$]' 			 >> $2 
		echo 'log_scale_x= no'						 >> $2
		echo 'log_scale_x_sub= no'					 >> $2

		echo 'x_range_min= 0'			 			 >> $2
		echo 'x_range_max= 6.74389'			 			 >> $2

		echo 'minor_ticks_x_space= 0.25'			 	>> $2
		echo 'minor_ticks_y_space= no'			 	>> $2
		echo 'xticks_minor= 0'					>> $2
	endif

	echo 'xticks_minor_offset= 0'				>> $2

	echo 'float_format_y= 3f'         				>> $2 
	#echo 'no_last_yticks= True'		 			 >> $2

	if ($yerr == '1') then
		echo 'error_bars_y= yes'					 >> $2
		echo 'error_bars_y_sub= yes'					 >> $2
	endif
	echo 'subplot_errorbars= True' 		 			 >> $2
	echo 'filled_between_subplot= True' 		 		 >> $2

	if ( $xIsz == 'True' ) then
		set add_x = 'yess'
	else
		set add_x = 'no'
	endif

	if ($key =~ plotXY_msta*) then
		echo 'y_title=  $M_*$ [$M_{\odot}$]'  >> $2	
		echo 'y_range_min= 4e9'			 	>> $2
		#echo 'y_range_min= 9e10'			 	>> $2  
		echo 'y_range_max= 3e11'			>> $2

		if ($sample == 300) then
				echo 'y_range_min= 4e9'			 	>> $2
				echo 'y_range_max= 1e13'			 >> $2
		endif

		echo 'log_scale_y= yes'						 >> $2
		echo 'log_scale_y_sub= yes'					 >> $2

	else if ($key == 'plotXY_sfr') then
		echo 'y_title=  SFR [$M_{\odot}$ $yr^{-1}$]'  >> $2
		echo 'y_range_min= 0.01'			 		 	 >> $2 
		echo 'y_range_max= 200'			 			 >> $2

		echo 'y_range_min= 0'			 		 	 >> $2 
		echo 'y_range_max= 25'			 			 >> $2
		echo 'minor_ticks_y_space= 1'			>> $2
		if ($sample == 300) then
				echo 'y_range_min= 0.001'			 		 	 >> $2
				echo 'y_range_max= 200'			 			 >> $2
		endif
		echo 'z_print_position_x= 0.28'		>> $2		
		echo 'z_print_position_y= 0.75'		>> $2
		echo 'log_scale_y= no'						 >> $2
		echo 'log_scale_y_sub= no'					 >> $2

	else if ($key == 'plotXYsumSFR') then
		echo 'y_title= $\sum$ SFR [$M_{\odot}$ $yr^{-1}$ $Mpc^{-3}$]'  >> $2
		echo 'y_range_min= 0.03'			 		 	 >> $2 
		echo 'y_range_max= 90'			 			 >> $2
		echo 'y_range_min= 1e-6'			 		 	 >> $2 
		echo 'y_range_max= 0.005'			 			 >> $2

		echo 'error_bars_y= no'					 >> $2
		echo 'error_bars_y_sub= no'					 >> $2

		echo 'log_scale_y= yes'						 >> $2
		echo 'log_scale_y_sub= yes'					 >> $2
		echo 'z_print_position_x= 0.20'		>> $2
		echo 'z_print_position_y= 0.75'		>> $2

	else if ($key == 'plotXY_ssfr') then
		echo 'y_title=  sSFR [$yr^{-1}$]'  >> $2
		echo 'y_range_min= 3e-13'			 		 	 >> $2
		echo 'y_range_max= 5e-9'			 			 >> $2
		echo #'y_range_max= 2e-11'			 			 >> $2
		echo 'log_scale_y= yes'						 >> $2
		echo 'log_scale_y_sub= yes'					 >> $2

		if ($sample == 300) then
				echo 'y_range_min= 1e-15'			 		 	 >> $2
				echo 'y_range_max= 1e-8'			 			 >> $2
		endif
		echo 'z_print_position_x= 0.28'		>> $2
		echo 'z_print_position_y= 0.75'		>> $2

	else if ($key == 'plotXY_mhalo') then
		echo 'y_title= $M_{vir}$ [$M_{\odot}$]'  >> $2
		echo 'y_range_min= 3e11'			 		 	 >> $2
		#echo 'y_range_min= 1e13'			 		 	 >> $2
		echo 'y_range_max= 1.5e14'			 			 >> $2
		echo 'log_scale_y= yes'						 >> $2
		echo 'log_scale_y_sub= yes'					 >> $2
		if ($sample == 300) then
			echo 'y_range_min= 5e10'			 		 	 >> $2
			echo 'y_range_max= 8e15'			 			 >> $2
		endif

	else if ($key == 'plotXY_mhalo_200c') then
		echo 'y_title= $M_{200c}$ [$M_{\odot}$]'  >> $2
		echo 'y_range_min= 3e11'			 		 	 >> $2
		echo 'y_range_max= 5e13'			 			 >> $2

		echo 'log_scale_y= yes'						 >> $2
		echo 'log_scale_y_sub= yes'					 >> $2

	else if ($key == 'plotXY_SHMF') then
 		echo 'y_title= $M_{*}$/$M_{vir}$'  >> $2
		echo 'y_range_min= 0.001'			 		 	 >> $2
		echo 'y_range_max= 0.02'			 			 >> $2

		#echo 'y_range_min= 0.01'			 		 	 >> $2
		#echo 'y_range_max= 0.1'	ssfr		 			 >> $2

		echo 'log_scale_y= yes'						 >> $2
		echo 'log_scale_y_sub= yes'					 >> $2
		echo 'z_print_position_x= 0.70'		>> $2
		echo 'z_print_position_y= 0.15'		>> $2

	else if ($key == 'plotXY_SHMF_200c') then
 		echo 'y_title= $M_{*}$/$M_{200c}$'  >> $2
		echo 'y_range_min= 0.001'			 		 	 >> $2
		echo 'y_range_max= 0.03'			 			 >> $2

		echo 'log_scale_y= yes'						 >> $2
		echo 'log_scale_y_sub= yes'					 >> $2

	else if ($key == 'plotXY_logSHMF') then
 		echo 'y_title= $\log_{10}$ $M_{*}$/$M_{vir}$'  >> $2
		echo 'y_range_min= -2.5'			 		 	 >> $2
		echo 'y_range_max= -1.6'			 			 >> $2

		echo 'minor_ticks_y_space= 0.05'			 	>> $2
		echo 'log_scale_y= no'						 >> $2
		echo 'log_scale_y_sub= no'					 >> $2
		echo 'no_last_yticks= False'		 			 >> $2

		echo 'z_print_position_x= 0.70'		>> $2
		echo 'z_print_position_y= 0.15'		>> $2

	else if ($key == 'plotXY_logSHMF_200c') then
 		echo 'y_title= $\log_{10}$ $M_{*}$/$M_{200c}$'  >> $2
		echo 'y_range_min= -2.2'			 		 	 >> $2
		echo 'y_range_max= -1.6'			 			 >> $2

		echo 'minor_ticks_y_space= 0.05'			 	>> $2
		echo 'log_scale_y= no'						 >> $2
		echo 'log_scale_y_sub= no'					 >> $2
		echo 'no_last_yticks= False'		 			 >> $2

 	else if ($key == 'mstar_z0') then
		echo 'y_title=  $M_*$/$M_*$$_{z=0.56}$  '  >> $2

	else if ($key == 'plotXY_zcold') then
		echo 'y_title= $Z_{cold}$'  >> $2
		echo 'y_range_min= 8.7'			 		 	 >> $2
		echo 'y_range_max= 10.5'			 			 >> $2

		if ($sample == 300) then
			echo 'y_range_max= 9.5'			 			 >> $2
		endif

		echo 'log_scale_y= no'						 >> $2
		echo 'log_scale_y_sub= no'					 >> $2
		echo 'minor_ticks_y_space= 0.05'			>> $2

		echo 'z_print_position_x= 0.28'		>> $2		
		echo 'z_print_position_y= 0.75'		>> $2
		#echo 'z_print_position_y= 0.12'		>> $2

		echo 'float_format_y= 1f'         >> $2

	else if ($key == 'plotXY_mcold') then
		echo 'y_title=  $M_{cold}$ [$M_{\odot}$]'  >> $2	
		echo 'y_range_min= 5e6'			 	>> $2
		#echo 'y_range_min= 1e9'			 	>> $2 
		echo 'y_range_max= 8e10'			>> $2

		if ($sample == 300) then
				echo 'y_range_min= 1e9'			 		 	 >> $2
				echo 'y_range_max= 1e12'			 			 >> $2
		endif

		echo 'log_scale_y= yes'						 >> $2
		echo 'log_scale_y_sub= yes'					 >> $2
		echo 'z_print_position_x= 0.70'		>> $2
		echo 'z_print_position_y= 0.75'		>> $2

	else if ($key == 'plotXY_Mzgas') then
		echo 'y_title=  $M_{z_{cold}}$ [$M_{\odot}$]'  >> $2	
		echo 'y_range_min= 1e6'			 	>> $2
		#echo 'y_range_min= 1e8'			 	>> $2  
		echo 'y_range_max= 8e9'			>> $2

		if ($sample == 300) then
				echo 'y_range_min= 1e7'			 		 	 >> $2
				echo 'y_range_max= 5e10'			 			 >> $2
		endif
		echo 'z_print_position_x= 0.70'		>> $2
		echo 'z_print_position_y= 0.75'		>> $2
		echo 'log_scale_y= yes'						 >> $2
		echo 'log_scale_y_sub= yes'					 >> $2

	else if ($key == 'plotXY_Tcons') then
		echo 'y_title=  $M_{cold}/SFR$ [$yr$]'  >> $2	
		echo 'y_range_min= 5e7'			 	>> $2
		#echo 'y_range_min= 1e9'			 	>> $2 
		echo 'y_range_max= 8e10'			>> $2

		if ($sample == 300) then
				echo 'y_range_min= 1e9'			 		 	 >> $2
				echo 'y_range_max= 1e12'			 			 >> $2
		endif

		echo 'log_scale_y= yes'						 >> $2
		echo 'log_scale_y_sub= yes'					 >> $2
		echo 'z_print_position_x= 0.75'		>> $2
		echo 'z_print_position_y= 0.80'		>> $2

	else if ($key == 'plotXY_mbh') then
		echo 'y_title=  $M_{BH}$ [$M_{\odot}$]'  >> $2	
		echo 'y_range_min= 1e4'			 	>> $2
		echo 'y_range_max= 8e8'			>> $2

		if ($sample == 300) then
				echo 'y_range_min= 1e6'			 		 	 >> $2
				echo 'y_range_max= 5e10'			 			 >> $2
		endif

		echo 'log_scale_y= yes'						 >> $2
		echo 'log_scale_y_sub= yes'					 >> $2

	else if ($key == 'plotXY_cgf') then
		echo 'y_title=  $M_{cold}$/$M_{*}$'  >> $2
		echo 'y_range_min= 0.00005'			 		 	 >> $2
		echo 'y_range_max= 0.05'			 			 >> $2
		echo 'log_scale_y= yes'						 >> $2
		echo 'log_scale_y_sub= yes'					 >> $2
		echo 'z_print_position_x= 0.2'		>> $2
		echo 'z_print_position_y= 0.75'		>> $2

	else if ($key == 'plotXY_g-i') then
		echo 'y_title= $g-i$'  >> $2
		echo 'y_range_min= -0.1'			 		 	 >> $2
		#echo 'y_range_min= 1'			 		 	 >> $2
		echo 'y_range_max= 3.5'			 			 >> $2
		echo 'minor_ticks_y_space= 0.1'			>> $2
		echo 'log_scale_y= no'						 >> $2
		echo 'log_scale_y_sub= no'					 >> $2
		echo 'float_format_y= 1f'         >> $2

	else if ($key == 'plotXY_r-i') then
		echo 'y_title= $r-i$'  >> $2
		echo 'y_range_min= -0.1'			 		 	 >> $2
		echo 'y_range_max= 1.9'			 			 >> $2

		echo 'minor_ticks_y_space= 0.05'			>> $2
		echo 'log_scale_y= no'						 >> $2
		echo 'log_scale_y_sub= no'					 >> $2

		echo 'z_print_position_x= 0.28'		>> $2		
		echo 'z_print_position_y= 0.75'		>> $2


	else if ($key == 'plotXY_rhalfmass') then

		echo 'y_title=  $r_{1/2}$ [$Mpc$]'  >> $2
		echo 'y_range_min= 3e-4'			 	>> $2
		echo 'y_range_max= 3'			>> $2

		if ($sample == 300) then
				echo 'y_range_min= 1e-3'			 		 	 >> $2
				echo 'y_range_max= 50'			 			 >> $2
		endif
		echo 'z_print_position_x= 0.70'		>> $2
		echo 'z_print_position_y= 0.75'		>> $2
		echo 'log_scale_y= yes'						 >> $2
		echo 'log_scale_y_sub= yes'					 >> $2

		echo 'y_range_min= 0'			 	>> $2
		echo 'y_range_max= 115'			>> $2
		echo 'log_scale_y= no'						 >> $2
		echo 'log_scale_y_sub= no'					 >> $2
		echo 'minor_ticks_y_space= 5'			>> $2

	else if ($key =~ plotXY_r*) then

		if ($key == 'plotXY_rdisk') then
			echo 'y_title=  $r_{disk}$ [$kpc$]'  >> $2
		else
			echo 'y_title=  $r_{bulge}$ [$kpc$]'  >> $2
		endif
		echo 'y_range_min= 0'			 	>> $2
		echo 'y_range_max= 115'			>> $2

		if ($sample == 300) then
				echo 'y_range_min= 1e-3'			 		 	 >> $2
				echo 'y_range_max= 50'			 			 >> $2
		endif
		echo 'z_print_position_x= 0.75'		>> $2
		echo 'z_print_position_y= 0.80'		>> $2
		echo 'log_scale_y= no'						 >> $2
		echo 'log_scale_y_sub= no'					 >> $2
		echo 'minor_ticks_y_space= 5'			>> $2

	else if ($key == 'plotXY_rbulgevsrdisk') then
		echo 'y_title=  $r_{bulge}/r_{disk}$'  >> $2	
		echo 'y_range_min= 0.5'			 	>> $2
		echo 'y_range_max= 7e3'			>> $2

		echo 'log_scale_y= yes'						 >> $2
		echo 'log_scale_y_sub= yes'					 >> $2


	else if ($key == 'plotXY_vmax') then
		echo 'y_title=  $V_{max}$ [$kms^{-1}$]'  >> $2	
		echo 'y_range_min= 100'			 	>> $2
		echo 'y_range_max= 1e3'			>> $2

		echo 'log_scale_y= yes'						 >> $2
		echo 'log_scale_y_sub= yes'					 >> $2

	else if ($key == 'plotXY_vdisp') then
		echo 'y_title=  $V_{disp}$ [$kms^{-1}$]'  >> $2	
		echo 'y_range_min= 100'			 	>> $2
		echo 'y_range_max= 1e3'			>> $2

		echo 'log_scale_y= yes'						 >> $2
		echo 'log_scale_y_sub= yes'					 >> $2

	else if ($key == 'plotXY_mean_age_stars_disk') then
		echo 'y_title=  $age_{disk}$ [$Gyr$]'  >> $2	
		echo 'y_range_min= 0'			 	>> $2
		echo 'y_range_max= 7'			>> $2

		echo 'log_scale_y= no'						 >> $2
		echo 'log_scale_y_sub= no'					 >> $2

	else if ($key == 'plotXY_mean_age_stars_spheroid') then
		echo 'y_title=  $age_{bulge}$ [$Gyr$]'  >> $2	
		echo 'y_range_min= 0'			 	>> $2
		echo 'y_range_max= 7'			>> $2

		echo 'log_scale_y= no'						 >> $2
		echo 'log_scale_y_sub= no'					 >> $2

	else if ($key == 'con') then
		echo 'y_title= $C_{NFW}$'  >> $2

	else if ($key =~ *violi*) then
		echo 'z_print_position_x= 0.25'		>> $2
		echo 'z_print_position_y= 0.22'		>> $2
		echo 'text_fontsize= 30'				 	 >> $2
		if ($key =~ *mhalo) then 
			echo 'y_title= $\log_{10}$ ($M_{vir}$ [$M_{\odot}$])'				>> $2
			echo 'y_range_min= 9.5'			 		 	 >> $2 
			echo 'y_range_max= 15.5'			 			 >> $2
			echo 'yticks= 10,11,12,13,14,15' 			 >> $2
			echo 'z_print_position_x= 0.82'		>> $2
			echo 'z_print_position_y= 0.9'		>> $2

		else if ($key =~ *mstar) then 
			echo 'y_title= $\log_{10}$ ($M_{*}$ [$M_{\odot}$])'	>> $2
			echo 'y_range_min= 7.5'			 	 >> $2 
			echo 'y_range_max= 12.5'			 >> $2
			echo 'yticks= 8,9,10,11,12' 			 >> $2

		else if ($key =~ *mcold) then 
			echo 'y_title= $\log_{10}$ ($M_{Cold}$ [$M_{\odot}$])'	>> $2
			echo 'y_range_min= 4.5'			 	 >> $2 
			echo 'y_range_max= 11.5'			 >> $2
			echo 'yticks= 5,6,7,8,9,10,11'		 >> $2

		else if ($key =~ *Mzgas) then 
			echo 'y_title= $\log_{10}$ ($M_{Z_{gas}}$ [$M_{\odot}$])'	>> $2
			echo 'y_range_min= 1.5'			 	 >> $2 
			echo 'y_range_max= 10.5'			 >> $2
			echo 'yticks= 2,3,4,5,6,7,8,9,10'		 >> $2

		else if ($key =~ *g-i) then 
			echo 'y_title= $g-i$'	>> $2
			echo 'y_range_min= -0.5'			 	 >> $2 
			echo 'y_range_max= 4.5'			 >> $2
			echo 'yticks= 0,1,2,3,4'		 >> $2
			echo 'z_print_position_x= 0.25'		>> $2
			echo 'z_print_position_y= 0.9'		>> $2

		else if ($key =~ *zcold) then 
			echo 'y_title= $Z_{Cold}$'	>> $2
			echo 'y_range_min= 7.5'			 	 >> $2 
			echo 'y_range_max= 11.2'			 >> $2
			echo 'yticks= 8,9,10,11' 			 >> $2

		else if ($key =~ *_sfr) then 
			echo 'y_title= $\log_{10}$ (SFR [$M_{\odot}yr^{-1}$])'	>> $2
			echo 'y_range_min= -4.5'			 	 >> $2 
			echo 'y_range_max= 2.8'			 >> $2
			echo 'yticks= -4,-3,-2,-1,0,1,2' 			 >> $2

		else if ($key =~ *sfr) then 
			echo 'y_title= $\log_{10}$ (sSFR [$yr^{-1}$])'	>> $2
			echo 'y_range_min= -14.5'			 	 >> $2 
			echo 'y_range_max= -6.5'			 >> $2
			echo 'yticks= -14,-13,-12,-11,-10,-9,-8,-7' 			 >> $2

		else if ($key =~ *SHMR) then 
			echo 'y_title= $\log_{10}$ ($M_{vir}$/$M_*$)'	>> $2
			echo 'y_range_min= -3.5'			 	 >> $2 
			echo 'y_range_max= -0.5'			 >> $2
			echo 'yticks= -3,-2,-1' 			 >> $2

		endif

		echo 'adjust_left= 0.12'				 >> $2
		echo 'adjust_right= 0.99'					>> $2
		echo 'adjust_bottom= 0.18'					>> $2
		echo 'adjust_top= 0.99'					>> $2

		echo 'float_format_y= 0f'         				>> $2

		set use_cb_colors = 5

	else if ($key =~ plotXY_g* || $key =~ plotXY_sta*) then


		if ($key == plotXY_gr) then
			echo 'hline_ypos= 0.5'					 >> $2
			echo 'y_range_min= 0.0'			 		 	 >> $2 
			echo 'y_range_max= 1.08'			 			 >> $2

			echo 'y_title= X$_{z}$/X$_{z_0}$'				>> $2
			echo 'no_last_yticks= False'		 			 >> $2

		else if ($key == plotXY_gr_res) then
			echo 'hline_ypos= -0.5,0,0.5'					 >> $2
			echo 'y_range_min= -1.05'			 		 	 >> $2 
			echo 'y_range_max= 1.05'			 			 >> $2
			echo 'yticks= -1.0,-0.5,0.0,0.5,1.0' 			 >> $2

			echo 'y_title= gr$_{M1}$/gr$_{M2}-1$'				>> $2
			echo 'no_last_yticks= True'		 			 >> $2

		else if ($key == plotXY_stats_N) then
			echo 'y_range_min= 0.0'			 		 	 >> $2 
			echo 'y_range_max= 1.0'			 			 >> $2
			set plot = 'stats_N'
			echo 'y_title= '				>> $2
			echo 'cut_part2= '$plot				 >> $2
			echo 'yticks= 0.25,0.5,0.75,1' 			 >> $2
			echo 'hline_ypos= 0.25,0.5,0.75'			 >> $2

			echo 'print_redshift= False'				 	 >> $2

		else if ($key == plotXY_stats_xbar_M1found_M1) then
			set plot = 'stats_xbar_M1found_M1'
			echo 'y_title= '				>> $2
			echo 'cut_part2= '$plot					 >> $2
			echo 'no_first_yticks= True'		 			 >> $2

		else if ($key == plotXY_stats_xbar_M1found_M2) then
			set plot = 'stats_xbar_M1found_M2'
			echo 'y_title= '				>> $2
			echo 'cut_part2= '$plot					 >> $2
			echo 'no_first_yticks= True'		 			 >> $2

		else if ($key == plotXY_stats_xbar_M1_M2) then
			set plot = 'stats_xbar_M1_M2'
			echo 'y_title= '				>> $2
			echo 'cut_part2= '$plot					 >> $2

		else if ($key == plotXY_stats_xbar_M1-M2_M2) then
			set plot = 'stats_xbar_M1-M2_M2'
			echo 'y_title= '				>> $2
			echo 'cut_part2= '$plot					 >> $2
			echo 'no_first_yticks= True'		 			 >> $2
		else if ($key == plotXY_stats_xbar_MAD) then
			set plot = 'stats_xbar_MAD'
			echo 'y_title= '				>> $2
			echo 'cut_part2= '$plot					 >> $2
			echo 'no_first_yticks= False'		 			 >> $2
			echo 'y_range_min= -1.0'			 		 	 >> $2 
			echo 'y_range_max= 1.0'			 			 >> $2

			echo 'yticks= -1,-0.5,0,0.5,1' 			 >> $2
			echo 'hline_ypos= -0.5,0.0,0.5'			 >> $2
		else
			echo 'y_range_min= 0.0'			 		 	 >> $2 
			echo 'y_range_max= 1.15'			 			 >> $2
			echo 'y_title= mass growth of modelled LRGs'				>> $2
			echo 'y_title= growth'				>> $2

		endif

		echo 'ylabel_pad= 0.03'	         			>> $2
		#echo 'float_format_y= 3f'         				>> $2 

		echo 'error_bars_y= no'					 >> $2
		echo 'error_bars_y_sub= no'					 >> $2

		echo 'log_scale_y= no'						 >> $2
		echo 'log_scale_y_sub= no'					 >> $2
		echo 'minor_ticks_y_space= 0.05'			>> $2
		if ($small_plot == 'True') then
			if ($key =~ plotXY_stats_xbar*) then
				echo 'print_redshift= True'				 >> $2
				echo 'z_print_position_x= 0.2'		>> $2
				echo 'z_print_position_y= 0.83'		>> $2
			else
				echo 'print_redshift= False'				 >> $2
			endif

			echo 'plot_legend= no'				 	 >> $2
			echo 'size_x= 13'				 		 >> $2
			echo 'size_y= 5'				 		 >> $2
			echo 'adjust_left= 0.17'				 >> $2
			echo 'adjust_right= 0.99'					>> $2
			echo 'adjust_bottom= 0.20'					>> $2
			echo 'adjust_top= 0.95'					>> $2
		endif

	else if ($key =~ plotXY_w*) then

		set key2 = 'noyd'
		set key1 = 'noxd'

		echo 'z_print_position_x= 0.22'		>> $2
		echo 'z_print_position_y= 0.90'		>> $2

		echo 'use_loc_co_x= 0.17'		>> $2
		echo 'use_loc_co_y= 0.60'		>> $2

		echo 'log_scale_x= yes'						 >> $2
		echo 'log_scale_x_sub= yes'					 >> $2

		echo 'adjust_bottom= 0.13'					>> $2
		echo 'adjust_left= 0.14'					>> $2
		#echo 'ylabel_pos= 0.45'	         			>> $2

		echo 'size_x= 12'				 		 >> $2
		echo 'size_y= 8.5'				 		 >> $2

		if ($key1 == 'nox') then
			echo 'x_title= ' >> $2
			echo 'no_xticks= True'        					>> $2
		else
			echo 'x_title= $r_p$ $[Mpc]$'	 	>> $2
		endif
		if ($key2 == 'noy') then
			echo 'y_title= ' 				 >> $2
			echo 'no_yticks= True'        					>> $2
			echo 'no_last_yticks= False'        					>> $2
		else
			echo 'y_title= $\log_{10}$ ($r_p w_p$ $[Mpc^2]$)'	>> $2
		endif

		#echo 'minor_ticks_x_space= no'			 	>> $2
		echo 'xticks= 1,2,5,10,25,50,75'	>> $2
		echo 'xticks_minor= 0'					>> $2
		echo 'xticks_minor_offset= 0'				>> $2
		echo 'minor_ticks_y_space= 0.05'			>> $2
		echo 'x_range_min= 0.9'			 		 >> $2
		echo 'x_range_max= 90'			 		 >> $2
		echo 'y_range_min= 1.35'		 			 >> $2
		echo 'y_range_max= 2.95'			 		 >> $2

		echo 'y_range_min= 1.65'		 			 >> $2
		echo 'y_range_max= 3.35'			 		 >> $2

		echo 'size_x= 10'			 		 >> $2
		echo 'adjust_left= 0.18'					>> $2
	else if ($key =~ plotXY_envr*) then

		echo 'y_range_min= 0'			 			 >> $2
		echo 'y_range_max= 1'			 			 >> $2
		echo 'yticks= 0,1'							 >> $2
		echo 'log_scale_y= no'						 >> $2
		echo 'log_scale_y_sub= no'					 >> $2
		echo 'minor_ticks_y_space= 1'			>> $2
		echo 'no_first_last_yticks= True'			>> $2

	else
		echo 'y_title= cSFR density [$M_{\odot}$ $yr^{-1}$ $Mpc^{-3}$]'  >> $2
		#echo 'x_range_min= 0'			 			 >> $2
		#echo 'x_range_max= 7'			 			 >> $2
		#echo 'xticks= 0,1,2,3,4,5,6,7'		>> $2
		#echo 'xticks_minor_offset= 1'				>> $2
		echo 'y_range_min= 3e-7'			 		 	 >> $2 
		echo 'y_range_max=  0.003'			 			 >> $2


		echo 'z_print_position_x= 0.28'		>> $2
		echo 'z_print_position_y= 0.75'		>> $2
	endif

	if ($add_x == 'yes') then
		#add axis x
		echo '0= yes'		>> $22
		echo 'add_axis= yes'				>> $2
		echo 'add_which_col_x= 24'				>> $2
		echo 'log_scale_add_x= no'			>> $2
		echo 'add_axis_plot_exact_ticks= False'				>> $2
		if ($no_top == True) then
			echo 'title_add_x0= '		 	 >> $2
		else
			echo 'title_add_x0= lookback time [$Gyr$]'		 	 >> $2
		endif
		echo 'minor_ticks_x_set_top= off'				>> $2
	endif

	if ($share_y == 'yes') then

		echo 'size_x= 14'				 		 >> $2
		echo 'size_y= 9'				 		 >> $2
		#share axis y
		echo 'share_axis_y= yes'				>> $2
		echo 'share_which_col_y= 6'				>> $2
		echo 'log_scale_share_y= no'			>> $2
		if ($no_right == True) then
			echo 'title_share_y= '		 	 >> $2
		else
			echo 'title_share_y= $M_{cold}/M_*$'		 	 >> $2
		endif
		#echo 'minor_ticks_y_set_right= off'				>> $2
		echo 'minor_ticks_y_share_space= 2'			 	>> $2
		echo 'adjust_right= 0.85'					>> $2

		echo 'share_y_range_min= 0.1'					>> $2
		echo 'share_y_range_max= 52'					>> $2
	endif

	if ($no_y == 'noy') then
		echo 'y_title= ' >> $2
		echo 'no_yticks= True'        					>> $2
		echo 'no_last_yticks= False'        					>> $2
		echo 'plot_legend= no'			>> $2
			
	endif

	echo 'float_format_y= 3f'         >> $2

	if ($key =~ *gr_one) then
		set use_cb_colors = 7
		echo 'legend_ncols= 7'		 			 >> $2
		echo 'use_loc_co_y= 0.83'		>> $2
		echo 'use_loc_co_x= 0.01'		>> $2

		echo 'z_print_position_x= 0.15'		>> $2
		echo 'z_print_position_y= 0.86'		>> $2

		echo 'z_print_position_x= 0.3'		>> $2
		echo 'z_print_position_y= 0.25'		>> $2
		echo 'z_print_position_x= 0.75'		>> $2
		echo 'z_print_position_y= 0.57'		>> $2
		echo 'adjust_top= 0.99'					>> $2
		echo 'adjust_bottom= 0.17'					>> $2

		echo 'hline_ypos= 0.5,1.0'					 >> $2

		echo 'ylabel_pos= 0.5'	         			>> $2    
		echo 'size_y= 6'	         			>> $2
		echo 'size_x= 12'	         			>> $2
		echo 'keyword_change1= plot_legend0,$M_{vir}$' 		>> $2
		echo 'keyword_change2= plot_legend1,$M_{*}$' 		>> $2
		echo 'keyword_change3= plot_legend2,$M_{cold}$' 	>> $2
		echo 'keyword_change4= plot_legend3,$M_{z_{cold}}$' 	>> $2
		echo 'keyword_change5= plot_legend4,$M_{BH}$' 		>> $2
		echo 'keyword_change6= plot_legend5,$r_{1/2}$' 		>> $2
		echo 'keyword_change7= plot_legend6,$Z_{cold}$' 		>> $2

	else if ($key =~ *wp_one) then
		#['0.56','0.59', '0.63', '0.66', '0.7','0.74','0.78','0.82','0.86','0.9','0.94','0.99',\
		set use_cb_colors = 6
		echo 'legend_ncols= 2'			>> $2
		echo 'use_loc_co_x= 0.48'		>> $2
		echo 'use_loc_co_y= 0.8'		>> $2

		echo 'z_print_position_x= 0.22'		>> $2
		echo 'z_print_position_y= 0.90'		>> $2

		if ($3 =~ *400*) then
			if ($sample == 1) then
				echo 'keyword_change1= plot_legend0,z=0.55' 	>> $2
				echo 'keyword_change2= plot_legend1,z=0.59' 	>> $2
				echo 'keyword_change3= plot_legend2,z=0.64' 	>> $2
				echo 'keyword_change4= plot_legend3,z=0.67' 	>> $2
				echo 'keyword_change5= plot_legend4,z=0.71' 		>> $2
				echo 'keyword_change6= plot_legend5,z=0.74' 	>> $2
				echo 'keyword_change7= plot_legend6,z=0.78' 	>> $2
				echo 'keyword_change8= plot_legend7,z=0.8' 	>> $2
				echo 'keyword_change9= plot_legend8,z=0.82' 	>> $2
				echo 'keyword_change10= plot_legend9,z=0.84' 	>> $2
				echo 'keyword_change11= plot_legend10,z=0.89' 	>> $2
				echo 'keyword_change12= plot_legend11,z=1.0' 	>> $2
				echo 'keyword_change13= plot_legend12,z=1.0' 	>> $2
				echo 'keyword_change14= plot_legend13,z=1.0' 	>> $2
			else
				echo 'keyword_change1= plot_legend0,z=0.55' 	>> $2
				echo 'keyword_change2= plot_legend1,z=0.59' 	>> $2
				echo 'keyword_change3= plot_legend2,z=0.64' 	>> $2
				echo 'keyword_change4= plot_legend3,z=0.71' 	>> $2
				echo 'keyword_change5= plot_legend4,z=0.78' 		>> $2
				echo 'keyword_change6= plot_legend5,z=0.8' 	>> $2
				echo 'keyword_change7= plot_legend6,z=0.84' 	>> $2
				echo 'keyword_change8= plot_legend7,z=0.89' 	>> $2
				echo 'keyword_change9= plot_legend8,z=1.0' 	>> $2
				echo 'keyword_change10= plot_legend9,z=1.26' 	>> $2
				echo 'keyword_change11= plot_legend10,z=1.56' 	>> $2
				echo 'keyword_change12= plot_legend11,z=2.02' 	>> $2
				echo 'keyword_change13= plot_legend12,z=1.0' 	>> $2
				echo 'keyword_change14= plot_legend13,z=1.0' 	>> $2
			endif

		else
			if ($sample == 1) then
				echo 'keyword_change1= plot_legend0,z=0.56' 		>> $2
				echo 'keyword_change2= plot_legend1,z=0.59' 		>> $2
				echo 'keyword_change3= plot_legend2,z=0.63' 	>> $2
				echo 'keyword_change4= plot_legend3,z=0.66' 	>> $2
				echo 'keyword_change5= plot_legend4,z=0.7' 		>> $2
				echo 'keyword_change6= plot_legend5,z=0.74' 		>> $2
				echo 'keyword_change7= plot_legend6,z=0.78' 	>> $2
				echo 'keyword_change8= plot_legend7,z=0.82' 	>> $2
				echo 'keyword_change9= plot_legend8,z=0.86' 	>> $2
				echo 'keyword_change10= plot_legend9,z=0.9' 	>> $2
				echo 'keyword_change11= plot_legend10,z=0.94' 	>> $2
				echo 'keyword_change12= plot_legend11,z=0.99' 	>> $2
				echo 'keyword_change13= plot_legend12,z=1.0' 	>> $2
				echo 'keyword_change14= plot_legend13,z=1.0' 	>> $2

			else if ($sample == 'more') then
				set use_cb_colors = 2
				echo 'legend_ncols= 2'			>> $2
				echo 'use_loc_co_x= 0.25'		>> $2
				echo 'use_loc_co_y= 0.15'		>> $2
				echo 'keyword_change1= plot_legend0,red' 		>> $2
				echo 'keyword_change1= plot_legend0,low-$Z_{cold}$' 		>> $2
				echo 'keyword_change2= plot_legend1,high-$Z_{cold}$' 	>> $2
			else if ($sample == 'pop') then
				echo 'print_redshift2=  low-$Z_{cold}$:         high-$Z_{cold}$:'				>> $2
				echo 'z_print_position_x2= 0.25'		>> $2
				echo 'z_print_position_y2= 0.28'		>> $2

				echo 'use_loc_co_x= 0.25'		>> $2
				echo 'use_loc_co_y= 0.15'		>> $2
				echo 'keyword_change1= plot_legend0,PopB+k' 	>> $2
				echo 'keyword_change2= plot_legend1,PopB+f' 	>> $2
				echo 'keyword_change3= plot_legend2,PopA+k' 		>> $2
				echo 'keyword_change4= plot_legend3,PopA+f' 		>> $2
			else if ($sample == 'envr') then
				set use_cb_colors = 6
				echo 'print_redshift2=  low-$Z_{cold}$:           high-$Z_{cold}$:'				>> $2
				echo 'z_print_position_x2= 0.25'		>> $2
				echo 'z_print_position_y2= 0.28'		>> $2

				echo 'use_loc_co_x= 0.25'		>> $2
				echo 'use_loc_co_y= 0.15'		>> $2
				echo 'keyword_change1= plot_legend0,knots' 		>> $2
				echo 'keyword_change2= plot_legend1,filaments' 		>> $2
				echo 'keyword_change3= plot_legend2,knots' 	>> $2
				echo 'keyword_change4= plot_legend3,filaments' 	>> $2
				echo 'keyword_change5= plot_legend4,knots' 		>> $2
				echo 'keyword_change6= plot_legend5,filaments' 		>> $2
			else
				echo 'keyword_change1= plot_legend0,z=0.56' 		>> $2
				echo 'keyword_change2= plot_legend1,z=0.59' 		>> $2
				echo 'keyword_change3= plot_legend2,z=0.63' 	>> $2
				echo 'keyword_change4= plot_legend3,z=0.7' 	>> $2
				echo 'keyword_change5= plot_legend4,z=0.78' 		>> $2
				echo 'keyword_change6= plot_legend5,z=0.86' 		>> $2
				echo 'keyword_change7= plot_legend6,z=0.9' 	>> $2
				echo 'keyword_change8= plot_legend7,z=0.94' 	>> $2
				echo 'keyword_change9= plot_legend8,z=1.03' 	>> $2
				echo 'keyword_change10= plot_legend9,z=1.22' 	>> $2
				echo 'keyword_change11= plot_legend10,z=1.54' 	>> $2
				echo 'keyword_change12= plot_legend11,z=2.03' 	>> $2
				echo 'keyword_change13= plot_legend12,z=1.0' 	>> $2
				echo 'keyword_change14= plot_legend13,z=1.0' 	>> $2
			endif
		endif
	else

		#echo 'use_loc_co_y= 0.82'		>> $2
		#echo 'use_loc_co_x= 0.01'		>> $2

		#echo 'z_print_position_x= 0.73'		>> $2
		#echo 'z_print_position_y= 0.77'		>> $2

		echo 'use_loc_co_x= 0.02'		>> $2
		echo 'use_loc_co_y= 0.58'		>> $2

		echo 'legend_ncols= 1'			>> $2

		if ($sample == 1b) then
			set use_cb_colors = 14
		else
			set use_cb_colors = 9
		endif

		if ($sample == 1 || $sample == 1a || $sample == 1b) then
			echo 'keyword_change1= plot_legend0,low-SFR' 	>> $2
			echo 'keyword_change2= plot_legend1,low-SFR' 	>> $2
			echo 'keyword_change3= plot_legend2,high-SFR' 	>> $2
			echo 'keyword_change4= plot_legend3,passive' 	>> $2
			echo 'keyword_change5= plot_legend4,active' 	>> $2
			echo 'keyword_change6= plot_legend5,red' 	>> $2
			echo 'keyword_change7= plot_legend6,blue' 	>> $2
			echo 'keyword_change8= plot_legend7,low-$Z_{Cold}$' 	>> $2
			echo 'keyword_change9= plot_legend8,high-$Z_{Cold}$' 	>> $2
			echo 'keyword_change10= plot_legend9,$M_*>11$' 	>> $2
			echo 'keyword_change11= plot_legend10,$Z_{Cold}<9$' 	>> $2
			echo 'keyword_change12= plot_legend11,$Z_{Cold}>9$' 	>> $2
			echo 'keyword_change13= plot_legend12,red+$M_{vir}\sim13$' 	>> $2
			echo 'keyword_change14= plot_legend13,low-SFR+$M_{vir}\sim13$' 	>> $2

			if ($reduced == True) then
				echo '0= yes'		>> $21
				echo '1= no'		>> $21
				echo '2= no'		>> $21
				echo '3= no'		>> $21
				echo '4= no'		>> $21
				echo '5= yes'		>> $21
				echo '6= yes'		>> $21
				echo '7= yes'		>> $21
				echo '8= yes'		>> $21
				echo '9= no'		>> $21
				echo '10= no'		>> $21
				echo '11= no'		>> $21
				echo '12= yes'		>> $21
				echo '13= no'		>> $21
			endif

			if ($reduced == True) then
				echo '0= yes'		>> $21
				echo '1= no'		>> $21
				echo '2= no'		>> $21
				echo '3= yes'		>> $21
				echo '4= no'		>> $21
				echo '5= yes'		>> $21
				echo '6= no'		>> $21
				echo '7= yes'		>> $21
				echo '8= yes'		>> $21
				echo '9= no'		>> $21
				echo '10= no'		>> $21
				echo '11= no'		>> $21
				echo '12= no'		>> $21
				echo '13= no'		>> $21
			endif


		else if ($sample == final) then
			set use_cb_colors = 5
			echo 'keyword_change1= plot_legend0,low-SFR' 	>> $2
			echo 'keyword_change2= plot_legend1,passive' 	>> $2
			echo 'keyword_change3= plot_legend2,red' 	>> $2
			echo 'keyword_change4= plot_legend3,low-$Z_{Cold}$' 	>> $2
			echo 'keyword_change5= plot_legend4,high-$Z_{Cold}$' 	>> $2

		else if ($sample == 2a) then
			set use_cb_colors = 11
			echo 'use_loc_co_y= 0.82'		>> $2
			echo 'keyword_change1= plot_legend0,all' 	>> $2
			echo 'keyword_change2= plot_legend1,$M_*>11$' 	>> $2
			echo 'keyword_change3= plot_legend2,$M_{vir}>12$' 	>> $2
			echo 'keyword_change4= plot_legend3,$M_*<11$' 	>> $2
			echo 'keyword_change5= plot_legend4,$M_{vir}<12$' 	>> $2
			echo 'keyword_change6= plot_legend5,$Z_{Cold}<9$' 	>> $2
			echo 'keyword_change7= plot_legend6,$Z_{Cold}>9$' 	>> $2
			echo 'keyword_change8= plot_legend7,red+$M_*\sim11$' 	>> $2
			echo 'keyword_change9= plot_legend8,red+$M_{vir}\sim13$' 	>> $2
			echo 'keyword_change10= plot_legend9,low-SFR+$M_*\sim11$' 	>> $2
			echo 'keyword_change11= plot_legend10,low-SFR+$M_{vir}\sim13$' 	>> $2

		else if ($sample == 2b) then
			#reduced sample for paper
			set use_cb_colors = 7
			echo 'use_loc_co_y= 0.82'		>> $2
			echo 'keyword_change1= plot_legend0,all' 	>> $2
			echo 'keyword_change2= plot_legend1,$M_*>11$' 	>> $2
			echo 'keyword_change3= plot_legend2,$M_{vir}>12$' 	>> $2
			echo 'keyword_change4= plot_legend3,$Z_{Cold}<9$' 	>> $2
			echo 'keyword_change5= plot_legend4,$Z_{Cold}>9$' 	>> $2
			echo 'keyword_change6= plot_legend5,red+$M_{vir}\sim13$' 	>> $2
			echo 'keyword_change7= plot_legend6,low-SFR+$M_{vir}\sim13$' 	>> $2
		else if ($sample == 3) then
			set use_cb_colors = 6
			echo 'keyword_change1= plot_legend0,high $M_{vir}$' 	>> $2
			echo 'keyword_change2= plot_legend1,high $M_*$' 	>> $2
			echo 'keyword_change3= plot_legend2,high $M_{Z_{gas}}$' 	>> $2
			echo 'keyword_change4= plot_legend3,low $M_{vir}$' 	>> $2
			echo 'keyword_change5= plot_legend4,low $M_*$' 	>> $2
			echo 'keyword_change6= plot_legend5,low $M_{Z_{gas}}$' 	>> $2
		else if ($sample == 300) then
			set use_cb_colors = 6
			echo 'keyword_change1= plot_legend0,Gal-r0001' 	>> $2
			echo 'keyword_change2= plot_legend1,high $M_*$' 	>> $2
			echo 'keyword_change3= plot_legend2,high $M_{Z_{gas}}$' 	>> $2
			echo 'keyword_change4= plot_legend3,low $M_{vir}$' 	>> $2
			echo 'keyword_change5= plot_legend4,low $M_*$' 	>> $2
			echo 'keyword_change6= plot_legend5,low $M_{Z_{gas}}$' 	>> $2
		else
			echo 'keyword_change1= plot_legend0,red' 	>> $2
			echo 'keyword_change2= plot_legend1,red' 	>> $2
			echo 'keyword_change3= plot_legend2,red' 	>> $2
			echo 'keyword_change4= plot_legend3,red' 	>> $2
			echo 'keyword_change5= plot_legend4,blue' 	>> $2
			echo 'keyword_change6= plot_legend5,blue' 	>> $2
			echo 'keyword_change7= plot_legend6,blue' 	>> $2
			echo 'keyword_change8= plot_legend7,blue' 	>> $2
		endif
	endif

else if ($1 == plotXY && $plotXY_key == 'zevol') then
	set no_x = nox
	set no_y = noyx
	set no_top = False
	set xIsz = 'True'

	set sample = 2
	set key = $24

	set yerr = 3

	echo '0= yes'		>> $21
	echo '1= no'		>> $21
	echo '2= no'		>> $21
	echo '3= yes'		>> $21
	echo '4= no'		>> $21
	echo '5= no'		>> $21
	echo '6= yes'		>> $21
	echo '7= no'		>> $21
	echo '8= no'		>> $21
	echo '9= sep'		>> $21
	echo '10= sep'		>> $21
	echo '11= sep'		>> $21
	echo '12= yes'		>> $21
	echo '13= yes'		>> $21
	echo '14= yes'		>> $21
	echo '15= yes'		>> $21
	echo '16= yes'		>> $21

	echo 'legend_sep_xpos= 0.80'	>> $2
	echo 'legend_sep_ypos= 0.78'	>> $2
	echo 'legend_sep_offset= 0.04'	>> $2
	echo 'plot_legend= yes'			>> $2

	echo 'lw_offset= 0'				 		 >> $2
	echo 'text_fontsize= 26'					>> $2
	echo 'legend_fontsize= 20'				 	 >> $2
	echo 'print_redshift= yes'				 	 >> $2

	set small_plot = 'False'

	if ( $small_plot == 'True') then
		echo 'size_x= 13'				 		 >> $2
		echo 'size_y= 4.5'				 		 >> $2
		echo 'adjust_left= 0.16'					>> $2
		echo 'adjust_right= 0.99'					>> $2
		echo 'adjust_bottom= 0.20'					>> $2
		echo 'adjust_top= 0.80'					>> $2

		echo 'z_print_position_x= 0.005'		>> $2
		echo 'z_print_position_y= 0.72'		>> $2

		echo 'use_loc_co_x= 0.45'		>> $2
		echo 'use_loc_co_y= 0.55'		>> $2

	else
		echo 'size_x= 10'				 		 >> $2
		echo 'size_y= 6'				 		 >> $2
		echo 'adjust_left= 0.18'					>> $2
		echo 'adjust_right= 0.99'					>> $2
		echo 'adjust_bottom= 0.16'					>> $2
		echo 'adjust_top= 0.99'					>> $2

		echo 'z_print_position_x= 0.23'		>> $2
		echo 'z_print_position_y= 0.87'		>> $2

		echo 'use_loc_co_x= 0.22'		>> $2
		echo 'use_loc_co_y= 0.18'		>> $2
	endif

	echo 'ylabel_pos= 0.55'	         			>> $2    

	if ( $xIsz == 'True') then
		echo 'x_title= z'				 			 >> $2 
		echo 'minor_ticks_x_space= no'			 	>> $2
		echo 'minor_ticks_y_space= no'			 	>> $2
		echo 'xticks_minor= 0.1'					>> $2
		echo 'x_range_min= -0.15'			 			 >> $2
		echo 'x_range_max= 3.15'			 			 >> $2
		#echo 'xticks= 0,0.1,0.5,1,2,3'				>> $2
		if ($no_x == 'nox') then
			echo 'x_title= ' >> $2
			echo 'no_xticks= True'        					>> $2	
		endif		

	else if ($xIsz = 'Gyr')
		echo 'x_title= lookback time [$Gyr$]' 			 >> $2 
		echo 'log_scale_x= no'						 >> $2
		echo 'log_scale_x_sub= no'					 >> $2

		echo 'x_range_min= 0'			 			 >> $2
		echo 'x_range_max= 6.74389'			 			 >> $2

		echo 'minor_ticks_x_space= 0.25'			 	>> $2
		echo 'minor_ticks_y_space= no'			 	>> $2
		echo 'xticks_minor= 0'					>> $2
	else
		echo 'log_scale_x= yes'						 >> $2
		echo 'log_scale_x_sub= yes'					 >> $2
		echo 'minor_ticks_y_space= 0.25'			 	>> $2
	endif

	echo 'xticks_minor_offset= 0'				>> $2
	#echo 'float_format_y= 0f'         				>> $2 
	echo 'no_first_last_yticks= True'		 			 >> $2

	echo 'subplot_errorbars= True' 		 			 >> $2
	echo 'filled_between_subplot= True' 		 		 >> $2

	if ( $xIsz == 'True' ) then
		set add_x = 'yess'
	else
		set add_x = 'no'
	endif

	if ($key == 'plotXY_mstar') then
		echo 'y_title=   $\log_{10}$ ($M_{\star}$ [$M_{\odot}$])'  >> $2	
		echo 'y_range_min= 8.7'			 	>> $2  
		echo 'y_range_max= 11.8'			>> $2

	else if ($key == 'plotXY_sfr') then
		echo 'y_title=   $\log_{10}$ (SFR [$M_{\odot}$ $yr^{-1}$])'  >> $2
		echo 'y_range_min= -4'			 		 	 >> $2 
		echo 'y_range_max= 2'			 			 >> $2

	else if ($key == 'plotXY_ssfr') then
		echo 'y_title=   $\log_{10}$ (sSFR [$yr^{-1}$])'  >> $2
		echo 'y_range_min= -13.5'			 		 	 >> $2
		echo 'y_range_max= -7.5'			 			 >> $2

	else if ($key == 'plotXY_mhalo') then
		echo 'y_title=  $\log_{10}$ ($M_{200c}$ [$M_{\odot}$])'  >> $2
		echo 'y_range_min= 10'			 		 	 >> $2
		echo 'y_range_max= 14'			 			 >> $2
	endif

	if ($no_y == 'noy') then
		echo 'y_title= ' >> $2
		echo 'no_yticks= True'        					>> $2
		echo 'no_first_last_yticks= False'		 			 >> $2	
	endif	

	if ($yerr == '1') then
		echo 'error_bars_y= yes'					 >> $2
		echo 'error_bars_y_sub= no'					 >> $2
	endif

	echo 'keyword_change1= plot_legend0,Gal-MD1000' 	>> $2
	echo 'keyword_change2= plot_legend1,Gal-MD1000 CUT2' 	>> $2
	echo 'keyword_change3= plot_legend2,Gal-MD1000 CUT3' 	>> $2
	echo 'keyword_change4= plot_legend3,Gal-SMD400' 	>> $2
	echo 'keyword_change5= plot_legend4,Gal-SMD400 CUT2' 	>> $2
	echo 'keyword_change6= plot_legend5,Gal-SMD400 CUT3' 	>> $2
	echo 'keyword_change7= plot_legend6,Ill-TNG300' 	>> $2
	echo 'keyword_change8= plot_legend7,Ill-TNG300 CUT2' 	>> $2
	echo 'keyword_change9= plot_legend8,Ill-TNG300 CUT3' 	>> $2

else if ($1 == plotXY && $plotXY_key == 'wp_mstar') then
	echo 'mylw= 6'				 		 	 >> $2
	set plot = 'wp'
	set key2 = 'noyb'
	set key1 = 'noxb'

	echo 'log_scale_x= yes'						 >> $2
	echo 'log_scale_x_sub= yes'					 >> $2

	if ($key1 == 'nox') then
		echo 'x_title= ' >> $2
		echo 'no_xticks= True'        					>> $2
		echo 'no_first_yticks= True'        					>> $2
	else
		echo 'x_title= $r_p$ $[Mpc]$'	 	>> $2
	endif

	if ($key2 == 'noy') then
		echo 'y_title= ' >> $2
		echo 'no_yticks= True'        					>> $2
	else
		echo 'y_title= $\log_{10}$ ($r_p w_p$ $[Mpc^2]$)' >> $2
	endif

	#echo 'adjust_left= 0.12'					>> $2
	#echo 'adjust_right= 0.99'					>> $2
	#echo 'adjust_bottom= 0.12'					>> $2
	#echo 'adjust_top= 0.99'					>> $2

	echo 'x_range_min= 0.5'			 		 >> $2
	echo 'x_range_max= 87'			 		 >> $2
	echo 'y_range_min= 2.3'			 		 >> $2
	echo 'y_range_max= 3.55'			 		 >> $2

	echo 'error_bars_y= no'				 >> $2

	echo 'use_loc_co_x= 0.82'			>> $2
	echo 'use_loc_co_y= 0.76'			>> $2

	echo 'z_print_position_x= 0.16'		>> $2
	echo 'z_print_position_y= 0.78'		>> $2

	echo 'minor_ticks_x_space= no'			 	>> $2
	echo 'xticks= 1,10'	>> $2
	echo 'xticks_minor= default'					>> $2



	#echo 'xticks_minor= 0'					>> $2

	echo $plot'x_range_min= 0.5'			 		 >> $2
	echo $plot'x_range_max= 87'			 		 >> $2

	echo 'print_redshift_sep= CMASS DR12:' >> $2
	echo 'legend_sep_xpos= 0.63'	>> $2
	echo 'legend_sep_ypos= 0.16'	>> $2
	echo 'legend_sep_offset= 0.035'	>> $2

	echo 'text_fontsize= 24'					>> $2
	echo 'legend_fontsize= 18'				 	 >> $2
	echo 'print_redshift= yes'				 	 >> $2

	if ($use_cb_colors == 5) then
		echo 'keyword_change1= plot_legend0,cut1' 	>> $2
		echo 'keyword_change2= plot_legend1,cut2' 	>> $2
		echo 'keyword_change3= plot_legend2,cut3' 	>> $2
		echo 'keyword_change4= plot_legend3,cut4' 	>> $2
		echo 'keyword_change5= plot_legend4,cut5' 	>> $2
	else
		echo 'keyword_change1= plot_legend0,bin1' 	>> $2
		echo 'keyword_change2= plot_legend1,bin2' 	>> $2
		echo 'keyword_change3= plot_legend2,bin3' 	>> $2
	endif

	echo '0= yes'		>> $21
	echo '1= yes'		>> $21
	echo '2= yes'		>> $21
	echo '3= yes'		>> $21
	echo '4= yes'		>> $21
	echo '5= sep'		>> $21
	echo '6= sep'		>> $21
	echo '7= sep'		>> $21
	echo '8= sep'		>> $21
	echo '9= sep'		>> $21
	echo '10= sep'		>> $21

	#fixed n_z
	#=======================
	#Chabrier
	#echo 'keyword_change1= plot_legend0,$11.21<\log_{10}(M_*$ $[M_{\odot}])<11.30$' 	>> $2
	#echo 'keyword_change2= plot_legend1,$11.32<\log_{10}(M_*$ $[M_{\odot}])<11.48$' 	>> $2
	#echo 'keyword_change3= plot_legend2,$\log_{10}(M_*$ $[M_{\odot}])>11.39$' 	>> $2

	#Kroupa
	#echo 'keyword_change1= plot_legend0,$\log_{10}(M_*$ $[M_{\odot}])>11.27$' 	>> $2
	#echo 'keyword_change2= plot_legend1,$\log_{10}(M_*$ $[M_{\odot}])>11.34$' 	>> $2
	#echo 'keyword_change3= plot_legend2,$\log_{10}(M_*$ $[M_{\odot}])>11.43$' 	>> $2
	#echo 'keyword_change4= plot_legend3,$\log_{10}(M_*$ $[M_{\odot}])>11.52$' 	>> $2
	#echo 'keyword_change5= plot_legend4,$\log_{10}(M_*$ $[M_{\odot}])>11.63$' 	>> $2

	#echo 'keyword_change1= plot_legend0,$11.28<\log_{10}(M_*$ $[M_{\odot}])<11.36$' 	>> $2
	#echo 'keyword_change2= plot_legend1,$11.36<\log_{10}(M_*$ $[M_{\odot}])<11.48$' 	>> $2
	#echo 'keyword_change3= plot_legend2,$\log_{10}(M_*$ $[M_{\odot}])>11.46$' 	>> $2

#####################################################################################################################################################
else if ($1 == plotXY && $plotXY_key == 'wp') then
	echo 'mylw= 8'				 		 	 >> $2
	set plot = 'wp'
	set key2 = 'noyb'
	set key1 = 'noxb'
	set CUT = 'CUT33'
	set sel = 'CMASS_wp'

	echo 'log_scale_x= no'						 >> $2
	echo 'log_scale_x_sub= no'					 >> $2

	#echo 'adjust_left= 0.14'					>> $2
	#echo 'adjust_right= 0.99'					>> $2
	#echo 'adjust_bottom= 0.12'					>> $2
	#echo 'adjust_top= 0.99'					>> $2

	if ($key1 == 'nox') then
		echo 'x_title= ' >> $2
		echo 'no_xticks= True'        					>> $2
	else
		echo 'x_title= $r$ $[Mpc]$'	 	 >> $2
	endif

	if ($key2 == 'noy') then
		echo 'y_title= ' >> $2
		echo 'no_yticks= True'        					>> $2
	else
		echo 'y_title= $r^{2}$ $\xi(r)$ $[Mpc^2]$' >> $2
	endif

	echo 'x_range_min= 0.5'			 		 >> $2
	echo 'x_range_max= 150'			 		 >> $2
	echo 'y_range_min= 2.05'			 		 >> $2
	echo 'y_range_max= 3.6'			 		 >> $2

	echo 'error_bars_y= yes'				 >> $2
	echo 'error_bars_y_sub= yes'				 >> $2
	echo 'nr_plots_x= 1'			 		 >> $2
	echo 'nr_plots_y= 1'			 		 >> $2

	echo 'hratio1= 3'			 		 >> $2
	echo 'hratio2= 1'			 		 >> $2
	echo 'tight_layout= False'		 		 >> $2


	#middle
	echo 'use_loc_co_x= 0.53'			>> $2
	#three curves
	#echo 'use_loc_co_y= 0.38'			>> $2

	#top right
	#echo 'use_loc_co_x= 0.72'			>> $2
	#two curves
	echo 'use_loc_co_y= 0.81'			>> $2

	echo 'use_loc_co_x= 0.61'			>> $2
	echo 'use_loc_co_y= 0.78'			>> $2

	#xi
	if ($key2 == 'noy') then
		echo $plot'y_title= ' >> $2
		echo 'no_yticks= True'        					>> $2
	else
		echo $plot'y_title= $r^{2}$ $\xi(r)$ $[Mpc^2]$' >> $2
	endif

	echo $plot'cut_part1= '					 >> $2
	echo $plot'cut_part2= '					 >> $2

	echo $plot'x_range_min= 0.5'			 		 >> $2
	echo $plot'x_range_max= 150'			 		 >> $2
	echo $plot'y_range_min= 0.4'		 			 >> $2
	echo $plot'y_range_max= 2.6'			 		 >> $2

	echo $plot'vline_xpos= 2.0'				 	 >> $2
	echo $plot'hline_ypos= False'					 >> $2

	echo $plot'print_redshift= yes'				 	 >> $2
	echo $plot'plot_legend= yes'	 				 >> $2

	if ($sel == 'CMASS_xi') then

		echo 'use_loc_co_x= 0.50'			>> $2
		echo 'use_loc_co_y= 0.62'			>> $2

		echo 'minor_ticks_y_space= 5'			 	>> $2
		echo 'z_print_position_x= 5'		>> $2
		echo 'z_print_position_y= -30'		>> $2

		echo 'xticks= 5,25,50,75,100,125,150,175,200'		>> $2
		echo 'xticks_minor= 5'					>> $2
		echo 'xticks_minor_offset= 0'				>> $2

		echo $plot'x_range_min= 1'			 		 >> $2
		echo $plot'x_range_max= 205'			 		 >> $2
		echo $plot'y_range_min= -40'		 			 >> $2
		echo $plot'y_range_max= 150'			 		 >> $2
		echo $plot'vline_xpos= 100.0'				 	 >> $2
	endif


	echo $plot'size_x= 12'				 		 >> $2
	echo $plot'size_y= 10'				 		 >> $2


	set plot = 'xi_ref'
	#xi_ref
	echo $plot'x_title= $r$ $[Mpc]$'	 	 >> $2
	echo $plot'cut_part1= '					 >> $2

	if ($key1 == 'nox') then
		echo $plot'x_title= ' >> $2
		echo 'no_xticks= True'        					>> $2
	else
		echo $plot'x_title= $r$ $[Mpc]$'	 	 >> $2
	endif
	if ($key2 == 'noy') then
		echo $plot'y_title= ' 				 >> $2
		echo $plot'cut_part2= '				 >> $2
		echo 'no_yticks= True'   			 >> $2
	else
		echo $plot'y_title= $\xi(r)$/' 				 >> $2
		echo $plot'cut_part2= bar'				 >> $2
	endif

	echo $plot'y_range_min= -0.5'		 			 >> $2
	echo $plot'y_range_max= 0.5'			 		 >> $2
	echo $plot'x_range_min= -0.5'			 		 >> $2
	echo $plot'x_range_max= 2.3'			 		 >> $2

	echo $plot'hline_ypos= 0.0'					 >> $2
	echo $plot'vline_xpos= 2.0'				 	 >> $2
	echo $plot'print_redshift= '				 	 >> $2
	echo $plot'plot_legend= no'	 				 >> $2

	if ($sel == 'CMASS_xi') then
		echo $plot'x_range_min= 1'			 		 >> $2
		echo $plot'x_range_max= 205'			 		 >> $2
		echo $plot'y_range_min= -1'		 			 >> $2
		echo $plot'y_range_max= 1'			 		 >> $2
		echo $plot'vline_xpos= 100.0'				 	 >> $2
	endif

	echo $plot'size_x= 12'				 		 >> $2
	echo $plot'size_y= 10'				 		 >> $2

	set plot = 'wp'
#wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp
	echo $plot'cut_part1= '					 >> $2
	echo $plot'cut_part2= '					 >> $2
	if ($key1 == 'nox') then
		echo $plot'x_title= ' >> $2
		echo 'no_xticks= True'        					>> $2
	else
		echo $plot'x_title= $\log_{10}$ ($r_p$ $[Mpc]$)'	 	>> $2
	endif
	if ($key2 == 'noy') then
		echo $plot'y_title= ' 				 >> $2
	else
		echo $plot'y_title= $\log_{10}$ ($r_p w_p$ $[Mpc^2]$)'	>> $2
	endif

	echo $plot'x_range_min= -0.8'			 		 >> $2
	echo $plot'x_range_max= 1.7'			 		 >> $2
	echo $plot'y_range_min= 0'		 			 >> $2
	echo $plot'y_range_max= 3.5'			 		 >> $2

	if ($sel == 'CMASS_wp') then
		echo 'z_print_position_x= 0.65'		>> $2
		echo 'z_print_position_y= 2.74'		>> $2
		echo 'z_print_position_y= 2.6'		>> $2
		#echo 'xticks= -0.25,0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0'	>> $2
		echo 'minor_ticks_x_space= no'			 	>> $2
		echo 'xticks= 1,10'	>> $2
		echo 'xticks_minor= 0'					>> $2
		echo 'xticks_minor_offset= 0'				>> $2

		echo $plot'x_range_min= 0.5'			 		 >> $2
		echo $plot'x_range_max= 87'			 		 >> $2
		echo $plot'y_range_min= 2.15'		 			 >> $2
		echo $plot'y_range_max= 3.05'			 		 >> $2
	else
		echo 'z_print_position_x= 0.7'		>> $2
		echo 'z_print_position_y= 0.8'		>> $2
		echo 'xticks= 0.5,1,2,5,10,50'	>> $2
		#echo 'xticks= -0.5,0.0,0.5,1.0,1.5'	>> $2
		echo 'minor_ticks_x_space= no'			 	>> $2
		echo 'xticks_minor= 0'					>> $2
		echo 'xticks_minor_offset= 0'				>> $2

		echo 'x_range_min= 0.5'			 		 >> $2
		echo 'x_range_max= 60'			 		 >> $2
		echo 'y_range_min= 1.7'		 			 >> $2
		echo 'y_range_max= 3.1'			 		 >> $2

		echo 'use_loc_co_x= 0.16'			>> $2
		echo 'use_loc_co_y= 0.70'			>> $2
		echo 'no_last_yticks= True'        					>> $2

	endif

	echo $plot'hline_ypos= False'					 >> $2
	echo $plot'vline_xpos= False'					 >> $2
	echo $plot'print_redshift= '				 	 >> $2
	echo $plot'plot_legend= yes'	 				 >> $2

	set plot = 'wp_ref'	
	#wp_ref
	echo $plot'cut_part1= '					 >> $2
	if ($key1 == 'nox') then
		echo $plot'x_title= ' >> $2
		echo 'no_xticks= True'        					>> $2
	else
		echo $plot'x_title= $r_p$ $[Mpc]$'	 	>> $2
	endif
	if ($key2 == 'noy') then
		echo $plot'y_title= ' 				 >> $2
		echo $plot'cut_part2= '					 >>	$2
	else
		echo $plot'y_title= $w_p$/' 				 >> $2
		echo $plot'cut_part2= data_res'					 >> $2
	endif

	echo $plot'y_range_min= -0.58'		 			 >> $2
	echo $plot'y_range_max= 0.58'			 		 >> $2
	echo $plot'x_range_min= -0.8'			 		 >> $2
	echo $plot'x_range_max= 1.7'			 		 >> $2

	echo 'log_scale_x= yes'						 >> $2
	echo 'log_scale_x_sub= yes'					 >> $2

	if ($sel == 'CMASS_wp') then
		echo $plot'x_range_min= 0.5'			 		 >> $2
		echo $plot'x_range_max= 87'			 		 >> $2
	endif

	echo $plot'hline_ypos= 0.0'					 >> $2
	echo $plot'vline_xpos= False'					 >> $2
	echo $plot'print_redshift= '				 	 >> $2
	echo $plot'plot_legend= no'	 				 >> $2

	#DEFAULT
	echo 'x_title= $r_p$ $[Mpc]$' >> $2
	#echo 'y_title= $\log_{10}$ $w(r_p)$' >> $2

#wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp wp


	#echo 'text_fontsize= 28'					>> $2
	echo 'legend_fontsize= 30'				 	 >> $2

	#echo 'keyword_change1= plot_legend0,'$prefix'-cols' 	>> $2
	#echo 'keyword_change2= plot_legend1,'$prefix'-dens'  	>> $2
	#echo 'keyword_change3= plot_legend2,'$prefix'-mass' 	>> $2
	#echo 'keyword_change4= plot_legend3,Gal-all'  	>> $2
	#echo 'keyword_change4= plot_legend3,Rodriguez+16 BigMD HAM' 	>> $2
	#echo 'keyword_change4= plot_legend3,CMASS DR12' 	>> $2

	#echo 'keyword_change1= plot_legend0,Gal-cols $M_r>-21.5$' 	>> $2
	#echo 'keyword_change1= plot_legend0,$n_{CMASS}$' 	>> $2
	#echo 'keyword_change2= plot_legend1,Guo+13'  	>> $2
	#echo 'keyword_change3= plot_legend2,$g-i>2.35$'  	>> $2
	
	#echo 'keyword_change1= plot_legend0,Gal-dens' 	>> $2
	#echo 'keyword_change2= plot_legend1,Gal-cross'  	>> $2
	#echo 'keyword_change3= plot_legend2,Gal-cross2'  	>> $2

	echo 'keyword_change1= plot_legend0,Gal-dens' 	>> $2
	echo 'keyword_change2= plot_legend1,Gal400-dens'  	>> $2
	echo 'keyword_change3= plot_legend2,Gal2-dens mstar'  	>> $2

	set sign = '<'
	set scale = 'd'
	if ($scale == 'l') then
		echo 'keyword_change1= plot_legend0,$all$' 	>> $2
		echo 'keyword_change2= plot_legend1,$>11.25$' 	>> $2
		echo 'keyword_change3= plot_legend2,$>11.35$' 	>> $2
		echo 'keyword_change4= plot_legend3,$>11.45$' 	>> $2
		echo 'keyword_change5= plot_legend4,$>11.55$' 	>> $2
		echo 'keyword_change6= plot_legend5,$>11.65$' 	>> $2
		echo 'keyword_change7= plot_legend6,CMASS DR12' 	>> $2

	else if ($scale == 's') then
		echo 'keyword_change1= plot_legend0,$\log_{10}(M_*$ $[M_{\odot}])<11.18$' 	>> $2
		echo 'keyword_change2= plot_legend1,$\log_{10}(M_*$ $[M_{\odot}])<11.24$' 	>> $2
		echo 'keyword_change3= plot_legend2,$\log_{10}(M_*$ $[M_{\odot}])<11.3$' 	>> $2
		echo 'keyword_change4= plot_legend3,$\log_{10}(M_*$ $[M_{\odot}])<11.35$' 	>> $2
		echo 'keyword_change5= plot_legend4,$\log_{10}(M_*$ $[M_{\odot}])<11.48$' 	>> $2
		echo 'keyword_change6= plot_legend5,Rodriguez+16 BigMD HAM' 	>> $2
		echo 'keyword_change7= plot_legend6,CMASS DR12' 	>> $2

	else if ($sign == 'bin1') then
		#1.5e11-2e11, 1.73e11-3e11, 2.24e11-4e11, 3e11-4e11, 1.05e11-3.45e11, 3.45e11-1e13
		echo 'keyword_change1= plot_legend0,$11.00 < \log_{10}(M_*$ $[M_{\odot}]) < 11.54$' 	>> $2
		echo 'keyword_change2= plot_legend1,$11.18 < \log_{10}(M_*$ $[M_{\odot}]) < 11.30$' 	>> $2
		echo 'keyword_change3= plot_legend2,$11.24 < \log_{10}(M_*$ $[M_{\odot}]) < 11.48$' 	>> $2
		#echo 'keyword_change3= plot_legend2,$11.35 < \log_{10}(M_*$ $[M_{\odot}]) < 11.60$' 	>> $2
		#echo 'keyword_change4= plot_legend3,$11.48 < \log_{10}(M_*$ $[M_{\odot}]) < 11.60$' 	>> $2

		echo 'keyword_change6= plot_legend5,$11.30 > \log_{10}(M_*$ $[M_{\odot}])$' 	>> $2
	else if ($sign == 'bin2') then
		#1.5e11-1.73e11, 1.73e11-2e11, 2e11-2.24e11, 2.24e11-3e11, 3e11-3.5e11, 3.5e11-4e11
		echo 'keyword_change1= plot_legend0,$11.18 < \log_{10}(M_*$ $[M_{\odot}]) < 11.24$' 	>> $2
		echo 'keyword_change2= plot_legend1,$11.24 < \log_{10}(M_*$ $[M_{\odot}]) < 11.3$' 	>> $2
		echo 'keyword_change3= plot_legend2,$11.3 < \log_{10}(M_*$ $[M_{\odot}]) < 11.35$' 	>> $2
		echo 'keyword_change4= plot_legend3,$11.35 < \log_{10}(M_*$ $[M_{\odot}]) < 11.48$' 	>> $2
		echo 'keyword_change5= plot_legend4,$11.48 < \log_{10}(M_*$ $[M_{\odot}]) < 11.54$' 	>> $2
		echo 'keyword_change6= plot_legend5,$11.54 < \log_{10}(M_*$ $[M_{\odot}]) < 11.6$' 	>> $2
		#echo 'keyword_change7= plot_legend6,$11.6 < \log_{10}(M_*$ $[M_{\odot}]) < 11.7$' 	>> $2

	else if ($scale == 'sm') then
		if ($use_cb_colors == 3) then

			echo 'use_loc_co_x= 0.57'			>> $2
			echo 'use_loc_co_y= 0.77'			>> $2
			#2.24e11, 3e11, 3.5e11, 4e11, 1e15
			echo 'keyword_change1= plot_legend0,$\log_{10}(M_*$ $[M_{\odot}])<11.18$' 	>> $2 
			echo 'keyword_change2= plot_legend1,$\log_{10}(M_*$ $[M_{\odot}])<11.35$' 	>> $2 
			#echo 'keyword_change1= plot_legend0,$\log_{10}(M_*$ $[M_{\odot}])'$sign'11.35$' 	>> $2
			#echo 'keyword_change2= plot_legend1,$\log_{10}(M_*$ $[M_{\odot}])'$sign'11.48$' 	>> $2
			#echo 'keyword_change3= plot_legend2,$\log_{10}(M_*$ $[M_{\odot}])'$sign'11.54$' 	>> $2
			#echo 'keyword_change3= plot_legend2,$\log_{10}(M_*$ $[M_{\odot}])'$sign'11.6$' 	>> $2
			echo 'keyword_change3= plot_legend2,all' 	>> $2
		else

			#1.5e11, 1.73e11, 2e11, 2.24e11, 3e11, 3.5e11, 4e11, 1e15 
			echo 'keyword_change1= plot_legend0,$\log_{10}(M_*$ $[M_{\odot}])'$sign'11.18$' 	>> $2
			echo 'keyword_change2= plot_legend1,$\log_{10}(M_*$ $[M_{\odot}])'$sign'11.24$' 	>> $2
			echo 'keyword_change3= plot_legend2,$\log_{10}(M_*$ $[M_{\odot}])'$sign'11.3$' 	>> $2
			echo 'keyword_change4= plot_legend3,$\log_{10}(M_*$ $[M_{\odot}])'$sign'11.35$' 	>> $2
			echo 'keyword_change5= plot_legend4,$\log_{10}(M_*$ $[M_{\odot}])'$sign'11.48$' 	>> $2
			echo 'keyword_change6= plot_legend5,$\log_{10}(M_*$ $[M_{\odot}])'$sign'11.54$' 	>> $2
			echo 'keyword_change7= plot_legend6,$\log_{10}(M_*$ $[M_{\odot}])'$sign'11.6$' 	>> $2
			echo 'keyword_change8= plot_legend7,all' 	>> $2
		endif
	endif

	#echo 'keyword_change1= plot_legend0,$-19<M_r<-18$' 	>> $2
	#echo 'keyword_change2= plot_legend1,$-20<M_r<-19$' 	>> $2
	#echo 'keyword_change3= plot_legend2,$-21<M_r<-20$' 	>> $2
	#echo 'keyword_change4= plot_legend3,$-22<M_r<-21$' 	>> $2
	#echo 'keyword_change5= plot_legend4,CUT1 $M_*$'  	>> $2
	#echo 'keyword_change6= plot_legend5,CUT2 $M_*$'  	>> $2
	#echo 'keyword_change7= plot_legend6,CUT3 $M_*$'  	>> $2
	#echo 'keyword_change8= plot_legend7,CUT1 $M_{Cold}$'  	>> $2
	#echo 'keyword_change9= plot_legend8,CUT2 $M_{Cold}$'  	>> $2
	#echo 'keyword_change10= plot_legend9,CUT3 $M_{Cold}$'  	>> $2
	#echo 'keyword_change11= plot_legend10,CUT1 SFR'  	>> $2
	#echo 'keyword_change12= plot_legend11,CUT2 SFR'  	>> $2
	#echo 'keyword_change13= plot_legend12,CUT3 SFR'  	>> $2

	#echo 'keyword_change8= plot_legend7,CUT3 $M_*$ cents'  	>> $2
	#echo 'keyword_change9= plot_legend8,CUT3 $M_*$ no'  	>> $2
	#echo 'keyword_change10= plot_legend9,CUT3 $M_{Cold}$ cents'  	>> $2
	#echo 'keyword_change11= plot_legend10,CUT3 $M_{Cold}$ no'  	>> $2
	#echo 'keyword_change12= plot_legend11,CUT3 SFR cents'  	>> $2
	#echo 'keyword_change13= plot_legend12,CUT3 SFR no'  	>> $2

	#Galacticus optimal
	#echo 'keyword_change1= plot_legend0,$M_*$ CUT2 no' 	>> $2
	#echo 'keyword_change2= plot_legend1,$M_{Cold}$ CUT2 no'  	>> $2
	#echo 'keyword_change3= plot_legend2,$M_{Cold}$ CUT3 no' 	>> $2
	#echo 'keyword_change4= plot_legend3,SFR CUT2 no' 	>> $2
	#echo 'keyword_change5= plot_legend4,SFR CUT1 no' 	>> $2
	#echo 'keyword_change6= plot_legend5,$M_r$ CUT3 no'  	>> $2

	#Galacticus Mcold Zcold-Subsamples
	#echo 'keyword_change1= plot_legend0,high $Z_{Cold}$-high $M_*$' 	>> $2
	#echo 'keyword_change2= plot_legend1,high $Z_{Cold}$-low $M_*$'  	>> $2
	#echo 'keyword_change3= plot_legend2,low $Z_{Cold}$-high $M_*$' 	>> $2
	#echo 'keyword_change4= plot_legend3,low $Z_{Cold}$-low $M_*$'	>> $2
	#echo 'keyword_change5= plot_legend4,all no' 	>> $2
	#echo 'keyword_change6= plot_legend5,all centrals'  	>> $2

	#SAG optimal
	#echo 'keyword_change1= plot_legend0,$M_*$ CUT3 no' 	>> $2
	#echo 'keyword_change2= plot_legend1,$M_*$ CUT3 centrals' 	>> $2
	#echo 'keyword_change3= plot_legend2,$M_{Cold}$ CUT3' 	>> $2
	#echo 'keyword_change4= plot_legend3,$M_*$ CUT1 centrals'  	>> $2
	#echo 'keyword_change5= plot_legend4,SFR CUT3' 	>> $2
	#echo 'keyword_change6= plot_legend5,$M_r$ CUT3'  	>> $2

	#SAGE optimal
	#echo 'keyword_change1= plot_legend0,$M_*$ CUT3' 	>> $2
	#echo 'keyword_change2= plot_legend1,$M_*$ CUT1 centrals' 	>> $2
	#echo 'keyword_change3= plot_legend2,$M_{Cold}$ CUT3' 	>> $2
	#echo 'keyword_change4= plot_legend3,SFR CUT3 centrals' 	>> $2
	#echo 'keyword_change5= plot_legend4,$M_*$ CUT3 centrals'  	>> $2
	#echo 'keyword_change6= plot_legend5,$M_r$ CUT3'  	>> $2

	#clustering test plots with Mr
	#echo 'keyword_change1= plot_legend0,test1' 	>> $2
	#echo 'keyword_change2= plot_legend1,test2' 	>> $2
	#echo 'keyword_change3= plot_legend2,test3' 	>> $2
	#echo 'keyword_change4= plot_legend3,test4' 	>> $2
	#echo 'keyword_change5= plot_legend4,test5'  	>> $2
	#echo 'keyword_change6= plot_legend5,$M_r$ CUT3'  	>> $2

endif

#DEFAULT configurations for contour plots
if ($1 =~ analyse_tarsel_*sfr_* || $1 =~ analyse_tarsel_*zcold*) then
	echo 'plot_type= hexbins'				>> $2
	echo 'contour_log= no'			>> $2
	echo 'colorbar_title= '       >> $2 
	echo 'error_bars_y= yes'				 >> $2
	echo 'error_bars_y_sub= yes'				 >> $2
	echo 'legend_ncols= 2'			>> $2

	if ($1 =~ *ssfr*) then
		echo 'y_range_min= -14'		 			 	 >> $2
		echo 'y_range_max= -10.4'				 		 >> $2
		echo 'float_format_y= 0f'         				>> $2
		echo 'yticks= -14,-13,-12,-11'		>> $2
		echo 'yticks_minor_offset= 0'				>> $2
	else
		echo 'y_range_min= -2'		 			 	 >> $2
		echo 'y_range_max= 1.1'				 		 >> $2
	endif

	echo 'plot_legend= yes'		 			>> $2

	echo 'hline_ypos= 2.35'					 >> $2

	echo 'use_loc_co_x= 0.48'				>> $2
	echo 'use_loc_co_y= 0.15'				>> $2
	echo 'z_print_position_x= 0.65'		>> $2
	echo 'z_print_position_y= 0.68'		>> $2

	echo 'xticks_minor_offset= 0'				>> $2

	if ($1 =~ *_ssfr_sfr*) then
		echo 'x_range_min= -2.7'			 		 >> $2
		echo 'x_range_max= 1.1'			 		 >> $2
		echo 'xticks= -2.5,-1.5,-0.5,0.5'		>> $2
		echo 'histo_num_density_max= 0.12'				>> $2

		echo 'use_loc_co_x= 0.57'				>> $2
		echo 'use_loc_co_y= 0.16'				>> $2

	else if ($1 =~ *_ssfr_mhal*) then
		echo 'x_range_min= 12'			 		 >> $2
		echo 'x_range_max= 15'			 		 >> $2
		echo 'xticks= 12,12.5,13,13.5,14,14.5'		>> $2
		echo 'histo_num_density_max= 0.12'				>> $2
		echo 'use_loc_co_x= 0.57'				>> $2
		echo 'use_loc_co_y= 0.66'				>> $2

	else if ($1 =~ *_ssfr_i) then
		echo 'x_range_min= 18.3'			 		 >> $2
		echo 'x_range_max= 21.5'			 		 >> $2
		echo 'xticks= 19,20,21'		>> $2
		echo 'histo_num_density_max= 0.12'				>> $2
		echo 'use_loc_co_x= 0.57'				>> $2
		echo 'use_loc_co_y= 0.66'				>> $2

	else if ($1 =~ *_ssfr_r) then
		echo 'x_range_min= 17'			 		 >> $2
		echo 'x_range_max= 22'			 		 >> $2
		#echo 'xticks= 19,20,21'		>> $2
		echo 'histo_num_density_max= 0.12'				>> $2
		echo 'use_loc_co_x= 0.57'				>> $2
		echo 'use_loc_co_y= 0.66'				>> $2

	else if ($1 =~ *_ssfr_zcold*) then
		echo 'x_range_min= 8.3'			 		 >> $2
		echo 'x_range_max= 10.5'			 		 >> $2
		#echo 'xticks= 12,12.5,13,13.5,14,14.5'		>> $2
		echo 'histo_num_density_max= 0.10'				>> $2
		echo 'use_loc_co_x= 0.16'				>> $2
		echo 'use_loc_co_y= 0.14'				>> $2
		echo 'y_title=  '	 	 	 	>> $2

	else if ($1 =~ *_ssfr_mcold*) then
		echo 'x_range_min= -4'			 		 >> $2
		echo 'x_range_max= -0.2'			 		 >> $2
		echo 'xticks= -4,-3,-2,-1'		>> $2
		echo 'histo_num_density_max= 0.12'				>> $2
		echo 'use_loc_co_x= 0.57'				>> $2
		echo 'use_loc_co_y= 0.66'				>> $2

	else if ($1 =~ *_zcold_mhal*) then
		echo 'x_range_min= 12'			 		 >> $2
		echo 'x_range_max= 15'			 		 >> $2
		echo 'xticks= 12,12.5,13,13.5,14,14.5'		>> $2
		echo 'histo_num_density_max= 0.12'				>> $2
		echo 'use_loc_co_x= 0.57'				>> $2
		echo 'use_loc_co_y= 0.66'				>> $2

		echo 'y_range_min= 8.3'		 			 	 >> $2
		echo 'y_range_max= 10.8'				 		 >> $2
		echo 'float_format_y= 1f'         				>> $2
		echo 'yticks= 8.5,9.0,9.5,10.0,10.5'		>> $2
		echo 'yticks_minor_offset= 0'				>> $2
	else if ($1 =~ *_zcold_*sfr*) then
		echo 'x_range_min= -3.5'			 		 >> $2
		echo 'x_range_max= 1.3'			 		 >> $2
		echo 'xticks= -3,-2,-1,0,1'		>> $2
		echo 'histo_num_density_max= 0.12'				>> $2
		echo 'use_loc_co_x= 0.57'				>> $2
		echo 'use_loc_co_y= 0.66'				>> $2

		echo 'y_range_min= 8.3'		 			 	 >> $2
		echo 'y_range_max= 10.8'				 		 >> $2
		echo 'float_format_y= 1f'         				>> $2
		echo 'yticks= 8.5,9.0,9.5,10.0,10.5'		>> $2
		echo 'yticks_minor_offset= 0'				>> $2
	else
		echo 'x_range_min= 10.65'			 		 >> $2
		echo 'x_range_max= 12.1'			 		 >> $2
		echo 'minor_ticks_x_space= 0.05'			 >> $2
		echo 'xticks= 11,11.5,12'		>> $2
	endif
	
	echo 'histo_num_density_max= 0.20'				>> $2
	echo 'histo_panels= yes'			>> $2
	#echo 'no_yticks= False'			>> $2
	echo 'minor_ticks_x_space= 0.1'			 >> $2
	echo 'minor_ticks_y_space= 0.1'			 >> $2

	#Violetas Paper:
	#----------------------	
	#ssfr = 0.3/t_H(z) Gyr-1
	#echo 'hline_ypos= -1.66'					 >> $2
	#ssfr = 1/t_H(z) Gyr-1
	#echo 'hline_ypos= -1.14'					 >> $2
	#echo 'y_range_min= -6'		 			 	 >> $2
	#echo 'y_range_max= 1'				 		 >> $2
	#echo 'no_last_xticks= True'         			>> $2
	#echo 'no_last_yticks= True'         			>> $2

else if ($1 =~ analyse_tarsel_*g-*_ssfr) then
	echo 'plot_type= hexbins'				>> $2
	echo 'contour_log= no'			>> $2
	echo 'colorbar_title= '       >> $2
	echo 'histo_num_density_max= 0.12'				>> $2
	echo 'y_range_min= 0.95'			 		 >> $2
	echo 'y_range_max= 1.85'			 		 >> $2

	echo 'x_range_min= -13.5'		 			 	 >> $2
	echo 'x_range_max= -10.0'				 		 >> $2

	echo 'plot_legend= yes'		 			>> $2

	#echo 'hline_ypos= -11'					 >> $2
	echo 'filled_between= False'					 >> $2

	echo 'histo_panels= yes'			>> $2
	echo 'use_loc_co_x= 0.18'				>> $2
	echo 'use_loc_co_y= 0.15'				>> $2
	echo 'z_print_position_x= 0.18'		>> $2
	echo 'z_print_position_y= 0.75'		>> $2

	echo 'minor_ticks_x_space= 0.1'			 	>> $2
	echo 'minor_ticks_y_space= 0.05'			 	>> $2

	echo 'xticks= -13.5,-12.5,-11.5,-10.5'		>> $2
	echo 'xticks_minor_offset= 0'				>> $2

	#echo 'no_last_xticks= True'         			>> $2
	echo 'no_last_yticks= True'         			>> $2

else if ($1 =~ analyse_tarsel_*dmesa*_ssfr) then
	echo 'plot_type= hexbins'				>> $2
	echo 'contour_log= no'			>> $2
	echo 'colorbar_title= '       >> $2
	echo 'histo_num_density_max= 0.15'				>> $2
	echo 'y_range_min= 0.55'			 		 >> $2
	echo 'y_range_max= 0.93'			 		 >> $2
	echo 'float_format_y= 2f'         				>> $2
	echo 'x_range_min= -13.5'		 			 	 >> $2
	echo 'x_range_max= -10'				 		 >> $2

	echo 'plot_legend= yes'		 			>> $2

	#echo 'hline_ypos= -11'					 >> $2
	echo 'filled_between= False'					 >> $2

	echo 'histo_panels= yes'			>> $2
	echo 'use_loc_co_x= 0.16'				>> $2
	echo 'use_loc_co_y= 0.14'				>> $2
	echo 'z_print_position_x= 0.18'		>> $2
	echo 'z_print_position_y= 0.73'		>> $2

	echo 'minor_ticks_x_space= 0.1'			 	>> $2
	echo 'minor_ticks_y_space= 0.01'			 	>> $2

	echo 'xticks= -13.5,-12.5,-11.5,-10.5'		>> $2
	echo 'xticks_minor_offset= 0'				>> $2
	echo 'yticks= 0.6,0.7,0.8,0.9'		>> $2
	echo 'yticks_minor_offset= 0'				>> 
	#echo 'no_last_xticks= True'         			>> $2
	echo 'no_last_yticks= True'         			>> $2

else if ($1 =~ analyse_tarsel_*rdisk* || $1 =~ analyse_tarsel_*rhalf*) then
	echo 'plot_type= hexbins'				>> $2
	echo 'contour_log= no'			>> $2
	#echo 'contour_log_wich_axis= both'		>> $2
	echo 'colorbar_title= '       >> $2

	set key = 'half'

	if ($key == 'half') then
		echo 'minor_ticks_x_space= 0.05'			 	>> $2
		echo 'minor_ticks_y_space= 0.2'			 	>> $2

		echo 'x_range_min= 10.5'			 		 >> $2
		echo 'x_range_max= 12.2'			 		 >> $2

		echo 'y_range_min= -2'		 			 	 >> $2
		echo 'y_range_max= 4'				 		 >> $2
	else
		echo 'minor_ticks_x_space= 0.1'			 	>> $2
		echo 'minor_ticks_y_space= 0.1'			 	>> $2
		echo 'x_range_min= -2'			 		 >> $2
		echo 'x_range_max= 3'			 		 >> $2

		echo 'y_range_min= -2'		 			 	 >> $2
		echo 'y_range_max= 3'				 		 >> $2
	endif

	echo 'plot_legend= no'		 			>> $2
	echo 'use_loc_co_x= 0.9'			>> $2
	echo 'use_loc_co_y= 0.17'			>> $2

	echo 'z_print_position_x= 0.17'		>> $2
	echo 'z_print_position_y= 0.9'		>> $2

else if ($1 =~ analyse_tarsel_mstar_mhalo) then
	echo 'plot_type= hexbins'				>> $2
	echo 'contour_log= yes'			>> $2
	echo 'contour_log_wich_axis= both'		>> $2
	echo 'colorbar_title= '       >> $2 
	echo 'x_range_min= 12'			 		 >> $2
	echo 'x_range_max= 15'			 		 >> $2
	echo 'y_range_min= 10.6'		 			 	 >> $2
	echo 'y_range_max= 12.0'				 		 >> $2
	echo 'plot_legend= yes'		 			>> $2
	echo 'histo_num_density_max= 0.15'				>> $2
	echo 'xticks= 12,12.5,13,13.5,14,14.5'		>> $2
	
	echo 'yticks= 10.6,11.0,11.4,11.8'		>> $2
	echo 'yticks_minor_offset= 0'				>> $2

	echo 'histo_panels= yes'			>> $2
	echo 'use_loc_co_x= 0.57'				>> $2
	echo 'use_loc_co_y= 0.16'				>> $2
	echo 'z_print_position_x= 0.18'		>> $2
	echo 'z_print_position_y= 0.73'		>> $2

	echo 'no_yticks= False'					>> $2
	echo 'minor_ticks_x_space= 0.1'			 >> $2
	echo 'minor_ticks_y_space= 0.05'			 >> $2


else if ($1 =~ analyse_tarsel_mbh_mstar*) then
	echo 'plot_type= hexbins'				>> $2
	echo 'contour_log= yes'			>> $2
	echo 'contour_log_wich_axis= both'		>> $2
	echo 'colorbar_title= '       >> $2 
	echo 'x_range_min= 10.6'			 		 >> $2
	echo 'x_range_max= 12.2'			 		 >> $2
	echo 'y_range_min= 7'		 			 	 >> $2
	echo 'y_range_max= 9.75'				 		 >> $2
	echo 'plot_legend= yes'		 			>> $2
	echo 'use_loc_co_x= 0.55'			>> $2
	echo 'use_loc_co_y= 0.16'			>> $2								

	echo 'xticks= 10.6,11.0,11.4,11.8'		>> $2
	echo 'xticks_minor_offset= 0'				>> $2
	echo 'histo_panels= yes'			>> $2
	echo 'minor_ticks_x_space= 0.1'			 	>> $2
	echo 'minor_ticks_y_space= 0.25'			 	>> $2
	echo 'no_yticks= False'					>> $2



else if ($1 =~ *zgas*mcold) then
	if ($1 =~ *tarsel*) then
		echo 'HERE: 547'
		echo 'plot_type= hexbins'				>> $2
		echo 'contour_log= yes'			>> $2
		echo 'contour_log_wich_axis= x'		>> $2
		echo 'plot_legend= yes'		 			>> $2
		#echo 'y_title= False'	         				>> $2
		echo 'use_loc_co_x= 0.55'			>> $2
		echo 'use_loc_co_y= 0.66'			>> $2								

		#echo 'xticks= 10.6,11.0,11.4,11.8'		>> $2
		echo 'xticks_minor_offset= 0'				>> $2
		echo 'histo_panels= yes'			>> $2
		echo 'minor_ticks_x_space= 0.1'			 	>> $2
		echo 'minor_ticks_y_space= 0.1'			 	>> $2
		echo 'no_yticks= False'					>> $2

	else
		echo 'x_title= $\log_{10}$ ($M_{Cold}/M_{_*}$)' 	        >> $2	
		echo 'y_title= $Z_{Cold}$'       >> $2
		echo 'y_title= '       >> $2
		echo 'plot_legend= yes'		 			>> $2
		echo 'add_which_axis_sub= x'				>> $2
		echo 'add_which_col_sub= 6'					>> $2
		echo 'log_scale_add_sub= no'				>> $2 
		echo 'add_axis_sub= yes'					>> $2

		echo 'title_add0= $\log_{10}$ ($M_{_*}$ $[M_{\odot}]$) Galacticus'	>> $2
		echo 'title_add1= $\log_{10}$ ($M_{_*}$ $[M_{\odot}]$) SAG'	>> $2
		echo 'title_add2= $\log_{10}$ ($M_{_*}$ $[M_{\odot}]$) SAGE'	>> $2
	endif

	#echo 'keyword_change1= plot_legend0,$8.7 < \log M_{_*} < 9.0$' 	>> $2
	echo 'colorbar_title= '       >> $2 
	echo 'x_range_min= -4.0'			 		 >> $2
	echo 'x_range_max= -0.25'			 		 >> $2
	echo 'y_range_min= 8'		 			 	 >> $2
	echo 'y_range_max= 10.75'				 		 >> $2
	echo 'no_last_xticks= True'					>> $2
	#echo 'vline_xpos= 0.0'					 >> $2
	#echo 'hline_ypos= 8.69'					 >> $2
	echo 'histo_num_density_max= 0.15'				>> $2

else if ($1 =~ analyse_tarsel_Mzgas*mstar) then
	echo 'plot_type= hexbins'				>> $2
	echo 'contour_log= yes'			>> $2
	echo 'contour_log_wich_axis= x'		>> $2
	echo 'colorbar_title= '       >> $2 
	echo 'x_range_min= 10.6'			 		 >> $2
	echo 'x_range_max= 12.2'			 		 >> $2
	echo 'y_range_min= 8'		 			 	 >> $2
	echo 'y_range_max= 10.75'				 		 >> $2
	echo 'plot_legend= yes'		 			>> $2
	echo 'use_loc_co_x= 0.55'			>> $2
	echo 'use_loc_co_y= 0.66'			>> $2								

	echo 'xticks= 10.6,11.0,11.4,11.8'		>> $2
	echo 'xticks_minor_offset= 0'				>> $2
	echo 'histo_panels= yes'			>> $2
	echo 'minor_ticks_x_space= 0.1'			 	>> $2
	echo 'minor_ticks_y_space= 0.25'			 	>> $2
	echo 'no_yticks= False'					>> $2
	echo 'minor_ticks_x_space= 0.1'			 >> $2
	echo 'minor_ticks_y_space= 0.1'			 >> $2

else if ($1 =~ analyse_tarsel_ang*mbar) then
	echo 'plot_type= hexbins'				>> $2
	echo 'contour_log= yes'			>> $2
	echo 'contour_log_wich_axis= both'		>> $2
	echo 'colorbar_title= '       >> $2 
	echo 'x_range_min= 9.4'			 		 >> $2
	echo 'x_range_max= 11.6'			 		 >> $2
	echo 'y_range_min= 1'		 			 	 >> $2
	echo 'y_range_max= 4.5'				 		 >> $2
	echo 'plot_legend= yes'		 			>> $2
	echo 'use_loc_co_x= 0.55'			>> $2
	echo 'use_loc_co_y= 0.14'			>> $2								

	echo 'xticks= 9.4,9.8,10.2,10.6,11.0,11.4'		>> $2
	echo 'xticks_minor_offset= 0'				>> $2
	echo 'histo_panels= yes'			>> $2
	echo 'minor_ticks_x_space= 0.1'			 	>> $2
	echo 'minor_ticks_y_space= 0.25'			 	>> $2
	echo 'no_yticks= False'					>> $2
	#echo 'no_last_xticks= True'					>> $2
	echo 'no_last_yticks= True'					>> $2
	echo 'minor_ticks_x_space= 0.1'			 >> $2
	echo 'minor_ticks_y_space= 0.1'			 >> $2

	echo 'histo_num_density_max= 0.12'				>> $2

else if ($1 =~ analyse_tarsel_bdisk*) then
	echo 'plot_type= hexbins'				>> $2
	echo 'contour_log= no'			>> $2
	echo 'colorbar_title= '       >> $2 
	echo 'x_range_min= -5.25'			 		 >> $2
	echo 'x_range_max= -2.75'			 		 >> $2
	echo 'y_range_min= -5'		 			 	 >> $2
	echo 'y_range_max= -1'				 		 >> $2
	echo 'plot_legend= no'		 			>> $2
	echo 'use_loc_co_x= 0.55'			>> $2
	echo 'use_loc_co_y= 0.14'			>> $2								

	echo 'format_xticks= False'	         			>> $2
	#echo 'xticks= 1,2,3,4,5'		>> $2
	echo 'xticks_minor_offset= 0'				>> $2
	echo 'histo_panels= yes'			>> $2
	echo 'minor_ticks_x_space= 0.25'			 	>> $2
	echo 'minor_ticks_y_space= 0.25'			 	>> $2
	echo 'no_yticks= False'					>> $2
	echo 'no_last_xticks= True'					>> $2
	echo 'no_last_yticks= True'					>> $2
	echo 'minor_ticks_x_space= 0.1'			 >> $2
	echo 'minor_ticks_y_space= 0.1'			 >> $2

	echo 'histo_num_density_max= 0.12'				>> $2
								
else if ($1 =~ analyse_tarsel_mcold*) then
	echo 'plot_type= hexbins'				>> $2
	echo 'contour_log= yes'			>> $2
	echo 'contour_log_wich_axis= both'		>> $2
	echo 'histo_panels= yes'			>> $2
	echo 'colorbar_title= '       >> $2 
	echo 'x_range_min= 8.5'			 		 >> $2
	echo 'x_range_max= 11.5'			 		 >> $2
	echo 'y_range_min= -3'		 			 	 >> $2
	echo 'y_range_max= 2'				 		 >> $2
	echo 'plot_legend= yes'		 			>> $2
	echo 'use_loc_co_x= 0.42'			>> $2
	echo 'use_loc_co_y= 0.68'			>> $2								

	#echo 'xticks= 10.6,11.0,11.4,11.8'		>> $2
	echo 'xticks_minor_offset= 0'				>> $2

	echo 'minor_ticks_x_space= 0.1'			 	>> $2
	echo 'minor_ticks_y_space= 0.25'			 	>> $2
	echo 'no_yticks= False'					>> $2

	echo 'no_last_yticks= True'					>> $2
	echo 'no_last_xticks= True'					>> $2
	echo 'histo_num_density_max= 0.10'				>> $2

	echo 'z_print_position_x= 0.18'		>> $2
	echo 'z_print_position_y= 0.16'		>> $2

	echo 'error_bars_y_sub= yes'					 >> $2

#contour u-r vs r
else if ($1 =~ *u-*r) then

	echo 'plot_type= hexbins'				>> $2
	echo 'contour_log= no'		>> $2
	echo 'x_range_min= -22.5'		 			 >> $2
	echo 'x_range_max= -18'			 		 >> $2
	echo 'y_range_min= 0.5'		 			 >> $2
	echo 'y_range_max= 5'				 		 >> $2
	echo 'z_print_position_x= 0.62'		>> $2
	echo 'z_print_position_y= 0.2'		>> $2
	echo 'format_xticks= False'	         			>> $2
	echo 'plot_legend= yes'		 			>> $2
	echo 'filled_between= False'					 >> $2

	echo 'histo_panels= yes'			>> $2
	echo 'use_loc_co_x= 0.55'				>> $2
	echo 'use_loc_co_y= 0.15'				>> $2
	echo 'z_print_position_x= 0.18'		>> $2
	echo 'z_print_position_y= 0.73'		>> $2

	#echo 'no_last_xticks= True'         			>> $2
	echo 'no_last_yticks= True'         			>> $2

	#echo 'xticks= -23.5,-22.5,-21.5,-20.5,-19.5,-18.5'		>> $2
	echo 'xticks= -22.5,-21.5,-20.5,-19.5,-18.5'		>> $2
	echo 'xticks_minor= 0.25'					>> $2
	echo 'xticks_minor_offset= 0'				>> $2
	echo 'yticks_minor= 0.25'					>> $2
	echo 'histo_num_density_max= 0.18'				>> $2
	echo 'minor_ticks_x_space= 0.1'			 	>> $2
	echo 'minor_ticks_y_space= 0.1'			 	>> $2

#contour r-i vs i
else if ($1 =~ *-i_i*) then
	echo 'tight_layout= False'		 		 	>> $2
	echo 'plot_type= hexbins'				>> $2
	echo 'contour_log= no'		>> $2
	echo 'minor_ticks_y_space= 0.05'			 	>> $2
	set mag = 'M'
	set band = 'i'
	if ($mag == M) then
		echo 'x_range_min= -23.3'		 			 >> $2
		echo 'x_range_max= -19.6'			 		 >> $2
		echo 'xticks= -23,-22,-21,-20'		>> $2
	else
		echo 'x_range_min= 23.0'		 			 >> $2
		echo 'x_range_max= 17.0'		 		 >> $2

	endif

	if ($band == g) then
		echo 'y_range_min= 0.5'		 			 >> $2
		echo 'y_range_max= 3.5'				 	 >> $2
	else
		echo 'y_range_min= 0.65'		 			 >> $2
		echo 'y_range_max= 1.15'				 	 >> $2
	endif
	echo 'histo_num_density_max= 0.18'				>> $2
	echo 'histo_panels= yes'				>> $2
	echo 'error_bars_y_sub= no'					 >> $2
	echo 'z_print_position_x= 0.16'				>> $2
	echo 'z_print_position_y= 0.9'				>> $2
	echo 'use_loc_co_x= 0.16'				>> $2
	echo 'use_loc_co_y= 0.15'				>> $2
	echo 'plot_legend= yes'		 			>> $2
	echo 'no_last_xticks= False'         			>> $2
	echo 'no_last_yticks= True'         			>> $2
#contour r-i vs mstar
else if ($1 == analyse_tarsel_r-i_mstar || $1 == analyse_tarsel_g-i_mstar || $1 == analyse_tarsel_g-i_mhalo) then
	echo 'plot_type= hexbins'				>> $2
	echo 'tight_layout= False'		 		 	>> $2
	echo 'contour_log= no'			>> $2
	echo 'contour_log_wich_axis= x'		>> $2

	if ($1 =~ analyse_tarsel_*mhalo) then
		echo 'x_range_min= 12.0'		 			 >> $2
		echo 'x_range_max= 14.9'			 		 >> $2

		echo 'xticks= 12.5,13.5,14.5'		>> $2
	else
		echo 'x_range_min= 10.0'		 			 >> $2
		echo 'x_range_max= 11.7'			 		 >> $2

		echo 'xticks= 10.0,10.5,11.0,11.5'		>> $2
	endif

	echo 'plot_legend= yes'		 			>> $2
	echo 'histo_num_density_max= 0.2'				>> $2
	echo 'legend_sep_xpos= 0.55'	>> $2
	
	set band = 'g'
	if ($band == g) then
		echo 'y_range_min= 1.8'		 			 >> $2
		echo 'y_range_max= 2.9'				 	 >> $2

		echo 'yticks= 1.8,2.0,2.2,2.4,2.6,2.8'		>> $2
		echo 'hline_ypos= 2.35'					 >> $2
	else
		echo 'y_range_min= 0.0'		 			 >> $2
		echo 'y_range_max= 1.2'				 	 >> $2
	endif

	set histo_panels = 'yes'
	echo 'histo_panels= '$histo_panels				>> $2
	if ($histo_panels == 'yes') then
		echo 'use_loc_co_x= 0.6'				>> $2
		echo 'use_loc_co_y= 0.15'				>> $2
		echo 'z_print_position_x= 0.62'		>> $2
		echo 'z_print_position_y= 0.38'		>> $2
	else
		echo 'use_loc_co_x= 0.60'				>> $2
		echo 'use_loc_co_y= 0.85'				>> $2
		echo 'z_print_position_x= 0.16'		>> $2
		echo 'z_print_position_y= 0.9'		>> $2
	endif	

	echo 'yticks_minor= 0.01'				>> $2
	echo 'colorbar_anchor_right= 0.9'		 		 >> $2
	echo 'no_last_xticks= False'         			>> $2
	echo 'no_last_yticks= False'         			>> $2

#contour g-r vs r-i
else if ($1 =~ *r-i*g-*) then
	echo 'contour_log= no'		>> $2
	echo 'no_last_yticks= False'         			>> $2

	echo 'x_range_min= 0.0'		 			 >> $2
	echo 'x_range_max= 1'			 		 >> $2
	echo 'y_range_min= 0.0'		 			 >> $2
	echo 'y_range_max= 0.8'				 		 >> $2
	echo 'histo_num_density_max= 0.15'				>> $2

	echo 'legend_sep_xpos= 0.16'	>> $2
	echo 'legend_sep_ypos= 0.63'	>> $2

	set zoom = 'no'
	if ($zoom == 'yes') then
		echo 'histo_num_density_max= 0.12'				>> $2
		echo 'float_format_y= 2f'         				>> $2
		echo 'x_range_min= 1.0'		 			 >> $2
		echo 'x_range_max= 1.8'			 		 >> $2
		echo 'y_range_min= 0.8'		 			 >> $2
		echo 'y_range_max= 1.15'				 	 >> $2
	endif

	echo 'minor_ticks_y_space= 0.1'			 	>> $2
	echo 'error_bars_y_sub= no'					 >> $2
	echo 'colorbar_anchor_right= 0.9'		 		 >> $2
	echo 'plot_legend= no'		 			>> $2

	set histo_panels = 'yes'
	echo 'histo_panels= '$histo_panels				>> $2
	if ($histo_panels == 'yes') then
		echo 'use_loc_co_x= 0.53'				>> $2
		echo 'use_loc_co_y= 0.15'				>> $2
		echo 'z_print_position_x= 0.16'		>> $2
		echo 'z_print_position_y= 0.75'		>> $2
	else
		echo 'use_loc_co_x= 0.60'				>> $2
		echo 'use_loc_co_y= 0.85'				>> $2
		echo 'z_print_position_x= 0.16'		>> $2
		echo 'z_print_position_y= 0.9'		>> $2
	endif
	echo 'no_last_xticks= False'         			>> $2
	echo 'no_last_yticks= True'         			>> $2

	#echo 'xticks= 0,0.5,1,1.5,2'		>> $2

#contour dmesa vs i-band
else if ($1 == analyse_tarsel_dmesa_i) then
	echo 'contour_log= no'					>> $2 
	echo 'plot_legend= yes'		 			>> $2
	echo 'minor_ticks_x_space= 0.1'			 >> $2
	echo 'minor_ticks_y_space= 0.05'			 >> $2
	echo 'histo_num_density_max= 0.18'				>> $2
	echo 'legend_sep_xpos= 0.57'	>> $2

	set zoom = 'no'
	if ($zoom == 'yes') then
		#zoomed in
		echo 'x_range_min= 18'		 	>> $2
		echo 'x_range_max= 21.5'			 >> $2
		echo 'x_range_min= 18'		 	>> $2
		echo 'x_range_max= 21.7'			 >> $2
		echo 'y_range_min= 0.45'		 		>> $2
		echo 'y_range_max= 1.0'				 >> $2
		echo 'z_print_position_x= 0.16'		>> $2
	else
		echo 'x_range_min= 17.5'		 	>> $2
		echo 'x_range_max= 25'			 >> $2
		echo 'y_range_min= 0.45'		 	>> $2
		echo 'y_range_max= 1.1'				 >> $2
		echo 'x_range_min= 17.5'		 	>> $2
		echo 'x_range_max= 25'			 >> $2
		echo 'y_range_min= 0.2'		 	>> $2
		echo 'y_range_max= 1.1'				 >> $2
	endif

	echo 'use_loc_co_x= 0.55'				>> $2
	echo 'use_loc_co_y= 0.65'				>> $2
	echo 'error_bars_y_sub= no'					 >> $2

	echo 'xticks= 18,19,20,21,22,23,24'		>> $2

	set histo_panels = 'yes'
	if ($histo_panels == 'yes') then
		echo 'histo_panels= '$histo_panels			>> $2
		echo 'use_loc_co_x= 0.49'				>> $2
		echo 'use_loc_co_y= 0.65'				>> $2
		echo 'z_print_position_x= 0.18'		>> $2
		echo 'z_print_position_y= 0.75'		>> $2
	else
		echo 'use_loc_co_x= 0.60'				>> $2
		echo 'use_loc_co_y= 0.85'				>> $2
		echo 'z_print_position_x= 0.55'		>> $2
		echo 'z_print_position_y= 0.2'		>> $2
	endif	

	echo 'no_last_yticks= True'         			>> $2

#contour dmesa vs mstar
else if ($1 == analyse_tarsel_dmesa_mstar) then
	echo 'plot_type= hexbins'				>> $2
	echo 'contour_log= yes'					>> $2
	echo 'contour_log_wich_axis= x'		>> $2 2
	echo 'plot_legend= yes'		 			>> $2
	echo 'minor_ticks_x_space= 0.05'			 >> $2
	echo 'minor_ticks_y_space= 0.05'			 >> $2
	echo 'histo_panels= yes'				>> $2
	echo 'z_print_position_y= 0.9'		>> $2
	set zoom = 'no'
	if ($zoom == 'yes') then
		#zoomed in
		echo 'x_range_min= 11'		 	>> $2
		echo 'x_range_max= 12.5'			 >> $2
		echo 'y_range_min= 0.50'		 		>> $2
		echo 'y_range_max= 1.2'				 >> $2
		echo 'z_print_position_x= 0.16'		>> $2
	else
		echo 'x_range_min= 11'		 	>> $2
		echo 'x_range_max= 12.5'			 >> $2
		echo 'y_range_min= 0.55'		 	>> $2
		echo 'y_range_max= 1.0'				 >> $2
	endif

	echo 'use_loc_co_x= 0.16'				>> $2
	echo 'use_loc_co_y= 0.14'				>> $2
	echo 'z_print_position_x= 0.16'		>> $2
	echo 'z_print_position_y= 0.75'		>> $2

	echo 'error_bars_y_sub= no'					 >> $2
	echo 'no_last_yticks= True'         			>> $2
	echo 'no_last_xticks= True'         			>> $2

#histo
else if ($1 == analyse_tarsel_histo) then
	echo 'plot_type= barhisto'				>> $2
	echo 'x_range_min= 10.8'			 		 >> $2
	echo 'x_range_max= 12.2'			 		 >> $2
	echo 'y_range_min= 0'		 			 	 >> $2
	echo 'y_range_max= 0.35'				 		 >> $2
	echo 'legend_position= upper right'			 	 >> $2
	echo 'keyword_change1= plot_legend0,BOSS CMASS 0.50<z<0.60'	 >> $2
	echo 'keyword_change2= plot_legend1,Galacticus CMASS  z=0.56' 	 >> $2

else if ($1 =~ analyse*) then
	echo 'plot_type= hexbins'				>> $2
	echo 'contour_log= yes'					>> $2
	echo 'contour_log_wich_axis= x'		>> $2 
	echo 'plot_legend= yes'		 			>> $2
	echo 'minor_ticks_x_space= 0.1'			 >> $2
	echo 'minor_ticks_y_space= 0.1'			 >> $2
	echo 'histo_panels= yes'				>> $2
	echo 'z_print_position_y= 0.9'		>> $2

	echo 'x_range_min= 10.5'		 	>> $2
	echo 'x_range_max= 12.2'			 >> $2
	echo 'y_range_min= 17.5'		 	>> $2
	echo 'y_range_max= 22'				 >> $2

	echo 'use_loc_co_x= 0.16'				>> $2
	echo 'use_loc_co_y= 0.14'				>> $2
	echo 'z_print_position_x= 0.16'		>> $2
	echo 'z_print_position_y= 0.75'		>> $2

	echo 'xticks= 10.5,11.0,11.5,12.0'		>> $2
	#echo 'xticks= -23,-22,-21,-20'		>> $2
	echo 'xticks_minor= 0.25'					>> $2
	echo 'xticks_minor_offset= 0'				>> $2

	echo 'error_bars_y_sub= no'					 >> $2
	echo 'no_last_yticks= True'         			>> $2
	echo 'no_last_xticks= True'         			>> $2
								
endif
#####################################################################################################################################################
#####################################################################################################################################################

if ($1 == plotOnly) then

	set mylog_scale = 'no'
	set x_range_min = 'min'
	set x_range_max = 'max'
	set y_range_min = -8
	set y_range_max = 0

	if ($5 == mstar) then
		echo $5
		set mylog_scale = 'yes'
		set x_range_min = '1e10'
		set x_range_max = '1e13'
		set y_range_min = '1e-8'
		set y_range_max = '0.1'
	else if ($5 == mAB_dA_total_g) then
		echo $5
		set mylog_scale = 'no'
		set x_range_min = '18.5'
		set x_range_max = '23'
	else if ($5 == mAB_dA_total_r) then
		echo $5
		set mylog_scale = 'no'
		set x_range_min = '18'
		set x_range_max = '21'
	else if ($5 =~ mAB*total_i) then
		echo $5
		set mylog_scale = 'no'
		set x_range_min = -10.5
		set x_range_max = 7
	else if ($5 =~ MAB*total_i) then
		echo $5
		set mylog_scale = 'no'
		set x_range_min = -40
		set x_range_max = -15
	else if ($5 == mAB_dA_total_z) then
		echo $5
		set mylog_scale = 'no'
		set x_range_min = '16'
		set x_range_max = '20'
	endif

	echo 'title= MDPL2' >> $2
	echo 'x_title= $M_{_*}$ [$M_{\odot}$]' 	         >> $2
	#echo 'x_title= $MAB_{i}$' 	         >> $2
	echo 'y_title=  $\log_{10}$ ($\Phi$ [$Mpc^{-3}$])'       >> $2  
	echo 'x_range_min= '$x_range_min			 			 >> $2
	echo 'x_range_max= '$x_range_max		 		 >> $2
	echo 'y_range_min= '$y_range_min		 			 >> $2
	echo 'y_range_max= '$y_range_max				 		 >> $2
	echo 'legend_position= upper right'			 	 >> $2
	echo 'add_axis= no'						 >> $2
	echo 'add_which_axis= no'					>> $2
	echo 'title_add= $number$ $count$ '$6		 	 >> $2
	echo 'range_min_add= min'					 >> $2
	echo 'range_max_add= max'					 >> $2
	echo 'title_add2= $number$ $count$ '$7	 	 >> $2

	echo 'linestyle= --'			 			 >> $2
	echo 'linestyle_change= '		 			 >> $2
	echo 'linestyle_sub= --'	 			 >> $2
	echo 'linestyle_sub_change= -.'		 			 >> $2
	echo 'keyword_change1= plot_legend0,'$6' z=0.56 Mstar>1.2e11h'		 >> $2
	echo 'keyword_change2= plot_legend1,Galacticus z=0.56 BOSS-CMASS'		 			>> $2
	echo 'keyword_change3= default'		 			>> $2
	echo 'log_scale_x= '$mylog_scale				 >> $2
	echo 'log_scale_y= yes'						 >> $2
	echo 'log_scale_x_sub= '$mylog_scale					 >> $2
	echo 'log_scale_y_sub= yes'					 >> $2
	echo 'error_bars_x= no'					 	 >> $2
	echo 'error_bars_y= yes'					 >> $2
	echo 'error_bars_x_sub= no'					 >> $2
	echo 'error_bars_y_sub= yes'					>> $2

								
endif


if ($1 == SMF || $1 == SFRF || $1 == sSFRF) then

	set use_cb_colors = 2
	echo 'use_cb_colors= 2'	>> $2
	echo 'error_bars_x= no'					 	 >> $2
	echo 'error_bars_y= yes'					 >> $2
	echo 'error_bars_x_sub= no'					 >> $2
	echo 'error_bars_y_sub= yes'					 >> $2

	set plot_CUT = 'CMASS'
	echo 'plot_CUT= '$plot_CUT >> $2
	#echo 'print_redshift= z=3.0' >> $2
	#echo 'hline_ypos= -2.01,-2.51,-3.01'					 >> $2

	echo 'keyword_change1= plot_legend0,Gal-MD1000'	 		 >> $2
	echo 'keyword_change2= plot_legend1,Gal-SMD400'	 		 >> $2
	echo 'keyword_change3= plot_legend2,Ill-TNG300'			 >> $2


	if ($plot_CUT == 'mcold') then
		echo 'x_title= $\log_{10}$ ($M_{Cold}$ [$M_{\odot}$])'	         >> $2
		echo 'y_title= '  >> $2
		echo 'x_range_min= 7.5'			 		 	 >> $2 
		echo 'x_range_max= 12.5'			 		 >> $2
		echo 'y_range_min= -7.2'		 			 >> $2
		echo 'y_range_max= 0'				 		 >> $2
		#echo 'hline_ypos= -1.93'					 >> $2
		echo 'tight_layout= False'		 		 	>> $2

	else if ($plot_CUT == 'sfr') then
		echo 'x_title= $\log_{10}$ (SFR [$M_{\odot} yr^{-1}$])' 	 >> $2
		echo 'y_title= $\log_{10}$ ($n_{gal}$ (> SFR $[M_{\odot} yr^{-1}$]) [$Mpc^{-3}$])'  >> $2
		echo 'x_range_min= -5'			 		 >> $2 
		echo 'x_range_max= 3'   					 >> $2
		echo 'y_range_min= -7.2'		 			 >> $2
		echo 'y_range_max= -1'				 		 >> $2

	else if ($plot_CUT == 'ssfr') then
		echo 'x_title= $\log_{10}$ (sSFR [$yr^{-1}$])' 	 >> $2
		echo 'y_title= $\log_{10}$ ($n_{gal}$ (> sSFR [$yr^{-1}$]) [$Mpc^{-3}$])'  >> $2
		echo 'x_range_min= -15'			 		 >> $2 
		echo 'x_range_max= -7'   					 >> $2
		echo 'y_range_min= -7.2'		 			 >> $2
		echo 'y_range_max= -1'				 		 >> $2

	else if ($plot_CUT == 'mstar') then
		echo 'x_title= $\log_{10}$ ($M_{\star}$ [$M_{\odot}$])'	         >> $2
		echo 'y_title= $\log_{10}$ ($n_{gal}$ (> X) [$Mpc^{-3}$])'  >> $2
		echo 'x_range_min= 8.7'			 		 >> $2 
		echo 'x_range_max= 13.2'			 		 >> $2
		echo 'xticks= 9,10,11,12,13'
		echo 'y_range_min= -7.2'		 			 >> $2
		echo 'y_range_max= -1'				 		 >> $2
		echo 'tight_layout= False'		 		 	>> $2


	else if ($plot_CUT == 'CMASS') then
		echo 'print_redshift= z$\sim0.55$' >> $2
		echo 'legend_fontsize= 30'					>> $2
		echo 'adjust_right= 0.98'					>> $2
		echo 'adjust_top= 0.98'					>> $2
		#echo 'set_thin_line= 3'				 	>> $2			 	 
		echo 'x_range_min= 10.7'			 			 >> $2
		echo 'x_range_max= 12.3'			 		 	 >> $2
		echo 'y_range_min= -7.2'			 		 	 >> $2
		echo 'y_range_max= -3'			 		 >> $2
		echo 'x_title= $\log_{10}$ ($M_{\star}$ [$M_{\odot}$])'	         >> $2
		echo 'y_title= $\log_{10}$ ($\Phi$ [$Mpc^{-3}$ $dex^{-1}$])  '       >> $2
		echo 'size_x= 12'				 		 >> $2
		echo 'size_y= 8'				 		 >> $2

		echo 'keyword_change1= plot_legend0,Gal-dens'	 		 >> $2
		echo 'keyword_change2= plot_legend1,Gal400-dens'	 		 >> $2

	else

		#OII
		#echo 'x_range_min= 9'			 			 >> $2
		#echo 'x_range_max= 13'			 		 	 >> $2
		#echo 'y_range_min= -7.5'			 		 	 >> $2
		#echo 'y_range_max= -1'			 		 >> $2
		echo 'float_format_y= 0f'         				>> $2
		echo 'adjust_right= 0.99'					>> $2
		echo 'adjust_top= 0.99'					>> $2
		echo 'adjust_left= 0.16'					>> $2
		echo 'adjust_bottom= 0.18'					>> $2
		echo 'size_x= 12'				 		 >> $2
		echo 'size_y= 8'				 		 >> $2
		#LG				 	 
		#echo 'x_range_min= 10.35'			 			 >> $2
		#echo 'x_range_max= 12.05'			 		 	 >> $2
		#echo 'y_range_min= -6'			 		 	 >> $2
		#echo 'y_range_max= -1.8'			 		 >> $2

	endif

	echo 'adjust_right= 0.99'					>> $2
	echo 'adjust_top= 0.99'					>> $2
	echo 'adjust_left= 0.15'					>> $2
	echo 'adjust_bottom= 0.12'					>> $2
	echo 'adjust_bottom= 0.17'					>> $2

	echo 'minor_ticks_x_space= 0.1'			 	>> $2
	echo 'minor_ticks_y_space= 0.1'			 	>> $2

	echo 'mylw= 8'				 		 	 >> $2
	echo 'error_bars_y= yes'					 >> $2
	echo 'error_bars_y_sub= yes'					 >> $2
	#echo 'no_yticks= False'					>> $2
	echo 'xticks_minor_offset= 0'				>> $2
	
	#echo 'x_title= '	         >> $2
	#echo 'no_xticks= True'         			>> $2
	#echo 'no_first_xticks= True'         			>> $2
	echo 'no_first_last_xticks= True'         			>> $2

	#echo 'y_title= $\log_{10}$ ($n_{gal}$ (> X) [$Mpc^{-3}$])'  >> $2
	#echo 'y_title= $\log_{10}$ ($\Phi$/Mpc$^{-3}$ dlog$_{10}$M$_{\star}$)  '       >> $2
	#echo 'y_title= $\log_{10}$ ($\Phi$ [$Mpc^{-3}$ $dex^{-1}$])  '       >> $2
	echo 'no_last_yticks= True'         			>> $2
	#echo 'y_title= '  >> $2
	#echo 'no_yticks= True'         			>> $2

	
	echo 'z_print_position_x= 0.8'		>> $2
	echo 'z_print_position_y= 0.9'		>> $2

	echo 'use_loc_co_x= 0.2'		>> $2
	echo 'use_loc_co_y= 0.2'		>> $2

								
else if ($1 == dSFRF) then
	echo 'print_redshift= z=0.14' >> $2
	echo 'title= SFR Function MDPL '$3' '$4				 >> $2
	echo 'x_title= $\log_{10}$ (SFR [$M_{\odot} yr^{-1}$])' 	 >> $2
	#echo 'x_title= ' 	 >> $2
	echo 'no_first_yticks= False'         			>> $2
	echo 'no_xticks= False'         			>> $2
	echo 'y_title= $\log_{10}$ ($\Phi$/Mpc$^{-3}$ dlog$_{10}$SFR)  '       >> $2
	echo 'no_yticks= False'					>> $2
	#echo 'markersize= 4'			 			 >> $2
	echo 'x_range_min= -1'			 			 >> $2
	echo 'x_range_max= 3.5'			 		 	 >> $2
	echo 'y_range_min= -6'			 			 >> $2
	echo 'y_range_max= -1'			 			 >> $2
	echo 'use_loc_co= yes'			>> $2

	echo 'float_format_y= 0f'         				>> $2
	echo 'adjust_right= 0.96'					>> $2
	echo 'adjust_top= 0.97'					>> $2
	echo 'adjust_left= 0.12'					>> $2
	echo 'adjust_bottom= 0.16'					>> $2
	echo 'size_x= 14'				 		 >> $2
	echo 'size_y= 7.5'				 		 >> $2

	echo 'minor_ticks_x_space= 0.1'			 	>> $2
	echo 'minor_ticks_y_space= 0.1'			 	>> $2

	echo 'mylw= 6'				 		 	 >> $2
	echo 'error_bars_y= no'					 >> $2
	echo 'error_bars_y_sub= yes'					 >> $2

	echo 'xticks_minor_offset= 0'				>> $2

	echo 'z_print_position_x= 0.20'		>> $2
	echo 'z_print_position_y= 0.4'		>> $2
	echo 'use_loc_co_x= 0.14'		>> $2
	echo 'use_loc_co_y= 0.2'		>> $2



	echo 'keyword_change1= plot_legend0,Gal'	 	>> $2
	echo 'keyword_change2= plot_legend1,Gal2'	 >> $2
	echo 'keyword_change3= plot_legend2,SAGv2'	 >> $2
	echo 'keyword_change4= plot_legend3,SAGE'	 >> $2

else if ($1 == dsSFRF) then
	echo 'title= specific SFR MDPL' $3' '$4	 		 >> $2
	echo 'x_title= sSFR [$yr^{-1}$]'			 >> $2
	echo 'y_title= $\Phi$  [$Mpc^{-3}$ $dex^{-1}$] '        >> $2 

	echo 'x_range_min= -17'			 			 >> $2
	echo 'x_range_max= -7'			 		 	 >> $2
	echo 'y_range_min= -9'			 			 >> $2
	echo 'y_range_max= 0.0'			 			 >> $2
			
endif

if ($1 == cSFRD) then
	echo 'x_title= z'				 			 >> $2 
	echo 'y_title=  $\log_{10}$ (cSFR density [$M_{\odot}$ $yr^{-1}$ $Mpc^{-3}$])  '  >> $2
	echo 'cut_part1= '					 >> $2
	echo 'cut_part2= '					 >> $2
	#if x-axis is log x_ranges must be min/max -- trust me!!!

	#echo 'x_range_max= max'			 		 	 >> $2 
	#OII-sfr-paper
	#SFH cSFRD 	
	echo 'y_range_min= -2.3'			 		 	 >> $2 
	echo 'y_range_max= -0.5'			 			 >> $2
	#lower left
	echo 'use_loc_co_x= 0.50'			>> $2
	echo 'use_loc_co_y= 0.14'			>> $2
	#echo 'y_range_min= min'			 		 	 >> $2 
	#echo 'y_range_max= max'			 			 >> $2 


	echo 'minor_ticks_x_space= no'			 	>> $2
	echo 'xticks_minor= 0'					>> $2
	echo 'minor_ticks_x_space= no'			 	>> $2
	echo 'xticks= 0,0.1,0.25,0.5,0.75,1,1.25,1.5,2'		>> $2
	#echo 'xticks= 0.5,0.75,1,1.5,2,3,4,6'		>> $2
	echo 'xticks_minor= 0'					>> $2
	echo 'xticks_minor_offset= 1'				>> $2

	echo 'x_range_min= min'			 		 	 >> $2 
	echo 'x_range_max= max'			 			 >> $2 

	echo 'log_scale_x= yes'						 >> $2
	echo 'log_scale_x_sub= yes'					 >> $2

	echo 'error_bars_x= no'					 >> $2
	echo 'error_bars_x_sub= yes'					 >> $2
	echo 'error_bars_y= no'					 >> $2
	echo 'error_bars_y_sub= yes'					 >> $2

	echo 'no_last_yticks= True'			 	>> $2

	echo 'keyword_change1= plot_legend0,Galacticus'	 		 >> $2
	echo 'keyword_change2= plot_legend1,SAG'	 		 >> $2
	echo 'keyword_change3= plot_legend2,SAGE'			 >> $2
	echo 'keyword_change4= plot_legend3,SAGE'			 >> $2
	 
endif


if ($1 == sfr2mstar) then
	echo 'title= sSFR/Mstar MDPL '$3' '$4	 >> $2
	echo 'x_title= $\log_{10}$ ($M_{_*}$ [$M_{\odot}$])'	 >> $2
	echo 'y_title= $\log_{10}$ (SFR [$Gyr^{-1}$])  '         >> $2 

	echo 'x_range_min= 8.5'			 			 >> $2
	echo 'x_range_max= 12'			 		 	 >> $2
	echo 'y_range_min= 6'			 		 	 >> $2
	echo 'y_range_max= 12'			 		 >> $2
	echo 'legend_position= upper right'			 	 >> $2
	echo 'linestyle= -'			 			 >> $2
	echo 'linestyle_change=  '		 			 >> $2
	echo 'linestyle_sub= '			 			 >> $2
	echo 'linestyle_sub_change= '		 			 >> $2
	#echo 'keyword_change1= plot_legend0,'$6' cal01 125 Mpc/h z='$8		 >> $2
	#echo 'keyword_change2= plot_legend1,'$7' 400 Mpc/h z='$9	 >> $2
	#echo 'keyword_change3= default'		 			 >> $2
	echo 'log_scale_x= no'						 >> $2
	echo 'log_scale_y= no'						 >> $2
	echo 'log_scale_x_sub= no'					 >> $2
	echo 'log_scale_y_sub= no'					 >> $2
	echo 'error_bars_x= no'					 >> $2
	echo 'error_bars_y= yes'					 >> $2
	echo 'error_bars_x_sub= no'					 >> $2
	echo 'error_bars_y_sub= yes'					 >> $2

endif

if ($1 == oh2mstar || $1 == zgas2mstar || $plotXY_key == zcold) then
	echo 'print_redshift= MDPL2 SAGE' >> $2
	echo 'print_redshift2= sSFR cut: Henriques+20 A1' >> $2
	echo 'title= (O/H) vs. Mstar MDPL '$3' '$4	 >> $2
	echo 'x_range_min= 9'			 			 >> $2
	echo 'x_range_max= 11.5'			 		 	 >> $2
	echo 'y_range_min= 8'			 		 	 >> $2
	echo 'y_range_max= 9.6'			 		 >> $2
	
	echo 'add_axis= yes'						 >> $2
	echo 'add_which_axis= x'					>> $2
	echo 'add_which_col= 8'						>> $2
	echo 'log_scale_add= no'						>> $2
	echo 'y_title= 12+ $\log_{10}$ (0/H)'		 	 >> $2
	echo 'x_title= $\log_{10}$ ($M_{_*}$ [$M_{\odot}$])'	 	 >> $2
	echo 'title_add= $\log_{10}$ ($M_{Cold}/M_{_*}$)'		 	 >> $2


	echo 'use_loc_co= yes'			>> $2

	echo 'minor_ticks_x_space= 0.1'			 	>> $2
	echo 'minor_ticks_y_space= 0.1'			 	>> $2

	#echo 'xticks= 8.5,9.5,10.5,11.5'		>> $2

	echo 'error_bars_y= yes'					 >> $2
	echo 'error_bars_y_sub= yes'					 >> $2
	echo 'no_yticks= False'					>> $2

	echo 'z_print_position_x= 0.18'		>> $2
	echo 'z_print_position_y= 0.92'		>> $2

	echo 'z_print_position_x2= 0.18'		>> $2
	echo 'z_print_position_y2= 0.18'		>> $2

	echo 'use_loc_co_x= 0.75'		>> $2
	echo 'use_loc_co_y= 0.16'		>> $2

	echo 'use_loc_co_x= 0.55'		>> $2
	echo 'use_loc_co_y= 0.72'		>> $2
	echo 'no_last_yticks= True'         			>> $2
	echo 'no_last_xticks= True'         			>> $2

	echo 'keyword_change1= plot_legend0,Galacticus'	 	>> $2
	#echo 'keyword_change2= plot_legend1,Galacticus 2'	 >> $2
	echo 'keyword_change2= plot_legend1,SAG'	 >> $2
	echo 'keyword_change3= plot_legend2,SAGE'	 >> $2

	echo 'keyword_change1= plot_legend0,z=0.0'	 	>> $2
	echo 'keyword_change2= plot_legend1,z=0.5'	 >> $2
	echo 'keyword_change3= plot_legend2,z=1.0'	 >> $2
	echo 'keyword_change4= plot_legend3,z=2.0'	 >> $2
	echo 'keyword_change5= plot_legend4,z=3.0'	 >> $2
	echo 'keyword_change6= plot_legend5,z=4.0'	 >> $2
	echo 'keyword_change7= plot_legend6,z=5.0'	 >> $2

	set add_x = 'yess'
	if ($add_x == 'yes') then
			echo 'add_axis= yes'				>> $2
			echo 'add_which_col_x= 7'				>> $2
			echo 'log_scale_add= no'						>> $2
			echo 'title_add_x0= $\log_{10}$ ($M_{Cold}/M_{_*}$) Galacticus'		 	 >> $2
			#echo 'title_add_x0= '		 	 >> $2
			echo 'title_add_x1= $\log_{10}$ ($M_{Cold}/M_{_*}$) SAG'		 	 >> $2
			echo 'title_add_x2= $\log_{10}$ ($M_{Cold}/M_{_*}$) SAGE'		 	 >> $2		

			echo 'size_y= 12'				 		 >> $2
			echo 'z_print_position_y= 0.65'		>> $2
			echo 'minor_ticks_x_set_top= off'				>> $2

			echo 'adjust_left= 0.14'					>> $2
			echo 'adjust_right= 0.99'					>> $2
			echo 'adjust_bottom= 0.12'					>> $2
			echo 'adjust_top= 0.72'					>> $2

			echo 'ylabel_pos= 0.4'	         			>> $2 

		endif

endif

echo 'ssfr2mstar_ref= False'					 >> $2
if ($1 == ssfr2mstar) then
	echo 'title= sSFR/Mstar MDPL '$3' '$4	 >> $2
	echo 'x_title= $\log_{10}$ ($M_{_*}$ [$M_{\odot}$])'	 >> $2
	echo 'y_title= $\log_{10}$ (sSFR [$yr^{-1}$])  '         >> $2 
	echo 'x_range_min= 8.5'			 			 >> $2
	echo 'x_range_max= 12.0'			 		 	 >> $2
	echo 'y_range_min= -11.25'			 		 	 >> $2
	echo 'y_range_max= -9'			 		 >> $2
	echo 'legend_position= upper right'			 	 >> $2
	echo 'ssfr2mstar_ref= Elbaz11'					 >> $2

	#star forming Paper I 
	echo 'hline_ypos= -11'					 >> $2

	#star forming OII-sfr paper 2.51e-11 (0.3/t_hubble(z=0.0)) --> for Paper I
	echo 'hline_ypos= -10.67'					 >> $2	

	#star forming OII-sfr paper 3.547e-11 (0.3/t_hubble(z=0.1)) --> for SFR-OII Paper
	#echo 'hline_ypos= -10.45'					 >> $2	
endif

if ($1 == HMFc || $1 == HMF_noc ) then
	echo 'title= Halo Mass Function (incl. orphans) '$3' '$4	 		 >> $2
	echo 'x_title= $M_{AB_i}$'	         >> $2
	echo 'y_title= $\log_{10}$ ($\Phi$ [$Mpc^{-3}$] $mag^{-1}$)  '       >> $2

	echo 'x_range_min= -27'			 		 >> $2
	echo 'x_range_max= -15'			 		 >> $2
	echo 'y_range_min= -7'			 			 >> $2
	echo 'y_range_max=  0'			 			 >> $2
	echo 'keyword_change1= plot_legend0,Galacticus region0001'	 		 >> $2
	echo 'keyword_change2= plot_legend1,SAG'	 		 >> $2
	echo 'keyword_change3= plot_legend2,SAGE'			 >> $2
		
endif

if ($1 == HMF|| $1 == HMF_no ) then
	echo 'title= Halo Mass Function (incl. orphans) '$3' '$4	 		 >> $2
	echo 'y_title= $\log_{10}$ ($\Phi$ [$h^3Mpc^{-3}dex^{-1}$])  '       >> $2
	echo 'x_title= $\log_{10}$ ($M_{Halo}$ [$h^{-1}M_{\odot}$])'       >> $2

	echo 'keyword_change1= plot_legend0,Galacticus'	 		 >> $2
	echo 'keyword_change2= plot_legend1,SAG'	 		 >> $2
	echo 'keyword_change3= plot_legend2,SAGE'			 >> $2
	#echo 'y_title= $\log_{10}$ ($n_{gal}$ (> $M_{Halo}$ [$h^{-1}M_{\odot}$]) [$h^{3}Mpc^{-3}$])'  >> $2
	echo 'x_range_min= 9'			 		 >> $2
	echo 'x_range_max= 15'			 		 >> $2
	echo 'y_range_min= -2'			 			 >> $2
	echo 'y_range_max= 2'			 			 >> $2

	echo 'print_redshift= region0021, r=1 [h-1Mpc]' 		>> $2
		
endif


if ($1 == mstar2mhalo || $1 == mstar2mhalo_no || $1 == mstar2mhalovsSFR || $plotXY_key == 'HODsc') then
	echo 'title= Mstar/Mhalo  MDPL '$3' '$4	 			 >> $2

	if ($1 == mstar2mhalovsSFR) then 
		echo 'x_title= $\log_{10}$ (sSFR $[yr^{-1}]$)'	 	>> $2
		echo 'x_range_min= -15'			 		 	>> $2
		echo 'x_range_max= -7'			 		 	>> $2
	else
		echo 'x_title= $\log_{10}$ ($M_{200c}$ $[M_{\odot}]$)'	 	>> $2
		echo 'x_range_min= 10.0'			 		 	>> $2
		echo 'x_range_max= 15.4'			 		 	>> $2
	endif
	echo 'y_title= $\log_{10}$ ($M_{\star}/M_{Halo}$)'	 	 >> $2
	#echo 'x_title=  '	 	 	 	>> $2
	#echo 'no_xticks= True'         				>> $2
	#echo 'y_title=  '	 	 	 	>> $2
	#echo 'no_yticks= True'         				>> $2

	set use_cb_colors = 3
	echo 'use_cb_colors= 3'	>> $2

	echo 'y_range_min= -3.0'			 		 >> $2
	echo 'y_range_max= -0.5'		 			 	 >> $2

	set add_x = 'yes'
	set add_y = 'no'

	echo 'minor_ticks_x_space= 0.1'			 	>> $2
	echo 'minor_ticks_x_set_top= off'				>> $2
	echo 'minor_ticks_y_space= 0.1'			 	>> $2
	echo 'minor_ticks_y_set_right= on'				>> $2
	echo 'xticks= default'			 			>> $2
	echo 'yticks= default'			 			>> $2
	echo 'plot_legend= yes'		 			>> $2
	echo 'legend_fontsize= 24'				 	 >> $2
	echo 'text_fontsize= 24'				 	 >> $2

	echo 'use_loc_co_x= 0.70'			>> $2
	echo 'use_loc_co_y= 0.45'			>> $2
	echo 'print_redshift= z=0.0'				 	 >> $2
	#echo 'print_redshift2= $M_{Halo}$: before infall' >> $2
	#echo 'print_redshift2= $M_{Halo}$: bound mass' >> $2

	echo 'z_print_position_x= 0.50'		>> $2
	echo 'z_print_position_y= 0.65'		>> $2
	echo 'z_print_position_x2= 0.50'		>> $2
	echo 'z_print_position_y2= 0.60'		>> $2
	#echo '3= #c5c1c2' >> $16
	echo 'size_y= 8'				 		 >> $2
	echo 'adjust_bottom= 0.11'					>> $2
	echo 'error_bars_y= yes'					 >> $2

	if ($add_x == 'yes') then
		echo 'add_axis= yes'				>> $2
		echo 'add_which_col_x= 7'				>> $2
		echo 'log_scale_add= no'						>> $2
		echo 'title_add_x0= $\log_{10}$ ($M_{\star}$ $[M_{\odot}]$) Gal-MD1000'		 	 >> $2
		#echo 'title_add_x0= '		 	 >> $2
		echo 'title_add_x1= $\log_{10}$ ($M_{\star}$ $[M_{\odot}]$) Gal-SMD400'		 	 >> $2
		echo 'title_add_x2= $\log_{10}$ ($M_{\star}$ $[M_{\odot}]$) Ill-TNG300'		 	 >> $2	
		echo 'adjust_top= 0.72'					>> $2
		echo 'size_y= 11'				 		 >> $2

		#echo 'minor_ticks_x_set_top= off'				>> $2
		#echo 'ylabel_pos= 0.4'	         			>> $2   
	endif

	if ($add_y == 'yes') then

		echo 'add_axis= yes'				>> $2	
		echo 'adjust_right= 0.72'					>> $2
		echo 'adjust_left= 0.11'					>> $2
		echo 'z_print_position_x= 0.55'		>> $2
		echo 'z_print_position_y= 0.90'		>> $2
		echo 'z_print_position_x2= 0.55'		>> $2
		echo 'z_print_position_y2= 0.85'		>> $2

		#echo 'minor_ticks_y_set_right= off'				>> $2

		echo 'ylabel_pos= 0.42'				>> $2
		echo 'xlabel_pos= 0.42'				>> $2
		echo 'add_which_col_y= 6'				>> $2
		echo 'xlabel_pad= 0.02'					>> $2

		echo 'title_add_y0= $\log_{10}$ (SFR $[M_{\odot}yr^{-1}]$) Gal-MD1000'		 	 >> $2
		#echo 'title_add_x0= '		 	 >> $2
		echo 'title_add_y1= $\log_{10}$ (SFR $[M_{\odot}yr^{-1}]$) Gal-SMD400'		 	 >> $2
		echo 'title_add_y2= $\log_{10}$ (SFR $[M_{\odot}yr^{-1}]$) Ill-TNG300'		 	 >> $2	
		echo 'size_x= 16' >> $2
		#echo 'title_add_y0= '		 	 >> $2

		echo 'use_loc_co_x= 0.75'			>> $2
		echo 'use_loc_co_y= 0.73'			>> $2

	endif

	set key = 'mstarr'
	if ($key == 'mstar') then
		echo 'y_title= $\log_{10}$ ($M_{_*}$ $[M_{\odot}]$)'		 	 >> $2
		echo 'title_add0= $\log_{10}$ ($M_{_*}/M_{200c}$) '$prefix'-cols'	 	 >> $2
		#echo 'x_range_max= 15'			 		 	>> $2
		#echo 'x_title=  '	 	 	 	>> $2
		#echo 'no_xticks= True'         				>> $2
		#echo 'y_title=  '	 	 	 	>> $2
		#echo 'no_yticks= True'         				>> $2

		echo 'use_loc_co_x= 0.72'			>> $2
		echo 'use_loc_co_y= 0.25'			>> $2

		echo 'y_range_min= 8.5'			 		 >> $2
		echo 'y_range_max= 13'		 			 >> $2
		echo 'minor_ticks_x_space= 0.1'			 	>> $2
		echo 'minor_ticks_y_space= 0.25'			 	>> $2
		echo 'minor_ticks_x_set_top= on'				>> $2
		echo 'add_axis= no'						 >> $2

		echo 'size_y= 6'				 		 >> $2
		echo 'z_print_position_x= 0.19'		>> $2
		echo 'z_print_position_y= 0.90'		>> $2
		echo 'z_print_position_x2= 0.19'		>> $2
		echo 'z_print_position_y2= 0.85'		>> $2
		#echo 'print_redshift2= ' >> $2
		echo 'plot_legend= no'		 			>> $2		

		echo 'adjust_top= 0.99'					>> $2
		echo 'adjust_bottom= 0.18'					>> $2
		#echo 'adjust_left= 0.15'					>> $2
		#echo 'minor_ticks_y_set_right= on'				>> $2
		#echo 'adjust_right= 0.99'					>> $2
		echo 'print_redshift= '				 	 >> $2

		echo 'error_bars_y= yes'					 >> $2
		echo 'hline_ypos= 9.5,10.5,11.5'					 >> $2	
		echo 'vline_xpos= 11,12,13.5'					 >> $2
		echo 'legend_fancy= True'					>> $2
	
	else if ($plotXY_key == 'HODsc' && $1 !~ mstar2mhal* ) then

		echo 'y_title= $\sigma_{\log_{10}M_{_*}}$ [dex]'		 	 >> $2
		#echo 'y_title=  '	 	 	 	>> $2
		#echo 'no_yticks= True'         				>> $2
		#echo 'no_first_yticks= True'         				>> $2

		echo 'use_loc_co_x= 0.6'			>> $2
		echo 'use_loc_co_y= 0.3'			>> $2
		#echo 'x_range_max= 15.0'			 		 	>> $2
		echo 'y_range_min= 0.0'			 		 >> $2
		echo 'y_range_max= 0.5'		 			 >> $2

		echo 'print_redshift2= ' >> $2

		echo 'plot_legend= no'		 			>> $2

		echo 'add_axis= no'						 >> $2
		echo 'no_first_last_yticks= True'         				>> $2
		echo 'float_format_y= 1f'         				>> $2
		echo 'adjust_top= 0.99'					>> $2
		echo 'adjust_bottom= 0.28'					>> $2

		echo 'xlabel_pad= 0.03'					>> $2
		echo 'ylabel_pad= 0.005'				>> $2
		#echo 'xlabel_pos= 0.45'				>> $2
		echo 'ylabel_pos= 0.6'				>> $2

		echo 'major_ticks_x_space= 0.25'			 	>> $2
		echo 'minor_ticks_x_space= 0.1'			 	>> $2
		echo 'minor_ticks_y_space= 0.05'			 	>> $2
		echo 'minor_ticks_y_set_right= on'				>> $2
		echo 'minor_ticks_x_set_top= on'				>> $2

		echo 'size_y= 4'				 		 >> $2
		#echo 'size_4= 16'				 		 >> $2
		echo 'z_print_position_x= 0.55'		>> $2
		echo 'z_print_position_y= 0.38'		>> $2
		#echo 'adjust_right= 0.72'					>> $2
		echo 'print_redshift= False'				 	 >> $2

		echo 'error_bars_y= no'					 >> $2
	endif

	echo 'keyword_change1= plot_legend0,'$prefix'-cols'	 		 >> $2
	echo 'keyword_change2= plot_legend1,'$prefix'-dens $n_{CMASS}$'			 >> $2
	echo 'keyword_change3= plot_legend2,'$prefix'-dens g-i>2.35'			 >> $2
	echo 'keyword_change4= plot_legend3,'$prefix'-dens Guo+13'			 >> $2
	echo 'keyword_change5= plot_legend4,'$prefix'-mass'			 >> $2

	if ($catname =~ 'Ga'*) then
		echo 'keyword_change1= plot_legend0,'$catname''	 		 >> $2
		echo 'keyword_change2= plot_legend1,'$catname'2'			 >> $2
		echo 'keyword_change3= plot_legend2,'$catname'-mass'			 >> $2
	else

		echo 'keyword_change1= plot_legend0,'$catname'-dens'	 		 >> $2
		echo 'keyword_change2= plot_legend1,'$catname'-mass'			 >> $2
		echo 'keyword_change3= plot_legend2,'$catname'-mass'			 >> $2
	endif

	if ($key != 'mstar') then
		echo 'z_print_position_x= 0.73'		>> $2
		echo 'z_print_position_y= 0.65'		>> $2
		echo 'z_print_position_x2= 0.73'		>> $2
		echo 'z_print_position_y2= 0.61'		>> $2
		echo 'z_print_position_x3= 0.73'		>> $2
		echo 'z_print_position_y3= 0.57'		>> $2
	endif

	echo 'keyword_change1= plot_legend0,Gal-MD1000'	 		 >> $2
	echo 'keyword_change2= plot_legend1,Gal-SMD400'	 		 >> $2
	echo 'keyword_change3= plot_legend2,Ill-TNG300'			 >> $2
	
	#echo 'print_redshift= z=0.0, centrals'				 	 >> $2
	echo 'print_redshift2= '				 	 >> $2
	#echo 'print_redshift2= $\log_{10}$ ($M_{200c}$ $[M_{\odot}]$) $\leq$11'				 	 >> $2
	#echo 'print_redshift2= 11$\leq$ $\log_{10}$ ($M_{Halo}$ $[M_{\odot}]$) $\leq$12'				 	 >> $2
	#echo 'print_redshift2= 12$\leq$ $\log_{10}$ ($M_{Halo}$ $[M_{\odot}]$) $\leq$13.5'				 	 >> $2
	#echo 'print_redshift2= $\log_{10}$ ($M_{200c}$ $[M_{\odot}]$) $\geq$13.5'				 	 >> $2
	echo 'print_redshift3= '				 	 >> $2	
	#echo 'print_redshift3= $\log_{10}$ ($M_{\star}$ $[M_{\odot}]$) $\leq$9.5'				 	 >> $2
	#echo 'print_redshift3= 9.5$\leq$ $\log_{10}$ ($M_{\star}$ $[M_{\odot}]$) $\leq$10.5'				 	 >> $2
	#echo 'print_redshift3= 10.5$\leq$ $\log_{10}$ ($M_{\star}$ $[M_{\odot}]$) $\leq$11.5'				 	 >> $2
	#echo 'print_redshift3= $\log_{10}$ ($M_{\star}$ $[M_{\odot}]$) $\geq$11.5'				 	 >> $2
	echo 'vline_xpos= 11,12,13.5'					 >> $2
	echo 'legend_fancy= True'					>> $2
	#echo 'x_title=  '	 	 	 	>> $2
	#echo 'no_xticks= True'         				>> $2	

	#echo 'y_title=  '	 	 	 	>> $2
	echo 'no_last_yticks= True'         				>> $2

endif

if ($1 == zgas2mstarr) then
	echo 'title= $Z_{Gas}$ to $M_{_*}$  '$3' '$4	 	>> $2
	echo 'x_title= $M_{_*}$ [$M_{\odot}$]'			>> $2
	echo 'y_title= $Z_{Cold}$'	>> $2 

	echo 'x_range_min= 8.5'				 	 	 >> $2
	echo 'x_range_max=  12.5'				 	 >> $2
	echo 'y_range_min= 7'			 			 >> $2
	echo 'y_range_max=  10'			 			 >> $2
endif

if ($1 == mcold2mstar) then
	echo 'title= Cold Gas Mass $M_{Cold}$ to Mstar $M_{_*}$ '$3' '$4 >> $2
	echo 'x_title= $\log_{10}$ ($M_{_*}$ [$M_{\odot}$])'			         >> $2
	echo 'y_title= $\log_{10}$ ($M_{Cold}/M_{_*}$)'		         >> $2 
	echo 'legend_position= upper right'			 	 >> $2
	echo 'use_loc_co_x= 0.56'			>> $2
	echo 'use_loc_co_y= 0.72'			>> $2
	echo 'x_range_min= 8.5'				 	 	 >> $2
	echo 'x_range_max= 12.5'				 	 >> $2
	echo 'y_range_min= -4'			 			 >> $2
	echo 'y_range_max= 2'			 			 >> $2
	echo 'print_redshift= 0.0'				 	 >> $2
		
endif

if ($1 == mbh2mstarsph) then
	echo 'title= Black hole mass $M_{bh}$ to bulge mass '$3' '$4 >> $2
	echo 'x_title= $\log_{10}$ ($M_{Bulge}$ [$M_{\odot}$])'	 >> $2
	echo 'y_title= $\log_{10}$ ($M_{BH}$ [$M_{\odot}$])'		 >> $2 
	echo 'legend_position= lower right'			 	 >> $2
	echo 'use_loc_co_x= 0.60'			>> $2
	echo 'use_loc_co_y= 0.05'			>> $2
	echo 'x_range_min= 9'				 	 	 >> $2
	echo 'x_range_max= 13'				 	 >> $2
	echo 'y_range_min= 5'			 			 >> $2
	echo 'y_range_max= 11'			 			 >> $2

	echo 'error_bars_x= no'					 	 >> $2
	echo 'error_bars_y= yes'					 >> $2
	echo 'error_bars_x_sub= no'					 >> $2
	echo 'error_bars_y_sub= no'					 >> $2
	echo 'print_redshift= 0.0'				 	 >> $2

endif

#if ($1 == SFRF || $1 == sSFRF || $1 =~ *22* ) then
#	echo 'keyword_change1= plot_legend0,Galacticus' 	>> $2
#	echo 'keyword_change2= plot_legend1,SAG'		>> $2
#	echo 'keyword_change3= plot_legend2,SAGE'		>> $2
#endif

#OII-sfr paper
#echo 'keyword_change1= plot_legend0,Galacticus' 	>> $2
#echo 'keyword_change2= plot_legend1,SAG v3'		>> $2
#echo 'keyword_change3= plot_legend2,SAG v2'		>> $2
#echo 'keyword_change4= plot_legend3,SAGE'	>> $2

#echo 'keyword_change1= plot_legend0,Gal-all'	 		 >> $2
#echo 'keyword_change1= plot_legend0,Gal-cols'	 		 >> $2
#echo 'keyword_change2= plot_legend1,Gal-dens'			 >> $2
#echo 'keyword_change3= plot_legend2,Gal-mass'			 >> $2

#echo 'keyword_change1= plot_legend0,Guo+13' 		 >> $2
#echo 'keyword_change2= plot_legend1,$n_{CMASS}$'	 		 >> $2
#echo 'keyword_change3= plot_legend2,g-i>2.35'			 >> $2
#echo 'keyword_change4= plot_legend3,Gal-all'			 >> $2

#echo 'keyword_change1= plot_legend0,z=0.0' 	>> $2

echo 'caption1= ' >> $2
echo 'caption2= ' >> $2
echo 'caption3= ' >> $2

if ($catname == 'SAGg') then

	echo '0= --'	>> $13
	echo '1= :'	>> $13
	echo '2= '	>> $13

	echo '2= o'	>> $15
	echo '3= ^'	>> $15
	echo '4= s'	>> $15

	echo '2= w'	>> $16
	echo '3= #ffe200' >> $16
	echo '4= #c5c1c2' >> $16

endif	

#DEFAULT colour style for contour plots
if ($1 =~ analyse* || $1 =~ plotX*) then

	#linestyle table: 	MATPLOT_LINESTYLE
	echo '0= -'	>> $13

	#linestyle colour table:MATPLOT_COL
	echo '0= w' >> $14 
	echo '1= k' >> $14  
	echo '3= r' >> $14
	echo '2= #ffab00' >> $14
	#Butterfly green
	echo '5= #4baf0c' >> $14
	echo '4= #AA4499' >> $14
	echo '6= #4da6ff' >> $14
	echo '7= k' >> $14
	echo '8= k' >> $14
	echo '9= k' >> $14
	echo '10= k' >> $14
	echo '5= k' >> $14
	#echo '3= #f3bb1d' >> $14
	#echo '3= k' >> $14
	#echo '6= k' >> $14

	#greysafecols = ['#809BC8', '#FF6666', '#FFCC66', '#64C204']f

	#Choose nr of colours-code which should be color-blind (cb) friendly created: choose 0, if you do not want any cb friendly color
	#echo 'use_cb_colors= '$use_cb_colors	>> $2
	echo 'use_cb_colors= 0'	>> $2

	if ($1 =~ analyse*) then
		echo 'legend_fontsize= 18'				 	 >> $2
		echo 'size_x= 12'				 		 >> $2
		echo 'size_y= 10'				 		 >> $2

		if ($1 =~ analyse*_sfr_mstarr*) then 
			echo '1= --'	>> $13
			echo '2= '	>> $13
			echo '3= '	>> $13

			#markerstyle table: 	MATPLOT_MARKERSTYLE
			echo '0= '	>> $15
			echo '1= '	>> $15
			echo '2= s'	>> $15
			echo '3= o'	>> $15

			#marker colour table: 	MATPLOT_MARKERCOL
			echo '0= w'	>> $16
			echo '1= w'	>> $16
			echo '2= #ffe200'	>> $16
			echo '3= k'	>> $16

		else if ($1 =~ analyse*dmesa* || $1 =~ analyse*r-i* || $1 =~ analyse*rdisk* || $1 =~ analyse*rhalf* || $1 =~ analyse*mstar* || $1 =~ analyse*mbar* || $1 =~ analyse*bdisk* || $1 =~ analyse*u-*r* || $1 =~ analyse*zgas*  || $1 =~ analyse*sfr* || $1 =~ analyse*mhalo*  ) then
			echo '1= :'	>> $13
			echo '2= -'	>> $13
			echo '3= --'	>> $13
			echo '4= -'	>> $13
			echo '5= -'	>> $13
			echo '6= :'	>> $13
			echo '7= '	>> $13
			echo '8= --'	>> $13
			echo '9= --'	>> $13

			#markerstyle table: 	MATPLOT_MARKERSTYLE
			echo '0= '	>> $15
			echo '1= '	>> $15
			echo '2= '	>> $15
			echo '3= '	>> $15
			echo '4= '	>> $15
			echo '5= '	>> $15
			echo '6= '	>> $15
			echo '7= o'	>> $15
			echo '8= '	>> $15
			echo '9= '	>> $15

			#marker colour table: 	MATPLOT_MARKERCOL
			echo '0= w'	>> $16
			echo '1= w'	>> $16
			echo '2= #ffab00'	>> $16
			echo '3= w'	>> $16
			echo '4= #AA4499'	>> $16
			echo '5= None'	>> $16
			echo '6= None'	>> $16
			#echo '6= k'	>> $16
			echo '7= #ffe200'	>> $16
			echo '8= #ffe200'	>> $16
			echo '9= #ffe200'	>> $16

		else if ($1 =~ analyse*) then

			#linestyle colour table:MATPLOT_COL
			echo '0= k' >> $14
			echo '1= #225588' >> $14  
			echo '2= k' >> $14  
			echo '3= w' >> $14
			echo '4= #4baf0c' >> $14
			echo '5= m' >> $14
			echo '6= k' >> $14
			echo '6= k' >> $14
			echo '7= k' >> $14
			echo '8= m' >> $14


			echo '1= -'	>> $13
			echo '2= :'	>> $13
			echo '3= --'	>> $13
			echo '4= '	>> $13
			echo '5= --'	>> $13
			echo '6= -'	>> $13
			echo '7= -'	>> $13
			echo '8= '	>> $13
			echo '9= -'	>> $13

			#markerstyle table: 	MATPLOT_MARKERSTYLE
			echo '0= '	>> $15
			echo '1= '	>> $15
			echo '2= '	>> $15
			echo '3= '	>> $15
			echo '4= o'	>> $15
			echo '5= '	>> $15
			echo '6= '	>> $15
			echo '7= '	>> $15
			echo '8= o'	>> $15
			echo '9= '	>> $15

			#marker colour table: 	MATPLOT_MARKERCOL
			echo '0= w'	>> $16
			echo '1= w'	>> $16
			echo '6= w'	>> $16
			echo '2= w'	>> $16
			echo '3= #ffe200'	>> $16
			echo '5= m'	>> $16
			echo '4= #ffe200'	>> $16
			echo '8= #ffe200'	>> $16
			echo '7= w'	>> $16
			echo '9= w'	>> $16

		else
			echo '1= -'	>> $13
			echo '2= --'	>> $13
			echo '3= '	>> $13
			echo '4= '	>> $13
			echo '5= -'	>> $13
			echo '6= --'	>> $13

			#markerstyle table: 	MATPLOT_MARKERSTYLE
			echo '0= '	>> $15
			echo '1= '	>> $15
			echo '2= '	>> $15
			echo '3= '	>> $15
			echo '4= o'	>> $15
			echo '5= '	>> $15
			echo '6= '	>> $15

			#marker colour table: 	MATPLOT_MARKERCOL
			echo '0= w'	>> $16
			echo '1= w'	>> $16
			echo '2= w'	>> $16
			echo '3= w'	>> $16
			echo '4= w'	>> $16
			echo '6= w'	>> $16
			echo '5= w'	>> $16
		endif

	else if ($24 =~ *isto* ) then

		#Choose nr of colours-code which should be color-blind (cb) friendly created: choose 0, if you do not want any cb friendly color
		echo 'use_cb_colors= 0'					>> $2
		#linestyle table: 	MATPLOT_LINESTYLE
		echo '0= -'	>> $13
		echo '1= -'	>> $13
		echo '2= '	>> $13
		echo '3= --'	>> $13
		echo '4= '	>> $13
		echo '5= '	>> $13

		#linestyle colour table: MATPLOT_COL
		echo '0= k' >> $14 
		echo '1= #ffe200' >> $14 
		echo '2= k' >> $14
		echo '3= r' >> $14
		echo '4= w' >> $14
		echo '5= w' >> $14

		#marker colour table: 	MATPLOT_MARKERCOL
		echo '0= k'	>> $16
		echo '1= w'	>> $16
		echo '2= w'	>> $16
		echo '3= w'	>> $16
		echo '4= r'	>> $16
		echo '5= k'	>> $16

		#markerstyle table: 	MATPLOT_MARKERSTYLE --> here hatch
		echo '0= '	>> $15
		echo '1= /'	>> $15
		echo '2= .'	>> $15
		echo '3= '	>> $15
		echo '4= '	>> $15
		echo '5= '	>> $15

		#marker alpha: 	MATPLOT_ALPHA
		echo '0= 1'		>> $20
		echo '1= 1'		>> $20
		echo '2= 0.7'		>> $20
		echo '3= 1'		>> $20
		echo '4= 0.2'		>> $20
		echo '5= 0.1'		>> $20

	else if ($plotXY_key == 'zevol') then

		echo '0= #225588' >> $14 
		echo '1= #225588' >> $14  
		echo '2= #225588' >> $14

		echo '3= #CC6677' >> $14 
		echo '4= #CC6677' >> $14  
		echo '5= #CC6677' >> $14

		echo '6= #DDCC77' >> $14 
		echo '7= #DDCC77' >> $14  
		echo '8= #DDCC77' >> $14
		echo '9= k' >> $14 
		echo '10= k' >> $14  
		echo '11= k' >> $14

		echo '0= -'	>> $13						
		echo '1= --'	>> $13
		echo '2= :'	>> $13
		echo '3= -'	>> $13
		echo '4= --'	>> $13
		echo '5= :'	>> $13
		echo '6= -'	>> $13
		echo '7= --'	>> $13	
		echo '8= :'	>> $13
		echo '9= -'	>> $13
		echo '10= --'	>> $13
		echo '11= :'	>> $13
		echo '12= -'	>> $13
		echo '13= --'	>> $13
		echo '14= :'	>> $13

		#markerstyle table: 	MATPLOT_MARKERSTYLE
		echo '0= .'	>> $15
		echo '1= '	>> $15
		echo '2= '	>> $15
		echo '3= '	>> $15
		echo '4= '	>> $15
		echo '5= '	>> $15
		echo '6= '	>> $15
		echo '7= '	>> $15
		echo '8= '	>> $15
		echo '9= '	>> $15
		echo '10= '	>> $15
		echo '11= '	>> $15

	else if ($use_cb_colors == 2 || $use_cb_colors == 3 ) then

		echo '0= -'	>> $13						
		echo '1= -'	>> $13
		echo '2= c1'	>> $13
		echo '3= --'	>> $13
		echo '4= --'	>> $13
		echo '5= --'	>> $13

		echo '3= #225588' >> $14
		echo '4= #CC6677' >> $14  
		echo '5= #DDCC77' >> $14

		
		#markerstyle table: 	MATPLOT_MARKERSTYLE
		echo '0= '	>> $15
		echo '1= '	>> $15
		echo '2= '	>> $15
		echo '3= '	>> $15
		echo '4= '	>> $15
		echo '5= '	>> $15

		#marker colour table: 	MATPLOT_MARKERCOL
		echo '0= w'	>> $16
		echo '1= r'	>> $16
		echo '2= w'	>> $16
		echo '2= k'	>> $16
		echo '3= #225588'	>> $16
		echo '3= w'	>> $16
		echo '4= #CC6677'	>> $16
		echo '5= #DDCC77'	>> $16

	else if ($use_cb_colors == 4) then

		echo '0= -'	>> $13						
		echo '1= --'	>> $13
		echo '2= -'	>> $13
		echo '3= -.'	>> $13
		echo '4= '	>> $13
		echo '5= '	>> $13
		
		#markerstyle table: 	MATPLOT_MARKERSTYLE
		echo '0= '	>> $15
		echo '1= '	>> $15
		echo '2= '	>> $15
		echo '3= '	>> $15
		echo '4= <'	>> $15
		echo '5= h'	>> $15

		#marker colour table: 	MATPLOT_MARKERCOL
		echo '0= w'	>> $16
		echo '1= r'	>> $16
		echo '2= w'	>> $16
		echo '3= #225588'	>> $16
		echo '3= w'	>> $16
		echo '4= #CC6677'	>> $16
		echo '5= #DDCC77'	>> $16

	else if ($use_cb_colors == 5) then

		echo '0= -'	>> $13						
		echo '1= -'	>> $13
		echo '2= --'	>> $13
		echo '3= -'	>> $13
		echo '4= -'	>> $13
		
		#markerstyle table: 	MATPLOT_MARKERSTYLE
		echo '0= '	>> $15
		echo '1= '	>> $15
		echo '2= '	>> $15
		echo '3= '	>> $15
		echo '4= '	>> $15

		#marker colour table: 	MATPLOT_MARKERCOL
		echo '0= w'	>> $16
		echo '1= r'	>> $16
		echo '2= w'	>> $16
		echo '3= #225588'	>> $16
		echo '3= w'	>> $16
		echo '4= #CC6677'	>> $16

	else if ($use_cb_colors > 5) then

		#linestyle colour table:MATPLOT_COL
		echo '0= k' >> $14
		echo '1= k' >> $14 
		echo '2= #225588' >> $14 
		echo '3= r' >> $14  
		echo '4= #ffe200' >> $14
		echo '5= r' >> $14
		echo '6= k' >> $14
		echo '7= r' >> $14
		echo '8= k' >> $14
		echo '9= #225588' >> $14
		echo '10= #225588' >> $14
		echo '11= r' >> $14
		echo '12= r' >> $14
		echo '14= #4d4d4d' >> $14
		echo '15= #4d4d4d' >> $14
		#Galacticus Blue
		echo '16= #225588' >> $14


		echo '0= -'	>> $13						
		echo '1= -'	>> $13
		echo '2= -.'	>> $13
		echo '3= -'	>> $13
		echo '4= --'	>> $13
		echo '5= -'	>> $13
		echo '6= :'	>> $13
		#echo '4= '	>> $13
		#echo '5= '	>> $13
		#echo '6= '	>> $13
		echo '7= '	>> $13
		echo '8= '	>> $13
		echo '9= '	>> $13
		echo '10= '	>> $13
		echo '11= '	>> $13
		echo '12= '	>> $13
		echo '13= '	>> $13
		echo '14= '	>> $13
		echo '15= '	>> $13
		echo '16= '	>> $13
		if ($use_cb_colors == 5) then
			echo '5= '	>> $13
		else if ($use_cb_colors == 6) then
			echo '0= -'	>> $13						
			echo '1= :'	>> $13
			echo '2= -'	>> $13
			echo '3= --'	>> $13
			echo '4= -'	>> $13
			echo '5= c4'	>> $13
			echo '6= -'	>> $13
			echo '7= -'	>> $13
		else if ($use_cb_colors == 7) then
			echo '7= '	>> $13
		else if ($use_cb_colors == 8) then
			echo '7= -'	>> $13
		else if ($use_cb_colors == 9) then
			echo '7= -'	>> $13
			echo '8= -'	>> $13
		else if ($use_cb_colors == 11) then
			echo '7= -'	>> $13
			echo '8= -.'	>> $13
			echo '9= -'	>> $13
			echo '10= --'	>> $13
		else if ($use_cb_colors > 11) then
			echo '7= -'	>> $13
			echo '8= c1'	>> $13
			echo '9= c4'	>> $13
			echo '10= --'	>> $13
			echo '11= c3'	>> $13
			echo '12= c2'	>> $13
			echo '13= -'	>> $13
			echo '14= '	>> $13
			echo '15= '	>> $13
			echo '16= '	>> $13
		endif

		if ($reduced == True) then
			if ($key =~ plotXY_ww*) then
				echo '0= -'	>> $13						
				echo '1= '	>> $13
				echo '2= '	>> $13
				echo '3= '	>> $13
				echo '4= '	>> $13
				echo '5= -'	>> $13
				echo '6= '	>> $13
				echo '7= -'	>> $13
				echo '8= '	>> $13
				echo '9= '	>> $13
				echo '10= '	>> $13
				echo '11= '	>> $13
				echo '12= '	>> $13
				echo '13= '	>> $13

			else
				echo '0= -'	>> $13						
				echo '1= ' >> $13
				echo '2= '	>> $13
				echo '3= -'	>> $13
				echo '4= '	>> $13
				echo '5= --'	>> $13
				echo '6= '	>> $13
				echo '7= -'	>> $13
				echo '8= c1'	>> $13
				echo '9= '	>> $13
				echo '10= '	>> $13
				echo '11= '	>> $13
				echo '12= '	>> $13
				echo '13= '	>> $13
		endif


		#markerstyle table: 	MATPLOT_MARKERSTYLE
		echo '0= '	>> $15
		echo '1= '	>> $15
		echo '2= '	>> $15
		echo '3= '	>> $15
		echo '4= '	>> $15
		echo '5= '	>> $15
		echo '6= '	>> $15
		echo '7= '	>> $15
		echo '8= '	>> $15
		echo '9= '	>> $15
		echo '10= '	>> $15
		echo '11= '	>> $15

		if ($use_cb_colors == 5) then
			echo '5= o'	>> $15
		else if ($use_cb_colors == 6) then
			echo '6= '	>> $15
			echo '7= .'	>> $15
			echo '7= ^'	>> $15
			echo '8= v'	>> $15
			echo '9= d'	>> $15

		else if ($use_cb_colors == 7) then
			echo '7= '	>> $15
		else if ($use_cb_colors == 8) then
			echo '7= '	>> $15
			echo '8= v'	>> $15
			echo '9= d'	>> $15
			echo '10= o'	>> $15
			echo '11= ^'	>> $15
		else if ($use_cb_colors == 9) then
			echo '8= '	>> $15
			echo '7= '	>> $15
			echo '8= '	>> $15
			echo '9= '	>> $15
			echo '10= '	>> $15
			echo '11= o'	>> $15
			echo '12= ^'	>> $15
			echo '13= d'	>> $15
			echo '14= v'	>> $15
		else if ($use_cb_colors > 9) then
			echo '8= '	>> $15
			echo '9= '	>> $15
			echo '10= '	>> $15
			echo '11= '	>> $15
			echo '12= '	>> $15
			echo '13= '	>> $15
			echo '14= '	>> $15
			echo '15= '	>> $15
			echo '16= '	>> $15
		endif

		#marker colour table: 	MATPLOT_MARKERCOL
		echo '0= w'	>> $16
		echo '1= w'	>> $16
		echo '3= #797979'	>> $16
		echo '4= k'	>> $16
		echo '2= #ffe200'	>> $16
		echo '5= #332288'	>> $16
		echo '6= w'	>> $16
		echo '7= #ffe200'	>> $16
		echo '8= r'	>> $16
		echo '9= k'	>> $16
		echo '10= r'	>> $16
		echo '11= k'	>> $16
		echo '12= r'	>> $16
		echo '13= w'	>> $16
		echo '14= #ffe200'	>> $16
		echo '15= r'	>> $16
		echo '16= k'	>> $16

	else if ($plotXY_key == 'HOD') then
		#linestyle colour table:MATPLOT_COL
		echo '0= k' >> $14
		echo '1= r' >> $14 
		echo '0= #225588' >> $14  
		echo '7= #ffe200' >> $14
		echo '5= r' >> $14
		echo '2= #4d4d4d' >> $14
		echo '15= #4d4d4d' >> $14
		#Galacticus Blue
		echo '13= #225588' >> $14

		#orange: #ff6600
		echo '0= #225588' >> $14 
		echo '1= #225588' >> $14  
		echo '2= #225588' >> $14

		echo '3= #CC6677' >> $14 
		echo '4= #CC6677' >> $14  
		echo '5= #CC6677' >> $14

		if ($hod_plot_key == 'single' || $hod_plot_key == 'CMASS') then

			echo '6= k' >> $14 
			echo '7= k' >> $14  
			echo '8= k' >> $14

		else
			echo '6= #DDCC77' >> $14 
			echo '7= #DDCC77' >> $14  
			echo '8= #DDCC77' >> $14
			echo '9= k' >> $14 
			echo '10= k' >> $14  
			echo '11= k' >> $14
		endif


		echo '0= -'	>> $13						
		echo '1= --'	>> $13
		echo '2= :'	>> $13
		echo '3= -'	>> $13
		echo '4= --'	>> $13
		echo '5= :'	>> $13
		echo '6= -'	>> $13
		echo '7= --'	>> $13
		if ($use_cb_colors == 5) then
			echo '5= '	>> $13
		else if ($use_cb_colors == 6) then
			echo '6= '	>> $13
		else if ($use_cb_colors == 7) then
			echo '7= '	>> $13
		else if ($use_cb_colors == 8) then
			echo '7= -'	>> $13
		endif		
		echo '8= :'	>> $13
		echo '9= -'	>> $13
		echo '10= --'	>> $13
		echo '11= :'	>> $13
		echo '12= -'	>> $13
		echo '13= --'	>> $13
		echo '14= :'	>> $13

		#markerstyle table: 	MATPLOT_MARKERSTYLE
		echo '0= '	>> $15
		echo '1= '	>> $15
		echo '2= '	>> $15
		echo '3= '	>> $15
		echo '4= '	>> $15
		echo '5= '	>> $15
		echo '6= o'	>> $15
		echo '7= .'	>> $15
		if ($use_cb_colors == 5) then
			echo '5= o'	>> $15
		else if ($use_cb_colors == 6) then
			echo '6= o'	>> $15
		else if ($use_cb_colors == 7) then
			echo '7= o'	>> $15
		else if ($use_cb_colors == 8) then
			echo '7= '	>> $15
		endif
		echo '8= .'	>> $15
		echo '9= '	>> $15
		echo '10= '	>> $15
		echo '11= '	>> $15

			echo '6= '	>> $15
			echo '7= '	>> $15
			echo '8= '	>> $15
		endif

		if ($hod_plot_key == 'CMASS') then
			echo '6= .'	>> $15
			echo '7= o'	>> $15
			echo '8= .'	>> $15
		endif

		#marker colour table: 	MATPLOT_MARKERCOL
		echo '0= w'	>> $16
		echo '1= r'	>> $16
		echo '3= #797979'	>> $16
		echo '2= w'	>> $16
		echo '7= #ffe200'	>> $16
		echo '5= #ffe200'	>> $16
		echo '3= #ffe200'	>> $16
		echo '4= #ffe200'	>> $16
		echo '6= #ffe200'	>> $16
		echo '8= #ffe200'	>> $16
		echo '9= #ffe200'	>> $16
		echo '10= #ffe200'	>> $16
		echo '11= #ffe200'	>> $16

	else if ($plotXY_key == 'wp' || $plotXY_key == 'HODsc') then

		set colors = $use_cb_colors
		if ($colors == 2) then
			echo '0= --'	>> $13
			echo '1= :'	>> $13
			echo '2= '	>> $13
			echo '3= '	>> $13
			echo '4= -'	>> $13
			echo '5= '	>> $13
			echo '6= '	>> $13
 
			echo '0= r' >> $14
			echo '1= k' >> $14
			echo '2= #44AA99' >> $14  
			echo '3= #ffe200' >> $14
			echo '4= #ffe200' >> $14
			echo '5= #ffe200' >> $14

			echo '2= o'	>> $15
		else if ($colors == 3) then
			echo '0= -'	>> $13
			echo '1= --'	>> $13
			echo '2= :'	>> $13
			echo '3= '	>> $13
			echo '4= '	>> $13
			echo '5= -'	>> $13
			echo '6= '	>> $13
			echo '7= '	>> $13

			echo '0= #225588' >> $14   
			echo '1= r' >> $14
			echo '2= k' >> $14
			echo '3= #44AA99' >> $14  
			echo '4= #ffe200' >> $14
			echo '5= #ffe200' >> $14
			echo '6= #ffe200' >> $14

			echo '3= o'	>> $15


		else
			echo '0= #225588' >> $14  
			echo '1= #7B369A' >> $14
			echo '2= k' >> $14  
			echo '3= r' >> $14
			echo '4= k' >> $14
			echo '5= #44AA99' >> $14  
			echo '6= #ffe200' >> $14
			echo '7= #ffe200' >> $14

			echo '0= -'	>> $13
			echo '1= -.'	>> $13
			echo '2= -'	>> $13
			echo '3= --'	>> $13
			echo '4= :'	>> $13
			echo '5= -'	>> $13
			echo '6= '	>> $13
			echo '7= -'	>> $13

			echo '3= '	>> $15
		endif

		echo '8= '	>> $13
		echo '9= '	>> $13

		echo '0= '	>> $15
		echo '1= '	>> $15
		echo '2= '	>> $15
		echo '4= '	>> $15
		echo '5= '	>> $15
		echo '6= o'	>> $15
		echo '7= '	>> $15

		echo '0= k'	>> $16
		echo '1= r'	>> $16
		echo '2= w'	>> $16
		echo '3= #ffe200'	>> $16
		echo '4= #ffe200'	>> $16
		echo '5= #ffe200'	>> $16
		echo '6= #ffe200'	>> $16
		echo '7= #ffe200'	>> $16
	else
		set colors = 2
		if ($colors == 2) then
			echo '0= --'	>> $13
			echo '1= -'	>> $13
			echo '2= '	>> $13
			echo '3= '	>> $13
			echo '4= -'	>> $13
			echo '5= '	>> $13
			echo '6= '	>> $13

			echo '0= r' >> $14
			echo '1= k' >> $14
			echo '2= k' >> $14  
			echo '3= #ffe200' >> $14
			echo '4= #ffe200' >> $14
			echo '5= #ffe200' >> $14

			echo '2= o'	>> $15		
			echo '3= o'	>> $15
			echo '3= #ffe200'	>> $16
		else if ($colors == 3) then

			echo '0= -'	>> $13
			echo '1= -'	>> $13
			echo '2= --'	>> $13
			echo '3= :'	>> $13
			echo '4= '	>> $13
			echo '5= -'	>> $13
			echo '6= '	>> $13
			echo '7= '	>> $13

			echo '0= #ffe200' >> $14 
			echo '3= k' >> $14
			echo '1= r' >> $14
			echo '2= k' >> $14  
			echo '4= #ffe200' >> $14
			echo '5= #ffe200' >> $14
			echo '6= #ffe200' >> $14

			echo '4= o'	>> $15
			echo '3= '	>> $15
			echo '2= '	>> $15
			echo '5= .'	>> $15
			echo '4= #ffe200'	>> $16
			echo '5= k'	>> $16
		else
			#Line plots
			#/////////////////////////
			#linestyle colour table:MATPLOT_COL
			echo '0= #225588' >> $14  
			echo '1= r' >> $14
			echo '2= k' >> $14  
			echo '3= k' >> $14
			echo '4= #4baf0c' >> $14
			echo '5= m' >> $14
			echo '6= k' >> $14
			echo '7= k' >> $14
			echo '8= k' >> $14
			echo '9= k' >> $14

			echo '0= -'	>> $13
			echo '1= --'	>> $13
			echo '2= :'	>> $13
			echo '3= '	>> $13
			echo '4= '	>> $13
			echo '5= '	>> $13
			echo '6= '	>> $13
			echo '7= '	>> $13
			echo '8= '	>> $13
			echo '9= '	>> $13
			echo '10= '	>> $13

			#markerstyle table: 	MATPLOT_MARKERSTYLE
			echo '0= '	>> $15
			echo '1= '	>> $15
			echo '2= '	>> $15
			echo '3= o'	>> $15
			echo '4= o'	>> $15
			echo '5= o'	>> $15
			echo '6= o'	>> $15
			echo '7= v'	>> $15
			echo '8= .'	>> $15
			echo '9= o'	>> $15
			echo '10= *'	>> $15

			#marker colour table: 	MATPLOT_MARKERCOL
			echo '0= w'	>> $16
			echo '1= w'	>> $16
			echo '2= w'	>> $16
			echo '3= #ffe200'	>> $16
			echo '4= #ffe200'	>> $16
			echo '6= r'	>> $16
			echo '5= #ffe200'	>> $16
			echo '7= g'	>> $16
			echo '8= k'	>> $16
			echo '9= b'	>> $16
			echo '10= r'	>> $16
		endif


	endif

else if ($OII_ls_set == 'True') then

	#Choose nr of colours-code which should be color-blind (cb) friendly created: choose 0, if you do not want any cb friendly color
	echo 'use_cb_colors= '$use_cb_colors					>> $2
	#linestyle table: 	MATPLOT_LINESTYLE
	echo '0= -.'	>> $13
	echo '1= --'	>> $13
	echo '2= -'	>> $13
	echo '3= '	>> $13
	echo '4= '	>> $13
	echo '5= '	>> $13
	echo '6= '	>> $13
	echo '7= --'	>> $13
	echo '8= -.'	>> $13
	echo '9= :'	>> $13
	echo '10= :'	>> $13

	#linestyle colour table: MATPLOT_COL
	echo '0= w' >> $14 
	echo '1= #c5c1c2' >> $14  
	echo '2= r' >> $14
	echo '4= k' >> $14
	echo '3= g' >> $14

	echo '5= k' >> $14
	echo '6= r' >> $14

	echo '14= b' >> $14
	echo '8= b' >> $14
	echo '9= k' >> $14
	echo '3= k' >> $14
	#dark grey
	echo '9= #212121'	>> $14
	#middle grey
	echo '2= #4d4d4d'	>> $14
	#light grey
	echo '10= #797979'	>> $14

	#marker style table: 	MATPLOT_MARKERSTYLE
	echo '0= '	>> $15
	echo '1= '	>> $15
	echo '2= '	>> $15
	echo '3= o'	>> $15
	echo '4= o'	>> $15
	echo '5= o'	>> $15
	echo '6= ^'	>> $15
	echo '7= '	>> $15
	echo '8= '	>> $15
	echo '9= '	>> $15
	echo '10= '	>> $15


	#marker colour table: 	MATPLOT_MARKERCOL
	echo '0= w'	>> $16
	echo '1= k'	>> $16
	echo '2= w'	>> $16
	echo '3= k'	>> $16
	#echo '3= #ffe200'	>> $16
	echo '4= w'	>> $16
	echo '5= k'	>> $16
	echo '6= r'	>> $16
	echo '7= w'	>> $16
	echo '8= w' >> $16
	echo '9= w'	>> $16
	echo '10= w'	>> $16
	#very light grey
	echo '11= #c5c1c2' >> $16


else

	#Choose nr of colours-code which should be color-blind (cb) friendly created: choose 0, if you do not want any cb friendly color
	#echo 'use_cb_colors= 3'					>> $2
	#linestyle table: 	MATPLOT_LINESTYLE

	if ($catname =~ 'Ga'*) then
		set colors = 3
	else
		set colors = 2
	endif		
	if ($colors == 2) then
		echo '0= :'	>> $13
		echo '1= --'	>> $13
		echo '2= '	>> $13
		echo '3= '	>> $13
		echo '4= -'	>> $13
		echo '5= '	>> $13
		echo '6= '	>> $13

		echo '0= #CC6677' >> $14
		echo '1= k' >> $14
		echo '2= #44AA99' >> $14
		echo '2= k' >> $14    
		echo '3= #ffe200' >> $14
		echo '4= #ffe200' >> $14
		echo '5= #ffe200' >> $14

		echo '0= '	>> $15
		echo '1= '	>> $15
		echo '2= '	>> $15
		echo '3= o'	>> $15
		echo '4= '	>> $15

		echo '2= #ffe200'	>> $16	
		echo '3= k'	>> $16
		echo '4= w'		>> $16

	else if ($colors == 3) then
		echo '0= -'	>> $13
		echo '1= -'	>> $13
		echo '2= -'	>> $13
		echo '3= '	>> $13
		echo '4= '	>> $13
		echo '5= -'	>> $13
		echo '6= '	>> $13
		echo '7= '	>> $13

		echo '0= #225588' >> $14   
		echo '1= r' >> $14
		echo '2= k' >> $14
		echo '3= #ffe200' >> $14  
		echo '4= #ffe200' >> $14
		echo '5= #ffe200' >> $14
		echo '6= #ffe200' >> $14

		echo '4= '	>> $15
		echo '3= '	>> $15
		echo '2= .'	>> $15
		echo '1= '	>> $15
		echo '0= '	>> $15

		echo '2= k'	>> $16	
		echo '3= #ffe200'	>> $16
		echo '5= #ffe200'	>> $16

	else
		echo '0= #225588' >> $14  
		echo '1= #7B369A' >> $14
		echo '2= k' >> $14  
		echo '3= r' >> $14
		echo '4= k' >> $14
		echo '5= #44AA99' >> $14  
		echo '6= #ffe200' >> $14
		echo '7= #ffe200' >> $14

		echo '0= -'	>> $13
		echo '1= -.'	>> $13
		echo '2= -'	>> $13
		echo '3= --'	>> $13
		echo '4= :'	>> $13
		echo '5= -'	>> $13
		echo '6= '	>> $13
		echo '7= -'	>> $13
		echo '8= '	>> $13
		echo '9= '	>> $13

		#linestyle colour table: MATPLOT_COL
		echo '7= #ffe200' >> $14
		echo '8= k' >> $14
	endif

endif

if ($demo == 'True') then
		echo 'use_cb_colors= 0'	>> $2

		echo '0= #CC6677' >> $14
		echo '1= #CC6677' >> $14
		echo '2= #CC6677' >> $14
		echo '3= #225588' >> $14
		echo '4= #225588' >> $14
		echo '5= #225588' >> $14

		echo '0= --'	>> $13						
		echo '1= -'	>> $13
		echo '2= c1'	>> $13
		echo '3= --'	>> $13
		echo '4= -'	>> $13
		echo '5= c1'	>> $13

else
#print polka-dots on Gal-dens all galaxie sample
#echo '0= '	>> $15
#echo '0= w'	>> $16

#echo '1= .'	>> $15
#echo '1= w'	>> $16
echo 'Plot config .... DONE!'
