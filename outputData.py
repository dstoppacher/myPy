import config as conf
import myLib as mL
import numpy as np
import h5py
import time
ts = time.time()
import datetime 
date_time_stamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H:%M:%S')

import matplotlib as mpl
#set this or there will be a backend error for 'Tkinter'! see: http://matplotlib.org/faq/usage_faq.html#what-is-a-backend 17
mpl.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import NullFormatter
nullfmt = NullFormatter()

class OutputData:

    def __init__(
            self,
            config):
        self.obs_color =  '#1e1e28' #ff9900'
        #self.obs_color =  '#ff9900'
        
           
    def multiPlot(self,
                mydata,
                multipanel=False,
                mydir='mydir/',
                myfilename='myfile',
                data_block_offset=0,
                my_add_subplot=[],
                nr_added_subplots=1,
                legend_added_subplots=[],
                add_subplot2=[],
                nr_added_subplots2=1,
                legend_added_subplots2=[],
                data_block_subplot_offset=0,
                config_path='mypath',
                myconfig_datafile='mydata',
                loops=1,
                mykeyword_change1=['keyword', 'default'],
                mykeyword_change2=['keyword', 'default'],
                mykeyword_change3=['keyword', 'default'],
                mykeyword_change4=['keyword', 'default'],
                mykeyword_change5=['keyword', 'default'],
                mykeyword_change6=['keyword', 'default'],
                mykeyword_change7=['keyword', 'default'],
                mykeyword_change8=['keyword', 'default'],
                mykeyword_change9=['keyword', 'default'],
                mykeyword_change10=['keyword', 'default'],
                mykeyword_change11=['keyword', 'default'],
                mykeyword_change12=['keyword', 'default'],
                mykeyword_change13=['keyword', 'default'],
                mykeyword_change14=['keyword', 'default'],                 
                add_plot_par1=['keyword', 'default'],
                add_plot_par2=['keyword', 'default'],
                add_plot_par3=['keyword', 'default'],
                add_plot_par4=['keyword', 'default'],
                add_plot_par5=['keyword', 'default'],
                add_plot_par6=['keyword', 'default'],
                add_plot_par7=['keyword', 'default'],
                add_plot_par8=['keyword', 'default'],
                add_plot_par9=['keyword', 'default'],
                add_plot_par10=['keyword', 'default'],
                add_plot_par11=['keyword', 'default'],
                add_plot_par12=['keyword', 'default'],
                add_plot_par13=['keyword', 'default'],
                add_plot_par14=['keyword', 'default'],                  
                format_y_axis="0.1f",
                minor_ticks_lenght=8,
                major_ticks_lenght=12,
                minor_ticks_width=2,
                major_ticks_width=3,
                plot_key=False,
                catname=False,
                catname_sub=False,
                print_redshift=False,
                print_redshift2=False,
                print_redshift3=False,
                custom_plot_key=None):

        def find_latex_code_times(prefix):
            if myConfig.plot_config_array[prefix+'cut_part2']=='times':
                latex_code_times=r'$\times$'
            elif myConfig.plot_config_array[prefix+'cut_part2']=='halo':    
                latex_code_times=r'$\newmoon$'
            elif myConfig.plot_config_array['cut_part2'].find('stats')!=-1:
                latex_code_times=''
                if myConfig.plot_config_array['cut_part2']=='stats_N':
                    myConfig.plot_config_array['y_title']=r'N$_{\rm M1_{\rm found}}$/N$_{\rm M1_{\rm total}}$-1'
                elif myConfig.plot_config_array['cut_part2']=='stats_xbar_M1found_M1':
                    myConfig.plot_config_array['y_title']=r'$\overline{x}_{\rm M1_{\rm found}}$/$\overline{x}_{\rm M1}}$-1'                
                elif myConfig.plot_config_array['cut_part2']=='stats_xbar_M1found_M2':
                    myConfig.plot_config_array['y_title']=r'$\overline{x}_{\rm M1_{\rm found}}$/$\overline{x}_{\rm M2}}$-1'                
                elif myConfig.plot_config_array['cut_part2']=='stats_xbar_M1_M2':
                    myConfig.plot_config_array['y_title']=r'$\overline{x}_{\rm M1}$/$\overline{x}_{\rm M2}}$-1'
                elif myConfig.plot_config_array['cut_part2']=='stats_xbar_M1-M2_M2':
                    myConfig.plot_config_array['y_title']=r'$\overline{x}_{\rm M1}$/$\overline{x}_{\rm M2}-1$'
                    myConfig.plot_config_array['y_title']=r'$\overline{x}_{\rm M1_{found}}$/$\overline{x}_{\rm M1}-1$'
                    myConfig.plot_config_array['y_title']=r'$\overline{x}_{\rm M1_{found}}$/$\overline{x}_{\rm M2}-1$'
                    myConfig.plot_config_array['y_title']=r'($\overline{x}_{\rm M1}$-$\overline{x}_{\rm M1_{\rm found}}$)/$\overline{x}_{\rm M2}$'
            elif myConfig.plot_config_array[prefix+'cut_part2']!='':
                latex_code_times=''
                if prefix.startswith('xi_'):
                    myConfig.plot_config_array[prefix+'cut_part2']=r'$\overline{\xi}(r) - 1$'
                elif prefix.startswith('wp_'):
                    if myConfig.plot_config_array[prefix+'y_title'].find('Phi')!=-1:
                        myConfig.plot_config_array[prefix+'cut_part2']=r'${\Phi_{\rm DR12}} - 1$'
                    elif myConfig.plot_config_array[prefix+'cut_part2']=='data_res':
                        myConfig.plot_config_array[prefix+'cut_part2']=r'${w_{p_{\rm DR12}}} - 1$'                                               
                    else:
                        myConfig.plot_config_array[prefix+'cut_part2']=r'$\overline{w}(r_p) - 1$'
            
            else:
                latex_code_times=''
            
            print 'here:', latex_code_times
            return latex_code_times


        def plot(myax,
                 mydata,
                 prefix='',
                 add_obs=False,
                 loops=1,
                 print_redshift=False,
                 catname=False):

            def add_subplot(add_subplot,i,j):
                                         
                key_sub1 = myax.twinx()
                key_sub2 = 'axis' + str(i)
                key_sub1.set_zorder(25)             
                plt.setp(key_sub1.get_yticklabels(), visible=False)

                mymarker_code, mycolor_code, mylinestyle, mycolor_map, mymarkerfacecolor, mydashes, colorbar, myalpha, myplotlegend = get_styles(i,i) 
    
                #print add_subplot, legend_added_subplots[str(j)]
                if legend_added_subplots[str(j)].find('MD')!=-1:
                    #print 'here', legend_added_subplots[str(j)]
                    myConfig.plot_config_array['error_bars_y_sub']='yes'
                    myConfig.plot_config_array['subplot_errorbars']='True'
                    myConfig.plot_config_array['filled_between_subplot'] = 'False'
#                else:                 
#                    myConfig.plot_config_array['subplot_errorbars']='True'
#                    myConfig.plot_config_array['filled_between_subplot'] = 'False'
                    
#                if i==9:# or i==6 or i==3:
#                    myConfig.plot_config_array['lw_offset']=float(myConfig.plot_config_array['lw_offset'])-1.5                    
                    
                if legend_added_subplots[str(j)].startswith('Shan'): myConfig.plot_config_array['lw_offset']=4
               
                key_sub2, = key_sub1.plot(add_subplot[:,0], 
                                          add_subplot[:,1], 
                                          color=mycolor_code, 
                                          ls=mylinestyle, 
                                          lw=float(myConfig.plot_config_array['mylw'])+float(myConfig.plot_config_array['lw_offset']))
                
                                
#                if i==8:
#                    myConfig.plot_config_array['mymarkeredgewidth']=float(myConfig.plot_config_array['mymarkeredgewidth'])-1.0
#                    myConfig.plot_config_array['markersize'] = float(myConfig.plot_config_array['markersize'])-3.0 
                
                if mymarker_code!=' ' or mymarker_code!='no':
                    key_sub2 = add_marker(key_sub2, 
                                          key_sub1, 
                                          add_subplot[:,0],
                                          add_subplot[:,1], 
                                          mymarker_code, 
                                          mycolor_code, 
                                          mylinestyle, 
                                          mymarkerfacecolor, 
                                          float(myConfig.plot_config_array['mylw'])+float(myConfig.plot_config_array['lw_offset']),
                                          size=float(myConfig.plot_config_array['markersize']),
                                          width=float(myConfig.plot_config_array['mymarkeredgewidth']))
    
#                if legend_added_subplots[str(j)].startswith('Kormendy') or legend_added_subplots[str(j)].startswith('McC') or legend_added_subplots[str(j)].startswith('Bower') or legend_added_subplots[str(j)].startswith('Maraston') or legend_added_subplots[str(j)].find('BigMD')!=-1:
#                    myConfig.plot_config_array['error_bars_y_sub']='no'
#                else:
#                    myConfig.plot_config_array['error_bars_y_sub']='yes'                    


                if legend_added_subplots[str(j)].startswith('Mad'):
                    myConfig.plot_config_array['error_bars_x_sub']='yes'
                else:
                    myConfig.plot_config_array['error_bars_x_sub']='no'                    

                if legend_added_subplots[str(j)].startswith('Sridhar'):
                    myConfig.plot_config_array['error_bars_y_sub']='yes'
                    
                if myConfig.plot_config_array['subplot_use_cat_colors']=='True':
                    if myConfig.plot_config_array['subplot_is_obs']=='True':
                        mymarker_code = 'o'
                        mycolor_code = self.obs_color 
                        mylinestyle = ' '

                if legend_added_subplots[str(j)].find('g-i')!=-1:
                    myConfig.plot_config_array['subplot_errorbars']='True'
                    myConfig.plot_config_array['filled_between_subplot'] = 'False'  

                if legend_added_subplots[str(j)].startswith('Guo') or \
                    legend_added_subplots[str(j)].startswith('Beh') or \
                    legend_added_subplots[str(j)].startswith('Shan') or \
                    legend_added_subplots[str(j)].startswith('Montero'):
                    #myConfig.plot_config_array['subplot_errorbars']='no'
                    #myConfig.plot_config_array['filled_between_subplot'] = 'True'
                    myConfig.plot_config_array['error_bars_x_sub']='no'
                    
                elif legend_added_subplots[str(j)].startswith('Reid') or \
                     legend_added_subplots[str(j)].startswith('Bigd'):
                         
                    myConfig.plot_config_array['error_bars_y_sub']='no'                         
                else:
                    myConfig.plot_config_array['subplot_errorbars']='True'
                    myConfig.plot_config_array['filled_between_subplot'] = 'False'                
                        
#                if legend_added_subplots[str(j)]=='best fit':
#                    myConfig.plot_config_array['error_bars_y_sub']='no'
                #Errorbars
                add_errorbars(key_sub1,
                              myConfig.plot_config_array['subplot_errorbars'],
                              myConfig.plot_config_array['error_bars_x_sub'], 
                              myConfig.plot_config_array['error_bars_y_sub'],
                              myConfig.plot_config_array['filled_between_subplot'],
                              add_subplot[:,0],
                              add_subplot[:,1],
                              add_subplot[:,2],
                              add_subplot[:,4],
                              xerr_data2=add_subplot[:,3],
                              yerr_data2=add_subplot[:,5],
                              color=mycolor_code,
                              fcolor=mycolor_code,
                              alpha=float(myConfig.plot_config_array['alpha']))

                print 'ADD PLOT SUB j:', j, 'i:', i, 'marker:', mymarker_code, 'mylw:', float(myConfig.plot_config_array['mylw'])+float(myConfig.plot_config_array['lw_offset']), \
                      'fcol:', mymarkerfacecolor, 'myls:', mylinestyle, 'size:', myConfig.plot_config_array['markersize_subplot'], \
                      'errorbars:', myConfig.plot_config_array['subplot_errorbars'], 'filled:', myConfig.plot_config_array['filled_between_subplot']
    
                add_tick_params(key_sub1, 'both', 'minor', length=minor_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']))
                add_tick_params(key_sub1, 'both', 'major', length=major_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']))
                
                if myConfig.plot_config_array['add_axis']=='yes':
                    add_tick_params(key_sub1, 'both', 'major', length=0, width=float(myConfig.plot_config_array['myframe_lw']), set_visible=False)
                else:
                    add_tick_params(key_sub1, 'both', 'major', length=major_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), set_visible=False)  

                try:
                    if myConfig.plot_config_array['add_axis_sub']=='yes' and myConfig.plot_config_array['add_axis_sub_which_ref'+str(j)]=='True':
                        #print count_add_axis, int(myConfig.plot_config_array['add_axis_steps_sub'])
                        #ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS                
                        add_axis(add_subplot, 
                                 'key_add_sub'+str(j), 
                                 key0, 
                                 count_add_axis, 
                                 int(myConfig.plot_config_array['add_which_col_sub']), 
                                 int(myConfig.plot_config_array['add_axis_steps_sub'])) 
                except:
                    print 'sub add axis: key_add_sub'+str(j), 'has no valid axis to add!'                       

                if i>1:
                    try:
                        plt.setp(key1.get_yticklabels(), visible=False) 
                    except:
                        plt.setp(key_sub1.get_yticklabels(), visible=False)
                        
                return key_sub1, key_sub2, legend_added_subplots[str(j)], myplotlegend
 
            def get_styles(a,b):
                dashes=''
                if myConfig.plot_config_array['use_cat_colors']=='True':
                    colorbar='no'
                    if plot_key.find('analyseTargetSelection')!=-1:
                        marker_code       = myConfig.plot_marker_color_map_array[catname+'_marker_code']
                        color_code        = myConfig.plot_marker_color_map_array[catname+'_color_code']
                        linestyle         = myConfig.plot_marker_color_map_array[catname+'_linestyle_code']
                        color_map         = myConfig.plot_marker_color_map_array[catname+'_color_map']
                        markerfacecolor   = myConfig.plot_marker_color_map_array[catname+'_marker_facecolor_code']
                        alpha             = myConfig.plot_marker_color_map_array[catname+'_alpha']                        
    
                    else:
                        marker_code       = myConfig.plot_marker_color_map_array[myConfig.plot_marker_color_map_array['catname'+str(a)]+'_marker_code']
                        color_code        = myConfig.plot_marker_color_map_array[myConfig.plot_marker_color_map_array['catname'+str(a)]+'_color_code']
                        linestyle         = myConfig.plot_marker_color_map_array[myConfig.plot_marker_color_map_array['catname'+str(a)]+'_linestyle_code']
                        color_map         = myConfig.plot_marker_color_map_array[myConfig.plot_marker_color_map_array['catname'+str(a)]+'_color_map']
                        markerfacecolor   = myConfig.plot_marker_color_map_array[myConfig.plot_marker_color_map_array['catname'+str(a)]+'_marker_facecolor_code']
                        alpha             = myConfig.plot_marker_color_map_array[myConfig.plot_marker_color_map_array['catname'+str(a)]+'_alpha']                        
                        
                else:
                    #print 'a:', a, 'b:', b
                    marker_code       = myMarkerTable[str(b)]
                    color_code        = myColourTable[str(b)]
                    linestyle         = myLinestyleTable[str(b)]
                    color_map         = myColormapTable[str(b)]
                    markerfacecolor   = myMarkerfacecolorTable[str(b)]
                    colorbar          = myColorbarTable[str(b)]
                    alpha             = float(myAlphaTable[str(b)])
                    plotlegend        = myPlotLegendTable[str(b)]                    

                #print 'linestyle:', linestyle, 'a:', a, 'b:', b
                if linestyle.startswith('c'):
                    if linestyle.find('1')!=-1:
                        linestyle=(3,(10,3,3,3))
                    elif linestyle.find('2')!=-1:
                        linestyle=(0,(7,2,1,7))
                    elif linestyle.find('3')!=-1:
                        linestyle=(0,(12,2,12,2))
                    else:
                        linestyle=(0,(2,2,2,2))                        
                            
                return marker_code, color_code, linestyle, color_map, markerfacecolor, dashes, colorbar, alpha, plotlegend

            def add_axis(add_which_axis, data2add, key_add, parent_axis, axis_id, count, col, steps):

                print 'add_axis:', key_add, 'parent:', parent_axis, add_which_axis, 'which col:', col, 'axis_id', axis_id, 'axis_count:', count
                mymarker_code, mycolor_code, mylinestyle, mycolor_map, mymarkerfacecolor, mydashes, colorbar, alpha, plotlegend = get_styles(axis_id,axis_id)  

                if add_which_axis=='add_x':
                    #data_add_x= data2add[:,1]#+data2add[:,0]
                    #print data2add[:,0]
                    data2add=data2add[np.where(np.isfinite(data2add[:, int(myConfig.plot_config_array['add_which_col_x'])])==True)[0][:]]                   
                    data_add_x=data2add[:, int(myConfig.plot_config_array['add_which_col_x'])]
                    key_add_x = parent_axis.twiny()

                    key_add_x.spines["top"].set_position(("axes", 1+0.16*count))                    
                    key_add_x.spines["top"].set_linewidth(float(myConfig.plot_config_array['myframe_lw']))
                    key_add_x.set_xscale('linear')
                    

                    xmin_add=np.nanmin(data_add_x)
                    xmax_add=np.nanmax(data_add_x)
                    #key_add_x.plot(data2add[:,0], data_add_x, ls=' ', marker=' ')
                    key_add_x.plot(data2add[:,0], data2add[:,1], ls=' ', marker=' ')

                    #print data_add_x, data2add[:,1]
                    #print 'xmin:', xmin_add, 'xmax:', xmax_add

                else:
                    data2add=data2add[np.where(np.isfinite(data2add[:, int(myConfig.plot_config_array['add_which_col_y'])])==True)[0][:]] 
                    data_add_y=data2add[:, int(myConfig.plot_config_array['add_which_col_y'])]
                    key_add_y = parent_axis.twinx()
                    key_add_y.spines["right"].set_position(("axes", 1+0.14*count))
                    key_add_y.spines["right"].set_linewidth(float(myConfig.plot_config_array['myframe_lw']))
                    key_add_y.set_yscale('linear')                                    
 
                    ymin_add=min(data_add_y)
                    ymax_add=max(data_add_y)               
                    key_add_y.plot(data2add[:,0], data2add[:,1],  ls=' ', marker=' ')

                if add_which_axis=='add_x':               
                    key_add_x.set_xlabel(myConfig.plot_config_array['title_add_x'+str(axis_id)], color='k', fontsize=int(myConfig.plot_config_array['axis_label_fontsize'])-10,labelpad=10)
                    #key_add_x.set_xlim(xmin,xmax)
                    key_add_x.set_xlim(parent_axis.get_xlim())
                    #key_add_x.set_ylim(ymin,ymax)
                else:
                    key_add_y.set_ylabel(myConfig.plot_config_array['title_add_y'+str(axis_id)], color='k', fontsize=int(myConfig.plot_config_array['axis_label_fontsize'])-10, rotation=270, labelpad=25+count*6)
                    #key_add_y.set_ylim(ymin,ymax) 
                    key_add_y.set_ylim(parent_axis.get_ylim())
                
                ticks=[]
                data=[]
                fig.canvas.draw()
                #print data2add[:,col]
                if myConfig.plot_config_array['add_axis_plot_exact_ticks']=='True':
                    a=0
                    while a<data2add[:,0].size:
                        if add_which_axis=='add_x':
                            #print 'a:', a, 'data_ass_x:', data_add_x[a], 'data2add:', data2add[a,0]
                            ticks+=[float("{0:.2f}".format(np.round(data_add_x[a],2)))]
                            #plt.axvline(x=float("{0:.2f}".format(np.round(mydata[a,c*data_block_offset+int(myConfig.plot_config_array['add_which_col_x'])],1))), linewidth=1, color=mycolor_code, ls=mylinestyle)
                            #plt.axvline(x=float("{0:.2f}".format(np.round(mydata[a,c*data_block_offset],1))), linewidth=1, color=mycolor_code, ls=mylinestyle)
                            data+=[float("{0:.2f}".format(np.round(data2add[a,0],2)))]
                            a+=steps
                        else:
                            ticks+=[float("{0:.2f}".format(np.round(data_add_y[a],1)))]
                            #plt.axvline(x=float("{0:.2f}".format(np.round(mydata[a,c*data_block_offset+int(myConfig.plot_config_array['add_which_col_y'])],1))), linewidth=1, color=mycolor_code, ls=mylinestyle)
                            #plt.axvline(x=float("{0:.2f}".format(np.round(mydata[a,c*data_block_offset],1))), linewidth=1, color=mycolor_code, ls=mylinestyle)
                            data+=[float("{0:.2f}".format(np.round(data2add[a,1],1)))]                        
                            a+=steps
                     
                    #print 'ticks:', ticks
                    #print 'data:', data
                
                if add_which_axis=='add_x':
                    
                    key_add_x.set_xticks(data)
                    if ticks!=[]:
                        key_add_x.set_xticklabels(ticks, color='k', fontsize=int(myConfig.plot_config_array['axis_label_fontsize']))                     

#                    for item in key_add_x.get_xticklabels():
#                        print item
#                    for item in key_add_x.get_xticks():
#                        print item

                    if myConfig.plot_config_array['add_axis_plot_exact_ticks']=='True':                        
                        add_tick_params(key_add_x, 
                                        'both', 
                                        'minor', 
                                        bottom='off', 
                                        top='off', 
                                        length=minor_ticks_lenght, 
                                        width=float(myConfig.plot_config_array['myframe_lw']), 
                                        set_visible=True, 
                                        label_size_offset=8, 
                                        pad=6)
                        add_tick_params(key_add_x, 
                                        'both', 
                                        'major', 
                                        bottom='off', 
                                        top='on', 
                                        length=major_ticks_lenght,
                                        width=float(myConfig.plot_config_array['myframe_lw']),
                                        set_visible=True,
                                        label_size_offset=8,
                                        pad=6) 

                    else:
                        add_tick_params(key_add_x,
                                        'both', 
                                        'minor', 
                                        isadd_x=True,
                                        xmin_add=xmin_add, 
                                        xmax_add=xmax_add, 
                                        bottom='off', 
                                        top='off', 
                                        length=minor_ticks_lenght,
                                        width=float(myConfig.plot_config_array['myframe_lw']), 
                                        set_visible=True, 
                                        label_size_offset=8, 
                                        pad=6)
                        add_tick_params(key_add_x, 
                                        'both', 
                                        'major', 
                                        isadd_x=True,
                                        xmin_add=xmin_add, 
                                        xmax_add=xmax_add, 
                                        bottom='off', 
                                        top='on', 
                                        length=major_ticks_lenght, 
                                        width=float(myConfig.plot_config_array['myframe_lw']), 
                                        set_visible=True, 
                                        label_size_offset=8, 
                                        pad=6) 
        
                else:
                    key_add_y.set_yticks(data)
                    if ticks!=[]:
                        key_add_y.set_yticklabels(ticks, color='k', fontsize=int(myConfig.plot_config_array['axis_label_fontsize']))

#                    for item in key_add_y.get_yticklabels():
#                        print item
#                    for item in key_add_y.get_yticks():
#                        print item

                    if myConfig.plot_config_array['add_axis_plot_exact_ticks']=='True':
                        add_tick_params(key_add_y, 
                                        'both', 
                                        'minor', 
                                        left='off', 
                                        right='off', 
                                        length=minor_ticks_lenght, 
                                        width=float(myConfig.plot_config_array['myframe_lw']), 
                                        set_visible=True, 
                                        label_size_offset=8, 
                                        pad=6)
                        add_tick_params(key_add_y, 
                                        'both', 
                                        'major', 
                                        left='off', 
                                        right='on', 
                                        length=major_ticks_lenght, 
                                        width=float(myConfig.plot_config_array['myframe_lw']), 
                                        set_visible=True, 
                                        label_size_offset=8, 
                                        pad=6)
                       
                    else:
                        add_tick_params(key_add_y, 
                                        'both', 
                                        'minor', 
                                        isadd_y=True, 
                                        ymin_add=ymin_add, 
                                        ymax_add=ymax_add, 
                                        left='off', 
                                        right='off', 
                                        length=minor_ticks_lenght, 
                                        width=float(myConfig.plot_config_array['myframe_lw']), 
                                        set_visible=True, 
                                        label_size_offset=8, 
                                        pad=6)
                        add_tick_params(key_add_y, 
                                        'both', 
                                        'major', 
                                        isadd_y=True, 
                                        ymin_add=ymin_add, 
                                        ymax_add=ymax_add, 
                                        left='off', 
                                        right='on', 
                                        length=major_ticks_lenght, 
                                        width=float(myConfig.plot_config_array['myframe_lw']), 
                                        set_visible=True, 
                                        label_size_offset=8, 
                                        pad=6)
                                  
                      
            def add_marker(key, 
                           ax, 
                           x_data, 
                           y_data, 
                           marker, 
                           color, 
                           ls, 
                           mfc, 
                           mylw, 
                           size=float(myConfig.plot_config_array['markersize']), 
                           width=float(myConfig.plot_config_array['mymarkeredgewidth']),
                           alpha=1):
                          
                key, = ax.plot(x_data, 
                              y_data, 
                              ls=ls, 
                              color=color,
                              marker=marker, 
                              lw=mylw, 
                              markersize=size,
                              markeredgewidth=width,
                              markerfacecolor=mfc,
                              alpha=alpha)
                return key
    
    
            def add_errorbars(ax, 
                              errorbars, 
                              x_bar, y_bar, 
                              filled_between, 
                              x_data, 
                              y_data, 
                              xerr_data1, 
                              yerr_data1, 
                              xerr_data2=[], 
                              yerr_data2=[], 
                              color='k', 
                              fcolor='k', 
                              alpha=0.3,
                              ls='-'):               
                print 'add_errorbars():', errorbars, 'x_bar:', x_bar, 'y_bar:', y_bar, 'filled_between:', filled_between

                if x_bar=='yes':
                    if errorbars=='True':
                        ax.errorbar(x_data,
                                    y_data, 
                                    xerr=xerr_data1, 
                                    color='k',
                                    ls='')                                      
                if y_bar=='yes':
                    if errorbars=='True': 
                        ax.errorbar(x_data,
                                    y_data,
                                    yerr=[yerr_data1,yerr_data2], 
                                    color='k',
                                    ls='')

                #import mpl.transforms as mtransforms
                trans = mpl.transforms.blended_transform_factory(ax.transData, ax.transAxes)
                                    
                if filled_between=='True' and x_bar=='yes':

                    # use the data coordinates for the x-axis and the axes coordinates for the y-axis
                    #import mpl.transforms as mtransforms
                    trans = mpl.transforms.blended_transform_factory(ax.transData, ax.transAxes)

                    ax.fill_between(np.linspace(xerr_data1[0],xerr_data2[0], len(y_data)), 
                                    y_data.min(),
                                    y_data.max(),
                                    color=color, 
                                    facecolor=fcolor, 
                                    alpha=alpha, 
                                    transform=trans,
                                    interpolate=True,
                                    lw=0.0,
                                    linestyle=ls)                     

                if filled_between=='True' and y_bar=='yes':
                    if myConfig.plot_config_array['fill_no_facecolor']=='True':
                        ax.fill_between(x_data, 
                                        yerr_data2, 
                                        yerr_data1, 
                                        color=color, 
                                        facecolor='none', 
                                        alpha=0.5,
                                        interpolate=False,
                                        lw=3.0,
                                        linestyle=ls)
                    else:
                        ax.fill_between(x_data, 
                                        yerr_data2, 
                                        yerr_data1, 
                                        color=color, 
                                        facecolor=fcolor, 
                                        alpha=alpha,
                                        interpolate=False,
                                        lw=1.5,
                                        linestyle='-')                    
                    
    
            def add_tick_params(ax, 
                                which_axis, 
                                which_ticks, 
                                length=12, 
                                width=float(myConfig.plot_config_array['myframe_lw']), 
                                bottom='on', 
                                top='on', 
                                right='on', 
                                left='on',
                                pad=10, 
                                set_visible=False, 
                                label_size_offset=0, 
                                isadd_x=False, 
                                xmin_add=None, 
                                xmax_add=None, 
                                isadd_y=False, 
                                ymin_add=None, 
                                ymax_add=None,
                                isshare_x=False, 
                                xmin_share=None, 
                                xmax_share=None, 
                                isshare_y=False, 
                                ymin_share=None, 
                                ymax_share=None,
                                color='k',
                                parent_axis=None):


                #print 'ax:', ax, 'which_axis:', which_axis, 'which_ticks:', which_ticks, 'bottom:', bottom, 'top:', top, 'left:', left, 'right:', right, 'set_visible:', set_visible

                ax.patch.set_visible(False) 
                ax.tick_params(axis=which_axis,
                               which=which_ticks,
                               length=length, 
                               width=width,
                               bottom=bottom,
                               right=right,
                               left=left,
                               top=top, 
                               labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize'])-label_size_offset, 
                               pad=pad,
                               color=color)
                
                if myConfig.plot_config_array['plot_grid']=='yes':        
                    plt.grid(True)
                else:
                    plt.grid(False)                           

                if myConfig.plot_config_array['log_scale_x']=='symlog' and isadd_x==False and isshare_x==False:
                    ax.set_xscale('symlog')
                elif myConfig.plot_config_array['log_scale_x']=='yes' and isadd_x==False and isshare_x==False:
                    ax.set_xscale('log')
                elif myConfig.plot_config_array['log_scale_share_x']=='yes' and isadd_x==False and isshare_x==True:
                    ax.set_xscale('log')                      
#                    
                if myConfig.plot_config_array['log_scale_y']=='yes' and isadd_y==False and isshare_y==False:
                    ax.set_yscale('log')

                elif myConfig.plot_config_array['log_scale_share_y']=='yes' and isadd_y==False and isshare_y==True:
                    ax.set_yscale('log')
                    if myConfig.plot_config_array['log_scale_y']=='yes':
                        myax.set_yscale('log')
                    else:
                        myax.set_yscale('linear')

                #print 'xmin:', xmin, 'xmax:', xmax, 'ymin:', ymin, 'ymax:', ymax
                #print 'isadd_x:', isadd_x, 'isadd_y:', isadd_y
                #print 'xmin_share:', xmin_share, 'xmax_share:', xmax_share, 'ymin_share:', ymin_share, 'ymax_share:', ymax_share
                #print 'isshare_x:', isshare_x, 'isshare_y:', isshare_y
                if isadd_x==True:
                    plt.xlim(xmin_add,xmax_add)
                    plt.ylim(ymin,ymax) 
                elif isadd_y==True:
                    plt.xlim(xmin,xmax)
                    plt.ylim(ymin_add,ymax_add)
                elif isshare_x==True:
                    plt.xlim(xmin_share,xmax_share)
                    plt.ylim(ymin,ymax) 
                elif isshare_y==True:
                    plt.xlim(xmin,xmax)
                    plt.ylim(ymin_share,ymax_share)                     
                else:
                    plt.xlim(xmin,xmax)
                    plt.ylim(ymin,ymax)            
                try:               
                    plt.zlim(myConfig.plot_config_array[prefix+'y_range_min'],myConfig.plot_config_array[prefix+'y_range_max'])
                except:
                    pass  

                plt.setp(ax.get_xticklabels(), visible=set_visible)                
                plt.setp(ax.get_yticklabels(), visible=set_visible)
      

            def plot_lines_min_max_population(H,X,Y,
                                              binsize_x, 
                                              binsize_y,
                                              plot_min=False,
                                              color='k', 
                                              linestyle='-', 
                                              alpha=0.3):

                print 'plot line at vertical/horizontal position at the',                    
                if plot_min!=False:
                    print 'MINIMUM',
                    x,y = np.unravel_index(np.argmin(H),H.shape) 
                else:
                    print 'MAXIMUM',
                    x,y = np.unravel_index(np.argmax(H),H.shape)                         
                print 'of the population!'
                print 'coordinates to located:\nH[x,y]:', H[x][y], 'x/y coordinate:', x, '/', y, 'x/y values:', X[x,y], '/', Y[x,y]

                line_arr =np.zeros((10,7), dtype=np.float32)              
                line_arr[:,0]=X[x,y]
                line_arr[:,1]=Y[x,y]
                line_arr[:,2]=X[x,y]-binsize_x/2.
                line_arr[:,3]=X[x,y]+binsize_x/2. 
                line_arr[:,4]=Y[x,y]-binsize_y/2.  
                line_arr[:,5]=Y[x,y]+binsize_y/2.                 
                line_arr[:,6]=np.linspace(-100.0,100.0,10)

                myax.plot(line_arr[:,0], line_arr[:,6], linestyle=linestyle, color=color, lw=1.5)                  
                myax.plot(line_arr[:,6], line_arr[:,1], linestyle=linestyle, color=color, lw=1.5)                

                add_errorbars(myax,
                              'False',                          
                              'yes', 
                              'no', 
                              'True',
                              line_arr[:,0],
                              line_arr[:,6],
                              line_arr[:,2],                                  
                              line_arr[:,2],
                              xerr_data2=line_arr[:,3],
                              yerr_data2=line_arr[:,3],                              
                              color=color,
                              fcolor=color,
                              alpha=alpha)

                add_errorbars(myax,
                              'False',                          
                              'no', 
                              'yes', 
                              'True',
                              line_arr[:,6],
                              line_arr[:,1],
                              line_arr[:,2],                                  
                              line_arr[:,4],
                              xerr_data2=line_arr[:,3],
                              yerr_data2=line_arr[:,5],                              
                              color=color,
                              fcolor=color,
                              alpha=alpha)
 
            def histo_number_density(data_x,
                                     data_y,
                                     weights,
                                     orientation='vertical',
                                     mycolor='k',
                                     myfacecolor='k',                                     
                                     myls='-',
                                     myalpha=0.5,
                                     mycolor_map='no'):
                print 'HISTO NUMER DENSITY:',                
                # definitions for the axes
                left = 0.16
                width=0.8-left
                bottom=0.12
                height=0.80-bottom
                bottom_h = left_h = left + width
                
                rect_histx = [left, bottom_h, width, 0.18]
                rect_histy = [left_h, bottom, 0.18, height]
                                
                axHistx = plt.axes(rect_histx)
                axHisty = plt.axes(rect_histy)
                
                for axis in ['top','bottom','left','right']:
                    axHistx.spines[axis].set_linewidth(float(myConfig.plot_config_array['myframe_lw']))
                    axHisty.spines[axis].set_linewidth(float(myConfig.plot_config_array['myframe_lw']))
                    axHistx.spines[axis].set_linewidth(float(myConfig.plot_config_array['myframe_lw']))
                    axHisty.spines[axis].set_linewidth(float(myConfig.plot_config_array['myframe_lw']))                    

                if mycolor_map!='no':
                    mycmap = mpl.cm.get_cmap(mycolor_map)
                    myfacecolor = mycmap(0.1)
                
                # no labels
                axHistx.xaxis.set_major_formatter(nullfmt)
                axHisty.yaxis.set_major_formatter(nullfmt)
                                
                print 'x min/max:', xmin, '/', xmax, 'y min/max:', ymin, '/', ymax, 'mycolor:', mycolor, 'myls:', myls, 'myalpha:', myalpha,

                myhisto_ticks=np.linspace(float(myConfig.plot_config_array['histo_num_density_min']),float(myConfig.plot_config_array['histo_num_density_max']),int(myConfig.plot_config_array['histo_num_density_number_major_ticks']))                

                #Top panel
                #------------------------------------------------#

                binsx = np.arange(xmin, xmax + (xmax-xmin)/int(myConfig.plot_config_array['histo_num_density_nbins_top']), (xmax-xmin)/int(myConfig.plot_config_array['histo_num_density_nbins_top']))
                #print 'binsx:', binsx, 'norm:', float(len(data_x))

                myhisto_weights = weights/float(len(data_x))
                
                axis = axHistx.hist(data_x, bins=binsx, alpha=myalpha, edgecolor=mycolor, facecolor=myfacecolor, ls=myls, lw=2.5, normed=None, weights=myhisto_weights)

                axHistx.set_xlim(xmin, xmax)
                axHistx.set_ylim(float(myhisto_ticks[0]), float(myhisto_ticks[-1]))

                print 'myhisto_ticks min/max:', float(myhisto_ticks[0]), float(myhisto_ticks[-1])              

                xticks=[]
                if myConfig.plot_config_array['xticks']!='':
                    myticks=mL.multicolTestAlgorithm(myConfig.plot_config_array['xticks'])
                    k=0
                    while k<len(myticks):
                        xticks+=[float(myticks[k])+float(myConfig.plot_config_array['xticks_minor_offset'])]
                        k+=1                               

                    axHistx.set_xticks(xticks)
                    
                axHistx.set_yticks(myhisto_ticks)        
                axHistx.set_yticklabels(myhisto_ticks, fontsize=14)
                #plt.setp(axHistx.get_yticklabels(), visible=False)
                #plt.setp(axHistx.get_yticklabels()[-1], visible=False)
                              
                axHistx.tick_params(axis='y', which='minor', direction='in', length=minor_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), zorder=20)
                axHistx.tick_params(axis='y', which='major', direction='in', length=major_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), zorder=20) 
                axHistx.tick_params(axis='x', which='minor', bottom='off', top='on', direction='in', length=minor_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), zorder=20)
                axHistx.tick_params(axis='x', which='major', bottom='off', top='on', direction='in', length=major_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), zorder=20)                


                minor_x = MultipleLocator(float(myConfig.plot_config_array['minor_ticks_x_space']))
                axHistx.xaxis.set_minor_locator(minor_x)
             
                minor_y = MultipleLocator(float(myhisto_ticks[1]/2))
                axHistx.yaxis.set_minor_locator(minor_y)
                axHistx.set_zorder(20)

                axHistx.set_ylabel(r'$f_{N_{\rm gal}}$', fontsize=int(myConfig.plot_config_array['axis_label_fontsize'])-7)
                if myhisto_ticks[0]>1e3:                
                    axHistx.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{0:.1e}".format(x)))
                else:
                    axHistx.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{0:.2f}".format(x)))

                #ssfr border
                #axHistx.axvline(x=-10.8, ymin=-100, ymax=100, color='k', ls='--', lw=2.0)
                #axHistx.axvline(x=-11.8, ymin=-100, ymax=100, color='k', ls='--', lw=2.0) 
               
                
                #Left panel                
                #------------------------------------------------#

                binsy = np.arange(ymin, ymax + (ymax-ymin)/int(myConfig.plot_config_array['histo_num_density_nbins_left']), (ymax-ymin)/int(myConfig.plot_config_array['histo_num_density_nbins_left']))

                myhisto_weights = weights/float(len(data_y))
                
                axHisty.hist(data_y, bins=binsy, orientation='horizontal', alpha=myalpha, edgecolor=mycolor, facecolor=myfacecolor, ls=myls, lw=2.5, normed=None, weights=myhisto_weights)
                axHisty.set_ylim(ymin, ymax)
                axHisty.set_xlim(myhisto_ticks[0], myhisto_ticks[-1])

                axHisty.set_xticks(myhisto_ticks)        
                axHisty.set_xticklabels(myhisto_ticks, fontsize=14)
                axHisty.tick_params(axis='x', which='minor', direction='in', length=minor_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), zorder=20)
                axHisty.tick_params(axis='x', which='major', direction='in', length=major_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), zorder=20)                 
                axHisty.tick_params(axis='y', which='minor', left='off', right='on', direction='in', length=minor_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), zorder=20)
                axHisty.tick_params(axis='y', which='major', left='off', right='on', direction='in', length=major_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), zorder=20)                               

                minor_x = MultipleLocator(float(myhisto_ticks[1]/2))
                axHisty.xaxis.set_minor_locator(minor_x)

                yticks=[]
                if myConfig.plot_config_array['yticks']!='':
                    myticks=mL.multicolTestAlgorithm(myConfig.plot_config_array['yticks'])
                    k=0
                    while k<len(myticks):
                        yticks+=[float(myticks[k])+float(myConfig.plot_config_array['yticks_minor_offset'])]
                        k+=1                               

                    axHisty.set_yticks(yticks)
 
                
                #plt.setp(axHisty.get_xticklabels()[-1], visible=False)

               
                minor_y = MultipleLocator(float(myConfig.plot_config_array['minor_ticks_y_space']))
                axHisty.yaxis.set_minor_locator(minor_y)
                axHisty.set_zorder(20)
                
                plt.setp(axHisty.get_xticklabels(), visible=True, rotation=-90)               
                #plt.setp(axHisty.get_xticklabels(), visible=False)
                axHisty.set_xlabel(r'$f_{N_{\rm gal}}$', fontsize=int(myConfig.plot_config_array['axis_label_fontsize'])-7)
                if myhisto_ticks[0]>1e3:                
                    axHisty.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{0:.1e}".format(x)))
                else:
                    axHisty.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{0:.2f}".format(x)))                    
                    
                #axHisty.axhline(y=-10.8, xmin=-100, xmax=100, color='k', ls='--', lw=2.0)
                #axHisty.axhline(y=-11.8, xmin=-100, xmax=100, color='k', ls='--', lw=2.0)
                    
            def plot_colorbar(data,
                              offset,
                              myticks_min,
                              myticks_max,
                              myticks_steps=1):

                myticks=np.arange(myticks_min,myticks_max+myticks_steps, myticks_steps)
                
                cbaxis = fig.colorbar(data, use_gridspec=True, pad=-0.03, aspect=10, orientation='horizontal', shrink=1, format='%i', ticks=myticks)
                cbaxis.ax.minorticks_on()
                cbaxis.ax.tick_params(labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize'])-5, which='minor', length=minor_ticks_lenght-1.5, width=float(myConfig.plot_config_array['myframe_lw'])-0.5)
                cbaxis.ax.tick_params(labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize'])-5, which='major', length=minor_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']))
                
                cbaxis.ax.set_position([float(myConfig.plot_config_array['colorbar_anchor_right'])-0.04*len(myticks),0.15+0.1*offset,0.04*len(myticks),0.05])
                cbaxis.ax.set_zorder(30)                

            def Surface3DPlot(i,j):

                N=int(myConfig.plot_config_array['contour_histo_nbins']) 
                if len(mydata['name'+str(j)+'_data'][mydata['col_name_x']])<5e6:
                    N=70
                
                mymarker_code, mycolor_code, mylinestyle, mycolor_map, mymarkerfacecolor, mydashes, colorbar, myalpha, myplotlegend = get_styles(i,j) 
 
                print 'SURFACE 3D-PLOT!'               
                print 'j:', j, 'i:', i, 'marker:', mymarker_code, 'mylw:', float(myConfig.plot_config_array['mylw'])+float(myConfig.plot_config_array['lw_offset']), \
                      'fcol:', mymarkerfacecolor, 'myls:', mylinestyle, 'size:', myConfig.plot_config_array['markersize_subplot'], \
                      'errorbars:', myConfig.plot_config_array['subplot_errorbars'], 'filled:', myConfig.plot_config_array['filled_between_subplot'], 'nbins:', N

                print mydata['col_name_x'], mydata['col_name_y']
                print mydata['name'+str(c)+'_data'][mydata['col_name_x']], mydata['name'+str(c)+'_data'][mydata['col_name_y']]
                x=mydata['name'+str(c)+'_data'][mydata['col_name_x']]
                y=mydata['name'+str(c)+'_data'][mydata['col_name_z']]
                #x = np.linspace(min(mydata['name'+str(j)+'_data'][mydata['col_name_x']]), max(mydata['name'+str(j)+'_data'][mydata['col_name_x']]), N)
                #y = np.linspace(min(mydata['name'+str(j)+'_data'][mydata['col_name_z']]), max(mydata['name'+str(j)+'_data'][mydata['col_name_z']]), N)
               
                #indexing='ij' to use matrix indexing in order to get a meshgrid which is corresponding to a Matrix of (X,Y) --> see numpy documentatcion about meshgrids
                #X,Y = np.meshgrid(x,y, indexing='ij')
                
                #X,Y = np.meshgrid(mydata['name'+str(j)+'_data'][mydata['col_name_x']],mydata['name'+str(j)+'_data'][mydata['col_name_z']], indexing='ij')
                Z = mydata['name'+str(j)+'_data'][mydata['col_name_z']] - mydata['name'+str(j)+'_data'][mydata['col_name_x']]
                #Z = np.linspace(min(Z), max(Z), N)
                #H, xedges, yedges = np.histogram2d(mydata['name'+str(j)+'_data'][mydata['col_name_x']], mydata['name'+str(j)+'_data'][mydata['col_name_z']], bins=N)     
         
                # Plot the surface.
                fig=plt.figure()
                ax0=fig.gca(projection='3d')             
                surf= ax0.plot_surface(mydata['name'+str(j)+'_data'][mydata['col_name_x']], mydata['name'+str(j)+'_data'][mydata['col_name_z']], Z, cmap='coolwarm', linewidth=0, antialiased=True)



                ax0 = plt.axes(projection='3d')
                ax0.plot_trisurf(x, y, Z, cmap='viridis', edgecolor='none');



                fig.colorbar(surf, shrink=0.5, aspect=5)                
                # Customize the z axis.
                #ax0.set_xlim(min(x), max(x))
                #ax0.set_ylim(min(y), max(y))
                #ax0.set_zlim(Z.min(), Z.max())                
                #sp3D.zaxis.set_major_locator(LinearLocator(1))
                #sp3D.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))                           

                pp = PdfPages(mydir+myfilename+'.pdf')
                if myConfig.plot_config_array['tight_layout']=='True':
                    mybbox_inches='tight'
                else:
                    mybbox_inches=None
                plt.savefig(pp, format='pdf', rasterized=True, dpi=float(myConfig.plot_config_array['quality_dpi']), pad_inches=0.05, bbox_inches=mybbox_inches, transparent=True)
                pp.close()
                
                return myax,  myax.plot([], color=mycolor_code, linewidth=2.0, alpha=1, ls=mylinestyle)                

            def ContourPlot(i,j):
            #CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR CONTOUR                    

                mymarker_code, mycolor_code, mylinestyle, mycolor_map, mymarkerfacecolor, mydashes, colorbar, myalpha, myplotlegend = get_styles(i,j) 

                N=int(myConfig.plot_config_array['contour_histo_nbins']) 
                if len(mydata['name'+str(j)+'_data'][mydata['col_name_x']])<5e6:
                    N=70
  
                print 'CONTOUR-PLOT!'               
                print 'j:', j, 'i:', i, 'marker:', mymarker_code, 'mylw:', float(myConfig.plot_config_array['mylw'])+float(myConfig.plot_config_array['lw_offset']), \
                      'fcol:', mymarkerfacecolor, 'myls:', mylinestyle, 'size:', myConfig.plot_config_array['markersize_subplot'], \
                      'errorbars:', myConfig.plot_config_array['subplot_errorbars'], 'filled:', myConfig.plot_config_array['filled_between_subplot'], 'nbins:', N

                #print mydata['col_name_x'], mydata['col_name_y']
                #print mydata['name'+str(j)+'_data'][mydata['col_name_x']], mydata['name'+str(j)+'_data'][mydata['col_name_y']]
                print mydata['name'+str(j)+'_data'][mydata['col_name_weights']]      

                if j==2:
                    x = np.linspace(min(mydata['name'+str(j)+'_data'][mydata['col_name_x']]), max(mydata['name'+str(j)+'_data'][mydata['col_name_x']]), N)
                else:
                    x = np.linspace(min(mydata['name'+str(j)+'_data'][mydata['col_name_x']]), max(mydata['name'+str(j)+'_data'][mydata['col_name_x']]), N)

                y = np.linspace(min(mydata['name'+str(j)+'_data'][mydata['col_name_y']]), max(mydata['name'+str(j)+'_data'][mydata['col_name_y']]), N)
                #print x
                #print y
                H, xedges, yedges = np.histogram2d(mydata['name'+str(j)+'_data'][mydata['col_name_x']], 
                                                   mydata['name'+str(j)+'_data'][mydata['col_name_y']], 
                                                   bins=N, 
                                                   weights=mydata['name'+str(j)+'_data'][mydata['col_name_weights']])       

                binsize_x = (xedges.max()-xedges.min())/N
                binsize_y = (yedges.max()-yedges.min())/N                
                #indexing='ij' to use matrix indexing in order to get a meshgrid which is corresponding to a Matrix of (X,Y) --> see numpy documentatcion about meshgrids
                X,Y = np.meshgrid(x,y, indexing='ij')
                                                    
                #confidence levels / contour levels
                if len(mydata['name'+str(j)+'_data'][mydata['col_name_x']])<1e6:
                    mylevels=[#H.max()/100.0*0.1,
                              H.max()/100.0*13.6,
                              H.max()/100.0*31.74,
                              H.max()/100.0*68.26,
                              H.max()/100.0*95.0,
                              H.max()/100.0*99.7]
                    
                elif len(mydata['name'+str(j)+'_data'][mydata['col_name_x']])>1e6:
                    mylevels=[H.max()/100.0*2.1,
                              H.max()/100.0*13.6,
                              H.max()/100.0*31.74,
                              H.max()/100.0*68.26,
                              H.max()/100.0*90.0,
                              H.max()/100.0*95.0,
                              H.max()/100.0*99.7]
                              
                else: 
                    
                    mylevels=[H.max()/100.0*2.1,
                              H.max()/100.0*13.6,
                              H.max()/100.0*31.74,
                              H.max()/100.0*68.26,
                              H.max()/100.0*95.0,
                              H.max()/100.0*99.7]

                    
                    #level according to 68-95-99.7 rule          
#                    mylevels=[H.max()/100.0*68,
#                              H.max()/100.0*95,
#                              H.max()/100.0*99.7]                              

                print 'max:', H.max(), 'min:', H.min()#, 'mylevels:', mylevels
                print 'levels:', mylevels
                #print xedges
                #print yedges
         
                CS=key0.contour(X,Y,H, len(mylevels), linewidths=3.0, linestyles=mylinestyle, colors=mycolor_code, levels=mylevels)
                #plt.clabel(CS, colors = 'k', fmt = '%0.1f', fontsize=12)
                if mycolor_map!='no':
                    CS=key0.contourf(X,Y,H, len(mylevels), levels=mylevels, cmap=mycolor_map, alpha=myalpha)     


                if myConfig.plot_config_array['plot_lines_min_max_population']=='yes':
                    plot_lines_min_max_population(H,X,Y, binsize_x, binsize_y, color=mycolor_code)
                
                histo_number_density(mydata['name'+str(j)+'_data'][mydata['col_name_x']],
                                     mydata['name'+str(j)+'_data'][mydata['col_name_y']],
                                     mydata['name'+str(j)+'_data'][mydata['col_name_weights']],                                     
                                     mycolor=mycolor_code,
                                     myfacecolor=mymarkerfacecolor,
                                     myls=mylinestyle,
                                     myalpha=myalpha,
                                     mycolor_map=mycolor_map)
                
                if myConfig.plot_config_array['add_scalebar']=='True':
                    scalb.add_scalebar(key0, 
                                       sizex=float(myConfig.plot_config_array['scalebar_width_x']), 
                                       sizey=float(myConfig.plot_config_array['scalebar_width_y']), 
                                       hidex=False, hidey=False,
                                       matchx=False, matchy=False,
                                       loc=int(myConfig.plot_config_array['scalebar_loc_code']), 
                                       labelx=myConfig.plot_config_array['scalebar_width_x'], 
                                       labely=myConfig.plot_config_array['scalebar_width_y'],
                                       pad=0.1, borderpad=0.5, sep=5)
                
                axis0, = myax.plot([], color=mycolor_code, linewidth=3.0, alpha=1, ls=mylinestyle)
                return key0, axis0, mydata['name'+str(j)+'_legend'], myplotlegend 
                       
                
            def BarHistogram(i,j):
            #BARPLOT BARPLOT BARPLOT BARPLOT BARPLOT BARPLOT BARPLOT BARPLOT BARPLOT BARPLOT BARPLOT BARPLOT BARPLOT BARPLOT BARPLOT BARPLOT BARPLOT BARPLOT BARPLOT BARPLOT BARPLOT

                mymarker_code, mycolor_code, mylinestyle, mycolor_map, mymarkerfacecolor, mydashes, colorbar, myalpha, myplotlegend = get_styles(i,j)

                N=int(myConfig.plot_config_array['contour_histo_nbins']) 
                if len(mydata['name'+str(j)+'_data'][mydata['col_name_x']])<5e6:
                    N=70                
                
                print 'BAR HISTOGRAMM PLOT!'
                print 'j:', j, 'i:', i, 'marker:', mymarker_code, 'mylw:', float(myConfig.plot_config_array['mylw'])+float(myConfig.plot_config_array['lw_offset']), \
                      'fcol:', mymarkerfacecolor, 'myls:', mylinestyle, 'size:', myConfig.plot_config_array['markersize_subplot'], \
                      'errorbars:', myConfig.plot_config_array['subplot_errorbars'], 'filled:', myConfig.plot_config_array['filled_between_subplot'], 'nbins:', N                
                myConfig.plot_config_array['y_title']='fraction of galaxies'
                if c==0: 
                    edgecol='r'
                    myalpha=0.8
                else:
                    edgecol='k'#'#003964'
                    myalpha=1.0
                    #mydata=10**mydata
                    #mydata[:,0] = np.log10(mydata[:,0] + 0.24) - 0.1739  #Charbier to Kroupa Violeta
                    #mydata[:,0]= np.log10(1.2 * mydata[:,0])             #Charbier to Kroupa Croton+16
            
                myweights=np.ones_like(mydata[:,0])/len(mydata['name'+str(c)+'_data'][mydata['col_name_x']])
                print myweights, 'myalpha:', myalpha
               
                n, bins, patches=myax.hist(mydata[:,0], 
                                              int(myConfig.plot_config_array['contour_histo_nbins']),
                                              range=[10.6,12.2],
                                              edgecolor=edgecol, 
                                              facecolor=mycolor_code, 
                                              alpha=float(myConfig.plot_config_array['alpha']), 
                                              lw=2.0,
                                              histtype='step',
                                              weights=myweights,
                                              normed=False,
                                              stacked=False,
                                              ls=mylinestyle)                                                 

                
                print 'n:', n
                print 'bins:', bins
                print 'sum:', np.sum(n)
#                    
#                    axis0, = cb_axis.plot(bins[0:20], 
#                                      n/len(mydata[:,0]), 
#                                      color=mycolor_code, 
#                                      ls=mylinestyle, 
#                                      lw=float(myConfig.plot_config_array['mylw'])) 
                axis0, = myax.plot([], color=edgecol, linewidth=float(myConfig.plot_config_array['myframe_lw']), alpha=1, ls=mylinestyle)

                return key0, axis0

            def ScatterPlot(i, j):
            #SCATTER-PLOT SCATTER-PLOT SCATTER-PLOT SCATTER-PLOT SCATTER-PLOT SCATTER-PLOT SCATTER-PLOT SCATTER-PLOT SCATTER-PLOT SCATTER-PLOT SCATTER-PLOT SCATTER-PLOT SCATTER-PLOT SCATTER-PLOT                     
               
                mymarker_code, mycolor_code, mylinestyle, mycolor_map, mymarkerfacecolor, mydashes, colorbar, myalpha, myplotlegend = get_styles(i,j)
                N=int(myConfig.plot_config_array['contour_histo_nbins'])                
                if len(mydata['name'+str(j)+'_data'][mydata['col_name_x']])<5e6:
                    N=70
                    
                print 'SCATTER-PLOT!'               
                print 'j:', j, 'i:', i, 'marker:', mymarker_code, 'mylw:', float(myConfig.plot_config_array['mylw'])+float(myConfig.plot_config_array['lw_offset']), \
                      'fcol:', mymarkerfacecolor, 'myls:', mylinestyle, 'size:', myConfig.plot_config_array['markersize_subplot'], \
                      'errorbars:', myConfig.plot_config_array['subplot_errorbars'], 'filled:', myConfig.plot_config_array['filled_between_subplot'], 'nbins:', N               
                #key0=myax
                print 'max X:', np.max(mydata['name'+str(j)+'_data'][mydata['col_name_x']]), 'min X:', np.min(mydata['name'+str(j)+'_data'][mydata['col_name_x']])
                print 'max Y:', np.max(mydata['name'+str(j)+'_data'][mydata['col_name_y']]), 'min Y:', np.min(mydata['name'+str(j)+'_data'][mydata['col_name_y']])
                try:
                    print 'max Z:', np.max(mydata['name'+str(j)+'_data'][mydata['col_name_z']]), 'min Z:', np.min(mydata['name'+str(j)+'_data'][mydata['col_name_z']])
                except:
                    pass

                if plot_dim=='2d':

                    print mydata['name'+str(j)+'_data']
                    
                    key0.scatter(mydata['name'+str(j)+'_data'][mydata['col_name_x']],
                                mydata['name'+str(j)+'_data'][mydata['col_name_y']],
                                s=50,
                                alpha=1,
                                marker='x',
                                color='k')
                                
                    axis0, = key0.plot([],[], color='w', linewidth=float(myConfig.plot_config_array['myframe_lw']), alpha=1, ls=mylinestyle)
                else:
                    myax = fig.add_subplot(111, projection='3d')                   
                    myax.scatter(mydata['name'+str(j)+'_data'][mydata['col_name_x']],
                                mydata['name'+str(j)+'_data'][mydata['col_name_y']],
                                mydata['name'+str(j)+'_data'][mydata['col_name_z']],
                                s=5,
                                alpha=1,
                                marker='o',
                                zdir='z',
                                color='g')
                            
                    axis0, = key0.plot([],[],[], color='w', linewidth=float(myConfig.plot_config_array['myframe_lw']), alpha=1, ls=mylinestyle)
                    
                return key0, axis0, None, myplotlegend
                
            def HexBinPlot(i,j):              
                print 'HEXBIN-PLOT!'
                print 'j:', j, 'i:', i, 'marker:', mymarker_code, 'mylw:', float(myConfig.plot_config_array['mylw'])+float(myConfig.plot_config_array['lw_offset']), \
                      'fcol:', mymarkerfacecolor, 'myls:', mylinestyle, 'size:', myConfig.plot_config_array['markersize_subplot'], \
                      'errorbars:', myConfig.plot_config_array['subplot_errorbars'], 'filled:', myConfig.plot_config_array['filled_between_subplot'], \
                      'nbins:', myConfig.plot_config_array['bins']
                      
                key0 = myax
                hb = key0.hexbin(mydata['name'+str(c)+'_data'][mydata['col_name_x']],
                                 mydata['name'+str(c)+'_data'][mydata['col_name_y']],
                                 vmin=float(myConfig.plot_config_array['cb_min']),
                                 vmax=float(myConfig.plot_config_array['cb_max']),
                                 bins=myConfig.plot_config_array['bins'],
                                 mincnt=float(myConfig.plot_config_array['min_count']),
                                 gridsize=int(myConfig.plot_config_array['gridsize']),
                                 cmap='binary')

                if colorbar=='yes':
                    plot_colorbar(hb, j, int(myConfig.plot_config_array['cb_min']), int(myConfig.plot_config_array['cb_max']), int(myConfig.plot_config_array['cb_steps']))
                    
                return key0, None, None, myplotlegend            

            def Default():
            #START DEFAULT LOOP START DEFAULT LOOP START DEFAULT LOOP START DEFAULT LOOP START DEFAULT LOOP START DEFAULT LOOP
                
                mymarker_code, mycolor_code, mylinestyle, mycolor_map, mymarkerfacecolor, mydashes, colorbar, myalpha, myplotlegend = get_styles(i,i)             

                print 'MAIN plot: i:', i, 'myls:', mylinestyle, 'marker:', mymarker_code, 'mylw:', float(myConfig.plot_config_array['mylw']), \
                      'mycolor_code:', mycolor_code, 'errorbars:', myConfig.plot_config_array['error_bars_y'], \
                      'filled_between:', myConfig.plot_config_array['filled_between']

 
                #myax.set_zorder(20)
                key0 = myax
                #print mydata[0:5,[0,1,4,5]]
                key0.set_zorder(100)

                if myPlottypeTable[str(i)]=='default':                      
                    axis0, = key0.plot(mydata[:,0], 
                                      mydata[:,1],
                                      color=mycolor_code, 
                                      ls=mylinestyle, 
                                      lw=float(myConfig.plot_config_array['mylw']))
                    if mymarker_code!=' ' and mymarker_code!='no':
                        axis0 = add_marker(axis0, 
                                           key0,
                                           mydata[:,0], 
                                           mydata[:,1], 
                                           mymarker_code, 
                                           mycolor_code, 
                                           mylinestyle, 
                                           mymarkerfacecolor, 
                                           float(myConfig.plot_config_array['mylw']))                  
               
                    try:
                        add_errorbars(key0,
                                      myConfig.plot_config_array['errorbars'],                          
                                      myConfig.plot_config_array['error_bars_x'], 
                                      myConfig.plot_config_array['error_bars_y'], 
                                      myConfig.plot_config_array['filled_between'],
                                      mydata[:,0],
                                      mydata[:,1],
                                      mydata[:,2],                                  
                                      mydata[:,4],
                                      xerr_data2=mydata[:,3],
                                      yerr_data2=mydata[:,5],
                                      color=mycolor_code,
                                      fcolor=mycolor_code,
                                      alpha=float(myConfig.plot_config_array['alpha']))
                    except:
                        print 'error bars could not been set!'

                    add_tick_params(key0, 'both', 'minor', length=minor_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), set_visible=False)
                    add_tick_params(key0, 'both', 'major', length=major_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), set_visible=False) 
    
                    add_tick_params(myax, 'both', 'minor', length=minor_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), set_visible=True)
                    add_tick_params(myax, 'both', 'major', length=major_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), set_visible=True)

                    if myConfig.plot_config_array['share_axis_x']=='yes' or myConfig.plot_config_array['share_axis_y']=='yes':
                        
                        if myConfig.plot_config_array['share_axis_x']=='yes':
                            print 'TODO: share x'
                            
                        else:
                            key_share_y0=myax.twinx()
                            
                            axis_share_y0, = key_share_y0.plot(mydata[:,0], 
                                                              mydata[:,int(myConfig.plot_config_array['share_which_col_y'])],
                                                              color=mycolor_code, 
                                                              ls=mylinestyle, 
                                                              lw=float(myConfig.plot_config_array['mylw'])-float(myConfig.plot_config_array['mylw_offset']),
                                                              alpha=float(myConfig.plot_config_array['alpha'])-float(myConfig.plot_config_array['alpha_offset']))
                            
                            key_share_y0.set_ylabel(myConfig.plot_config_array['title_share_y'],
                                                    color='k',
                                                    fontsize=int(myConfig.plot_config_array['axis_label_fontsize']),
                                                    labelpad=50,
                                                    rotation=270)

                            if mymarker_code!=' ' and mymarker_code!='no':
                                axis_share_y0 = add_marker(axis_share_y0, 
                                                           key_share_y0,
                                                           mydata[:,0], 
                                                           mydata[:,int(myConfig.plot_config_array['share_which_col_y'])], 
                                                           mymarker_code, 
                                                           mycolor_code, 
                                                           mylinestyle, 
                                                           mymarkerfacecolor, 
                                                           float(myConfig.plot_config_array['mylw'])-float(myConfig.plot_config_array['mylw_offset']),
                                                           size=float(myConfig.plot_config_array['markersize'])-float(myConfig.plot_config_array['markersize_offset']),
                                                           alpha=1.0-float(myConfig.plot_config_array['alpha_offset']))                  
                       
                            try:
                                print 'add errorbars!'
                                add_errorbars(key_share_y0,
                                              myConfig.plot_config_array['errorbars'],                          
                                              myConfig.plot_config_array['error_bars_x'], 
                                              myConfig.plot_config_array['error_bars_y'], 
                                              myConfig.plot_config_array['filled_between'],
                                              mydata[:,0],
                                              mydata[:,int(myConfig.plot_config_array['share_which_col_y'])],
                                              mydata[:,2],                                  
                                              mydata[:,7],
                                              xerr_data2=mydata[:,3],
                                              yerr_data2=mydata[:,8],
                                              color=mycolor_code,
                                              fcolor=mycolor_code,
                                              alpha=float(myConfig.plot_config_array['alpha']))
                            except:
                                print 'error bars could not been set!'
                            
                            print 'ymin/ymax_share:', ymin_share, '/', ymax_share 
                            
                            add_tick_params(key_share_y0, 
                                            'y', 
                                            'minor', 
                                            isshare_y=True, 
                                            ymin_share=ymin_share, 
                                            ymax_share=ymax_share,
                                            top='off',
                                            bottom='off',                                            
                                            left='off', 
                                            right='on',
                                            length=minor_ticks_lenght, 
                                            width=float(myConfig.plot_config_array['myframe_lw']), 
                                            set_visible=True, 
                                            label_size_offset=0, 
                                            pad=6,
                                            color='k')
                            add_tick_params(key_share_y0, 
                                            'y', 
                                            'major',                                           
                                            isshare_y=True, 
                                            ymin_share=ymin_share, 
                                            ymax_share=ymax_share,
                                            top='off',
                                            bottom='off', 
                                            left='off', 
                                            right='on',
                                            length=major_ticks_lenght, 
                                            width=float(myConfig.plot_config_array['myframe_lw']), 
                                            set_visible=True, 
                                            label_size_offset=0, 
                                            pad=6,
                                            color='k')
                            
                            if myConfig.plot_config_array['minor_ticks_y_share_space']!='no':
                                key_share_y0.tick_params(axis='y',
                                which='minor', 
                                left='off', 
                                right=myConfig.plot_config_array['minor_ticks_y_set_right'], 
                                pad=10, 
                                labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), 
                                length=minor_ticks_lenght, 
                                width=float(myConfig.plot_config_array['myframe_lw']), 
                                direction='in', zorder=20) 
                                minor_y = MultipleLocator(float(myConfig.plot_config_array['minor_ticks_y_share_space']))
                                key_share_y0.yaxis.set_minor_locator(minor_y)
                                key_share_y0.set_zorder(50)
      
                else:
                    print 'BAR PLOT!'

                    myConfig.plot_config_array['y_title']='fraction of galaxies'

                    if mylinestyle=='':
                        mylinewidth=0
                        mylinestyle='-'
                    else:
                        mylinewidth=3
                        mymarkerfacecolor='None'
                 
                    axis0=key0.bar(mydata[:,0], 
                                   mydata[:,1]/sum(mydata[:,1]),
                                   (mydata[1,0]-mydata[0,0]),
                                   bottom=0,
                                   align='center',
                                   edgecolor=mycolor_code, 
                                   facecolor=mymarkerfacecolor, 
                                   alpha=myalpha, 
                                   linewidth=mylinewidth,
                                   ls=mylinestyle,
                                   hatch=mymarker_code)                                             

                myax.set_zorder(50)
                            
                return key0, axis0, None, myplotlegend           
            #END DEFAULT LOOP END DEFAULT LOOP END DEFAULT LOOP END DEFAULT LOOP END DEFAULT LOOP END DEFAULT LOOP END DEFAULT LOOP END DEFAULT LOOP END DEFAULT LOOP END
 
########################################################################################################################################################################               
#START START START START START START START START START START START START START START START START START START START START START START START START START START START START

            #arrange additional or change keywords
            add_plot_pars=[add_plot_par1,add_plot_par2,add_plot_par3,add_plot_par4,add_plot_par5,add_plot_par6,add_plot_par7,add_plot_par8,add_plot_par9,add_plot_par10,add_plot_par11,add_plot_par12,add_plot_par13,add_plot_par14]
            mykeyword_changes=[mykeyword_change1,mykeyword_change2,mykeyword_change3,mykeyword_change4,mykeyword_change5,mykeyword_change6,mykeyword_change7,mykeyword_change8,mykeyword_change9,mykeyword_change10,mykeyword_change11,mykeyword_change12,mykeyword_change13,mykeyword_change14]                    

            i=0                
            for entries in mykeyword_changes:
                #print entries
                if entries[0]!='keyword':
                    myConfig.plot_config_array[entries[0]] = entries[1]
                i+=1

            i=0                
            for entries in add_plot_pars:
                #print entries
                if entries[0]!='default':
                    myConfig.plot_config_array[entries[0]] = entries[1]  
                i+=1
            i=1
            while i<len(add_plot_pars)+1:
                if myConfig.plot_config_array['keyword_change'+str(i)]!='default':
                    changed_keyword=mL.multicolTestAlgorithm(myConfig.plot_config_array['keyword_change'+str(i)])  
                    myConfig.plot_config_array[changed_keyword[0]] = changed_keyword[1]
                i+=1
            
    
            plt.rcParams['axes.linewidth'] = 2
    
            if myConfig.plot_config_array['2D-plot']=='yes':
                plot_dim='2d'
            else:
                plot_dim='3d'
           
            if print_redshift==False: 
                print_redshift=myConfig.plot_config_array['print_redshift']                
                try:
                    float(print_redshift)
                    print_redshift='z='+str(print_redshift)
                except:
                    pass    


#PLOT MAIN LOOP PLOT MAIN LOOP PLOT MAIN LOOP PLOT MAIN LOOP PLOT MAIN LOOP PLOT MAIN LOOP PLOT MAIN LOOP PLOT MAIN LOOP PLOT MAIN LOOP PLOT MAIN LOOP PLOT MAIN LOOP PLOT MAIN LOOP                    
            i=0
            c=0           
            count_add_axis=0
            #print 'catname:', catname, 'plot_key:', plot_key, 'x/y_title', myConfig.plot_config_array['x_title'], '/', myConfig.plot_config_array['y_title'], \
                  #'use_cat_colors:', myConfig.plot_config_array['use_cat_colors']       


########################################################################################################################################################################
#CHOOSE A PLOT TYPE OR GO TO DEFAULT Main() 

            for f in [prefix+'x_range_min', prefix+'x_range_max', prefix+'y_range_min', prefix+'y_range_max', prefix+'share_x_range_min', prefix+'share_x_range_max', prefix+'share_y_range_min', prefix+'share_y_range_max']:
                if f.find('min') and myConfig.plot_config_array[f]=='min':                    
                    myConfig.plot_config_array[f]=min(mydata[:,i*data_block_offset])
                elif f.find('max')  and myConfig.plot_config_array[f]=='max':
                    myConfig.plot_config_array[f]=max(mydata[:,i*data_block_offset])
                else:
                    myConfig.plot_config_array[f]=myConfig.plot_config_array[f].astype(np.float)
                        #print 'f:', f, myConfig.plot_config_array[f]

            xmin = myConfig.plot_config_array[prefix+'x_range_min']
            xmax = myConfig.plot_config_array[prefix+'x_range_max']

            ymin = myConfig.plot_config_array[prefix+'y_range_min']
            ymax = myConfig.plot_config_array[prefix+'y_range_max'] 
            
            xmin_share = myConfig.plot_config_array[prefix+'share_x_range_min']
            xmax_share = myConfig.plot_config_array[prefix+'share_x_range_max']

            ymin_share = myConfig.plot_config_array[prefix+'share_y_range_min']
            ymax_share = myConfig.plot_config_array[prefix+'share_y_range_max']            

            #print myPlottypeTable[str(i)]            
            if myPlottypeTable[str(i)]=='default' or myPlottypeTable[str(i)]=='barplot':

                key0, axis0, mylegend, myplotlegend = Default()
                
                try:                
                    mylegend=myConfig.plot_config_array['plot_legend'+str(i)]
                except:
                    mylegend=mydata['name'+str(i)+'_legend']

                if i==0:
                    axis = [axis0]
                    legend= [mylegend]  

            else:

                
                if myConfig.plot_config_array['x_title']!='False': 
                    myConfig.plot_config_array['x_title']=mydata['name_x']
                if myConfig.plot_config_array['y_title']!='False':                     
                    myConfig.plot_config_array['y_title']=mydata['name_y']
                if myConfig.plot_config_array['z_title']!='False':                     
                    try:
                        myConfig.plot_config_array['z_title']=mydata['name_z']
                    except:
                        print 'no "z_title" found!'
                        myConfig.plot_config_array['z_title']=''

                #print 'mydata[name_x/name_y]', mydata['name_x'], '/', mydata['name_y']

                while c<loops:

                    mymarker_code, mycolor_code, mylinestyle, mycolor_map, mymarkerfacecolor, mydashes, colorbar, myalpha, myplotlegend = get_styles(i,c) 
                                                                    
                    catname = mydata['catname'+str(c)]
                    try:                
                        mylegend=myConfig.plot_config_array['plot_legend'+str(c)] #mydata['name'+str(c)+'_legend'] 
                    except:
                        mylegend=mydata['name'+str(c)+'_legend']          
                    
                    print 'i:', i, 'c:', c, 'catname:', catname, 'mylegend:', mylegend, 'myls:', mylinestyle, 'mycolor_code:', mycolor_code, 'color_map:', mycolor_map, 
                    print 'col_name_x:', mydata['col_name_x'], 'col_name_y:', mydata['col_name_y']
                    
                    if myConfig.plot_config_array['contour_log']=='yes':
    
                        if myConfig.plot_config_array['contour_log_wich_axis']=='x':
                            mydata['name'+str(c)+'_data'][mydata['col_name_x']]=np.log10(mydata['name'+str(c)+'_data'][mydata['col_name_x']])
                        elif myConfig.plot_config_array['contour_log_wich_axis']=='y':
                            mydata['name'+str(c)+'_data'][mydata['col_name_y']]=np.log10(mydata['name'+str(c)+'_data'][mydata['col_name_y']]) 
                        elif myConfig.plot_config_array['contour_log_wich_axis']=='both':     
                            mydata['name'+str(c)+'_data'][mydata['col_name_x']]=np.log10(mydata['name'+str(c)+'_data'][mydata['col_name_x']])
                            mydata['name'+str(c)+'_data'][mydata['col_name_y']]=np.log10(mydata['name'+str(c)+'_data'][mydata['col_name_y']]) 

                    def caseSwitcher(plot_type):
                    
                        choose = {
                            'contour': ContourPlot,
                            'scatter': ScatterPlot,
                            'barhisto': BarHistogram,
                            'hexbins': HexBinPlot,
                            'surface3D': Surface3DPlot
                            }
                            
                        func = choose.get(plot_type)
                        return func(i,c)
                               
                    key0, axis0, mylegend, myplotlegend = caseSwitcher(myPlottypeTable[str(c)]) 
                    
                    if myConfig.plot_config_array['x_title']!='False':
                        fig.text(0.5, 0.03, myConfig.plot_config_array['x_title'], ha='center', fontsize=float(myConfig.plot_config_array['axis_label_fontsize']))
                    if myConfig.plot_config_array['y_title']!='False':
                        fig.text(0.01, 0.5, myConfig.plot_config_array['y_title'], va='center', rotation='vertical', fontsize=float(myConfig.plot_config_array['axis_label_fontsize']))                           
                    
                    xmin = myConfig.plot_config_array[prefix+'x_range_min']
                    xmax = myConfig.plot_config_array[prefix+'x_range_max']
        
                    ymin = myConfig.plot_config_array[prefix+'y_range_min']
                    ymax = myConfig.plot_config_array[prefix+'y_range_max']                                

                    xmin_share = myConfig.plot_config_array[prefix+'share_x_range_min']
                    xmax_share = myConfig.plot_config_array[prefix+'share_x_range_max']

                    ymin_share = myConfig.plot_config_array[prefix+'share_y_range_min']
                    ymax_share = myConfig.plot_config_array[prefix+'share_y_range_max'] 

                    #add_tick_params(key0, 'both', 'minor', length=minor_ticks_lenght, width=minor_ticks_width, set_visible=False)
                    #add_tick_params(key0, 'both', 'major', length=major_ticks_lenght, width=major_ticks_width, set_visible=False)    
       
                    if c==0 and myplotlegend!='no':                                   
                        axis = [axis0]
                        legend= [mylegend]
                    elif c==0 and myplotlegend=='no':
                        axis=[]
                        legend=[]
                      
                    c+=1                       

                add_tick_params(myax, 'both', 'minor', color='r', length=minor_ticks_lenght, width=minor_ticks_width, set_visible=True)
                add_tick_params(myax, 'both', 'major', color='r', length=major_ticks_lenght, width=major_ticks_width, set_visible=True)
                
                if myConfig.plot_config_array['add_axis']=='yes':
                    add_tick_params(myax, 'both', 'major', length=0, width=major_ticks_width, set_visible=False)
                else:
                    add_tick_params(myax, 'both', 'major', length=major_ticks_lenght, width=major_ticks_width, set_visible=False)
            
                loops=0
                i=c

########################################################################################################################################################################
#ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS ADD-AXIS                
            p=0
            count_add_axis_x=0
            count_add_axis_y=0
            while p<loops:
                if add_xaxisTable[str(p)]=='yes':
                    add_axis('add_x', 
                             mydata[:,data_block_offset*p : data_block_offset*p + data_block_offset], 
                             'key_add'+str(p), 
                             key0, 
                             p, 
                             count_add_axis_x, 
                             int(myConfig.plot_config_array['add_which_col_x']), 
                             int(myConfig.plot_config_array['add_axis_steps_x']))   
                    count_add_axis_x+=1
                if add_yaxisTable[str(p)]=='yes':
                    add_axis('add_y', mydata[:,data_block_offset*p : data_block_offset*p + data_block_offset], 'key_add'+str(p), key0, p, count_add_axis_y, int(myConfig.plot_config_array['add_which_col_y']), int(myConfig.plot_config_array['add_axis_steps_y']))   
                    count_add_axis_y+=1
                p+=1     

            axis_sep=[]
            legend_sep=[]

            i=1    
            if (myPlottypeTable['0']=='default' or myPlottypeTable['0']=='barplot') and myPlotLegendTable['0']!='no':
                axis = [axis0]
                legend = [myConfig.plot_config_array['plot_legend0']]
            elif myplotlegend=='no':
                axis=[]
                legend=[]
            else:
                i=c

########################################################################################################################################################################
#SUBPLOT LOOP
            while i<loops:
                mymarker_code, mycolor_code, mylinestyle, mycolor_map, mymarkerfacecolor, mydashes, colorbar, myalpha, myplotlegend = get_styles(i,i) 
                print 'LOOP plot: i:', i, 'c:', c, 'myls:', mylinestyle, 'marker:', mymarker_code, 'mylw:', float(myConfig.plot_config_array['mylw']), \
                       'mycolor_code:', mycolor_code, 'fill between?', myConfig.plot_config_array['filled_between']           

                #print mydata[:,[i*data_block_offset,i*data_block_offset+1]]
#                myConfig.plot_config_array['error_bars_y']='no'
#                print myConfig.plot_config_array['plot_legend'+str(i)]
#                if myConfig.plot_config_array['plot_legend'+str(i)].find('Shan+17')!=-1:# and add_obs!=False:
#                    myConfig.plot_config_array['error_bars_y']='yes'
#                    myConfig.plot_config_array['errorbars']='False'
#                    myConfig.plot_config_array['filled_between'] = 'True'
#                elif myConfig.plot_config_array['plot_legend'+str(i)].find('g-i')!=-1:
#                    myConfig.plot_config_array['errorbars']='True'
#                    myConfig.plot_config_array['filled_between'] = 'False'                    
                    
                key1 = myax.twinx()
                key2 = 'axis' + str(i)           
                key1.set_zorder(20)
                if i==int(myConfig.plot_config_array['set_thin_line']):
                    myConfig.plot_config_array['mylw']=float(myConfig.plot_config_array['mylw'])-int(myConfig.plot_config_array['thin_line'])
                elif i==int(myConfig.plot_config_array['set_thin_line2']):
                    myConfig.plot_config_array['mylw']=float(myConfig.plot_config_array['mylw'])-int(myConfig.plot_config_array['thin_line2'])
                
                if myPlottypeTable[str(i)]=='default':
                    key2, = key1.plot(mydata[:,i*data_block_offset], 
                                      mydata[:,i*data_block_offset+1],
                                      color=mycolor_code, 
                                      ls=mylinestyle, 
                                      lw=float(myConfig.plot_config_array['mylw']))
                    if mymarker_code!=' ' and mymarker_code!='no':
                        key2 = add_marker(key2, 
                                          key1, 
                                          mydata[:,i*data_block_offset],
                                          mydata[:,i*data_block_offset+1], 
                                          mymarker_code, 
                                          mycolor_code, 
                                          mylinestyle, 
                                          mymarkerfacecolor, 
                                          float(myConfig.plot_config_array['mylw']))
                                  
                    #Errorbars
                    try:
                        if i==3 or i==5 or i==7 or i==8:
                            add_errorbars(key1,
                                          myConfig.plot_config_array['errorbars'],                              
                                          myConfig.plot_config_array['error_bars_x'], 
                                          myConfig.plot_config_array['error_bars_y'], 
                                          myConfig.plot_config_array['filled_between'],
                                          mydata[:,i*data_block_offset],
                                          mydata[:,i*data_block_offset+1],
                                          mydata[:,i*data_block_offset+2],                             
                                          mydata[:,i*data_block_offset+4],
                                          xerr_data2=mydata[:,i*data_block_offset+3],
                                          yerr_data2=mydata[:,i*data_block_offset+5],
                                          color=mycolor_code,
                                          fcolor=mycolor_code,
                                          alpha=float(myConfig.plot_config_array['alpha']),
                                          ls=mylinestyle)
                    except:
                        pass


                    if myConfig.plot_config_array['share_axis_x']=='yes' or myConfig.plot_config_array['share_axis_y']=='yes':
                        
                        if myConfig.plot_config_array['share_axis_x']=='yes':
                            print 'TODO: share x'
                            
                        if myConfig.plot_config_array['share_axis_y']=='yes':
                            key_share_y=myax.twinx()
                            print mydata[:,i*data_block_offset][0:10,], mydata[:,i*data_block_offset+int(myConfig.plot_config_array['share_which_col_y'])][0:10]
                            axis_share_y, = key_share_y.plot(mydata[:,i*data_block_offset], 
                                                              mydata[:,i*data_block_offset+int(myConfig.plot_config_array['share_which_col_y'])],
                                                              color=mycolor_code, 
                                                              ls=mylinestyle, 
                                                              lw=float(myConfig.plot_config_array['mylw'])-float(myConfig.plot_config_array['mylw_offset']),
                                                              alpha=float(myConfig.plot_config_array['alpha'])-float(myConfig.plot_config_array['alpha_offset']))

                            if mymarker_code!=' ' and mymarker_code!='no':
                                axis_share_y = add_marker(axis_share_y, 
                                                           key_share_y,
                                                           mydata[:,i*data_block_offset], 
                                                           mydata[:,i*data_block_offset+int(myConfig.plot_config_array['share_which_col_y'])], 
                                                           mymarker_code, 
                                                           mycolor_code, 
                                                           mylinestyle, 
                                                           mymarkerfacecolor, 
                                                           float(myConfig.plot_config_array['mylw'])-float(myConfig.plot_config_array['mylw_offset']),
                                                           size=float(myConfig.plot_config_array['markersize'])-float(myConfig.plot_config_array['markersize_offset']),
                                                           alpha=1.0-float(myConfig.plot_config_array['alpha_offset']))                  
                       
                            try:
                                print 'add errorbars!'
                                add_errorbars(key_share_y,
                                              myConfig.plot_config_array['errorbars'],                          
                                              myConfig.plot_config_array['error_bars_x'], 
                                              myConfig.plot_config_array['error_bars_y'], 
                                              myConfig.plot_config_array['filled_between'],
                                              mydata[:,i*data_block_offset],
                                              mydata[:,i*data_block_offset+int(myConfig.plot_config_array['share_which_col_y'])],
                                              mydata[:,i*data_block_offset+2],                                  
                                              mydata[:,i*data_block_offset+7],
                                              xerr_data2=mydata[:,i*data_block_offset+3],
                                              yerr_data2=mydata[:,i*data_block_offset+8],
                                              color=mycolor_code,
                                              fcolor=mycolor_code,
                                              alpha=float(myConfig.plot_config_array['alpha']))
                            except:
                                print 'error bars could not been set!'

                            add_tick_params(key_share_y, 'both', 'minor', color='r', length=0, width=major_ticks_width, set_visible=False)                                
                            add_tick_params(key_share_y, 'both', 'major', color='r', length=0, width=major_ticks_width, set_visible=False)                                
                                
                else:
                    print 'BAR PLOT!'

                    myConfig.plot_config_array['y_title']='fraction of galaxies'
                    
                    if mylinestyle=='':
                        mylinewidth=0
                        mylinestyle='-'
                    else:
                        mylinewidth=3
                        mymarkerfacecolor='None'
                        
                    key2=key1.bar(mydata[:,i*data_block_offset], 
                                   mydata[:,i*data_block_offset+1]/sum(mydata[:,i*data_block_offset+1]),
                                   (mydata[1,i*data_block_offset]-mydata[0,i*data_block_offset]),
                                   bottom=0,
                                   align='edge',
                                   edgecolor=mycolor_code, 
                                   facecolor=mymarkerfacecolor,
                                   alpha=myalpha, 
                                   linewidth=mylinewidth,
                                   ls=mylinestyle,
                                   hatch=mymarker_code)   



                if myConfig.plot_config_array['share_axis_y']=='yes':

                    add_tick_params(key1, 'x', 'minor', length=minor_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), set_visible=True)
                    add_tick_params(key1, 'x', 'major', length=major_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), set_visible=True)
                                        
                    add_tick_params(key1, 'y', 'minor', right='off', length=minor_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), set_visible=True)
                    add_tick_params(key1, 'y', 'major', right='off', length=major_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), set_visible=True)

                elif myConfig.plot_config_array['add_axis']=='yes':
                    add_tick_params(key1, 'both', 'major', length=0, width=major_ticks_width, set_visible=False)
                else:
                    add_tick_params(key1, 'both', 'minor', length=minor_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), set_visible=False)
                    add_tick_params(key1, 'both', 'major', length=major_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), set_visible=False)                   
                
                if i==int(myConfig.plot_config_array['set_thin_line']):
                    myConfig.plot_config_array['mylw']=float(myConfig.plot_config_array['mylw'])+int(myConfig.plot_config_array['thin_line'])
                elif i==int(myConfig.plot_config_array['set_thin_line2']):
                    myConfig.plot_config_array['mylw']=float(myConfig.plot_config_array['mylw'])+int(myConfig.plot_config_array['thin_line2'])
    
                if myplotlegend!='no':
                    axis.extend([key2])
                    legend.extend([myConfig.plot_config_array['plot_legend'+str(i)]])
                i+=1                           

                #CMASS density+mass SMF cut
                #plt.plot([np.log10(1.15741161e+11),np.log10(1.15741161e+11)],[np.log10(2.49150233e-04),-10], color='#225588',marker='', ls='-', lw=4.0, zorder=zorder_cut)  
                #plt.plot([np.log10(1.80972295e+11),np.log10(1.80972295e+11)],[np.log10(1.82572068e-04),-10], color='k',marker='', ls=':', lw=4.0, zorder=zorder_cut) 

########################################################################################################################################################################
#ADD OBSERVATIONS
            if my_add_subplot!=[] and nr_added_subplots!=0 and add_obs!=False:
                p=0
                for f in add_obs:

                    def caseSwitcher(plot_type):
                    
                        choose = {
                            'contour': ContourPlot,
                            'barhisto': BarHistogram,
                            'scatter': ScatterPlot,
                            'hexbins': HexBinPlot,
                            'default': Default_sub,
                            'barplot': Default_sub
                            }
                            
                        func = choose.get(plot_type)
                        if myPlottypeTable[str(0)]=='default':
                            return func(int(i)+int(p),int(p))
                        else:
                            return func(int(i)+int(p),int(i)+int(p))                            
                             
                    def Default_sub(i,j):
                        #print 'HERE! Default_sub j:', j, 'i:', i, my_add_subplot[str(j)][0]
                        return add_subplot(my_add_subplot[str(j)],i,j)
                       
                    #print 'i:', i, 'c:', c, 'p:', p, 'f:', f#, 'data:', my_add_subplot[str(f)][0]
                    key_sub1, key_sub2, mylegend, myplotlegend = caseSwitcher(myPlottypeTable[str(int(i)+int(p))])
                    print 'legend:', mylegend, 'myplotlegend:', myplotlegend#'axis:', key_sub2
#                    if mylegend.find('Big')!=-1:
#                        key_sub2, = myax.plot([], color=mycolor_code, linewidth=8.0, alpha=0.3, ls='-')
                    if mylegend!=None and myplotlegend!='no':
                        if myplotlegend=='sep':
                            axis_sep.extend([key_sub2])
                            legend_sep.extend([mylegend])
                        else:
                            axis.extend([key_sub2])
                            legend.extend([mylegend])
                    p+=1
    

########################################################################################################################################################################
#FURTHER SPECIFICATIONS

            #fig.text(0.15, 0.35, r'${\mathrm{\pi_{max}}=40}$ [Mpc]', fontsize=float(myConfig.plot_config_array['text_fontsize'])-5)             
            #fig.text(0.6, 0.35, r'${\mathrm{\pi_{max}}=150}$ [Mpc]', fontsize=float(myConfig.plot_config_array['text_fontsize'])-5)
            #fig.text(0.75, 0.92, r'${n \times 10^{-2}} {[h^3 Mpc^{-3}]}$', fontsize=float(myConfig.plot_config_array['text_fontsize'])-8)
            #fig.text(0.85, 0.86, '$n$=3.160', fontsize=float(myConfig.plot_config_array['text_fontsize'])-8)
            #fig.text(0.85, 0.79, '$n$=1.000', fontsize=float(myConfig.plot_config_array['text_fontsize'])-8)
            #fig.text(0.88, 0.72, '$n$=0.316', fontsize=float(myConfig.plot_config_array['text_fontsize'])-8)
            #plt.axvline(x=1.028, ymin=-100, ymax=100, color='k', ls='-', lw=2.0, zorder=30) 

            #ADD VERTICAL OR HORIZONTAL LINES:
            if myConfig.plot_config_array[prefix+'hline_ypos']!='False':
                my_hlines=mL.multicolTestAlgorithm(myConfig.plot_config_array[prefix+'hline_ypos'])
                a=0
                while a<len(my_hlines):
                    myax.axhline(y=float(my_hlines[a]), xmin=-100, xmax=100, color='k', ls='--', lw=2.0, zorder=100)
                    a+=1                

            if myConfig.plot_config_array[prefix+'vline_xpos']!='False':
                my_vlines=mL.multicolTestAlgorithm(myConfig.plot_config_array[prefix+'vline_xpos'])
                a=0
                while a<len(my_vlines):
                    myax.axvline(x=float(my_vlines[a]), ymin=-100, ymax=100, color='k', ls='--', lw=2.0, zorder=100) 
                    a+=1                  

            #SET TICKS                   
            yticks=[]
            if myConfig.plot_config_array['yticks']!='default':
                my_yticks=mL.multicolTestAlgorithm(myConfig.plot_config_array['yticks'])
                a=0
                while a<len(my_yticks):
                    yticks+=[float(my_yticks[a])+float(myConfig.plot_config_array['yticks_minor_offset'])]
                    a+=1

                myax.set_yticks(yticks)         
                myax.set_yticklabels(my_yticks)
                myax.tick_params(axis='y',
                                 which='minor', 
                                 bottom='off', 
                                 top='off', 
                                 pad=10, 
                                 labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), 
                                 direction='in', 
                                 width=float(myConfig.plot_config_array['myframe_lw']), 
                                 zorder=20)          
                
                if myConfig.plot_config_array['yticks_minor']=='default':
                     myax.tick_params(axis='y', 
                                      which='minor', 
                                      top='on', 
                                      bottom='on', 
                                      pad=10, 
                                      labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), 
                                      length=minor_ticks_lenght, 
                                      width=float(myConfig.plot_config_array['myframe_lw']),
                                      direction='in',
                                      zorder=20)                    
                elif float(myConfig.plot_config_array['yticks_minor'])!=0.0:
                    minor_y= MultipleLocator(float(myConfig.plot_config_array['yticks_minor']))
                    myax.yaxis.set_minor_locator(minor_y)   


            xticks=[]
            if myConfig.plot_config_array['xticks']!='default':
                my_xticks=mL.multicolTestAlgorithm(myConfig.plot_config_array['xticks'])
                a=0
                while a<len(my_xticks):
                    xticks+=[float(my_xticks[a])+float(myConfig.plot_config_array['xticks_minor_offset'])]
                    a+=1

                myax.set_xticks(xticks)         
                myax.set_xticklabels(my_xticks)
                myax.tick_params(axis='x',
                                 which='minor',
                                 bottom='off',
                                 top='off',
                                 pad=10,
                                 labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']),
                                 direction='in', width=float(myConfig.plot_config_array['myframe_lw']),
                                 zorder=20)          
                
                if myConfig.plot_config_array['xticks_minor']=='default':
                     myax.tick_params(axis='x',
                                      which='minor',
                                      top='on', 
                                      bottom='on', 
                                      pad=10, 
                                      labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), 
                                      length=minor_ticks_lenght, 
                                      width=float(myConfig.plot_config_array['myframe_lw']), 
                                      direction='in', 
                                      zorder=20)                    
                elif float(myConfig.plot_config_array['xticks_minor'])!=0.0:
                    minor_x = MultipleLocator(float(myConfig.plot_config_array['xticks_minor']))
                    myax.xaxis.set_minor_locator(minor_x)     
                
                if myConfig.plot_config_array['share_axis_y']!='yes':                    
                    myax.tick_params(axis='y', which='both', left='on', right='on', pad=10, labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), width=float(myConfig.plot_config_array['myframe_lw']), direction='in', zorder=20) 

            else:
                if myConfig.plot_config_array['add_axis']=='yes':                   
                    if any(x=='yes' for x in add_xaxisTable.values()) and any(x=='yes' for x in add_yaxisTable.values()):
                        
                        myax.tick_params(axis='x', which='minor', bottom='on', top='off', pad=10, labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), width=float(myConfig.plot_config_array['myframe_lw']), direction='in', zorder=20) 
                        myax.tick_params(axis='x', which='major', bottom='on', top='off', pad=10, labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), width=float(myConfig.plot_config_array['myframe_lw']), direction='in', zorder=20) 
                                               
                        myax.tick_params(axis='y', which='minor', left='on', right='off', pad=10, labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), width=float(myConfig.plot_config_array['myframe_lw']), direction='in', zorder=20)                    
                        myax.tick_params(axis='y', which='major', left='on', right='off', pad=10, labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), width=float(myConfig.plot_config_array['myframe_lw']), direction='in', zorder=20)                    

                    elif any(x=='yes' for x in add_xaxisTable.values()):
                        
                        myax.tick_params(axis='x', which='minor', bottom='on', top='off', pad=10, labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), width=float(myConfig.plot_config_array['myframe_lw']), direction='in', zorder=20) 
                        myax.tick_params(axis='x', which='major', bottom='on', top='off', pad=10, labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), width=float(myConfig.plot_config_array['myframe_lw']), direction='in', zorder=20) 

                        myax.tick_params(axis='y', which='both', left='on', right='on', pad=10, labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), width=float(myConfig.plot_config_array['myframe_lw']), direction='in', zorder=20) 
                                         
                    elif any(x=='yes' for x in add_yaxisTable.values()):
                        
                        myax.tick_params(axis='x', which='both', bottom='on', top='on', pad=10, labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), width=float(myConfig.plot_config_array['myframe_lw']), direction='in', zorder=20) 

                        myax.tick_params(axis='y', which='minor', left='on', right='off', pad=10, labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), width=float(myConfig.plot_config_array['myframe_lw']), direction='in', zorder=20)                    
                        myax.tick_params(axis='y', which='major', left='on', right='off', pad=10, labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), width=float(myConfig.plot_config_array['myframe_lw']), direction='in', zorder=20)                    
 
                else:
                    myax.tick_params(axis='x', which='both', bottom='on', top='on', pad=10, labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), width=float(myConfig.plot_config_array['myframe_lw']), direction='in', zorder=20) 
                    myax.tick_params(axis='y', which='both', left='on', right='on', pad=10, labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), width=float(myConfig.plot_config_array['myframe_lw']), direction='in', zorder=20) 

            #SET AXIS SPECIFICATIONS
            if myConfig.plot_config_array['minor_ticks_x_space']!='no':
                myax.tick_params(axis='x', 
                                 which='minor',
                                 bottom=myConfig.plot_config_array['minor_ticks_x_set_bottom'], 
                                 top=myConfig.plot_config_array['minor_ticks_x_set_top'],
                                 pad=10, 
                                 labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']),
                                 length=minor_ticks_lenght, 
                                 width=float(myConfig.plot_config_array['myframe_lw']), 
                                 direction='in', 
                                 zorder=20)          
                minor_x = MultipleLocator(float(myConfig.plot_config_array['minor_ticks_x_space']))
                myax.xaxis.set_minor_locator(minor_x)
                myax.set_zorder(50)
                
            if myConfig.plot_config_array['minor_ticks_y_space']!='no':
                myax.tick_params(axis='y',
                                 which='minor', 
                                 left=myConfig.plot_config_array['minor_ticks_y_set_left'], 
                                 right=myConfig.plot_config_array['minor_ticks_y_set_right'], 
                                 pad=10, 
                                 labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), 
                                 length=minor_ticks_lenght, 
                                 width=float(myConfig.plot_config_array['myframe_lw']), 
                                 direction='in', zorder=20) 
                minor_y = MultipleLocator(float(myConfig.plot_config_array['minor_ticks_y_space']))
                myax.yaxis.set_minor_locator(minor_y)
                myax.set_zorder(50)
                            

            
            #LEGEND          
            if myConfig.plot_config_array['use_loc_co']=='yes':
                myloc=[float(myConfig.plot_config_array['use_loc_co_x']), float(myConfig.plot_config_array['use_loc_co_y'])]
            else:
                myloc=myConfig.plot_config_array['legend_position']

   
            if myConfig.plot_config_array[prefix+'plot_legend']=='yes':
 
                if myConfig.plot_config_array['legend_fancy']=='True':
                    myfancybox=True
                    myshadow=True
                    myframeon=True
                else:
                    myfancybox=False
                    myshadow=False
                    myframeon=False               
                
                mylegend = fig.legend(axis,
                           legend,
                           loc=myloc,
                           ncol=int(myConfig.plot_config_array['legend_ncols']),
                           fontsize=float(myConfig.plot_config_array['legend_fontsize']),
                           shadow=myshadow, 
                           fancybox=myfancybox,
                           borderaxespad=1,
                           facecolor='w',
                           numpoints=1, 
                           frameon=myframeon,
                           labelspacing=0.3,
                           handlelength=1.6,
                           columnspacing=0.95,
                           handletextpad=0.5)#float(myConfig.plot_config_array['legend_handletextpad']))

                mylegend.set_zorder(30)

            if legend_sep!=[]:
                #print myConfig.plot_config_array['print_redshift_sep']
                if myConfig.plot_config_array['print_redshift_sep']!='':
                    fig.text(float(myConfig.plot_config_array['legend_sep_xpos']), float(myConfig.plot_config_array['legend_sep_ypos'])+0.01+float(myConfig.plot_config_array['legend_sep_offset'])*len(legend_sep), myConfig.plot_config_array['print_redshift_sep'], fontsize=float(myConfig.plot_config_array['legend_fontsize'])-2)      
                
                mylegend = fig.legend(axis_sep,
                           legend_sep,
                           loc=[float(myConfig.plot_config_array['legend_sep_xpos']),float(myConfig.plot_config_array['legend_sep_ypos'])],
                           fontsize=float(myConfig.plot_config_array['legend_fontsize'])-2,
                           shadow=False, 
                           fancybox=False,
                           borderaxespad=3,
                           facecolor='w',
                           numpoints=1, 
                           frameon=False,
                           labelspacing=0.2)                
                
                mylegend.set_zorder(30)
              

                   
########################################################################################################################################################################

        
#END END END END END END END END END END END END END END END END END END END END END END END END END END END END END END END END END END END END END END END END END END END END END 
########################################################################################################################################################################

        def plot_axis_ticks(axis):
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            #WARNING! THIS FUNCTION ONLY WORKS IF X- OR YTICKS ARE SET MANUALLY VIA!
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if myConfig.plot_config_array['no_first_last_xticks']=='True':  
                plt.setp(axis.get_xticklabels()[0], visible=False)    
                plt.setp(axis.get_xticklabels()[-1], visible=False)
            elif myConfig.plot_config_array['no_first_xticks']=='True':
                plt.setp(axis.get_xticklabels()[0], visible=False)
            elif myConfig.plot_config_array['no_last_xticks']=='True':
                plt.setp(axis.get_xticklabels()[-1], visible=False)
            elif myConfig.plot_config_array['no_xticks']=='True':
                 plt.setp(axis.get_xticklabels(), visible=False)
                  
            if myConfig.plot_config_array['no_first_last_yticks']=='True':  
                plt.setp(axis.get_yticklabels()[0], visible=False)    
                plt.setp(axis.get_yticklabels()[-1], visible=False)
            elif myConfig.plot_config_array['no_first_yticks']=='True':
                print 'HERE! 2062', axis.get_yticklabels()[0]
                plt.setp(axis.get_yticklabels()[0], visible=False)
            elif myConfig.plot_config_array['no_last_yticks']=='True': 
                plt.setp(axis.get_yticklabels()[-1], visible=False)
            elif myConfig.plot_config_array['no_yticks']=='True':
                 plt.setp(axis.get_yticklabels(), visible=False)

########################################################################################################################################################################
#MAIN PLOT

        myConfig = conf.Configuration()
        
        mpl.style.use('classic')            
        mpl.mathtext.use_cm = False
        mpl.rcParams['mathtext.fontset'] = 'custom'
        mpl.rcParams['mathtext.tt'] = 'Typewriter'
        mpl.mathtext.fallback_to_cm = True

        myConfig.plotMarkerColorKewords() 
        if custom_plot_key==None:
            myConfig.plotKeywords(myconfig_datafile=config_path+myconfig_datafile)
        else:
            print config_path+myconfig_datafile[0:len(myconfig_datafile)-11]+custom_plot_key+'_config.txt'
            myConfig.plotKeywords(myconfig_datafile=config_path+myconfig_datafile[0:len(myconfig_datafile)-11]+custom_plot_key+'_config.txt')  
            
        myLinestyleTable, myColourTable, myMarkerTable, myMarkerfacecolorTable, myColormapTable, myColorbarTable, myPlottypeTable, myAlphaTable, myPlotLegendTable, add_xaxisTable, add_yaxisTable = myConfig.load_matplot_stylefiles(myConfig.plot_config_array['plot_legend'+str(0)], ncolours=int(myConfig.plot_config_array['use_cb_colors']))
        
        nr_plots_x  = int(myConfig.plot_config_array['nr_plots_x'])   
        nr_plots_y  = int(myConfig.plot_config_array['nr_plots_y'])

        k=0
        l=0
        if nr_plots_x==1 and nr_plots_y==1:
            fig = plt.figure(figsize=(float(myConfig.plot_config_array['size_x']),float(myConfig.plot_config_array['size_y'])), dpi=300)
            
            if myPlottypeTable=='surface3D':
                #ax0=fig.gca(projection='3d')
                pass
            else:
                fig_x_pos=float(myConfig.plot_config_array['xlabel_pos'])
                fig_y_pos=float(myConfig.plot_config_array['ylabel_pos'])

                ax0 = fig.add_subplot(mpl.gridspec.GridSpec(nr_plots_y, nr_plots_x)[0])
                if myPlottypeTable['0'].find('hexbins')!=-1 and myConfig.plot_config_array['histo_panels']!='yes':
                    fig.subplots_adjust(hspace=0.0, wspace=0.0, left=0.16, bottom=0.01, right=0.95, top=0.95)
                    x_pos=0.05
                    y_pos=0.05
                elif myConfig.plot_config_array['histo_panels']=='yes':
                    fig.subplots_adjust(hspace=0.0, wspace=0.0, left=0.16, bottom=0.12, right=0.80, top=0.80)
                    x_pos=0.05
                    y_pos=0.0
                    fig_x_pos=0.5
                    fig_y_pos=0.5                      
                else:
                    fig.subplots_adjust(hspace=0.0, wspace=0.0, 
                                        left=float(myConfig.plot_config_array['adjust_left']), 
                                        bottom=float(myConfig.plot_config_array['adjust_bottom']), 
                                        right=float(myConfig.plot_config_array['adjust_right']), 
                                        top=float(myConfig.plot_config_array['adjust_top']))                      
                    
                    x_pos=float(myConfig.plot_config_array['xlabel_pad'])
                    y_pos=float(myConfig.plot_config_array['ylabel_pad'])                 

            #ajust x- or y-title if there are special symbols found
            latex_code_times=find_latex_code_times('')
                
            if myConfig.plot_config_array['subtitle']!='False':        
                fig.suptitle(myConfig.plot_config_array['title'], fontsize=float(myConfig.plot_config_array['axis_label_fontsize']), color='k', y=0.94, fontname='Bitstream Vera Sans Mono')
            if myConfig.plot_config_array['x_title']!='False' and plot_key.find('tarsel')==-1:                   
                fig.text(fig_x_pos, x_pos, myConfig.plot_config_array['x_title'], ha='center', fontsize=float(myConfig.plot_config_array['axis_label_fontsize']))
            if myConfig.plot_config_array['y_title']!='False' and plot_key.find('tarsel')==-1:             
                fig.text(y_pos, fig_y_pos, myConfig.plot_config_array['y_title'], va='center', rotation='vertical', fontsize=float(myConfig.plot_config_array['axis_label_fontsize']))
            if print_redshift==False:
                print_redshift=myConfig.plot_config_array['print_redshift']
            if print_redshift2==False:
                print_redshift2=myConfig.plot_config_array['print_redshift2']
            if print_redshift3==False:
                print_redshift3=myConfig.plot_config_array['print_redshift3']                
            if myConfig.plot_config_array['print_redshift']!='False':  
                fig.text(float(myConfig.plot_config_array['z_print_position_x']), float(myConfig.plot_config_array['z_print_position_y']), print_redshift, fontsize=float(myConfig.plot_config_array['text_fontsize'])) #+r'$\times 10^{-3}$ $h^3$ $Mpc^{-3}$'      
            if myConfig.plot_config_array['print_redshift2']!='False':  
                fig.text(float(myConfig.plot_config_array['z_print_position_x2']), float(myConfig.plot_config_array['z_print_position_y2']), print_redshift2, fontsize=float(myConfig.plot_config_array['text_fontsize'])-5) #+r'$\times 10^{-3}$ $h^3$ $Mpc^{-3}$'      
            if myConfig.plot_config_array['print_redshift3']!='False':  
                fig.text(float(myConfig.plot_config_array['z_print_position_x3']), float(myConfig.plot_config_array['z_print_position_y3']), print_redshift3, fontsize=float(myConfig.plot_config_array['text_fontsize'])-5) #+r'$\times 10^{-3}$ $h^3$ $Mpc^{-3}$'      
 

             
            plot(ax0, 
                 mydata, 
                 prefix='',
                 add_obs=my_add_subplot,                
                 loops=loops)
          
            #plot user defined axis ticks                
            plt.setp(ax0.get_xticklabels(), visible=True)
            plt.setp(ax0.get_yticklabels(), visible=True)
            plot_axis_ticks(ax0)

            for axis in ['top','bottom','left','right']:
                ax0.spines[axis].set_linewidth(float(myConfig.plot_config_array['myframe_lw'])+1)                

        else:

            fig = plt.figure(figsize=(float(myConfig.plot_config_array['size_x']),float(myConfig.plot_config_array['size_y'])), dpi=300)
            if myConfig.plot_config_array['subtitle']!='False':        
                fig.suptitle(myConfig.plot_config_array['title'], fontsize=float(myConfig.plot_config_array['axis_label_fontsize']), color='k', y=0.94, fontname='Bitstream Vera Sans Mono')

            fig.subplots_adjust(hspace=0.0, wspace=0.0, 
                    left=float(myConfig.plot_config_array['adjust_left']), 
                    bottom=float(myConfig.plot_config_array['adjust_bottom']), 
                    right=float(myConfig.plot_config_array['adjust_right']), 
                    top=float(myConfig.plot_config_array['adjust_top']))   
 
            gs = mpl.gridspec.GridSpec(nr_plots_y, nr_plots_x, height_ratios=(int(myConfig.plot_config_array['hratio1']),int(myConfig.plot_config_array['hratio2'])))            


            while k<nr_plots_y:
                l=0
                while l<nr_plots_x:
                    print ''
                    print 'MAIN START!+++++++++++++++++++++++++++++++++++++++++++++++'                
                   
                    key = fig.add_subplot(gs[k])
                    for axis in ['top','bottom','left','right']:
                        key.spines[axis].set_linewidth(float(myConfig.plot_config_array['myframe_lw']))                     
                    
                    print 'k:', k, 'l:', l, 'axis ID:', key, 'name:', mydata[k]['name'], 'add_obs:', mydata[k]['add_obs'], 'print_redshift:', mydata[k]['print_redshift'], 'plot_legend:', mydata[k]['plot_legend']                     
                    
                    if mydata[k]['print_redshift']!=False:
                        key.text(float(myConfig.plot_config_array['z_print_position_x']), float(myConfig.plot_config_array['z_print_position_y']), mydata[k]['print_redshift'], fontsize=float(myConfig.plot_config_array['text_fontsize'])) #+r'$\times 10^{-3}$ $h^3$ $Mpc^{-3}$'                        
                  
                    #print mydata[k]['add_obs']

                    plot(key, 
                         mydata[k]['data'], 
                         prefix=mydata[k]['name'],
                         add_obs=mydata[k]['add_obs'],
                         loops=loops)

                    latex_code_times=find_latex_code_times(mydata[k]['name'])                    
                    if k==0 and l==0:
                        plt.setp(key.get_xticklabels(), visible=False)    
                    #print 'HERE: 998',  str(mydata[l]['print_redshift'])
                    labelpad_xaxis=6
                    if k==0: 
                        labelpad=18
                    else:
                        labelpad=0
                        
                    if str(mydata[l]['print_redshift']).find('*')!=-1:
                        key.set_ylabel(myConfig.plot_config_array[mydata[k]['name']+'y_title']+myConfig.plot_config_array[mydata[k]['name']+'cut_part1']+latex_code_times+myConfig.plot_config_array[mydata[k]['name']+'cut_part2'], color='k', fontsize=int(myConfig.plot_config_array['axis_label_fontsize']), labelpad=labelpad)
                        plt.setp(key.get_yticklabels(), visible=True)

                    else:                        
                        #plot user defined axis ticks                
                        plot_axis_ticks(key)
                        key.set_ylabel(myConfig.plot_config_array[mydata[k]['name']+'y_title']+myConfig.plot_config_array[mydata[k]['name']+'cut_part1']+latex_code_times+myConfig.plot_config_array[mydata[k]['name']+'cut_part2'], color='k', fontsize=int(myConfig.plot_config_array['axis_label_fontsize']), labelpad=labelpad)                           

                    if k==nr_plots_y-1 and str(mydata[l]['print_redshift']).find('centralss')!=-1:
                        key.set_xlabel(myConfig.plot_config_array[mydata[k]['name']+'x_title'], color='k', fontsize=int(myConfig.plot_config_array['axis_label_fontsize']))        
                        plt.setp(key.get_xticklabels(), visible=True)
                    else:
                        #plot user defined axis ticks                
                        plot_axis_ticks(key)
                        if k==nr_plots_y-1:
                            key.set_xlabel(myConfig.plot_config_array[mydata[k]['name']+'x_title'], color='k', fontsize=int(myConfig.plot_config_array['axis_label_fontsize']), labelpad=labelpad_xaxis)        

#                    if str(mydata[l]['print_redshift']).find('M_r')!=-1: 
#                        if str(mydata[l]['print_redshift']).find('(a)')!=-1:
#                            key.set_ylabel(myConfig.plot_config_array[mydata[k]['name']+'y_title']+myConfig.plot_config_array[mydata[k]['name']+'cut_part1']+latex_code_times+myConfig.plot_config_array[mydata[k]['name']+'cut_part2'], color='k', fontsize=int(myConfig.plot_config_array['axis_label_fontsize']), labelpad=labelpad)                                                                                                                                
#                            if l==0 and k>0:
#                                plt.setp(key.get_xticklabels(), visible=False) 
#                                plt.setp(key.get_yticklabels(), visible=True)
#        
#                            else:
#                                plt.setp(key.get_xticklabels(), visible=False)                        
#                                plt.setp(key.get_yticklabels(), visible=True)
#                                
#                        elif str(mydata[l]['print_redshift']).find('(b)')!=-1:                                                                 
#                            plt.setp(key.get_xticklabels(), visible=False)                        
#                            plt.setp(key.get_yticklabels(), visible=False)
#    
#                        elif str(mydata[l]['print_redshift']).find('(c)')!=-1:
#                            key.set_ylabel(myConfig.plot_config_array[mydata[k]['name']+'y_title']+myConfig.plot_config_array[mydata[k]['name']+'cut_part1']+latex_code_times+myConfig.plot_config_array[mydata[k]['name']+'cut_part2'], color='k', fontsize=int(myConfig.plot_config_array['axis_label_fontsize']), labelpad=labelpad)                                                                
#                            if l==0 and k>0:
#                                key.set_xlabel(myConfig.plot_config_array[mydata[k]['name']+'x_title'], color='k', fontsize=int(myConfig.plot_config_array['axis_label_fontsize']))        
#                                
#                                plt.setp(key.get_xticklabels(), visible=True)
#                                plt.setp(key.get_yticklabels(), visible=True)
#        
#                            else:
#                                plt.setp(key.get_xticklabels(), visible=False)                        
#                                plt.setp(key.get_yticklabels(), visible=True)
#    
#                        elif str(mydata[l]['print_redshift']).find('(d)')!=-1:
#                            if l==0 and k>0:
#                                key.set_xlabel(myConfig.plot_config_array[mydata[k]['name']+'x_title'], color='k', fontsize=int(myConfig.plot_config_array['axis_label_fontsize']), labelpad=labelpad)        
#                                
#                                plt.setp(key.get_xticklabels(), visible=True)
#                                plt.setp(key.get_yticklabels(), visible=False)
#        
#                            else:
#                                plt.setp(key.get_xticklabels(), visible=False)                        
#                                plt.setp(key.get_yticklabels(), visible=False)                       

                    try: 
                        key.axhline(y=float(myConfig.plot_config_array[mydata[k]['name']+'hline_ypos']), xmin=-100, xmax=100, color='k', ls='--', lw=2.0, zorder=30)
                        print 'plotted axhline!'
                    except:
                        pass
                    
                    try:
                        key.axvline(x=float(myConfig.plot_config_array[mydata[k]['name']+'vline_xpos']), ymin=-100, ymax=100, color='k', ls='--', lw=2.0, zorder=0)
                        print 'plotted axvline!'
                    except:
                        pass                                                 

                    #plt.setp(key.get_yticklabels()[0], visible=False)    
                    #plt.setp(key.get_yticklabels()[-1], visible=False)                        
                        
                    #key.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{0:.1f}".format(x)))
                    #key.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{0:.2f}".format(x)))

                    key.tick_params(axis='x', which='major', top='on', bottom='on', pad=10, labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), length=major_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), direction='in', zorder=20) 
                    key.tick_params(axis='x', which='minor', top='on', bottom='on', pad=10, labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), length=minor_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), direction='in', zorder=20) 

                    key.tick_params(axis='y', which='minor', left='on', right='on', pad=10, labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), length=minor_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), direction='in', zorder=20) 
                    key.tick_params(axis='y', which='major', left='on', right='on', pad=10, labelsize=float(myConfig.plot_config_array['axis_ticks_fontsize']), length=major_ticks_lenght, width=float(myConfig.plot_config_array['myframe_lw']), direction='in', zorder=20) 
                    

                    print '+++++++++++++++++++++++++++++++++++++++++++++++MAIN END!'
                    print ''   
 
                    l+=1                  
                k+=1

        try:
            if myConfig.plot_config_array['float_format_x']=='1f':            
                ax0.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{0:.1f}".format(x)))
            elif myConfig.plot_config_array['float_format_x']=='2f':
                ax0.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{0:.2f}".format(x))) 
            elif myConfig.plot_config_array['float_format_x']=='3f':
                ax0.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{0:.3f}".format(x)))
            elif myConfig.plot_config_array['float_format_x']=='4f':
                ax0.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{0:.4f}".format(x)))
    
            if myConfig.plot_config_array['float_format_y']=='1f':            
                ax0.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{0:.1f}".format(x)))
            elif myConfig.plot_config_array['float_format_y']=='2f':
                ax0.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{0:.2f}".format(x)))
            elif myConfig.plot_config_array['float_format_y']=='3f':
                ax0.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{0:.3f}".format(x)))
            elif myConfig.plot_config_array['float_format_y']=='4f':
                ax0.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{0:.5f}".format(x)))
        except:
            pass

######################################################################################################################################################################## 
#PRINT OUTPUT
                
        if myConfig.plot_config_array['print_plot'] =='yes':
            print '... print into file!'
            if myConfig.plot_config_array['print_pdf'] =='yes':
                pp = PdfPages(mydir+myfilename+'.pdf')
                if myConfig.plot_config_array['tight_layout']=='True':
                    mybbox_inches='tight'
                else:
                    mybbox_inches=None
                plt.savefig(pp, format='pdf', rasterized=True, dpi=float(myConfig.plot_config_array['quality_dpi']), pad_inches=0.05, bbox_inches=mybbox_inches, transparent=True)
                pp.close()

            else:
                plt.savefig(mydir+myfilename+'.png', format='png', rasterized=True, dpi=500, bbox_inches='tight', pad_inches=0.05)   

        if myConfig.plot_config_array['show_plot'] =='yes': 
            plt.show()

        plt.close(fig)
        
    def writeIntoFile(
            self,
            filename,
            data2write,
            mydelimiter='',
            myheader='',
            data_format='%.18e',
            append_mytext=False,
            data_is_string=False):
        
        print 'wirteINotFile:', filename      
        #print 'myheader:', myheader
        #print data_format
                
        if append_mytext!=True:
            f_handle = file(filename, 'w')  
            if data_is_string==True:
                f_handle.write(data2write)
            else:
                np.savetxt(filename, data2write, fmt=data_format, header=myheader, delimiter=mydelimiter)
        else:
            f_handle = open(filename, 'a+')
            if data_is_string==True:
                f_handle.write(data2write)
            else:
                np.savetxt(f_handle, data2write, fmt=data_format, delimiter=mydelimiter)
        f_handle.close()
            
    def write2HDF5(self,
                   data,
                   myconds_array,
                   myconfig_array,
                   filename,
                   catname,
                   volume=False,
                   redshift=False,
                   scale_factor=False,
                   write_attributes=True):

        print 'wirte2HDF5:', filename 

#        if data[:,0].size>50000:
#            data = mL.choose_random_sample(data, 50000)
        
        with h5py.File(filename, 'w') as hf:

            if write_attributes==True:
                try:
                    hf.attrs['simName'] = myconfig_array[catname+'_simulation_name']
                    hf.attrs['catName'] = catname
                    hf.attrs['redshift'] = redshift
                    hf.attrs['scaleFactor'] = scale_factor
                    hf.attrs['boxSize'] = myconfig_array[catname+'_box_size']
                    hf.attrs['boxSizeUnit'] = 'h-1Mpc'
                    #hf.attrs['physcVolume'] = volume
                    #hf.attrs['physcVolumeUnit'] = 'Mpc3'
                    #hf.attrs['comvVolume'] = float(myconfig_array[catname+'_box_size'])**3
                    #hf.attrs['comvVolumeUnit'] = 'Mpc3h-3'
                    hf.attrs['hubblePar'] = myconfig_array[catname+'_hubble_par']
                    hf.attrs['sampleInfo'] = myconfig_array[catname+'_sample_info']
                    hf.attrs['cosmology'] = myconfig_array[catname+'_cosmology']
                    hf.attrs['annotation'] = myconfig_array[catname+'_my_annotation']
                        
                except:
                    pass
                hf.attrs['dateTimeStamp'] = date_time_stamp
                hf.attrs['author'] = 'DS'
              

            print 'Check attributes written to HDF5...'            
            i=0
            while i<len(hf.attrs.keys()):
                print hf.attrs.keys()[i], hf.attrs.values()[i]
                i+=1
            #print 'HERE:', myconds_array
            #exit()
            i=0
            while i<myconds_array['nr_entries']:
                
                #print 'name:', myconds_array['name'+str(i)], 'col_id:', myconds_array[myconds_array['name'+str(i)]+'_col_id'], myconds_array[myconds_array['name'+str(i)]+'_output_col_name']

                if myconds_array[myconds_array['name'+str(i)]+'_exclude']!='yes':
                    if myconds_array['name'+str(i)].startswith('mag')!=-1: myconds_array['name'+str(i)]+'_no_K-corr'

                    if myconds_array[myconds_array['name'+str(i)]+'_output_col_name']!='False':
                        col_name=myconds_array[myconds_array['name'+str(i)]+'_output_col_name']
                    else:
                        col_name=myconds_array['name'+str(i)]
                    
                    dset = hf.create_dataset(col_name, (len(data[myconds_array['name'+str(i)]]),), dtype=data[myconds_array['name'+str(i)]].dtype)
                    
                    dset[0::] = data[myconds_array['name'+str(i)]]
                    
                    dset.attrs['colName'] = col_name
                        
                    try:
                        dset.attrs['unit'] = myconds_array[myconds_array['name'+str(i)]+'_unit']
                        dset.attrs['cutMinValue'] = myconds_array[myconds_array['name'+str(i)]+'_min']
                        dset.attrs['cutMaxValue'] = myconds_array[myconds_array['name'+str(i)]+'_max']
                    except:
                        pass
                
                i+=1
        

#The later functions are form Stackoverflow Question: "Text box with line wrapping in matplotlib?" edited Oct 8 '14 at 1:41

    def wrapText(self, text, margin=4):

    # Text Wrapping
    # Defines wrapText which will attach an event to a given mpl.text object,
    # wrapping it within the parent axes object.  Also defines a the convenience
    # function textBox() which effectively converts an axes to a text box.
        """ Attaches an on-draw event to a given mpl.text object which will
            automatically wrap its string wthin the parent axes object.
    
            The margin argument controls the gap between the text and axes frame
            in points.
        """
        ax = text.get_axes()
        margin = margin / 72 * ax.figure.get_dpi()
    
        def _wrap(event):
            """Wraps text within its parent axes."""
            def _width(s):
                """Gets the length of a string in pixels."""
                text.set_text(s)
                return text.get_window_extent().width
    
            # Find available space
            clip = ax.get_window_extent()
            x0, y0 = text.get_transform().transform(text.get_position())
            if text.get_horizontalalignment() == 'left':
                width = clip.x1 - x0 - margin
            elif text.get_horizontalalignment() == 'right':
                width = x0 - clip.x0 - margin
            else:
                width = (min(clip.x1 - x0, x0 - clip.x0) - margin) * 2
    
            # Wrap the text string
            words = [''] + _splitText(text.get_text())[::-1]
            wrapped = []
    
            line = words.pop()
            while words:
                line = line if line else words.pop()
                lastLine = line
    
                while _width(line) <= width:
                    if words:
                        lastLine = line
                        line += words.pop()
                        # Add in any whitespace since it will not affect redraw width
                        while words and (words[-1].strip() == ''):
                            line += words.pop()
                    else:
                        lastLine = line
                        break
    
                wrapped.append(lastLine)
                line = line[len(lastLine):]
                if not words and line:
                    wrapped.append(line)
    
            text.set_text('\n'.join(wrapped))
    
            # Draw wrapped string after disabling events to prevent recursion
            handles = ax.figure.canvas.callbacks.callbacks[event.name]
            ax.figure.canvas.callbacks.callbacks[event.name] = {}
            ax.figure.canvas.draw()
            ax.figure.canvas.callbacks.callbacks[event.name] = handles
    
        ax.figure.canvas.mpl_connect('draw_event', _wrap)
    
    def _splitText(text):
        """ Splits a string into its underlying chucks for wordwrapping.  This
            mostly relies on the textwrap library but has some additional logic to
            avoid splitting latex/mathtext segments.
        """
        import textwrap
        import re
        math_re = re.compile(r'(?<!\\)\$')
        textWrapper = textwrap.TextWrapper()
    
        if len(math_re.findall(text)) <= 1:
            return textWrapper._split(text)
        else:
            chunks = []
            for n, segment in enumerate(math_re.split(text)):
                if segment and (n % 2):
                    # Mathtext
                    chunks.append('${}$'.format(segment))
                else:
                    chunks += textWrapper._split(segment)
            return chunks
    
    def textBox(self, text, axes, ha='left', fontsize=12, margin=None, frame=True, **kwargs):
        """ Converts an axes to a text box by removing its ticks and creating a
            wrapped annotation.
        """
        if margin is None:
            margin = 6 if frame else 0
        axes.set_xticks([])
        axes.set_yticks([])
        axes.set_frame_on(frame)
    
        an = axes.annotate(text, fontsize=fontsize, xy=({'left':0, 'right':1, 'center':0.5}[ha], 1), ha=ha, va='top',
                           xytext=(margin, -margin), xycoords='axes fraction', textcoords='offset points', **kwargs)
        wrapText(an, margin=margin)
        return an