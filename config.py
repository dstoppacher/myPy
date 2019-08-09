#configure a structure array where configuration details are stored
#--------------------------------------------
import os
system_info = os.getcwd()
start = system_info.find('anaconda')
mycomp = system_info[:start]

import arangeData as aD
import numpy as np
import myLib as mL

class Configuration:
    
    def __init__(self):
        
        myData = aD.ArangeData()
        data_array = myData.readAnyFormat(config=False, mydtype=np.str_, mypath=mycomp+'anaconda/pro/myRun/run_config/mypath_handler.txt', data_shape='shaped', comment='#')

        self.mypath_handler ={}
        j=0
        while j<data_array[:,0].size:
            self.mypath_handler[str(data_array[j,0])] = data_array[j,1]                     
            j+=1        
    
    def readMyRunConfig(
            self,
            read_filenr_sequence=False,
            myfile_format='.txt'):
 
       
        data = {}
        self.snapid_array = {}
        
        myData = aD.ArangeData()
 
        data_array = myData.readAnyFormat(config=False, mydtype=np.str_, mypath=self.mypath_handler['RUNFILE'], data_shape='shaped', comment='#')

        my_count=np.char.count(data_array, ':')
        data['nr_cats'] = np.sum(my_count)

        j=0
        count=0
        while j<np.sum(my_count):
            
            i=0 # i=redshift
            a=0 # a=filenr
            #print 'j:', j, 'i:', i, 'count:', count
            start = data_array[i+count].find(':')
            key_catname='catname'+str(j)
            data[key_catname]                     = data_array[i+count][start+1::]
            i+=1
            

            key = data[key_catname]+'_mysoftlink_dir'
            data[key] = data_array[i+count]            
            #print 'data[key]:', data[key]
            i+=1
            
            key = data[key_catname]+'_mydir'
            data[key] = data_array[i+count]
            #print 'data[key]:', data[key]
            
            i+=1
            cat_test=-1
            
            data = self.myConfigArray(data, data[key]+data[key_catname]+'_config.txt', data[key_catname])

            self.snapIdzRed(self.snapid_array, data[key]+'snapidzred.txt', data[key_catname], data[data[key_catname]+'_inputfilename_part1'], data[data[key_catname]+'_snapid0'], data[data[key_catname]+'_inputfilename_part2'], myfile_format, data[data[key_catname]+'_use_snapidz_mapping'], data)

            if read_filenr_sequence==True:
                read_filenr_from= int(data_array[i+count][len(data[data[key_catname]+'_inputfilename_part1'])::])
                read_filenr_until= int(data_array[i+count+1][len(data[data[key_catname]+'_inputfilename_part1'])::])
                
                data['read_filenr_from'] = read_filenr_from
                data['read_filenr_until'] = read_filenr_until
                
                #print 'from:', read_filenr_from, 'until:', read_filenr_until
            
            while cat_test <= -1 and i+count!= data_array.size:
                #print 'CAT_test! i:', i, 'a:', a
                    
                if read_filenr_sequence==True:
                    
                    b=read_filenr_from
                    while b<read_filenr_until+1:
                        key = data[key_catname] + '_filename' + str(b)                      
                        data[key]=data[key_catname] + '_' + data[data[key_catname]+'_inputfilename_part1']+ str(b)

                        b+=1

                else:
                    #print 'key:', key
                    key = data[key_catname] + '_filename' +str(a)                   
                    data[key]=data[key_catname] + '_' + data_array[i+count]
                    #print 'data[key]:', data[key]
                                        
                    if key.find('DIV')!=-1:
                        self.snapIdzRed(self.snapid_array, data[key]+'snapidzred.txt', data[key_catname], data[data[key_catname]+'_inputfilename_part1'], data[data[key_catname]+'_snapid'], data[data[key_catname]+'_inputfilename_part2'], myfile_format, data[data[key_catname]+'_use_snapidz_mapping'], key_filename=data[key])

                    
                i+=1
                a+=1


                if i+count!=data_array.size: 
                    #print 'i+count:', i+count, 'mydata size', data_array.size
                    cat_test = data_array[i+count].find(':')
                    #print 'cat_test:', cat_test
            count=count+i
            j+=1
        
        #Count nr of Redshifts of
        data['nr_zs']= a    
        #print 'snapid_array:', self.snapid_array
        self.config_array = data
        #print self.config_array
 
    def myConfigArray(self,
                      data,
                      myconfig_datafile,
                      catkey):
               
        data_array = mL.myMultiColDetector(myconfig_datafile)

        myData = aD.ArangeData()
        unit_corr_array= myData.readAnyFormat(config=False, mydtype=np.str_, mypath=self.mypath_handler['MY_UNIT_CORRECT_FILE'], data_shape='shaped', skiprow=0, delim=' ')
        
        if np.prod(unit_corr_array.shape)==3:
            unit_corr_array = np.expand_dims(unit_corr_array, axis=0)
            
#        analyse_tarsel_array= myData.readAnyFormat(config=False, mydtype=np.str_, mypath=self.mypath_handler['ANALYSE_TARSEL_CONFIGFILE'], data_shape='shaped', skiprow=0, delim=' ')
#        
#        if np.prod(unit_corr_array.shape)==2:
#            analyse_tarsel_array = np.expand_dims(analyse_tarsel_array, axis=0)
                       
        unit_map_array={}        
        i=0
        while i<unit_corr_array[:,0].size:
            
            multicol_test = mL.multicolTestAlgorithm(unit_corr_array[i,1])

            c=0
            while c<multicol_test.size:              
                unit_map_array[catkey+'_unit_corr_'+unit_corr_array[i,0]+str(c)] = multicol_test[c]
                unit_map_array[catkey+'_data_type_'+unit_corr_array[i,0]] = unit_corr_array[i,2]
                unit_map_array[catkey+'_unit_corr_'+unit_corr_array[i,0]+'_nr_corrs'] = c+1
                
                #print 'i:', i, 'name:', unit_corr_array[i,0], 'c:', c,'multicol_test[c]:', multicol_test[c], '[i,1]:', unit_corr_array[i,1],'[i,2]:', unit_corr_array[i,2], unit_map_array[catkey+'_unit_corr_'+unit_corr_array[i,0]+str(c)]
                #print 'data_type:', unit_map_array[catkey+'_data_type_'+unit_corr_array[i,0]], 'unit_corr:', unit_map_array[catkey+'_unit_corr_'+unit_corr_array[i,0]+str(c)], 'nr_corrs:',  unit_map_array[catkey+'_unit_corr_'+unit_corr_array[i,0]+'_nr_corrs']
                #print '------------'
                #print ' '
                c+=1
            i+=1
            #unit_map_array[catkey+'_unit_corr_'+unit_corr_array[i,0]] = unit_corr_array[i+c,1]
            #unit_map_array[catkey+'_data_type_'+unit_corr_array[i,0]] = unit_corr_array[i+c,2]
            #i+=1
        data[catkey+'_id_col_array'] = {}
        data[catkey+'_halo_id_col_array'] = {}
        data[catkey+'_analyse_tarsel_id_col_array'] = {}


        i=0
        x=0 
        y=0
        z=0
        nr_cols2read_col=0
        nr_cols2read_halo_col=0
        nr_cols2read_tarsel_col=0
        while i<data_array[:,0].size:
            
            #print 'i:', i, 'a:', data_array[i,1], data_array[i,:]
            #print 'name', data_array[i,2], 'value:', data_array[i,0]
            #print data_array[i,2]

            data[catkey+'_'+data_array[i,2]] = mL.check_datatype(data_array[i,0])

            if data_array[i,2].startswith('col_'):
                #print 'here: assembly id_col_array!'
                if data_array[i,0]!=str(99):
                    nr_cols2read_col+=1
                    p=0
                    while p<unit_map_array[catkey+'_unit_corr_'+data_array[i,2][4::]+'_nr_corrs']:
                        #print 'x:', x, data_array[i,2][4::], 'p:', 'nr_corrs:', unit_map_array[catkey+'_unit_corr_'+data_array[i,2][4::]+'_nr_corrs']
                        #print unit_map_array[catkey+'_unit_corr_'+data_array[i,2][4::]+str(p)]
                        #print data[catkey+'_unit_corr_'+unit_map_array[catkey+'_unit_corr_'+data_array[i,2][4::]+str(p)]]
                        data[catkey+'_id_col_array'].update({'name'+str(x): data_array[i,2][4::], data_array[i,2][4::]+'_col_id': data[catkey+'_'+data_array[int(data_array[i,1]),2]], str(data[catkey+'_'+data_array[int(data_array[i,1]),2]])+'_name'+str(x): data_array[i,2][4::],  'data_type'+str(x): unit_map_array[catkey+'_data_type_'+data_array[i,2][4::]], 'col_id'+str(x): data[catkey+'_'+data_array[int(data_array[i,1]),2]], 'corr_type'+str(x): unit_map_array[catkey+'_unit_corr_'+data_array[i,2][4::]+str(p)], 'unit_corr'+str(x): data[catkey+'_unit_corr_'+unit_map_array[catkey+'_unit_corr_'+data_array[i,2][4::]+str(p)]], 'conv_to_AB_mag'+str(x): data[catkey+'_conv_to_AB_mag']})
                        p+=1                        
                        x+=1
                                          
            if data_array[i,2].startswith('halo_col_'):
                #print 'here: assembly id_col_array!'
                nr_cols2read_halo_col+=1
                if data_array[i,0]!=str(99):
                    
                    #print 'y:', y, data_array[i,2][9::]
                    #print unit_map_array[catkey+'_unit_corr_'+data_array[i,2][9::]]
                    #print data[catkey+'_unit_corr_'+unit_map_array[catkey+'_unit_corr_'+data_array[i,2][9::]]]
                                  
                    data[catkey+'_halo_id_col_array'].update({'name'+str(y): data_array[i,2][4::], data_array[i,2][9::]+'_col_id': data[catkey+'_'+data_array[int(data_array[i,1]),2]],  'data_type'+str(y): unit_map_array[catkey+'_data_type_'+data_array[i,2][9::]], 'col_id'+str(y): data[catkey+'_'+data_array[int(data_array[i,1]),2]], 'corr_type'+str(y): unit_map_array[catkey+'_unit_corr_'+data_array[i,2][9::]], 'unit_corr'+str(y): data[catkey+'_unit_corr_'+unit_map_array[catkey+'_unit_corr_'+data_array[i,2][9::]]], 'dummy'+str(y): 'dummy'})
                    y+=1
                

            if data_array[i,2].startswith('analyse_tarsel_'):
                #print 'here: assembly id_col_array!'
                                                
                nr_cols2read_tarsel_col+=1
                name_map={}
                k=1
                find_name_before=len('analyse_tarsel_')

                if data_array[i,2].find('hist')!=-1:
                    #print 'analyse_tarsel_HIST!!!'                 
                    name_map.update({'col_name1': data[catkey+'_id_col_array']['name0'], 'id_col_1': 0, 'col_name2': data[catkey+'_id_col_array']['name0'], 'id_col_2': 0})
                
                else:
                    
                    while k<3:
                        #print 'name', data_array[i,2], 'value:', data_array[i,0], 'k:', k, 'find_name_before:', find_name_before, 'len of name to find:', len(data_array[i,2])
                        find_name = data_array[i,2][find_name_before::].find('-')
    
                        if find_name!=-1:
                            try:
                                #print data_array[i,2][find_name_before+find_name-1:find_name_before+find_name+2]
                                find_band_name = data_array[i,2][find_name_before+find_name-1:find_name_before+find_name+2]
                                find_minus_sign=find_band_name.find('-')
                                #correct the minus sign
                                #print 'minus:', find_minus_sign
                                col_name='mAB_dA_total_cut_'+find_band_name[find_minus_sign-1:find_minus_sign]+'_'+find_band_name[find_minus_sign+1::]
                                #print 'col_name:', col_name

                                try:
                                    id_col=data[catkey+'_id_col_array'][col_name+'_col_id']
                                except:
                                    try:
                                        #print '-: try2:',
                                        col_name='mAB_dA_total_'+find_band_name[0]
                                        id_col=data[catkey+'_id_col_array'][col_name+'_col_id']
                                        #print 'col_name:', col_name
                                    except:
                                        col_name='MAB_dA_total_'+find_band_name[0]
                                        id_col=data[catkey+'_id_col_array'][col_name+'_col_id']
                                #print 'col_name:', col_name
                            except:
                                id_col = 99
                                data[catkey+'_id_col_array'].update({'name': col_name, col_name+'_col_id': id_col})
                           
                            find_name=3
                        else:
                            #print 'find_name_before:', find_name_before
                            find_name = data_array[i,2][find_name_before::].find('_')
                            #print 'find_name else:', find_name
                            if find_name==0: 
                                find_name_before+=1
                                find_name=len(data_array[i,2])-find_name_before
                            #print 'find_name+1:', find_name_before, 'find_name:', find_name
                            try:
                                try:
                                    #print 'col name --> mAB!'
                                    col_name='mAB_dA_total_'+data_array[i,2][find_name_before:find_name_before+find_name]
                                    print 'try:', data_array[i,2][find_name_before:find_name_before+find_name],  'col_name:', col_name,data[catkey+'_id_col_array'][col_name+'_col_id']
                                except:
                                    try:
                                        print 'try2:',
                                        col_name='mAB_dA_total_'+find_band_name[0]
                                        id_col=data[catkey+'_id_col_array'][col_name+'_col_id']
                                        #print 'col_name:', col_name
                                    except:
                                        col_name='MAB_dA_total_'+find_band_name[0]
                                        id_col=data[catkey+'_id_col_array'][col_name+'_col_id']
                                #print 'col_name:', col_name
                            except:
                                try:
                                    if data_array[i,2][find_name_before:find_name_before+find_name].find('I')!=-1:
                                        #print 'col name --> MAB!'
                                        col_name='MAB_total_'+data_array[i,2][find_name_before+1:find_name_before+1+find_name-1]
                                    else:
                                        col_name='mAB_total_'+data_array[i,2][find_name_before:find_name_before+find_name]
                                    print 'try3:', data[catkey+'_id_col_array'][col_name+'_col_id']
                                                            
                                except:
                                    try:                                
                                        col_name='mAB_dA_total_cut_'+data_array[i,2][find_name_before:find_name_before+find_name]
                                        print data[catkey+'_id_col_array'][col_name+'_col_id']
                                    
                                    except:
                                        if col_name.find('rhalfmass')!=-1:
                                            col_name='rhalf_mass'
                                        elif col_name.find('rhalfd')!=-1:
                                            col_name='rhalf_disk'
                                        elif col_name.find('rhalfb')!=-1:
                                            col_name='rhalf_bulge'
                                        elif col_name.find('NFW')!=-1:
                                            col_name='NFW_con'
                                        elif col_name.find('mhalo200c')!=-1:
                                            col_name='mhalo_200c'
                                        elif col_name.find('mbasic200c')!=-1:
                                            col_name='mbasic_200c'
                                        elif col_name.find('angMdisk')!=-1:
                                            col_name='angM_disk'
                                        elif col_name.find('angMspheroid')!=-1:
                                            col_name='angM_spheroid'
                                        elif col_name.find('mbar')!=-1:
                                            col_name='mstar'
                                        elif col_name.find('bdisk')!=-1:
                                            col_name='angM_disk'                                             
                                        elif col_name.find('bbulge')!=-1:
                                            col_name='angM_spheroid'                                                 
                                        else:
                                            col_name=data_array[i,2][find_name_before:find_name_before+find_name]
                                            
                                        print 'except:', col_name
                                        try:
                                            #print 'except2:', col_name, 
                                            if catkey=='SAGE_1Gpc' and col_name=='mcold': col_name='mcold_disk'
                                            #if catkey=='Galacticus_1Gpc' and col_name=='Mzgas': col_name='zgas_spheroid'
                                            #print data[catkey+'_id_col_array'][col_name+'_col_id']
                                            
                                        except:
                                            #print 'col_name:', col_name, 'total name:', data_array[i,2][find_name_before::],
                                            find_name=data_array[i,2][find_name_before+len(col_name)+1::].find('_')

                                            #print 'except3:', 'find name:', find_name, 
                                            col_name=data_array[i,2][find_name_before:find_name_before+len(col_name)+1+find_name]
                                            find_name=len(col_name)
                                            #print 'new col_name:', col_name, 'new find_name_before:', find_name_before,                                                                              
                                            #print data[catkey+'_id_col_array'][col_name+'_col_id']
                        
                        find_name_before+=find_name
                        #name_map.update({'col_name'+str(k): col_name, 'id_col_'+str(k): id_col})
                        name_map.update({'col_name'+str(k): col_name})                        
                        
                        k+=1
                    
                print 'name_map:', name_map

                #data[catkey+'_analyse_tarsel_id_col_array'].update({'name'+str(z): data_array[i,2], data_array[i,2]+'_col_name1': name_map['col_name1'],   data_array[i,2]+'_col_name_id1': data[catkey+'_id_col_array'][name_map['col_name1']+'_col_id'],  data_array[i,2]+'_col_name2': name_map['col_name2'],   data_array[i,2]+'_col_name_id2': data[catkey+'_id_col_array'][name_map['col_name2']+'_col_id']})
                data[catkey+'_analyse_tarsel_id_col_array'].update({'name'+str(z): data_array[i,2], data_array[i,2]+'_col_name1': name_map['col_name1'],  data_array[i,2]+'_col_name2': name_map['col_name2']})
                
                z+=1
            
            if data_array[i,2].startswith('snapid'):
                
                snapidID=int(data_array[i,1])       
                data[catkey+'_snapid0']  = data_array[i,0]
                #print 'i:', i, 'snapid'+str(0), data[catkey+'_snapid'+str(0)]
                
                a=i+1
                count=1
                while a<data_array[:,0].size:
                    if data_array[a,1]==str(snapidID):
                        data[catkey+'_snapid'+str(count)]=data_array[a,0]
                        #print 'i:', i, 'a:', a, data_array[a,:], 'snapid'+str(count), data[catkey+'_snapid'+str(count)]
                        count+=1
                    a+=1
                i+=count-1

                data[catkey+'_nr_zs_count'] = count
                

            i+=1
  
        data[catkey+'_id_col_array'].update({'nr_entries': x})
        data[catkey+'_id_col_array'].update({'nr_cols2read': nr_cols2read_col})
        data[catkey+'_halo_id_col_array'].update({'nr_entries': y})
        data[catkey+'_halo_id_col_array'].update({'nr_cols2read': nr_cols2read_halo_col})
        data[catkey+'_analyse_tarsel_id_col_array'].update({'nr_entries': z})

        return data

    def plotMappingArray(self):

        myData = aD.ArangeData()
        data_array = myData.readAnyFormat(config=False, mydtype=np.str_, mypath=self.mypath_handler['MAP_PLOT_FILE'], data_shape='shaped', skiprow=0, delim=' ')
       
        if np.prod(data_array.shape)==1 or np.prod(data_array.shape)==5:
            data_array = np.expand_dims(data_array, axis=0)
        #print data_array
        self.plot_map_array = {}
        j=0
        self.plot_map_array['run_nr'] = data_array[j,1]
        while j<data_array[:,0].size:
            self.plot_map_array['plot_map_id'+str(j)]           = data_array[j,0]
            self.plot_map_array[data_array[j,0]+'_config']      = data_array[j,1]           
            self.plot_map_array['nbins_'+data_array[j,0]]       = data_array[j,2]
            self.plot_map_array['data_offset_'+data_array[j,0]] = data_array[j,3]
            self.plot_map_array['load_obs_'+data_array[j,0]]    = data_array[j,4]
            j+=1

        self.plot_map_array['nr_plot_keys'] = j         


    def plotMarkerColorKewords(self):
       

        myData = aD.ArangeData()
        data_array = myData.readAnyFormat(config=False, mydtype=np.str_, mypath=self.mypath_handler['MARKER_COLOR_MAP_FILE'], data_shape='shaped', skiprow=0, delim=' ')
 
        if data_array.size==6:
            data_array = np.expand_dims(data_array, axis=0)            
            
        self.plot_marker_color_map_array = {}
        j=0
        while j<data_array[:,0].size:
            self.plot_marker_color_map_array['catname'+str(j)]                   = data_array[j,0]
            self.plot_marker_color_map_array[data_array[j,0]+'_marker_code']     = data_array[j,1]
            if len(data_array[j,2])==1:
                self.plot_marker_color_map_array[data_array[j,0]+'_color_code']      = data_array[j,2]
            else:
                self.plot_marker_color_map_array[data_array[j,0]+'_color_code']      = '#'+data_array[j,2]
            self.plot_marker_color_map_array[data_array[j,0]+'_linestyle_code']  = data_array[j,3] 
            self.plot_marker_color_map_array[data_array[j,0]+'_color_map']  = data_array[j,4]
            self.plot_marker_color_map_array[data_array[j,0]+'_marker_facecolor_code']  = data_array[j,5] 
            j+=1

       
    def plotKeywords(
                self,
                myconfig_datafile='myconfig.txt'):
        
        self.plot_config_array = {}
        myData = aD.ArangeData()
        data_array = myData.readAnyFormat(config=False, mydtype=np.str_, mypath=myconfig_datafile, data_shape='shaped', comment='#', delim='= ')

        i=0
        while i<data_array[:,0].size:
            self.plot_config_array[data_array[i,0]] = data_array[i,1]
            i+=1

        #print self.plot_config_array

    def load_matplot_stylefiles(self,
                                catname,
                                ncolours=0):
    
        color = {}
        linestyle ={}
        marker = {}
        markercol = {}
        colormap = {}
        colorbar = {}
        plottype = {}
        alpha = {}
        plotlegend = {}
        add_xaxis = {}
        add_yaxis = {}
        
        myData = aD.ArangeData()
        
        mypaths = {'0': self.mypath_handler['MATPLOT_LINESTYLE'], 'dic0': linestyle,
                   '1': self.mypath_handler['MATPLOT_COL'], 'dic1': color,
                   '2': self.mypath_handler['MATPLOT_MARKERSTYLE'], 'dic2': marker,
                   '3': self.mypath_handler['MATPLOT_MARKERCOL'], 'dic3': markercol,        
                   '4': self.mypath_handler['MATPLOT_COLORMAP'], 'dic4': colormap,
                   '5': self.mypath_handler['MATPLOT_COLORBAR'], 'dic5': colorbar,
                   '6': self.mypath_handler['MATPLOT_PLOTTYPE'], 'dic6': plottype,
                   '7': self.mypath_handler['MATPLOT_ALPHA'], 'dic7': alpha,
                   '8': self.mypath_handler['MATPLOT_PLOTLEGEND'], 'dic8': plotlegend,
                   '9': self.mypath_handler['MATPLOT_ADD_XAXIS'], 'dic9': add_xaxis,
                   '10': self.mypath_handler['MATPLOT_ADD_YAXIS'], 'dic10': add_yaxis}          
        a=0
        while a<=10:
            data_array = myData.readAnyFormat(config=False, mydtype=np.str_, mypath=mypaths[str(a)], data_shape='shaped', comment='', delim='= ')
            #print 'mypath:', mypaths[str(a)]            
            #print 'data_array:', data_array
            if data_array.size==2:
                data_array = np.expand_dims(data_array, axis=0)
                
            i=0
            while i<data_array[:,0].size:
                mypaths['dic'+str(a)].update({data_array[i,0]: data_array[i,1]})
                i+=1
            a+=1
            
        if ncolours!=0:
            import distinct_colours as cb_col 
            #print 'use color-blind friendly colors!'
           # print catname
            cb_colours = cb_col.get_distinct(ncolours)           
            if ncolours==1:
                if catname.find('Galacticus')!=-1:
                    #print 'Gal'
                    color['0']='#225588'
                elif catname.find('SAG_')!=-1 or str(catname).find('SAG_1Gpc_v2')!=-1:
                    #print 'SAG'
                    color['0']='#CC6677'
                elif catname.find('SAGE')!=-1:
                    #print 'SAGE'
                    color['0']='#DDCC77'
                else:
                    #print 'else'
                    color['0']=cb_colours[0]
            else:
                i=0
                while i<ncolours:
                    color[str(i)]=cb_colours[i]
                    i+=1
            
            
        return linestyle, color, marker, markercol, colormap, colorbar, plottype, alpha, plotlegend, add_xaxis, add_yaxis

    def physicsSpecs(self):
    
        self.physics_specs = {}
        myData = aD.ArangeData()
        data_array = myData.readAnyFormat(config=False, mydtype=np.str_, mypath=self.mypath_handler['MY_PHYSICS_SPECS'], data_shape='shaped', comment='#', delim='= ')

        if data_array.size==2:
            data_array = np.expand_dims(data_array, axis=0)
            
        i=0
        while i<data_array[:,0].size:
            self.physics_specs[data_array[i,0]] = data_array[i,1]
            i+=1

    def histoConfigs(self):
    
        self.histo_config_array = {}
        myData = aD.ArangeData()
        data_array = myData.readAnyFormat(config=False, mydtype=np.str_, mypath=self.mypath_handler['HISTO_CONFIGFILE'], data_shape='shaped', comment='#', delim=': ')

        if data_array.size==2:
            data_array = np.expand_dims(data_array, axis=0)
            
        i=0
        while i<data_array[:,0].size:
            self.histo_config_array[data_array[i,0]] = data_array[i,1]
            i+=1

    def nameConvMap(self):
    
        self.name_conv_map = {}
        myData = aD.ArangeData()
        data_array = myData.readAnyFormat(config=False, mydtype=np.str_, mypath=self.mypath_handler['MY_NAME_CONV_MAP_FILE'], data_shape='shaped', comment='#', delim=' ')

        if data_array.size==2:
            data_array = np.expand_dims(data_array, axis=0)
            
        i=0
        while i<data_array[:,0].size:
            #print 'i:', i, data_array[i,0], '-->', data_array[i,1]
            self.name_conv_map[data_array[i,0]] = data_array[i,1]
            i+=1 
    
    def snapIdzRed(
            self,
            data,
            myconfig_datafile,
            catkey,
            inputfilename_part1,
            snapid,
            inputfilename_part2,
            file_format,
            use_snapidz_mapping,
            config_array):
        #print 'HERE!', use_snapidz_mapping        
        if use_snapidz_mapping=='True':
            mySnap = aD.ArangeData()
            if catkey.startswith('sussing')==-1:
                snap_array = mySnap.readAnyFormat(config=False, mypath=myconfig_datafile, data_shape='unshaped', nr_col=6, mydat_off=0, nr_rows=120, skiprow=1)
            else:
                snap_array = mySnap.readAnyFormat(config=False, mypath=myconfig_datafile, data_shape='shaped', comment='#', delim='\t') 
                
            #print 'snap_array'
          
            i=0
            while i<snap_array[:,0].size:
                
                if int(snap_array[i,0])<10:
                    if len(snapid)==3: #061
                         name_correction = '00'
                    else:
                         name_correction = '0'
                
                elif int(snap_array[i,0])>=100:
                        name_correction = ''
                
                else:
                   if len(snapid)==3:
                         name_correction = '0'
                   else:
                         name_correction = ''
                
                if config_array[catkey+'_data_format']=='HDF5':
                    key = inputfilename_part1 +'_'+ name_correction + str(int(snap_array[i,0]))
                else:
                    key = inputfilename_part1 + name_correction + str(int(snap_array[i,0])) +  inputfilename_part2 + file_format
                
                a=0
                while a<config_array[catkey+'_nr_zs_count']:
                    if catkey.startswith('sussing')==-1:
                        self.snapid_array[catkey+'_'+key+'_snapid'+str(a)] = {'id': key, 'a': snap_array[i,1], 'z': snap_array[i,2], 't(t0)': snap_array[i,3], 't(year)': snap_array[i,4]}
                    else:
                        self.snapid_array[catkey+'_'+catkey+'_snapid'+str(a)] = {'id': key, 'a': snap_array[i,2], 'z': snap_array[i,1], 't(year)': snap_array[i,3]}
                        
                    a+=1                
                    
                i+=1

        else:
                i=0
                while i<config_array[catkey+'_nr_zs_count']:
                    #print 'key_filename:', inputfilename_part1+snapid+inputfilename_part2

                    if config_array[catkey+'_create_subcat']=='True':
                        #print 'create subcat!'
                        #print 'scale_factor_map:', self.SAM_scale_factor_map
                        my_z=self.SAM_scale_factor_map[catkey+'_redshift'+str(i)]
                    
                    elif config_array[catkey+'_use_snapidz_mapping']=='False' and config_array[catkey+'_load_from_file']=='False' and config_array[catkey+'_load_subcat']!='True':                  
                        #print 'snapidz_mapping False!'
                        my_z=config_array[catkey+'_snapid'+str(i)]
      
                    elif config_array[catkey+'_load_from_file']=='True' or config_array[catkey+'_load_subcat']=='True':
                        #print 'manual z:', config_array[catkey+'_manual_input_redshift']
                        my_z=config_array[catkey+'_manual_input_redshift']
                             
                    self.snapid_array[catkey+'_'+inputfilename_part1+'_snapid'+str(i)] = {'z': my_z}
                    i+=1               
        #print self.snapid_array
        
    def filterDataConfig(self):

        self.mycond_config_array = {}
        myData = aD.ArangeData()
        data_array = myData.readAnyFormat(config=False, mydtype=np.str_, mypath=self.mypath_handler['MYCUT_VALUES_CONFIGFILE'], data_shape='shaped', comment='#', delim=' ')

        if data_array.size==10:
            data_array = np.expand_dims(data_array, axis=0)
        #print data_array
        i=0
        while i<data_array[:,0].size:
            self.mycond_config_array['name'+str(i)] = data_array[i,0]
            self.mycond_config_array[data_array[i,0]+'_min'] = data_array[i,1]
            self.mycond_config_array[data_array[i,0]+'_max'] = data_array[i,2]
            self.mycond_config_array[data_array[i,0]+'_unit'] = data_array[i,3]
            self.mycond_config_array[data_array[i,0]+'_data_type'] = data_array[i,4]
            self.mycond_config_array[data_array[i,0]+'_format'] = data_array[i,5]
            self.mycond_config_array[data_array[i,0]+'_col_id'] = data_array[i,6]
            self.mycond_config_array[data_array[i,6]+'_name']= data_array[i,0]
            self.mycond_config_array[data_array[i,0]+'_name_in_plot'] = data_array[i,7]
            self.mycond_config_array[data_array[i,0]+'_exclude'] = data_array[i,8]
            self.mycond_config_array[data_array[i,0]+'_output_col_name'] = data_array[i,9]
            i+=1
        
        self.mycond_config_array['nr_entries'] = i
        #print self.mycond_config_array
        return self.mycond_config_array

    def getCUTEParameterFile(self):

        self.CUTE_params = {}
        myData = aD.ArangeData()
        data_array = myData.readAnyFormat(config=False, mydtype=np.str_, mypath=self.mypath_handler['MYCUTE_PARAMFILE'], data_shape='shaped', comment='#', delim='= ')
       
        i=0
        while i<data_array[:,0].size:
            self.CUTE_params['name'+str(i)] = data_array[i,0]
            self.CUTE_params[data_array[i,0]+'_param'] = data_array[i,1] 
            i+=1
            
        self.CUTE_params['nr_entries'] = i
        
        return self.CUTE_params


    def SAMScaleFactorMapping(self):

        self.SAM_scale_factor_map = {}
        myData = aD.ArangeData()
        data_array = myData.readAnyFormat(config=False, mydtype=np.str_, mypath=self.mypath_handler['SAM_SCALE_FACTOR_MAPPING'], data_shape='shaped', comment='#', delim=' ')
        #print data_array
        def check_data(data_array):
        
            try:
                if data_array.size==4:
                    data_array = np.expand_dims(data_array, axis=0)
                    
                if data_array[0,2]!='None':
                    i=0
                    while i<data_array[:,0].size:
                        self.SAM_scale_factor_map[data_array[i,0]+'_redshift'+str(i)] = float(data_array[i,2])
                        self.SAM_scale_factor_map[data_array[i,0]+'_snapid'+str(i)] = data_array[i,1]
                        self.SAM_scale_factor_map[data_array[i,0]+'_scale_factor'+str(i)] = float(data_array[i,3])
                        i+=1
            except:
                print 'data_array=None'

        #print 'HERE:', data_array
        

        check_data(data_array)
           
        #print self.SAM_scale_factor_map 
        return self.SAM_scale_factor_map  
                  