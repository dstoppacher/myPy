# Load packages 
from __future__ import print_function
import numpy as np
import sys
import config as conf
import outputData as oD
import arangeData as aD
import myFuncs as mF
import myLib as mL

#from cosmolopy import cd, fidcosmo
#import loadObs as lO
import subprocess as subs

myOutput = oD.OutputData(config=False)
myData = aD.ArangeData()
myConfig = conf.Configuration()
myFuncs = mF.MyFunctions()

import os
system_info = os.getcwd()
start = system_info.find('anaconda')
mycomp = system_info[:start]

import time
ts = time.time()
import datetime 
date_time_stamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H:%M:%S')

from time import time
        
class Pipeline:
    
    def __init__(
            self,
            myconfig_array,
            mysnap_array,
            mycond_array,
            myplot_map_array,
            SAM_scale_factor_map,
            name_conv_map):
        
        self.histos = {}
        self.myconfig_array = myconfig_array
        self.mysnap_array = mysnap_array
        self.myplot_map_array = myplot_map_array
        self.mycord_array = {}
        self.units_array = {}
        self.mycond_array = mycond_array
        self.SAM_scale_factor_map = SAM_scale_factor_map
        self.name_conv_map = name_conv_map
        
        #print(SAM_scale_factor_map)
                     
    def readData(self,
                 myfilename=False,
                 myhalocat_code=False,
                 preprocessing_only=False):
             
        if myfilename!=False:
            filename=myfilename
        else:
            filename=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)][len(self.myconfig_array['catname'+str(self.a)])+1::]

        print('\n########################################################################################################\n#\n')
        print('#     PROGRESS STATUS: readData() catname:', self.myconfig_array['catname'+str(self.a)], 'self.i:', self.i, 'self.a:', self.a)
        print('\n#\n#######################################################################################################')

            
        data2process = myData.readAnyFormat( data_shape=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_format_info'],
                                             data_format=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format'],
                                             mypath=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_mydir'],
                                             myfilename=filename,
                                             mypath_softlink=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_mysoftlink_dir'],                            
                                             nr_rows=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_nr_rows'],
                                             nr_col=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_nr_col'],                 
                                             skiprow=1,
                                             delim=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_delimiter'],
                                             mydtype_after_read=np.float64,
                                             id_col_array=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array'],
                                             config_array=self.myconfig_array,
                                             name_conv_map=self.name_conv_map,
                                             catname=self.myconfig_array['catname'+str(self.a)],
                                             snapid=str(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_snapid'+str(self.i)]),
                                             value_is_redshift=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_use_snapidz_mapping'],
                                             file_count=self.i,
                                             halocat_code=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_halocat_code'],
                                             halo_id_col_array=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_halo_id_col_array'],
                                             start_fileID=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_start_fileID'],
                                             nr_files_snapshot=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_nr_files_snapshot'],
                                             end_fileID=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_end_fileID'],
                                             scale_factor_map=self.SAM_scale_factor_map)

        self.myData2Process = data2process
        self.output_filename_code_space='' 
        self.tarsel_code_space= ''
        
        #create spaces '_' in the filename:
        if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_output_filename_code']!='':
             self.output_filename_code_space='_'
        if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_tarsel_code']!='':
             self.tarsel_code_space='_'
       
        if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5' or\
           self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format'].find('BINARY')!=-1 or\
           self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='CATASCII' or\
           self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='ROCKSTAR_ASCII' or\
           self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format'].find('EAGLE_ASCII')!=-1 or\
           self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='CHOLLAHDF5' or\
           self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='HYDROHDF5' or\
           self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='MDHDF5':            
                self.correctAndConvertUnits(halocat_code=myhalocat_code,
                                            preprocessing_only=preprocessing_only)
        else:
            try:
                self.volume = myData.cat_attributes['volume']
            except:
                print('volume attribute not exciting! ...')
                try:
                    print('Box size:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'], end='')
                    self.volume=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size']**3
                    print('Volume set!')
                except:
                    print('Box size does not exists!')
                    self.volume=False
        
                                                                                
    def correctAndConvertUnits(self,
                               halocat_code=False,
                               preprocessing_only=False):
        
        #print('here 130 self.mysnap_array:', self.mysnap_array, myData.redshift)
                
        try:
            print('test format of [z] in mysnap_array ...', end=' ')
            print(float(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'])+1.0, end=' ')
            print('is float! CORRECT!')
        except:
            try:
                print('filename correct? --> ', end='')
                self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z']=float(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'])
                print('YES!')
                myData.redshift=float(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'])
                myData.scale_factor=mL.redshift_to_expfactor(myData.redshift)
                print('z+1=', float(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'])+1.0, '--> worked!')
            except:
                try:
                    print('NO! Exclude filename! --> z=', float(self.mysnap_array[self.myconfig_array['catname'+str(self.a)]+'_snapid'+str(self.i)]['z']), end='')
                    myData.redshift=float(self.mysnap_array[self.myconfig_array['catname'+str(self.a)]+'_snapid'+str(self.i)]['z'])
                    myData.scale_factor=mL.redshift_to_expfactor(myData.redshift)
                    print('--> worked!')
                except:
                    try:
                        print('Try with scale factor map ...', end='')
                        myData.redshift=self.SAM_scale_factor_map[self.myconfig_array['catname'+str(self.a)]+'_redshift'+str(self.i)]
                        myData.scale_factor=self.SAM_scale_factor_map[self.myconfig_array['catname'+str(self.a)]+'_scale_factor'+str(self.i)]
                        print('--> worked!')
                    except:
                        print('did not work! Convert [z] in mysnap_array to float or check cat- and filenames!')
            
            
        if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5' or\
           self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='BINARY_SAGE' or\
           self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='HYDROHDF5' or\
           self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='MDHDF5':
        
            #print('here: correctAndConvertUnits() --> set redshift: '
            #print('here 166:', self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'], myData.redshift)
            #self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'] = myData.redshift
        
            try:
                myData.redshift=self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z']
                myData.scale_factor=mL.redshift_to_expfactor(myData.redshift)
            except:
                try:
                    print('--> did not work, except!', end='')
                    myData.redshift=float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_manual_input_redshift'])
                    self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z']=myData.redshift
                    myData.scale_factor=mL.redshift_to_expfactor(myData.redshift)
                    print(myData.redshift+1.0, myData.scale_factor+1.0, 'redshift and scale factor set correctly!!!')
                except:
                    print('--> did not work, try to set scale factor!', end='')
                    myData.scale_factor= float(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['a'])           
                    try:
                        print('---> did not work, second except!', end='')
                        myData.scale_factor=mL.redshift_to_expfactor(myData.redshift)
                        print(myData.redshift+1.0, myData.scale_factor+1.0, 'redshift and scale factor set correctly!!!')   
                    except:
                        print('WARNING REDSHIFT WAS NOT CORRECTLY SET!!!')
       

        print('myData.redshift:', myData.redshift, 'myData.scale_factor:', myData.scale_factor)       
        #try:            
        self.volume = mL.survey_VolumeSqDeg_to_Mpc(box_size=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'], 
                                               redshift=myData.redshift,
                                               h=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_hubble_par'],
                                               scale_factor=myData.scale_factor)
#        except:
#            self.volume=False
                                   

        if halocat_code==False:
            id_col_array_code='_id_col_array'
        else:
            id_col_array_code='_halo_id_col_array' 

 
        try:           
            self.myData2Process = mL.correct_units(self.myconfig_array['catname'+str(self.a)],
                                                   self.myData2Process,
                                                   self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+id_col_array_code],
                                                   self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_correct_little_h'],
                                                   myData.scale_factor)
        except:
            print('WARNING! Property correction to standard values not done!')
            
            
        #print('here: 216 postproc', preprocessing_only)
        if preprocessing_only==True:
            mstar_lower=1e6
            if self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z']<1.0:
                mstar_lower = 1e9
    
            elif self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z']<2.0:
                mstar_lower = 1e8
                
            elif self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z']<3.0:            
                mstar_lower = 1e7
                
            try:
                self.myData2Process = self.myData2Process[np.where((self.myData2Process['mstar_disk']+self.myData2Process['mstar_spheroid'])>mstar_lower)[:][0]]                    
            except:
                self.myData2Process = self.myData2Process[np.where(self.myData2Process['mstar']>mstar_lower)[:][0]]
                
            print('z=', self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'],\
                  ': low mass galaxies with mstar<', mstar_lower, '[Msun] are excluded:', self.myData2Process.shape)


        if self.myconfig_array['catname'+str(self.a)].find('Galacticus')!=-1:                                                            
            self.myData2Process = mL.convert_units_Galacticus(self.myconfig_array['catname'+str(self.a)],
                                                  self.myData2Process,
                                                  self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+id_col_array_code],
                                                  myData.redshift,
                                                  telescope_name=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_telescope_name'],
                                                  apply_k_corr_app=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_apply_k_corr_app'],
                                                  apply_k_corr_abs=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_apply_k_corr_abs'],                                              
                                                  apply_z_boost=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_apply_z_boost'],
                                                  hubble_par=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_hubble_par'],
                                                  unit_code=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_unit_code'],
                                                  use_kcorrect=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_which_kcorrect'])
        else:
            self.myData2Process = mL.convert_units(self.myconfig_array['catname'+str(self.a)],
                                                  self.myData2Process,
                                                  self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+id_col_array_code],
                                                  myData.redshift,
                                                  telescope_name=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_telescope_name'],
                                                  apply_k_corr_app=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_apply_k_corr_app'],
                                                  apply_k_corr_abs=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_apply_k_corr_abs'],                                              
                                                  apply_z_boost=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_apply_z_boost'],
                                                  hubble_par=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_hubble_par'],
                                                  unit_code=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_unit_code'],
                                                  use_kcorrect=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_which_kcorrect'])
   
    def filterAndExtractData(self,
                             myconds_array,
                             data=False,
                             myheader=False,
                             mylog=False,
                             myselection_name='filtered subcatalogue - default settings!',
                             redshift=False,
                             scale_factor=False,
                             subcatsel_code=False,
                             preprocessing_only=False):

                                      
        def check_data(data):
        
            try:
                if data==False:
                    #print('check data=False!'
                    return self.myData2Process
            except:
                return data


        def check_redshift():
            """Check if redshift is set correctly, do so if not"""
            try:
                print('myData.redshift:', myData.redshift)
                redshift=myData.redshift
            except:
                print('Could not find REDSHIFT! --> Taken from SCALE_FACTOR_MAPPING instead!\nz:', end='')
                print(self.SAM_scale_factor_map)
                redshift=self.SAM_scale_factor_map[self.myconfig_array['catname'+str(self.a)]+'_redshift'+str(self.i)]
                myData.redshift=redshift
                print(redshift)
            return redshift
                
        def check_scale_factor():
            """Check if scale_factor is set correctly, do so if not"""              
            try:
                print('myData.scale_factor:', myData.scale_factor)
                scale_factor=myData.scale_factor
            except:
                print('Could not find SCALE FACTOR! --> Taken from SCALE_FACTOR_MAPPING instead!\na:')
                scale_factor=self.SAM_scale_factor_map[self.myconfig_array['catname'+str(self.a)]+'_scale_factor'+str(self.i)]
                myData.scale_factor=scale_factor
                print(scale_factor)
            return scale_factor
        
        
        def filterDensity(myheader, mylog):

            def get_ngal(ndens):
                print('number density:', ndens, 'box size:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'], '[h-1Mpc] little-h:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_hubble_par'], end='')
                if str(ndens).find('h3')!=-1:                    
                    ndens_unit=ndens[ndens.find('h3')::]
                    ndens=ndens[0:ndens.find('h3')]
                print('ndens:', ndens, 'unit:', ndens_unit)
                if ndens_unit.find('h3')!=-1:
                    print('here unit with littel-h absorbed!\nngal=ndens*box^3')
                    ngal=int(float(ndens)*float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'])**3)
                else:
                    print('here unit with littel-h factored out!\nngal=ndens*(box/h)^3')
                    ngal=int(float(ndens)*(float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'])/float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_hubble_par']))**3)
                print('----> ngal:', ngal)
                return ngal           

            d=0
            while d<10:              
                try:
                    ngal = self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filter_density_ngal'+str(d)]
                    print('ngal:', ngal)
                    try:
                        print('ngal: ', int(ngal), end='')
                    except:
                        print('--> not found ...!')
                        ngal=get_ngal(ngal)
                    cut_name = self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filter_density_cut_name'+str(d)]
                except:
                    print('END OF DENISTY FILTERING!')
                    exit()  

                
                c=0
                while c<100:
                    try:
                        name = self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filter_density'+str(c)]
                    except:
                        if c==0:
                            name = self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filter_density'] 
                        else:
                            print('END OF NAME LOOP!\n-------------------\n')
                            break
                  
                    mylog_filter = mylog    
                    print('\nfilter catalog to a certain number density!\n-----------------------------------------------')
                    print('--> which property? ', name, 'ngal:', ngal, 'col id of the sorted column:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array'][name+'_col_id'])
                    print('d:', d, 'c:', c, 'cut_name:', cut_name)
    
                    data_sorted = data[np.argsort(data[name])]
                    
                    mylog_filter+='\nselected for number density property '+name+': '+str(int(ngal)/float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'])**3)+' h3Mpc-3\n'
                    print('data after sort:', data.shape, 'mynumber_density_ngal', int(ngal), end='')
                                    
                    size=len(data_sorted[name])
                    if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filter_density_select_highest']=='True':
                        start_size=size-int(ngal)

                        if start_size<0: start_size=0
                        print('--> select sorted data set: data[', start_size, ':', size, ']')                          
                        data_sorted=data_sorted[start_size:size]
                        
                    else:
                        print('--> select smallest ...')
                        data_sorted=data_sorted[0:int(ngal)]

  
                        
                    print('check size of the filtered data_array:', len(data_sorted[name]), data_sorted.shape)
                    if mycomp.find('z')!=-1:
                        myhdf5_filename='/store/erebos/doris/'+self.myconfig_array['catname'+str(self.a)]+'_z_'+str(float("{0:.2f}".format(redshift)))+'_'+cut_name+'_'+name+self.output_filename_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_output_filename_code']+'.hdf5'                
                    else:
                        myhdf5_filename=mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/'+self.myconfig_array['catname'+str(self.a)]+'_z_'+str(float("{0:.2f}".format(redshift)))+'_'+cut_name+'_'+name+self.output_filename_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_output_filename_code']+'.hdf5'                
                        
                    mylog_filter+='\nngalaxies in the file: '+str(len(data_sorted[name]))+'\n\n'+'link and filename of catalouge: '+myhdf5_filename+'\n\nafter density filtering was applied:\n'                    
        
                    data_sorted, myheader, mydataformat, mylog_filter, check_size, sel_col_list = mL.filter_data_before_write2file(data_sorted,
                                                                                                                            myconds_array,
                                                                                                                            myconds_array,
                                                                                                                            myheader,
                                                                                                                            mylog_filter)

                    if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filter_density_write_coordinates']=='True':
                        writeCoordinates(data_sorted, myhdf5_filename, ngal)
                    if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_OUTPUT_ASCII']!='False':				
                        writeASCII(data_sorted[sel_col_list], myhdf5_filename[0:len(myhdf5_filename)-5]+'.txt', myheader, mydataformat)                                                
                    writeHDF5(data_sorted, myhdf5_filename)
                    writeLOGFILE(mylog_filter, myhdf5_filename[0:len(myhdf5_filename)-5]+'_log.txt')
                    
                    mL.give_galaxy_type_info(data_sorted['orphan'])                       
                    c+=1
                d+=1

            print('END OF FILTERING DENSITY!\n-----------------------------')
            exit()

        def selectRegion(data, myheader, mylog, mylog_region):

            print('filterData() --> selectRegion() ...\n-------------------\n')
                         
            #read file with information about the cutted region
            cut_info_array = myData.readAnyFormat(config=False, mypath=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_path_to_data']+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_filename'], data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=0)
            #print(cut_info_array
            
            cut_id_array={'0': 'x_pos', '0_id_selReg': self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_col_id_x_pos'], 'name0': 'X [h-1Mpc]',
                          '1': 'y_pos', '1_id_selReg': self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_col_id_y_pos'], 'name1': 'Y [h-1Mpc]',
                          '2': 'z_pos', '2_id_selReg': self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_col_id_z_pos'], 'name2': 'Z [h-1Mpc]',
                          }
            myunit='h-1Mpc'              
            if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_volume_type']=='SPHERE':
                size='radius'
            else:
                size='size'

            radius=cut_info_array[0,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_col_id_radius']]
            box_size=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size']
                
            if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_periodic_conds']=='True':
                print('set periodic boundary conditions ...')

                print('box_size', box_size, 'radius:', radius)
                pos=['x_pos','y_pos','z_pos']
               
                i=0
                for k in pos:
                    print('k:', k, 'i:', i, end='')
                    data_new = data[np.where(data[k]<radius)]
                    data_new[k]+=float(box_size)
                    print('len1:', len(data_new), end='')
                    if i==0:
                        data_add=data_new
                    else:
                        data_add = np.concatenate((data_add, data_new), axis=0)                    
                    data_new = data[np.where(data[k]>box_size-radius)]
                    data_new[k]-=float(box_size)
                    print('len2:', len(data_new), end='')
                    data_add = np.concatenate((data_add, data_new), axis=0)
                    i+=1
                    print('data_add:', len(data_add))
            
                data = np.concatenate((data, data_add), axis=0)
                
                print('data after periodic boundary conditions!', data.shape)
            
            mylog_region+='name of region'.ljust(16)+('X ['+myunit+']').ljust(12)+('Y ['+myunit+']').ljust(12)+('Z ['+myunit+']').ljust(12)+(str(size)+' ['+myunit+']').ljust(18)+'total ngal'.ljust(12)+'centrals'.ljust(12)+'satellites'.ljust(12)+'orphans'.ljust(12)+'mhalo'.ljust(10)+'X/Y/Z'.ljust(24)+'MainHaloID'.ljust(20)+'offset: X/Y/Z'.ljust(24)+'\n'
            i=0
            while i<cut_info_array[:,0].size:
                #cut_info_array[i,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_col_id_radius']]=5.0                                             
                mylog_filter=mylog

                print('--------------\ncheck: original data size:', data.shape, 'i:', i, 'region name', end='')

                if i<4:
                    region_name='VoidMDPL_'+str(int(cut_info_array[i,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_region_name']])).zfill(4)
                else:
                    region_name='region'+str(int(cut_info_array[i,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_region_name']])).zfill(4)

                print(region_name)

                if str(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_units'])!=myunit:
                    print('adjusted units to [h-1Mpc] --> ', cut_info_array[i,:], end='')
                    
                    cut_info_array[i,[cut_id_array['0_id_selReg'], cut_id_array['1_id_selReg'], cut_id_array['2_id_selReg'], self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_col_id_radius']]]*=0.6777                    

                mylog_filter+='\nselect a specific subvolume \n---------------------------------------------- \nregion number: '+region_name+', radius ['+myunit+']: '+str(cut_info_array[i,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_col_id_radius']]).zfill(4)+', coordinates: X ['+myunit+']: '+str(cut_info_array[i,cut_id_array['0_id_selReg']])+', Y ['+myunit+']: '+str(cut_info_array[i,cut_id_array['1_id_selReg']])+', Z ['+myunit+']: '+str(cut_info_array[i,cut_id_array['2_id_selReg']])+'\n\n'
 
                if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_volume_type']=='SPHERE':
                    print('Volume: SPHERE! selection coordinate:', cut_info_array[i,cut_id_array['0_id_selReg']],cut_info_array[i,cut_id_array['1_id_selReg']],cut_info_array[i,cut_id_array['2_id_selReg']],'radius: ', cut_info_array[i,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_col_id_radius']])
                   
                    def find_cluster_members(data_pos):
                        a=0
                        while a<3:
                            print('a:', a, cut_id_array[str(a)], cut_info_array[i,cut_id_array[str(a)+'_id_selReg']],  end='')
                            data_pos[cut_id_array[str(a)]]+=(500.0-cut_info_array[i,cut_id_array[str(a)+'_id_selReg']])
    
                            print('ngal>boxsize:', len(data_pos[cut_id_array[str(a)]][np.where(data_pos[cut_id_array[str(a)]]>box_size)]), end='')
                            print('ngal<0.0    :', len(data_pos[cut_id_array[str(a)]][np.where(data_pos[cut_id_array[str(a)]]<0.0)]))
    
                            data_pos[cut_id_array[str(a)]][np.where(data_pos[cut_id_array[str(a)]]>box_size)]-=box_size
                            data_pos[cut_id_array[str(a)]][np.where(data_pos[cut_id_array[str(a)]]<0.0)]+=box_size                       
                            a+=1
                            
                        return np.where(((data_pos['x_pos']-500.0)**2+(data_pos['y_pos']-500.0)**2+(data_pos['z_pos']-500.0)**2)**0.5 < cut_info_array[i,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_col_id_radius']])
                        
                    mask = find_cluster_members(data[['x_pos', 'y_pos', 'z_pos']])
                    
                    data_sorted=data[mask[:][0]]

                print('after selection -->', data_sorted.shape, len(data_sorted))
                try:
                    print(min(data_sorted['x_pos']), '/', max(data_sorted['x_pos']), ',', min(data_sorted['y_pos']), '/', max(data_sorted['y_pos']), ',', min(data_sorted['z_pos']), '/', max(data_sorted['z_pos']), end='')  
                except:
                    pass
                                    

                if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_output_filename_code']=='': 
                    output_space=''
                else:
                    output_space='_'

               
                myfilename=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_output_filename']+\
                           self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_simulation_name']+'_'+\
                           self.myconfig_array['catname'+str(self.a)][0:self.myconfig_array['catname'+str(self.a)].find('_')]+\
                           '_z'+ str(format(float(redshift), '.2f'))+'_'+\
                           region_name+output_space+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_output_filename_code']+'.txt'                           
                
              
                
                #try:
                num_centrals=len(data_sorted['orphan'][np.where(data_sorted['orphan']==0)])                        
                num_sats=len(data_sorted['orphan'][np.where(data_sorted['orphan']==1)])     
                num_orphans=len(data_sorted['orphan'][np.where(data_sorted['orphan']==2)])
                
                print('centrals:', num_centrals, 'sats:', num_sats, 'orphans:', num_orphans, 'total:', num_centrals+num_sats+num_orphans)
                print('data_sorted[[mhalo,x,y,z]]', data_sorted[['mhalo','x_pos', 'y_pos', 'z_pos']])
                max_mhalo= max(data_sorted['mhalo'])
                
                offset_x = cut_info_array[i,cut_id_array['0_id_selReg']] - data_sorted['x_pos'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0]
                offset_y = cut_info_array[i,cut_id_array['1_id_selReg']]- data_sorted['y_pos'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0] 
                offset_z = cut_info_array[i,cut_id_array['2_id_selReg']] - data_sorted['z_pos'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0] 
             
                print('max mhalo:', max_mhalo, 'where(max_mhalo)[0]:', np.where(data_sorted['mhalo']==max_mhalo)[0], end='')
                print('max(data_sorted[[mhalo,x,y,z])', data_sorted['mhalo'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0],\
                                                        data_sorted['x_pos'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0], \
                                                        data_sorted['y_pos'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0],\
                                                        data_sorted['z_pos'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0], \
                                                        data_sorted['haloid'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0])
                mylog_region+=region_name.ljust(16)+str(cut_info_array[i,cut_id_array['0_id_selReg']]).ljust(12)+\
                                                    str(cut_info_array[i,cut_id_array['1_id_selReg']]).ljust(12)+\
                                                    str(cut_info_array[i,cut_id_array['2_id_selReg']]).ljust(12)+\
                                                    str(cut_info_array[i,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_col_id_radius']]).ljust(18)+\
                                                    str(len(data_sorted)).ljust(12)+\
                                                    str(num_centrals).ljust(12)+\
                                                    str(num_sats).ljust(12)+\
                                                    str(num_orphans).ljust(12)+\
                                                    str("{0:.2e}".format(max(data_sorted['mhalo']))).ljust(10)+\
                                                    (str("{0:.2f}".format(data_sorted['x_pos'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0]))+\
                                                     '/'+str("{0:.2f}".format(data_sorted['y_pos'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0]))+'/'+\
                                                    str("{0:.2f}".format(data_sorted['z_pos'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0]))).ljust(24)+\
                                                    (str(data_sorted['haloid'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0])).ljust(20)+\
                                                    (str("{0:.2f}".format(offset_x))+\
                                                     '/'+str("{0:.2f}".format(offset_y))+'/'+\
                                                    str("{0:.2f}".format(offset_z))).ljust(24)+'\n' 
     
                #print('outputfilename:', myfilename
              
                
                data_sorted, header, mydataformat, mylog_filter, check_size, sel_col_list = mL.filter_data_before_write2file(data_sorted,
                                                                                                                        myconds_array,
                                                                                                                        myconds_array,
                                                                                                                        myheader,
                                                                                                                        mylog_filter)                
                                                                                                                        
                mylog_filter+='\nngalaxies in the file: '+str(data_sorted.size)+'\n\n'+'link and filename of catalogue: '+myfilename

                #print(sel_col_list
                #print(mydataformat
            
                writeHDF5(data_sorted, myfilename[0:len(myfilename)-3]+'hdf5')
                writeASCII(data_sorted[sel_col_list], myfilename, header, mydataformat)
                     
                writeLOGFILE(mylog_filter, myfilename[0:len(myfilename)-4]+'_log.txt')

                i+=1               
                #except:
                #    print('no galaxies found in', region_name
            writeLOGFILE(mylog_region, myfilename[0: myfilename.find('reg')]+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_output_filename_code']+output_space+'subsamples_log.txt')
            print('END OF SELECTREGION!\n-----------------------')
            exit()
            
        def writeCoordinates(data,
                             filename,
                             ngal):   
           

            print('Write coordinates in file!')
                 
            myOutput.writeIntoFile(filename[0:len(filename)-5]+'_coordinates.txt',
                                   data[['orphan','x_pos','y_pos','z_pos']],
                                   myheader='',
                                   data_format="%i %0.6e %0.6e %0.6e",
                                   mydelimiter=' ')                         

        def faster_central_pos_for_orphans(data):
            """
            Assigns the centrals spatial positions to all the hosted orphans
        
            """

            central_ind   = (np.where(data['orphan'] == 0))[0] 
            data_centrals = data[central_ind]
            orphan_ind    = (np.where(data['orphan'] == 2))[0]
            data_orphans  = data[orphan_ind]
            if len(orphan_ind) == 0:
                return data
            
            test = np.in1d(data_centrals['haloid'], data_orphans['haloid'])
            centrals=data_centrals[np.where(test==True)]            
            
            print('norphans:', len(data_orphans), 'ncentrals', len(centrals))

            orphan_haloid_unique, orphan_counts = np.unique(data['haloid'][orphan_ind],
                                                            return_counts=True)

            central_haloid_unique, idx_first_haloid_in_central = np.unique(centrals['haloid'],
                                                                              return_index=True)
      
            sort_mask = np.argsort(data['haloid'][orphan_ind])

            orig_idx_haloid_to_orig_array_ind = orphan_ind[sort_mask]

            centrals=centrals[idx_first_haloid_in_central]
                
            curr=0
            for (host_idx, norphans) in zip(np.arange(len(centrals)), orphan_counts):
                dest_sel = np.s_[curr:curr+norphans]

                orphan_indices_this_host = orig_idx_haloid_to_orig_array_ind[dest_sel]

                #print('host_idx:', host_idx, 'haloid centrals:', centrals['haloid'][host_idx], 'haloid_orphan:', data['haloid'][orphan_indices_this_host]
                #print('dest_sel:', dest_sel, 'counts:', norphans, 'orphan_ind this host:', orphan_indices_this_host
                #print('before:', data['x_pos'][orphan_indices_this_host], centrals['x_pos'][host_idx],  
                
                #data['x_pos'][orphan_indices_this_host]=centrals['x_pos'][host_idx]
                for f in ['x_pos', 'y_pos', 'z_pos']:
                    data[f][orphan_indices_this_host] = centrals[f][host_idx]

                #print('--> after:', data['x_pos'][orphan_indices_this_host] 
                #print('-------------'
                #print(' '
                curr += norphans

            return data

        def faster_central_mhalo_for_sats(data, 
                                          id_name='haloid'):
            """
            Assigns the centrals mhalos to all the hosted satellites
            haloid:     ID of the central galaxy/halo
        
            """
          
            central_ind   = (np.where(data['orphan'] == 0))[0] 
            data_centrals = data[central_ind]
            orphan_ind    = (np.where(data['orphan'] > 0))[0]
            data_orphans  = data[orphan_ind]
            if len(orphan_ind) == 0:
                return data
            
            test = np.in1d(data_centrals[id_name], data_orphans[id_name])
            centrals=data_centrals[np.where(test==True)]            
            
            print('norphans:', len(data_orphans), 'ncentrals', len(centrals)) 

            orphan_haloid_unique, orphan_counts = np.unique(data[id_name][orphan_ind],
                                                            return_counts=True)

            central_haloid_unique, idx_first_haloid_in_central = np.unique(centrals[id_name],
                                                                              return_index=True)
      
            sort_mask = np.argsort(data[id_name][orphan_ind])

            orig_idx_haloid_to_orig_array_ind = orphan_ind[sort_mask]

            centrals=centrals[idx_first_haloid_in_central]
                
            curr=0
            for (host_idx, norphans) in zip(np.arange(len(centrals)), orphan_counts):
                dest_sel = np.s_[curr:curr+norphans]

                orphan_indices_this_host = orig_idx_haloid_to_orig_array_ind[dest_sel]

                #print('host_idx:', host_idx, 'haloid centrals:', centrals['haloid'][host_idx], 'haloid_orphan:', data['haloid'][orphan_indices_this_host]
                #print('dest_sel:', dest_sel, 'counts:', norphans, 'orphan_ind this host:', orphan_indices_this_host
                #print('before:', data['x_pos'][orphan_indices_this_host], centrals['x_pos'][host_idx],  
                
                #data['x_pos'][orphan_indices_this_host]=centrals['x_pos'][host_idx]
                for f in ['mhalo_cents']:
                    data[f][orphan_indices_this_host] = centrals[f][host_idx]

#                print('--> after:', data['x_pos'][orphan_indices_this_host] 
#                print('-------------'
#                print(' '
                curr += norphans

            return data

        def writeHDF5(data, filename):                  
            myOutput.write2HDF5(data,
                                myconds_array,
                                self.myconfig_array,
                                filename,
                                self.myconfig_array['catname'+str(self.a)],
                                volume=self.volume,
                                redshift=redshift,
                                scale_factor=scale_factor)

            
        def writeASCII(data, filename, header, dataformat):
            myOutput.writeIntoFile(filename,
                                   data,
                                   myheader=header,
                                   data_format=dataformat)

        def writeLOGFILE(log, filename):            
            myOutput.writeIntoFile(filename,
                                   log,
                                   myheader=myheader,
                                   data_format="%s",
                                   append_mytext=False,
                                   data_is_string=True)

        def units2log():            
            for item in self.name_conv_map:
                if self.mycond_array[item+'_unit'].find('Msun')!=-1 or\
                    self.mycond_array[item+'_unit'].startswith('yr-'):
                    data[item]=np.log10(data[item])
                    print(item, 'unit:', self.mycond_array[item+'_unit'],'--> log10')
                    unit_old=self.mycond_array[item+'_unit']
                    self.mycond_array[item+'_unit']='log10('+unit_old+')'




        print('\n########################################################################################################\n#\n')  
        print('#     PROGRESS STATUS: filterAndExtractData() catname:', self.myconfig_array['catname'+str(self.a)], 'self.i:', self.i, 'self.a:', self.a)
        print('\n#\n########################################################################################################')


        data = check_data(data)
        redshift=check_redshift()
        scale_factor=check_scale_factor()

        i=0
        count_Mzstar=0
        count_Mzgas=0
        count_mstar=0
        count_mstar_IC=0
        count_mhalo=0
        count_sfr=0
        count_mag=0
           
        while i<self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['nr_entries']:
                     
            if self.myconfig_array['catname'+str(self.a)].startswith('Galacticus') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)]=='hostid' and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5':
                #print('GALACTICUS correct HOSTID of satellites and orphans!'                
                mask = np.where(data['orphan']!=0)               
                data['hostid'][mask[:][0]]=data['satelliteNodeIndex'][mask[:][0]]

            if self.myconfig_array['catname'+str(self.a)].startswith('Galacticus') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)]=='haloid' and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5':
                #print('GALACTICUS: correct HALOID of satellites and orphans!'                    
                mask = np.where(data['orphan']!=0)
                data['haloid'][mask[:][0]]=data['parentIndex'][mask[:][0]]

            if self.myconfig_array['catname'+str(self.a)].startswith('Galacticus') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)]=='mhalo' and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5' and count_mhalo==0:
                #print('GALACTICUS: set mhalo of the satellites to sattelliteBoundMass aka mhalo_sat!!!!'                
                mask = np.where(data['orphan']!=0)
                data['mhalo'][mask[:][0]]=data['mhalo_sat'][mask[:][0]]
                mask = np.where(data['orphan']!=0)              

            if self.myconfig_array['catname'+str(self.a)].startswith('Galacticus') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)]=='Mzgas' and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5' and count_Mzgas==0:
                #print('GALACTICUS: MZgas = MZgas_spheroid + MZgas_disk!!!!'
                count_Mzgas+=1
                data['Mzgas']=data['Mzgas_spheroid']+data['Mzgas_disk']

            if self.myconfig_array['catname'+str(self.a)].startswith('Galacticus') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)]=='Mzstar' and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5' and count_Mzstar==0:
                #print('GALACTICUS: MZstar = MZstar_spheroid + MZstar_disk!!!!'
                count_Mzstar+=1
                data['Mzstar']=data['Mzstar_spheroid']+data['Mzstar_disk']



            if self.myconfig_array['catname'+str(self.a)].startswith('SAG_') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)]=='sfr_disk' and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5' and count_sfr==0:
                print('SAG: set sfr_disk!!!!')
                count_sfr+=1
                data['sfr_disk']-=data['sfr_spheroid']
            
            if self.myconfig_array['catname'+str(self.a)].startswith('SAGE') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)]=='mstar_disk' and (self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='BINARY_SAGE' or self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5')  and count_mstar==0:
                print('SAGE: set mstar_disk!!!!')
                count_mstar+=1
                data['mstar_disk']-=data['mstar_spheroid']
                
            if self.myconfig_array['catname'+str(self.a)].startswith('SAGE') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)]=='mstar_disk' and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='BINARY_SAGE' and count_mstar_IC==0:
                count_mstar_IC+=1
                try:
                    data['mstar+IC']+=data['mstar_IC']
                    print('SAGE: set total mstar! Mstar+Mstar_IC!!!!')
                except:
                    pass
            
            if self.myconfig_array['catname'+str(self.a)].startswith('SAGE') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)]=='Mzstar_disk' and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='BINARY_SAGE' and count_Mzstar==0:
                print('SAGE: set metals mstar_disk!!!!' )
                data['Mzstar_disk']-=data['Mzstar_spheroid']
                count_Mzstar+=1

            # if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)].find('AB')!=-1 and count_mag==0:
            #     if self.myconfig_array['catname'+str(self.a)].find('Galacticus')!=-1:
            #         data = mL.convert_units_Galacticus(
            #                               self.myconfig_array['catname'+str(self.a)],
            #                               data,
            #                               self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array'],
            #                               redshift,
            #                               telescope_name=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_telescope_name'],
            #                               apply_k_corr_app=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_apply_k_corr_app'],
            #                               apply_k_corr_abs=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_apply_k_corr_abs'],                                              
            #                               apply_z_boost=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_apply_z_boost'],
            #                               hubble_par=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_hubble_par'],
            #                               unit_code=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_unit_code'],
            #                               use_kcorrect=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_which_kcorrect'],
            #                               cosmology=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_cosmology'])
            #     else:
            #         data = mL.convert_units(
            #                               self.myconfig_array['catname'+str(self.a)],
            #                               data,
            #                               self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array'],
            #                               redshift,
            #                               telescope_name=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_telescope_name'],
            #                               apply_k_corr_app=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_apply_k_corr_app'],
            #                               apply_k_corr_abs=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_apply_k_corr_abs'],                                              
            #                               apply_z_boost=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_apply_z_boost'],
            #                               hubble_par=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_hubble_par'],
            #                               unit_code=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_unit_code'],
            #                               use_kcorrect=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_which_kcorrect'],
            #                               cosmology=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_cosmology'])
            #     count_mag+=1
                

            # if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)].find('cut')!=-1:       
            #     data = mL.calc_colour_cut_parameter(data, 
            #                                         self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array'],
            #                                         self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)])
            i+=1

            
        if myheader==False:
            myheader=self.myconfig_array['catname'+str(self.a)]+' z='+str(float("{0:.2f}".format(myData.redshift)))+'\n'

        #print(self.myconfig_array
            
        if mylog==False:
            mylog='Captains log-file: Stardate '+str(date_time_stamp)+'\nname: '+self.myconfig_array['catname'+str(self.a)]+\
                '\nredshift: '+str(float("{0:.2f}".format(myData.redshift)))+\
                '\nsimulation name: '+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_simulation_name']+\
                '\ncosmology: '+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_cosmology']+\
                '\nnote: '+str(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_my_annotation'])+\
                '\n++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n'
        else:
            mylog+='\n'+myselection_name+'\n\n'

        mylog_region=mylog

        if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_set_register_name']=='Skies':
            filename_part2='SkiesANDUniverses/MDPL2_'\
                            +self.myconfig_array['catname'+str(self.a)][0:len(self.myconfig_array['catname'+str(self.a)])-5] \
                            +'_z_'+str(format(float(redshift), '.2f'))+'.hdf5' 
            myhdf5_filename='/store/erebos/doris/'+filename_part2
        elif self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_set_register_name']=='MD-Galaxies':
            filename_part2='MultiDark-Galaxies/MDPL2_'\
                            +self.myconfig_array['catname'+str(self.a)][0:len(self.myconfig_array['catname'+str(self.a)])-5] \
                            +'_z_'+str(format(float(redshift), '.2f'))+'.hdf5' 
            myhdf5_filename='/store/erebos/doris/'+filename_part2
        elif self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_set_register_name']=='OII':
            filename_part2='SAG_OII/'\
                            +self.myconfig_array['catname'+str(self.a)]\
                            +'_z_'+str(float("{0:.2f}".format(redshift)))+'_tarsel'+self.tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]\
                            +'_tarsel_code']+self.output_filename_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_output_filename_code']\
                            +'.hdf5'
            myhdf5_filename='/store/erebos/doris/'+filename_part2
        elif self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_set_register_name']=='CMASS':
            filename_part2=self.myconfig_array['catname'+str(self.a)]\
                            +self.tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]\
                            +'_tarsel_code']+self.output_filename_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_output_filename_code']\
                            +'.hdf5'
            myhdf5_filename='/store/erebos/doris/'+filename_part2              
        else:
            filename_part2=self.myconfig_array['catname'+str(self.a)]\
                            +'_z_'+str(float("{0:.2f}".format(redshift)))+'_tarsel'+self.tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]\
                            +'_tarsel_code']+self.output_filename_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_output_filename_code']\
                            +'.hdf5'
            try:
                myhdf5_filename=mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/'+filename_part2
            except:
                myhdf5_filename='/store/erebos/doris/'+filename_part2

        if preprocessing_only==False: 
            if self.myconfig_array['catname'+str(self.a)].startswith('Galacticus') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5':      
                try:
                    data=faster_central_pos_for_orphans(data)
                    print('Correct Galacitucs orphans positions!!!! --> DONE!')
                except:
                    pass
    
            if self.myconfig_array['catname'+str(self.a)].startswith('Galacticus') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5':
                try:
                    data=faster_central_mhalo_for_sats(data)
                    print('Set central mhalo (HOD calculation)!!!! --> DONE!')
                except:
                    pass               



        if self.myconfig_array['catname'+str(self.a)].startswith('EAGLE'):
            data = mL.handle_HF_EAGLE_data(data, key='crossMatch')
            print('Line 1037: after crossmatching', data.size)
            #data = mL.handle_HF_EAGLE_data(data, key='print2Table')
            #data = mL.handle_HF_EAGLE_data(data, key='filterData')
            #data['ropt'] = mL.handle_HF_EAGLE_data(data, key='correctRopt')
            #print('Line 1041: correct radii', data.size)

        if self.myconfig_array['catname'+str(self.a)].startswith('Galacticus'):
            #300
           # data = mL.handle_Galacticus_300_data(data, key='calcClusterProps')
            data = mL.handle_Galacticus_300_data(data, key='crossMatch')
            
        #exit()
        calc_props=False
        
        if calc_props==True:
            try:
                data['mstar'] = data['mstar_disk']+data['mstar_spheroid']
                #print('mstar=mstar_disk+mstar_spheroid'
            except:
                 pass
            try:
                data['sfr'] = data['sfr_disk']+data['sfr_spheroid']
                #print('sfr=sfr_disk+sfr_spheroid'            
            except:
                pass            
            try:
                data['mcold'] = data['mcold_disk']+data['mcold_spheroid']
                #print('mcold=mcold_disk+mcold_spheroid'            
            except:
                pass
            try:
                data['Mzstar'] = data['Mzstar_disk']+data['Mzstar_spheroid']
                #print('Mzstar=Mzstar_disk+Mzstar_spheroid'            
            except:
                pass        
            try:
                data['Mzgas'] = data['Mzgas_disk']+data['Mzgas_spheroid']
                #print('Mzgas=Mzgas_disk+Mzgas_spheroid'            
            except:
                pass
            try:
                data['MzgasvsMzstar'] = data['Mzgas']/data['Mzstar']
                #print('Mzstar=Mzstar_disk+Mzstar_spheroid'            
            except:
                pass          
            try:
                data['ssfr'] = data['sfr']/data['mstar']
                print('ssfr=sfr/mstar\n min/max:', min(data['ssfr']), '/', max(data['ssfr']))            
            except:
                pass 
            try:
                data['ssfr_30kpc'] = data['sfr_30kpc']/data['mstar_30kpc']
                print('ssfr_30kpc=sfr_30kpc/mstar_30kpc\n min/max:', min(data['ssfr_30kpc']), '/', max(data['ssfr_30kpc']), '\n')           
            except:
                pass
            try:
                data['ssfr_gasSF'] = data['sfr_gasSF']/data['mPart_gasSF']
                print('ssfr_gasSF=sfr_gasSF/mPart_gasSF\n min/max:', min(data['ssfr_gasSF']), '/', max(data['ssfr_gasSF']), '\n')            
            except:
                pass        
            
            try:
                data['zcold']=8.69+np.log10(data['Mzgas']/(data['mcold']*0.0134))
                print('zcold=8.69+log10(Mzgas/Mcold/0.0134)\n')
            except:
                pass
            
            
            try:
                data['zcold_gasSF']=12.0+np.log10(data['zgasSF_O']/data['zgasSF_H']/16.0)
                print('zcold_gasSF=12+log10(O/H)\n')
            except:
                pass      
            
            try:
                data['zstar']=8.69+np.log10(data['Mzstar']/(data['mstar']*0.0134))
                print('zstar=8.69+log10(Mzstar/mstar/0.0134)\n')            
            except:
                pass
            try:
                data['zcold_zstar']=data['zcold']-data['zstar']
                print('zcold-zstar\n')           
            except:
                pass 
    
            try:
                data['OH_gas_30kpc'] = 12+np.log10(data['zgasSF_O']/data['zgasSF_H']/16.0)
            except:
                pass
        
            try:
                data['Tcons']=data['mcold']/data['sfr']/1e9
                print('gas depletion time: mcold/sfr/1e9 [Gyr]\n')         
            except:
                pass

            try:
                data['Tcons_30kpc']=data['mgas_30kpc']/data['sfr_30kpc']/1e9
                print('gas depletion time 30kpc: Tcons_30kpc=mgas_30kpc/sfr_30kpc/1e9 [Gyr]\n')         
            except:
                pass
            try:
                data['Tcons_1.5ropt']=data['mgas_1.5ropt']/data['sfr_1.5ropt']/1e9
                print('gas depletion time: Tcons_1.5ropt=mgas_1.5ropt/sfr_1.5ropt/1e9 [Gyr]\n')       
            except:
                pass
            try:
                data['Tcons_HI_30kpc']=data['Mgas_HI_30kpc_GK11']/data['sfr_30kpc']/1e9
                print('gas depletion time Tcons_HIH2_30kpc_GK11: Mgas_H2_30kpc_GK11/sfr_30kpc/1e9 [Gyr]\n')        
            except:
                pass
            try:
                data['Tcons_HIH2_30kpc']=(data['Mgas_HI_30kpc_GK11']+data['Mgas_H2_30kpc_GK11'])/data['sfr_30kpc']/1e9
                print('gas depletion time Tcons_HIH2_30kpc_GK11: (Mgas_H2_30kpc_GK11+Mgas_H2_30kpc_GK11)/sfr_30kpc/1e9 [Gyr]\n')        
            except:
                pass
            try:
                data['Tcons_gasSF']=data['mPart_gasSF']/data['sfr_gasSF']/1e9
                print('gas depletion time: Tcons_gasSF=mPart_gasSF/sfr_gasSF/1e9 [Gyr]\n')         
            except:
                pass
            
            try:
                data['fbar']=data['mcold']/(data['mcold']+data['mstar'])
                print('baryon fraction fbar=mcold/(mcold+mstar)\n')       
            except:
                pass
            
            try:
                data['fbar_30kpc']=data['mgas_30kpc']/(data['mgas_30kpc']+data['mstar_30kpc'])
                print('baryon fraction fbar_30kpc: mgas_30kpc/(mgas_30kpc+mstar_30kpc)\n')       
            except:
                pass
            try:
                data['fbar_1.5ropt']=data['mgas_1.5ropt']/(data['mgas_1.5ropt']+data['mstar_1.5ropt'])
                print('baryon fraction fbar_1.5ropt: mgas_1.5ropt/(mgas_1.5ropt+mstar_1.5ropt)\n')       
            except:
                pass
            try:
                data['fbar_HI_30kpc']=data['Mgas_HI_30kpc_GK11']/(data['Mgas_HI_30kpc_GK11']+data['Mgas_H2_30kpc_GK11']+data['mstar_30kpc'])
                print('baryon fraction fbar_HIH2_30kpc: Mgas_HI_30kpc_GK11/(Mgas_HI_30kpc_GK11+Mgas_H2_30kpc_GK11+mstar_30kpc)\n')       
            except:
                pass
            try:
                data['fbar_HIH2_30kpc']=(data['Mgas_HI_30kpc_GK11']+data['Mgas_H2_30kpc_GK11'])/(data['Mgas_HI_30kpc_GK11']+data['Mgas_H2_30kpc_GK11']+data['mstar_30kpc'])
                print('baryon fraction fbar_HIH2_30kpc: (Mgas_HI_30kpc_GK11+Mgas_H2_30kpc_GK11)/(Mgas_HI_30kpc_GK11+Mgas_H2_30kpc_GK11+mstar_30kpc)\n')       
            except:
                pass        
            try:
                data['fbar_gasSF']=data['mPart_gasSF']/(data['mPart_gasSF']+data['mPart_stars'])
                print('baryon fraction fbar_gasSF: mPart_gasSF/(mPart_gasSF+mPart_star)\n')       
            except:
                pass
    
            try:
                #atomic gas fraction from Rosas-Guevara+22 (page 5) and Obreschkow+16  
                data['fatom_30kpc_GK11']=1.35*data['Mgas_HI_30kpc_GK11']/(data['mstar_30kpc']+1.35*(data['Mgas_HI_30kpc_GK11']+data['Mgas_H2_30kpc_GK11']))
                print('atomic gas fraction fatmo: 1.35*Mgas_HI_30kpc_GK11/(mstar+ 1.35(Mgas_HI_30kpc_GK11+Mgas_H2_30kpc_GK11))\n')       
            except:
                pass
            try:
                #atomic gas fraction from Rosas-Guevara+22 (page 5) and Obreschkow+16  
                data['fmol_30kpc_GK11']=1.35*data['Mgas_H2_30kpc_GK11']/(data['mstar_30kpc']+1.35*(data['Mgas_HI_30kpc_GK11']+data['Mgas_H2_30kpc_GK11']))
                print('molecular gas fraction fmol: 1.35*Mgas_H2_30kpc_GK11/(mstar+ 1.35(Mgas_HI_30kpc_GK11+Mgas_H2_30kpc_GK11))\n')       
            except:
                pass        
    
            try:
                data['cgf']=data['mcold']/data['mstar']
                print('cold gas fraction cgf: mcold/mstar\n')       
            except:
                pass

            try:
                data['cgf_30kpc']=data['mgas_30kpc']/data['mstar_30kpc']
                print('cold gas fraction within 30kpc cgf_30kpc: mgas_30kpc/mstar_30kpc\n')       
            except:
                pass
            try:
                data['cgf_1.5ropt']=data['mgas_1.5ropt']/data['mstar_1.5ropt']
                print('cold gas fraction within 1.5ropt cgf_1.5ropt: mgas_1.5ropt/mstar_1.5ropt\n')       
            except:
                pass
            try:
                data['cgf_HI_30kpc']=data['Mgas_HI_30kpc_GK11']/data['mstar_30kpc']
                print('cold gas fraction of HI gas within 30kpc cgf_HI_30kpc: Mgas_HI_30kpc_GK11/mstar_30kpc\n')       
            except:
                pass          
            try:
                data['cgf_HIH2_30kpc']=(data['Mgas_HI_30kpc_GK11']+data['Mgas_H2_30kpc_GK11'])/data['mstar_30kpc']
                print('cold gas fraction of HI+H2 gas within 30kpc cgf_HIH2_30kpc: (Mgas_HI_30kpc_GK11+Mgas_H2_30kpc_GK11)/mstar_30kpc\n')       
            except:
                pass        
            try:
                data['cgf_gasSF']=data['mPart_gasSF']/data['mPart_stars']
                print('cold gas fraction of SF-particles cgf_gasSF: mPart_gasSF/mPart_stars\n')       
            except:
                pass       
     
            try:
                data['Sigma_gas_reff_1.5ropt']=(data['mgas_1.5ropt']/2.0)/(2*np.pi*(data['reff_gas_1.5ropt']*1000)**2)
                print('Sigma_gas_reff_1.5ropt=mgas_half_reff_gas_disk/2*Pi*(reff_gas_1.5ropt*1000)**2 [Msunpc-2]\n')
            except:
                pass       
     
            try:
                data['Sigma_gas_reff_disk_1.5ropt']=data['mgas_half_reff_gas_disk']/(2*np.pi*(data['reff_gas_disk_1.5ropt']*1000)**2)
                print('Sigma_gas_reff_disk_1.5roptt=mgas_half_reff_gas_disk/2*Pi*(reff_gas_disk_1.5ropt*1000)**2 [Msunpc-2]\n')
            except:
                pass
            try:
                data['Sigma_gas_1.5ropt']=data['mgas_1.5ropt']/(2*np.pi*(1.5*data['ropt']*1000)**2)
                print('Sigma_gas_1.5ropt=mgas_1.5ropt/2*Pi*(1.5ropt*1000)**2 [Msunpc-2]\n')
            except:
                pass
            try:
                data['Sigma_sfr_1.5ropt']=data['sfr_1.5ropt']/(2*np.pi*(1.5*data['ropt'])**2)
                print('Sigma_sfr_1.5ropt=sfr_1.5ropt/2*Pi*(1.5ropt)**2 [Msunyr-1kpc-2]\n')
            except:
                pass
            try:
                data['Sigma_sfr_30kpc']=data['sfr_30kpc']/(2*np.pi*(data['rhalf_stars_30kpc'])**2)
                print('Sigma_sfr_30kpc=(sfr_30kpc)/2*Pi*(rhalf_stars_30kpc)**2 [Msunyr-1kpc-2]\n')
            except:
                pass
            try:
                data['Sigma_HI_30kpc']=data['Mgas_HI_30kpc_GK11']/2.0/(2*np.pi*(data['rhalf_stars_30kpc']*1000)**2)
                print('Sigma_HI_30kpc=Mgas_HI_30kpc_GK11/2*Pi*(rhalf_stars_30kpc_2D*1000)**2 [Msunpc-2]\n')
            except:
                pass
            try:
                data['Sigma_H2_30kpc']=data['Mgas_H2_30kpc_GK11']/2.0/(2*np.pi*(data['rhalf_stars_30kpc']*1000)**2)
                print('Sigma_H2_30kpc=Mgas_H2_30kpc_GK11/2*Pi*(rhalf_stars_30kpc_2D*1000)**2 [Msunpc-2]\n')
            except:
                pass        
            try:
                 data['Sigma_HIH2_30kpc']=(data['Mgas_HI_30kpc_GK11']+data['Mgas_H2_30kpc_GK11'])/2.0/(2*np.pi*(data['rhalf_stars_30kpc']*1000)**2)
                 print('Sigma_HIH2_30kpc=(Mgas_H2_30kpc_GK11+Mgas_H2_30kpc_GK11)/2*Pi*(rhalf_stars_30kpc_2D*1000)**2 [Msunpc-2]\n')
            except:
                 pass       
            try:
                data['Sigma_stars_30kpc']=(data['mstar_30kpc']/2.0)/(2*np.pi*(data['rhalf_stars_30kpc'])**2)
                print('Sigma_stars_30kpc=(mstar_30kpc/2)/2*Pi*(rhalf_stars_30kpc_2D)**2 [Msunkpc-2]\n')
            except:
                pass
            try:
                data['Sigma_stars_1.5ropt']=(data['mstar_1.5ropt']/2.0)/(2*np.pi*(data['rhalf_stars_1.5ropt'])**2)
                print('Sigma_star_1.5ropt=mstar_1.5ropt/2.0/2*Pi*(rhalf_star_1.5ropt)**2 [Msunkpc-2]\n')
            except:
                pass
            try:
                data['Dgas_disk_1.5ropt']=2.0*data['reff_gas_disk_1.5ropt']
            except:
                pass       
          
            try:
                data['BvT']=data['mstar_spheroid']/data['mstar']
                print('B/T=mstar_spheroid/mstar --> bulge mass to total stellar mass')
            except:
                pass
            try:
                data['rbulgevsrdisk']=data['rbulge']/data['rdisk']
                print('rbulgevsrdsik=rbulge/rdisk --> bulge to disk radius')
            except:
                pass      
            try:
                data['SHMR']=data['mstar']/data['mhalo']
                print('SHMR=mstar/mhalo --> Galaxy halo relation')
            except:
                pass
            try:
                data['SHMR']=data['mstar']/data['mhalo_200c']
                print('SHMR=mstar/mhalo --> Galaxy halo relation')
            except:
                pass 
            try:
                data['SHMR']=data['mPart_stars']/data['mhalo_200c']
                print('SHMR=mPart_stars/mhalo --> Galaxy halo relation')
            except:
                pass       
            try:
                data['SHMR_30kpc']=data['mstar_30kpc']/data['mhalo_200c']
                print('SHMR_30kpc=mstar_30kpc/mhalo_200c --> Galaxy halo relation')
            except:
                pass
            try:
                data['SHMR_1.5ropt']=data['mstar_1.5ropt']/data['mhalo_200c']
                print('SHMR=mstar_1.5ropt/mhalo_200c --> Galaxy halo relation\n')
            except:
                pass        
            
            try:
                data['sfe_HI_30kpc']=data['Sigma_sfr_30kpc']/(data['Sigma_HI_30kpc']*1000.0**2)
                print('total HI star formation efficiency sfe_HI_30kpc: Sigma_sfr_30kpc/(Sigma_HI_30kpc*1000^2) [yr-1]\n')
            except:
                pass
    
            try:
                data['sfe_HIH2_30kpc']=data['Sigma_sfr_30kpc']/(data['Sigma_HIH2_30kpc']*1000.0**2)
                print('star formation efficiency of total HI+H2  sfe_HIH2_30kpc: Sigma_sfr_30kpc/(Sigma_HIH2_30kpc*1000^2) [yr-1]\n')
            except:
                pass        
     
            try:
                data['sfe_gas_reff_1.5ropt']=data['Sigma_sfr_1.5ropt']/(data['Sigma_gas_reff_1.5ropt']*1000.0**2)
                print('star formation efficiency within gas disk and 1.5ropt sfe_gas_reff_1.5ropt=Sigma_sfr_1.5ropt/(Sigma_gas_reff_1.5ropt*1000^2) [yr-1]\n')
            except:
                pass
            try:
                data['sfe_gas_reff_disk_1.5ropt']=data['Sigma_sfr_1.5ropt']/(data['Sigma_gas_reff_disk_1.5ropt']*1000.0**2)
                print('star formation efficiency within gas disk and 1.5ropt sfe_gas_reff_disk_1.5ropt=Sigma_sfr_1.5ropt/(Sigma_gas_disk_1.5ropt*1000^2) [yr-1]\n')
            except:
                pass 
     
            try:
                data['sfe_gas_1.5ropt']=data['Sigma_sfr_1.5ropt']/(data['Sigma_gas_1.5ropt']*1000.0**2)
                print('star formation efficiency of the gas within and 1.5ropt sfe_gas_1.5ropt=Sigma_sfr_1.5ropt/(Sigma_gas_1.5ropt*1000^2) [yr-1]\n')
            except:
                pass            

            
            for filter_band in ['B', 'V']:
                try:
                    data['MAB_dA_Johnson_'+filter_band] = mL.convert_SDSS_to_Johnson_Lupton2005(filter_band,
                                                                                            data['MAB_dA_SDSS_g'], 
                                                                                            data['MAB_dA_SDSS_r'])
                    print('filter_band:', filter_band, 'successfully converted using Lupton (2005)!')
                except:
                    print('filter_band:', filter_band, 'Conversion FAILED!')  
            
            try: 
                radius = mL.conv_radius_pc_to_arcsec(data['ropt']*1000.0)
                data['SB_mu_ropt_B']=data['MAB_dA_Johnson_B']+ 2.5*np.log10(2.0*np.pi*(radius)**2.0)
                #data = data[np.where(data['MAB_dA_Johnson_B']!=0.0)[:][0]]

                print('optical surface brigthness Johnson B-band within Ropt [mag arcsec-2]\n')

            except:
                pass
            try: 
                radius = mL.conv_radius_pc_to_arcsec(1.5*data['ropt']*1000.0)
                data['SB_mu_1.5ropt_B']=data['MAB_dA_Johnson_B']+ 2.5*np.log10(2.0*np.pi*(radius)**2.0)
                #data = data[np.where(data['MAB_dA_Johnson_B']!=0.0)[:][0]]
                print('optical surface brigthness Johnson B-band within 1.5Ropt [mag arcsec-2]\n')
                print('Line 1386: after MB!=0', data.size)
            except:
                pass
            try: 
                radius = mL.conv_radius_pc_to_arcsec(data['ropt']*1000.0)
                data['SB_mu_corr_ropt_B']=data['MAB_dA_Johnson_B']+ 2.5*np.log10(2.0*np.pi*(radius)**2.0)
                data = data[np.where(data['MAB_dA_Johnson_B']!=0.0)[:][0]]
                print('corrected optical surface brigthness Johnson B-band within the corrected Ropt [mag arcsec-2]\n')
                print('Line 1386: after MB!=0', data.size)
            except:
                pass            
            try:
                radius = mL.conv_radius_pc_to_arcsec(data['reff_gas_disk_1.5ropt']*1000.0)
                data['SB_mu_eff_gas_disk_B']=data['MAB_dA_Johnson_B']+2.5*np.log10(2.0*np.pi*(radius)**2.0)
                print('effective surface brigthness SDSS B-band within gas disk radius [mag arcsec-2]\n')
            except:
                pass
            try:
                radius = mL.conv_radius_pc_to_arcsec(data['rhalf_stars_1.5ropt']*1000.0)
                data['SB_mu_eff_stars_1.5ropt_r']=data['MAB_dA_SDSS_r']+2.5*np.log10(2.0*np.pi*(radius)**2.0)
                #data = data[np.where(data['MAB_dA_SDSS_r']!=0.0)[:][0]]
                print('stellar surface brigthness SDSS r-band within half mass radius stars [mag arcsec-2]\n')
            except:
                pass        
            try:
                radius = mL.conv_radius_pc_to_arcsec(data['rhalf_stars_30kpc']*1000.0)
                data['SB_mu_eff_stars_30kpc_r']=data['MAB_dA_SDSS_r']+2.5*np.log10(2.0*np.pi*(radius)**2.0)
                #data = data[np.where(data['MAB_dA_SDSS_r']!=0.0)[:][0]]
                print('stellar surface brigthness r-band within 30kpc aperture [mag arcsec-2]\n')
            except:
                pass         
    
            try: 
                radius = mL.conv_radius_pc_to_arcsec(data['ropt']*1000.0)
                data['SB_mu_ropt_r']=data['MAB_dA_SDSS_r']+ 2.5*np.log10(2.0*np.pi*(radius)**2.0)
                data = data[np.where(data['MAB_dA_SDSS_r']!=0.0)[:][0]]
                print('optical surface brigthness SDSS r-band within Ropt [mag arcsec-2]\n')
    
            except:
                pass
            try:
                radius = mL.conv_radius_pc_to_arcsec(data['reff_gas_disk_1.5ropt']*1000.0)
                data['SB_mu_eff_gas_disk_r']=data['MAB_dA_SDSS_r']+2.5*np.log10(2.0*np.pi*(radius)**2.0)
                #data = data[np.where(data['MAB_dA_SDSS_r']!=0.0)[:][0]]
                print('effective surface brigthness SDSS r-band within gas disk radius [mag arcsec-2]\n')
            except:
                pass
            try:
                radius = mL.conv_radius_pc_to_arcsec(data['rhalf_stars_1.5ropt']*1000.0)
                data['SB_mu_eff_stars_1.5ropt_B']=data['MAB_dA_Johnson_B']+2.5*np.log10(2.0*np.pi*(radius)**2.0)
                print('stellar surface brigthness Johnson B-band within half mass radius stars [mag arcsec-2]\n')
            except:
                pass        
            try:
                radius = mL.conv_radius_pc_to_arcsec(data['rhalf_stars_30kpc']*1000.0)
                data['SB_mu_eff_stars_30kpc_B']=data['MAB_dA_Johnson_B']+2.5*np.log10(2.0*np.pi*(radius)**2.0)
                #data = data[np.where(data['MAB_dA_Johnson_B']!=0.0)[:][0]]
                print('stellar surface brigthness Johnson B-band within 30kpc aperture [mag arcsec-2]\n')
            except:
                pass 
    
            try: 
                radius = mL.conv_radius_pc_to_arcsec(data['ropt']*1000.0)
                data['SB_mu_ropt_V']=data['MAB_dA_Johnson_V']+ 2.5*np.log10(2.0*np.pi*(radius)**2.0)
                #data = data[np.where(data['MAB_dA_Johnson_V']!=0.0)[:][0]]
                print('optical surface brigthness Johnson V-band within Ropt [mag arcsec-2]\n')
            except:
                pass
                 
            try:
                radius = mL.conv_radius_pc_to_arcsec(data['reff_gas_disk_1.5ropt']*1000.0)
                data['SB_mu_eff_gas_disk_V']=data['MAB_dA_Johnson_V']+2.5*np.log10(2.0*np.pi*(radius)**2.0)
                #data = data[np.where(data['MAB_dA_Johnson_V']!=0.0)[:][0]]
                print('effective surface brigthness Johnson V-band within gas disk radius [mag arcsec-2]\n')
            except:
                pass
            try:
                radius = mL.conv_radius_pc_to_arcsec(data['rhalf_stars_1.5ropt']*1000.0)
                data['SB_mu_eff_stars_1.5ropt_V']=data['MAB_dA_Johnson_V']+2.5*np.log10(2.0*np.pi*(radius)**2.0)
                #data = data[np.where(data['MAB_dA_Johnson_V']!=0.0)[:][0]]
                print('stellar surface brigthness Johnson V-band within half mass radius stars [mag arcsec-2]\n')
            except:
                pass        
            try:
                radius = mL.conv_radius_pc_to_arcsec(data['rhalf_stars_30kpc']*1000.0)
                data['SB_mu_eff_stars_30kpc_V']=data['MAB_dA_Johnson_V']+2.5*np.log10(2.0*np.pi*(radius)**2.0)
                #data = data[np.where(data['MAB_dA_Johnson_V']!=0.0)[:][0]]
                print('stellar surface brigthness Johnson V-band within 30kpc aperture [mag arcsec-2]\n')
            except:
                pass
            try: 
                radius = mL.conv_radius_pc_to_arcsec(data['ropt']*1000.0)
                data['SB_mu_ropt_g']=data['MAB_dA_SDSS_g']+ 2.5*np.log10(2.0*np.pi*(radius)**2.0)
                #data = data[np.where(data['MAB_dA_SDSS_g']!=0.0)[:][0]]
                print('optical surface brigthness SDSS g-band within Ropt [mag arcsec-2]\n')
    
            except:
                pass
            try:
                radius = mL.conv_radius_pc_to_arcsec(data['reff_gas_disk_1.5ropt']*1000.0)
                data['SB_mu_eff_gas_disk_g']=data['MAB_dA_SDSS_g']+2.5*np.log10(2.0*np.pi*(radius)**2.0)
                #data = data[np.where(data['MAB_dA_SDSS_g']!=0.0)[:][0]]
                print('effective surface brigthness SDSS g-band within gas disk radius [mag arcsec-2]\n')
            except:
                pass
            try:
                radius = mL.conv_radius_pc_to_arcsec(data['rhalf_stars_1.5ropt']*1000.0)
                data['SB_mu_eff_stars_1.5ropt_g']=data['MAB_dA_SDSS_g']+2.5*np.log10(2.0*np.pi*(radius)**2.0)
                #data = data[np.where(data['MAB_dA_SDSS_g']!=0.0)[:][0]]
                print('stellar surface brigthness SDSS g-band within half mass radius stars [mag arcsec-2]\n')
            except:
                pass        
            try:
                radius = mL.conv_radius_pc_to_arcsec(data['rhalf_stars_30kpc']*1000.0)
                data['SB_mu_eff_stars_30kpc_g']=data['MAB_dA_SDSS_g']+2.5*np.log10(2.0*np.pi*(radius)**2.0)
                #data = data[np.where(data['MAB_dA_SDSS_g']!=0.0)[:][0]]
                print('stellar surface brigthness SDSS g-band within 30kpc aperture [mag arcsec-2]\n')
            except:
                pass        
    
            try:
                data['bheff']=data['mbh']/data['mhalo']
                print('bheff=mbh/mhalo --> black hole efficiency')
            except:
                try:
                    data['bheff']=data['mbh']/data['mhalo_200c']
                    print('bheff=mbh/mhalo_200c --> black hole efficiency')
                except:
                    pass         
    
            try:
                data['mAB_dA_total_cut_r_i'] = data['MAB_dA_SDSS_r']-data['MAB_dA_SDSS_i']
                data['mAB_dA_total_cut_g_r'] = data['MAB_dA_SDSS_g']-data['MAB_dA_SDSS_r'] 
                data['mAB_dA_total_cut_g_i'] = data['MAB_dA_SDSS_g']-data['MAB_dA_SDSS_i']
            except:
                pass
            try:
                data['mAB_dA_total_cut_r_i'] = data['MAB_dA_total_r']-data['MAB_dA_total_i']
                data['mAB_dA_total_cut_g_r'] = data['MAB_dA_total_g']-data['MAB_dA_total_r'] 
                data['mAB_dA_total_cut_g_i'] = data['MAB_dA_total_g']-data['MAB_dA_total_i']          
            except:
                pass        
            try:
                data['mAB_dA_total_cut_r_i'] = data['MAB_dA_SDSS_r']-data['MAB_dA_SDSS_i']
                data['mAB_dA_total_cut_g_r'] = data['MAB_dA_SDSS_g']-data['MAB_dA_SDSS_r'] 
                data['mAB_dA_total_cut_g_i'] = data['MAB_dA_SDSS_g']-data['MAB_dA_SDSS_i']          
            except:
                pass       
            try:
                data['mAB_total_cut_r_i'] = data['MAB_total_r']-data['MAB_total_i']
                data['mAB_total_cut_g_r'] = data['MAB_total_g']-data['MAB_total_r'] 
                data['mAB_total_cut_g_i'] = data['MAB_total_g']-data['MAB_total_i']          
            except:
                pass
            try:
                data['mAB_dA_total_cut_B_V'] = data['MAB_dA_Johnson_B']-data['MAB_dA_Johnson_V']
                data['mAB_dA_total_cut_U_V'] = data['MAB_dA_Johnson_U']-data['MAB_dA_Johnson_V']        
            except:
                pass
            
            try:
                data['Mgas_HIH2_30kpc_GK11'] = data['Mgas_HI_30kpc_GK11']+data['Mgas_H2_30kpc_GK11']        
            except:
                pass       
           
            try:
                from cosmolopy import cparam, cd,cc
                fidcosmo = cparam.PlanckMD(flat=True, extras=False)    
                data['mean_age_stars_disk'] = cd.age(redshift, **fidcosmo)/cc.Gyr_s-(data['age_sfr_int_disk']/data['sfr_int_disk']/1e9)
                
                print('z=', redshift, 'age:', format(cd.age(redshift, **fidcosmo)/cc.Gyr_s, '0.2f'), 'mean_age_stars_disk=age_Universe(z)-age_sfr_int_disk/sfr_int_disk [Gyr]', end='')
                print('min/max:', format(min(data['mean_age_stars_disk']), '0.2f'), '/', format(max(data['mean_age_stars_disk']), '0.2f'))           
            except:
                pass
    
            try:
                from cosmolopy import cparam, cd,cc
                fidcosmo = cparam.PlanckMD(flat=True, extras=False)    
                data['mean_age_stars_spheroid'] = cd.age(redshift, **fidcosmo)/cc.Gyr_s-(data['age_sfr_int_spheroid']/data['sfr_int_spheroid']/1e9)
                
                print('z=', redshift, 'age:', format(cd.age(redshift, **fidcosmo)/cc.Gyr_s, '0.2f'), 'mean_age_stars_spheroid=age_Universe(z)-age_sfr_int_spheroid/sfr_int_spheroid [Gyr]', end='')
                print('min/max:', format(min(data['mean_age_stars_spheroid']), '0.2f'), '/', format(max(data['mean_age_stars_spheroid']), '0.2f'))           
            except:
                pass
            
            try:
                 data['mstar_IC']=data['mstar_IC']+data['mstar_disk']+data['mstar_spheroid']
                 print('Mstar_IC+Mstar')            
            except:
                 pass
             
            try:
                data['jbar'] =  (data['angM_disk']+ data['angM_spheroid'])/(data['mstar'] + data['mcold'])
                print('jbar --> DONE!\n')            
            except:
                 pass
            try:
                data['angM_norm_1.5ropt_bar'] =  data['angM_norm_1.5ropt_gas']+ data['angM_norm_1.5ropt_stars']
                print('jbar_1.5ropt --> DONE!\n')            
            except:
                  pass             
             
            try:
                 data['jdisk'] =  data['angM_disk']/data['mstar_disk']
                 print('jdisk --> DONE!\n')            
            except:
                  pass
    
            try:
                data['jbulge'] =  data['angM_spheroid']/data['mstar_spheroid']
                print('jbulge --> DONE!\n')            
            except:
                 pass
             
            try:
                data['jhotHalo'] =  data['angM_hotHalo']/data['mhot']
                print('jhotHalo --> DONE!\n')            
            except:
                 pass
             
            try:
                data['joutHotHalo'] =  data['angM_outflowHotHalo']/data['mhot_outflow']
                print('joutHotHalo --> DONE!\n')            
            except:
                 pass
             
            try:
                data['delta_age_stars_rband'] =  data['age_stars_rband_r502D'] - data['age_stars_rband_2r502D']
                print('delta_age_stars_rband --> DONE!\n')            
            except:
                 pass
        
        
            try:
                if np.all(data['v200c']==-99.0):
                    GravConst_in_Mpc = 4.3009e-9
                    data['v200c'] =  (GravConst_in_Mpc*data['mhalo_200c']/data['r200c'])
                    print('v200c --> DONE!\n')            
            except:
                 pass
            try:
                data['Vvir'] =  (GravConst_in_Mpc*data['mhalo']/data['rvir'])
                print('Vvir --> DONE!\n')            
            except:
                 pass       
    
          
         
     
            try:
                if np.all(data['NFW_con']==0.0) or np.all(data['NFW_con']==-99.0):
                    data['NFW_con'] = mL.calculate_NFW_con_with_fit(data['mhalo_200c'], 
                                                           float(redshift),
                                                           cosmology='Planck', 
                                                           overdens='200c')
                    print('Calculate NFW_con for 200c and Planck cosmology --> DONE!\n')
            except:
                pass
            
            try:
                NFW_con = mL.calculate_NFW_con_with_fit(data['mhalo_cents'], 
                                                       float(redshift),
                                                       cosmology='Planck', 
                                                       overdens='vir')
                
                data['mhalo_cents_200c']=mL.convert_halo_mass(data['mhalo_cents'], 
                                                                    NFW_con,
                                                                    data['orphan'])
          
                print('--> DONE!\n')
            except:
                pass 
            
            try:
                data['reff_gas_disk_1.5ropt'][np.where(data['np_gas_1.5ropt']<600)[:][0]] = -99.0
                data['reff_gas_1.5ropt'][np.where(data['np_gas_1.5ropt']<600)[:][0]] = -99.0
                data['mgas_1.5ropt'][np.where(data['np_gas_1.5ropt']<600)[:][0]] = -99.0
            except:
                pass
    
    #
    #        try:
    #            if np.all(data['weight_tot']==-99.0):
    #                print('\n+++++++++++++++++++\napply weights to CMASS DR12',
    #                data['weight_tot']=data['weight_systot']*(data['weight_noz']+data['weight_cp']-1)
    #                totgal = sum(data['weight_noz']+data['weight_cp']-1)
    #                weight_max = max(np.cumsum(data['weight_tot']))
    #                
    #                print('ngal weight normal:', sum(data['weight_tot']), 'ngal totgal:', totgal, 'ngal weight_max:', weight_max, '---> ngal special weight:',           
    #                print(sum(data['weight_tot']*totgal/weight_max)
    #                print('--> DONE!\n'                
    #
    #        except:
            try:
                if np.all(data['weight_tot']==-99.0):
                    print('set weights to 1!', end='')
                    data['weight_tot']=np.ones((data.size,), dtype=np.int8)
                    print('--> DONE!\n')
            except:
                pass
        else:
            print('Skipping caluclating properties!!!')
                
        preprocessing_only=False
        
        # data = mL.assign_sats(data,
        #                       name_ID_central='fofID')

        #data = mL.count_nr_sats(data)

        if self.myconfig_array['catname'+str(self.a)].startswith('EAGLE'):
            print('Line 1692: after prop calc', data.size)
            #data = mL.handle_HF_EAGLE_data(data, key='filterData')
            data = mL.handle_HF_EAGLE_data(data, key='calcAssemblyHistory')
            #data = mL.handle_HF_EAGLE_data(data, key='assignSubsamples')
            #data = mL.handle_HF_EAGLE_data(data, key='calcProps')
            #data = mL.handle_HF_EAGLE_data(data, key='generateBoxPlot')

            #data = mL.handle_HF_EAGLE_data(data, key='print2Table')            
            #data = mL.handle_HF_EAGLE_data(data, key='print2TableBins')
            #units2log()
 
        #print(data[['sfe_gas_reff_1.5ropt', 'sfe_gas_reff_disk_1.5ropt', 'sfe_gas_1.5ropt', 'sfe_HIH2_30kpc']]

        # if self.myconfig_array['catname'+str(self.a)].startswith('Galacticus'):
        #     print('Galacticus!')
            #SFH
            #data = mL.handle_SFH_Galacticus_data(data, key='assignSubsamples')
            #data = mL.handle_SFH_Galacticus_data(data, key='print2Table') 
            #mL.handle_SFH_Galacticus_data(data, key='calcStats')
            

        #exit()    
            
        # data['fbar_30kpc']=(data['mstar_30kpc']+data['mgas_30kpc']+data['mbh'])/data['mhalo_200c']
        
        # data['angM_norm_SFgas2stars'] = data['angM_norm_SFgas']/data['angM_norm_stars']
        # data['angM_norm_NSFgas2stars'] = data['angM_norm_NSFgas']/data['angM_norm_stars']
        # data['angM_norm_gas2stars'] = (data['angM_norm_SFgas']+data['angM_norm_NSFgas'])/data['angM_norm_stars']
        
        # data['delta_age_stars_rband']=data['age_stars_rband_2r502D']-data['age_stars_rband_r502D']
        # data['delta_lambda_edgeOn']=data['lambda_2r502D_edgeOn']-data['lambda_r502D_edgeOn']
        # data['delta_vdisp_edgeOn']=data['vdisp_r502D_edgeOn']-data['vdisp_2r502D_edgeOn']
        # data['delta_VvS_edgeOn']=data['VvS_2r502D_edgeOn']-data['VvS_r502D_edgeOn']
        
        # for ap in ['', '2']:
        #     print('ap:', ap)
        #     for prop, delta in zip(['age_stars_rband_'+ap+'r502D','lambda_'+ap+'r502D_edgeOn','vdisp_'+ap+'r502D_edgeOn', 'VvS_'+ap+'r502D_edgeOn'],
        #                            ['delta_age_stars_rband','delta_lambda_edgeOn','delta_vdisp_edgeOn', 'delta_VvS_edgeOn']):
        #         print('delta:', delta, 'prop:', prop, '-->', data[np.where(data[prop]==-99.0)[:][0]].size, 'found!')
        #         data[delta][np.where(data[prop]==-99.0)[:][0]]=-99.0
                
                #print(data[np.where(data[prop]==-99.0)[:][0]][[prop,delta]]
                
                
        #units2log()
        #exit()
                

###########################################################################################################################################################

        #data=data[np.where(data['mhalo_200c']>4.46e10)[:][0]]  

        # print(data['OH_gas_disk_bulge']
        # print(min(data['OH_gas_disk_bulge']), max(data['OH_gas_disk_bulge'])
        
        # exit()

#        import h5py as hdf5
#        path = '/data/256_hydro_50Mpc_halo_tests/65_particles.h5.3'
#        f=hdf5.File(path, "r")
#        
#        print(f['particle_IDs'][0:100]
#
#        ha_lib.scan_file_format_Cholla(f,path)  
        #Constructing main progenitor merger trees from ROCKSTAR halo catalog
        #1) assign predIndex
#        try:
#            redshift_before=float(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i-1)]+'_snapid'+str(self.i-1)]['z'])
#        except:
#            redshift_before=redshift
        #data = ha.assign_progenitor_indices(data,redshift_before)
        
        #2) construct the merger trees     
      
        #data=ha.construct_merger_tree(data)

        
        #3) test the output by comparing with the particle informaiton in the hydro simulation files
        #myHaloAnalysis.MAIN()
 
        #work with Gadget and Wfirst data from Nicole
        #ha.load_binary_GADGET()
        
        #work with Amiga Halo Finder and data from Roman (former WFIRST) from Nicole
        #ha.load_AHF()

        #mL.generate_snapidzred_file(258,end_snapid=0)
        #exit()

      
#        mL.property_stats(mydata)
#        exit()        

##########################################################################################################################################################        

        if preprocessing_only==False:        
            if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion']!='False':
                selectRegion(data, myheader, mylog, mylog_region)                                                                                                       
    
            if redshift=='False':
                try:
                    redshift=float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_manual_input_redshift'])
                    print('redshift reset to:', redshift)
                except:
                    print('redshift could not be reset! Current redshift value is:', redshift)
                    
    
            if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_convert_sph_to_cart_coords']!='False':
                try:
                    data=mL.convert_coordinates(data)
                    print('data new coords:', 'x_pos','y_pos','z_pos','RA','DEC','Z')
                    print(data[0:10][['x_pos','y_pos','z_pos','RA','DEC','Z']])
                except:
                    print('Converstion from spherical to euclidian coordinates is not possible!\n')


            #Reset MD-output name & units conventions to DS
            if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_UNIT_CODE']=='MD2DS':
                for item in self.name_conv_map: 
                    #print(myconds_array[item+'_output_col_name'], item
                    myconds_array[item+'_output_col_name']=item

            data, myheader, mydataformat, mylog, check_size, sel_col_list = mL.filter_data_before_write2file(data,
                                                                                                            myconds_array,
                                                                                                            myconds_array,
                                                                                                            myheader,
                                                                                                            mylog,
                                                                                                            set_header_Topcat_format=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_set_header_Topcat_format'] )
    
            try:
                name = self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filter_density0']
            except:
                name = self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filter_density']  
                       
            if name!='False':
                filterDensity(myheader, mylog)            
                  
            
            if subcatsel_code!=False and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_output_filename_code']=='':
                self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_output_filename_code'] = '_subcat'
    
            if check_size==True:          
                myfilename=myhdf5_filename[0:len(myhdf5_filename)-5]+'.txt'
    
                #print('myheader:', myheader
                #print('mydataformat:', str(mydataformat)
            else:
                print('Sorry, but we ran out of galaxies!!!\nThe selections of your data sample has zero entries! Try again with a different selection!!!\nProgramm exiting ...')
                exit()
    
    
            #select randomly a subsample of data if the catalouge has more galaxies as your threshold!    
            if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_sample_info'].find('random')!=-1:            
                mylog+='\nCatalogue is randomly selected: Yes! ngalaxies selected? '+str(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_ngal_random_sample'])           
                data= mL.choose_random_sample(data, self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_ngal_random_sample'])
               
                    
            mylog+='\nngalaxies in the file: '+str(len(data[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name0']]))+'\nlink and filename of catalogue: '+myhdf5_filename
            print(myhdf5_filename)
       
            #writeCoordinates(data, myhdf5_filename)
            if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_OUTPUT_ASCII']!='False':				
                writeASCII(data[sel_col_list], myhdf5_filename[0:len(myhdf5_filename)-5]+'.txt', myheader, mydataformat)
            writeHDF5(data, myhdf5_filename)
            writeLOGFILE(mylog, myhdf5_filename[0:len(myhdf5_filename)-5]+'_log.txt')    
        else:
            print('ONLY PREPROZESSING DONE!')
            self.myData2Process=data
            
                  
    def analyseTargetSelection(self,
                                redshift,
                                filename,
                                myconds_array,
                                plot_key):

    
        print('analysisTargetSelection()', 'self.i:', self.i, 'self.a:', self.a, '++++++++++++++++++++++++++++++++++++++++++++++\n')

        def rawSAMHDF5Catalouge():
            print('rawSAMHDF5Catalouge()\n+++++++++++++++++++++')
            
            data = myData.readAnyFormat(config=True, 
                                        mypath=filename,
                                        data_format='HDF5', 
                                        data_shape='shaped', 
                                        nr_rows=500000000, 
                                        nr_col=myconds_array['nr_entries'], 
                                        delim='  ', 
                                        mydtype=np.float,
                                        id_col_array=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array'],
                                        skiprow=2)
              
                
            self.volume = mL.survey_VolumeSqDeg_to_Mpc(box_size=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'], 
                                   redshift=redshift,
                                   h=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_hubble_par'],
                                   scale_factor=self.SAM_scale_factor_map[self.myconfig_array['catname'+str(self.a)]+'_'+str(redshift)])
                                   
            cat_info = {'volume': self.volume, 'redshift': redshift, 'scaleFactor': self.SAM_scale_factor_map[self.myconfig_array['catname'+str(self.a)]+'_'+str(redshift)]}
          
            return data, cat_info
                        
        def correctedCatalouge():
            print('correctedCatalouge()\n+++++++++++++++++++++')
            self.readData(myfilename=filename)
            
            return self.myData2Process, myData.cat_attributes

        
        def caseSwitcher(fileformat):
        
            choose = {
                'SAMHDF5': rawSAMHDF5Catalouge,
                'HDF5': correctedCatalouge
                }
                
            func = choose.get(fileformat)
            return func()
        
        data, cat_info = caseSwitcher(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format'])        

             
       
        #Filter galaxy cat before intersecting with halo cat to exclude galaxy previously in order to save calc time!!!
        myheader=self.myconfig_array['catname'+str(self.a)]+' z='+str(redshift)+'\n'  
        mylog='Captains log-file: Stardate '+str(date_time_stamp)+'\n'+myheader+'++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n'+'box volume: '+str(self.volume)+' [Mpc3]\n\nPreselection before matching galaxy cat with halo cat:\n\n'

        if self.myconfig_array['catname'+str(self.a)].startswith('SAGE') and self.myconfig_array['catname'+str(self.a)].find('run')!=-1: 
            self.correctAndConvertUnits()           
 
        
        if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filterSample']=='True':        

            data[np.argsort(data[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['haloid_col_id']])]            
            
        
            data, myfilename, myheader, mydataformat, mylogfilename, mylog, sel_col_list = self.filterAndExtractData(myconds_array,
                                                                                                       data=data, 
                                                                                                       myheader=myheader,
                                                                                                       myselection_name=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selection_name'],
                                                                                                       redshift=cat_info['redshift'],        
                                                                                                       scale_factor=cat_info['scaleFactor'])
                                                                                                       
            exit()                                                                                              
        
        if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_calcHistoSample']=='True':         

            #Calculate a histo from all coluumns!

            def caseSwitcher(histo_method):
            
                choose = {
                    'fast':  self.calcFastHisto,               
                    }
                    
                func = choose.get(histo_method)
                return func()
            
            caseSwitcher('fast')

                       
        if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_plotSample']=='True':

            print('analyseTargetSelection() --> plot sample ...\n-------------------')
           
            i=0           
            while i<self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_analyse_tarsel_id_col_array']['nr_entries']:
                tarsel_plot_name = str(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_analyse_tarsel_id_col_array']['name'+str(i)])
                              
                name_y=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_analyse_tarsel_id_col_array'][tarsel_plot_name+'_col_name1']
                name_x=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_analyse_tarsel_id_col_array'][tarsel_plot_name+'_col_name2']
                name_z=False
                if name_z==False:
                    name_z=name_x

                #data=data[np.where(data['mstar']>10**10.8)[:][0]]
                #data=data[np.where(data['orphan']>0)[:][0]]
                name_weights='weight_tot'

                #data = data[(np.where(data['MAB_dA_total_g']-data['MAB_dA_total_i']>2.35))[:][0]]
                #x=data['mAB_dA_total_g']-data['mAB_dA_total_r']                        
                #y=data['mAB_dA_total_r']-data['mAB_dA_total_i']
                #prop=(y+1.05)/1.3-x
                #print(min(prop), max(prop)
                #data = data[np.where(prop<0.0)[:][0]]                
#                print('name_x', name_x, 'name_y', name_y, 'name_z:', name_z

#                data = myData.selectData2Compute(data, 
#                                                 selected_col='orphan', 
#                                                 operator='==', 
#                                                 condition=2)                                
                
                try:
                    if self.myconfig_array['catname'+str(self.a)].find('SAdG_')!=-1:
                        data['mstar']=data['mstar_disk']+data['mstar_spheroid']
                        data['mcold']=data['mcold_disk']+data['mcold_spheroid']
                        data['Mzgas']=data['Mzgas_disk']+data['Mzgas_spheroid']

                    mstar_min=5e8
                    
                    #data = mL.choose_random_sample(data,int(data.size*0.10))
                    
                    #mstar_min=10**(11.29-0.039+np.log10((0.73/0.6777)**2))                        
                    print('Mstar > ', mstar_min, 'selection:', data.shape, '-->', end='')
                    data = myData.selectData2Compute(data, 
                                                     selected_col='mstar', 
                                                     operator='>', 
                                                     condition=mstar_min)
                    print(data.shape)
                                 
                except:
                    try:
                        print('Mbulge >', mstar_min, 'selection:', data.shape, '-->', end='')
                        data = myData.selectData2Compute(data, 
                                                         selected_col='mstar_spheroid', 
                                                         operator='>', 
                                                         condition=mstar_min)
                        print(data.shape)
                    except:
                        pass


                if tarsel_plot_name.find('oh_mstar')!=-1:
                    print('oh vs. mstar! name_y:', name_y)
                    data[name_y]=12+np.log10(0.0356*data['Mzgas']/data['mcold'])
                        
                    myconds_array[name_y+'_name_in_plot'] = r'12+ $\log_{10}$(OH)'  

                if tarsel_plot_name.find('_sfr_mstar')!=-1 or tarsel_plot_name.find('mcold_mstar')!=-1:
                    print('sfr vs. mstar or mcold vs. mstar! name_x:', name_x, 'name_y:', name_y)
                    myconds_array[name_y+'_name_in_plot']=r'$\log_{10}$ ($M_{Cold}/M_{_*}$)'

                    data = myData.selectData2Compute(data, 
                                                    selected_col=name_y, 
                                                    operator='>', 
                                                    condition=0.0)
                    
                    if tarsel_plot_name.find('mcold_mstar')!=-1:
                        data[name_y]/=data[name_x]

                    data[name_x]=np.log10(data[name_x])
                    data[name_y]=np.log10(data[name_y])

                    data = myData.selectData2Compute(data, 
                                                    selected_col=name_y, 
                                                    operator='>', 
                                                    condition=-16)
                    print('ssfr > 1e-16 -->', data.shape)                      
                    
                                                    
                if tarsel_plot_name.find('mbh_mstar')!=-1:
                    print('mbh vs mbulge!')
                    try:
                        data['mstar_spheroid']+=data['mcold_spheroid']                    
                    except:
                        print('no mcold bulge found!')
                        
                    data = myData.selectData2Compute(data, 
                                                     selected_col=name_y, 
                                                     operator='>', 
                                                     condition=0.0)                       

                if tarsel_plot_name.find('mhalo_mstar')!=-1:
                    print('mstar vs. mhalo!')
                    data = myData.selectData2Compute(data, 
                                                    selected_col='orphan', 
                                                    operator='<', 
                                                    condition=1)
                    print('data after central selection!', data.shape)
                
                if tarsel_plot_name.find('zgas')!=-1:                           

#                    data = myData.selectData2Compute(data, 
#                                                     selected_col=name_x, 
#                                                     operator='>', 
#                                                     condition=1e11)
                    
#                    if self.myconfig_array['catname'+str(self.a)]=='Galacticus_1Gpc':
#                        print('name_y:', name_y
#                        print('Galacticus: zgas vs. mstar/cold! --> Zgas = Zgas_spheroid+Zgas_disk/Mcold (column Zgas_spheroid and Zgas_disk are wrongly named in the catalog it should be Mzgas_spheroid and Mzgas_disk)'
#                        data[name_y]+=data['zgas_disk']

                    data = myData.selectData2Compute(data, 
                                                    selected_col=name_y, 
                                                    operator='>', 
                                                    condition=0.0)
                    print('Mzgas > 0.0 -->', data.shape)

                    try:
                        data = myData.selectData2Compute(data, 
                                                    selected_col='mcold', 
                                                    operator='>', 
                                                    condition=0.0)
                        data[name_y]=data[name_y]/(data['mcold']*0.0134)
                        
                    except:
                        data = myData.selectData2Compute(data, 
                                                    selected_col='mcold_disk', 
                                                    operator='>', 
                                                    condition=0.0)
                        data[name_y]=data[name_y]/(data['mcold_disk']*0.0134)                    
                    

                                                    
                    data[name_y]=8.69+np.log10(data[name_y])
                    myconds_array[name_y+'_name_in_plot'] = r'$Z_{Cold}$'

                    if tarsel_plot_name.find('zgas_mcold')!=-1:                  
                        data[name_x]/=data['mstar']
                        myconds_array[name_x+'_name_in_plot'] = r'$\log_{10}$ ($M_{Cold}/M_{*}$)'   
                        
                        data = myData.selectData2Compute(data, 
                                                         selected_col=name_x, 
                                                         operator='>', 
                                                         condition=0.001)                        

                        data = myData.selectData2Compute(data, 
                                                         selected_col=name_x, 
                                                         operator='<', 
                                                         condition=100)

                if tarsel_plot_name.find('-')!=-1: 
                    print('color vs color! name_x:', name_x, 'name_y:', name_y)
                    print(tarsel_plot_name)
                    
#                    if name_x.find('MAB')==-1:
#                        data = myData.selectData2Compute(data, 
#                                                         selected_col=name_x, 
#                                                         operator='>', 
#                                                         condition=0.0)
#                    if name_y.find('MAB')==-1:
#                        data = myData.selectData2Compute(data, 
#                                                         selected_col=name_y, 
#                                                         operator='>', 
#                                                         condition=0.0)
#                        
#                    print('> 0.0 -->', data.shape
                                           
                    try:
                        print('version:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['mAB_dA_total_r_col_id'], 'mAB_dA')
                        version = 'mAB_dA'
                    except:
                        try:
                            print('version:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['MAB_dA_total_r_col_id'], 'MAB_dA')
                            version = 'MAB_dA'
                        except:
                            print('version:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['mAB_total_r_col_id'], 'mAB')
                            version = 'mAB'                    
                    #print(tarsel_plot_name, version   
                    if tarsel_plot_name.find('r-i')!=-1 and tarsel_plot_name.find('g-r')==-1:
                        print('HERE: r-i')
                        name_y=version+'_total_r'
                        name_x=version+'_total_i'
                        myconds_array[name_y+'_name_in_plot'] = '$r-i$'
                        data[name_y]-=data[name_x]
                        if tarsel_plot_name.find('mstar')!=-1:                        
                            name_x='mstar'                                              
                            data['mstar']= np.log10(data['mstar'])
                        elif tarsel_plot_name.find('i_i')!=-1:
                            data[name_x]-=5*np.log10(0.6777)
                    if tarsel_plot_name.find('g-r')!=-1 and tarsel_plot_name.find('u-g')!=-1:
                        print('HERE: g-r vs u-g')
                        name_y=version+'_total_g'
                        name_x=version+'_total_u' 
                        data[name_x]=data[version+'_total_u']-data[version+'_total_g']
                        data[name_y]=data[version+'_total_g']-data[version+'_total_r']
            
                    if tarsel_plot_name.find('g-r')!=-1 and tarsel_plot_name.find('r-i')!=-1:
                        print('HERE: g-r vs r-i')
                        name_x=version+'_total_g' 
                        name_y=version+'_total_r'  
                        data[name_x]-=data[version+'_total_r']
                        data[name_y]-=data[version+'_total_i']
                     

                    if tarsel_plot_name.find('g-r')!=-1 and tarsel_plot_name.find('g-i')!=-1:
                        print('HERE 2043: g-r vs g-i')
                        name_x=version+'_total_i' 
                        name_y=version+'_total_g'
                        data[name_x]=data[version+'_total_g']-data[version+'_total_i']
                        data[name_y]-=data[version+'_total_r']

                    if tarsel_plot_name.find('g-i')!=-1 and tarsel_plot_name.find('g-r')==-1:
                        print('HERE: g-i')
                        name_y=version+'_total_g'
                        name_x=version+'_total_i'
                        data[name_y]-=data[version+'_total_i']
                        myconds_array[name_y+'_name_in_plot'] = '$g-i$'
                        
                        if tarsel_plot_name.find('mstar')!=-1:
                            data[name_x]=np.log10(data['mstar'])
                        
                    if tarsel_plot_name.find('u-r')!=-1:
                        print('HERE: u-r', end='')
                        name_y=version+'_total_u'
                        name_x=version+'_total_r' 
                        myconds_array[name_y+'_name_in_plot'] = '$u-r$'
                        print('name_x:', name_x, 'name_y:', name_y)
                        data[name_y]-=data[name_x]
                        data[name_x]-=5*np.log10(0.6777)
                        print('check!')

                    if tarsel_plot_name.find('g-u')!=-1:
                        print('HERE: g-u')
                        name_y=version+'_total_g'
                        name_x=version+'_total_u' 
                        myconds_array[name_y+'_name_in_plot'] = '$g-u$'                        

                if tarsel_plot_name.find('ssfr')!=-1 or tarsel_plot_name.find('zcold')!=-1 or tarsel_plot_name.find('mhalo')!=-1:
                    print('name_x:', end='')
                    if tarsel_plot_name.find('ssfr')!=-1 and tarsel_plot_name.find('mstar')==-1:
                        try:
                            name_x='ssfr'
                            data[name_x]=np.log10(data['ssfr'])
                        except:
                            name_x='sfr'
                            data[name_x]=np.log10(data['sfr']/data['mstar'])                            
                        print( name_x, 'name_y:', name_y)
    
                        data = myData.selectData2Compute(data, 
                                                        selected_col=name_x, 
                                                        operator='>', 
                                                        condition=-16)
                        print('ssfr > 1e-16 -->', data.shape)                        
                  
                    if tarsel_plot_name.find('g-r')!=-1:
                        print('HERE: g-r vs ssfr! name_x', name_x, 'name_y', name_y)
                        data[name_y]=data[version+'_total_g']-data[version+'_total_r']
                        
                    if tarsel_plot_name.find('dmesa')!=-1:
                        print('HERE: dmesa vs ssfr! name_x:', name_x, 'name_y:', name_y)
                    
                    if tarsel_plot_name.find('mhalo')!=-1:
                        name_x='mhalo'
                        print('HERE: ssfr vs mhalo! or g-i vs mhalo, name_x:', name_x, 'name_y:', name_y)
                        print('only centrals!')
                        data = myData.selectData2Compute(data, 
                                                        selected_col='orphan', 
                                                        operator='==', 
                                                        condition=0)                           
                    
#                        print('check environment!'
#                        data = myData.selectData2Compute(data, 
#                                                        selected_col='env_1024', 
#                                                        operator='==', 
#                                                        condition=2)

                        if tarsel_plot_name.find('mstar')!=-1:
                            data['mstar']= np.log10(data['mstar'])
                              
                        try:
                            data[name_x]= np.log10(data['mhalo_200c'])
                        except:
                            data[name_x]= np.log10(data['mhalo'])

                    if tarsel_plot_name.find('mstar')!=-1 and tarsel_plot_name.find('mhalo')==-1:
                        name_x='mstar'                                              
                        data['mstar']= np.log10(data['mstar'])

                    if tarsel_plot_name.find('ssfr_mstar')!=-1:                       
                        name_y='ssfr'                                              
                        try:
                            data[name_y]=np.log10(data['ssfr'])
                        except:
                            name_y='sfr'
                            data[name_y]=np.log10(data['sfr']/data['mstar'])                        

                        name_x='mstar'                                              
                        #data['mstar']=np.log10(data['mstar']) 
                        
                    if tarsel_plot_name.find('_sfr')!=-1:
                        name_x='sfr'
                        data[name_x]= np.log10(data['sfr'])                                             

                    if tarsel_plot_name.find('ssfr_zcold')!=-1:
                        name_x='zcold'
                        name_y='ssfr'   

                    if tarsel_plot_name.find('ssfr_i')!=-1:
                        name_x='mAB_dA_total_i'
                        name_y='ssfr'

                    if tarsel_plot_name.find('ssfr_r')!=-1:
                        name_x='mAB_dA_total_r'
                        name_y='ssfr'                                           

                    if tarsel_plot_name.find('ssfr_mcold')!=-1:
                        name_x='mcold'
                        name_y='ssfr'                  

                        data[name_x]/=data['mstar']
                        data[name_x]=np.log10(data[name_x])
                        myconds_array[name_x+'_name_in_plot'] = r'$\log_{10}$ ($M_{Cold}/M_{*}$)'  

                if tarsel_plot_name.find('dmesa')!=-1:
                    print('HERE: dmesa vs.', end='') 
                    if tarsel_plot_name.find('dmesa_i')!=-1:
                        print('i!')
                    elif tarsel_plot_name.find('dmesa_mstar')!=-1:
                        print('mstar!')                    
                    myconds_array[name_y+'_name_in_plot'] = r'$d_{\perp}$'

                if tarsel_plot_name.find('_mstar')!=-1 and tarsel_plot_name.find('_mstar_sph')==-1:
                    print('set x-axis to mstar! x_name:', name_x)
                    myconds_array[name_x+'_name_in_plot'] = r'$\log_{10}$ ($M_{_*}$ [$M_{\odot}$])'
                    #data[name_x] = np.log10(data['mstar'])

                if tarsel_plot_name.find('_rdisk')!=-1:
                    print('rdisk in kpc and log10!')
                    data[name_y]= np.log10(data['rdisk']*1e3)

                if tarsel_plot_name.find('_rbulge')!=-1:
                    print('rbulge in kpc and log10!')
                    
                    data[name_x]= np.log10(data['rbulge']*1e3)                   
                    
                if tarsel_plot_name.find('_rdisk')!=-1 and tarsel_plot_name.find('_rbulge')!=-1:
                    print('rdisk vs rbulge in kpc! x_name:', name_x, 'name_y:', name_y)
                    
                    data[name_x]= np.log10(data['rbulge']*1e3)
                    data[name_y]= np.log10(data['rdisk']*1e3)   

                elif tarsel_plot_name.find('rdisk')!=-1 and tarsel_plot_name.find('mstar')!=-1:
                    print('rdisk in kpc vs mstar! x_name:', name_x, 'name_y:', name_y)
                    data[name_x]= np.log10(data['mstar'])
                    data[name_y]= np.log10(data['rdisk']*1e3)

                elif tarsel_plot_name.find('rhalfdisk')!=-1 and tarsel_plot_name.find('mstar')!=-1:
                    print('rhalf_disk in kpc vs mstar! x_name:', name_x, 'name_y:', name_y)
                    data[name_x]= np.log10(data['mstar'])
                    data[name_y]= np.log10(data['rhalf_disk']*1e3)

                elif tarsel_plot_name.find('rbulge_')!=-1 and tarsel_plot_name.find('_rhalf')!=-1:
                    print('rhalf vs rbulge in kpc! x_name:', name_x, 'name_y:', name_y)                                  
                    data[name_x]= np.log10(data['rbulge']*1e3)
                    data[name_y]= np.log10(data['rhalf_mass']*1e3)
                    
                elif tarsel_plot_name.find('rhalfbulge')!=-1 and tarsel_plot_name.find('mstar')!=-1:
                    print('rhalf_bulge in kpc vs mstar! x_name:', name_x, 'name_y:', name_y)
                    data[name_x]= np.log10(data['mstar'])
                    data[name_y]= np.log10(data['rhalf_bulge']*1e3)                    

                elif tarsel_plot_name.find('rhalfmass')!=-1 and tarsel_plot_name.find('mstar')!=-1:
                    print('rhalf_mass in kpc vs mstar! x_name:', name_x, 'name_y:', name_y)
                    data[name_x]= np.log10(data['mstar'])
                    data[name_y]= np.log10(data['rhalf_mass']*1e3)

                elif tarsel_plot_name.find('bdisk_bbulge')!=-1:
                    name_x='angM_spheroid'
                    name_y='angM_disk'
                    print('b_disk vs b_bulge! x_name:', name_x, 'name_y:', name_y)
                    data[name_y]=np.log10(data['angM_disk']/(data['mstar_disk']+data['mcold_disk']))-2/3.0*np.log10(data['mstar_disk']+data['mcold_disk'])
                    data[name_x]=np.log10(data['angM_spheroid']/(data['mstar_spheroid']+data['mcold_spheroid']))-2/3.0*np.log10(data['mstar_spheroid']+data['mcold_spheroid'])        

                elif tarsel_plot_name.find('angM_mstar')!=-1:
                    name_x='mstar'
                    name_y='angM_disk'
                    print('j vs mstar! x_name:', name_x, 'name_y:', name_y)
                    data[name_y]= np.log10((data['angM_disk']+data['angM_spheroid'])/data['mstar'])
                    data[name_x]= np.log10(data['mstar'])

                elif tarsel_plot_name.find('angM_mdisk')!=-1:
                    name_x='mstar'
                    name_y='angM_disk'
                    print('j_disk vs mdisk! x_name:', name_x, 'name_y:', name_y)
                    data[name_y]= np.log10(data['angM_disk']/data['mstar_disk'])
                    data[name_x]= np.log10(data['mstar_disk'])                    

                elif tarsel_plot_name.find('angM_mbulge')!=-1:
                    name_x='mstar'
                    name_y='angM_spheroid'
                    print('j_bulge vs mbulge! x_name:', name_x, 'name_y:', name_y)
                    data[name_y]= np.log10(data['angM_spheroid']/data['mstar_spheroid'])
                    data[name_x]= np.log10(data['mstar_spheroid'])                     

                elif tarsel_plot_name.find('angM_mbar')!=-1:
                    name_x='mstar'
                    name_y='angM_disk'
                    print('jbar vs mbar! x_name:', name_x, 'name_y:', name_y)
                    data[name_y]= np.log10((data['angM_disk']+data['angM_spheroid'])/(data['mstar']+data['mcold']))
                    data[name_x]= np.log10(data['mstar']+data['mcold'])


                elif tarsel_plot_name.find('angMdisk')!=-1 and tarsel_plot_name.find('mbar')!=-1:
                    name_x='mstar'
                    name_y='angM_disk'
                    print('jdisk vs mbar disk! x_name:', name_x, 'name_y:', name_y)
                    data[name_y]= np.log10(data['angM_disk']/(data['mstar_disk']+data['mcold_disk']))
                    data[name_x]= np.log10(data['mstar_disk']+data['mcold_disk'])

                elif str(tarsel_plot_name).find('angMspheroid')!=-1 and tarsel_plot_name.find('mbar')!=-1:
                    name_x='mstar'
                    name_y='angM_spheroid'
                    print('jbluge vs mbar bulge! x_name:', name_x, 'name_y:', name_y)
                    data[name_y]= np.log10(data['angM_spheroid']/(data['mstar_spheroid']+data['mcold_spheroid']))
                    data[name_x]= np.log10(data['mstar_spheroid']+data['mcold_spheroid'])                
                    
                for name in [name_x,name_y,name_z]:
                    print('set axis labels! name:', name)
                    if name=='mstar':
                        if tarsel_plot_name.find('mbar')!=-1:
                            myconds_array[name+'_name_in_plot']=r'$\log_{10}$ ($M_{bar}$ [$M_{\odot}$])'
                        else:
                            myconds_array[name+'_name_in_plot']=r'$\log_{10}$ ($M_{_*}$ [$M_{\odot}$])'
                            
                        if tarsel_plot_name.find('angMdisk_mbar')!=-1:
                            myconds_array[name+'_name_in_plot']=r'$\log_{10}$ ($M^{disk}_{bar}$ [$M_{\odot}$])'
                        elif tarsel_plot_name.find('angMspheroid_mbar')!=-1:
                            myconds_array[name+'_name_in_plot']=r'$\log_{10}$ ($M^{bulge}_{bar}$ [$M_{\odot}$])'                            

                    elif name=='Mzstar':
                        myconds_array[name+'_name_in_plot']=r'$\log_{10}$ ($mass metals M_{_*}$ [$M_{\odot}$])'                        
                    elif name=='mhalo':
                        myconds_array[name+'_name_in_plot']=r'$\log_{10}$ ($M_{vir}$ [$M_{\odot}$])'
                    elif name=='mhalo_200c':
                        myconds_array[name+'_name_in_plot']=r'$\log_{10}$ ($M_{200c}$ [$M_{\odot}$])'                        
                    elif name=='mbh':
                        myconds_array[name+'_name_in_plot']=r'$\log_{10}$ ($M_{BH}$ [$M_{\odot}$])'
                    elif name=='mstar_spheroid':
                        #myconds_array[name+'_name_in_plot']=r'$\log_{10}$ ($M_{Spheroid}$+$M_{Gas_{Bulge}}$ [$M_{\odot}$])'
                        myconds_array[name+'_name_in_plot']=r'$\log_{10}$ ($M_{Bulge}$ [$M_{\odot}$])'                          
                    elif name=='MAB_dA_total_i':                       
                        myconds_array[name+'_name_in_plot']='$M_{AB_i}$'
                    elif name=='MAB_dA_total_r':                       
                        myconds_array[name+'_name_in_plot']='$M_{AB_r}$'
                    elif name=='mAB_dA_total_i':                       
                        myconds_array[name+'_name_in_plot']='$m_{AB_i}$'                        
                    elif name=='mAB_dA_total_r':                       
                        myconds_array[name+'_name_in_plot']='$m_{AB_r}$'
                    elif name=='sfr' and str(tarsel_plot_name).find('ssfr_mstar')==-1:
                        myconds_array[name+'_name_in_plot'] = r'$\log_{10}$ ($SFR$ [$M_{\odot}$ $yr^{-1}$])'
                    elif name=='ssfr' or str(tarsel_plot_name).find('ssfr_mstar')!=-1:                 
                        myconds_array[name+'_name_in_plot'] = r'$\log_{10}$ ($sSFR$ [$yr^{-1}$])'                            
                    elif name=='rdisk':
                        myconds_array[name+'_name_in_plot'] = r'$\log_{10}$ ($r_{disk}$ [$kpc$])'
                        myconds_array[name+'_name_in_plot'] = r'$\log_{10}$ (R [$kpc$])'                         
                    elif name=='rbulge':
                        myconds_array[name+'_name_in_plot'] = r'$\log_{10}$ ($r_{bulge}$ [$kpc$])'
                    elif name=='rhalf_mass':
                        myconds_array[name+'_name_in_plot'] = r'$\log_{10}$ ($r_{1/2}$ [$kpc$])'                         
                    elif name=='rhalf_disk':
                        myconds_array[name+'_name_in_plot'] = r'$\log_{10}$ ($r^{disk}_{1/2}$ [$kpc$])'  
                    elif name=='rhalf_bulge':
                        myconds_array[name+'_name_in_plot'] = r'$\log_{10}$ ($r^{bulge}_{1/2}$ [$kpc$])'
                    elif name=='spinParameter':
                        myconds_array[name+'_name_in_plot'] = r'$\log_{10}$ (spin parameter [-])'
                    elif name=='zcold':
                        myconds_array[name+'_name_in_plot'] = '$Z_{Cold}$'
                    elif name=='angM_disk':
                        if tarsel_plot_name.find('angM_mbar')!=-1:
                            myconds_array[name+'_name_in_plot'] = '$j_{bar}$ [$kpc$ $km$ $s^{-1}$]'
                        elif tarsel_plot_name.find('angM_mdisk')!=-1:
                            myconds_array[name+'_name_in_plot'] = '$j_{disk}$ [$kpc$ $km$ $s^{-1}$]'
                        elif tarsel_plot_name.find('angM_mstar')!=-1:
                            myconds_array[name+'_name_in_plot'] = '$j_{*}$ [$kpc$ $km$ $s^{-1}$]'                               
                        else:
                            myconds_array[name+'_name_in_plot'] = '$j^{disk}_{bar}$ [$kpc$ $km$ $s^{-1}$]'
                            if tarsel_plot_name.find('bdisk')!=-1:
                                myconds_array[name+'_name_in_plot'] = '$b^{disk}_{bar}$' 
                    elif name=='angM_spheroid':
                        myconds_array[name+'_name_in_plot'] = '$j^{bulge}_{bar}$ [$kpc$ $km$ $s^{-1}$]'
                        if tarsel_plot_name.find('angM_mbulge')!=-1:
                            myconds_array[name+'_name_in_plot'] = '$j_{bulge}$ [$kpc$ $km$ $s^{-1}$]'
                        elif tarsel_plot_name.find('bbulge')!=-1:
                            myconds_array[name+'_name_in_plot'] = '$b^{bulge}_{bar}$'                          
                        
#                    if tarsel_plot_name.find('_mstar')!=-1 and tarsel_plot_name.find('_mstar_sph')==-1:
#                        myconds_array[name_x+'_name_in_plot'] = r'$\log_{10}$ ($M_{_*}$ [$M_{\odot}$])'
                    if tarsel_plot_name.find('r-i')!=-1:
                        myconds_array[name_y+'_name_in_plot'] = '$r-i$'
                    if tarsel_plot_name.find('g-r')!=-1:
                        myconds_array[name_y+'_name_in_plot'] = '$g-r$'                         
                    if tarsel_plot_name.find('g-r')!=-1 and tarsel_plot_name.find('u-g')!=-1:
                        myconds_array[name_y+'_name_in_plot'] = '$g-r$'
                        myconds_array[name_x+'_name_in_plot'] = '$u-g$'
                    if tarsel_plot_name.find('r-i')!=-1 and tarsel_plot_name.find('g-r')!=-1:
                        myconds_array[name_y+'_name_in_plot'] = '$r-i$'
                        myconds_array[name_x+'_name_in_plot'] = '$g-r$'
                    if tarsel_plot_name.find('g-r')!=-1 and tarsel_plot_name.find('g-i')!=-1:
                        myconds_array[name_y+'_name_in_plot'] = '$g-r$'
                        myconds_array[name_x+'_name_in_plot'] = '$g-i$'                         
#                    elif name=='zgas_spheroid':
#                        myconds_array[name_y+'_name_in_plot'] = r'$\log$ ($Z_{Gas}$)'                         
                    #a+=1

#                y=np.log10(data['ssfr'])
#                x=np.log10(data['sfr'])
#                prop=(y+11.16)/1.12-x
#                data = data[np.where(prop>0.0)[:][0]] 

                mask=np.where(np.isfinite(data[name_x]))
                data=data[mask[:][0]]
                mask=np.where(np.isfinite(data[name_y]))
                data=data[mask[:][0]]                
                
                print('min/max x:', min(data[name_x]), max(data[name_x]))          
                print('min/max y:', min(data[name_y]), max(data[name_y]))
                print('name y/x:', name_y, name_x, 'name_in_plot y/x:', myconds_array[name_y+'_name_in_plot'], myconds_array[name_x+'_name_in_plot'], 'data.shape:', data.shape)
                #if len(data[0])>1e6:
                    #data = mL.choose_random_sample(data, 100000)
                if np.all(data[name_weights]==-99.0) or np.all(data[name_weights]==0.0):
                    print('set weights to 1!')
                    data[name_weights]=np.ones((data.size,), dtype=np.int8)
   
                    #print(data[name_weights]
                    

                i+=1
                
        print('CHECK:', self.myconfig_array['catname'+str(self.a)]+' z='+str(redshift),  myconds_array[name_x+'_name_in_plot'], myconds_array[name_y+'_name_in_plot'])
        try:
            return data[[name_x, name_y, name_z, name_weights]], self.myconfig_array['catname'+str(self.a)]+' z='+str(redshift), name_x, name_y, name_z, myconds_array[name_x+'_name_in_plot'], myconds_array[name_y+'_name_in_plot'], myconds_array[name_z+'_name_in_plot']
        except:
            #print('except:'
            return data[[name_x, name_y, name_weights]], self.myconfig_array['catname'+str(self.a)]+' z='+str(redshift), name_x, name_y, name_z, name_weights, myconds_array[name_x+'_name_in_plot'], myconds_array[name_y+'_name_in_plot'], myconds_array[name_z+'_name_in_plot'], myconds_array[name_weights+'_name_in_plot']

    def plotXY(self,
                plot_key,
                filename):
       
        print('\nplotXY()', 'self.i:', self.i, 'self.a:', self.a,'\n++++++++++++++++++++++++++++++++++++++\n')

        def test_300():
            import pandas as pd
                     
            for n_cl in list(range(1,325)):
                print('cluster number:', n_cl, end='')
                    
                filename=mycomp+'anaconda/pro/data/300/Galacticus/MDPL2_Galacticus_z0.00_region'+str(n_cl).zfill(4)+'.txt'
                print(filename)
                mydata_names=['ID','hostid','haloid','orphan','x_pos','y_pos','z_pos','x_vel','y_vel','z_vel','mstar','mstar_disk',\
                              'mcold_spheroid','mcold_disk','mhot','mbh','sfr','sfr_spheroid','sfr_disk','mean_age_stars','mhalo','rest',\
                              'Vmax','Vpeak','NFW_con','spin','MZstar_spheroid','Mzstar_disk', 'Mzgas_spheroid', 'Mzgas_disk','mzhot_halo',\
                              'L_SDSS_dA_total_u','L_SDSS_dA_total_g','L_SDSS_dA_total_r','L_SDSS_dA_total_i','L_SDSS_dA_total_z',\
                              'MAB_dA_total_u','MAB_dA_total_g','MAB_dA_total_r','MAB_dA_total_i','MAB_dA_total_z',\
                              'rdisk','rbulge','rhalf_mass','mhot_outflow']
                
                df=pd.read_csv(filename, skiprows=2, names=mydata_names, sep='  ')
                df['ID']=n_cl
                df['ID']=df['ID'].apply(lambda x: '{0:0>4}'.format(x))
                
                data = mL.df_to_sarray(df)                
                
                
                data=data[np.where(data['mhalo']==max(data['mhalo']))[:][0]]


                print(data['ID'][0])
                data['mstar']+=data['mstar_disk']
                #print('data:', data[['haloid','hostid','orphan','mhalo','mstar']]
                
                
                
                if n_cl==1:
                    mydata=data
                else:
                    mydata = np.append(mydata, data, axis=0)
                    
                

                 
            print('size:', mydata.size)
            #print(mydata
            
            mydata['mhalo']/=0.6777
            mydata['mstar']/=0.6777

            mydata[::-1].sort(order=['mhalo'], axis=0)

            myOutput.writeIntoFile(mycomp+'anaconda/pro/data/300/Galacticus/MDPL2_Galacticus_z0.00_cluster_regions_haloids.txt',
                       mydata[['ID','haloid','hostid','orphan','mhalo','mstar']],
                       myheader= 'z=0.0 \n(1) cluster (2) haloid (3) hostid (4) orphan (5) mhalo [Msun] (6) mstar [Msun]',
                       data_format="%s\t%i\t%i\t%i\t%0.8e\t%0.8e",
                       mydelimiter='\t')                
                

            
        
        def caseSwitcher(myplotXY):
        
            choose = {
                'test_300': test_300                 
                }
                
            func = choose.get(myplotXY)
            return func()
        
        caseSwitcher('test_300')

    def cSFRD(
            self,
            data,
            histo_data_min=0.0,
            histo_data_max=0.0,
            mycond_min=0.0,
            mycond_max=0.0,
            data_offset=0,
            ssfr_cut=None,
            histo_output_offset=0,
            log10_xaxis=False,
            log10_yaxis=False,
            norm_z0=False,
            sfr_cut_min='min',
            sfr_cut_max='max'):
                                       
        print('\n\nSFR vs Redshift(): -->  self.i:', self.i, 'self.a:', self.a,'++++++++++++++++++++++++++++++++++++++++++++++')   
        
        self.histo_data_ZvsSFR = data
        print('mstar min/max:', mycond_min, '/', mycond_max, '\t', end='')
        myfiltered_data = mL.filter_data_before_analysis(self.myData2Process, 
                                                         mycond_min, 
                                                         mycond_max, 
                                                         'mstar')
        
        print('sfr min/max:', sfr_cut_min, '/', sfr_cut_max, '\t', end='')
        myfiltered_data = mL.filter_data_before_analysis(myfiltered_data, 
                                                         sfr_cut_min, 
                                                         sfr_cut_max, 
                                                         'sfr')
        print('ngal after selection -->', myfiltered_data.shape)
        
        print('mstar min/max:', format(min(myfiltered_data['mstar']),'0.4e') ,'/', format(max(myfiltered_data['mstar']),'0.4e'))
        print('sfr min/max:', format(min(myfiltered_data['sfr']),'0.4e') ,'/', format(max(myfiltered_data['sfr']),'0.4e'))        
        
        if ssfr_cut!=None:
            import numpy.lib.recfunctions as rcfuncs
            myfiltered_data = rcfuncs.append_fields([myfiltered_data], ['ssfr'] ,[myfiltered_data['sfr']], usemask=False)
            #print(myfiltered_data[0:5]
            print('expand -->', np.info(myfiltered_data), end='')
            
            myfiltered_data['ssfr']/=myfiltered_data['mstar']
            #print(myfiltered_data[0:5]
    
            #MD-paper ssfr >1e-11
            myfiltered_data = myfiltered_data[np.where(myfiltered_data['ssfr']>1e-11)[0]]
            #MD-paper ssfr >2.17e-11
            #myfiltered_data = myfiltered_data[np.where(myfiltered_data['ssfr']>2.17e-11)[0]]        
            
            print('ssfr max:', max(myfiltered_data['ssfr']), 'min:', min(myfiltered_data['ssfr']))
    
            print('after selection:', myfiltered_data.shape)
        #Attention!
        #This routine is only running if z=0 is calculated first. Although the normalisation of the volume in wrong!
        if self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(0)]+'_snapid'+str(0)]['z']!=0.0:
            print('Self construction sequence initiated ..... you have 5 sec to overthink your choices of redshifts!\n')
            print('... in order to abort self construct choose z=0 as the first redshift to calculate cosmic star formation history!')
            print('This routine is only running if z=0 is calculated first or otherwise the normalisation of the volume will be wrong!')

        self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +4] = self.volume
        print('self.volume:', format(self.volume, '0.2e'), 'sum sfr:', format(np.sum(myfiltered_data['sfr']),'0.2e'), ' / volume z=0:', format(self.histo_data_ZvsSFR[0, data_offset*(self.a) +4], '0.2e'), '=', format(np.sum(myfiltered_data['sfr'])/self.histo_data_ZvsSFR[0, data_offset*(self.a) +4], '0.2e'))  
        
        self.histo_data_ZvsSFR[self.i, data_offset*(self.a)] = self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z']
        self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +1] = np.sum(myfiltered_data['sfr'])/self.histo_data_ZvsSFR[0, data_offset*(self.a) +4]    
        self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +2] = np.sum(myfiltered_data['sfr']) 
    
        #normalise sfr to sfr(z=0)
        self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +3] =  self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +1] / self.histo_data_ZvsSFR[0, data_offset*(self.a) +1]
             
        print('z:', format(self.histo_data_ZvsSFR[self.i, data_offset*(self.a)], '0.3f'))
        print('sum sfr / volume(z=0) [Msun yr-1 Mpc-3]', format(self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +1], '0.2e'))
        print('sum sfr [Msun yr-1]', format(self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +2], '0.2e'))
        print('sum sfr/sfr(z=0):', format(self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +3], '0.2e'))
        print('volume [Mpc3]:', format(self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +4], '0.2e'))
        
    def SFR2Z(
            self,
            data,
            filename,
            SFH_key_count,
            error_count,
            mycond_array,
            histo_data_min=0.0,
            histo_data_max=0.0,
            mycond_min=0.0,
            mycond_max=0.0,
            data_offset=0,
            histo_output_offset=0,
            log10_xaxis=False,
            log10_yaxis=False,
            norm_z0=False,
            SFH_key='',
            method='',
            env_name='',
            env_code=0):
                    
                  
        print('SFR vs Redshift(): --> Selection:', SFH_key,  'method:', method, 'self.i:', self.i, 'self.a:', self.a,'\n++++++++++++++++++++++++++++++++++++++++++++++++++')  

        def get_output_properties_string(props):
            """Produced a string of various galaxy and halo properties which can be written to a textfile. Use that if you want to 'append text'
            to an already exciting file"""
            
            prop_dict = mL.get_property_dict()
            
            format_string=''
            header_string=''
            
            for i, item in enumerate(props):                
                
                header_string+='('+str(i+1)+')'+item+'['+prop_dict[item]['unit']+'] '
                format_string+=prop_dict[item]['format']
                if i!=len(myprops)-1:
                    format_string+=' '

            return format_string, header_string

        def write_progenitors(data,
                              filename):

            myOutput.writeIntoFile(filename,
                       data[['parentIndex','nodeIndex']],
                       myheader= 'z='+str(format(redshift, '0.2f'))+'\n(1)parentIndex (2)nodeIndex',
                       data_format="%i\t%i",
                       mydelimiter='\t')
            
        def write_progenitors_300(data,
                              filename):

            myOutput.writeIntoFile(filename,
                       data[['parentIndex','nodeIndex','regionID']],
                       myheader= 'z='+str(format(redshift, '0.2f'))+'\n(1)parentIndex (2)nodeIndex (3)regionID',
                       data_format="%i\t%i\t%i",
                       mydelimiter='\t')
            
        def write_properties(data,
                             filename,
                             props):
            
            myformat_string, myheader_string = get_output_properties_string(props)
            
            #print(myformat_string, r'\n', myheader_string)           

            myOutput.writeIntoFile(filename[:-4]+'_props.txt',
                   data[props],
                   myheader= 'MDPL2 Galacticus TTH z='+str(format(redshift, '0.2f'))+'\n'+myheader_string,
                   data_format=myformat_string,
                   mydelimiter=' ')  
                        
        def select_300Clusters():
           
            import pandas as pd
                       
            self.myData2Process = self.myData2Process[np.where(self.myData2Process['orphan']!=2)[:][0]]
            print('after central selection:' ,self.myData2Process.shape)
            #self.myData2Process[::-1].sort(order=['mhalo'], axis=0)
            #print(self.myData2Process['haloid'][0:5])

            #print(self.myData2Process[0:2])

            # for prop in self.myData2Process.dtype.names:
            #     prop_dtype = self.myData2Process[prop].dtype
            #     print(prop, prop_dtype, '[0]:', self.myData2Process[prop][0] )
            #exit()
            
            
            #filename_parent_sample=mycomp+'/anaconda/pro/data/300/MDPL2_Galacticus_z0.00_cluster_regions_haloids.txt'
            filename_parent_sample=mycomp+'/anaconda/pro/data/300/MDPL_Galacticus_324_CCGs.txt' 
            #filename_parent_sample=mycomp+'/anaconda/pro/data/300/MDPL2_Galacticus_z0.00_cluster_regions_haloids_local_test.txt'
            #(1)regionID[ID] (2)haloid[ID] (3)nsats[Mstar>1e9Msun,count] (4)x_pos[comMpc] (5)y_pos[comMpc] (6)z_pos[comMpc] (7)mhalo[Msun] (8)mstar[Msun]            
            mydata_names_props= ['regionID','haloid','nsats','x_pos','y_pos','z_pos','mhalo', 'mstar']
                        
            parent_data = mL.df_to_sarray(pd.read_csv(filename_parent_sample, skiprows=2, names=mydata_names_props, sep='  '))
            #print(parent_data)
            #test_haloids = [12569515903, 12569210835, 12569696351, 12569842360, 12569770239]
            indices = np.in1d(self.myData2Process['haloid'], parent_data['haloid'])

            self.myData2Process = self.myData2Process[np.where(indices==True)[:][0]]
            #print('here 2777:', self.myData2Process.size)
            
            parentIndices, index1, index2 = np.intersect1d(self.myData2Process['haloid'], parent_data['haloid'], return_indices=True) 
            
            self.myData2Process[['regionID','nsats']]=-99
            for prop in ['regionID','nsats']:
                #print('property:', prop'
                self.myData2Process[prop][index1]=parent_data[prop][index2]
            
            
            self.myData2Process[::-1].sort(order=['haloid','mhalo'], axis=0)
            #print(self.myData2Process[['haloid','mhalo','orphan']][0:10])            

            data_2nd_CGs = mL.get_second_MGC(self.myData2Process, parent_data['haloid'])        
            
            print('shape of data after 300 selection:', self.myData2Process.shape)
                        
            return self.myData2Process, data_2nd_CGs

        def select_CMASS_z056():
            if self.myconfig_array['catname'+str(self.a)].find('Gala')!=-1:
                #CMASS_IDs = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/data/Galacticus_400Mpc/Galacticus_400Mpc_z_0.55_CMASS_down_sample3_haloid.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.uint64, skiprow=2)                         
                #CMASS_IDs = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/data/Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3_haloid.txt', data_format='ASCII', data_shape='shaped', delim='  ', mydtype=np.uint64, skiprow=2)
                #CMASS_IDs = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/data/Galacticus_1Gpc/Galacticus_1Gpc_z_0.09_tarsel_CUT3_Contreras+13_mcold_haloid.txt', data_format='ASCII', data_shape='shaped', delim='  ', mydtype=np.uint64, skiprow=2)
                #print(CMASS_IDs
                #exit()
                
                import pandas as pd
                
                filename_CMASS=mycomp+'anaconda/pro/data/Galacticus_1Gpc/SFH/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3.txt'
                mydata_names_props= ['haloid','hostid','mstar','mhalo','orphan','mcold','Mzgas','Mzstar','mbh','mstar_spheroid','mstar_disk','mcold_disk',\
                                      'mcold_spheroid','mhot','sfr_spheroid','sfr_disk','sfr','Mzgas_spheroid','Mzgas_disk','Mzstar_spheroid','Mzstar_disk',\
                                      'Mzhot_halo','x_pos','y_pos','z_pos','x_vel','y_vel','z_vel','L_SDSS_dA_total_u','L_SDSS_dA_total_g','L_SDSS_dA_total_r',\
                                      'L_SDSS_dA_total_i','L_SDSS_dA_total_z','rdisk','rbulge','rhalf_mass','spinParameter','MAB_dA_total_u','MAB_dA_total_g',\
                                      'MAB_dA_total_r','MAB_dA_total_i','MAB_dA_total_z','mAB_dA_total_u','mAB_dA_total_g','mAB_dA_total_r','mAB_dA_total_i',\
                                      'mAB_dA_total_z','mAB_dA_total_cut_r_i','mAB_dA_total_cut_dmesa','mAB_dA_total_cut_i_lt_dmesa','mhalo_sat','nodeIndex',\
                                      'parentIndex','satelliteNodeIndex','satelliteIndex','siblingIndex','satelliteMergeTime','isolated','timeLastIsolated',\
                                      'mhalo_200c','zcold','cgf','mbasic','mbasic_200c','NFW_con','ssfr','weight_tot','env_512','env_1024','pop']
                
                CMASS_data = mL.df_to_sarray(pd.read_csv(filename_CMASS, skiprows=2, names=mydata_names_props, sep='  ', engine='python'))
                #print(np.info(CMASS_data)
                indices = np.in1d(self.myData2Process['hostid'], CMASS_data['hostid'])
                #print(np.info(self.myData2Process)
                self.myData2Process = self.myData2Process[np.where(indices==True)[:][0]]
                print('shape of data after CMASS_ID selection:', self.myData2Process.shape, end='')
                
                #Select only centrals
                self.myData2Process = self.myData2Process[np.where(self.myData2Process['orphan']==0)[:][0]]
                print('only cents! final ngal:', self.myData2Process.size) 
    
                print('method: ', method, '--> select samples and follow progenitors! SFH_key:', SFH_key)
                if method=='M2':
                    data_sample=mL.sample_selection(self.myData2Process, sample_key=SFH_key)
                    
                    if env_name!='':
                        CMASS_IDs_env=CMASS_data['hostid'][np.where(CMASS_data['env_1024']==env_code)[0][:]]
                        
                        indices_env = np.in1d(data_sample['hostid'], CMASS_IDs_env)

                        data_sample = data_sample[np.where(indices_env==True)[:][0]]
                        
                        print('shape of data after environment! Name:', env_name, 'Code:', env_code, ' selection:', data_sample.shape)
                              
                    write_properties(data_sample,                                     
                                     mycomp+'/anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index_'+str(format(redshift, '0.2f'))+'_'+SFH_key+'_props_'+method+env_space+env_name+'.txt')

                    
                    return data_sample
                
                else:
                     write_properties(self.myData2Process,
                                      mycomp+'/anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index_'+str(format(redshift, '0.2f'))+'_props_'+method+env_space+env_name+'.txt')
                 
                     return self.myData2Process               

            else:
                self.myData2Process = self.myData2Process[np.where(self.myData2Process['hostid']>-2)[:][0]]
                
                return self.myData2Process

        def calc_extra_props():
            print('calc_estra_props --->', end=' ')
            
            self.myData2Process = rcfuncs.append_fields([self.myData2Process], ['color', 'z'] ,[self.myData2Process[mag_prefix+'i'],self.myData2Process[mag_prefix+'i']], usemask=False)
 
            print('expand array!', end=' ')               

            self.myData2Process['color']=self.myData2Process['MAB_dA_total_g']-self.myData2Process['MAB_dA_total_i']
                
            self.myData2Process['z']=redshift

            print('--> DONE!')


####### CONFIGURATIONS    ##################################################################    

        #print(self.myData2Process.shape
        #print(np.info(self.myData2Process)
        redshift = self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z']
        ext_redshift = '0.0'
        mag_prefix='mAB_dA_total_'
        mag_prefix_rest='mAB_rst_total_'
        
        run_interrupted=True
        
        if env_name!='':
            env_space='_'
        else:
            env_space=''
            
        
        import numpy.lib.recfunctions as rcfuncs
        from cosmolopy import cparam, cd, cc
        fidcosmo = cparam.PlanckMD(flat=True, extras=True)
        #print('lum_dist:', cd.luminosity_distance(redshift, **fidcosmo), 'z=', redshift)
        print('SFH_key:', SFH_key, 'method:', method, 'SFH_key_count', SFH_key_count, end=' ')
        
        if method=='M2':
            filename_indices=mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index'+SFH_key+env_space+env_name+'.txt'
            filename_properties_redshift=mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index_z_'+str(format(redshift, '0.2f'))+'_'+SFH_key+env_space+env_name+'.txt'

        elif method=='300':
            filename_indices=mycomp+'anaconda/pro/data/300/SFH_Indices/SFH_Index_300'+SFH_key+env_space+env_name+'.txt'
            filename_indices_redshift=mycomp+'anaconda/pro/data/300/SFH_Properties/SFH_Properties_300_z_'+str(format(redshift, '0.2f'))+'_'+SFH_key+env_space+env_name+'.txt'             
        else:
            filename_indices=mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index.txt'
            
        print('\n++++++++++++++++++++++++++++++++++++++++++++++++++\nfilename_indices:', filename_indices, '\n++++++++++++++++++++++++++++++++++++++++++++++++++\n')      
                
        calc_extra_props()
        z_test = str(format(redshift, '0.2f'))

####### SAMPLE SELECTION     ##################################################################

        if method=='300':
            
            myprops=['regionID', 'haloid', 'hostid', 'parentIndex', 'orphan', 'nsats', 'x_pos','y_pos','z_pos', 'mhalo', 'mstar', 'SHMR', 'mcold', 'mhot', 'mbh', 'mhot_outflow', 'Mzgas', 'Mzstar', 'Mzhot_halo', 'Mzhot_outflowHalo',\
                     'zcold', 'zstar', 'sfr', 'ssfr', 'Tcons', 'bheff','rhalf_mass', 'rhotHalo', 'cgf', 'jbar', 'jhotHalo', 'joutHotHalo', \
                     'L_SDSS_dA_total_g','L_SDSS_dA_total_r','L_SDSS_dA_total_i', 'mAB_dA_total_g', 'mAB_dA_total_r', 'mAB_dA_total_i',\
                     'MAB_dA_total_g', 'MAB_dA_total_r', 'MAB_dA_total_i', 'mAB_dA_total_cut_r_i', 'mAB_dA_total_cut_g_r', 'mAB_dA_total_cut_g_i']

            if self.i==0 and run_interrupted==False:            
            
                print('i:', self.i, 'z=', redshift, 'inital sample selection for 300-Clusters!')
                
                myfiltered_data, myfiltered_data_2nd_CGs = select_300Clusters()            

                #write properties off all galaxies    
                write_properties(myfiltered_data_2nd_CGs, filename_indices_redshift, myprops)
                
                myfiltered_data = myfiltered_data[np.where(myfiltered_data['orphan']==0)[:][0]]
                
                #Write progenitors for the inital selection (because in sample_selection() the keyword "sample_key"!= 'main').
                write_progenitors_300(myfiltered_data,
                                      filename_indices)

            # elif self.i==0 and run_interrupted==True:
                
            #     print('i:', self.i, 'z=', redshift, 'continue after interruption!')
                    
            #     myfiltered_data, myfiltered_data_2nd_CGs =mL.sample_selection(self.myData2Process,
            #                                                                     sample_key='main',
            #                                                                     filename_indices=filename_indices,
            #                                                                     redshift=redshift,
            #                                                                     include_2nd_massive_CGs=True)
            #     #write properties off all galaxies    
            #     write_properties(myfiltered_data_2nd_CGs, filename_indices_redshift, myprops)
                

            #     write_progenitors_300(myfiltered_data,
            #                           filename_indices)
                               
            else:
                print('i:', self.i, 'z=', format(redshift, '0.2f'), 'follow merger trees ...')
                
                myfiltered_data, myfiltered_data_2nd_CGs =mL.sample_selection(self.myData2Process,
                                                                                sample_key='main',
                                                                                filename_indices=filename_indices,
                                                                                redshift=redshift,
                                                                                include_2nd_massive_CGs=True)
                
                #write properties off all galaxies    
                write_properties(myfiltered_data_2nd_CGs, filename_indices_redshift, myprops)
                                                                    
                write_progenitors_300(myfiltered_data,
                                      filename_indices)
                                   
                              
        
        elif self.i==0 and SFH_key_count==0 and method!='M2':
                        
            myfiltered_data = select_CMASS_z056()            

            if self.myconfig_array['catname'+str(self.a)].find('Gala')!=-1:
                #Write progenitors for the inital selection (because in sample_selection() the keyword "sample_key"!= 'main').              
                write_progenitors(myfiltered_data,
                                  filename_indices)
       
        elif method=='M2':
            if self.i>0:
                print('i>0:')

                if z_test==ext_redshift and (SFH_key.find('st')!=-1 or SFH_key.find('gt')!=-1):
                    myfiltered_data=mL.sample_selection(self.myData2Process,
                                                        sample_key=SFH_key,
                                                        redshift=redshift)
                    #Write progenitors for the inital selection (because in sample_selection() the keyword "sample_key"!= 'main').    
                    write_progenitors(myfiltered_data,
                                      filename_indices)                        
                        
                else:
                    if self.myconfig_array['catname'+str(self.a)].find('Gala')!=-1:
                        myfiltered_data=mL.sample_selection(self.myData2Process,
                                                            sample_key='main',
                                                            filename_indices=filename_indices,
                                                            redshift=redshift)
                    else:
                        myfiltered_data=mL.sample_selection(self.myData2Process,
                                                            sample_key='valid',
                                                            filename_indices=filename_indices,
                                                            redshift=redshift)
                        
                               
            else:
                print('i=0', self.i)
                #Write progenitors for the inital selection (because in sample_selection() the keyword "sample_key"!= 'main').
                myfiltered_data = select_CMASS_z056()
                write_progenitors(myfiltered_data,
                                  filename_indices)                
            
        elif self.i>0 and SFH_key_count==0:
            if self.myconfig_array['catname'+str(self.a)].find('Gala')!=-1:
                myfiltered_data=mL.sample_selection(self.myData2Process,
                                                 sample_key='main',
                                                 filename_indices=filename_indices,
                                                 redshift=redshift)
                 
            else:
                myfiltered_data=mL.sample_selection(self.myData2Process,
                                                 sample_key='valid',
                                                 filename_indices=filename_indices,
                                                 redshift=redshift)
          


        print('after selection:', myfiltered_data.shape, end='')
        if method=='M1':
                     
            myfiltered_data=mL.sample_selection(self.myData2Process, 
                                             sample_key=SFH_key,
                                             redshift=redshift)                
            print('ngal after selection:', myfiltered_data.shape,'\n')#, 'check original:', self.myData2Process.shape
         
        
####### CALUCLATE PROPERTIES & RETURN TO MAIN   ##################################################################
        percent_low = 16.0
        percent_high = 84.0

        #print(np.info(data))
        #try:
        self.histo_data_ZvsSFR = data
   
        self.histo_data_ZvsSFR[self.i, data_offset*(self.a)]    = redshift

        #lookback time
        self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +1] = cd.lookback_time(self.histo_data_ZvsSFR[self.i, data_offset*(self.a)], z0=self.histo_data_ZvsSFR[0, data_offset*(self.a)], **fidcosmo)/3.1536e+16 #seconds to Gyrs

        data_cents=myfiltered_data[np.where(myfiltered_data['orphan']==0)[0][:]]    
        data_sats=myfiltered_data[np.where(myfiltered_data['orphan']==1)[0][:]]
        data_os=myfiltered_data[np.where(myfiltered_data['orphan']==2)[0][:]]

        #galaxy type count
        self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +2] = data_cents.size
        self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +3] = data_sats.size
        self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +4] = data_os.size

        #number density    
        self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +5] = myfiltered_data.size/float(self.volume)/1e-4
        try:
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +6] = np.nansum(myfiltered_data['sfr'])
        except:
            pass
        try:
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +7] = np.nansum(myfiltered_data['mstar'])
        except:
            pass
        
        try:
            #bluge to total mass stars (BvT)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +8]= np.nanmedian(myfiltered_data['mstar_spheroid']/myfiltered_data['mstar'])
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +9]= np.nanpercentile((myfiltered_data['mstar_spheroid']/myfiltered_data['mstar']),percent_low)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +10]= np.nanpercentile((myfiltered_data['mstar_spheroid']/myfiltered_data['mstar']),percent_high)        
        except:
            pass
        try:
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +11]= np.nanmedian(myfiltered_data['MAB_dA_total_g']-myfiltered_data['MAB_dA_total_i'])
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +12]= np.nanpercentile((myfiltered_data['MAB_dA_total_g']-myfiltered_data['MAB_dA_total_i']),percent_low)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +13]= np.nanpercentile((myfiltered_data['MAB_dA_total_g']-myfiltered_data['MAB_dA_total_i']),percent_high)
        except:
            print('g-i color calculation faild!')

        try:
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +14]= np.nanmedian(myfiltered_data['MAB_dA_total_g']-myfiltered_data['MAB_dA_total_r']) 
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +15]= np.nanpercentile((myfiltered_data['MAB_dA_total_g']-myfiltered_data['MAB_dA_total_r']),percent_low)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +16]= np.nanpercentile((myfiltered_data['MAB_dA_total_g']-myfiltered_data['MAB_dA_total_r']),percent_high)
        except:
            print('g-r color calculation faild!') 
 
        try:
            #black hole to halo mass ratio (black hole efficiency)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +17]= np.nanmedian(myfiltered_data['mbh']/myfiltered_data['mhalo'])
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +18]= np.nanpercentile(myfiltered_data['mbh']/myfiltered_data['mhalo'],percent_low)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +19]= np.nanpercentile(myfiltered_data['mbh']/myfiltered_data['mhalo'],percent_high)            
        except:
            pass
        
        try:
            #gas deplection time (Tcons) in Gyr
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +20]= np.nanmedian(myfiltered_data['mcold']/myfiltered_data['sfr']/1e9)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +21]= np.nanpercentile(myfiltered_data['mcold']/myfiltered_data['sfr']/1e9,percent_low)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +22]= np.nanpercentile(myfiltered_data['mcold']/myfiltered_data['sfr']/1e9,percent_high)            
        except:
            pass

        try:
            #mbar
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +23]= np.nanmedian(myfiltered_data['mstar']+myfiltered_data['mcold'])
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +24]= np.nanpercentile(myfiltered_data['mstar']+myfiltered_data['mcold'],percent_low)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +25]= np.nanpercentile(myfiltered_data['mstar']+myfiltered_data['mcold'],percent_high)            
        except:
            pass


        try:
            #jbar
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +26]= np.nanmedian((myfiltered_data['angM_disk']+myfiltered_data['angM_spheroid'])/(myfiltered_data['mstar']+myfiltered_data['mcold']))
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +27]= np.nanpercentile((myfiltered_data['angM_disk']+myfiltered_data['angM_spheroid'])/(myfiltered_data['mstar']+myfiltered_data['mcold']),percent_low)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +28]= np.nanpercentile((myfiltered_data['angM_disk']+myfiltered_data['angM_spheroid'])/(myfiltered_data['mstar']+myfiltered_data['mcold']),percent_high)            
        except:
            pass
        

        
        try:
            #jbulge
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +29]= np.nanmedian(myfiltered_data['angM_spheroid']/(myfiltered_data['mstar_spheroid']+myfiltered_data['mcold_spheroid']))
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +30]= np.nanpercentile(myfiltered_data['angM_spheroid']/(myfiltered_data['mstar_spheroid']+myfiltered_data['mcold_spheroid']),percent_low)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +31]= np.nanpercentile(myfiltered_data['angM_spheroid']/(myfiltered_data['mstar_spheroid']+myfiltered_data['mcold_spheroid']),percent_high)            
        except:
            pass
        
        try:
            #jdisk
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +32]= np.nanmedian(myfiltered_data['angM_disk']/(myfiltered_data['mstar_disk']+myfiltered_data['mcold_disk']))
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +33]= np.nanpercentile(myfiltered_data['angM_disk']/(myfiltered_data['mstar_disk']+myfiltered_data['mcold_disk']),percent_low)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +34]= np.nanpercentile(myfiltered_data['angM_disk']/(myfiltered_data['mstar_disk']+myfiltered_data['mcold_disk']),percent_high)            
        except:
            pass 

        try:
            #cold gas to stellar mass ratio (cgf)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +35]= np.nanmedian(myfiltered_data['mcold']/myfiltered_data['mstar'])
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +36]= np.nanpercentile(myfiltered_data['mcold']/myfiltered_data['mstar'],percent_low)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +37]= np.nanpercentile(myfiltered_data['mcold']/myfiltered_data['mstar'],percent_high)            
        except:
            pass
        
        try:
            #baryon fraction (fbar)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +38]= np.nanmedian(myfiltered_data['mcold']/(myfiltered_data['mcold']+myfiltered_data['mstar']))
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +39]= np.nanpercentile(myfiltered_data['mcold']/(myfiltered_data['mcold']+myfiltered_data['mstar']),percent_low)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +40]= np.nanpercentile(myfiltered_data['mcold']/(myfiltered_data['mcold']+myfiltered_data['mstar']),percent_high)            
        except:
            pass            

        try:
            #bluge to disk mass (BvT)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +41]= np.nanmedian(myfiltered_data['rbulge']/myfiltered_data['rdisk'])
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +42]= np.nanpercentile((myfiltered_data['rbulge']/myfiltered_data['rdisk']),percent_low)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +43]= np.nanpercentile((myfiltered_data['rbulge']/myfiltered_data['rdisk']),percent_high)        
        except:
            pass
        
        prop=['z','t','n_cents','n_sats','n_orphans','ndensity',\
              'sum_sfr','sum_mstar','BvT','g-i','g-r','bheff',\
              'Tcons','mbar','jbar','jbulge','jdisk','cgf','fbar','rbulgevsrdisk']
        
        add_props = ['sfr','ssfr','mstar','mhalo', 'mbh','SHMR','mcold','Mzgas','zcold','rdisk','rbulge','rhalf_mass',\
                     'mean_age_stars_disk','mean_age_stars_spheroid','vmax','vdisp','vdisk','vbulge',\
                     'L_SDSS_dA_total_g', 'L_SDSS_dA_total_r', 'L_SDSS_dA_total_i',\
                     mag_prefix+'g', mag_prefix+'r', mag_prefix+'i', 'MAB_dA_total_g', 'MAB_dA_total_r','MAB_dA_total_i'] 
            
        count=44
        for item in add_props:
            #print('count:', count, 'item:', item
            try: 
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +count]= np.nanmedian(myfiltered_data[item])
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +count+1]= np.nanpercentile(myfiltered_data[item],percent_low)
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +count+2]= np.nanpercentile(myfiltered_data[item],percent_high)            
            except:
                #print('--> not found!'
                pass
            count+=3

        if myfiltered_data.size>25:
            try:
                self.calcFastHisto(myfiltered_data['mstar'], filename[0:len(filename)-4]+'_mstar_histo.txt', 'mstar', 'Msun')
            except:
                pass
            try:
                self.calcFastHisto(myfiltered_data['sfr'], filename[0:len(filename)-4]+'_sfr_histo.txt', 'sfr', 'Msun yr-1')
            except:
                pass
            try:
                self.calcFastHisto(myfiltered_data['ssfr'], filename[0:len(filename)-4]+'_ssfr_histo.txt', 'ssfr', 'yr-1')
            except:
                pass
            try:
                self.calcFastHisto(myfiltered_data['rbulge']/myfiltered_data['rdisk'], filename[0:len(filename)-4]+'_rbulgevsrdisk_histo.txt', 'rbulgevsrdisk', '-')        
            except:
                pass
            try:
                self.calcFastHisto(myfiltered_data['rhalf_mass'], filename[0:len(filename)-4]+'_rhalf_mass_histo.txt', 'rhalf_mass', '-')        
            except:
                pass        
            try:
                self.calcFastHisto(myfiltered_data['mhalo'], filename[0:len(filename)-4]+'_mhalo_histo.txt', 'mhalo', 'Msun')        
            except:
                pass        
            try:
                self.calcFastHisto(myfiltered_data['MAB_dA_total_g']-myfiltered_data['MAB_dA_total_i'], filename[0:len(filename)-4]+'_g-i_histo.txt', 'g-i', '-', binning='lin')        
            except:
                pass
            try:
                self.calcFastHisto(myfiltered_data['zcold'], filename[0:len(filename)-4]+'_zcold_histo.txt', 'zcold', '-', binning='lin')        
            except:
                pass       
    
            try:     
                self.calcFastHisto(myfiltered_data[['mhalo','mstar','orphan']], filename[0:len(filename)-4]+'_SHMF_histo.txt', 'mhalo', '-', binup2D=True)        
            except:
                pass
            try:     
                self.calcFastHisto(myfiltered_data['mbh'], filename[0:len(filename)-4]+'_mbh_histo.txt', 'mbh', 'Msun')        
            except:
                pass        
            try:     
                self.calcFastHisto(myfiltered_data['mcold']+myfiltered_data['mstar'], filename[0:len(filename)-4]+'_mbar_histo.txt', 'mbar', 'Msun')        
            except:
                pass 
            try:     
                self.calcFastHisto(myfiltered_data['mean_age_stars_disk'], filename[0:len(filename)-4]+'_mean_age_stars_disk_histo.txt', 'mean_age_stars_disk', 'Gyr', binning='lin')        
            except:
                pass
            try:     
                self.calcFastHisto(myfiltered_data['mean_age_stars_spheroid'], filename[0:len(filename)-4]+'_mean_age_stars_spheroid_histo.txt', 'mean_age_stars_spheroid', 'Gyr', binning='lin')        
            except:
                pass
            try:     
                self.calcFastHisto(myfiltered_data['vmax'], filename[0:len(filename)-4]+'_vmax_histo.txt', 'vmax', 'kms-1')        
            except:
                pass
            try:     
                self.calcFastHisto(myfiltered_data['vdisp'], filename[0:len(filename)-4]+'_vdisp_histo.txt', 'vdisp', 'kms-1')        
            except:
                pass         
            try:     
                self.calcFastHisto(myfiltered_data['mcold']/myfiltered_data['mstar'], filename[0:len(filename)-4]+'_cgf_histo.txt', 'cgf', '-')        
            except:
                pass
            try:     
                self.calcFastHisto(myfiltered_data['mbh']/myfiltered_data['mhalo'], filename[0:len(filename)-4]+'_bheff_histo.txt', 'bheff', '-')        
            except:
                pass
            try:     
                self.calcFastHisto(myfiltered_data['mcold']/(myfiltered_data['mcold']+myfiltered_data['mstar']), filename[0:len(filename)-4]+'_fbar_histo.txt', 'fbar', '-')        
            except:
                pass         
            try:     
                self.calcFastHisto(myfiltered_data['mcold']/myfiltered_data['sfr']/1e9, filename[0:len(filename)-4]+'_Tcons_histo.txt', 'Tcons', 'Gyr-1')        
            except:
                pass
            try:     
                self.calcFastHisto((myfiltered_data['angM_disk']+myfiltered_data['angM_spheroid'])/(myfiltered_data['mstar']+myfiltered_data['mcold']), filename[0:len(filename)-4]+'_jbar_histo.txt', 'jbar', 'kpc kms-1')        
            except:
                pass
            try:     
                self.calcFastHisto(myfiltered_data['angM_disk']/(myfiltered_data['mstar_disk']+myfiltered_data['mcold_disk']), filename[0:len(filename)-4]+'_jdisk_histo.txt', 'jdisk', 'kpc kms-1')        
            except:
                pass         
            try:     
                self.calcFastHisto(myfiltered_data['angM_spheroid']/(myfiltered_data['mstar_spheroid']+data['mcold_spheroid']), filename[0:len(filename)-4]+'_jbulge_histo.txt', 'jbulge', 'kpc kms-1')        
            except:
                pass         
        #        
    #        try:
    #            self.HODFunction(myfiltered_data, filename, which_halomass='mhalo_200c')
    #        except:
    #            pass
       
            try:
                #print('mycond array:', mycond_array
                self.TwoPCF(myfiltered_data, filename, mycond_array)
            except:
                print('CLUSTERING CALCULACTION NOT POSSIBLE!')
            
        else:
            print('no halos found at this redshift OR data sample is too small!')
     
 
        #except:
            #print('no halos found at this redshift ...!')

        return error_count, prop+add_props, str(format(percent_high, '.0f')),  str(format(percent_low, '.0f'))
    
    def binAndFrac2D(
                    self,
                    data,
                    nbins,
                    name_x,
                    name_y=False,
                    histo_data_min='min',
                    histo_data_max='max',
                    mycond_min_x='min',
                    mycond_max_x='max',
                    mycond_min_y='min',
                    mycond_max_y='max',                    
                    data_offset=0,
                    histo_output_offset=0,
                    log10bin=True,
                    binUp2D=False,
                    div_y2x=False,
                    add_axis=False,
                    plot_key=False,
                    cumulative=False,
                    normalise=False,
                    starforming_cut=1e-11,
                    volume_units='Mpc-3',
                    weight_function=False):

        def caseSwitcher_plot_key(plot_key, myfiltered_data):
        
            choose = {
                'SMF':          SMF,
                'SFRF':         SFRF,
                'sSFRF':        sSFRF,
                'oh2mstar':     oh2mstar,
                'mstar2mhalo':  mstar2mhalo,
                'mstar2mhalo_no': mstar2mhalo,
                'mstar2mhalovsSFR': mstar2mhalo,
                'HMF':          HMF,
                'HMF_no':       HMF,
                'zgas2mstar':   zgas2mstar,
                'Mzgas2mcold':  Mzgas2mcold
                }
                
            func = choose.get(plot_key)
            
            return func(plot_key, myfiltered_data)

        def SMF(plot_key, myfiltered_data):

            if self.myconfig_array['catname'+str(self.a)]=='CMASS_SPALL':
                if self.myconfig_array['catname'+str(self.a)].find('SPALL')!=-1:
                    print('correct IMF --> Kroupa to Chabrier')
                    myfiltered_data['mstar']/=10**(0.03925)
                    myfiltered_data['mstar']/=10**(0.2)
    
                    #print('correct little-h for', self.myconfig_array['catname'+str(self.a)]
                    #correction_factor=0.7**2            #Msun --> Msunh-2 WMAP9
                    #correction_factor=(0.7/0.6777)**2    #Msunh-2 WMAP --> Msunh-2 Planck
                    #correction_factor/=0.6777           #Msunh-2 Planck --> Msunh-1 WMAP9                     
                    #self.myData2Process['mstar']*=(0.7/0.6777)
    
            try:
                myfiltered_data['mstar']=myfiltered_data['mstar_spheroid']+myfiltered_data['mstar_disk']
            except:
                pass
            
     #        myfiltered_data = myfiltered_data[np.argsort(myfiltered_data['ssfr'])]
    #        myfiltered_data=myfiltered_data[int(myfiltered_data['mstar'].size/100*40)::]                
    #        print('selection: --> highest', myfiltered_data.shape        
            
    #        Blue-Red separation
            #myfiltered_data=myfiltered_data[np.where(myfiltered_data['mAB_total_g']-myfiltered_data['mAB_total_i']>2.35)[0][:]] 
    #        Guo+13 cut
    #        (r-i)>0.679 -0.082(Mi+20)
#            cut=0.679-0.082*(myfiltered_data['MAB_dA_total_i']-5*np.log10(0.6777)+20.0)
#            myfiltered_data=myfiltered_data[np.where(myfiltered_data['mAB_dA_total_r']-myfiltered_data['mAB_dA_total_i']>cut)[0][:]]         
#            print(name_x, 'selection: --> r-i > ', myfiltered_data.shape
             
            #myfiltered_data=myfiltered_data[np.where(myfiltered_data['sfr']/myfiltered_data['mstar']>1e-11)[0][:]]
            #myfiltered_data=myfiltered_data[np.where(myfiltered_data['sfr']>1e-4)[0][:]]         
    #        
            #print(name_x, 'selection: --> g-i ', myfiltered_data.shape
            
            return myfiltered_data
    
        def SFRF(plot_key, myfiltered_data):
            return myfiltered_data
                       
            
        def sSFRF(plot_key, myfiltered_data):
            if np.all(myfiltered_data['ssfr']==-99.0):
                myfiltered_data['ssfr']=myfiltered_data['sfr']/myfiltered_data['mstar']
                
            return myfiltered_data 
 

        def ssfr2mstar(plot_key, myfiltered_data):
            
            myfiltered_data['sfr']/=myfiltered_data['mstar']           
            print('starforming cut:', starforming_cut)
            myfiltered_data = myData.selectData2Compute(myfiltered_data, 
                                                        selected_col='sfr', 
                                                        operator='>', 
                                                        condition=starforming_cut)

            return myfiltered_data
 
        def HMF(plot_key, myfiltered_data):

            print('HMF_no\n+++++++++++++++++')
            print('shape before non-orphan selection!', myfiltered_data.shape, end='')
         
            try:
                print('--> select non-orphans', end='')
                myfiltered_data = myData.selectData2Compute(myfiltered_data, 
                                                            selected_col='orphan', 
                                                            operator='<', 
                                                            condition=2)
                                                                                                   
                print('shape after non-isolated selection!', myfiltered_data.shape)                                           
            except:
                print('no col "orphan" found! Could not do selection ...')
                
            return myfiltered_data
                                                               
        def mstar2mhalo(plot_key, myfiltered_data):
            print('mstar2mhalo\n++++++++++++++++++')
            #approximated 200c overdensity converstion

            if plot_key.find('vsSFR')!=-1:
                myfiltered_data[name_y]=myfiltered_data['mstar']/myfiltered_data[name_y]
 
            return myfiltered_data                                                              

        def Mzgas2mcold(plot_key, myfiltered_data):
            print('Mzgas2mstar\n++++++++++++++++++')
            
            print('Zgas=8.69+log10(Mzgas/(Mcold*0.0134)) ref: Allende-Pieto+01, Asplund+09')
            if self.myconfig_array['catname'+str(self.a)].startswith('Galacticus') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5':
                #column Zgas_spheroid and Zgas_disk are wrongly named in the catalog it should be Mzgas_spheroid and Mzgas_disk
                print('Galaticus: zgas_disk+zgas_spheroid = Mzgas')
                myfiltered_data['zgas_spheroid']+=myfiltered_data['zgas_disk']

            myfiltered_data=myfiltered_data[np.where(myfiltered_data[add_axis]>0.0)[0][:]]

            myfiltered_data[name_y]=myfiltered_data[name_y]/(myfiltered_data[add_axis]*0.0134)                
            
            myfiltered_data=myfiltered_data[np.where(myfiltered_data[name_y]>0.0)[0][:]]
                                            
            myfiltered_data[name_y]=8.69+np.log10(myfiltered_data[name_y]) 

            myfiltered_data = myData.selectData2Compute(myfiltered_data, 
                                                        selected_col=name_y, 
                                                        operator='!=', 
                                                        condition='nan')

                        
            if name_x.find('mcold')!=-1:

                myfiltered_data=myfiltered_data[np.where(myfiltered_data[name_x]>1e3)[0][:]]                    
                
                
                myfiltered_data = mL.filter_data_before_analysis(myfiltered_data, 
                                                                   5e10, 
                                                                   1e11, 
                                                                   'mstar')                
                print('mstar selection: -->', myfiltered_data.shape, end='')
                
                myfiltered_data[name_x]/=myfiltered_data['mstar'] 

                myfiltered_data = mL.filter_data_before_analysis(myfiltered_data, 
                                                                   1e-4, 
                                                                   1e4, 
                                                                   name_x)
                print('Zcold selection: -->', myfiltered_data.shape, end='')                                                    
                myfiltered_data = mL.filter_data_before_analysis(myfiltered_data, 
                                                                   6, 
                                                                   11, 
                                                                   name_y)                                                                   
                print('Mcold/Mstar selection: -->', myfiltered_data.shape)


            return myfiltered_data                                                 



        def zgas2mstar(plot_key, data):
            print('\nzgas2mstar\n+++++++++++++++++')

            try:
                data['cgf']=data['mcold']/data['mstar']
                
                print('CGF min/max', format(np.log10(min(data['cgf'])), '0.3f'), '/', format(np.log10(max(data['cgf'])), '0.3f'))
            except:
                try:
                    data['cgf']=data['mcold_disk']/data['mstar']
                    
                    print('SAGE: CGF min/max', format(np.log10(min(data['cgf'])), '0.3f'), '/', format(np.log10(max(data['cgf'])), '0.3f'))
                except:
                    pass
            

            try:
                data['sfr']=data['sfr_disk']+data['sfr_spheroid']
                data['ssfr']=data['sfr']/(data['mstar'])
                print('DONE set ssfr!')
            except:
                try:
                    data['ssfr']=data['sfr']/data['mstar']                    
                except:
                    pass
            try:
                print('starforming cut?', starforming_cut, end='')
                
                data=data[np.where(data['ssfr']>starforming_cut)[:][0]]
                                                            
                print('min', min(data['ssfr']), 'max:', max(data['ssfr']), '-->',data.shape)
            except:
                print('starforming cut FAILED!') 

            try:
                data['Mzgas']=data['Mzgas_disk']+data['Mzgas_spheroid']
                print('Mzgas=Mzgas_disk+Mzgas_spheroid', '-->',data.shape)
            except:
                pass

            try:
                data=data[np.where(data['mcold']>0.0)[:][0]]
                print('mcold>0 -->',data.shape)
            except:
                pass
                  

            # Converstion value 0.0365507 comes from ...
            # (O/H)_gal / (O/H)_sun = Z_gal / Z_sun
            # (O/H)_gal = (O/H)_sun * Z_gal / Z_sun
            # (O/H)_sun = 10**(8.38-12) --> from 12+log10(0/H)_sun = 8.69 (Allende-Prieto+01)
            # Z_sun=0.0134 --> from (Asplund+09)
            # (O/H)_gal = M_metals / M_baryons --> MZ_gas total / Mcold gas                                                      
            # (O/H)_sun / Z_sun = 0.0365507
            # IF O/H available use left side of equation, if other properties available, the right side!

            # 12+np.log10(0.0365507*O/H) == 8.69+np.log10(Mzgas_cold/Mcold/0.0134)

            # myfiltered_data[name_y]=12+np.log10(0.0365507*myfiltered_data[name_y]) == myfiltered_data[name_y]=8.69+np.log10(myfiltered_data[add_axis]/myfiltered_data[name_y]/0.0134)

            try:
                data[name_y]=8.69+np.log10(data['Mzgas']/data['mcold']/0.0134)
            except:
                try:
                    print('calculation failed! Try with disk properties only!', end='')
                    data=data[np.where(data['mcold_disk']>0.0)[:][0]]
                    print('mcold>0 -->', data.shape)
                    data[name_y]=8.69+np.log10(data['Mzgas_disk']/data['mcold_disk']/0.0134)
                    print('--> DONE!')
                except:
                    print('--> could not calculate Zcold!')

                
            data=data[np.where(np.isfinite(data[name_y]))[:][0]]
            print( 'after isfinite -->', data.shape)
            try:
                data=data[np.where(np.isfinite(data[add_axis]))[:][0]]
            except:
                pass

#            try:
#                data = data[np.where(data['mstar_disk']/data['mstar']>0.7)[:][0]]
#                print('alternative starforming cut: mstar_disk/mstar>0.7 -->', data.size)
#            except:
#                pass
            
            #myfiltered_data[name_y]=12+np.log10(0.0365507*myfiltered_data['zcold'])
            print('Zcold min/max', min(data[name_y]), '/', max(data[name_y]), '-->', data.shape)
            return data                                               


                    
        def oh2mstar(plot_key, data):
            print('\noh2mstar\n++++++++++++++++++')

            try:
                data['cgf']=data['mcold']/data['mstar']
                
                print('CGF min/max', format(np.log10(min(data['cgf'])), '0.3f'), '/', format(np.log10(max(data['cgf'])), '0.3f'))
            except:
                try:
                    data['cgf']=data['mcold_disk']/data['mstar']
                    
                    print('SAGE: CGF min/max', format(np.log10(min(data['cgf'])), '0.3f'), '/', format(np.log10(max(data['cgf'])), '0.3f'))
                except:
                    pass

            
            try:
                data['sfr']=data['sfr_disk']+data['sfr_spheroid']
                data['ssfr']=data['sfr']/(data['mstar'])
                print('DONE set ssfr!')
            except:
                try:
                    data['ssfr']=data['sfr']/data['mstar']                    
                except:
                    pass
            try:
                print('starforming cut?', starforming_cut, end='')
                
                data=data[np.where(data['ssfr']>starforming_cut)[:][0]]
                                                            
                print('min', min(data['ssfr']), 'max:', max(data['ssfr']), '-->',data.shape)
            except:
                print('starforming cut FAILED!')

            try:
                data = data[np.where(data['mstar_disk']/data['mstar']>0.7)[:][0]]
                print('alternative starforming cut: mstar_disk/mstar>0.7 -->', data.size)
            except:
                pass 

           
            try:
                print('min', min(data['zgas_disk']), 'max:', max(data['zgas_disk']), '-->', data.shape)
                data=data[np.where(data['zgas_disk']>0.0)[:][0]] 
                data=data[np.where(data['zgas_disk']<=1.0)[:][0]] 
            except:
                data=data[np.where(data['OH_gas_disk_bulge']>0.0)[:][0]] 
                data=data[np.where(data['OH_gas_disk_bulge']<=1.0)[:][0]]
                
            print('check 0.0 > zgas >= 1.0 (its a metal abundance!) -->',data.shape)
            #data=data[np.where(data[str(add_axis)]>0.0)[:][0]]
            if self.myconfig_array['catname'+str(self.a)].startswith('SAGE'):
                print('SAGE!', end='')
                                                             
                data[name_y]=8.69+np.log10(data['zgas_disk']/0.0134)
            else:
                 
                data[name_y]=12.0+np.log10(data['OH_gas_disk_bulge']/16.0)
                

            data=data[np.where(np.isfinite(data[name_y]))[:][0]]
            
            print('oh min/max', min(data[name_y]), '/', max(data[name_y]))
            try:
                print('add axis min/max', add_axis, min(data[add_axis]), '/', max(data[add_axis]))
            except:
                pass
            return data

               
        print('\n\nbinAndFrac2D()', 'self.i:', self.i, 'self.a:', self.a, '\n++++++++++++++++++++++++++++++++++++++++++++++\n')

        self.ratio_data_binAndFrac2D = data        

        print('CONFIGURATIONS:')
        print('~~~~~~~~~~~~~~~~')
        print('PLOTKEY:', plot_key)
        print('binUp2D:', binUp2D)
        print('div_y2x:', div_y2x)
        print('add_axis:', add_axis)
        print('log10bin:', log10bin)
        print('cumulative:', cumulative)
        print('normalise:', normalise, '\n')
        
        print('name_x:', name_x, 'col id:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_'+name_x], 'mycond x-axis: min/max', mycond_min_x, '/', mycond_max_x)

        try:
            self.myData2Process['mstar']=self.myData2Process['mstar_spheroid']+self.myData2Process['mstar_disk']
        except:
            pass

#        for item in ['mstar', 'mstar_disk', 'mstar_spheroid']:
#            print('item:', item
#            self.myData2Process[item]=self.myData2Process[item]/0.6778

        #print(np.info(self.myData2Process)
        myfiltered_data = mL.filter_data_before_analysis(self.myData2Process, 
                                                           mycond_min_x, 
                                                           mycond_max_x, 
                                                           name_x)
        print('\t', name_x, 'selection: -->', myfiltered_data.shape)
                                                                                                                                           

        if name_y!=False:
            print('name_y:', name_y, 'col id:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_'+name_y], 'mycond y-axis: min/max', mycond_min_y, '/', mycond_max_y)
            myfiltered_data = mL.filter_data_before_analysis(myfiltered_data, 
                                                               mycond_min_y, 
                                                               mycond_max_y, 
                                                               name_y)                                                       
            print('\t', name_y, 'selection: -->', myfiltered_data.shape)

        try: 
            print('add_axis:', add_axis, 'col id', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_'+add_axis], '\n\t', add_axis, 'selection: -->', myfiltered_data.shape)
        except:
            print('\tno axis to add ... ')


        if plot_key.find('_no')!=-1 or plot_key.find('vsSFR')!=-1:                          
            print('exclude orphan galaxies:', end='')
            myfiltered_data=myfiltered_data[np.where(myfiltered_data['orphan']==0)[0][:]] 
            print('\tselection: -->', myfiltered_data.shape)
           
        
        data2bin = caseSwitcher_plot_key(plot_key, myfiltered_data)
    
        print('check data:\n+++++++++++++\nname_x:', name_x ,'min/max:', min(data2bin[name_x]), '/', max(data2bin[name_x]))
        if name_y!=False: 
            print('name_y:', name_y ,'min/max:', min(data2bin[name_y]), '/', max(data2bin[name_y]))
        try:                
            print('add_axis:', add_axis,'min/max:', min(data2bin[add_axis]), '/', max(data2bin[add_axis]))
        except:
            pass

        def calcHisto(mydata):
            #rint 'calculate!'
            return myFuncs.binUp(mydata,
                                nbins,
                                histo_min=histo_data_min,
                                histo_max=histo_data_max,
                                log10bin=log10bin,
                                binup_2D=binUp2D,
                                div_y2x=div_y2x,
                                use_MAD=True,
                                add_axis=add_axis,
                                cumulative=cumulative,
                                weights=weight_function)
            
        if float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'])==0.0:
            weight_function=True
            print('weight function --> OBERVATIONAL DATA SET [name_x,weight] --> col name weights:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_name_weights'])
            histo, binsize = calcHisto(data2bin[[name_x,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_name_weights']]])          
        else:   
                
            try:
                print('[name_x,name_y,add_axis] -->', end='')
                histo, binsize = calcHisto(data2bin[[name_x,name_y,add_axis]])
            except:
                try:
                    print('[name_x,name_y] -->', end='')
                    histo, binsize = calcHisto(data2bin[[name_x,name_y]])
                except:
                    print('[name_x] -->', end='')
                    histo, binsize = calcHisto(data2bin[name_x])

        if float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'])==0.0:
            print('... calculate survey volume. skycoverage:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_skycoverage'], '[deg2]', 'zmax:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_zmax'], 'zmin:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_zmin'], '-->', end='')
            mynorm_y, myunit_volume = mL.calc_survey_volume(float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_skycoverage']),
                                             float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_zmin']),
                                             float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_zmax']),
                                             little_h_out=False)
            print('volume:', mynorm_y, myunit_volume)
            
        elif volume_units.find('h3')!=-1:
            print(volume_units, '--> found h3!')
            mynorm_y=(float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size']))**3
        else:
            print(volume_units, '--> default!')
            mynorm_y=(float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'])/float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_hubble_par']))**3                            

        if normalise==True:
            print(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'], self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_hubble_par'])
            print('normalise histogram: by', (float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'])/float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_hubble_par']))**3)                                                                                                                   
                    
            self.ratio_data_binAndFrac2D[:,data_offset*self.a:data_offset*self.a+data_offset] = myFuncs.normaliseHisto(histo[:,[0,1,2,3,4,5,6]],
                                                                                                                       norm_y=mynorm_y,
                                                                                                                       x_error=binsize)
                                                                                    
        else:
            self.ratio_data_binAndFrac2D[:,data_offset*self.a:data_offset*self.a+data_offset]=histo                                                                                                      
   
               
    def TwoPCF(self,
               mydata,
               filename,
               myconds_array):
        """https://pypi.python.org/pypi/Corrfunc or https://manodeep.github.io/Corrfunc/ by Manodeep Sinha (manodeep@gmail.com) version: 1.1.0 (2016-06-08)

            Clustering Measures on a Cosmological box:
            -------------------------------------------
            
                All codes that work on cosmological boxes with co-moving positions are located in the xi_theory directory. The various clustering measures are:
        
                xi_of_r     – Measures auto/cross-correlations between two boxes. The boxes do not need to be cubes.
                xi          – Measures 3-d auto-correlation in a cubic cosmological box. Assumes PERIODIC boundary conditions.
                wp          – Measures auto 2-d point projected correlation function in a cubic cosmological box. Assumes PERIODIC boundary conditions.
                xi_rp_pi    – Measures the auto/cross correlation function between two boxes. The boxes do not need to be cubes.
                vpf         – Measures the void probability function + counts-in-cells.              


            Clustering measures on a Mock:
            --------------------------------
            
                All codes that work on mock catalogs (RA, DEC, CZ) are located in the xi_mocks directory. The various clustering measures are:
        
                DDrppi      – The standard auto/cross correlation between two data sets. The outputs, DD, DR and RR can be combined using wprp to produce the Landy-Szalay estimator for \(w_p(r_p)\).
                wtheta      – Computes angular correlation function between two data sets. The outputs from DDtheta_mocks need to be combined with wtheta to get the full \(\omega(\theta)\)
                vpf         – Computes the void probability function on mocks.       
                            
            Science options:
            ----------------------
            
                If you plan to use the command-line, then you will have to specify the code runtime options at compile-time. For theory routines, these options are in the file theory.options while for the mocks, these options are in file mocks.options.
            
                Note All options can be specified at runtime if you use the python interface or the static libraries. Each one of the following Makefile option has a corresponding entry for the runtime libraries.
                Theory (in theory.options)
            
                PERIODIC (ignored in case of wp/xi) -- switches periodic boundary conditions on/off. Enabled by default.           
                OUTPUT_RPAVG                        -- switches on output of <rp> in each rp bin. Can be a massive performance hit (~ 2.2x in case of wp). Disabled by default.            
                DOUBLE_PREC                         -- switches on calculations in double precision. Disabled by default (i.e., calculations are performed in single precision by default).
            
            Mocks (in mocks.options)
            
                OUTPUT_RPAVG    -- switches on output of <rp> in each rp bin for DDrppi_mocks. Enabled by default.            
                OUTPUT_THETAAVG -- switches on output of in each theta bin. Can be extremely slow (~5x) depending on compiler, and CPU capabilities. Disabled by default.            
                DOUBLE_PREC     -- switches on calculations in double precision. Disabled by default (i.e., calculations are performed in single precision by default).          
                LINK_IN_DEC     -- creates binning in declination for DDtheta. Please check that for your desired limits \theta, this binning does not produce incorrect results (due to numerical precision). Generally speaking, if your \thetamax (the max. \theta to consider pairs within) is too small (probaly less than 1 degree), then you should check with and without this option. Errors are typically sub-percent level.           
                LINK_IN_RA      -- creates binning in RA once binning in DEC has been enabled. Same numerical issues as LINK_IN_DEC            
                FAST_DIVIDE     -- Disabled by default. Divisions are slow but required DD(r_p,\pi). Enabling this option, replaces the divisions with a reciprocal followed by a Newton-Raphson. The code will run ~20% faster at the expense of some numerical precision. Please check that the loss of precision is not important for your use-case.            
                FAST_ACOS       -- Relevant only when OUTPUT_THETAAVG is enabled. Disabled by default. An arccos is required to calculate <\theta>. In absence of vectorized arccos (intel compiler, icc provides one via intel Short Vector Math Library), this calculation is extremely slow. However, we can approximate arccos using polynomials (with Remez Algorithm <https://en.wikipedia.org/wiki/Remez_algorithm>). The approximations are taken from implementations released by Geometric Tools <http://geometrictools.com/>. Depending on the level of accuracy desired, this implementation of fast acos can be tweaked in the file utils/fast_acos.h <utils/fast_acos.h>__. An alternate, less accurate implementation is already present in that file. Please check that the loss of precision is not important for your use-case.            
                COMOVING_DIST   -- Currently there is no support in Corrfunc for different cosmologies. However, for the mocks routines like, DDrppi_mocks and vpf_mocks, cosmology parameters are required to convert between redshift and co-moving distance. Both DDrppi_mocks and vpf_mocks expects to receive a redshift array as input; however, with this option enabled, the redshift array will be assumed to contain already converted co-moving distances. So, if you have redshifts and want to use an arbitrary cosmology, then convert the redshifts into co-moving distances, enable this option, and pass the co-moving distance array into the routines.
            """

        #Clustering 3D: Auto-correlation on periodic, cosmological boxes using the Natural Estimator: 
        #Pairs which are separated by less than the r bins in 3-D real space.
        def xi():       return cf.theory.xi(boxsize, nthreads, rbins, data['x_pos'], data['y_pos'], data['z_pos'], output_ravg=True) 

        #PAIR COUNT: Calculate the 3-D pair-counts corresponding to the real-space correlation function, ξ(r).
        #--> convert to ξ(r) via cf.utils.convert_3d_counts_to_cf
        def xi_of_r():  return cf.theory.DD(1, nthreads, rbins, data['x_pos'],data['y_pos'], data['z_pos'], boxsize=boxsize,  periodic=True, output_ravg=True) 



        def DD2xi(): 

            rand_X, rand_Y, rand_Z = create_random_positions(N, rand_N)
            DD_counts = cf.theory.DD(1, nthreads, rbins, data['x_pos'], data['y_pos'], data['z_pos'], boxsize=boxsize)
            DR_counts = cf.theory.DD(0, nthreads, rbins, data['x_pos'], data['y_pos'], data['z_pos'], X2=rand_X, Y2=rand_Y, Z2=rand_Z, boxsize=boxsize)
            RR_counts = cf.theory.DD(1, nthreads, rbins, rand_X, rand_Y, rand_Z, boxsize=boxsize)
            
            #Converts raw pair counts to a correlation function.  
            return cf.utils.convert_3d_counts_to_cf(N, N, rand_N, rand_N,
                                                    DD_counts, DR_counts, DR_counts, RR_counts) 

        #projected correlation function wp(rp) in a periodic cosmological box. Pairs which are separated by less than the rp bins in the X-Y plane,
        #and less than pimax in the Z-dimension are counted.
        def wp(): return cf.theory.wp(boxsize, pimax, nthreads, rbins, data['x_pos'], data['y_pos'], data['z_pos'], output_rpavg=True) 

        
        #PAIR COUNT in 2D for (auto or cross) correlations for a periodic COSMOLOGICAL BOX with separation rp
        #--> convert to ξ(rp,π) via cf.utils.convert_3d_counts_to_cf or directly use function DDrppi2xi OR
        #--> convert to wp(rp) via cf.utils.convert_rp_pi_counts_to_wp or directly use function DDrppi2wp
        def xi_rp_pi(): return cf.theory.DDrppi(1, nthreads, pimax, rbins, data['x_pos'], data['y_pos'], data['z_pos'], boxsize=boxsize,  periodic=True, output_rpavg=True)

        def create_DDrppi_pair_counts():
            #Calculate the 3-D pair-counts corresponding to the real-space correlation function, ξ(rp,π). Pairs which are separated by less than the rp bins in the X-Y plane,
            #and less than pimax in the Z-dimension are counted.
            rand_X, rand_Y, rand_Z = create_random_positions(N, rand_N)

            DD_counts = cf.theory.DDrppi(1, nthreads, pimax, rbins, data['x_pos'], data['y_pos'], data['z_pos'], boxsize=boxsize, periodic=True, output_rpavg=True)
            DR_counts = cf.theory.DDrppi(0, nthreads, pimax, rbins, data['x_pos'], data['y_pos'], data['z_pos'], X2=rand_X, Y2=rand_Y, Z2=rand_Z, boxsize=boxsize,  periodic=True, output_rpavg=True)
            RR_counts = cf.theory.DDrppi(1, nthreads, pimax, rbins, rand_X, rand_Y, rand_Z, boxsize=boxsize, periodic=True, output_rpavg=True)  
            
            return DD_counts, DR_counts, RR_counts

        def DDrppi2xi():  

            DD_counts, DR_counts, RR_counts = create_DDrppi_pair_counts()
            
            #Converting (rp,π) pairs into a real-space correlation function ξ(rp,π)
            return cf.utils.convert_3d_counts_to_cf(N, N, rand_N, rand_N,
                                                    DD_counts, DR_counts, DR_counts, RR_counts)  


        def DDrppi2wp():  

            DD_counts, DR_counts, RR_counts = create_DDrppi_pair_counts()
            
            #Converting (rp,π) pairs into a projected correlation function wp(rp)
            return cf.utils.convert_rp_pi_counts_to_wp(N, N, rand_N, rand_N,
                                                          DD_counts, DR_counts, DR_counts, RR_counts, 
                                                          nbins, pimax) 


        #Pair counts ONLY in 2D real space for (auto or cross) correlations for a periodic COSMOLOGICAL BOX with seperation s
        #returns pair counts and not the actual correlation function ξ(s,μ)
        #--> use cf.utils.convert_3d_counts_to_cf for computing ξ(s,μ) from the pair counts.
        def xi_s_pi(): return cf.theory.DDsmu(1, nthreads, rbins, mu_max, nbins_mu, data['x_pos'],data['y_pos'], data['z_pos'], output_savg=True, boxsize=boxsize, periodic=True)         

        def create_DDsmu_pair_counts():
            #Calculate the 2-D pair-counts corresponding to the redshift-space correlation function, ξ(s,μ) Pairs which are separated by less than the s bins (see xi_s_pi)

            rand_X, rand_Y, rand_Z = create_random_positions(N, rand_N)

            DD_counts = cf.theory.DDsmu(1, nthreads, rbins, mu_max, nbins_mu, data['x_pos'],data['y_pos'], data['z_pos'], output_savg=True, boxsize=boxsize, periodic=True)
            DR_counts = cf.theory.DDsmu(0, nthreads, rbins, mu_max, nbins_mu, data['x_pos'], data['y_pos'], data['z_pos'], X2=rand_X, Y2=rand_Y, Z2=rand_Z, output_savg=True, boxsize=boxsize, periodic=True)
            RR_counts = cf.theory.DDsmu(1, nthreads, rbins, mu_max, nbins_mu, rand_X, rand_Y, rand_Z, output_savg=True, boxsize=boxsize, periodic=True)
            
            return DD_counts, DR_counts, RR_counts

        def DDsmu2xi():
            
            DD_counts, DR_counts, RR_counts = create_DDsmu_pair_counts()
            
            return cf.utils.convert_3d_counts_to_cf(N, N, rand_N, rand_N,
                                                    DD_counts, DR_counts, DR_counts, RR_counts) 

        def DDsmu2wp():
            
            DD_counts, DR_counts, RR_counts = create_DDsmu_pair_counts()

            return cf.utils.convert_rp_pi_counts_to_wp(N, N, rand_N, rand_N,
                                                          DD_counts, DR_counts, DR_counts, RR_counts, 
                                                          nbins, pimax)             
            
        def create_random_positions(N, rand_N):

            seed = 42
            np.random.seed(seed)

            #return random positon [X,Y,Z]
            return np.random.uniform(0, boxsize, rand_N), np.random.uniform(0, boxsize, rand_N), np.random.uniform(0, boxsize, rand_N)

        def writeIntoFile(results_CF,
                          name):
                        
            myheader=self.myconfig_array['catname'+str(self.a)]+' z='+str(float("{0:.2f}".format(redshift)))+' r_max: '+str(rmax)+' ['+pos_unit+'] '+'r_min: '+str(rmin)+' ['+pos_unit+'], boxsize: '+str(boxsize)+' ['+pos_unit+'], '
            if name.find('wp')!=-1: 
                myheader+='ngal: '+str(len(data))+', pi_max: '+str(pimax)+' ['+pos_unit+']\n'
            elif name.find('_s_'):
                myheader+='ngal: '+str(len(data))+', pi_max: '+str(pimax)+' ['+pos_unit+'], mu_max: '+str(mu_max)+'\n'
            else:
                myheader+='ngal: '+str(len(data))+'\n'

            if name.startswith('xi') or name.startswith('wp'):
                data_format="%0.8e\t"*5+'%0.8e'
                prefix_name_cf='r'
                endfix_myheader=' (5) npairs [-] (6) weight [-]'
                if name.find('_s_')!=-1:
                    prefix_name_cf='s'
                elif name.find('_of_')!=-1:
                    data_format="%0.8e\t"*4+'%0.8e'
                    endfix_myheader=' [npairs count] (5) weight [-]'
                
                myheader+='(1) '+prefix_name_cf+'min ['+pos_unit+']  (2) '+prefix_name_cf+'max ['+pos_unit+']  (3) '+prefix_name_cf+'_avg ['+pos_unit+'] (4) '+str(name)+endfix_myheader  
                
            elif name.find('2')!=-1:

                if name.find('smu2xi')!=-1 or name.find('rppi2xi')!=-1:
                    results_CF = mL.dim_expander(results_CF,3)                    
                    data_format="%0.8e\t%0.8e\t%0.8e"
                    
                    if name.find('smu')!=-1:
                        z_dim_name='mu'
                        z_dim = mu_max
                        zbins = nbins_mu
                        array=np.arange(0,z_dim,float(z_dim)/zbins)
                    else:
                        z_dim_name='pi'
                        z_dim = pimax
                        array=np.arange(1.0,z_dim+1.0,1.0) #if pi_max=40 40 radial bins are created, starting with 1.0
                        
                    myheader+='(1) r_mean ['+pos_unit+']  (2) '+z_dim_name+' '+'['+pos_unit+'] (3) '+name+' ['+pos_unit+']'
                        
                    
                    size=len(array)
                    
                    for i, rbin in enumerate(rbins[:-1]): 
                        for k, zbin in enumerate(array):
                            print('i:', i, 'k:', k, i*size+k, 'rbin:', rbin, 'zbin:', zbin)
                            results_CF[i*size+k,0] = (rbins[i]*rbins[i+1])**0.5
                            results_CF[i*size+k,1] = zbin                            
                                               
                else:
                    results_CF = mL.dim_expander(results_CF,2)
                
                    myheader+='(1) r_mean ['+pos_unit+']  (2) '+name+' ['+pos_unit+']'
                    data_format="%0.8e\t%0.8e"
                    i=0
                    while i<results_CF[:,0].size:
                        results_CF[i,0]=(rbins[i]*rbins[i+1])**0.5
                        i+=1
                            
                print(results_CF)

            myheader+=' galaxy types cents: '+str(centrals)+' no-sats: '+str(sats)+' orphans: '+str(orphans)

            if cut_list_operators_dict[operator]=='_mstar_' or cut_list_operators_dict[operator]=='_kroup_mstar_':
                cut_label=cut_list_operators_dict[operator]+str(cut)+'_'+str(cut_list_max_dict[cut])
            elif cut_list_operators_dict[operator].find('mstar')!=-1:
                cut_label=cut_list_operators_dict[operator]+str(float("{0:.2f}".format(np.log10(cut))))                  
            elif cut_list_operators_dict[operator].find('Mr')!=-1:
                cut_label=cut_list_operators_dict[operator]
            elif operator.find('Pop')!=-1:
                 cut_label='_'+operator           
            elif operator=='filaments' or operator=='knots':
                cut_label='_'+operator               
            elif cut!='':
                cut_label=cut_list_operators_dict[operator]+str(cut)
            else:
                cut_label=''

            #print(cut_label
            #print(data_format
            myOutput.writeIntoFile(#filename[0:len(filename)-4]+'_'+name+'.txt',
                                   filename[0:len(filename)-4]+'_'+gtype+cut_label+'_'+str(rmin)+'_'+str(rmax)+'_'+str(pimax)+'_'+name+'_test.txt',
                                   results_CF,
                                   myheader=myheader,
                                   data_format=data_format)
            
                
        print('\n\nTwoPCF Corrfunc()', 'self.i:', self.i, 'self.a:', self.a,'++++++++++++++++++++++++++++++++++++++++++++++\n')
        #print('catname', self.myconfig_array['catname'+str(self.a)], 'which:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_which'], 'z:', self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'] ) 

        sys.path.append(mycomp+'anaconda/Corrfunc')
        import Corrfunc as cf
        from os.path import dirname, abspath, join as pjoin
        
        redshift  = float(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'])

        cosmology=2


        if myconds_array['x_pos_unit'].find('h-1')!=-1:
            pos_unit='h-1Mpc'
        else:
            pos_unit='Mpc'            
      
        gtype_list=['centrals']#, 'all']

        gtype_dict=          {'centrals': 0,    'no': 2, 'no-sats': 1, 'all': 10}       
        gtype_dict_operators={'centrals': '==', 'no': '<', 'no-sats': '==', 'all': '<'}

        #self.myData2Process= mL.choose_random_sample(self.myData2Process, 95683)

        if mydata==False:
            mydata=self.myData2Process
        #print(np.info(mydata)
#        mydata = mL.choose_random_sample(self.myData2Process,int(self.myData2Process.size*0.1))
#        print('randomly chosen 10 %', mydata.size
        #print('here 3858', mydata   
       
        for gtype in gtype_list:
            
            print('gtype:', gtype, 'operator:',  gtype_dict_operators[gtype], 'condition:', gtype_dict[gtype])
            cut_list=['']
            cut_list_max_dict={}
            #> or <
            #cut_list=[1e10]
    
            #cut_list_operators=['PopA+k', 'PopA+f', 'PopB+k', 'PopB+f']
            #cut_list_operators_dict={'PopA+k':'==', 'PopA+f':'==', 'PopB+k':'==', 'PopB+f':'=='}

            #cut_list_operators=['filaments', 'knots']
            #cut_list_operators_dict={'filaments': '2', 'knots': '3'}
            
            cut_list_operators=['']
            cut_list_operators_dict={'': ''}
            
            for cut in cut_list:
                for operator in cut_list_operators:    
                    print('second loop:\ncut:', cut, 'operator:', operator, 'cut_list_operators_dict:', cut_list_operators_dict[operator])
                    try:
                        data      = self.myData2Process[np.where(self.myData2Process['Z']>0.05)]
                        #data['Z'] = cd.comoving_distance(data['Z'], **fidcosmo)*0.6777
                        #print('cosmology:', fidcosmo
                        #print(data['Z'][0:20]
                    except:
 
                                             
                        data = myData.selectData2Compute(mydata, 
                                                        selected_col='orphan', 
                                                        operator=gtype_dict_operators[gtype], 
                                                        condition=gtype_dict[gtype])
                                 
                        if operator=='filaments' or operator=='knots':
                            
                            data = myData.selectData2Compute(data, 
                                                            selected_col='env_1024', 
                                                            operator='==', 
                                                            condition=int(cut_list_operators_dict[operator]))                                
                                     
                        print('size after cut:', data.size)
                        centrals=len(np.where(data['orphan']==0)[0])
                        sats=len(np.where(data['orphan']==1)[0])
                        orphans=len(np.where(data['orphan']==2)[0])
                        
                        #print('galaxy types in sample: cents:', centrals, 'no-sats:', sats, 'orphans:', orphans

                    if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_ngal_random_sample']!='False':
                        data= mL.choose_random_sample(data, int(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_ngal_random_sample']))
                        print('after random selection of', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_ngal_random_sample'], 'galaxies:', data.shape)
    
                    pi_max_list=[float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_pimax'])]

                    for pimax in pi_max_list:
                        
                        nthreads  = int(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_nthreads'])
                        print('pimax:', pimax)

                        #large scales:    
                        calc_list_rmin=[float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_rmin'])]
                        calc_dict_rmax={float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_rmin']): float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_rmax'])}
                        
                        #small scales:
                        #calc_list_rmin=[0.5, 1, 5, 20]
                        #calc_dict_rmax={0.5: 3, 1: 10, 5: 75, 20:200}
                        
                        #final scales:
                        #calc_list_rmin=[0.5, 10]
                        #calc_dict_rmax={0.5: 10, 10: 200}                        
                        
                        for r_min in calc_list_rmin:                   
                            # Setup the bins
                            rmin      = r_min
                            rmax      = calc_dict_rmax[rmin]
                            nbins     = int(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_nbins'])

                            N=data['x_pos'].size
                            rand_N=3*N
                            nbins_mu=40
                            mu_max=1.0
                            
                            # Create the bins
                            rbins     = np.logspace(np.log10(rmin), np.log10(rmax), nbins+1)
                    
                            #if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_which']=='BOX':
                            if pos_unit.find('h')==-1:
                                boxsize   = self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size']/0.6777
                            else:
                                 boxsize   = self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size']             
                #            if self.myconfig_array['catname'+str(self.a)].startswith('SAGE'):
                #                boxsize   = 350.0/0.6777+0.1
                                
                            #print('chosen BOXSIZE:', boxsize, 'rmin/rmax:', rmin, '/', rmax
                            
                            def caseSwitcher(name):
                        
                                choose = {

                                        'xi'            : xi,       #3D real-space correlation function ξ(r)
                                        'wp'            : wp,       #2D projected correlation function wp(rp)

                                        'xi_of_r'       : xi_of_r,  #3D real-space pair count! Convert to ξ(r)          via convert_3d_counts_to_cf
                                        'xi_rp_pi'      : xi_rp_pi, #2D pair count! Convert to ξ(r)                     via convert_3d_counts_to_cf
                                        'xi_s_pi'       : xi_s_pi,  #2D pair count in redshift space! Convert to ξ(s,μ) via convert_3d_count_to_cf
                                        
                                        'DD2xi'         : DD2xi,    #3D real-space correlation function ξ(r) from xi_of_r pair counts   
                                        'DDrppi2xi'     : DDrppi2xi,#3D real-space correlation function ξ(rp,π) from xi_rp_pi pair counts  
                                        'DDrppi2wp'     : DDrppi2wp,#2D projected correlation function wp(rp) from xi_rp_pi pair counts
                                        'DDsmu2xi'      : DDsmu2xi, #2D redshift-space correlation function ξ(s,μ) with line of side separation less than s from xi_s_pi pair counts
                                        'DDsmu2wp'      : DDsmu2wp 

                                        }
                                    
                                func = choose.get(name)
                                return func()
                                       

                            print('calculate:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_calculate'], end='')
                            results_CF = caseSwitcher(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_calculate'])
                            #print(results_CF)
                            writeIntoFile(results_CF,  self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_calculate'])
                            print('----------\n\n')
                  

    def HODFunction(self,
                    mydata,
                    filename,
                    which_halomass=''):
             
        print('HODFunction()', 'self.i:', self.i, 'self.a:', self.a,'\n++++++++++++++++++++++++++++++++++++++++++++++\n')
        #print('catname', self.myconfig_array['catname'+str(self.a)], 'z:', self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'], end='')  
        #print('filename:', filename)

        def GalaxyTypeStats(data):
    
            data = mL.filter_data_before_analysis(data, 
                                                   5e11, 
                                                   3e15, 
                                                   'mhalo')
    
            def find_galaxy_type(data):
                central_ind   = (np.where(data['orphan'] == 0))[0] 
                data_centrals = data[central_ind]
                sat_ind       = (np.where(data['orphan'] == 1))[0] 
                data_sats     = data[sat_ind]        
                orphan_ind    = (np.where(data['orphan'] == 2))[0]
                data_orphans  = data[orphan_ind]
                
                return data_centrals, data_sats, data_orphans
    
    
            def find_unique_objects2(data):
                
                 unique_array, index, count = np.unique(data,
                                                return_index=True,
                                                return_counts=True)
                                                
                 return unique_array, index, count
                                             
            all_centrals, all_sats, all_orphans = find_galaxy_type(data)
            print('all galaxies: centrals', len(all_centrals), 'sats:', len(all_sats), 'orphans:', len(all_orphans), '\nall haloid:\n-------------------')                  
            unique_haloid, index_haloid, count_haloid = find_unique_objects2(data['haloid'])
            print('haloid', len(unique_haloid), 'ngal:', np.sum(count_haloid))
    
            centrals, sats, orphans = find_galaxy_type(data[index_haloid])                  
            print('centrals', len(centrals), 'sats:', len(sats), 'orphans:', len(orphans), '\n')
            
            myhalo_data=data[['haloid', 'hostid', 'orphan', 'mhalo']][index_haloid]
            myhalo_data= mL.dim_expander_struct(myhalo_data, 'mhalo', 'ngal') 
            myhalo_data['ngal']=count_haloid
    
            print('all hostid:\n-------------------')        
            unique_hostid, index_hostid, count_hostid = find_unique_objects2(data['hostid'])
            print('hostid', len(unique_hostid), 'ngal:', np.sum(count_hostid))
    
            centrals, sats, orphans = find_galaxy_type(data[index_hostid])                  
            print('centrals', len(centrals), 'sats:', len(sats), 'orphans:', len(orphans), '\n')


        def find_unique_objects(data):
                                            
             return np.unique(data, return_index=False, return_counts=True)


        def crossmacht_catalogs(data1, data2):
                      
            test = np.in1d(data1['haloid'], data2['haloid'])
                       
            #print('test haloid:', test_haloid
            #print(data[np.where(data['haloid']==test_haloid)]
            data_to_check=data[np.where(test==True)[:][0]]
            print(data_to_check)
            
            print('ngal haloid==hostid', len(data_centrals), 'sats:', len(data_sats), 'orphans:', len(data_orphans), 'gal test=True:', data_to_check.size )             
            

        def check_parents(data):
               
            data2write_not_found=[]
            data2write_found=[]
            
            id_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/data/prova.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.int64, skiprow=1) 
   
            header='\nHostHaloID\tMainHaloID\tPOS_X/POS_Y/POS_Z [h-1Mpc] GalaxyType HaloMass [h-1Msun] Mstar [h-1Msun]'
            for haloid, hostid in zip(id_array[:,0], id_array[:,1]):
            
                print('satellite HostHaloID:', hostid, '--> find parent\t', haloid, end='')
                mytest_data = data[np.where(data['hostid']==hostid)[:][0]]
                
                #print('check from complete list! satellite hostid: MainHaloID:', mytest_data['haloid'], 'HostHaloID:', mytest_data['hostid']    
                #print('find this satellite:', hostid, '--> and its parent', haloid, '-->', 
                
                parent=data[np.where(data['haloid']==haloid)[:][0]]
                parent=parent[np.where(parent['orphan']==0)[:][0]]
                print('parent', parent['haloid'], end='')
                            
                if parent['haloid']==haloid:                   
                    print('\t\t\tfound!')                
                    data2write_found.append(str(hostid)+'\t'+str(haloid)+'\t'+str(mytest_data['x_pos'])+'\t'+str(mytest_data['y_pos'])+'\t'+str(mytest_data['z_pos'])+'\t'+str(mytest_data['orphan'])+'\t'+str(mytest_data['mhalo'])+'\t'+str(mytest_data['mstar']))
                else:
                    print('\tNOT FOUND!')
                    try:
                        data2write_not_found.append(str(hostid)+'\t'+str(haloid)+'\t'+str(mytest_data['x_pos'])+'\t'+str(mytest_data['y_pos'])+'\t'+str(mytest_data['z_pos'])+'\t'+str(mytest_data['orphan'])+'\t'+str(mytest_data['mhalo'])+'\t'+str(mytest_data['mstar']))
                    except:
                        data2write_not_found.append(str(hostid)+'\t'+str(haloid)+'\t'+str(mytest_data['x_pos'][0])+'\t'+str(mytest_data['y_pos'][0])+'\t'+str(mytest_data['z_pos'][0])+'\t'+str(mytest_data['orphan'][0])+'\t'+str(mytest_data['mhalo'][0])+'\t'+str(mytest_data['mstar'][0]))
    
            
            myOutput.writeIntoFile(mycomp+'anaconda/pro/data/SAG_z_0.09_prova_not_found.txt',
                       data2write_not_found,
                       myheader='not found parents '+self.myconfig_array['catname'+str(self.a)]+' z='+str(float("{0:.2f}".format(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'])))+header,
                       data_format='%s',
                       data_is_string=False)
    
            myOutput.writeIntoFile(mycomp+'anaconda/pro/data/SAG_z_0.09_prova_found.txt',
                       data2write_found,
                       myheader='found parents'+self.myconfig_array['catname'+str(self.a)]+' z='+str(float("{0:.2f}".format(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'])))+header,
                       data_format='%s',
                       data_is_string=False)
            
            print('total galaxies:', id_array[:,0].size, 'parents found', len(data2write_found), 'NOT FOUND:', len(data2write_not_found))
                
    
        hod=True
        if which_halomass=='':
            halomass='mhalo_cents_200c'
        else:
            halomass=which_halomass
            
        #central_id='fofID' #LGALAXIES        
        central_id='haloid' #MD
        #self.myData2Process=self.myData2Process[np.where(self.myData2Process['orphan']<2)[:][0]]
#        mL.give_galaxy_type_info(self.myData2Process['orphan'])
#        exit()
        #print(self.myData2Process.shape, 'ngal sats:', self.myData2Process[np.where(self.myData2Process['orphan']>0)[:][0]].size


        #choose random number of satellites!        
        #rand_sats = mL.choose_random_sample(self.myData2Process[np.where(self.myData2Process['orphan']>0)[:][0]], int(self.myData2Process[np.where(self.myData2Process['orphan']>0)[:][0]].size-26309))
        #test      = np.in1d(self.myData2Process['hostid'], rand_sats['hostid'])             
        #self.myData2Process      = self.myData2Process[np.where(test==False)[0][:]]

        redshift=float(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'])
        print('redshift', redshift)
        which_halocat=''
        sim='MD'
        if sim=='MD':
            #path_to_halos=mycomp+'anaconda/pro/OBS/hod_MD_0.56.txt'
            #path_to_halos=mycomp+'anaconda/pro/OBS/HMF_MPDL2_FOF_mass_Msun.txt'
            path_to_halos=mycomp+'anaconda/pro/OBS/HMF_MDPL2_Rockstar_z_'+str(redshift)+'_M200c_Msun_Nhalos.txt'
            #path_to_halos=mycomp+'anaconda/pro/OBS/HMF_MDPL2_FOF_z_'+str(redshift)+'_mass_Msun_Nhalos.txt'
            data_halo = myData.readAnyFormat(config=False, mypath=path_to_halos, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=11) 
            which_halocat='from MDPL2 Rockstar cosmosim.org'
            bins=data_halo[:,0]
            Nhalos=data_halo[:,1]        
        elif sim=='SMD':
            #path_to_halos=mycomp+'anaconda/pro/OBS/HMF_SMPDL_FOF_mass_Msun.txt' 
            path_to_halos=mycomp+'anaconda/pro/OBS/HMF_SMDPL_Rockstar_z_'+str(redshift)+'_M200c_Msun_Nhalos.txt'
            #path_to_halos=mycomp+'anaconda/pro/OBS/HMF_SMDPL_FOF_z_'+str(redshift)+'_mass_Msun_Nhalos.txt'
            data_halo = myData.readAnyFormat(config=False, mypath=path_to_halos, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=11) 
            which_halocat='from SMDP2 M200c cosmosim.org'
            bins=data_halo[:,0]
            Nhalos=data_halo[:,1]
        elif sim=='MSI':
            data_halo = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/MSI_LGALAXIES_500Mpc_z_0.56_Nhalos.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)            
            bins=data_halo[:,0]
            
        elif sim=='TNG':
            path_to_halos=mycomp+'anaconda/pro/OBS/HMF_IllustrisTNG300_z_'+str(redshift)+'_tarsel_Nhalos.txt'
            data_halo = myData.readAnyFormat(config=False, mypath=path_to_halos, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)             
            which_halocat='from TNG300-1 cats from Celeste'
            bins=np.log10(data_halo[:,0])
            #print('bins:', bins
            Nhalos=data_halo[:,6]           
        else:
            binsize=0.103
            mhalo_min=12.1690
            mhalo_max=15.271
            
            binsize      = abs(mhalo_max-mhalo_min)
            bins         = np.arange(data_halo[0,0]-binsize/2.,data_halo[-1,0]-binsize/2., binsize)
  
        print('Halos from --> ', path_to_halos)
        
        if hod==True:
            pass
            #print('<N> = sum hosthalos + 1 / nhalos per mhalo bin')


#        print(len(np.where(self.myData2Process['orphan']==0)[:][0]))
#        exit()

        
        data = mL.filter_data_before_analysis(self.myData2Process, 
                                               'min', 
                                               'max', 
                                               halomass)

        histo        = np.zeros((len(bins),7), dtype=np.float64)
        
        counts, edges= np.histogram(np.log10(data[halomass]), bins)
        histo[:,0]   = 10**bins

        for gtype in ['all','centrals', 'sats']:
              
            i=1
            while i<len(histo[:,0])+1:
                #try:
                print('i:', i, edges[i-1], end='')

                if i==len(histo[:,0]):
                    #print('here! MAX', edges[i-1], np.log10(histo[i-1,0]),
                    data_bin   = data[np.where(np.log10(data[halomass])<=max(np.log10(data[halomass])))[:][0]]
                else:               
                    data_bin    = data[np.where(np.log10(data[halomass])<edges[i])[:][0]]
                    
                data_bin    = data_bin[np.where(np.log10(data_bin[halomass])>edges[i-1])[:][0]]

                if gtype=='all':
                    data_sel=data[np.where(data['orphan']>-1)[:][0]] 
                elif gtype=='centrals':
                    data_sel=data[np.where(data['orphan']==0)[:][0]]    
                    
                else:
                    data_sel=data[np.where(data['orphan']!=0)[:][0]] 

                
                if hod==True:                                      
                    #print('data_bin', data_bin.size 
                    test = np.in1d(data_sel[central_id],data_bin[central_id])
                    data_to_check=data_sel[np.where(test==True)[:][0]]
                    
                    #print('N'+gtype+':', data_to_check.size
                    unique_haloid, count_haloid = find_unique_objects(data_to_check[central_id]) 
                    #print(data_to_check[['haloid','hostid','orphan']]
                    #print(count_haloid[np.where(count_haloid!=1)[:][0]]
 
                    print('Nhaloid:', sum(count_haloid), 'Nunique:', len(unique_haloid), 'Nhalos:', float(Nhalos[i-1]))
                        
                    try:
                        histo[i-1,1]=(sum(count_haloid))/float(Nhalos[i-1])                        
                    except:
                        histo[i-1,1]='NaN'
                        pass
                    
                    histo[i-1,4]=len(unique_haloid)
                    histo[i-1,5]=sum(count_haloid)
                    histo[i-1,6]=float(Nhalos[i-1])
                    
                else:
                    
                    if gtype=='all':
                        data_sel=data_bin[np.where(data_bin['orphan']>-1)[:][0]] 
                    elif gtype=='centrals':
                        data_sel=data_bin[np.where(data_bin['orphan']==0)[:][0]]    
                        
                    else:
                        data_sel=data_bin[np.where(data_bin['orphan']!=0)[:][0]] 
                    
                    print('calculate scatter: ', end='')
                    mean_mstar=np.mean(np.log10(data_sel['mstar']))
                    #scatter
                    histo[i-1,1]=np.std(np.log10(data_sel['mstar']))
                    histo[i-1,4]=histo[i-1,1]/(data_sel.size**0.5)
                    histo[i-1,5]=histo[i-1,4]
                    print(histo[i-1,1], np.log10(histo[i-1,0]))
                    
                    histo[i-1,6]=data_sel.size               
    
                i+=1


            #print(histo[:,[0,1]]
            if hod==True:
                header='\n(1) '+halomass+' [Msun] (2) <ngal> [-] (3) -dx \t(4) dx\t(5) N unique\t(6) N '+gtype+' (7) Nhalos '+which_halocat
                myOutput.writeIntoFile(filename[0:len(filename)-4]+'_'+halomass+'_'+gtype+'.txt',
                                       histo,
                                       myheader='HODFunction '+self.myconfig_array['catname'+str(self.a)]+' z='+str(float("{0:.2f}".format(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'])))+' cumulative: no'+header,
                                       data_format="%0.8e",
                                       mydelimiter='\t')
            else:
                header='\n(1) '+halomass+' [Msun] (2) sigma(log10 Mstar [Msun]) [-] (3) -dx\t(4) dx\t(5) -dy\t(6) dy\t(7) Ngal/bin'
                myOutput.writeIntoFile(filename[0:len(filename)-4]+'_'+halomass+'_sc_'+gtype+'.txt',
                                       histo,
                                       myheader='HOD scatter '+self.myconfig_array['catname'+str(self.a)]+' z='+str(float("{0:.2f}".format(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'])))+' cumulative: no'+header,
                                       data_format="%0.8e",
                                       mydelimiter='\t')  
                     
    def calcFastHisto(self,
                      data,
                      filename,
                      col_name,
                      col_unit,
                      binning='log',
                      nbins=25,
                      binup2D=False,
                      custom_min_max=True):

        if binup2D==True:

            if col_name.find('mhalo')!=-1:
                data = data[np.where(data['orphan']==0)[:][0]]
                data = data[np.where(data[col_name]>=1e11)[:][0]]
                data = data[np.where(data[col_name]<1e15)[:][0]]

                histo, binsize = myFuncs.binUp(data[[col_name,'mstar']],
                                                    nbins,
                                                    histo_min=1e11,
                                                    histo_max=1e15,
                                                    log10bin=True,
                                                    binup_2D=True,
                                                    div_y2x=True,
                                                    use_MAD=True,
                                                    equal_bins=True) 
                #print(histo[:,[0,1]]
            
        else:
                
            #print('fast histogramm -->', col_name, '[', col_unit, ']',
            if custom_min_max==True:
                if col_name=='sfr':
                    data_min = 1e-4
                    data_max = 1000
    
                elif col_name=='zcold':
                    data_min = 8
                    data_max = 11.5
    
                elif col_name=='g-i':
                    data_min = 0
                    data_max = 3.5
    
                elif col_name=='ssfr':
                    data_min = 1e-13
                    data_max = 1e-8
        
                elif col_name=='mstar':
                    data_min = 1e9
                    data_max = 1e12                                          
                elif col_name=='mbar':
                    data_min = 1e8
                    data_max = 1e12    
                elif col_name=='mbh':
                    data_min = 1e5
                    data_max = 1e9
                elif col_name=='rbulgevsrdisk':
                    data_min = 0.001
                    data_max = 10000    
                elif col_name.find('mhalo')!=-1:
                    data_min = 1e11
                    data_max = 1e15
                elif col_name.find('age')!=-1:
                    data_min = 0.1
                    data_max = 10
                elif col_name.find('rhalf_mass')!=-1:
                    data_min = 1e-3
                    data_max = 10
                elif col_name.find('rbulge')!=-1 or col_name.find('rdisk')!=-1:
                    data_min = 1e-4
                    data_max = 1                     
                elif col_name=='vmax' or col_name=='vdisp':
                    data_min = 10
                    data_max = 10000
                elif col_name=='cgf':
                    data_min = 1e-5
                    data_max = 10
                elif col_name=='fbar':
                    data_min = 1e-7
                    data_max = 1                    
                elif col_name=='Tcons':
                    data_min = 1e-3
                    data_max = 100
                elif col_name=='BHeff':
                    data_min = 1e-10
                    data_max = 0.001
                elif col_name.startswith('j'):
                    data_min = -1
                    data_max = 5                      
            else:
                data_min=min(data)
                data_max=max(data)
                
            data = data[np.where(data>=data_min)[:][0]]
            data = data[np.where(data<=data_max)[:][0]]

            if binning=='log':           
                data = np.log10(data)
                data_min=np.log10(data_min)
                data_max=np.log10(data_max)
                  
            data = data[np.where(np.isfinite(data))[:][0]]        
    
            #print('data_min:', data_min, 'data_max:', data_max)
            binsize = (data_max-data_min)/nbins
            
            bins = np.linspace(data_min, data_max, nbins+1)
            #print(bins)
           
            counts, edges = np.histogram(data, bins)
            #print(edges)
    
            histo=np.zeros((nbins, 7), dtype=np.float32)
           
            if binning =='log':
                histo[:,0]= (10**edges[1:]*10**edges[:-1])**0.5
            else:
                histo[:,0]= (edges[1:]+edges[:-1])/2
    
            if self.volume=='False' or self.volume==False:
                self.volume=(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size']/0.6777)**3
                print('box size:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size']/0.6777, '--> set volume:', self.volume)
                
            histo[:,1]=counts/binsize/self.volume
            histo[:,2]=binsize
            histo[:,4]=histo[:,1]/counts**0.5
            histo[:,5]=histo[:,4]
            histo[:,6]= counts
      
        header='\n(1) '+col_name+' ['+col_unit+'] (2) Phi [Mpc-3 dex-1] (3) dx \t(4) -dx (5) -dy\t(6) +dy\t(7) N count'    

        myOutput.writeIntoFile(filename,
                               histo,
                               myheader='FAST HISTOGRAMM! '+self.myconfig_array['catname'+str(self.a)]+' z='+str(float("{0:.2f}".format(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'])))+' cumulative: NO'+header,
                               data_format="%0.8e",
                               mydelimiter='\t')    
                