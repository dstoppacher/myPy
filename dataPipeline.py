# Load packages 
#from __future__ import print_function
import numpy as np
import config as conf
import outputData as oD
import arangeData as aD
import myFuncs as mF
import myLib as mL

from cosmolopy import cd, fidcosmo
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
                  
    def showConfig(self):

        print ' '
        print '*****************************************************'
        print '** CONFIGURATIONS:'
        print '** -----------------------------'
        print '** CATNAME:   :', self.myconfig_array['catname'+str(self.a)]
        print '** CAT HUBBLE :', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_hubble_par'], 'CORRECT HUBBLE IN GENERAL :', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_correct_little_h']       
        print '** REDSHIFT   :', self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z']
        print '** BOX SIZE   :', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size']
        print '** VOLUME     :', self.volume
        print '**'        
        print '** CORRECT UNITS'
        print '** -------------'
        print '** HUBBLE     :', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_unit_corr_h']
        print '** Mpc        :', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_unit_corr_Mpc']        
        print '** Gyr        :', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_unit_corr_Gyr']
        print '** LUM        :', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_unit_corr_lum']
        print '** COMV       :', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_unit_corr_COMV']
        print '** kms-1      :', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_unit_corr_kms-1']
        print '**'
        print '** CONVERT UNITS'
        print '** -------------'        
        print '** AB MAGS    :', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_conv_to_AB_mag']
        print '** K-CORR APP :', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_apply_k_corr_app']
        print '** K-CORR ABS :', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_apply_k_corr_abs']
        print '** Z-BOOST    :', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_apply_z_boost']
        print '**'               
        print '** FILENAME   :', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]
        print '** DIRNAME    :', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_mydir']
        print '** FORMAT INFO:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format'], self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_format_info']
        print '*****************************************************'
        print ' '
	
    
    def readData(self,
                 myfilename=False,
                 myhalocat_code=False):
             
        if myfilename!=False:
            filename=myfilename
        else:
            filename=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)][len(self.myconfig_array['catname'+str(self.a)])+1::]
        print ' '
        print '########################################################################################################'
        print '#                                                                                                      '
        print '#     PROGRESS STATUS: readData() catname:', self.myconfig_array['catname'+str(self.a)], 'self.i:', self.i, 'self.a:', self.a
        print '#                                                                                                      '
        print '########################################################################################################'
        print ' '
            
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
                                             end_fileID=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_end_fileID'])

        self.myData2Process = data2process
        self.output_filename_code_space='' 
        self.tarsel_code_space=''
        
        #create spaces '_' in the filename:
        if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_output_filename_code']!='':
             self.output_filename_code_space='_'
        if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_tarsel_code']!='':
             self.tarsel_code_space='_'
       
        if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5' or self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format'].find('BINARY')!=-1 or self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='CATASCII':            
                self.correctAndConvertUnits(halocat_code=myhalocat_code)
        else:
            try:
                self.volume = myData.cat_attributes['volume']
            except:
                print 'volume attribute not exciting! ...'
                self.volume=False
        
                                                                                
    def correctAndConvertUnits(self,
                               halocat_code=False):
        
        try:
            print 'test format of [z] in mysnap_array ...',
            print float(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'])+1.0
            print 'is float! CORRECT!'
        except:
            print 'did not work! Convert [z] in mysnap_array to float!',
            self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'] =float(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'])
            print float(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'])+1.0, '--> worked!'

        
        
        print 'here: correctAndConvertUnits() --> set redshift: '
        if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5' or self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='BINARY_SAGE':
            self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'] = myData.redshift
        
            try:
                myData.redshift=self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z']
                myData.scale_factor=mL.redshift_to_expfactor(myData.redshift)
            except:
                try:
                    print '--> did not work, except!',
                    myData.redshift=float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_manual_input_redshift'])
                    self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z']=myData.redshift
                    myData.scale_factor=mL.redshift_to_expfactor(myData.redshift)
                    print myData.redshift+1.0, myData.scale_factor+1.0, 'redshift and scale factor set correctly!!!'
                except:
                    print '--> did not work, try to set scale factor!',
                    myData.scale_factor= float(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['a'])           
                    try:
                        print '---> did not work, second except!',
                        myData.scale_factor=mL.redshift_to_expfactor(myData.redshift)
                        print myData.redshift+1.0, myData.scale_factor+1.0, 'redshift and scale factor set correctly!!!'   
                    except:
                        print 'WARNING REDSHIFT WAS NOT CORRECTLY SET!!!'
       
        print 'myData.redshift:', myData.redshift, 'z in mysnap_array:', self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'], 'myData.scale_factor:', myData.scale_factor       
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

 
            
        self.myData2Process = mL.correct_units(self.myconfig_array['catname'+str(self.a)],
                                               self.myData2Process,
                                               self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+id_col_array_code],
                                               self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_correct_little_h'],
                                               myData.scale_factor)
                                                            
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
    def matchHaloCatWithResultCat(self,
                                  myconds_array):
                    
        print ' '
        print ' '                 
        print 'matchHaloCatWithResultCat()', 'self.i:', self.i, 'self.a:', self.a
        print '++++++++++++++++++++++++++++++++++++++++++++++'
        print ' '
 
        #Filter galaxy cat before intersecting with halo cat to exclude galaxy previously in order to save calc time!!!
        myheader=self.myconfig_array['catname'+str(self.a)]+' z='+str(float("{0:.2f}".format(myData.redshift)))+'\n'  
        mylog='Captains log-file: Stardate '+str(date_time_stamp)+'\n'+myheader+'++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n'+\
            'box volume: '+str(self.volume)+' [Mpc3]\n\nPreselection before matching galaxy cat with halo cat:\n\n'
        
        self.myData2Process, myheader, mydataformat, mylog, check_size, sel_col_list = mL.filter_data_before_write2file(self.myData2Process,
                                                                                                                            myconds_array,
                                                                                                                            self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array'],
                                                                                                                            myheader,
                                                                                                                            mylog)          
        print 'after pre-selection: galaxy cat:', self.myData2Process[:,0].size 
        #Load Halocatalouge corrisponding to the snapshotfiel!
        mydata = self.readData(myhalocat_code=True)
      
        intersect = list(set(self.myData2Process[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['haloid_col_id']]).intersection(set(mydata[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_halo_id_col_array']['halo_haloid_col_id']])))      
        
        halomask=[]
        galmask=[]

        matched_cat = np.zeros((self.myData2Process[:,0].size, self.myData2Process[0,:].size+mydata[0,:].size), dtype=np.float32)
        
        print 'halocat.shape:', mydata.shape 
        print 'galcat.shape:', self.myData2Process.shape        
        print 'matched_cat.shape:', matched_cat.shape, 'intersect:', len(intersect)
        print ' '
        
        centralmask = np.where(self.myData2Process[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['orphan_col_id']]==0)       
        satmask = np.where(self.myData2Process[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['orphan_col_id']]==1)            
        orphanmask = np.where(self.myData2Process[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['orphan_col_id']]==2)       
        print 'centralmask:', len(centralmask[0]), 'orphanmask:', len(orphanmask[0]), 'satmask:', len(satmask[0])       
          
        count = 0
        i=0
        while i<len(intersect):
            
            halomask = np.where(mydata[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_halo_id_col_array']['halo_haloid_col_id']]==intersect[i])        
            galmask = np.where(self.myData2Process[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['haloid_col_id']]==intersect[i])
            #rint ' '            
            #print 'i:', i, 'gal:', galmask, 'count:', count
            #print 'halo:', mydata[halomask]
            #print 'gal:', self.myData2Process[galmask]   
            print 'Iam still alive!', len(intersect)-i, 'to go ...'
            
            a=0
            while a<len(galmask[0]):
                #print 'i:', i, 'a:', a
                #print 'a:', a, 'halomask:', halomask[:][i], 'galmask:', galmask[:][i], galmask[:][i][a], 'count:', count, 'count+a:', count+a
                #print 'halo:', mydata[halomask[:][i]] 
                #print 'gal:', galmask[0], 'galmask[0][a]:', galmask[0][a], 'len(galmask[0]):', len(galmask[0])
                    
                matched_cat[count+a,0:self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_halo_id_col_array']['nr_entries']] = mydata[halomask]            
                matched_cat[count+a,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_halo_id_col_array']['nr_entries']::] = self.myData2Process[galmask[0][a]] 

#                if galmask[:][i][a]==self.myData2Process[:,0].size-1: 'Max:', self.myData2Process[:,0].size-1
#                if galmask[:][i][a]==0: 'Min', 0
                
                a+=1
                
            count+=a
            #print '-------------------'
                
            #print 'matched_cat:'
            #print matched_cat[count-a:count,:]
            #print ' '
            #print ' '

            i+=1

        checkmask = np.where(matched_cat[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_halo_id_col_array']['halo_haloid_col_id']]!=matched_cat[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['haloid_col_id']])       

        if len(checkmask[0])>0:
            print 'matchHaloCatWithResultCat()! the Haloids are not matching!!!'
            print ' '
            print 'Are you sure that you did everything alright!?!  --> Check Halo-catalogue haloids!!!!'
            exit()

        mydata = mL.correct_units(self.myconfig_array['catname'+str(self.a)],
                                  mydata,
                                  self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_halo_id_col_array'],
                                  1)
        
        test_z = myData.selectData2Compute(mydata, 
                                       selected_col=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_halo_id_col_array']['halo_redshift_col_id'],
                                       operator='!=',
                                       condition=myData.redshift)
                                       
        print 'test_z:', test_z.size, 'checkmask:', len(checkmask[0]) 
        
        i=0
        while i<myconds_array['nr_entries']:
            if  myconds_array['name'+str(i)].find('halo_')==-1:
                #print myconds_array['name'+str(i)], myconds_array[myconds_array['name'+str(i)]+'_col_id']
                myconds_array[myconds_array['name'+str(i)]+'_col_id']=int(myconds_array[myconds_array['name'+str(i)]+'_col_id'])+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_halo_id_col_array']['nr_entries']
            i+=1
        
        print 'matched_cat.shape:', np.resize(matched_cat,(count,matched_cat[0,:].size)).shape, 'len(galmask):', len(galmask), 'len(halomask):', len(halomask), 'count:', count

        selection_name = 'Selection after galaxy and halo catalogue match:'
                                                      
        return self.filterAndExtractData(myconds_array, data=matched_cat, myheader=myheader, mylog=mylog, myselection_name=selection_name)
        
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
                    #print 'check data=False!'
                    return self.myData2Process
            except:
                return data

        def createProgenitorList(data,
                                 header,
                                 log):

            print ' '
            print ' '                 
            print 'createProgenitorList()', 'self.i:', self.i, 'self.a:', self.a
            print '++++++++++++++++++++++++++++++++++++++++++++++'
            print ' '
    
            data_past = myData.readAnyFormat( data_shape=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_format_info'],
                                                 data_format=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format'],
                                                 mypath=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_mydir'],
                                                 #myfilename=filename,
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
                                                 snapid=str(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_snapid'+str(self.i+1)]),
                                                 value_is_redshift=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_use_snapidz_mapping'],
                                                 file_count=self.i+1,
                                                 halocat_code=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_halocat_code'],
                                                 halo_id_col_array=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_halo_id_col_array'],
                                                 start_fileID=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_start_fileID'],
                                                 nr_files_snapshot=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_nr_files_snapshot'],
                                                 end_fileID=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_end_fileID'])

                               
#            data_past =  myData.readAnyFormat(data_format=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format'],                     
#                                                id_col_array=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array'],
#                                                config_array=self.myconfig_array,
#                                                name_conv_map=self.name_conv_map,
#                                                catname=self.myconfig_array['catname'+str(self.a)],
#                                                snapid=str(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_snapid'+str(self.i)])

            print 'current redshift:', myData.redshift, 'z one step into past:', self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i+1)]['z'] 
            print 'data:', np.info(data)
            print 'data_past:', np.info(data_past)
          
            dset={}
            curr=0L            
            nhalo=0            
            for (haloid, nprog) in zip(data['haloid'], np.arange(0,data['npros'][nhalo],1)):               

                dset.update({'haloid'+str(nhalo): haloid})
                
                pro_idx = np.where(data_past['haloid']==data['firstProgenitorID'][nhalo])[0]

                print 'pro_idx:', pro_idx, 'data[pro_idx]:', data['haloid'][pro_idx]
                
                proList={}
                
                
                for n in nprog:

                    sib_idx = np.where(data_past['haloid']==data_past['siblingIndex'][n])[0]                  
                    sib     = data_past[sib_idx]
                    print 'n:', n, 'sib:', sib, 'sib_idx:', sib_idx
                    
                    proList.update({sib})
                    print 'proList:', proList

                dset.update({haloid+'_proList': proList})
                curr+=nprog
                nhalo+=1

            return dset

        def filterDensity(myheader, mylog):

            def get_ngal(ndens):
                print 'number density:', ndens, 'box size:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'], '[h-1Mpc] little-h:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_hubble_par'],
                if str(ndens).find('h3')!=-1:                    
                    ndens_unit=ndens[ndens.find('h3')::]
                    ndens=ndens[0:ndens.find('h3')]
                print 'ndens:', ndens, 'unit:', ndens_unit
                if ndens_unit.find('h3')!=-1:
                    print 'here unit with littel-h absorbed!\nngal=ndens*box^3'
                    ngal=int(float(ndens)*float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'])**3)
                else:
                    print 'here unit with littel-h factored out!\nngal=ndens*(box/h)^3'
                    ngal=int(float(ndens)*(float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'])/float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_hubble_par']))**3)
                print '----> ngal:', ngal
                return ngal           

            d=0
            while d<10:              
                try:
                    ngal = self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filter_density_ngal'+str(d)]
                    print 'ngal:', ngal
                    try:
                        print 'ngal: ', int(ngal),
                    except:
                        print '--> not found ...!'
                        ngal=get_ngal(ngal)
                    cut_name = self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filter_density_cut_name'+str(d)]
                except:
                    print 'END OF DENISTY FILTERING!'
                    exit()  

                
                c=0
                while c<100:
                    try:
                        name = self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filter_density'+str(c)]
                    except:
                        if c==0:
                            name = self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filter_density'] 
                        else:
                            print 'END OF NAME LOOP!'
                            print '-------------------'
                            print ' '
                            break
                  
                    mylog_filter = mylog    
                    print ' '
                    print 'filter catalog to a certain number density!'
                    print '-----------------------------------------------'
                    print '--> which property? ', name, 'ngal:', ngal, 'col id of the sorted column:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array'][name+'_col_id']
                    print 'd:', d, 'c:', c, 'cut_name:', cut_name
    
                    data_sorted = data[np.argsort(data[name])]
                    
                    mylog_filter+='\nselected for number density property '+name+': '+str(int(ngal)/float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'])**3)+' h3Mpc-3\n'
                    print 'data after sort:', data.shape, 'mynumber_density_ngal', int(ngal),
                                    
                    size=len(data_sorted[name])
                    if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filter_density_select_highest']=='True':
                        start_size=size-int(ngal)

                        if start_size<0: start_size=0
                        print '--> select sorted data set: data[', start_size, ':', size, ']'                          
                        data_sorted=data_sorted[start_size:size]
                        
                    else:
                        print '--> select smallest ...'
                        data_sorted=data_sorted[0:int(ngal)]

  
                        
                    print 'check size of the filtered data_array:', len(data_sorted[name]), data_sorted.shape
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

            print 'END OF FILTERING DENSITY!'
            print '-----------------------------'
            exit()

        def selectRegion(data, myheader, mylog, mylog_region):

            print ' '
            print 'filterData() --> selectRegion() ...'
            print '-------------------'
            print ' '                            
            

            #read file with information about the cutted region
            cut_info_array = myData.readAnyFormat(config=False, mypath=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_path_to_data']+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_filename'], data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=0)
            #print cut_info_array
            
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
                print 'set periodic boundary conditions ...'

                print 'box_size', box_size, 'radius:', radius
                pos=['x_pos','y_pos','z_pos']
               
                i=0
                for k in pos:
                    print 'k:', k, 'i:', i,
                    data_new = data[np.where(data[k]<radius)]
                    data_new[k]+=float(box_size)
                    print 'len1:', len(data_new),
                    if i==0:
                        data_add=data_new
                    else:
                        data_add = np.concatenate((data_add, data_new), axis=0)                    
                    data_new = data[np.where(data[k]>box_size-radius)]
                    data_new[k]-=float(box_size)
                    print 'len2:', len(data_new),
                    data_add = np.concatenate((data_add, data_new), axis=0)
                    i+=1
                    print 'data_add:', len(data_add)
            
                data = np.concatenate((data, data_add), axis=0)
                
                print 'data after periodic boundary conditions!', data.shape
            
            mylog_region+='name of region'.ljust(16)+('X ['+myunit+']').ljust(12)+('Y ['+myunit+']').ljust(12)+('Z ['+myunit+']').ljust(12)+(str(size)+' ['+myunit+']').ljust(18)+'total ngal'.ljust(12)+'centrals'.ljust(12)+'satellites'.ljust(12)+'orphans'.ljust(12)+'mhalo'.ljust(10)+'X/Y/Z'.ljust(24)+'MainHaloID'.ljust(20)+'offset: X/Y/Z'.ljust(24)+'\n'
            i=0
            while i<cut_info_array[:,0].size:
                #cut_info_array[i,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_col_id_radius']]=5.0                                             
                mylog_filter=mylog
                print ' '
                print '--------------'
                print 'check: original data size:', data.shape
                print 'i:', i, 'region name',

                if i<4:
                    region_name='VoidMDPL_'+str(int(cut_info_array[i,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_region_name']])).zfill(4)
                else:
                    region_name='region'+str(int(cut_info_array[i,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_region_name']])).zfill(4)

                print region_name

                if str(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_units'])!=myunit:
                    print 'adjusted units to [h-1Mpc] --> ', cut_info_array[i,:],
                    
                    cut_info_array[i,[cut_id_array['0_id_selReg'], cut_id_array['1_id_selReg'], cut_id_array['2_id_selReg'], self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_col_id_radius']]]*=0.6777                    

                mylog_filter+='\nselect a specific subvolume \n---------------------------------------------- \nregion number: '+region_name+', radius ['+myunit+']: '+str(cut_info_array[i,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_col_id_radius']]).zfill(4)+', coordinates: X ['+myunit+']: '+str(cut_info_array[i,cut_id_array['0_id_selReg']])+', Y ['+myunit+']: '+str(cut_info_array[i,cut_id_array['1_id_selReg']])+', Z ['+myunit+']: '+str(cut_info_array[i,cut_id_array['2_id_selReg']])+'\n\n'
 
                if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_volume_type']=='SPHERE':
                    print 'Volume: SPHERE! selection coordinate:', cut_info_array[i,cut_id_array['0_id_selReg']],cut_info_array[i,cut_id_array['1_id_selReg']],cut_info_array[i,cut_id_array['2_id_selReg']],'radius: ', cut_info_array[i,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_col_id_radius']]
                   
                    def find_cluster_members(data_pos):
                        a=0
                        while a<3:
                            print 'a:', a, cut_id_array[str(a)], cut_info_array[i,cut_id_array[str(a)+'_id_selReg']], 
                            data_pos[cut_id_array[str(a)]]+=(500.0-cut_info_array[i,cut_id_array[str(a)+'_id_selReg']])
    
                            print 'ngal>boxsize:', len(data_pos[cut_id_array[str(a)]][np.where(data_pos[cut_id_array[str(a)]]>box_size)]),
                            print 'ngal<0.0    :', len(data_pos[cut_id_array[str(a)]][np.where(data_pos[cut_id_array[str(a)]]<0.0)])
    
                            data_pos[cut_id_array[str(a)]][np.where(data_pos[cut_id_array[str(a)]]>box_size)]-=box_size
                            data_pos[cut_id_array[str(a)]][np.where(data_pos[cut_id_array[str(a)]]<0.0)]+=box_size                       
                            a+=1
                            
                        return np.where(((data_pos['x_pos']-500.0)**2+(data_pos['y_pos']-500.0)**2+(data_pos['z_pos']-500.0)**2)**0.5 < cut_info_array[i,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion_col_id_radius']])
                        
                    mask = find_cluster_members(data[['x_pos', 'y_pos', 'z_pos']])
                    
                    data_sorted=data[mask[:][0]]

                print 'after selection -->', data_sorted.shape, len(data_sorted)
                try:
                    print min(data_sorted['x_pos']), '/', max(data_sorted['x_pos']), ',', min(data_sorted['y_pos']), '/', max(data_sorted['y_pos']), ',', min(data_sorted['z_pos']), '/', max(data_sorted['z_pos']),  
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
                
                print 'centrals:', num_centrals, 'sats:', num_sats, 'orphans:', num_orphans, 'total:', num_centrals+num_sats+num_orphans
                print 'data_sorted[[mhalo,x,y,z]]', data_sorted[['mhalo','x_pos', 'y_pos', 'z_pos']]
                max_mhalo= max(data_sorted['mhalo'])
                
                offset_x = cut_info_array[i,cut_id_array['0_id_selReg']] - data_sorted['x_pos'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0]
                offset_y = cut_info_array[i,cut_id_array['1_id_selReg']]- data_sorted['y_pos'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0] 
                offset_z = cut_info_array[i,cut_id_array['2_id_selReg']] - data_sorted['z_pos'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0] 
             
                print 'max mhalo:', max_mhalo, 'where(max_mhalo)[0]:', np.where(data_sorted['mhalo']==max_mhalo)[0],
                print 'max(data_sorted[[mhalo,x,y,z])', data_sorted['mhalo'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0],\
                                                        data_sorted['x_pos'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0], \
                                                        data_sorted['y_pos'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0],\
                                                        data_sorted['z_pos'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0], \
                                                        data_sorted['haloid'][np.where(data_sorted['mhalo']==max_mhalo)[0]][0]
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
     
                #print 'outputfilename:', myfilename
              
                
                data_sorted, header, mydataformat, mylog_filter, check_size, sel_col_list = mL.filter_data_before_write2file(data_sorted,
                                                                                                                        myconds_array,
                                                                                                                        myconds_array,
                                                                                                                        myheader,
                                                                                                                        mylog_filter)                
                                                                                                                        
                mylog_filter+='\nngalaxies in the file: '+str(data_sorted.size)+'\n\n'+'link and filename of catalogue: '+myfilename

                #print sel_col_list
                #print mydataformat
            
                writeHDF5(data_sorted, myfilename[0:len(myfilename)-3]+'hdf5')
                writeASCII(data_sorted[sel_col_list], myfilename, header, mydataformat)
                     
                writeLOGFILE(mylog_filter, myfilename[0:len(myfilename)-4]+'_log.txt')

                i+=1               
                #except:
                #    print 'no galaxies found in', region_name
            writeLOGFILE(mylog_region, myfilename[0: myfilename.find('reg')]+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_output_filename_code']+output_space+'subsamples_log.txt')
            print 'END OF SELECTREGION!'
            print '-----------------------'
            exit()
            
        def writeCoordinates(data,
                             filename,
                             ngal):   
           

            print 'Write coordinates in file!'
                 
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
            
            print 'norphans:', len(data_orphans), 'ncentrals', len(centrals) 

            orphan_haloid_unique, orphan_counts = np.unique(data['haloid'][orphan_ind],
                                                            return_counts=True)

            central_haloid_unique, idx_first_haloid_in_central = np.unique(centrals['haloid'],
                                                                              return_index=True)
      
            sort_mask = np.argsort(data['haloid'][orphan_ind])

            orig_idx_haloid_to_orig_array_ind = orphan_ind[sort_mask]

            centrals=centrals[idx_first_haloid_in_central]
                
            curr=0L
            for (host_idx, norphans) in zip(np.arange(len(centrals)), orphan_counts):
                dest_sel = np.s_[curr:curr+norphans]

                orphan_indices_this_host = orig_idx_haloid_to_orig_array_ind[dest_sel]

                #print 'host_idx:', host_idx, 'haloid centrals:', centrals['haloid'][host_idx], 'haloid_orphan:', data['haloid'][orphan_indices_this_host]
                #print 'dest_sel:', dest_sel, 'counts:', norphans, 'orphan_ind this host:', orphan_indices_this_host
                #print 'before:', data['x_pos'][orphan_indices_this_host], centrals['x_pos'][host_idx],  
                
                #data['x_pos'][orphan_indices_this_host]=centrals['x_pos'][host_idx]
                for f in ['x_pos', 'y_pos', 'z_pos']:
                    data[f][orphan_indices_this_host] = centrals[f][host_idx]

                #print '--> after:', data['x_pos'][orphan_indices_this_host] 
                #print '-------------'
                #print ' '
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
            
            print 'norphans:', len(data_orphans), 'ncentrals', len(centrals) 

            orphan_haloid_unique, orphan_counts = np.unique(data[id_name][orphan_ind],
                                                            return_counts=True)

            central_haloid_unique, idx_first_haloid_in_central = np.unique(centrals[id_name],
                                                                              return_index=True)
      
            sort_mask = np.argsort(data[id_name][orphan_ind])

            orig_idx_haloid_to_orig_array_ind = orphan_ind[sort_mask]

            centrals=centrals[idx_first_haloid_in_central]
                
            curr=0L
            for (host_idx, norphans) in zip(np.arange(len(centrals)), orphan_counts):
                dest_sel = np.s_[curr:curr+norphans]

                orphan_indices_this_host = orig_idx_haloid_to_orig_array_ind[dest_sel]

                #print 'host_idx:', host_idx, 'haloid centrals:', centrals['haloid'][host_idx], 'haloid_orphan:', data['haloid'][orphan_indices_this_host]
                #print 'dest_sel:', dest_sel, 'counts:', norphans, 'orphan_ind this host:', orphan_indices_this_host
                #print 'before:', data['x_pos'][orphan_indices_this_host], centrals['x_pos'][host_idx],  
                
                #data['x_pos'][orphan_indices_this_host]=centrals['x_pos'][host_idx]
                for f in ['mhalo_cents']:
                    data[f][orphan_indices_this_host] = centrals[f][host_idx]

#                print '--> after:', data['x_pos'][orphan_indices_this_host] 
#                print '-------------'
#                print ' '
                curr += norphans

            return data

        def extractHalomassSergio(data):
            """
            #Extract halomass for Sergio
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            - Halo Mass (of the group, not the individual galaxy, eg. the satellite halo mass = central halo mass)
            - Halo ID (again, from the group)
            - Stellar Mass
            - SFR
            - Type
            - Position
            
            In principle, I will need it for the full simulation, but for the first part of the project, for halo masses above 10^14 h^-1 Msun
            should be enough. Also, for now only at z=0.
            
            Galacticus: all subhalo masses are included in the mainhalomass
            
            
            """                
            data_centrals = data[np.where(data['orphan']==0)[0]]
            print 'nr centrals:', data_centrals.shape, 
            
            data_centrals = myData.selectData2Compute(data_centrals, 
                                             selected_col='mhalo', 
                                             operator='>', 
                                             condition=1e14)
             
            print 'after mhalo > 1e14', data.shape,

#            u=0
#            while u<10:
#            haloid=data_centrals['hostid'][u]

            test = np.in1d(data['haloid'],data_centrals['haloid'])
            print 'test.shape:', test.shape,
    
            data_selected=data[np.where(test==True)[0]]

            print 'output:', data.shape            
                                            
            all_centrals, all_sats, all_orphans = mL.return_lists_by_galaxy_type(data_selected)   

            print 'all galaxies: centrals', len(all_centrals), 'sats:', len(all_sats), 'orphans:', len(all_orphans), '\n' 
    
#                print data_selected
#                u+=1
#            exit()
            return data_selected

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



        print ' '
        print '########################################################################################################'
        print '#                                                                                                      '
        print '#     PROGRESS STATUS: filterAndExtractData() catname:', self.myconfig_array['catname'+str(self.a)], 'self.i:', self.i, 'self.a:', self.a
        print '#                                                                                                      '
        print '########################################################################################################'
        print ' ' 

        data = check_data(data)

        if redshift==False:
            redshift=myData.redshift
        if scale_factor==False:
            scale_factor=myData.scale_factor

        i=0
        count_Mzstar=0
        count_Mzgas=0
        count_mstar=0
        count_mstar_IC=0
        count_mhalo=0
        count_sfr=0
        count_mag=0

        if preprocessing_only==True:  
            if self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z']<1.0:
                try:
                    data = data[np.where((data['mstar_disk']+data['mstar_spheroid'])>1e8)[:][0]]                    
                except:
                    data = data[np.where(data['mstar']>1e8)[:][0]]
                print 'z<1: low mass galaxies with mstar<1e8 [Msun] are excluded:', data.shape
    
            elif self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z']<2.0:            
                try:
                    data = data[np.where((data['mstar_disk']+data['mstar_spheroid'])>1e7)[:][0]]                    
                except:
                    data = data[np.where(data['mstar']>1e7)[:][0]]
                print 'z<2: low mass galaxies with mstar<1e7 [Msun] are excluded:', data.shape 
                
            elif self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z']<3.0:            
                try:
                    data = data[np.where((data['mstar_disk']+data['mstar_spheroid'])>1e6)[:][0]]                    
                except:
                    data = data[np.where(data['mstar']>1e6)[:][0]]
                print 'z<3: low mass galaxies with mstar<1e6 [Msun] are excluded:', data.shape   

           
        while i<self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['nr_entries']:
                     
            if self.myconfig_array['catname'+str(self.a)].startswith('Galacticus') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)]=='hostid' and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5':
                #print 'GALACTICUS correct HOSTID of satellites and orphans!'                
                mask = np.where(data['orphan']!=0)               
                data['hostid'][mask[:][0]]=data['satelliteNodeIndex'][mask[:][0]]

            if self.myconfig_array['catname'+str(self.a)].startswith('Galacticus') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)]=='haloid' and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5':
                #print 'GALACTICUS: correct HALOID of satellites and orphans!'                    
                mask = np.where(data['orphan']!=0)
                data['haloid'][mask[:][0]]=data['parentIndex'][mask[:][0]]

            if self.myconfig_array['catname'+str(self.a)].startswith('Galacticus') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)]=='mhalo' and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5' and count_mhalo==0:
                #print 'GALACTICUS: set mhalo of the satellites to sattelliteBoundMass aka mhalo_sat!!!!'                
                mask = np.where(data['orphan']!=0)
                data['mhalo'][mask[:][0]]=data['mhalo_sat'][mask[:][0]]
                mask = np.where(data['orphan']!=0)              

            if self.myconfig_array['catname'+str(self.a)].startswith('Galacticus') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)]=='Mzgas' and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5' and count_Mzgas==0:
                #print 'GALACTICUS: MZgas = MZgas_spheroid + MZgas_disk!!!!'
                count_Mzgas+=1
                data['Mzgas']=data['Mzgas_spheroid']+data['Mzgas_disk']

            if self.myconfig_array['catname'+str(self.a)].startswith('Galacticus') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)]=='Mzstar' and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5' and count_Mzstar==0:
                #print 'GALACTICUS: MZstar = MZstar_spheroid + MZstar_disk!!!!'
                count_Mzstar+=1
                data['Mzstar']=data['Mzstar_spheroid']+data['Mzstar_disk']



            if self.myconfig_array['catname'+str(self.a)].startswith('SAG_') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)]=='sfr_disk' and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5' and count_sfr==0:
                print 'SAG: set sfr_disk!!!!'
                count_sfr+=1
                data['sfr_disk']-=data['sfr_spheroid']
            
            if self.myconfig_array['catname'+str(self.a)].startswith('SAGE') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)]=='mstar_disk' and (self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='BINARY_SAGE' or self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5')  and count_mstar==0:
                print 'SAGE: set mstar_disk!!!!'
                count_mstar+=1
                data['mstar_disk']-=data['mstar_spheroid']
                
            if self.myconfig_array['catname'+str(self.a)].startswith('SAGE') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)]=='mstar_disk' and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='BINARY_SAGE' and count_mstar_IC==0:
                print 'SAGE: set total mstar! Mstar+Mstar_IC!!!!'
                count_mstar_IC+=1
                try:
                    data['mstar+IC']+=data['mstar_IC']
                except:
                    pass
            
            if self.myconfig_array['catname'+str(self.a)].startswith('SAGE') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)]=='Mzstar_disk' and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='BINARY_SAGE' and count_Mzstar==0:
                print 'SAGE: set metals mstar_disk!!!!' 
                data['Mzstar_disk']-=data['Mzstar_spheroid']
                count_Mzstar+=1

            if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)].find('AB')!=-1 and count_mag==0:
                data = mL.convert_units_Galacticus(
                                      self.myconfig_array['catname'+str(self.a)],
                                      data,
                                      self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array'],
                                      redshift,
                                      telescope_name=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_telescope_name'],
                                      apply_k_corr_app=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_apply_k_corr_app'],
                                      apply_k_corr_abs=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_apply_k_corr_abs'],                                              
                                      apply_z_boost=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_apply_z_boost'],
                                      hubble_par=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_hubble_par'],
                                      unit_code=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_unit_code'],
                                      use_kcorrect=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_which_kcorrect'],
                                      cosmology=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_cosmology'])
                count_mag+=1

            if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)].find('cut')!=-1:       
                data = mL.calc_colour_cut_parameter(data, 
                                                    self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array'],
                                                    self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name'+str(i)])
            i+=1

            
        if myheader==False:
            myheader=self.myconfig_array['catname'+str(self.a)]+' z='+str(float("{0:.2f}".format(myData.redshift)))+'\n'

        #print self.myconfig_array
            
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
                    print 'Correct Galacitucs orphans positions!!!! --> DONE!'
                except:
                    pass
    
            if self.myconfig_array['catname'+str(self.a)].startswith('Galacticus') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5':
                try:
                    data=faster_central_mhalo_for_sats(data)
                    print 'Set central mhalo (HOD calculation)!!!! --> DONE!'
                except:
                    pass               
              
        try:
            data['mstar'] = data['mstar_disk']+data['mstar_spheroid']
            #print 'mstar=mstar_disk+mstar_spheroid'
        except:
             pass
        try:
            data['sfr'] = data['sfr_disk']+data['sfr_spheroid']
            #print 'sfr=sfr_disk+sfr_spheroid'            
        except:
            pass            
        try:
            data['mcold'] = data['mcold_disk']+data['mcold_spheroid']
            #print 'mcold=mcold_disk+mcold_spheroid'            
        except:
            pass
        try:
            data['Mzstar'] = data['Mzstar_disk']+data['Mzstar_spheroid']
            #print 'Mzstar=Mzstar_disk+Mzstar_spheroid'            
        except:
            pass        
        try:
            data['Mzgas'] = data['Mzgas_disk']+data['Mzgas_spheroid']
            #print 'Mzgas=Mzgas_disk+Mzgas_spheroid'            
        except:
            pass
        try:
            data['ssfr'] = data['sfr']/data['mstar']
            print 'ssfr=sfr/mstar'
            print 'ssfr min/max:', min(data['ssfr']), '/', max(data['ssfr'])            
        except:
            pass
        try:
            data['zcold']=8.69+np.log10(data['Mzgas']/(data['mcold']*0.0134))
            print 'zcold=8.69+log10(Mzgas/Mcold/0.0134)'
        except:
            pass
        try:
            data['cgf']=data['mcold']/data['mstar']
            print 'cgf=mcold/mstar'
        except:
            pass

        #print 'vmax: min/max:', format(min(data['vmax']), '0.2f'), '/', format(max(data['vmax']), '0.2f')           
        #print 'vdisp: min/max:', format(min(data['vdisp']), '0.2f'), '/', format(max(data['vdisp']), '0.2f')         

#        try:
#            data['mAB_dA_total_cut_r_i'] = data['MAB_dA_total_r']-data['MAB_dA_total_i']
#            data['mAB_dA_total_cut_g_r'] = data['MAB_dA_total_g']-data['MAB_dA_total_r'] 
#            data['mAB_dA_total_cut_g_i'] = data['MAB_dA_total_g']-data['MAB_dA_total_i']          
#        except:
#            pass
#
        
        try:
            from cosmolopy import cparam, cd,cc
            fidcosmo = cparam.Planck(flat=True, extras=False)    
            data['mean_age_stars_disk'] = cd.age(redshift, **fidcosmo)/cc.Gyr_s-(data['age_sfr_int_disk']/data['sfr_int_disk']/1e9)
            
            print 'z=', redshift, 'age:', format(cd.age(redshift, **fidcosmo)/cc.Gyr_s, '0.2f'), 'mean_age_stars_disk=age_Universe(z)-age_sfr_int_disk/sfr_int_disk [Gyr]',
            print 'min/max:', format(min(data['mean_age_stars_disk']), '0.2f'), '/', format(max(data['mean_age_stars_disk']), '0.2f')           
        except:
            pass

        try:
            from cosmolopy import cparam, cd,cc
            fidcosmo = cparam.Planck(flat=True, extras=False)    
            data['mean_age_stars_spheroid'] = cd.age(redshift, **fidcosmo)/cc.Gyr_s-(data['age_sfr_int_spheroid']/data['sfr_int_spheroid']/1e9)
            
            print 'z=', redshift, 'age:', format(cd.age(redshift, **fidcosmo)/cc.Gyr_s, '0.2f'), 'mean_age_stars_spheroid=age_Universe(z)-age_sfr_int_spheroid/sfr_int_spheroid [Gyr]',
            print 'min/max:', format(min(data['mean_age_stars_spheroid']), '0.2f'), '/', format(max(data['mean_age_stars_spheroid']), '0.2f')           
        except:
            pass
        
#        
        try:
            if np.all(data['pop']==-99):
                print 'divide sample!'
                y=np.log10(data['ssfr'])
                x=np.log10(data['sfr'])
                prop=(y+11.16)/1.12-x
                print prop
                data['pop'][np.where(prop<=0.0)]=2
                data['pop'][np.where(prop>0.0)]=1
    
                print 'min/max', min(data['pop']), max(data['pop'])               
        except:
            pass   
        
        
#
#        try:
#            if np.all(data['env_512']==-99.0) or np.all(data['env_1024']==-99.0):
#                print '\n+++++++++++++++++++\nset environment',
#                data2cross = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/data/Galacticus_1Gpc/Galac_with_Environments_new.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=2) 
#                                
#                try:
#                    data['env_512']=data2cross[:,3]
#                    print 'env_512 set!'
#                except:
#                    print 'env_512 not excisting ...'                        
#                try:
#                    data['env_1024']=data2cross[:,4]
#                    print 'env_1024 set!'
#                except:
#                    print 'env_1024 not excisting ...'
#                    data2cross = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/data/Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3_haloid.txt', data_format='ASCII', data_shape='shaped', delim='  ', mydtype=np.uint64, skiprow=2) 
#                    print data2cross[0:10,0]
#        
#                    test = np.in1d(data2cross[:,1],data['hostid'])
#                    data2cross=data2cross[test]
#                    print data2cross[:,0].size
#                    for i, hostid in enumerate(data2cross[:,1]):
#                        #print 'hostid', hostid, 'i:', i, 'check hostid:', data2cross[i,1], 'env_512:',  data2cross[i,3], 'env_1024:', data2cross[i,4]
#                        data['env_512'][np.where(data['hostid']==hostid)[:][0]]=data2cross[i,3]
#                    data['env_1024'][np.where(data['hostid']==hostid)[:][0]]=data2cross[i,4]
#                    print '--> DONE!\n'
#                    
#        except:
#            pass
#
#            
#        try:
#
#            if np.all(data['mbasic_200c']==-99.0):
#                print '\n+++++++++++++++++++\nconvert mbasic from viral to 200c',
#                data['mbasic_200c'], data['NFW_con']=mL.convert_halo_mass(data['mbasic'],
#                                                                         data['NFW_con'],
#                                                                         data['orphan'])
#                print '--> DONE!\n'
#
#        except:
#            pass
#
#
#        try:
#            if np.all(data['mbasic_cents_200c']==-99.0):
#                print '\n+++++++++++++++++++\ncreate coluumn for central mbasic masses in satellites',
#                data['mbasic_cents_200c'], data['NFW_con']=mL.convert_halo_mass(data['mbasic_cents'],
#                                                                         data['NFW_con'],
#                                                                         data['orphan'])
#                print '--> DONE!\n'
#        except:
#            pass
##
 
#        data_cents=data[np.where(data['orphan']==0)[0][:]]     
#        data_sats=data[np.where(data['orphan']==1)[0][:]]
#        data_os=data[np.where(data['orphan']==2)[0][:]]
#        print 'all:', data.size, 'centrals:', data_cents.size, 'sats:', data_sats.size+data_os.size, 'orphans:', data_os.size
#        print 'fractions cent:', data_cents.size/float(data.size), 'sats:', (data_sats.size+data_os.size)/float(data.size), 'orphans in sats:', data_os.size/float(data_sats.size)
#        
#        
#        data = mL.choose_random_sample(data,
#                         int(data.size*0.10))     
#        print np.info(data)
      

        print 'CHECK HERE: redshift=', redshift
        #redshift=0.0
        try:
            if np.all(data['mhalo_200c']==-99.0):
                data['NFW_con'] = mL.calculate_NFW_con_with_fit(data['mhalo'], 
                                                       float(redshift),
                                                       cosmology='Planck', 
                                                       overdens='vir')
            
                data['mhalo_200c']=mL.convert_halo_mass(data['mhalo'], 
                                                                    data['NFW_con'],
                                                                    data['orphan'])
            
                print '--> DONE!\n'
        except:
            pass
        
        try:
            if np.all(data['mhalo_cents_200c']==-99.0):
                NFW_con = mL.calculate_NFW_con_with_fit(data['mhalo_cents'], 
                                                       float(redshift),
                                                       cosmology='Planck', 
                                                       overdens='vir')
                
                data['mhalo_cents_200c']=mL.convert_halo_mass(data['mhalo_cents'], 
                                                                    NFW_con,
                                                                    data['orphan'])
          
                print '--> DONE!\n'
        except:
            pass 

#
#        try:
#            if np.all(data['weight_tot']==-99.0):
#                print '\n+++++++++++++++++++\napply weights to CMASS DR12',
#                data['weight_tot']=data['weight_systot']*(data['weight_noz']+data['weight_cp']-1)
#                totgal = sum(data['weight_noz']+data['weight_cp']-1)
#                weight_max = max(np.cumsum(data['weight_tot']))
#                
#                print 'ngal weight normal:', sum(data['weight_tot']), 'ngal totgal:', totgal, 'ngal weight_max:', weight_max, '---> ngal special weight:',           
#                print sum(data['weight_tot']*totgal/weight_max)
#                print '--> DONE!\n'                
#
#        except:
        try:
            if np.all(data['weight_tot']==-99.0):
                print 'set weights to 1!',
                data['weight_tot']=np.ones((data.size,), dtype=np.int8)
                print '--> DONE!\n'
        except:
            pass
#
          
#        import pandas as pd
#        
#        for n_cl in list(range(1,324)):
#            print 'cluster number:', n_cl,
#                
#            filename=mycomp+'anaconda/pro/data/300/Galacticus/MDPL2_Galacticus_z0.00_region'+str(n_cl).zfill(4)+'.txt'
#            print filename
#            mydata_names=['z','hostid','haloid','orphan','x_pos','y_pos','z_pos','x_vel','y_vel','z_vel','mstar_spheroid','mstar_disk',\
#                          'mcold_spheroid','mcold_disk','mhot','mbh','sfr','sfr_spheroid','sfr_disk','mean_age_stars','mhalo','rest',\
#                          'Vmax','Vpeak','NFW_con','spin','MZstar_spheroid','Mzstar_disk', 'Mzgas_spheroid', 'Mzgas_disk','mzhot_halo',\
#                          'L_SDSS_dA_total_u','L_SDSS_dA_total_g','L_SDSS_dA_total_r','L_SDSS_dA_total_i','L_SDSS_dA_total_z',\
#                          'MAB_dA_total_u','MAB_dA_total_g','MAB_dA_total_r','MAB_dA_total_i','MAB_dA_total_z',\
#                          'rdisk','rbulge','rhalf_mass','mhot_outflow']
#            
#            mydata = mL.df_to_sarray(pd.read_csv(filename, skiprows=2, names=mydata_names, sep='  '))                
#            
#            if n_cl==1:
#                data=mydata
#            else:
#                data = np.append(data, mydata, axis=0)
#                
#            print 'size:', data.size
        
#        try:
#            if np.all(data['zcold']==-99.0):
#                print '\n+++++++++++++++++++\ncalculate gas-phase metallicity',
#                data['zcold']=8.69+np.log10(data['Mzgas']/(data['mcold']*0.0134))
#                #data['OH_gas_disk_bulge']=12.0+np.log10(data['OH_gas_disk_bulge']/16.0)
#                print '--> DONE!\n'
#        except:
#            print 'SAGE:',
#            #data['zcold'] = 8.69+np.log10(data['Mzgas_disk']/(data['mcold_disk']*0.0134))
#            data['zcold'] = 8.69+np.log10(data['zgas_disk']/0.0134)
#            print '--> DONE\n'
###
#        try:
#            if np.all(data['cgf']==-99.0):
#                print '\n+++++++++++++++++++\ncalculate gas-phase metallicity',
#                data['cgf']=data['mcold']/data['mstar']
#                print '--> DONE!\n'
#        except:
#            print 'SAGE:',
#            data['cgf'] = (data['mcold_disk']/data['mstar'])
#            print '--> DONE'
#
#        try:
#            mags='mAB_dA_total_'
#            bands=[mags+'u', mags+'g', mags+'r', mags+'i', mags+'z']
#            
#            if np.all(data['kcorr_u']==-99.0):
#    
#                #data=data[np.where(data['mstar']>1e11)[0][:]]
#                #print 'MSTAR CUT ', data.shape            
#                
#                print '\n+++++++++++++++++++\ncalculate k-correction with approximate formula!',
#                k_corr_array = mL.kcorrect_approx(data,
#                                          redshift,
#                                          bands)
#    
#    
#                print np.info(k_corr_array)
#                for k, band in enumerate(bands):
#                    print 'band:', band, 'num:', k
#                    data['kcorr_'+band[-1::]] = k_corr_array[band]
#                    print 'kcorrection: median =', "{0:.2f}".format(np.median(data['kcorr_'+band[-1::]])), '16th/84th', "{0:.2f}".format(np.percentile(data['kcorr_'+band[-1::]], 16)), '/', "{0:.2f}".format(np.percentile(data['kcorr_'+band[-1::]], 84))
#    
#                print '--> DONE!\n'
#        except:
#           pass
                                        
#        try:
#            data['mstar_IC']=data['mstar_IC']+data['mstar_disk']+data['mstar_spheroid']
#            print 'Mstar_IC+Mstar'            
#        except:
#            pass
            

##########################################################################################################################################################
#        data = data[np.where(data['mstar']>5e8)[:][0]]
#        #data = data[np.where(data['mstar']<=9e11)[:][0]]
#        print data.size        

#        preprocessing_only=True
#
#        print 'sample name:', 
#        mL.give_galaxy_type_info(data['orphan'])
#
#
#        mycut=myhdf5_filename[myhdf5_filename.find('CUT'):len(myhdf5_filename)-5]
#        print 'mycut:', mycut
#
#        for item in ['mstar','sfr','ssfr','mhalo_cents_200c']:
#            mydata = data[item][np.where(data[item]>0.0)[:][0]]
#            print item, 'min/max:', format(np.log10(np.nanmin(mydata)), '0.2f') , '/', format(np.log10(np.nanmax(mydata)), '0.2f')
#            print '\t\t', "{0:.2f}".format(np.log10(np.nanmedian(mydata))), '\t', "{0:.2f}".format(np.log10(np.nanmedian(mydata))-np.log10(np.nanpercentile(mydata, 32))), '/', "{0:.2f}".format(np.log10(np.nanpercentile(mydata, 68))-np.log10(np.nanmedian(mydata)))
#
#            stats=str(redshift)+'\t'+str("{0:.2f}".format(np.log10(np.nanmedian(mydata))))+'\t'+str("{0:.2f}".format(np.log10(np.nanmedian(mydata))-np.log10(np.nanpercentile(mydata, 32))))+'\t'+str("{0:.2f}".format(np.log10(np.nanpercentile(mydata, 68))-np.log10(np.nanmedian(mydata))))
#            if redshift==0.0 or redshift=='0.0':
#                myOutput.writeIntoFile(mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/'+self.myconfig_array['catname'+str(self.a)]+'_HOD-SFgalaxies_z-evolution_'+mycut+'_'+item+'.txt',
#                           [stats],
#                           myheader='HOD-SFGalaxies project z-evolution of properties! cat: '+self.myconfig_array['catname'+str(self.a)]+' cut: '+mycut+' prop: '+item+'\n(1) z (2) log10('+item+') (3) -dy (4) +dy',
#                           append_mytext=False,
#                           data_is_string=False,
#                           data_format='%s')
#            else:
#                myOutput.writeIntoFile(mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/'+self.myconfig_array['catname'+str(self.a)]+'_HOD-SFgalaxies_z-evolution_'+mycut+'_'+item+'.txt',
#                           stats+'\n',
#                           append_mytext=True,
#                           data_is_string=True,
#                           data_format='%s')                
#            

        #exit()
#        data=mL.downsample_SMF(data)
#        mL.give_galaxy_type_info(data['orphan'])

        
        #ngal=int((self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size']/0.6778)**3*0.008)
        #number density: Gal400
        #data = mL.choose_random_sample(data,3288480)#n=0.016
        #data = mL.choose_random_sample(data,2466360)#n=0.012
        #data = mL.choose_random_sample(data,1664240)#n=0.008
     
        #data = mL.choose_random_sample(data,ngal)        
        #print 'randomly chosen:', data.size


#        data = mL.test_haloids(data)
###
#        #data = data[np.where(data['orphan']==1)[0][:]]
#        print np.info(data)
#        mL.property_stats(data)
#        myOutput.writeIntoFile(mycomp+'anaconda/pro/data/Galacticus_1Gpc_run2/Galacticus_1Gpc_run2_z_0.56_same_haloid_Gal-dens.txt',
#                   data[['haloid','hostid','orphan','mhalo','mstar', 'haloid_test_against', 'hostid_test_against', 'orphan_test_against', 'env_512', 'env_1024', 'pop']],
#                   myheader= 'Galacticus 1Gpc z='+str(format(redshift, '0.2f'))+' CMASS sample of Gal2-dens sample galaxies with same hostid as Gal-dens \n(1) haloid (2) hostid (3) orphan (4) mhalo_200c [Msun] (5) mstar [Msun] (6) haloid_Gal2 (7) hostid_Gal2 (8) orphan_Gal2 (9) env_512 (10) env_1024 (11) pop',
#                   data_format="%i\t%i\t%i\t%0.8e\t%0.8e\t%i\t%i\t%i\t%i\t%i\t%i",
#                   mydelimiter='\t')         
#
#        myOutput.writeIntoFile(mycomp+'anaconda/pro/data/Galacticus_1Gpc_run2/Galacticus_1Gpc_run2_z_0.56_CMASS_down_sample3_haloid_from_Gal-dens.txt',
#                   data[['haloid','hostid']],
#                   myheader= 'Galacticus 1Gpc z='+str(format(redshift, '0.2f'))+' CMASS down sample3 (color cut: Guo+13)\n(1) haloid (2) hostid',
#                   data_format="%i\t%i",
#                   mydelimiter='\t')  
        
#        myOutput.writeIntoFile(mycomp+'anaconda/pro/data/Galacticus_400Mpc/Galacticus_400Mpc_z_0.55_CMASS_down_sample3_haloid.txt',
#                   data[['haloid','hostid']],
#                   myheader= 'Galacticus 400Mpc z='+str(format(redshift, '0.2f'))+' CMASS down sample3\n(1) haloid (2) hostid',
#                   data_format="%i\t%i",
#                   mydelimiter='\t')  

#        myOutput.writeIntoFile(mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/'+self.myconfig_array['catname'+str(self.a)]+'_z_'+str(format(redshift, '0.2f'))+'_rand-sample_n_0.016.txt',
#                   data[['mhalo_200c','orphan','mstar','sfr','ssfr','mhalo_cents_200c', 'mbh', 'mAB_total_cut_g_i', 'rhalf_mass', 'rdisk']],
#                   myheader= self.myconfig_array['catname'+str(self.a)]+' z='+str(format(redshift, '0.2f'))+' rand sample n=0.016 [Mpc-1] as TNG300-1 (Mstar>3e8 [h-1Msun]) \n(1) mhalo_200c [Msun] (2) orphan (3) mstar [Msun] (4) sfr [Msunyr-1] (5) ssfr [yr-1] (6) mhalo_cents_200c [Msun] (7) mbh [Msun] (8) mAB_total_cut_g_i [-] (9) rhalf_mass [Mpc] (10) rdisk/rbulge [-]',
#                   data_format="%0.10e\t%i\t%0.10e\t%0.10e\t%0.10e\t%0.10e\t%0.10e",
#                   mydelimiter='\t') 
      
#        exit()

#        data.sort(order=['mcold'], axis=0)      
#        data=data[data.size-530000:data.size]
#        print 'size sorted arry:', data.size
#        data = data[np.where(data['orphan']==0)[0][:]] 
#        
#        data = data[np.where(data['zcold']<=8.5)[0][:]]
#        data = data[np.where(data['mstar']>4e10)[0][:]]
#        print 'lowZ:', data.size      
#        
#
#        try:
#            if np.all(data['mhalo_200c']==-99.0):
#                data['NFW_con'] = mL.calculate_NFW_con_with_fit(data['mhalo'], 
#                                                       redshift,
#                                                       cosmology='Planck', 
#                                                       overdens='vir')
#                
#                data['mhalo_200c']=mL.convert_halo_mass(data['mhalo'], 
#                                                                    data['NFW_con'],
#                                                                    data['orphan'])
#                print '--> DONE!\n'
#        except:
#            pass
#
#
#        print 'mcold CUT3 centrals', data.size        
#        mL.property_stats(data)
#
#                   
##        mydata = data[np.where(data['zcold']>8.5)[0][:]]
##        mydata = mydata[np.where(mydata['mstar']>1e11)[0][:]]
##        print 'highZ:', mydata.size     
##        mL.property_stats(mydata)
#        
#                            
#
#        
#        mL.property_stats(mydata)
#        exit()        

##########################################################################################################################################################        

        if preprocessing_only==False:        
            if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_selectRegion']!='False':
                selectRegion(data, myheader, mylog, mylog_region)                                                                                                       
    
            if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filter_halomass_Sergio']!='False':
                data=extractHalomassSergio(data) 
    
            if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_createProgenitorList']!='False':
                print 'HERE 911'
                data=createProgenitorList(data, myheader, mylog)                                                                                                   
    
    
    
            if redshift=='False' or redshift==False:
                try:
                    redshift=float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_manual_input_redshift'])
                    print 'redshift reset to:', redshift
                except:
                    print 'redshift could not be reset! Current redshift value is:', redshift
    
            if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_convert_sph_to_cart_coords']!='False':        
                try:
                    data=mL.convert_coordinates(data)
                    print 'data new coords:', 'x_pos','y_pos','z_pos','RA','DEC','Z'
                    print data[0:10][['x_pos','y_pos','z_pos','RA','DEC','Z']]
                except:
                    print 'Converstion from spherical to euclidian coordinates is not possible!\n'
    
         
            data, myheader, mydataformat, mylog, check_size, sel_col_list = mL.filter_data_before_write2file(data,
                                                                                                            myconds_array,
                                                                                                            myconds_array,
                                                                                                            myheader,
                                                                                                            mylog)
    
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
    
                #print 'myheader:', myheader
                #print 'mydataformat:', str(mydataformat)
            else:
                print 'Sorry, but we ran out of galaxies!!!'
                print ' '
                print 'The selections of your data sample has zero entries! Try again with a different selection!!!'
                print 'Programm exiting ...'
                exit()
    
    
            #select randomly a subsample of data if the catalouge has more galaxies as your threshold!    
            if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_sample_info'].find('random')!=-1:            
                mylog+='\nCatalogue is randomly selected: Yes! ngalaxies selected? '+str(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_ngal_random_sample'])           
                data= mL.choose_random_sample(data, self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_ngal_random_sample'])
               
                    
            mylog+='\nngalaxies in the file: '+str(len(data[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['name0']]))+'\nlink and filename of catalogue: '+myhdf5_filename
          
            #writeCoordinates(data, myhdf5_filename)
            if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_OUTPUT_ASCII']!='False':				
                writeASCII(data[sel_col_list], myhdf5_filename[0:len(myhdf5_filename)-5]+'.txt', myheader, mydataformat)
            writeHDF5(data, myhdf5_filename)
            writeLOGFILE(mylog, myhdf5_filename[0:len(myhdf5_filename)-5]+'_log.txt')    
        else:
            print 'ONLY PREPROZESSING DONE!'
            self.myData2Process=data
            #print np.info(data)
            #print np.info(self.myData2Process)
            
                  
    def analyseTargetSelection(self,
                                redshift,
                                filename,
                                myconds_array,
                                plot_key):

    
        print 'analysisTargetSelection()', 'self.i:', self.i, 'self.a:', self.a
        print '++++++++++++++++++++++++++++++++++++++++++++++'
        print ' ' 

        def rawSAMHDF5Catalouge():
            print 'rawSAMHDF5Catalouge()'
            print '+++++++++++++++++++++'
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
            print 'correctedCatalouge()'
            print '+++++++++++++++++++++'
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

            def lumHisto():      
                bands=['i']
                cuts=[2.35]#,1.6,1.4,1.2]           
                
                method='1'
                for cut in cuts:
                    try:
                        mydata=data
                        i=0
                        for band in bands:
                            print 'analyseTargetSelection() --> Calculate Histogramm ...'
                            print '-------------------'
                            print ' '
                            print 'i:', i, 'name:',  myconds_array['name'+str(i)], 'sel_col:', myconds_array[myconds_array['name'+str(i)]+'_col_id'], mydata.shape
            
    
                            name = 'MAB_dA_total_'+band
                            #name='mstar'
                            print 'name:', name
                            
                            lum_dist=cd.luminosity_distance(0.55, **fidcosmo)
                            print fidcosmo
                            DM = 5.0*np.log10(lum_dist*1e6/10.0)
                            K_corr = -2.5*np.log10(1.0+0.555)
                           
                            if method=='1':
                                print 'method 1:'
                                mydata['MAB_dA_total_'+band] = -2.5*np.log10(mydata['L_SDSS_dA_total_'+band]) -(-2.5*np.log10(1.0+0.557)) + K_corr
                                mydata['mAB_dA_total_i'] = -2.5*np.log10(mydata['L_SDSS_dA_total_i']) + DM + K_corr
                                mydata['mAB_dA_total_r'] = -2.5*np.log10(mydata['L_SDSS_dA_total_r']) + DM + K_corr
                                #mydata['mAB_dA_total_g'] = -2.5*np.log10(mydata['L_SDSS_dA_total_g']) + DM + K_corr                          
    
                            elif method=='CM':
                                mydata['MAB_dA_total_'+band]=mydata['mAB_dA_total_'+band]-DM+5.0*np.log10(0.7/0.6777)-(-2.5*np.log10(1.0+0.55))
                                #print 18.055-DM+5.0*np.log10(0.7/0.6777)-(-2.5*np.log10(1.0+0.55))
    #                            mydata = mydata[(np.where(mydata['mAB_dA_total_i']>18.4999))[:][0]]
    #                            mydata = mydata[(np.where(mydata['mAB_dA_total_i']<18.5001))[:][0]] 
    #                            print mydata['mAB_dA_total_'+cut]
    #                            print mydata['MAB_dA_total_'+cut]                            
                            else:
    #                            print 'method 2:'
                                lum=mydata['L_SDSS_dA_total_'+cut]/(4.0*np.pi*(lum_dist*1e6*3.0857e16*100)**2)
                                mydata['mAB_dA_total_'+cut] = -2.5*np.log10(lum*(1+0.56)/3631e-23)                                                
                                
                                mydata['MAB_dA_total_'+cut] = mydata['mAB_dA_total_'+cut] - DM -5*np.log10(0.6777) -(-2.5*np.log10(1.0+0.55)) 
    #
    #                        print mydata
                            #x=mydata['mAB_dA_total_g']-mydata['mAB_dA_total_r']                        
                            #y=mydata['mAB_dA_total_r']-mydata['mAB_dA_total_i']
    
                            #test1
                            #prop=-0.26+(y+0.7)/0.92857-x                        
                            #test2
                            #prop=-0.3+(y+0.7)/0.92857-x
                            #test3
                            #prop=(y+1.05)/1.3-x
                            #mydata = mydata[np.where(prop<0.0)[:][0]]
    
                            #mydata = mydata[np.where(mydata['orphan']<2)[:][0]]
                            if cut=='Guo13':
                                print 'cut:', cut
                                mycut=0.679-0.082*(data['MAB_dA_total_i']+20.0)
                                mydata=mydata[np.where(mydata['mAB_dA_total_r']-mydata['mAB_dA_total_i']>mycut)[0][:]]
                            else:
                                print 'g-i>', cut
                                mydata = mydata[np.where(mydata['mAB_dA_total_g']-mydata['mAB_dA_total_i']>float(cut))[:][0]]                            
    
                              
                            print 'after color_cut:', mydata.shape
                                                 
                            if name.find('MAB')!=-1:
                                cond_min=-25
                                cond_max=-19
            
                            else:
                                cond_min=16
                                cond_max=23
            
                            data_histo = myData.selectData2Compute(mydata, 
                                                            selected_col=name, 
                                                            operator='<', 
                                                            condition=cond_max)
                                                    
                                                            
            
                            data_histo = myData.selectData2Compute(data_histo, 
                                                            selected_col=name, 
                                                            operator='>', 
                                                            condition=cond_min)
                            print 'after max:', data_histo.shape         
                            #data['mstar']*=1e10     
                            nbins=30
            
                            if name.find('AB')!=-1:
                                data_histo[name]=abs(data_histo[name])             
                            else:
                                data_histo[name]=np.log10(data_histo[name])
                            data_min = min(data_histo[name])
                            data_max = max(data_histo[name])
            
                            binsize = (data_max-data_min)/nbins
                            print 'data min/max:', data_min, '/', data_max, 'binsize:', binsize,
                            counts, edges = np.histogram(data_histo[name],nbins,weights=data_histo['weight_tot'])
    
                            histo=np.zeros((nbins, 7), dtype=np.float64)
                            if name.find('AB')!=-1:
                                if name.startswith('MAB'):
                                    histo[:,0]= -(edges[1:]+edges[:-1])/2
                                else:
                                    histo[:,0]= (edges[1:]+edges[:-1])/2
                            else:
                                histo[:,0]= (10**edges[1:]*10**edges[:-1])**0.5
                            
                            if method=='CM':
                                Vol=4.147e9 #CMASS DR12 
                            else:
                                Vol  = (float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'])/0.6777)**3
                                print 'Vol:', Vol, 'box_size:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'], 'data_min:', data_min, 'data_max:', data_max, 'binsize:', binsize
    
                                                       
                            histo[:,1]= counts/(Vol*binsize)
                            histo[:,4]=histo[:,1]/counts**0.5
                            histo[:,5]=histo[:,4]
                            histo[:,6]= counts
                            k=0
                            for cols in histo[:,6]:
                                #print cols
                                if cols<25.0:
                                    histo[k,2]=1.0
                                k+=1
                            
                            #print_5logh=' 5logh'
                            print_5logh=''
                            #selection='_-Kcorr_047_i_g-i_'+str(cut)+'_'
                            #selection='_g-i_gt_'+str(cut)+'_'                
                            selection=''
                            #key='_350Mpc_run_1235'
                            key='_'+str(cut)
                            key='_g-i_gt_'+str(cut)
                            
                            #filename=mycomp+'anaconda/pro/myRun/histos/analyseTargetSelection/histo_'+self.myconfig_array['catname'+str(self.a)]+'_v3_z_0.1_'+name+'-5logh_h-units.txt'
                            filename=mycomp+'anaconda/pro/myRun/histos/analyseTargetSelection/histo_'+self.myconfig_array['catname'+str(self.a)]+'_z_'+str(redshift)+'_'+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_tarsel_code']+'_'+selection+name+key+'.txt'
            
                            print 'name:', name, 'z:', redshift, 'format:', myconds_array[name+'_format'], 'unit:', myconds_array[name+'_unit']
                            myOutput.writeIntoFile(filename,
                                           histo,
                                           myheader= self.myconfig_array['catname'+str(self.a)]+' z '+str(redshift)+' '+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_tarsel_code']+' '+name+' histo\n(1) '+name+print_5logh+' ['+myconds_array[name+'_unit']+']\t(2) Phi [Mpc-3 dex-1]\t(3) -dx\t(4) dx\t(5) -dy\t(6) +dy\t(7) N count',
                                           data_format='%0.8e',
                                           mydelimiter='\t')
                
                            print '--------------'
                            print ' '
                            i+=1
                    except:
                        print 'cut:', cut, 'not possible!'    

            def caseSwitcher(histo_method):
            
                choose = {
                    'lum':   lumHisto,
                    'fast':  self.calcFastHisto,               
                    }
                    
                func = choose.get(histo_method)
                return func()
            
            caseSwitcher('fast')

                       
        if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_plotSample']=='True':
            print ' '
            print 'analyseTargetSelection() --> plot sample ...'
            print '-------------------'
            print ' '
           
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
                #print min(prop), max(prop)
                #data = data[np.where(prop<0.0)[:][0]]                
#                print 'name_x', name_x, 'name_y', name_y, 'name_z:', name_z

#                data = myData.selectData2Compute(data, 
#                                                 selected_col='orphan', 
#                                                 operator='==', 
#                                                 condition=2)                                
                
                try:
                    if self.myconfig_array['catname'+str(self.a)].find('SAdG_')!=-1:
                        data['mstar']=data['mstar_disk']+data['mstar_spheroid']
                        data['mcold']=data['mcold_disk']+data['mcold_spheroid']
                        data['Mzgas']=data['Mzgas_disk']+data['Mzgas_spheroid']

                    mstar_min=1e10
                    #mstar_min=10**(11.29-0.039+np.log10((0.73/0.6777)**2))                        
                    print 'Mstar > ', mstar_min, 'selection:', data.shape, '-->',
                    data = myData.selectData2Compute(data, 
                                                     selected_col='mstar', 
                                                     operator='>', 
                                                     condition=mstar_min)
                    print data.shape
                                 
                except:
                    try:
                        print 'Mbulge >', mstar_min, 'selection:', data.shape, '-->',
                        data = myData.selectData2Compute(data, 
                                                         selected_col='mstar_spheroid', 
                                                         operator='>', 
                                                         condition=mstar_min)
                        print data.shape
                    except:
                        pass


                if tarsel_plot_name.find('oh_mstar')!=-1:
                    print 'oh vs. mstar! name_y:', name_y
                    data[name_y]=12+np.log10(0.0356*data['Mzgas']/data['mcold'])
                        
                    myconds_array[name_y+'_name_in_plot'] = '12+ $\log_{10}$(OH)'  

                if tarsel_plot_name.find('_sfr_mstar')!=-1 or tarsel_plot_name.find('mcold_mstar')!=-1:
                    print 'sfr vs. mstar or mcold vs. mstar! name_x:', name_x, 'name_y:', name_y
                    myconds_array[name_y+'_name_in_plot']='$\log_{10}$ ($M_{Cold}/M_{_*}$)'

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
                    print 'ssfr > 1e-16 -->', data.shape                      
                    
                                                    
                if tarsel_plot_name.find('mbh_mstar')!=-1:
                    print 'mbh vs mbulge!'
                    try:
                        data['mstar_spheroid']+=data['mcold_spheroid']                    
                    except:
                        print 'no mcold bulge found!'
                        
                    data = myData.selectData2Compute(data, 
                                                     selected_col=name_y, 
                                                     operator='>', 
                                                     condition=0.0)                       

                if tarsel_plot_name.find('mhalo_mstar')!=-1:
                    print 'mstar vs. mhalo!'
                    data = myData.selectData2Compute(data, 
                                                    selected_col='orphan', 
                                                    operator='<', 
                                                    condition=1)
                    print 'data after central selection!', data.shape
                
                if tarsel_plot_name.find('zgas')!=-1:                           

#                    data = myData.selectData2Compute(data, 
#                                                     selected_col=name_x, 
#                                                     operator='>', 
#                                                     condition=1e11)
                    
#                    if self.myconfig_array['catname'+str(self.a)]=='Galacticus_1Gpc':
#                        print 'name_y:', name_y
#                        print 'Galacticus: zgas vs. mstar/cold! --> Zgas = Zgas_spheroid+Zgas_disk/Mcold (column Zgas_spheroid and Zgas_disk are wrongly named in the catalog it should be Mzgas_spheroid and Mzgas_disk)'
#                        data[name_y]+=data['zgas_disk']

                    data = myData.selectData2Compute(data, 
                                                    selected_col=name_y, 
                                                    operator='>', 
                                                    condition=0.0)
                    print 'Mzgas > 0.0 -->', data.shape

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
                    myconds_array[name_y+'_name_in_plot'] = '$Z_{Cold}$'

                    if tarsel_plot_name.find('zgas_mcold')!=-1:                  
                        data[name_x]/=data['mstar']
                        myconds_array[name_x+'_name_in_plot'] = '$\log_{10}$ ($M_{Cold}/M_{*}$)'   
                        
                        data = myData.selectData2Compute(data, 
                                                         selected_col=name_x, 
                                                         operator='>', 
                                                         condition=0.001)                        

                        data = myData.selectData2Compute(data, 
                                                         selected_col=name_x, 
                                                         operator='<', 
                                                         condition=100)

                if tarsel_plot_name.find('-')!=-1: 
                    print 'color vs color! name_x:', name_x, 'name_y:', name_y
                    print tarsel_plot_name
                    
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
#                    print '> 0.0 -->', data.shape
                                           
                    try:
                        print 'version:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['mAB_dA_total_r_col_id'], 'mAB_dA'
                        version = 'mAB_dA'
                    except:
                        try:
                            print 'version:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['MAB_dA_total_r_col_id'], 'MAB_dA'
                            version = 'MAB_dA'
                        except:
                            print 'version:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_id_col_array']['mAB_total_r_col_id'], 'mAB'
                            version = 'mAB'                    
                    #print tarsel_plot_name, version   
                    if tarsel_plot_name.find('r-i')!=-1 and tarsel_plot_name.find('g-r')==-1:
                        print 'HERE: r-i'
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
                        print 'HERE: g-r vs u-g'
                        name_y=version+'_total_g'
                        name_x=version+'_total_u' 
                        data[name_x]=data[version+'_total_u']-data[version+'_total_g']
                        data[name_y]=data[version+'_total_g']-data[version+'_total_r']
            
                    if tarsel_plot_name.find('g-r')!=-1 and tarsel_plot_name.find('r-i')!=-1:
                        print 'HERE: g-r vs r-i'
                        name_x=version+'_total_g' 
                        name_y=version+'_total_r'
                        data[name_x]-=data[version+'_total_r']
                        data[name_y]-=data[version+'_total_i']
                        
                    if tarsel_plot_name.find('g-i')!=-1:
                        print 'HERE: g-i'
                        name_y=version+'_total_g'
                        name_x=version+'_total_i'
                        data[name_y]-=data[version+'_total_i']
                        myconds_array[name_y+'_name_in_plot'] = '$g-i$'
                        
                        if tarsel_plot_name.find('mstar')!=-1:
                            data[name_x]=np.log10(data['mstar'])
                        
                    if tarsel_plot_name.find('u-r')!=-1:
                        print 'HERE: u-r',
                        name_y=version+'_total_u'
                        name_x=version+'_total_r' 
                        myconds_array[name_y+'_name_in_plot'] = '$u-r$'
                        print 'name_x:', name_x, 'name_y:', name_y
                        data[name_y]-=data[name_x]
                        data[name_x]-=5*np.log10(0.6777)
                        print 'check!'

                    if tarsel_plot_name.find('g-u')!=-1:
                        print 'HERE: g-u'
                        name_y=version+'_total_g'
                        name_x=version+'_total_u' 
                        myconds_array[name_y+'_name_in_plot'] = '$g-u$'                        

                if tarsel_plot_name.find('ssfr')!=-1 or tarsel_plot_name.find('zcold')!=-1 or tarsel_plot_name.find('mhalo')!=-1:
                    print 'name_x:',
                    if tarsel_plot_name.find('ssfr')!=-1 and tarsel_plot_name.find('mstar')==-1:
                        try:
                            name_x='ssfr'
                            data[name_x]=np.log10(data['ssfr'])
                        except:
                            name_x='sfr'
                            data[name_x]=np.log10(data['sfr']/data['mstar'])                            
                        print  name_x, 'name_y:', name_y
    
                        data = myData.selectData2Compute(data, 
                                                        selected_col=name_x, 
                                                        operator='>', 
                                                        condition=-16)
                        print 'ssfr > 1e-16 -->', data.shape                        
                  
                    if tarsel_plot_name.find('g-r')!=-1:
                        print 'HERE: g-r vs ssfr! name_x', name_x, 'name_y', name_y
                        data[name_y]=data[version+'_total_g']-data[version+'_total_r']
                        
                    if tarsel_plot_name.find('dmesa')!=-1:
                        print 'HERE: dmesa vs ssfr! name_x:', name_x, 'name_y:', name_y
                    
                    if tarsel_plot_name.find('mhalo')!=-1:
                        name_x='mhalo'
                        print 'HERE: ssfr vs mhalo! name_x:', name_x, 'name_y:', name_y
                        print 'only centrals!'
                        data = myData.selectData2Compute(data, 
                                                        selected_col='orphan', 
                                                        operator='==', 
                                                        condition=0)                           
                    
#                        print 'check environment!'
#                        data = myData.selectData2Compute(data, 
#                                                        selected_col='env_1024', 
#                                                        operator='==', 
#                                                        condition=2)

                        if tarsel_plot_name.find('mstar')!=-1:
                            data['mstar']= np.log10(data['mstar'])                              
                        
                        data[name_x]= np.log10(data['mhalo_200c'])

                    if tarsel_plot_name.find('mstar')!=-1 and tarsel_plot_name.find('mhalo')==-1:
                        name_x='mstar'                                              
                        data['mstar']= np.log10(data['mstar'])

                    if tarsel_plot_name.find('ssfr_mstar')!=-1:                       
                        name_y='ssfr'                                              
                        try:
                            data[name_y]=np.log10(data['ssfr'])
                        except:
                            print 'except 1893'
                            name_y='sfr'
                            data[name_y]=np.log10(data['sfr']/data['mstar'])                        

                        name_x='mstar'                                              
                        data['mstar']=np.log10(data['mstar']) 
                        
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
                        myconds_array[name_x+'_name_in_plot'] = '$\log_{10}$ ($M_{Cold}/M_{*}$)'  

                if tarsel_plot_name.find('dmesa')!=-1:
                    print 'HERE: dmesa vs.', 
                    if tarsel_plot_name.find('dmesa_i')!=-1:
                        print 'i!'
                    elif tarsel_plot_name.find('dmesa_mstar')!=-1:
                        print 'mstar!'                     
                    myconds_array[name_y+'_name_in_plot'] = '$d_{\perp}$'

                if tarsel_plot_name.find('_mstar')!=-1 and tarsel_plot_name.find('_mstar_sph')==-1:
                    print 'set x-axis to mstar! x_name:', name_x
                    myconds_array[name_x+'_name_in_plot'] = '$\log_{10}$ ($M_{_*}$ [$M_{\odot}$])'
                    #data[name_x] = np.log10(data['mstar'])

                if tarsel_plot_name.find('_rdisk')!=-1 and tarsel_plot_name.find('_rbulge')!=-1:
                    print 'rdisk vs rbulge in kpc! x_name:', name_x, 'name_y:', name_y
                    
                    data[name_x]= np.log10(data['rbulge']*1e3)
                    data[name_y]= np.log10(data['rdisk']*1e3)

                elif tarsel_plot_name.find('rdisk')!=-1 and tarsel_plot_name.find('mstar')!=-1:
                    print 'rdisk in kpc vs mstar! x_name:', name_x, 'name_y:', name_y
                    data[name_x]= np.log10(data['mstar'])
                    data[name_y]= np.log10(data['rdisk']*1e3)

                elif tarsel_plot_name.find('rhalfdisk')!=-1 and tarsel_plot_name.find('mstar')!=-1:
                    print 'rhalf_disk in kpc vs mstar! x_name:', name_x, 'name_y:', name_y
                    data[name_x]= np.log10(data['mstar'])
                    data[name_y]= np.log10(data['rhalf_disk']*1e3)

                elif tarsel_plot_name.find('rbulge_')!=-1 and tarsel_plot_name.find('_rhalf')!=-1:
                    print 'rhalf vs rbulge in kpc! x_name:', name_x, 'name_y:', name_y                                  
                    data[name_x]= np.log10(data['rbulge']*1e3)
                    data[name_y]= np.log10(data['rhalf_mass']*1e3)
                    
                elif tarsel_plot_name.find('rhalfbulge')!=-1 and tarsel_plot_name.find('mstar')!=-1:
                    print 'rhalf_bulge in kpc vs mstar! x_name:', name_x, 'name_y:', name_y
                    data[name_x]= np.log10(data['mstar'])
                    data[name_y]= np.log10(data['rhalf_bulge']*1e3)                    

                elif tarsel_plot_name.find('rhalfmass')!=-1 and tarsel_plot_name.find('mstar')!=-1:
                    print 'rhalf_mass in kpc vs mstar! x_name:', name_x, 'name_y:', name_y
                    data[name_x]= np.log10(data['mstar'])
                    data[name_y]= np.log10(data['rhalf_mass']*1e3)

                elif tarsel_plot_name.find('bdisk_bbulge')!=-1:
                    name_x='angM_spheroid'
                    name_y='angM_disk'
                    print 'b_disk vs b_bulge! x_name:', name_x, 'name_y:', name_y
                    data[name_y]=np.log10(data['angM_disk']/(data['mstar_disk']+data['mcold_disk']))-2/3.0*np.log10(data['mstar_disk']+data['mcold_disk'])
                    data[name_x]=np.log10(data['angM_spheroid']/(data['mstar_spheroid']+data['mcold_spheroid']))-2/3.0*np.log10(data['mstar_spheroid']+data['mcold_spheroid'])        


                elif tarsel_plot_name.find('angM_mbar')!=-1:
                    name_x='mstar'
                    name_y='angM_disk'
                    print 'jbar vs mbar! x_name:', name_x, 'name_y:', name_y
                    data[name_y]= np.log10((data['angM_disk']+data['angM_spheroid'])/(data['mstar']+data['mcold']))
                    data[name_x]= np.log10(data['mstar']+data['mcold'])


                elif tarsel_plot_name.find('angMdisk')!=-1 and tarsel_plot_name.find('mbar')!=-1:
                    name_x='mstar'
                    print 'jdisk vs mbar disk! x_name:', name_x, 'name_y:', name_y
                    data[name_y]= np.log10(data['angM_disk']/(data['mstar_disk']+data['mcold_disk']))
                    data[name_x]= np.log10(data['mstar_disk']+data['mcold_disk'])

                elif str(tarsel_plot_name).find('angMspheroid')!=-1 and tarsel_plot_name.find('mbar')!=-1:
                    name_x='mstar'
                    print 'jbluge vs mbar bulge! x_name:', name_x, 'name_y:', name_y
                    data[name_y]= np.log10(data['angM_spheroid']/(data['mstar_spheroid']+data['mcold_spheroid']))
                    data[name_x]= np.log10(data['mstar_spheroid']+data['mcold_spheroid'])                
                    
                for name in [name_x,name_y,name_z]:
                    print 'set axis labels! name:', name
                    if name=='mstar':
                        if tarsel_plot_name.find('mbar')!=-1:
                            myconds_array[name+'_name_in_plot']='$\log_{10}$ ($M_{bar}$ [$M_{\odot}$])'
                        else:
                            myconds_array[name+'_name_in_plot']='$\log_{10}$ ($M_{_*}$ [$M_{\odot}$])'
                    elif name=='Mzstar':
                        myconds_array[name+'_name_in_plot']='$\log_{10}$ ($mass metals M_{_*}$ [$M_{\odot}$])'                        
                    elif name=='mhalo':
                        myconds_array[name+'_name_in_plot']='$\log_{10}$ ($M_{200c}$ [$M_{\odot}$])'
                    elif name=='mbh':
                        myconds_array[name+'_name_in_plot']='$\log_{10}$ ($M_{BH}$ [$M_{\odot}$])'
                    elif name=='mstar_spheroid':
                        #myconds_array[name+'_name_in_plot']='$\log_{10}$ ($M_{Spheroid}$+$M_{Gas_{Bulge}}$ [$M_{\odot}$])'
                        myconds_array[name+'_name_in_plot']='$\log_{10}$ ($M_{Bulge}$ [$M_{\odot}$])'                          
                    elif name=='MAB_dA_total_i':                       
                        myconds_array[name+'_name_in_plot']='$M_{AB_i}$'
                    elif name=='MAB_dA_total_r':                       
                        myconds_array[name+'_name_in_plot']='$M_{AB_r}$'
                    elif name=='mAB_dA_total_i':                       
                        myconds_array[name+'_name_in_plot']='$m_{AB_i}$'                        
                    elif name=='mAB_dA_total_r':                       
                        myconds_array[name+'_name_in_plot']='$m_{AB_r}$'
                    elif name=='sfr' and str(tarsel_plot_name).find('ssfr_mstar')==-1:
                        myconds_array[name+'_name_in_plot'] = '$\log_{10}$ ($SFR$ [$M_{\odot}$ $yr^{-1}$])'
                    elif name=='ssfr' or str(tarsel_plot_name).find('ssfr_mstar')!=-1:                 
                        myconds_array[name+'_name_in_plot'] = '$\log_{10}$ ($sSFR$ [$yr^{-1}$])'                            
                    elif name=='rdisk':
                        myconds_array[name+'_name_in_plot'] = '$\log_{10}$ ($r_{disk}$ [$kpc$])'                       
                    elif name=='rbulge':
                        myconds_array[name+'_name_in_plot'] = '$\log_{10}$ ($r_{bulge}$ [$kpc$])'
                    elif name=='rhalf_mass':
                        myconds_array[name+'_name_in_plot'] = '$\log_{10}$ ($r_{1/2}$ [$kpc$])'                         
                    elif name=='rhalf_disk':
                        myconds_array[name+'_name_in_plot'] = '$\log_{10}$ ($r^{disk}_{1/2}$ [$kpc$])'  
                    elif name=='rhalf_bulge':
                        myconds_array[name+'_name_in_plot'] = '$\log_{10}$ ($r^{bulge}_{1/2}$ [$kpc$])'
                    elif name=='spinParameter':
                        myconds_array[name+'_name_in_plot'] = '$\log_{10}$ (spin parameter [-])'
                    elif name=='zcold':
                        myconds_array[name+'_name_in_plot'] = '$Z_{Cold}$'
                    elif name=='angM_disk':
                        if tarsel_plot_name.find('angM_mbar')!=-1:
                            myconds_array[name+'_name_in_plot'] = '$j_{bar}$ [$kpc$ $km$ $s^{-1}$]'                          
                        else:
                            myconds_array[name+'_name_in_plot'] = '$j^{disk}_{bar}$ [$kpc$ $km$ $s^{-1}$]'
                            if tarsel_plot_name.find('bdisk')!=-1:
                                myconds_array[name+'_name_in_plot'] = '$b^{disk}_{bar}$' 
                    elif name=='angM_spheroid':
                        myconds_array[name+'_name_in_plot'] = '$j^{bulge}_{bar}$ [$kpc$ $km$ $s^{-1}$]'
                        if tarsel_plot_name.find('bbulge')!=-1:
                            myconds_array[name+'_name_in_plot'] = '$b^{bulge}_{bar}$'                          
                        
#                    if tarsel_plot_name.find('_mstar')!=-1 and tarsel_plot_name.find('_mstar_sph')==-1:
#                        myconds_array[name_x+'_name_in_plot'] = '$\log_{10}$ ($M_{_*}$ [$M_{\odot}$])'
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
#                    elif name=='zgas_spheroid':
#                        myconds_array[name_y+'_name_in_plot'] = '$\log$ ($Z_{Gas}$)'                         
                    #a+=1

#                y=np.log10(data['ssfr'])
#                x=np.log10(data['sfr'])
#                prop=(y+11.16)/1.12-x
#                data = data[np.where(prop>0.0)[:][0]] 

                mask=np.where(np.isfinite(data[name_x]))
                data=data[mask[:][0]]
                mask=np.where(np.isfinite(data[name_y]))
                data=data[mask[:][0]]                
                
                print 'min/max x:', min(data[name_x]), max(data[name_x])            
                print 'min/max y:', min(data[name_y]), max(data[name_y])
                print 'name y/x:', name_y, name_x, 'name_in_plot y/x:', myconds_array[name_y+'_name_in_plot'], myconds_array[name_x+'_name_in_plot'], 'data.shape:', data.shape
                #if len(data[0])>1e6:
                    #data = mL.choose_random_sample(data, 100000)
                if np.all(data[name_weights]==-99.0) or np.all(data[name_weights]==0.0):
                    print 'set weights to 1!'
                    data[name_weights]=np.ones((data.size,), dtype=np.int8)
   
                    #print data[name_weights]
                    

                i+=1
        
        print 'CHECK:', self.myconfig_array['catname'+str(self.a)]+' z='+str(redshift),  myconds_array[name_x+'_name_in_plot'], myconds_array[name_y+'_name_in_plot']
        try:
            return data[[name_x, name_y, name_z, name_weights]], self.myconfig_array['catname'+str(self.a)]+' z='+str(redshift), name_x, name_y, name_z, myconds_array[name_x+'_name_in_plot'], myconds_array[name_y+'_name_in_plot'], myconds_array[name_z+'_name_in_plot']
        except:
            return data[[name_x, name_y, name_weights]], self.myconfig_array['catname'+str(self.a)]+' z='+str(redshift), name_x, name_y, name_z, name_weights, myconds_array[name_x+'_name_in_plot'], myconds_array[name_y+'_name_in_plot'], myconds_array[name_z+'_name_in_plot'], myconds_array[name_weights+'_name_in_plot']




    def plotXY(self,
                plot_key,
                filename):
        print ' '        
        print 'plotXY()', 'self.i:', self.i, 'self.a:', self.a
        print '++++++++++++++++++++++++++++++++++++++++++++++'
        print ' '

        def calc_dmesa_cut():
            print 'filename:', filename
            
            
            

        def numFrac_short():
            print 'filename:', filename
            
            cuts = [0,1,2]
            
            data = myData.selectData2Compute(self.myData2Process, 
                                 selected_col='mstar', 
                                 operator='>', 
                                 condition=5e9)
            checknum=0
            for cut in cuts:
                print data.shape
                print 'select mstar > ', cut, ' -->',
                data_sorted = myData.selectData2Compute(data, 
                                                 selected_col='orphan', 
                                                 operator='==', 
                                                 condition=cut)
                print 'num: ', data_sorted.shape
                checknum+=data_sorted.size
                print 'total num:', checknum
                
#                data_sorted = myData.selectData2Compute(data_sorted, 
#                                                 selected_col='orphan', 
#                                                 operator='==', 
#                                                 condition=2)
#                print 'orphans:', data_sorted.shape

        def binup2props():

            data_filtered = myData.selectData2Compute(self.myData2Process, 
                                                         selected_col='sfr', 
                                                         operator='>', 
                                                         condition=0.0)            
            
            data_filtered['sfr_spheroid_inst']+=data_filtered['sfr_quies_inst']

            data_filtered = myData.selectData2Compute(data_filtered, 
                                                         selected_col='sfr_spheroid_inst', 
                                                         operator='>', 
                                                         condition=0.0)  

            print data_filtered.shape
            
            histo, binsiye = myFuncs.binUp(data_filtered[['sfr','sfr_spheroid_inst']],
                                           50,
                                           use_MAD=False,
                                           norm_by_binsize=False,
                                           log10bin=True)

            header='\n(1) sfr [Msunyr-1] (2) sfr inst (spheroid+quies) [Msunyr-1] (3) dx \t(4) -dy\t(5) +dy\t(6) N count'    

            myOutput.writeIntoFile(filename[0:len(filename)-5]+'_sfr_inst.txt',
                                   histo,
                                   myheader=self.myconfig_array['catname'+str(self.a)]+' z='+str(float("{0:.2f}".format(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'])))+' cumulative: NO'+header,
                                   data_format="%0.8e",
                                   mydelimiter='\t')   



        def numFrac_mags():
                          
            data_sorted = myData.selectData2Compute(self.myData2Process, 
                                             selected_col='MAB_dA_total_r', 
                                             operator='>', 
                                             condition=-22)
            
            data_sorted = myData.selectData2Compute(data_sorted, 
                                             selected_col='MAB_dA_total_r', 
                                             operator='<', 
                                             condition=-21)
            
            print self.myconfig_array['catname'+str(self.a)], 'total ngal:', len(data_sorted)
            print '-22<Mr<-21'
            data_type0=data_sorted[np.where(data_sorted['haloid']==data_sorted['hostid'])]
            
            print 'ngal          centrals:', len(data_type0), 'sats:', len(data_sorted)-len(data_type0) 
            print 'fraction [%]: centrals:', 100.0/len(data_sorted)*len(data_type0), 'sats', 100.0-100.0/len(data_sorted)*len(data_type0)
    
        def testHaloid():
            
            data_type0=self.myData2Process[np.where(self.myData2Process['orphan']==0)]
            data_sats=self.myData2Process[np.where(self.myData2Process['orphan']==1)]
            data_orphans=self.myData2Process[np.where(self.myData2Process['orphan']==2)]
            print 'ngal          centrals:', len(data_type0), 'sats+orphans:', len(data_sats)+len(data_orphans), 'sats:', len(data_sats), 'orphans:', len(data_orphans)
            print 'fraction [%]: centrals:', 100.0/len(self.myData2Process)*len(data_type0), 'sats+orphans:', 100.0/len(self.myData2Process)*(len(data_sats)+len(data_orphans)), 'sats:', 100.0/len(self.myData2Process)*len(data_sats), 'orphans:', 100.0/len(self.myData2Process)*len(data_orphans) 
        
            test = np.in1d(data_orphans['haloid'], data_centrals['haloid'])
    
            print 'len test:', len(test)
            print 'where test==False:', np.where(test==False)
            print data_orphans['haloid'][np.where(test==False)]
            orphan_error=data_orphans['haloid'][np.where(test==False)]
            print 'len orphan_error:', len(orphan_error)
            data_orphans=data_orphans[np.where(test==False)]
            print 'len orphans not set:', len(data_orphans)
            
            import numpy.lib.recfunctions as rcfuncs
            data = rcfuncs.append_fields([data], ['x_pos_centrals','y_pos_centrals','z_pos_centrals'], [data['x_pos'],data['y_pos'],data['z_pos']], usemask=False)

            myerror_log='ngal total: '+str(len(self.myData2Process))+' MainHalo==HostHalo: '+str(len(data_centrals))+' type0: '+str(len(data_type0))+' sats: '+str(len(data_sats))+' orphans: '+str(len(data_orphans))+'\n'#--------------------------------------------------------------------------------------\n\n'
            myerror_log+='all:                MIN/MAX values: mstar: '+str("{0:.2e}".format(min(data['mstar'])))+'/'+str("{0:.2e}".format(max(data['mstar'])))+' mcold: '+str("{0:.2e}".format(min(data['mcold'])))+'/'+str("{0:.2e}".format(max(data['mcold'])))+' sfr: '+str("{0:.2e}".format(min(data['sfr'])))+'/'+str("{0:.2e}".format(max(data['sfr'])))+'\n'            
            myerror_log+='MainHalo==HostHalo: MIN/MAX values: mstar: '+str("{0:.2e}".format(min(data_centrals['mstar'])))+'/'+str("{0:.2e}".format(max(data_centrals['mstar'])))+' mcold: '+str("{0:.2e}".format(min(data_centrals['mcold'])))+'/'+str("{0:.2e}".format(max(data_centrals['mcold'])))+' sfr: '+str("{0:.2e}".format(min(data_centrals['sfr'])))+'/'+str("{0:.2e}".format(max(data_centrals['sfr'])))+'\n'         
            myerror_log+='centrals type 0:    MIN/MAX values: mstar: '+str("{0:.2e}".format(min(data_type0['mstar'])))+'/'+str("{0:.2e}".format(max(data_type0['mstar'])))+' mcold: '+str("{0:.2e}".format(min(data_type0['mcold'])))+'/'+str("{0:.2e}".format(max(data_type0['mcold'])))+' sfr: '+str("{0:.2e}".format(min(data_type0['sfr'])))+'/'+str("{0:.2e}".format(max(data_type0['sfr'])))+'\n'           
            myerror_log+='sats type 1:        MIN/MAX values: mstar: '+str("{0:.2e}".format(min(data_sats['mstar'])))+'/'+str("{0:.2e}".format(max(data_sats['mstar'])))+' mcold: '+str("{0:.2e}".format(min(data_sats['mcold'])))+'/'+str("{0:.2e}".format(max(data_sats['mcold'])))+' sfr: '+str("{0:.2e}".format(min(data_sats['sfr'])))+'/'+str("{0:.2e}".format(max(data_sats['sfr'])))+'\n'            
            myerror_log+='orphans type 2:     MIN/MAX values: mstar: '+str("{0:.2e}".format(min(data_orphans['mstar'])))+'/'+str("{0:.2e}".format(max(data_orphans['mstar'])))+' mcold: '+str("{0:.2e}".format(min(data_orphans['mcold'])))+'/'+str("{0:.2e}".format(max(data_orphans['mcold'])))+' sfr: '+str("{0:.2e}".format(min(data_orphans['sfr'])))+'/'+str("{0:.2e}".format(max(data_orphans['sfr'])))
            myerror_log+='\n\n'+'ID in file'.ljust(12)+'MainHaloID'.ljust(15)+'HostHaloID'.ljust(15)+'type'.ljust(6)+'mstar'.ljust(12)+'mcold'.ljust(12)+'sfr'.ljust(12)+'X'.ljust(12)+'Y'.ljust(12)+'Z'.ljust(12)+'\n'            


        def numFrac():        
            cuts = [0, 1e9,1e10,1e11,1e12]
            
            for cut in cuts:
                print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'        
                print 'select mstar > 1e9 -->',
                data_sorted = myData.selectData2Compute(self.myData2Process, 
                                                 selected_col='mstar', 
                                                 operator='>', 
                                                 condition=cut)
                print data_sorted.shape
                
        
                print ' '
                print 'all:'
                print 'mstar min/max:',"{0:.2e}".format(min(data_sorted['mstar'])), "{0:.2e}".format(max(data_sorted['mstar']))
                print 'mcold min/max:',"{0:.2e}".format(min(data_sorted['mcold_disk'])), "{0:.2e}".format(max(data_sorted['mcold_disk']))
                print 'sfr   min/max:',"{0:.2e}".format(min(data_sorted['sfr'])), "{0:.2e}".format(max(data_sorted['sfr']))
                print ' '
                print 'centrals:'
                print 'mstar min/max:',"{0:.2e}".format(min(data_sorted[np.where(data_sorted['orphan']==0)]['mstar'])), "{0:.2e}".format(max(data_sorted[np.where(data_sorted['orphan']==0)]['mstar']))
                print 'mcold min/max:',"{0:.2e}".format(min(data_sorted[np.where(data_sorted['orphan']==0)]['mcold_disk'])), "{0:.2e}".format(max(data_sorted[np.where(data_sorted['orphan']==0)]['mcold_disk']))
                print 'sfr   min/max:',"{0:.2e}".format(min(data_sorted[np.where(data_sorted['orphan']==0)]['sfr'])), "{0:.2e}".format(max(data_sorted[np.where(data_sorted['orphan']==0)]['sfr']))        
                print ' '        
                print 'sats:'
                print 'mstar min/max:',"{0:.2e}".format(min(data_sorted[np.where(data_sorted['orphan']==1)]['mstar'])), "{0:.2e}".format(max(data_sorted[np.where(data_sorted['orphan']==1)]['mstar']))
                print 'mcold min/max:',"{0:.2e}".format(min(data_sorted[np.where(data_sorted['orphan']==1)]['mcold_disk'])), "{0:.2e}".format(max(data_sorted[np.where(data_sorted['orphan']==1)]['mcold_disk']))
                print 'sfr   min/max:',"{0:.2e}".format(min(data_sorted[np.where(data_sorted['orphan']==1)]['sfr'])), "{0:.2e}".format(max(data_sorted[np.where(data_sorted['orphan']==1)]['sfr']))             
                print ' '        
                print 'orphans:'
                print 'mstar min/max:',"{0:.2e}".format(min(data_sorted[np.where(data_sorted['orphan']==2)]['mstar'])), "{0:.2e}".format(max(data_sorted[np.where(data_sorted['orphan']==2)]['mstar']))
                print 'mcold min/max:',"{0:.2e}".format(min(data_sorted[np.where(data_sorted['orphan']==2)]['mcold'])), "{0:.2e}".format(max(data_sorted[np.where(data_sorted['orphan']==2)]['mcold']))
                print 'sfr   min/max:',"{0:.2e}".format(min(data_sorted[np.where(data_sorted['orphan']==2)]['sfr'])), "{0:.2e}".format(max(data_sorted[np.where(data_sorted['orphan']==2)]['sfr'])) 
                print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++' 
                print ' '
  

        def filterPortsmouth():
            #Filter Portsmouth BOSS D12 Starforming catalog
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            print max(self.myData2Process[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_Z']]), min(self.myData2Process[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_Z']])
    
            self.myData2Process = mL.filter_data_before_analysis(self.myData2Process, 
                                                             0.0, 
                                                             0.2, 
                                                             self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_Z'])
            
            print max(self.myData2Process[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_Z']]), min(self.myData2Process[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_Z']])
            print self.myData2Process.shape
            filename=mycomp+'anaconda/pro/data/SDSS_DR8/SDSS_DR8_z_0.1_Portsmouth_starforming_Salpeter.hdf5'
            myOutput.write2HDF5(self.myData2Process,
                                self.mycond_array,
                                self.myconfig_array,
                                filename,
                                self.myconfig_array['catname'+str(self.a)],
                                volume=None,
                                redshift='0.0-0.2',
                                scale_factor=None)

        def test_kcorrect_blanton():
            
            mL.kcorrect_blanton(0.1)
            
        def test_MDreadin():
            import read_MDGal as MDGal
            MDGal.readMDGal()

        def plot_scatter():

            import matplotlib.pyplot as plt
            
            fig = plt.figure(figsize=(12,10), dpi=150)
            myax = fig.add_subplot(111, projection='3d') 

            #plt.scatter3D(xs, ys, zs=0, zdir='z', s=20, c=None, depthshade=True, *args, **kwargs)                  
            myax.scatter3D(self.myData2Process['y_pos'],
                         self.myData2Process['z_pos'],
                         self.myData2Process['z_pos'],
                        s=0.01,
                        alpha=0.1,
                        marker='.',
                        zdir='z',
                        color='r')
            plt.show()
            plt.savefig(mycomp+'anaconda/pro/myRun/plots/plotXY/scatter_plot.png', format='png', rasterized=True, dpi=100, bbox_inches='tight', pad_inches=0.05) 
            #from matplotlib.backends.backend_pdf import PdfPages
            #pp = PdfPages(mycomp+'anaconda/pro/myRun/plots/plotXY/scatter_plot.pdf')
            #plt.savefig(pp, format='pdf', rasterized=True, dpi=20, pad_inches=0.05, transparent=True)
            #pp.close()
            
        def plot_only():

            import matplotlib.pyplot as plt

            histo_data=myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/SMF/SMF_CMASS_SPALL_z_0.5_0.6_spall_por_PS_DR12v4_compl_sample_wg_fixed_30bins.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)        
            #print histo_data[0,:]
            histo_model=myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/SMF/SMF_Galacticus_1Gpc_z_0.56_new_mags_fixed_30bins.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)        
      
            #print histo_model[0,:]
        
        
            frac_data=histo_data[:,1]/histo_model[:,1]

            
            fig = plt.figure(figsize=(12,10), dpi=150)
            myax = fig.add_subplot(111) 

            #plt.scatter3D(xs, ys, zs=0, zdir='z', s=20, c=None, depthshade=True, *args, **kwargs)                  
            myax.plot(np.log10(histo_data[:,0]),
                        frac_data,
                        lw=2.0,
                        marker='',
                        ls='-',
                        color='r')
            #plt.show()
            #plt.savefig(mycomp+'anaconda/pro/myRun/plots/plotXY/plot_only.png', format='png', rasterized=True, dpi=100, bbox_inches='tight', pad_inches=0.05) 
            from matplotlib.backends.backend_pdf import PdfPages
            pp = PdfPages(mycomp+'anaconda/pro/myRun/plots/plotXY/plot_oply.pdf')
            plt.savefig(pp, format='pdf', rasterized=True, dpi=20, pad_inches=0.05, transparent=True)
            pp.close()

        def test_300():
            import pandas as pd
            print 'here 300'
            

            
            for n_cl in list(range(1,325)):
                print 'cluster number:', n_cl,
                    
                filename=mycomp+'anaconda/pro/data/300/Galacticus/MDPL2_Galacticus_z0.00_region'+str(n_cl).zfill(4)+'.txt'
                print filename
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


                print data['ID'][0]
                data['mstar']+=data['mstar_disk']
                #print 'data:', data[['haloid','hostid','orphan','mhalo','mstar']]
                
                
                
                if n_cl==1:
                    mydata=data
                else:
                    mydata = np.append(mydata, data, axis=0)
                    
                

                 
            print 'size:', mydata.size
            #print mydata
            
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
                'calc_dmesa_cut':   calc_dmesa_cut,
                'numFrac':          numFrac,
                'numFrac_short':    numFrac_short,
                'filterPortsmouth': filterPortsmouth,
                'numFrac_mags':     numFrac_mags,
                'testHaloid':       testHaloid,
                'binup2props':      binup2props,
                'test_kcorrect_blanton': test_kcorrect_blanton,
                'test_MDreadin': test_MDreadin,
                'plot_scatter': plot_scatter,
                'plot_only': plot_only,
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
                    
        print ' ' 
        print ' '                    
        print 'SFR vs Redshift(): -->  self.i:', self.i, 'self.a:', self.a
        print '++++++++++++++++++++++++++++++++++++++++++++++'   
        
        self.histo_data_ZvsSFR = data
        print 'mstar min/max:', mycond_min, '/', mycond_max, '\t',
        myfiltered_data = mL.filter_data_before_analysis(self.myData2Process, 
                                                         mycond_min, 
                                                         mycond_max, 
                                                         'mstar')
        
        print 'sfr min/max:', sfr_cut_min, '/', sfr_cut_max, '\t',
        myfiltered_data = mL.filter_data_before_analysis(myfiltered_data, 
                                                         sfr_cut_min, 
                                                         sfr_cut_max, 
                                                         'sfr')
        print 'ngal after selection -->', myfiltered_data.shape
        
        print 'mstar min/max:', format(min(myfiltered_data['mstar']),'0.4e') ,'/', format(max(myfiltered_data['mstar']),'0.4e')
        print 'sfr min/max:', format(min(myfiltered_data['sfr']),'0.4e') ,'/', format(max(myfiltered_data['sfr']),'0.4e')        
        
        if ssfr_cut!=None:
            import numpy.lib.recfunctions as rcfuncs
            myfiltered_data = rcfuncs.append_fields([myfiltered_data], ['ssfr'] ,[myfiltered_data['sfr']], usemask=False)
            #print myfiltered_data[0:5]
            print 'expand -->', np.info(myfiltered_data),
            
            myfiltered_data['ssfr']/=myfiltered_data['mstar']
            #print myfiltered_data[0:5]
    
            #MD-paper ssfr >1e-11
            myfiltered_data = myfiltered_data[np.where(myfiltered_data['ssfr']>1e-11)[0]]
            #MD-paper ssfr >2.17e-11
            #myfiltered_data = myfiltered_data[np.where(myfiltered_data['ssfr']>2.17e-11)[0]]        
            
            print 'ssfr max:', max(myfiltered_data['ssfr']), 'min:', min(myfiltered_data['ssfr'])
    
            print 'after selection:', myfiltered_data.shape
        #Attention!
        #This routine is only running if z=0 is calculated first. Although the normalisation of the volume in wrong!
        if self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(0)]+'_snapid'+str(0)]['z']!=0.0:
            print 'Self construction sequence initiated ..... you have 5 sec to overthink your choices of redshifts!'
            print ' '
            print '... in order to abort self construct choose z=0 as the first redshift to calculate cosmic star formation history!'
            print 'This routine is only running if z=0 is calculated first or otherwise the normalisation of the volume will be wrong!'
            exit()

        self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +4] = self.volume
        print 'self.volume:', format(self.volume, '0.2e'), 'sum sfr:', format(np.sum(myfiltered_data['sfr']),'0.2e'), ' / volume z=0:', format(self.histo_data_ZvsSFR[0, data_offset*(self.a) +4], '0.2e'), '=', format(np.sum(myfiltered_data['sfr'])/self.histo_data_ZvsSFR[0, data_offset*(self.a) +4], '0.2e')  
        
        self.histo_data_ZvsSFR[self.i, data_offset*(self.a)] = self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z']
        self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +1] = np.sum(myfiltered_data['sfr'])/self.histo_data_ZvsSFR[0, data_offset*(self.a) +4]    
        self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +2] = np.sum(myfiltered_data['sfr']) 
    
        #normalise sfr to sfr(z=0)
        self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +3] =  self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +1] / self.histo_data_ZvsSFR[0, data_offset*(self.a) +1]
             
        print 'z:', format(self.histo_data_ZvsSFR[self.i, data_offset*(self.a)], '0.3f')
        print 'sum sfr / volume(z=0) [Msun yr-1 Mpc-3]', format(self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +1], '0.2e')
        print 'sum sfr [Msun yr-1]', format(self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +2], '0.2e')
        print 'sum sfr/sfr(z=0):', format(self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +3], '0.2e')
        print 'volume [Mpc3]:', format(self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +4], '0.2e')
        
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
            method=''):
                    
                  
        print 'SFR vs Redshift(): --> Selection:', SFH_key,  'method:', method, 'self.i:', self.i, 'self.a:', self.a
        print '++++++++++++++++++++++++++++++++++++++++++++++++++'   

        def sample_selection(data,
                             selection_key=None,
                             sample_key=None,
                             prop='sfr',
                             filename_indices=''): 


            
            print 'SELECTION',

            if sample_key=='passive2' or sample_key=='active2':
                Hubble_z=cd.hubble_z(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'], **fidcosmo)
                ssfr_cut=0.3*Hubble_z*31536000.0                                
                
                print 'ssfr_cut:', format(ssfr_cut, '.3e'), 'max:', format(max(data['ssfr']), '.3e'), 'min:', format(min(data['ssfr']), '.3e')
                if sample_key=='passive2':
                    try:
                        print 'select passive galaxies ssfr<=0.3*Hubble_z',
                        data = data[np.where(data['ssfr']<=ssfr_cut)[:][0]]                      
                    except:
                        print 'NO PASSIVE SELECTION POSSIBLE!'
                  
                elif sample_key=='active2':
                    print 'select active galaxies ssfr>0.3*Hubble_z'
                    try: 
                        data = data[np.where(data['ssfr']>ssfr_cut)[0][:]]        
                    except:
                        print 'NO ACTIVE SELECTION POSSIBLE!'
                        
                print 'ngalaxies selected: ', data.size

            elif sample_key=='passive' or sample_key=='active':
                print 'select galaxies by ssfr!',
                
                frac_to_select=int(np.floor(data.size/100.0*20.0)) 
                data = data[np.argsort(data['ssfr'])]
  
                if sample_key=='active':
                    print '--> select 20% most active galaxies >>'
                    start_size=data.size-frac_to_select                       
                    data=data[start_size:data.size]

                elif sample_key=='passive':
                    print '--> select 20% most passive galaxies <<'
                    data=data[0:int(frac_to_select)]

            elif sample_key=='red' or sample_key=='blue':
                print 'select galaxies by color!'
                
                frac_to_select=int(np.floor(data.size/100.0*20.0)) 
                data = data[np.argsort(data['color'])]
  
                if sample_key=='red':
                    print '--> select 20% reddest (g-i)>>'
                    start_size=data.size-frac_to_select                       
                    data=data[start_size:data.size]

                elif sample_key=='blue':
                    print '--> select 20% bluest (g-i)<<'
                    data=data[0:int(frac_to_select)]
               
    
            #sample selection
            elif sample_key=='low' or sample_key=='high':
                frac_to_select=int(np.floor(data.size/100.0*20.0)) 
                data = data[np.argsort(data[prop])]
                                   
                if sample_key=='high':
                    print '--> select 20% highest ', prop
                    start_size=data.size-frac_to_select
                    data=data[start_size:data.size]

                elif sample_key=='low':
                    print '--> select 20% lowest ', prop
                    data=data[0:int(frac_to_select)]
                
            elif sample_key=='SFR3gt':
                test = str(format(redshift, '0.2f'))
                if test==ext_redshift:
                    print sample_key, '--> SFR>=2 [Msun yr-1]',
                    data = data[np.where(data['sfr']>=2.0)[0][:]]
                    print 'size:', data.size
                else:
                    print 'DO NOTHING!'                    
     
            elif sample_key=='SFR3st':
                test = str(format(redshift, '0.2f'))
                if test==ext_redshift:
                #if test=='1.48':
                    print sample_key, 'select --> SFR<2 [Msun yr-1]',
                    data = data[np.where(data['sfr']<2.0)[0][:]]
                    print 'size:', data.size
                else:
                    print 'DO NOTHING!'

            elif sample_key=='mstar11gt':
                test = str(format(redshift, '0.2f'))
                if test==ext_redshift:
                    print sample_key, 'select --> log10 (Mstar [Msun]) >= 11',
                    data = data[np.where(data['mstar']>=10**11)[0][:]]
                    print 'size:', data.size
                else:
                    print 'DO NOTHING!'                    
     
            elif sample_key=='mstar11st':
                #old name was 'mstar10st'
                test = str(format(redshift, '0.2f'))
                if test==ext_redshift:
                    print sample_key, 'select --> log10 (Mstar [Msun]) < 11',
                    data = data[np.where(data['mstar']<10**11)[0][:]]
                    print 'size:', data.size
                else:
                    print 'DO NOTHING!'

            elif sample_key=='mhalo12gt':
                test = str(format(redshift, '0.2f'))
                if test==ext_redshift:
                    print sample_key, 'select --> log10 (Mvir [Msun]) >= 12.3',
                    data = data[np.where(data['mhalo']>=10**12.3)[0][:]]
                    print 'size:', data.size
                else:
                    print 'DO NOTHING!'                    
     
            elif sample_key=='mhalo12st':
                test = str(format(redshift, '0.2f'))
                if test==ext_redshift:
                    print sample_key, 'select --> log10 (Mvir [Msun]) < 12.3',
                    data = data[np.where(data['mhalo']<10**12.3)[0][:]]
                    print 'size:', data.size
                else:
                    print 'DO NOTHING!'

            elif sample_key=='high-zcold':
                frac_to_select=int(np.floor(data.size/100.0*20.0)) 
                data = data[np.argsort(data['zcold'])]               

                print '--> select 20% highest zcold metallicities >>'
                start_size=data.size-frac_to_select                       
                data=data[start_size:data.size]


            elif sample_key=='low-zcold':
                frac_to_select=int(np.floor(data.size/100.0*20.0)) 
                data = data[np.argsort(data['zcold'])]

                print '--> select 20% lowest zcold metallicities <<'
                data=data[0:int(frac_to_select)]


            elif sample_key=='zcold9gt':
                test = str(format(redshift, '0.2f'))
                if test==ext_redshift:
                    print sample_key, 'select --> zcold >= 9.6 (Gal-dens, Gal2-dens) >=9.1 (Gal400-dens)',
                    data = data[np.where(data['zcold']>=9.0)[0][:]]
                    print 'size:', data.size
                else:
                    print 'DO NOTHING!'                    
     
            elif sample_key=='zcold9st':
                test = str(format(redshift, '0.2f'))
                
                if test==ext_redshift:
                    print sample_key, 'select --> zcold < 9.6 (Gal-dens, Gal2-dens) <9.1 (Gal400-dens)',
                    data = data[np.where(data['zcold']<9.0)[0][:]]
                    print 'size:', data.size
                else:
                    print 'DO NOTHING!'

            elif sample_key=='redmstar11':
                test = str(format(redshift, '0.2f'))
                print sample_key, '1st step',
                frac_to_select=int(np.floor(data.size/100.0*20.0)) 
                data = data[np.argsort(data['color'])]

                print '--> select 20% reddest (g-i)>>',
                start_size=data.size-frac_to_select                       
                data=data[start_size:data.size]
                print 'size:', data.size
                
                if test==ext_redshift:
                    print  '2nd step --> red sample  10.6 >= log10 (Mstar [Msun]) < 11.3',
                    data = data[np.where(data['mstar']>=10**10.6)[0][:]]
                    data = data[np.where(data['mstar']<10**11.3)[0][:]]                    
                    print 'size:', data.size
                else:
                    print 'DO NOTHING!'                    
     

            elif sample_key=='redmhalo13':
                test = str(format(redshift, '0.2f'))
                print sample_key, '1st step',
                frac_to_select=int(np.floor(data.size/100.0*20.0)) 
                data = data[np.argsort(data['color'])]

                print '--> select 20% reddest (g-i)>>',
                start_size=data.size-frac_to_select                       
                data=data[start_size:data.size]
                print 'size:', data.size
                
                if test==ext_redshift:
                    print  '2nd step --> red sample  12.3 >= log10 (Mvir [Msun]) < 13.6',
                    data = data[np.where(data['mhalo']>=10**12.3)[0][:]]
                    data = data[np.where(data['mhalo']<10**13.6)[0][:]]                    
                    print 'size:', data.size
                else:
                    print 'DO NOTHING!'                      

            elif sample_key=='lowmstar11':
                test = str(format(redshift, '0.2f'))
                print sample_key, '1st step',
                frac_to_select=int(np.floor(data.size/100.0*20.0)) 
                data = data[np.argsort(data['sfr'])]

                print '--> select 20% lowest sfr',
                start_size=data.size-frac_to_select                       
                data=data[0:int(frac_to_select)]
                print 'size:', data.size
                
                if test==ext_redshift:
                    print  '2nd step --> low sfr sample  10.6 >= log10 (Mstar [Msun]) < 11.3',
                    data = data[np.where(data['mstar']>=10**10.6)[0][:]]
                    data = data[np.where(data['mstar']<10**11.3)[0][:]]                    
                    print 'size:', data.size
                else:
                    print 'DO NOTHING!'                    
     

            elif sample_key=='lowmhalo13':
                test = str(format(redshift, '0.2f'))
                print sample_key, '1st step',
                frac_to_select=int(np.floor(data.size/100.0*20.0)) 
                data = data[np.argsort(data['sfr'])]

                print '--> select 20% lowest sfr',
                start_size=data.size-frac_to_select                       
                data=data[0:int(frac_to_select)]
                print 'size:', data.size
                
                if test==ext_redshift:
                    print  '2nd step --> low sfr sample  12.3 >= log10 (Mvir [Msun]) < 13.6',
                    data = data[np.where(data['mhalo']>=10**12.3)[0][:]]
                    data = data[np.where(data['mhalo']<10**13.6)[0][:]]                    
                    print 'size:', data.size
                else:
                    print 'DO NOTHING!'  
     
            elif sample_key=='mhalo12st':
                test = str(format(redshift, '0.2f'))
                if test==ext_redshift:
                    print sample_key, 'select --> log10 (Mvir [Msun]) < 12.3',
                    data = data[np.where(data['mhalo']<10**12.3)[0][:]]
                    print 'size:', data.size
                else:
                    print 'DO NOTHING!'

            elif sample_key=='lowZcold-lowMstar':
                test = str(format(redshift, '0.2f'))
                if test==ext_redshift:
                    print sample_key, 'here--> low Zcold / low Mstar',
                    data = data[np.where(data['orphan']==0)[0][:]]
                    data = data[np.where(data['zcold']<=8.5)[0][:]]
                    data = data[np.where(data['mstar']<=4e10)[0][:]]
                    print 'size:', data.size
                else:
                    print 'DO NOTHING!'

            elif sample_key=='lowZcold-highMstar':
                test = str(format(redshift, '0.2f'))
                if test==ext_redshift:
                    print sample_key, 'here--> low Zcold / high Mstar',
                    data = data[np.where(data['orphan']==0)[0][:]]                    
                    data = data[np.where(data['zcold']<=8.5)[0][:]]
                    data = data[np.where(data['mstar']>4e10)[0][:]]
                    print 'size:', data.size
                else:
                    print 'DO NOTHING!'
                                    

            elif sample_key=='highZcold-lowMstar':
                test = str(format(redshift, '0.2f'))
                if test==ext_redshift:
                    print sample_key, 'here--> high Zcold / low Mstar',
                    data = data[np.where(data['orphan']==0)[0][:]]
                    data = data[np.where(data['zcold']>8.5)[0][:]]
                    data = data[np.where(data['mstar']>1.95e9)[0][:]]
                    data = data[np.where(data['mstar']<=1e11)[0][:]]
                    print 'size:', data.size
                else:
                    print 'DO NOTHING!'

            elif sample_key=='highZcold-highMstar':
                test = str(format(redshift, '0.2f'))
                if test==ext_redshift:
                    print sample_key, 'here--> high Zcold / high Mstar',
                    data = data[np.where(data['orphan']==0)[0][:]]                    
                    data = data[np.where(data['zcold']>8.5)[0][:]]
                    data = data[np.where(data['mstar']>1e11)[0][:]]
                    print 'size:', data.size
                else:
                    print 'DO NOTHING!'                    
                    

            elif sample_key=='main' or sample_key=='massive':
                
#                myOutput.writeIntoFile(mycomp+'/anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index_'+str(format(redshift, '0.2f'))+'_mstar0'+sample_key+'.txt',
#                           data[['haloid','parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar', 'SHMR']],
#                           myheader= 'z='+str(format(redshift, '0.2f'))+'\n(0) haloid  (2) parentIndex (1) hostid (3) nodeIndex (4) satelliteNodeIndex (5) orphan (6) Mvir [Msun] (7) Mstar [Msun] (8) SHMR [-]',
#                           data_format="%i\t%i\t%i\t%i\t%i\t%i\t%0.5e\t%0.5e\t%0.1f",
#                           mydelimiter='\t')     
                
                
                print 'filename_indices:', filename_indices
                indices = myData.readAnyFormat(config=False, mypath=filename_indices, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.uint64, skiprow=2) 
                #indices[:,0]-> parentIndex one snapshot before
                #indices[:,1]-> nodeIndex one snapshot before
                
                #expand dimension of the numpy-array if it has only 1 row!
                if np.sum(indices.size)==2:
                    indices = np.expand_dims(indices, axis=0)

                print 'nr indices:', indices[:,0].size
                data = data[np.where(data['orphan']==0)[:][0]] 
                #general selection
                #-------------------------
                #we intersect the current parentIndices with the nodeIndices a snapshot before, a unique array of parentIndices
                #which are present in both arrays are returned. Index1 gives the first occurance of the common IDs in parentIndex,
                #index2 of the same in the nodeIndex a snapshot before
                parentIndices, index1, index2 = np.intersect1d(data['parentIndex'], indices[:,1], return_indices=True)             

                #we test wheater the nodeIndex of this snapshot, in particular, those we already found in the first step 
                #(stored in index1) are linked to the a parentIndex. Thereby we check if those nodeIndices con be found in the
                #parentIndex at any position.          
                test = np.in1d(data['parentIndex'], data['nodeIndex'][index1])
                #if that is true we include those indices in 'index1' too. In other words we have found all satellites to a certain halo                    
                index1 = np.hstack((index1,np.where(test==True)[:][0]))

                #if we find galaxies with the same parentIndex, this means that there are two or more progenitors and we need to find
                #a correct 'main' progenitor.
                #-------------------------
                #we do not need to look through all of them, but only those who have been parents a snapshot before --> that is stored
                #parentIndex (indices[:,0]) a snapshot before. We test of parentIndices a snapshot before are still in the parentIndex
                #of this snapshot
                test_find_non_unique = np.in1d(data['parentIndex'], indices[:,0])
                
                #We indentfy the unique parentIndices among them and store the index of their first occurence (index) of the unique
                #parentIndices and how often they appear (count)
                non_unique, index, count = np.unique(data['parentIndex'][test_find_non_unique], return_index=True, return_counts=True)
                
                #we intersect the parentIndices which are non-unique we found at this snapshot (non-unique) and look for their indicies
                #in the original array data. We only compare those which occure more than one time (coun>1), because if the parentIndex
                #is unique it is automatically chosen since it is the only progenitor.
                #nodeIndices_nu, index_nu1, index_nu2 = np.intersect1d(data['parentIndex'], non_unique[np.where(count>1)[:][0]], return_indices=True)             
                nodeIndices_nu = np.intersect1d(data['parentIndex'], non_unique[np.where(count>1)[:][0]]) 
                #For the others (count>1) we have to choose a main progenitor. So we test where the non-unique nodeIndices are located in the
                #parentIndices of the total array of parentIndecies. 'index_nu' identifies the positions of the non-unique parentIndices
                #-------------------------              
                test_non_unique_indices = np.in1d(data['parentIndex'], nodeIndices_nu)                          
                index_nu = np.where(test_non_unique_indices==True)[:][0]
                #We test if the index of the non-unique (index_nu) is already represented in our index1 list (the list of all progenitors
                #we identified so far)
                test = np.in1d(index_nu, index1)

                #We include those which are not already included
                index1 = np.hstack((index1,index_nu[np.where(test==False)[:][0]]))
                #now we have all progenitors indices linked to the snapshot before and we can withdraw the rest of the IDs
                data=data[index1]
                
             
             
                #We reverse the ordering of the indices starting with the highest and then then first according to the haloid,
                #then parentIndex, and the hostid. This order the indices putting allways the central halo in fron
                
                # (0) haloid  (2) parentIndex (1) hostid (3) nodeIndex (4) satelliteNodeIndex (5) orphan (6) Mvir [Msun] (7) Mstar [Msun] (8) SHMR [-]
                #9416439760	9583218894	9416439760	9416439760	9416439760	0	4.41088e+12	6.51122e+10	67.7
                #9416439760	9416439760	9416439761	6987007766	9416439761	1	4.88697e+10	3.29770e+10	1.5
                #9416439760	9416439760	8394714535	6272408496	8394714535	2	2.44356e+10	4.07484e+09	6.0

                data[::-1].sort(order=['haloid','parentIndex','hostid'], axis=0)
              

#                myOutput.writeIntoFile(mycomp+'/anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index_'+str(format(redshift, '0.2f'))+'_step1_'+sample_key+'.txt',
#                           data[['haloid','parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar', 'SHMR']],
#                           myheader= 'z='+str(format(redshift, '0.2f'))+'--> step1: selecting all haloids linked to a progenitors and sorting remaining halos firstly to haloid, then parentIndex, and final hostid descending\n(0) haloid  (2) parentIndex (1) hostid (3) nodeIndex (4) satelliteNodeIndex (5) orphan (6) Mvir [Msun] (7) Mstar [Msun] (8) SHMR [-]',
#                           data_format="%i\t%i\t%i\t%i\t%i\t%i\t%0.5e\t%0.5e\t%0.1f",
#                           mydelimiter='\t')             
                
                #now we select all unique haloids again, which will atomatically select the central halos!
                unique_parent_id, index1, count_parentIndex = np.unique(data['haloid'], return_index=True, return_counts=True)
                data=data[index1][::-1]
 
#                myOutput.writeIntoFile(mycomp+'/anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index_'+str(format(redshift, '0.2f'))+'_step2_'+sample_key+'.txt',
#                           data[['haloid','parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar', 'SHMR']],
#                           myheader= 'z='+str(format(redshift, '0.2f'))+'--> step2: extracting unique haloids\n(0) haloid  (2) parentIndex (1) hostid (3) nodeIndex (4) satelliteNodeIndex (5) orphan (6) Mvir [Msun] (7) Mstar [Msun] (8) SHMR [-]',
#                           data_format="%i\t%i\t%i\t%i\t%i\t%i\t%0.5e\t%0.5e\t%0.1f",
#                           mydelimiter='\t')                  
#
#                if str(format(redshift, '0.2f'))=='0.59':
#                    data=data[np.where(data['nodeIndex']==8567381780)[:][0]]
                
                if sample_key=='main':
                    print 'SELECT MAIN PROGENITOR!'                               
                    #the most massive is always the one with the lowest ID!
                    
                    data.sort(order=['parentIndex','satelliteNodeIndex'], axis=0)
              
#                    myOutput.writeIntoFile(mycomp+'/anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index_'+str(format(redshift, '0.2f'))+'_step3_'+sample_key+'.txt',
#                           data[['haloid','parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar', 'SHMR']],
#                           myheader= 'z='+str(format(redshift, '0.2f'))+'--> step3: sorint first parentIndex and then satelliteNodeIndex descending\n(0) haloid  (2) parentIndex (1) hostid (3) nodeIndex (4) satelliteNodeIndex (5) orphan (6) Mvir [Msun] (7) Mstar [Msun] (8) SHMR [-]',
#                           data_format="%i\t%i\t%i\t%i\t%i\t%i\t%0.5e\t%0.5e\t%0.1f",
#                           mydelimiter='\t')                  
                  
                else:
                    print 'SELECT MOST MASSIVE!' 
 
                    data[::-1].sort(order=['parentIndex','mhalo'], axis=0)
              
#                    myOutput.writeIntoFile(mycomp+'/anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index_'+str(format(redshift, '0.2f'))+'_mstar2'+sample_key+'.txt',
#                           data[['haloid','parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar', 'SHMR']],
#                           myheader= 'z='+str(format(redshift, '0.2f'))+'\n(0) haloid  (2) parentIndex (1) hostid (3) nodeIndex (4) satelliteNodeIndex (5) orphan (6) Mvir [Msun] (7) Mstar [Msun] (8) SHMR [-]',
#                           data_format="%i\t%i\t%i\t%i\t%i\t%i\t%0.5e\t%0.5e\t%0.1f",
#                           mydelimiter='\t')                  
                    
                    
                unique_parent_id, index1, count_parentIndex = np.unique(data['parentIndex'], return_index=True, return_counts=True)
                data=data[index1][::-1]
                
                                              
#                myOutput.writeIntoFile(mycomp+'/anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index_'+str(format(redshift, '0.2f'))+'_step4_'+sample_key+'.txt',
#                           data[['haloid','parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar', 'SHMR']],
#                           myheader= 'z='+str(format(redshift, '0.2f'))+'--> step4: final selection of progenitors via unique parentIndex\n(0) haloid  (2) parentIndex (1) hostid (3) nodeIndex (4) satelliteNodeIndex (5) orphan (6) Mvir [Msun] (7) Mstar [Msun] (8) SHMR [-]',
#                           data_format="%i\t%i\t%i\t%i\t%i\t%i\t%0.5e\t%0.5e\t%0.1f",
#                           mydelimiter='\t')              
                     
                write_progenitors(data,
                                  filename_indices)
                
#                write_progenitors(data,
#                                  mycomp+'/anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index_'+str(format(redshift, '0.2f'))+'_final'+'.txt')                        

            elif sample_key=='valid':
                print 'check if hostids are valid!'
                data = data[np.where(data['hostid']>-2)[:][0]]
               
 
            else:
                print '--> NO SAMPLE SELECTION DONE!'
                
            return data


        def write_progenitors(data,
                              filename):

            myOutput.writeIntoFile(filename,
                       data[['parentIndex','nodeIndex']],
                       myheader= 'z='+str(format(redshift, '0.2f'))+'\n(1) parentIndex (2) nodeIndex',
                       data_format="%i\t%i",
                       mydelimiter='\t')
           

        def select_CMASS_z056():
            if self.myconfig_array['catname'+str(self.a)].find('Gala')!=-1:
                #CMASS_IDs = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/data/Galacticus_400Mpc/Galacticus_400Mpc_z_0.55_CMASS_down_sample3_haloid.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.uint64, skiprow=2)                         
                #CMASS_IDs = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/data/Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3_haloid.txt', data_format='ASCII', data_shape='shaped', delim='  ', mydtype=np.uint64, skiprow=2)
                #CMASS_IDs = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/data/Galacticus_1Gpc/Galacticus_1Gpc_z_0.09_tarsel_CUT3_Contreras+13_mcold_haloid.txt', data_format='ASCII', data_shape='shaped', delim='  ', mydtype=np.uint64, skiprow=2)
                CMASS_IDs = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/data/300/1h-1Mpc/Galacticus/log/00_MDPL2_Galacticus_z0.00_subsamples_r_1h-1Mpc_log.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.str_, skiprow=8)
                print CMASS_IDs
                exit()
                #indices = np.in1d(self.myData2Process['hostid'], CMASS_IDs[:,1])

                #self.myData2Process = self.myData2Process[np.where(indices==True)[:][0]]
                print 'shape of data after CMASS_ID selection:', self.myData2Process.shape
                #self.myData2Process = self.myData2Process[np.where(self.myData2Process['orphan']==0)[:][0]]
                #test sample                
                #self.myData2Process=self.myData2Process[0:107]
                #test1
                #self.myData2Process=self.myData2Process[np.where(self.myData2Process['haloid']==9583213112)[:][0]]
                #test2
                #self.myData2Process=self.myData2Process[np.where(self.myData2Process['haloid']==9583213304)[:][0]]
                #test3
                #self.myData2Process=self.myData2Process[np.where(self.myData2Process['haloid']==9583212851)[:][0]]
                #test4
                #self.myData2Process=self.myData2Process[np.where(self.myData2Process['haloid']==9583212725)[:][0]]
                #test5
                #self.myData2Process=self.myData2Process[np.where(self.myData2Process['haloid']==9583212324)[:][0]]                
                #test6
                #self.myData2Process=self.myData2Process[np.where(self.myData2Process['haloid']==9583210676)[:][0]]                 
                #test7
                #self.myData2Process=self.myData2Process[np.where(self.myData2Process['haloid']==9583213190)[:][0]]
                self.myData2Process=self.myData2Process[np.where(self.myData2Process['hostid']==12663401752)[:][0]]
  
                print 'method: ', method, '--> select samples and follow progenitors! SFH_key:', SFH_key
                if method=='M2':
                    data=sample_selection(self.myData2Process, sample_key=SFH_key)
                
                    myOutput.writeIntoFile(mycomp+'/anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index_'+str(format(redshift, '0.2f'))+'_'+SFH_key+'_props_'+method+'.txt',
                               data[['haloid','parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar', 'SHMR', 'mcold', 'Mzgas', 'zcold', 'color', 'sfr', 'ssfr']],
                               myheader= 'z='+str(format(redshift, '0.2f'))+'\n(0) haloid  (2) parentIndex (1) hostid (3) nodeIndex (4) satelliteNodeIndex (5) orphan (6) Mvir [Msun] (7) Mstar [Msun] (8) SHMR [-] (9) mcold [Msun] (10) Mzgas [Msun] (11) zcold [-] (12) g-i [-] (13) sfr [Msun yr-1] (14) ssfr [yr-1]',
                               data_format="%i\t%i\t%i\t%i\t%i\t%i\t%0.5e\t%0.5e\t%0.5f\t%0.8e\t%0.8e\t%0.5f\t%0.5f\t%0.8e\t%0.8e",
                               mydelimiter='\t')
                    
                    return data
                
                else:
                     myOutput.writeIntoFile(mycomp+'/anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index_'+str(format(redshift, '0.2f'))+'_props_'+method+'.txt',
                               self.myData2Process[['haloid','parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar', 'SHMR', 'mcold', 'Mzgas', 'zcold', 'color', 'sfr', 'ssfr']],
                               myheader= 'z='+str(format(redshift, '0.2f'))+'\n(0) haloid  (2) parentIndex (1) hostid (3) nodeIndex (4) satelliteNodeIndex (5) orphan (6) Mvir [Msun] (7) Mstar [Msun] (8) SHMR [-] (9) mcold [Msun] (10) Mzgas [Msun] (11) zcold [-] (12) g-i [-] (13) sfr [Msun yr-1] (14) ssfr [yr-1]',
                               data_format="%i\t%i\t%i\t%i\t%i\t%i\t%0.5e\t%0.5e\t%0.5f\t%0.8e\t%0.8e\t%0.5f\t%0.5f\t%0.8e\t%0.8e",
                               mydelimiter='\t')                   
                     return self.myData2Process               

            else:
                self.myData2Process = self.myData2Process[np.where(self.myData2Process['hostid']>-2)[:][0]]
                
                return self.myData2Process

            


            
            

        def calc_extra_props(): 
            
            try:
                #self.myData2Process = rcfuncs.append_fields([self.myData2Process], ['ssfr','SHMR','zcold', 'color', 'mhalo_200c', 'NFW_con'] ,[self.myData2Process['sfr'],self.myData2Process['mhalo'],self.myData2Process['spinParameter'],self.myData2Process[mag_prefix+'i'],self.myData2Process['mhalo'],self.myData2Process['sfr']], usemask=False)
    
                self.myData2Process = rcfuncs.append_fields([self.myData2Process], ['SHMR','color', 'z', 'nsats'] ,[self.myData2Process['mhalo'],self.myData2Process[mag_prefix+'i'],self.myData2Process[mag_prefix+'i'],self.myData2Process['orphan']], usemask=False)
     
                print 'expand array!',               
                #self.myData2Process['ssfr']=self.myData2Process['sfr']/self.myData2Process['mstar']
                self.myData2Process['SHMR']=self.myData2Process['mhalo']/self.myData2Process['mstar']
                #self.myData2Process['zcold']=8.69+np.log10(self.myData2Process['Mzgas']/(self.myData2Process['mcold']*0.0134))
                try:
                    self.myData2Process['color']=self.myData2Process[mag_prefix+'g']-self.myData2Process[mag_prefix+'i']
                except:
                    self.myData2Process['color']=self.myData2Process[mag_prefix+'g']-self.myData2Process[mag_prefix+'i'] 
                                  
                print '--> DONE!'
            except:
                pass
  
        #print self.myData2Process.shape
        #print np.info(self.myData2Process)
        redshift = self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z']
        ext_redshift = '0.0'
        mag_prefix='mAB_dA_total_'
        
        import numpy.lib.recfunctions as rcfuncs
        from cosmolopy import cparam, cd, cc
        fidcosmo = cparam.Planck(flat=True, extras=True)
        #print 'lum_dist:', cd.luminosity_distance(redshift, **fidcosmo), 'z=', redshift
        #print 'SFH_key:', SFH_key, 'method:', method,
        
        if method=='M2':
            filename_indices=mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index'+SFH_key+'.txt'
        else:
            filename_indices=mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index.txt'
            
        print '\n++++++++++++++++++++++++++++++++++++++++++++++++++\nfilename_indices:', filename_indices, '\n++++++++++++++++++++++++++++++++++++++++++++++++++\n'      
                
        calc_extra_props()
        z_test = str(format(redshift, '0.2f'))

        if method=='FT':

            self.myData2Process['z'] = redshift            
            IDs = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/data/ROCKSTAR/MDPL2_merger_tree_info_300-'+SFH_key+'.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.str_, skiprow=2)

            centralID=IDs[np.where(IDs[:,1]==z_test)[:][0]][:,2][0]
            print 'centralID:', centralID,

            
            myfiltered_data=self.myData2Process[np.where(self.myData2Process['haloid']==int(centralID))[:][0]]
            print 'HERE! size: ', myfiltered_data.size
            if myfiltered_data.size>0:               
                myfiltered_data['nsats']=myfiltered_data[np.where(myfiltered_data['orphan']!=0)[:][0]].size
                print 'nsats', myfiltered_data['nsats'][0]
                if myfiltered_data.size>1:
                    print 'more than one central!'
                    myfiltered_data=self.myData2Process[np.where(self.myData2Process['hostid']==int(centralID))[:][0]]
                
                print 'final galaxy found:\n', myfiltered_data
                if z_test==0.0 or z_test=='0.00' or z_test=='0.0':
                    myOutput.writeIntoFile(mycomp+'/anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index_300_test_'+method+'_'+SFH_key+'.txt',
                               myfiltered_data[['z', 'haloid','parentIndex', 'hostid', 'orphan', 'mhalo', 'mstar', 'SHMR', 'mcold', 'Mzgas', 'zcold', 'color', 'sfr', 'ssfr', 'x_pos', 'y_pos', 'z_pos', 'nsats']],
                               myheader= self.myconfig_array['catname'+str(self.a)]+'300 cluster region: '+SFH_key+'z='+str(format(redshift, '0.2f'))+'\n(1) z (2) haloid  (3) parentIndex (4) hostid (5) orphan (6) Mvir [Msun] (7) Mstar [Msun] (8) SHMR [-] (9) mcold [Msun] (10) Mzgas [Msun] (11) zcold [-] (12) g-i [-] (13) sfr [Msun yr-1] (14) ssfr [yr-1] (15) x [comvMpc] (16) y [comvMpc] (17) z [comvMpc] (18) nsats',
                               data_format="%0.2f\t%i\t%i\t%i\t%i\t%0.5e\t%0.5e\t%0.5f\t%0.8e\t%0.8e\t%0.5f\t%0.5f\t%0.8e\t%0.8e\t%0.6f\t%0.6f\t%0.6f\t%i",
                               mydelimiter='\t')
                else:                    
                    myOutput.writeIntoFile(mycomp+'/anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index_300_test_'+method+'_'+SFH_key+'.txt',
                               myfiltered_data[['z', 'haloid','parentIndex', 'hostid', 'orphan', 'mhalo', 'mstar', 'SHMR', 'mcold', 'Mzgas', 'zcold', 'color', 'sfr', 'ssfr', 'x_pos', 'y_pos', 'z_pos','nsats']],
                               myheader= 'z='+str(format(redshift, '0.2f'))+'\n(1) z (2) haloid  (3) parentIndex (4) hostid (5) orphan (6) Mvir [Msun] (7) Mstar [Msun] (8) SHMR [-] (9) mcold [Msun] (10) Mzgas [Msun] (11) zcold [-] (12) g-i [-] (13) sfr [Msun yr-1] (14) ssfr [yr-1] (15) x [comvMpc] (16) y [comvMpc] (17) z [comvMpc] (18) nsats',
                               data_format="%0.2f\t%i\t%i\t%i\t%i\t%0.5e\t%0.5e\t%0.5f\t%0.8e\t%0.8e\t%0.5f\t%0.5f\t%0.8e\t%0.8e\t%0.6f\t%0.6f\t%0.6f\t%i",
                               mydelimiter='\t',
                               append_mytext=True)
            else:
                return error_count

       
        elif self.i==0 and SFH_key_count==0 and method!='M2':
            
            select_CMASS_z056()            

            if self.myconfig_array['catname'+str(self.a)].find('Gala')!=-1:
                              
                write_progenitors(self.myData2Process,
                                  filename_indices)
       
        elif method=='M2':
            if self.i>0:
                print 'i>0:'

                if z_test==ext_redshift and (SFH_key.find('st')!=-1 or SFH_key.find('gt')!=-1):
                    myfiltered_data=sample_selection(self.myData2Process, sample_key=SFH_key)
                        
                    write_progenitors(myfiltered_data,
                                      mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index'+SFH_key+'.txt')                        
                        
                else:
                    if self.myconfig_array['catname'+str(self.a)].find('Gala')!=-1:
                        myfiltered_data=sample_selection(self.myData2Process, sample_key='main', filename_indices=filename_indices)
                    else:
                        myfiltered_data=sample_selection(self.myData2Process, sample_key='valid', filename_indices=filename_indices)
                        
                myOutput.writeIntoFile(mycomp+'/anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index_'+str(format(redshift, '0.2f'))+'_'+SFH_key+'_props_'+method+'.txt',
                       myfiltered_data[['haloid','parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar', 'SHMR', 'mcold', 'Mzgas', 'zcold', 'color', 'sfr', 'ssfr']],
                       myheader= 'z='+str(format(redshift, '0.2f'))+'\n(0) haloid  (2) parentIndex (1) hostid (3) nodeIndex (4) satelliteNodeIndex (5) orphan (6) Mvir [Msun] (7) Mstar [Msun] (8) SHMR [-] (9) mcold [Msun] (10) Mzgas [Msun] (11) zcold [-] (12) g-i [-] (13) sfr [Msun yr-1] (14) ssfr [yr-1]',
                       data_format="%i\t%i\t%i\t%i\t%i\t%i\t%0.5e\t%0.5e\t%0.5f\t%0.8e\t%0.8e\t%0.5f\t%0.5f\t%0.8e\t%0.8e",
                       mydelimiter='\t')
                               
            else:
                print 'i=0', self.i
                myfiltered_data = select_CMASS_z056()
                write_progenitors(myfiltered_data,
                                  filename_indices)                
            
        elif self.i>0 and SFH_key_count==0:
            if self.myconfig_array['catname'+str(self.a)].find('Gala')!=-1:
                self.myData2Process=sample_selection(self.myData2Process, sample_key='main', filename_indices=filename_indices)
            else:
                self.myData2Process=sample_selection(self.myData2Process, sample_key='valid', filename_indices=filename_indices)
          

####### SAMPLE SELECTION     ##################################################################
        #print 'after selection:', myfiltered_data.shape,
        if method=='M1':
                     
            myfiltered_data=sample_selection(self.myData2Process, sample_key=SFH_key)                
            print 'ngal after selection:', myfiltered_data.shape,'\n'#, 'check original:', self.myData2Process.shape

            myOutput.writeIntoFile(mycomp+'/anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)]+'/SFH_Index_'+str(format(redshift, '0.2f'))+'_'+SFH_key+'_props_'+method+'.txt',
                       myfiltered_data[['haloid','parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar', 'SHMR', 'mcold', 'Mzgas', 'zcold', 'color', 'sfr', 'ssfr']],
                       myheader= 'z='+str(format(redshift, '0.2f'))+'\n(0) haloid  (2) parentIndex (1) hostid (3) nodeIndex (4) satelliteNodeIndex (5) orphan (6) Mvir [Msun] (7) Mstar [Msun] (8) SHMR [-] (9) mcold [Msun] (10) Mzgas [Msun] (11) zcold [-] (12) g-i [-] (13) sfr [Msun yr-1] (14) ssfr [yr-1]',
                       data_format="%i\t%i\t%i\t%i\t%i\t%i\t%0.5e\t%0.5e\t%0.5f\t%0.8e\t%0.8e\t%0.5f\t%0.5f\t%0.8e\t%0.8e",
                       mydelimiter='\t')          
        
###############################################################################################

#        myfiltered_data['NFW_con'] = mL.calculate_NFW_con_with_fit(myfiltered_data['mhalo'], 
#                                                                       redshift,
#                                                                       cosmology='Planck', 
#                                                                       overdens='vir')
            
#        myfiltered_data['mhalo_200c']=mL.convert_halo_mass(myfiltered_data['mhalo'], 
#                                                            myfiltered_data['NFW_con'],
#                                                            myfiltered_data['orphan'])


#
        try:
            self.histo_data_ZvsSFR = data
    
            percent_low = 31.7
            percent_high = 68.3
            #calculate properties        
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a)]    = redshift
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +1] = np.nansum(myfiltered_data['sfr'])  
        
            #normalise sfr to sfr(z=0)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +2] =  self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +1] / self.histo_data_ZvsSFR[0, data_offset*(self.a) +1]
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +3] =  np.nanmedian(myfiltered_data['sfr']) 
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +4] =  np.nanpercentile(myfiltered_data['sfr'],percent_low)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +5] =  np.nanpercentile(myfiltered_data['sfr'],percent_high)
     
         
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +6] =  np.nansum(myfiltered_data['ssfr'])
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +7] =  self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +6]/self.histo_data_ZvsSFR[0, data_offset*(self.a) +6]
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +8] =  np.nanmedian(myfiltered_data['ssfr']) 
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +9] =  np.nanpercentile(myfiltered_data['ssfr'],percent_low)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +10] =  np.nanpercentile(myfiltered_data['ssfr'],percent_high)        
            
            try:
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +11]= np.nanmedian(myfiltered_data[mag_prefix+'g']-myfiltered_data[mag_prefix+'i'])
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +12]= np.nanpercentile((myfiltered_data[mag_prefix+'g']-myfiltered_data[mag_prefix+'i']),percent_low)
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +13]= np.nanpercentile((myfiltered_data[mag_prefix+'g']-myfiltered_data[mag_prefix+'i']),percent_high)
            except:
                print 'g-i color calculation faild!'
    
            try:
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +29]= np.nanmedian(myfiltered_data[mag_prefix+'r']-myfiltered_data[mag_prefix+'i']) 
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +30]= np.nanpercentile((myfiltered_data[mag_prefix+'r']-myfiltered_data[mag_prefix+'i']),percent_low)
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +31]= np.nanpercentile((myfiltered_data[mag_prefix+'r']-myfiltered_data[mag_prefix+'i']),percent_high)
            except:
                print 'r-i color calculation faild!' 
    
               
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +14]= np.nansum(myfiltered_data['mstar'])
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +15]= self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +14]/self.histo_data_ZvsSFR[0, data_offset*(self.a) +14]
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +16]= np.nanmedian(myfiltered_data['mstar'])
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +17]= np.nanpercentile(myfiltered_data['mstar'],percent_low)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +18]= np.nanpercentile(myfiltered_data['mstar'],percent_high)
    
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +19]= np.nanmedian(myfiltered_data['rdisk'])
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +20]= np.nanpercentile(myfiltered_data['rdisk'],percent_low)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +21]= np.nanmedian(myfiltered_data['rbulge']/myfiltered_data['rdisk'])
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +22]= np.nanpercentile((myfiltered_data['rbulge']/myfiltered_data['rdisk']),percent_low)
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +23]= np.nanpercentile((myfiltered_data['rbulge']/myfiltered_data['rdisk']),percent_high)        
    
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +24] = cd.lookback_time(self.histo_data_ZvsSFR[self.i, data_offset*(self.a)], z0=self.histo_data_ZvsSFR[0, data_offset*(self.a)], **fidcosmo)/3.1536e+16 #seconds to Gyrs
    
            data_cents=myfiltered_data[np.where(myfiltered_data['orphan']==0)[0][:]]    
            data_sats=myfiltered_data[np.where(myfiltered_data['orphan']==1)[0][:]]
            data_os=myfiltered_data[np.where(myfiltered_data['orphan']==2)[0][:]]
    
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +25] = format(data_cents.size/float(myfiltered_data.size), '0.3f')
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +26] = format(data_sats.size/float(myfiltered_data.size), '0.3f')
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +27] = format(data_os.size/float(myfiltered_data.size), '0.3f')
    
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +28] = format(myfiltered_data.size/float(self.volume)/1e-4,'0.3f')
    
            try:        
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +32]= np.nanmedian(myfiltered_data['mbh'])
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +33]= np.nanpercentile(myfiltered_data['mbh'],percent_low)
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +34]= np.nanpercentile(myfiltered_data['mbh'],percent_high)
            except:
                pass
            try:             
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +35]= np.nanmedian(myfiltered_data['rhalf_mass'])
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +36]= np.nanpercentile(myfiltered_data['rhalf_mass'],percent_low)
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +37]= np.nanpercentile(myfiltered_data['rhalf_mass'],percent_high)          
            except:
                pass
            try: 
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +38]= np.nansum(myfiltered_data['mcold'])
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +39]= self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +38]/self.histo_data_ZvsSFR[0, data_offset*(self.a) +38]
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +40]= np.nanmedian(myfiltered_data['mcold'])
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +41]= np.nanpercentile(myfiltered_data['mcold'],percent_low)
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +42]= np.nanpercentile(myfiltered_data['mcold'],percent_high)        
            except:
                pass    
            try: 
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +43]= np.nansum(myfiltered_data['Mzgas'])
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +44]= self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +43]/self.histo_data_ZvsSFR[0, data_offset*(self.a) +43]
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +45]= np.nanmedian(myfiltered_data['Mzgas'])
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +46]= np.nanpercentile(myfiltered_data['Mzgas'],percent_low)
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +47]= np.nanpercentile(myfiltered_data['Mzgas'],percent_high)
            except:
                pass    
            try: 
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +48]= np.nanmedian(myfiltered_data['zcold'])
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +49]= np.nanpercentile(myfiltered_data['zcold'],percent_low)
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +50]= np.nanpercentile(myfiltered_data['zcold'],percent_high)
            except:
                pass
            try: 
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +52]= np.nansum(myfiltered_data['mhalo'])
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +53]= self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +52]/self.histo_data_ZvsSFR[0, data_offset*(self.a) +52]
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +54]= np.nanmedian(myfiltered_data['mhalo'])
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +55]= np.nanpercentile(myfiltered_data['mhalo'],percent_low)
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +56]= np.nanpercentile(myfiltered_data['mhalo'],percent_high)
            except:
                pass
            try: 
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +57]= np.nanmedian(myfiltered_data['mstar']/myfiltered_data['mhalo'])
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +58]= np.nanpercentile(myfiltered_data['mstar']/myfiltered_data['mhalo'],percent_low)
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +59]= np.nanpercentile(myfiltered_data['mstar']/myfiltered_data['mhalo'],percent_high)            
            except:
                pass
    
            try: 
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +60]= np.nanmedian(myfiltered_data['mean_age_stars_disk'])
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +61]= np.nanpercentile(myfiltered_data['mean_age_stars_disk'],percent_low)
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +62]= np.nanpercentile(myfiltered_data['mean_age_stars_disk'],percent_high)            
            except:
                pass
            try: 
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +63]= np.nanmedian(myfiltered_data['mean_age_stars_spheroid'])
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +64]= np.nanpercentile(myfiltered_data['mean_age_stars_spheroid'],percent_low)
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +65]= np.nanpercentile(myfiltered_data['mean_age_stars_spheroid'],percent_high)            
            except:
                pass
            try: 
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +66]= np.nanmedian(myfiltered_data['vmax'])
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +67]= np.nanpercentile(myfiltered_data['vmax'],percent_low)
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +68]= np.nanpercentile(myfiltered_data['vmax'],percent_high)            
            except:
                pass
            try: 
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +69]= np.nanmedian(myfiltered_data['vdisp'])
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +70]= np.nanpercentile(myfiltered_data['vdisp'],percent_low)
                self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +71]= np.nanpercentile(myfiltered_data['vdisp'],percent_high)            
            except:
                pass
    
    
            self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +51]= self.histo_data_ZvsSFR[self.i, data_offset*(self.a) +1]/self.volume 
        except:
            print 'no halos found at this redshift ...!'
#        try:
#            self.calcFastHisto(myfiltered_data['mstar'], filename[0:len(filename)-4]+'_mstar_histo.txt', 'mstar', 'Msun')
#        except:
#            pass
#        try:
#            self.calcFastHisto(myfiltered_data['sfr'], filename[0:len(filename)-4]+'_sfr_histo.txt', 'sfr', 'Msun yr-1')
#        except:
#            pass
#        try:
#            self.calcFastHisto(myfiltered_data['ssfr'], filename[0:len(filename)-4]+'_ssfr_histo.txt', 'ssfr', 'yr-1')
#        except:
#            pass
#        try:
#            self.calcFastHisto(myfiltered_data['rbulge']/myfiltered_data['rdisk'], filename[0:len(filename)-4]+'_rbulgevsrdisk_histo.txt', 'rbulgevsrdisk', '-')        
#        except:
#            pass
#        try:
#            self.calcFastHisto(myfiltered_data['rhalfmass'], filename[0:len(filename)-4]+'_rhalfmass_histo.txt', 'rhalfmass', '-')        
#        except:
#            pass        
#        try:
#            self.calcFastHisto(myfiltered_data['mhalo'], filename[0:len(filename)-4]+'_mhalo_histo.txt', 'mhalo', 'Msun')        
#        except:
#            pass        
#        try:
#            self.calcFastHisto(myfiltered_data[mag_prefix+'g']-myfiltered_data[mag_prefix+'i'], filename[0:len(filename)-4]+'_g-i_histo.txt', 'g-i', '-', binning='lin')        
#        except:
#            pass
#        try:
#            self.calcFastHisto(myfiltered_data['zcold'], filename[0:len(filename)-4]+'_zcold_histo.txt', 'zcold', '-', binning='lin')        
#        except:
#            pass       
#
#        try:     
#            self.calcFastHisto(myfiltered_data[['mhalo','mstar','orphan']], filename[0:len(filename)-4]+'_SHMF_histo.txt', 'mhalo', '-', binup2D=True)        
#        except:
#            pass
#        try:     
#            self.calcFastHisto(myfiltered_data['mbh'], filename[0:len(filename)-4]+'_mbh_histo.txt', 'mbh', 'Msun')        
#        except:
#            pass        
#
#        try:     
#            self.calcFastHisto(myfiltered_data['mean_age_stars_disk'], filename[0:len(filename)-4]+'_mean_age_stars_disk_histo.txt', 'mean_age_stars_disk', 'Gyr', binning='lin')        
#        except:
#            pass
#        
#        try:     
#            self.calcFastHisto(myfiltered_data['mean_age_stars_spheroid'], filename[0:len(filename)-4]+'_mean_age_stars_spheroid_histo.txt', 'mean_age_stars_spheroid', 'Gyr', binning='lin')        
#        except:
#            pass
#        try:     
#            self.calcFastHisto(myfiltered_data['vmax'], filename[0:len(filename)-4]+'_vmax_histo.txt', 'vmax', 'kms-1')        
#        except:
#            pass
#        try:     
#            self.calcFastHisto(myfiltered_data['vdisp'], filename[0:len(filename)-4]+'_vdisp_histo.txt', 'vdisp', 'kms-1')        
#        except:
#            pass         
#        try:     
#            self.calcFastHisto(myfiltered_data['mcold']/myfiltered_data['mstar'], filename[0:len(filename)-4]+'_cgf_histo.txt', 'cgf', '-')        
#        except:
#            pass    
##        
##        try:
##            self.HODFunction(myfiltered_data, filename, which_halomass='mhalo_200c')
##        except:
##            pass
#       
#        #try:
#            #print 'mycond array:', mycond_array
#        self.TwoPCF(myfiltered_data, filename, mycond_array)
##        except:
##            print 'CLUSTERING CALCULACTION NOT POSSIBLE!'

#        if z_test=='0.56':
#
#            
#            NFW_con = mL.calculate_NFW_con_with_fit(myfiltered_data['mhalo'], 
#                                                                       redshift,
#                                                                       cosmology='Planck', 
#                                                                       overdens='vir')
#            
#            myfiltered_data['mhalo']=mL.convert_halo_mass(myfiltered_data['mhalo'], 
#                                                                NFW_con,
#                                                                myfiltered_data['orphan'])
#             
#            self.HODFunction(myfiltered_data, filename, which_halomass='mhalo')
#
#            myOutput.writeIntoFile(filename[0:len(filename)-4]+'_haloids.txt',
#                       myfiltered_data[['haloid','hostid','orphan','mhalo','mstar']],
#                       myheader= 'z='+str(format(redshift, '0.2f'))+'\n(1) haloid (2) hostid (3) orphan (4) mhalo_200c [Msun] (5) mstar [Msun]',
#                       data_format="%i\t%i\t%i\t%0.8e\t%0.8e",
#                       mydelimiter='\t')            
        
        return error_count
    
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
                    print 'correct IMF --> Kroupa to Chabrier'
                    myfiltered_data['mstar']/=10**(0.03925)
                    myfiltered_data['mstar']/=10**(0.2)
    
                    #print 'correct little-h for', self.myconfig_array['catname'+str(self.a)]
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
    #        print 'selection: --> highest', myfiltered_data.shape        
            
    #        Blue-Red separation
            #myfiltered_data=myfiltered_data[np.where(myfiltered_data['mAB_total_g']-myfiltered_data['mAB_total_i']>2.35)[0][:]] 
    #        Guo+13 cut
    #        (r-i)>0.679 -0.082(Mi+20)
#            cut=0.679-0.082*(myfiltered_data['MAB_dA_total_i']-5*np.log10(0.6777)+20.0)
#            myfiltered_data=myfiltered_data[np.where(myfiltered_data['mAB_dA_total_r']-myfiltered_data['mAB_dA_total_i']>cut)[0][:]]         
#            print name_x, 'selection: --> r-i > ', myfiltered_data.shape
             
            #myfiltered_data=myfiltered_data[np.where(myfiltered_data['sfr']/myfiltered_data['mstar']>1e-11)[0][:]]
            #myfiltered_data=myfiltered_data[np.where(myfiltered_data['sfr']>1e-4)[0][:]]         
    #        
            #print name_x, 'selection: --> g-i ', myfiltered_data.shape
            
            return myfiltered_data
    
        def SFRF(plot_key, myfiltered_data):
            return myfiltered_data
                       
            
        def sSFRF(plot_key, myfiltered_data):
            if np.all(myfiltered_data['ssfr']==-99.0):
                myfiltered_data['ssfr']=myfiltered_data['sfr']/myfiltered_data['mstar']
                
            return myfiltered_data 
 

        def ssfr2mstar(plot_key, myfiltered_data):
            
            myfiltered_data['sfr']/=myfiltered_data['mstar']           
            print 'starforming cut:', starforming_cut
            myfiltered_data = myData.selectData2Compute(myfiltered_data, 
                                                        selected_col='sfr', 
                                                        operator='>', 
                                                        condition=starforming_cut)

            return myfiltered_data
 
        def HMF(plot_key, myfiltered_data):
            print ' '
            print 'HMF_no' 
            print '+++++++++++++++++'
            print 'shape before non-orphan selection!', myfiltered_data.shape,
         
            try:
                print '--> select non-orphans',
                myfiltered_data = myData.selectData2Compute(myfiltered_data, 
                                                            selected_col='orphan', 
                                                            operator='<', 
                                                            condition=2)
                                                                                                   
                print 'shape after non-isolated selection!', myfiltered_data.shape                                            
            except:
                print 'no col "orphan" found! Could not do selection ...'
                
            return myfiltered_data
                                                               
        def mstar2mhalo(plot_key, myfiltered_data):
            print 'mstar2mhalo\n++++++++++++++++++'
            #approximated 200c overdensity converstion

            if plot_key.find('vsSFR')!=-1:
                myfiltered_data[name_y]=myfiltered_data['mstar']/myfiltered_data[name_y]
 
            return myfiltered_data                                                              

        def Mzgas2mcold(plot_key, myfiltered_data):
            print 'Mzgas2mstar\n++++++++++++++++++'
            
            print 'Zgas=8.69+log10(Mzgas/(Mcold*0.0134)) ref: Allende-Pieto+01, Asplund+09'
            if self.myconfig_array['catname'+str(self.a)].startswith('Galacticus') and self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_data_format']=='SAMHDF5':
                #column Zgas_spheroid and Zgas_disk are wrongly named in the catalog it should be Mzgas_spheroid and Mzgas_disk
                print 'Galaticus: zgas_disk+zgas_spheroid = Mzgas' 
                myfiltered_data['zgas_spheroid']+=myfiltered_data['zgas_disk']

            myfiltered_data=myfiltered_data[np.where(myfiltered_data[add_axis]>0.0)[0][:]]

            myfiltered_data[name_y]=myfiltered_data[name_y]/(myfiltered_data[add_axis]*0.0134)                
            
            myfiltered_data=myfiltered_data[np.where(myfiltered_data[name_y]>0.0)[0][:]]
                                            
            myfiltered_data[name_y]=8.69+np.log10(myfiltered_data[name_y]) 

            myfiltered_data = myData.selectData2Compute(myfiltered_data, 
                                                        selected_col=name_y, 
                                                        operator='!=', 
                                                        condition='nan')

            
            print myfiltered_data[name_y][0:10]
            
            if name_x.find('mcold')!=-1:

                myfiltered_data=myfiltered_data[np.where(myfiltered_data[name_x]>1e3)[0][:]]                    
                
                
                myfiltered_data = mL.filter_data_before_analysis(myfiltered_data, 
                                                                   5e10, 
                                                                   1e11, 
                                                                   'mstar')                
                print 'mstar selection: -->', myfiltered_data.shape,
                
                myfiltered_data[name_x]/=myfiltered_data['mstar'] 

                myfiltered_data = mL.filter_data_before_analysis(myfiltered_data, 
                                                                   1e-4, 
                                                                   1e4, 
                                                                   name_x)
                print 'Zcold selection: -->', myfiltered_data.shape,                                                    
                myfiltered_data = mL.filter_data_before_analysis(myfiltered_data, 
                                                                   6, 
                                                                   11, 
                                                                   name_y)                                                                   
                print 'Mcold/Mstar selection: -->', myfiltered_data.shape 


            return myfiltered_data                                                 



        def zgas2mstar(plot_key, data):
            print '\nzgas2mstar\n+++++++++++++++++'

            try:
                data['cgf']=data['mcold']/data['mstar']
                
                print 'CGF min/max', format(np.log10(min(data['cgf'])), '0.3f'), '/', format(np.log10(max(data['cgf'])), '0.3f')
            except:
                try:
                    data['cgf']=data['mcold_disk']/data['mstar']
                    
                    print 'SAGE: CGF min/max', format(np.log10(min(data['cgf'])), '0.3f'), '/', format(np.log10(max(data['cgf'])), '0.3f')
                except:
                    pass
            

            try:
                data['ssfr']=data['sfr']/data['mstar'] 
                print 'starforming cut?', starforming_cut,
                
                data=data[np.where(data['ssfr']>starforming_cut)[:][0]]
                                                            
                print 'min', min(data['ssfr']), 'max:', max(data['ssfr']), '-->',data.shape
            except:
                print 'starforming cut FAILED!' 

#            print 'min', min(data['Mzgas_disk']), 'max:', max(data['Mzgas_disk']), '-->',data.shape
#            print 'min', min(data['Mzgas_spheroid']), 'max:', max(data['Mzgas_spheroid']), '-->',data.shape


            try:
                data['Mzgas']=data['Mzgas_disk']+data['Mzgas_spheroid']
                print 'Mzgas=Mzgas_disk+Mzgas_spheroid', '-->',data.shape
            except:
                pass

            data=data[np.where(data['mcold']>0.0)[:][0]]
            print '-->',data.shape
            try:
                data=data[np.where(data[str(add_axis)]>0.1)[:][0]]
            except:
                pass                   

            # Converstion value 0365507 comes from ...
            # (O/H)_gal / (O/H)_sun = Z_gal / Z_sun
            # (O/H)_gal = (O/H)_sun * Z_gal / Z_sun
            # (O/H)_sun = 10**(8.38-12) --> from 12+log10(0/H)_sun = 8.69 (Allende-Prieto+01)
            # Z_sun=0.0134 --> from (Asplund+09)
            # (O/H)_gal = M_metals / M_baryons --> MZ_gas total / Mcold gas                                                      
            # (O/H)_sun / Z_sun = 0.0365507
            
            # myfiltered_data[name_y]=12+np.log10(0.0365507*myfiltered_data[name_y]) == myfiltered_data[name_y]=8.69+np.log10(myfiltered_data[add_axis]/myfiltered_data[name_y]/0.0134)

            try:
                data[name_y]=8.69+np.log10(data['Mzgas']/data['mcold']/0.0134)
            except:
                print 'calculation failed! Try with disk properties only!',
                data[name_y]=8.69+np.log10(data['Mzgas_disk']/data['mcold_disk']/0.0134)
                print '--> DONE!'

                
            data=data[np.where(np.isfinite(data[name_y]))[:][0]]
            print  'after isfinite -->', data.shape
            try:
                data=data[np.where(np.isfinite(data[add_axis]))[:][0]]
            except:
                pass
            
            #myfiltered_data[name_y]=12+np.log10(0.0365507*myfiltered_data['zcold'])
            print 'Zcold min/max', min(data[name_y]), '/', max(data[name_y]), '-->', data.shape
            return data                                               


                    
        def oh2mstar(plot_key, data):
            print '\noh2mstar\n++++++++++++++++++'

            try:
                data['cgf']=data['mcold']/data['mstar']
                
                print 'CGF min/max', format(np.log10(min(data['cgf'])), '0.3f'), '/', format(np.log10(max(data['cgf'])), '0.3f')
            except:
                try:
                    data['cgf']=data['mcold_disk']/data['mstar']
                    
                    print 'SAGE: CGF min/max', format(np.log10(min(data['cgf'])), '0.3f'), '/', format(np.log10(max(data['cgf'])), '0.3f')
                except:
                    pass

            
            try:
                data['ssfr']=data['sfr']/data['mstar'] 
                print 'starforming cut?', starforming_cut,
                
                data=myfiltered_data[np.where(data['ssfr']>starforming_cut)[:][0]]
                                                            
                print 'min', min(data['ssfr']), 'max:', max(data['ssfr']), '-->', data.shape
            except:
                pass
            
            try:
                data=data[np.where(data['zgas_disk']>0.0)[:][0]] 
                data=data[np.where(data['zgas_disk']<=1.0)[:][0]] 
            except:
                data=data[np.where(data['OH_gas_disk_bulge']>0.0)[:][0]] 
                data=data[np.where(data['OH_gas_disk_bulge']<=1.0)[:][0]]
                
            print 'check 0.0 > zgas >= 1.0 (its a metal abundance!) -->',data.shape
            data=data[np.where(data[str(add_axis)]>0.1)[:][0]]
            if self.myconfig_array['catname'+str(self.a)].startswith('SAGE'):
                print 'SAGE!',
                                                             
                data[name_y]=8.69+np.log10(data['zgas_disk']/0.0134)
            else:
                 
                data[name_y]=12.0+np.log10(data['OH_gas_disk_bulge']/16.0)                


            data=data[np.where(np.isfinite(data[name_y]))[:][0]]

            
            print 'oh min/max', min(data[name_y]), '/', max(data[name_y])
            try:
                print 'add axis min/max', add_axis, min(data[add_axis]), '/', max(data[add_axis])
            except:
                pass
            return data

          
        print ' '
        print ' '                 
        print 'binAndFrac2D()', 'self.i:', self.i, 'self.a:', self.a
        print '++++++++++++++++++++++++++++++++++++++++++++++'
        print ' '


        self.ratio_data_binAndFrac2D = data        
        print
        print 'CONFIGURATIONS:'
        print '~~~~~~~~~~~~~~~~'
        print 'PLOTKEY:', plot_key
        print 'binUp2D:', binUp2D
        print 'div_y2x:', div_y2x
        print 'add_axis:', add_axis
        print 'log10bin:', log10bin
        print 'cumulative:', cumulative
        print 'normalise:', normalise
        print ' '
        #self.correctAndConvertUnits()
        
#        try:
#            self.myData2Process['mhalo']=self.myData2Process['mhalo']*0.884
#            print 'doing 200c approx. convertion by multiply mhalo*0.884 ...' 
#        except:
#            pass
        
        print 'name_x:', name_x, 'col id:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_'+name_x], 'mycond x-axis: min/max', mycond_min_x, '/', mycond_max_x

        #print np.info(self.myData2Process)
        myfiltered_data = mL.filter_data_before_analysis(self.myData2Process, 
                                                           mycond_min_x, 
                                                           mycond_max_x, 
                                                           name_x)
        print '\t', name_x, 'selection: -->', myfiltered_data.shape
                                                                                                                                           

        if name_y!=False:
            print 'name_y:', name_y, 'col id:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_'+name_y], 'mycond y-axis: min/max', mycond_min_y, '/', mycond_max_y
            myfiltered_data = mL.filter_data_before_analysis(myfiltered_data, 
                                                               mycond_min_y, 
                                                               mycond_max_y, 
                                                               name_y)                                                       
            print '\t', name_y, 'selection: -->', myfiltered_data.shape

        try: 
            print 'add_axis:', add_axis, 'col id', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_'+add_axis]
            myfiltered_data = mL.filter_data_before_analysis(myfiltered_data, 
                                                               'min', 
                                                               3.16e9, 
                                                               'mstar')                                                       
            print '\t', add_axis, 'selection: -->', myfiltered_data.shape
        except:
            print '\tno axis to add ... '


        if plot_key.find('_no')!=-1 or plot_key.find('vsSFR')!=-1:                          
            print 'exclude orphan galaxies:',
            myfiltered_data=myfiltered_data[np.where(myfiltered_data['orphan']==0)[0][:]] 
            print '\tselection: -->', myfiltered_data.shape
           
        
        data2bin = caseSwitcher_plot_key(plot_key, myfiltered_data)

                    
        print ' '        
        print 'check data:'
        print '+++++++++++++'
                                    
        print 'name_x:', name_x ,'min/max:', min(data2bin[name_x]), '/', max(data2bin[name_x])
        if name_y!=False: 
            print 'name_y:', name_y ,'min/max:', min(data2bin[name_y]), '/', max(data2bin[name_y])
        try:                
            print 'add_axis:', add_axis,'min/max:', min(data2bin[add_axis]), '/', max(data2bin[add_axis])
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
            print 'weight function --> OBERVATIONAL DATA SET [name_x,weight] --> col name weights:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_name_weights']
            histo, binsize = calcHisto(data2bin[[name_x,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_name_weights']]])          
        else:   
                
            try:
                print '[name_x,name_y,add_axis] -->',
                histo, binsize = calcHisto(data2bin[[name_x,name_y,add_axis]])
            except:
                try:
                    print '[name_x,name_y] -->',
                    histo, binsize = calcHisto(data2bin[[name_x,name_y]])
                except:
                    print '[name_x] -->',
                    histo, binsize = calcHisto(data2bin[name_x])

        if float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'])==0.0:
            print '... calculate survey volume. skycoverage:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_skycoverage'], '[deg2]', 'zmax:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_zmax'], 'zmin:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_zmin'], '-->',
            mynorm_y, myunit_volume = mL.calc_survey_volume(float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_skycoverage']),
                                             float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_zmin']),
                                             float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_zmax']),
                                             little_h_out=False)
            print 'volume:', mynorm_y, myunit_volume
            
        elif volume_units.find('h3')!=-1:
            print volume_units, '--> found h3!'
            mynorm_y=(float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size']))**3
        else:
            print volume_units, '--> default!'
            mynorm_y=(float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'])/float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_hubble_par']))**3                            

        if normalise==True:
            print self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'], self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_hubble_par']
            print 'normalise histogram: by', (float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size'])/float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_hubble_par']))**3                                                                                                                   
                    
            self.ratio_data_binAndFrac2D[:,data_offset*self.a:data_offset*self.a+data_offset] = myFuncs.normaliseHisto(histo[:,[0,1,2,3,4,5,6]],
                                                                                                                       norm_y=mynorm_y,
                                                                                                                       x_error=binsize)
                                                                                    
        else:
            self.ratio_data_binAndFrac2D[:,data_offset*self.a:data_offset*self.a+data_offset]=histo                                                                                                      
   
               
    def TwoPCFCUTEBox(self,
                      myfilename):

        print '2PCF with CUTE box:'
        print '--------------------------------'
        myConfig.getCUTEParameterFile()
        
        print 'myfilename:', myfilename[len(mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)])+1::]
              
        myConfig.CUTE_params['output_filename_param']=mycomp+'anaconda/pro/myRun/histos/twoPCF/twoPCF_'+myfilename[len(mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.a)])+1::]
        myConfig.CUTE_params['box_size']='box_size'
        myConfig.CUTE_params['box_size_param']=str(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_CF_box_size'])
        myConfig.CUTE_params['corr_estimator_param']=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_corr_estimator']
        
        data_array = myData.readAnyFormat(config=False, 
                                          mypath=myfilename, 
                                          mydtype=np.float32, 
                                          data_format='ASCII', 
                                          data_shape='shaped',
                                          nr_col=3,
                                          nr_rows=50000,
                                          delim='  ')
                                          
        print 'data_array:', data_array.shape

        myOutput.writeIntoFile(mycomp+'anaconda/pro/myRun/histos/twoPCF/twoPCF_data.dat',
                               data_array,
                               data_format="%0.2f\t%0.2f\t%0.2f")
                               
        myConfig.CUTE_params['data_filename_param']=mycomp+'anaconda/pro/myRun/histos/twoPCF/twoPCF_data.dat'
        
        if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_skip_reading_data']!='yes':
            #Write galaxy position in CUTE readable file 'x y z' without any header!
            myOutput.writeIntoFile(myConfig.CUTE_params['data_filename_param'],
                                   self.myData2Process[:,[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_x_pos'], self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_y_pos'], self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_z_pos']]],
                                   mydelimiter='\t')
            redshift=myData.redshift
        else:
            redshift = self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_manual_input_redshift']
        
        myparams=''
        i=0
        while i<myConfig.CUTE_params['nr_entries']-1:
            #print 'i:', i, myConfig.CUTE_params['name'+str(i)], myConfig.CUTE_params[myConfig.CUTE_params['name'+str(i)]+'_param']
            myparams+=myConfig.CUTE_params['name'+str(i)]
            myparams+='= '
            myparams+=myConfig.CUTE_params[myConfig.CUTE_params['name'+str(i)]+'_param']
            myparams+='\n'                              
            i+=1

        myOutput.writeIntoFile(mycomp+myConfig.CUTE_params['path_to_CUTE_param'][0::]+'params.txt',
                               myparams,
                               data_is_string=True)             
            
        subs.call(mycomp+'anaconda/pro/myRun/invoke_CUTE'+self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_choose_CUTE']+'.sh', shell=True)
        
        data_array = myData.readAnyFormat(config=False, 
                                          mypath=myConfig.CUTE_params['output_filename_param'], 
                                          mydtype=np.float32, 
                                          data_format='ASCII', 
                                          data_shape='shaped',
                                          nr_col=4,
                                          nr_rows=self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_CF_NR_bins'],
                                          skiprow=2,
                                          delim=' ')     

        data = np.zeros((data_array[:,0].size, 6), dtype=np.float32)
        print 'data_array'
        data[:,[0,1,2,5]] = data_array[:,[0,1,2,3]]

        myOutput.writeIntoFile(myConfig.CUTE_params['output_filename_param'],
                               data,
                               myheader=self.myconfig_array['catname'+str(self.a)]+' z='+str(float("{0:.2f}".format(redshift)))+' Rmax: '+str(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_CF_R_max'])+' Mpc boxsize: '+str(myConfig.CUTE_params['box_size_param'])+' Mpc\n(1) r [Mpc]  (2) xi(r)  (3) error(r)  (4) - (5) - (6) DD(r)',
                               data_format="%0.4f\t%0.4f\t%0.4f\t%d\t%d\t%d",
                               mydelimiter='\t') 
        
        exit()

    def TwoPCFMonopol(
                self,
                data,
                nbins,
                non_orphan=False,
                histo_data_min='min',
                histo_data_max='max',
                mycond_min='min',
                mycond_max='max',
                mycond_min_boxsize='min',
                mycond_max_boxsize='max',
                mycond_ngal='max',
                data_offset=0,
                histo_output_offset=0,
                mylog10bin=True):
                    
        print ' '
        print ' '                 
        print '2PCC_Monopol()', 'self.i:', self.i, 'self.a:', self.a
        print '++++++++++++++++++++++++++++++++++++++++++++++'
        print ' '
        print 'catname', self.myconfig_array['catname'+str(self.a)]

        self.data_2PCF = data

#        print 'x min:', np.min(self.myData2Process[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_x_pos']]), 'x max:', np.max(self.myData2Process[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_x_pos']]), '[Mpc]'
#        print 'y min:', np.min(self.myData2Process[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_y_pos']]), 'y max:', np.max(self.myData2Process[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_y_pos']]), '[Mpc]' 
#        print 'z min:', np.min(self.myData2Process[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_z_pos']]), 'z max:', np.max(self.myData2Process[:,self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_z_pos']]), '[Mpc]'
#        
#                              
#        data_2PCF = mL.filter_data_before_analysis(self.myData2Process, 
#                                               mycond_min, 
#                                               mycond_max, 
#                                               self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_mstar'])                                                       
#
#        data_2PCF = mL.filter_data_before_analysis(data_2PCF, 
#                                                mycond_min_boxsize, 
#                                                mycond_max_boxsize, 
#                                                self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_x_pos'])
#
#        data_2PCF = mL.filter_data_before_analysis(data_2PCF, 
#                                                mycond_min_boxsize, 
#                                                mycond_max_boxsize, 
#                                                self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_y_pos'])
#
#        data_2PCF = mL.filter_data_before_analysis(data_2PCF, 
#                                                mycond_min_boxsize, 
#                                                mycond_max_boxsize, 
#                                                self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_z_pos'])


        data_2PCF = myData.readAnyFormat(config=False, 
                                    mypath=mycomp+'anaconda/pro/myRun/histos/twoPCF/data_twoPCF.txt', 
                                    data_format='ASCII', 
                                    data_shape='shaped', 
                                    delim='\t', 
                                    mydtype=np.float64, 
                                    skiprow=2)                      


        self.data_2PCF[:,data_offset*self.a:data_offset*self.a+2] = myFuncs.calc_twoPCF_monopol(data_2PCF[:, [self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_x_pos'], self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_y_pos'], self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_col_z_pos']] ], 
                                        nbins, 
                                        log10bin=mylog10bin,
                                        histo_min=histo_data_min,
                                        histo_max=histo_data_max,
                                        sub_box_size=float(mycond_max_boxsize)-float(mycond_min_boxsize),
                                        ngalaxies=mycond_ngal,
                                        use_random_cat=True)                                                                              

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



        #Auto-correlation on periodic, ξ(r) for a COSMOLOGICAL BOX
        def xi():       return cf.theory.xi(boxsize, nthreads, rbins, data['x_pos'], data['y_pos'], data['z_pos'], output_ravg=True, c_api_timer=True) 

        #Projected auto-correlation function, wp(rp) for a COSMOLOGICAL BOX
        def wp():       return cf.theory.wp(boxsize, pimax, nthreads, rbins, data['x_pos'], data['y_pos'], data['z_pos'], output_rpavg=True, c_api_timer=True) 

        #Clustering in 3-D Pair counts for (auto or cross) correlations for ξ(r)      
        def xi_of_r():  return cf.theory.xi_of_r(boxsize, pimax, nthreads, rbins, data['x_pos'],data['y_pos'], data['z_pos'], output_rpavg=True, c_api_timer=True) 
        
        #Pair counts (auto or cross) correlations for ξ(rp,π)
        def xi_rp_pi(): return cf.theory.xi_rp_pi(boxsize, pimax, nthreads, rbins, data['x_pos'],data['y_pos'], data['z_pos'], output_rpavg=True, c_api_timer=True) 
           
            

        # Auto pairs counts in data 
        def calc_DD(autocorr=1): return cf.mocks.DDrppi_mocks(autocorr, cosmology, nthreads, pimax, rbins, 
                                                              data['RA'], data['DEC'], data['Z'], 
                                                              is_comoving_dist=False)

        # Cross pair counts in data-random
        def calc_DR(autocorr=0): return cf.mocks.DDrppi_mocks(autocorr, cosmology, nthreads, pimax, rbins,
                                                              data['RA'], data['DEC'], data['Z'],
                                                              RA2=rand_RA, DEC2=rand_DEC, CZ2=rand_CZ,
                                                              is_comoving_dist=False)

        # Auto pairs counts in random                                      
        def calc_RR(autocorr=1): return cf.mocks.DDrppi_mocks(autocorr, cosmology, nthreads, pimax, rbins,
                                                              rand_RA, rand_DEC, rand_CZ)

        # All the pair counts are done, get the projected correlation function  
        def wp_mocks(): return convert_rp_pi_counts_to_wp(N, N, rand_N, rand_N,
                                                          DD_counts, DR_counts, DR_counts, RR_counts, 
                                                          nbins, pimax, c_api_timer=True) 
          


        def writeIntoFile(results_CF,
                          name):
     
            myheader=self.myconfig_array['catname'+str(self.a)]+' z='+str(float("{0:.2f}".format(redshift)))+' r_max: '+str(rmax)+' ['+pos_unit+'] '+'r_min: '+str(rmin)+' ['+pos_unit+'], boxsize: '+str(boxsize)+' ['+pos_unit+'], '
            if name.startswith('wp'): 
                myheader+='ngal: '+str(len(data))+', pi_max: '+str(pimax)+'\n'
            else:
                myheader+='ngal: '+str(len(data))+'\n'

            if name!='wp_mocks':
                data_format="%0.8e\t"*5+'%0.8e'                   
                myheader+='(1) rmin ['+pos_unit+']  (2) rmax ['+pos_unit+']  (3) r_avg ['+pos_unit+'] (4) '+str(name)+' (5) npairs [-] (6) weight [-]'   
            else:
                data_format="%0.8e\t%0.8e"
                myheader+='(1) r_mean ['+pos_unit+']  (2) wp ['+pos_unit+']'                    
                results_CF = mL.dim_expander(results_CF,2)

            
                i=0
                while i<results_CF[:,0].size:
                    results_CF[i,0]=(rbins[i]*rbins[i+1])**0.5
                    i+=1

            

            myheader+=' galaxy types cents: '+str(centrals)+' no-sats: '+str(sats)+' orphans: '+str(orphans)

            if cut_list_operators_dict[operator]=='_mstar_' or cut_list_operators_dict[operator]=='_kroup_mstar_':
                cut_label=cut_list_operators_dict[operator]+str(cut)+'_'+str(cut_list_max_dict[cut])
            elif cut_list_operators_dict[operator].find('mstar')!=-1:
                cut_label=cut_list_operators_dict[operator]+str(float("{0:.2f}".format(np.log10(cut))))                  
            elif cut_list_operators_dict[operator].find('Mr')!=-1:
                cut_label=cut_list_operators_dict[operator]
            elif cut!='':
                cut_label=cut_list_operators_dict[operator]+str(cut)
            else:
                cut_label=''

            #print cut_label
            #print data_format
            myOutput.writeIntoFile(#filename[0:len(filename)-4]+'_'+name+'.txt',
                                   filename[0:len(filename)-4]+'_'+gtype+cut_label+'_'+str(rmin)+'_'+str(rmax)+'_'+str(pimax)+'_'+name+'.txt',
                                   results_CF,
                                   myheader=myheader,
                                   data_format=data_format)
            


        print ' '
        print ' '                 
        print 'TwoPCF Corrfunc()', 'self.i:', self.i, 'self.a:', self.a
        print '++++++++++++++++++++++++++++++++++++++++++++++'
        print ' '
        #print 'catname', self.myconfig_array['catname'+str(self.a)], 'which:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_which'], 'z:', self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z']  
        

        import Corrfunc as cf
        #from cf.theory import wp 
        from os.path import dirname, abspath, join as pjoin
        from Corrfunc.io import read_catalog
        from Corrfunc.utils import convert_rp_pi_counts_to_wp     
        
        redshift  = float(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'])

        cosmology=2

        if myconds_array['x_pos_unit'].find('h-1')!=-1:
            pos_unit='h-1Mpc'
        else:
            pos_unit='Mpc'            
      
        gtype_list=['no', 'centrals']#, 'all']

        gtype_dict=          {'centrals': 0,    'no': 2, 'no-sats': 1, 'all': 10}       
        gtype_dict_operators={'centrals': '==', 'no': '<', 'no-sats': '==', 'all': '<'}

        #self.myData2Process= mL.choose_random_sample(self.myData2Process, 95683)
        
        #self.myData2Process['mstar']*=1.16 
        #print 'min/max mstar:', format(min(self.myData2Process['mstar']), '.3e'), '/', format(max(self.myData2Process['mstar']), '.3e')
 
        if mydata==False:
            mydata=self.myData2Process

#        mydata = mL.choose_random_sample(self.myData2Process,int(self.myData2Process.size*0.1))
#        print 'randomly chosen 10 %', mydata.size
        #print 'here 3858', mydata   
       
        for gtype in gtype_list:
            
            print 'gtype:', gtype, 'operator:',  gtype_dict_operators[gtype], 'condition:', gtype_dict[gtype]
            cut_list=['']
            cut_list_max_dict={}
            #> or <
            #cut_list=[1e10]

            #cut_list=[1.7783e11, 2.2909e11, 3.02e11]
            #cut_list_max_dict={1.7783e11: 2.2909e11, 2.2909e11: 3.02e11, 3.02e11:1e15}
            #mass bins
            #Msun
            #cut_list=[10**11.25, 10**11.35, 10**11.45, 10**11.55, 10**11.65]
            #Msun
            #cut_list=[10**11.08, 10**11.18, 10**11.28, 10**11.38, 10**11.48]
            #cut_list=[1.7783e11, 2.2387e11, 2.8184e11, 3.5481e11, 4.4668e11]
            #cut_list_max_dict={1.5e11: 2e11, 1.73e11: 3e11, 2e11: 3e11, 2.24e11: 4e11, 3e11: 4e11, 3.5e11: 5e11}
            #cut_list=[1.5e11, 1.73e11, 2e11, 2.24e11, 3e11, 3.5e11, 4e11, 1e15]

            #bin1
            #cut_list=[1.5e11, 1.73e11, 2.24e11, 3e11, 1.05e11, 2e11]
            #cut_list_max_dict={1.5e11: 2e11, 1.73e11: 3e11, 2.24e11: 4e11, 3e11: 4e11, 1.05e11: 3.45e11, 2e11: 1e13}

            #cut_list=[1.5e11, 1.73e11, 1.05e11]
            #cut_list_max_dict={1.5e11: 2e11, 1.73e11: 3e11, 1.05e11: 3.45e11}            
          
           
            cut_list_operators=['>']
            cut_list_operators_dict={'>':''}
            
            for cut in cut_list:
                for operator in cut_list_operators:    
                    #print 'second loop:\ncut:', cut, 'operator:', operator, 'cut_list_operators_dict:', cut_list_operators_dict[operator]
                    try:
                        data      = self.myData2Process[np.where(self.myData2Process['Z']>0.05)]
                        #data['Z'] = cd.comoving_distance(data['Z'], **fidcosmo)*0.6777
                        #print 'cosmology:', fidcosmo
                        #print data['Z'][0:20]
                    except:
 
                                             
                        data = myData.selectData2Compute(mydata, 
                                                        selected_col='orphan', 
                                                        operator=gtype_dict_operators[gtype], 
                                                        condition=gtype_dict[gtype])

                        #print 'min/max mstar:', format(min(self.myData2Process['mstar']), '.3e'), '/', format(max(self.myData2Process['mstar']), '.3e')

                        
                        #print 'after selection galaxy type:', data.shape  

                        if cut_list_operators_dict[operator].find('Mr')!=-1:
                            print 'here Mr > -21.5'
                            data = myData.selectData2Compute(data, 
                                                            selected_col='MAB_dA_total_r', 
                                                            operator='>', 
                                                            condition=-22.0)
#                            
                            data = myData.selectData2Compute(data, 
                                                            selected_col='MAB_dA_total_r', 
                                                            operator='<', 
                                                            condition=-21.0)
            
                        #small scales
            #            data = myData.selectData2Compute(data, 
            #                                            selected_col='mstar', 
            #                                            operator='<', 
            #                                            condition=1.73e11)
            
                        #large scales
                        if cut!='':
                            data = myData.selectData2Compute(data, 
                                                            selected_col='mstar', 
                                                            operator=operator, 
                                                            condition=cut)
                            
                            print operator, cut, ':after selection 1:', data.shape 
                            
                            if cut_list_operators_dict[operator]=='_mstar_' or cut_list_operators_dict[operator]=='_kroup_mstar_':
                                data = myData.selectData2Compute(data, 
                                                                selected_col='mstar', 
                                                                operator='<', 
                                                                condition=cut_list_max_dict[cut])
                                print '-->    <', cut_list_max_dict[cut], 'after selection 2:', data.shape
                        
                        if cut_list_operators_dict[operator]=='_BC_':
                            cut=2.35
                            print 'BC!'
                            data = data[(np.where(data['mAB_dA_total_g']-data['mAB_dA_total_i']<float(cut)))[:][0]]
                            
                        if cut_list_operators_dict[operator]=='_RS_':
                            print 'RS!'
                            cut=2.35
                            data = data[(np.where(data['mAB_dA_total_g']-data['mAB_dA_total_i']>float(cut)))[:][0]]                         #print 'after selection 2:', cut_list_max_dict[cut], data.shape                                                
            ##                                            
            ##
            #            data = myData.selectData2Compute(data, 
            #                                            selected_col='MAB_dA_total_r', 
            #                                            operator='>', 
            #                                            condition=-19)
            #            print 'after selection 2:', data.shape 
            
            
                        centrals=len(np.where(data['orphan']==0)[0])
                        sats=len(np.where(data['orphan']==1)[0])
                        orphans=len(np.where(data['orphan']==2)[0])
                        
                        #print 'galaxy types in sample: cents:', centrals, 'no-sats:', sats, 'orphans:', orphans
                        
                        #data      = self.myData2Process[np.where(self.myData2Process['orphan']>0)]
                        #print 'sats+orphans!'
                        #print 'select centrals! after selection:', data.shape,  'cents:', len(data[np.where(data['orphan']==0)]), 'sats', len(data[np.where(data['orphan']==1)]), 'orphans:',  len(data[np.where(data['orphan']==2)]), 'stats+orphans:',  len(data[np.where(data['orphan']>0)])
            
                    if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_ngal_random_sample']!='False':
                        data= mL.choose_random_sample(data, int(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_ngal_random_sample']))
                        print 'after random selection of', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_ngal_random_sample'], 'galaxies:', data.shape
            #        print np.info(data)
            #        print data['x_pos'][0:10]
            #        for f in ['x_pos', 'y_pos', 'z_pos']:
            #            data[f]/=0.6777
            #        print data['x_pos'][0:10]
            
                    #print 'units:', pos_unit, 'min/max: x,y,z', max(data['x_pos']), '/', max(data['x_pos']), max(data['y_pos']), '/', max(data['y_pos']), max(data['z_pos']), '/', max(data['z_pos'])
    
    
    
                    pi_max_list=[150]

                    for pi in pi_max_list:
                        
                        pimax     = float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_pimax'])
                        pimax=pi
                        nthreads  = int(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_nthreads'])


                        #large scales:    
                        calc_list_rmin=[0.5]
                        calc_dict_rmax={0.5: 150}
                        
                        #samle scales:
                        #calc_list_rmin=[0.5, 1, 5, 20]
                        #calc_dict_rmax={0.5: 3, 1: 10, 5: 75, 20:200}
                        
                        #final scales:
                        #calc_list_rmin=[0.5, 10]
                        #calc_dict_rmax={0.5: 10, 10: 200}                        
                        
                        for r_min in calc_list_rmin:                   
                            # Setup the bins
                            rmin      = float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_rmin'])
                            rmin=r_min
                            rmax      = float(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_rmax'])
                            rmax=calc_dict_rmax[rmin]
                            nbins     = int(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_nbins'])
                            
                            # Create the bins
                            rbins     = np.logspace(np.log10(rmin), np.log10(rmax), nbins)
                    
                            if self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_which']=='BOX':
                                if pos_unit.find('h')==-1:
                                    boxsize   = self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size']/0.6777
                                else:
                                     boxsize   = self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size']             
                    #            if self.myconfig_array['catname'+str(self.a)].startswith('SAGE'):
                    #                boxsize   = 350.0/0.6777+0.1
                                    
                                #print 'chosen BOXSIZE:', boxsize, 'rmin/rmax:', rmin, '/', rmax
                            else:
                                boxsize   = '-'
                                rand_RA, rand_DEC, rand_CZ = read_catalog(pjoin(dirname(abspath(cf.__file__)), "../mocks/tests/data", "Mr19_randoms_northonly.rdcz.ff"))
                                N         = len(data['RA'])            
                                rand_N    = len(rand_RA)
                                
                                DD_counts = calc_DD()
                                DR_counts = calc_DR()
                                RR_counts = calc_RR() 
                    
                           
                            
                            def caseSwitcher(name):
                        
                                choose = {
                                        'xi_of_r'       : xi_of_r,
                                        'xi'            : xi,
                                        'wp'            : wp,
                                        'xi_rp_pi'      : xi_rp_pi,
                                        'wp_mocks'      : wp_mocks,                   
                                        }
                                    
                                func = choose.get(name)
                                return func()
                                       
                    
                          
                            try:
                                #print 'calculate:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_calculate'],        
                                results_CF, calc_time = caseSwitcher(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_calculate'])           
                                writeIntoFile(results_CF,  self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_calculate'])
                                print '----------'
                                #print '... succsessful calculated! within:', calc_time, 'min /', calc_time, 'sec' 
                                
                            except:
                                pass
#                                a=0
#                                while a<10:
#                                    try:
#                                        print 'a:', a, 'calculate:', self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_calculate'+str(a)],
#                                        results_CF, calc_time = caseSwitcher(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_calculate'+str(a)])
#                                        writeIntoFile(results_CF,  self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_twoPCF_calculate'+str(a)])
#                                        print '... succsessful calculated! within:', calc_time/60.0, 'min /', calc_time, 'sec' 
#                                    except:
#                                        print 'no more functions to calculate ...'
#                                        break
#                                    a+=1                    

    def HODFunction(self,
                    mydata,
                    filename,
                    which_halomass=''):

        print ' '
        print ' '                 
        print 'HODFunction()', 'self.i:', self.i, 'self.a:', self.a
        print '++++++++++++++++++++++++++++++++++++++++++++++'
        print ' '
        #print 'catname', self.myconfig_array['catname'+str(self.a)], 'z:', self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'],  
        #print 'filename:', filename

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
            print 'all galaxies: centrals', len(all_centrals), 'sats:', len(all_sats), 'orphans:', len(all_orphans), '\n'
                    
            print 'all haloid:'
            print '-------------------'                  
            unique_haloid, index_haloid, count_haloid = find_unique_objects2(data['haloid'])
            print 'haloid', len(unique_haloid), 'ngal:', np.sum(count_haloid)
    
            centrals, sats, orphans = find_galaxy_type(data[index_haloid])                  
            print 'centrals', len(centrals), 'sats:', len(sats), 'orphans:', len(orphans), '\n'
            
            myhalo_data=data[['haloid', 'hostid', 'orphan', 'mhalo']][index_haloid]
            myhalo_data= mL.dim_expander_struct(myhalo_data, 'mhalo', 'ngal') 
            myhalo_data['ngal']=count_haloid
    
            print 'all hostid:'
            print '-------------------'        
            unique_hostid, index_hostid, count_hostid = find_unique_objects2(data['hostid'])
            print 'hostid', len(unique_hostid), 'ngal:', np.sum(count_hostid) 
    
            centrals, sats, orphans = find_galaxy_type(data[index_hostid])                  
            print 'centrals', len(centrals), 'sats:', len(sats), 'orphans:', len(orphans), '\n'


        def find_unique_objects(data):
                                            
             return np.unique(data, return_index=False, return_counts=True)


        def crossmacht_catalogs(data1, data2):
                      
            test = np.in1d(data1['haloid'], data2['haloid'])
                       
            #print 'test haloid:', test_haloid
            #print data[np.where(data['haloid']==test_haloid)]
            data_to_check=data[np.where(test==True)[:][0]]
            print data_to_check
            
            print 'ngal haloid==hostid', len(data_centrals), 'sats:', len(data_sats), 'orphans:', len(data_orphans), 'gal test=True:', data_to_check.size              
            

        def check_parents(data):
               
            data2write_not_found=[]
            data2write_found=[]
            
            id_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/data/prova.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.int64, skiprow=1) 
   
            header='\nHostHaloID\tMainHaloID\tPOS_X/POS_Y/POS_Z [h-1Mpc] GalaxyType HaloMass [h-1Msun] Mstar [h-1Msun]'
            for haloid, hostid in zip(id_array[:,0], id_array[:,1]):
            
                print 'satellite HostHaloID:', hostid, '--> find parent\t', haloid,
                mytest_data = data[np.where(data['hostid']==hostid)[:][0]]
                
                #print 'check from complete list! satellite hostid: MainHaloID:', mytest_data['haloid'], 'HostHaloID:', mytest_data['hostid']    
                #print 'find this satellite:', hostid, '--> and its parent', haloid, '-->', 
                
                parent=data[np.where(data['haloid']==haloid)[:][0]]
                parent=parent[np.where(parent['orphan']==0)[:][0]]
                print 'parent', parent['haloid'],
                            
                if parent['haloid']==haloid:                   
                    print '\t\t\tfound!'                
                    data2write_found.append(str(hostid)+'\t'+str(haloid)+'\t'+str(mytest_data['x_pos'])+'\t'+str(mytest_data['y_pos'])+'\t'+str(mytest_data['z_pos'])+'\t'+str(mytest_data['orphan'])+'\t'+str(mytest_data['mhalo'])+'\t'+str(mytest_data['mstar']))
                else:
                    print '\tNOT FOUND!'
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
            
            print 'total galaxies:', id_array[:,0].size, 'parents found', len(data2write_found), 'NOT FOUND:', len(data2write_not_found)
            
        
        #GalaxyTypeStats(self.myData2Process)
        
        #data= mL.choose_random_sample(self.myData2Process, 500000)        

#        mL.give_galaxy_type_info(self.myData2Process['orphan'])       
#        exit()        
    
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
        #print self.myData2Process.shape, 'ngal sats:', self.myData2Process[np.where(self.myData2Process['orphan']>0)[:][0]].size


        #choose random number of satellites!        
        #rand_sats = mL.choose_random_sample(self.myData2Process[np.where(self.myData2Process['orphan']>0)[:][0]], int(self.myData2Process[np.where(self.myData2Process['orphan']>0)[:][0]].size-26309))
        #test      = np.in1d(self.myData2Process['hostid'], rand_sats['hostid'])             
        #self.myData2Process      = self.myData2Process[np.where(test==False)[0][:]]

        redshift=float(self.mysnap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_filename'+str(self.i)]+'_snapid'+str(self.i)]['z'])
        print 'redshift', redshift
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
            #print 'bins:', bins
            Nhalos=data_halo[:,6]           
        else:
            binsize=0.103
            mhalo_min=12.1690
            mhalo_max=15.271
            
            binsize      = abs(mhalo_max-mhalo_min)
            bins         = np.arange(data_halo[0,0]-binsize/2.,data_halo[-1,0]-binsize/2., binsize)
  
        print 'Halos from --> ', path_to_halos
        
        if hod==True:
            pass
            #print '<N> = sum hosthalos + 1 / nhalos per mhalo bin'


#        print len(np.where(self.myData2Process['orphan']==0)[:][0])
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
                print 'i:', i, edges[i-1],

                if i==len(histo[:,0]):
                    #print 'here! MAX', edges[i-1], np.log10(histo[i-1,0]),
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
                    #print 'data_bin', data_bin.size 
                    test = np.in1d(data_sel[central_id],data_bin[central_id])
                    data_to_check=data_sel[np.where(test==True)[:][0]]
                    
                    #print 'N'+gtype+':', data_to_check.size
                    unique_haloid, count_haloid = find_unique_objects(data_to_check[central_id]) 
                    #print data_to_check[['haloid','hostid','orphan']]
                    #print count_haloid[np.where(count_haloid!=1)[:][0]]
 
                    print 'Nhaloid:', sum(count_haloid), 'Nunique:', len(unique_haloid), 'Nhalos:', float(Nhalos[i-1])
                        
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
                    
                    print 'calculate scatter: ', 
                    mean_mstar=np.mean(np.log10(data_sel['mstar']))
                    #scatter
                    histo[i-1,1]=np.std(np.log10(data_sel['mstar']))
                    histo[i-1,4]=histo[i-1,1]/(data_sel.size**0.5)
                    histo[i-1,5]=histo[i-1,4]
                    print histo[i-1,1], np.log10(histo[i-1,0])
                    
                    histo[i-1,6]=data_sel.size               
    
                i+=1


            #print histo[:,[0,1]]
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
                      binup2D=False):

        if binup2D==True:

            if col_name.find('mhalo')!=-1:
                data = data[np.where(data['orphan']==0)[:][0]]
                data = data[np.where(data[col_name]>=1e11)[:][0]]
                data = data[np.where(data[col_name]<1e15)[:][0]]

                histo, binsize = myFuncs.binUp(data[[col_name,'mstar']],
                                                    25,
                                                    histo_min=1e11,
                                                    histo_max=1e15,
                                                    log10bin=True,
                                                    binup_2D=True,
                                                    div_y2x=True,
                                                    use_MAD=True,
                                                    equal_bins=True) 
                #print histo[:,[0,1]]
            
        else:
                
            #print 'fast histogramm -->', col_name, '[', col_unit, ']',
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
            elif col_name.find('rhalfmass')!=-1:
                data_min = 1e-4
                data_max = 1                
            elif col_name=='vmax' or col_name=='vdisp':
                data_min = 10
                data_max = 10000
            elif col_name=='cgf':
                data_min = 1e-5
                data_max = 10                   
            else:
                data_min=min(data)
                data_max=max(data)
                
            data = data[np.where(data>=data_min)[:][0]]
            data = data[np.where(data<data_max)[:][0]]

            if binning=='log':           
                data = np.log10(data)
                data_min=np.log10(data_min)
                data_max=np.log10(data_max)
                  
            data = data[np.where(np.isfinite(data))[:][0]]        
    #        data = data[np.where(data>=np.percentile(data,0.01))[:][0]]
    #        data = data[np.where(data<np.percentile(data,99.9))[:][0]]
            #print data.shape                
            #print 'data min/max:', min(data), '/', max(data),
            nbins=25
    
            #print 'data_min:', data_min, 'data_max:', data_max
            binsize = (data_max-data_min)/nbins
            
            bins = np.linspace(data_min, data_max, nbins+1)
            #print bins
           
            counts, edges = np.histogram(data, bins)
    
            histo=np.zeros((nbins, 7), dtype=np.float32)
           
            if binning =='log':
                histo[:,0]= (10**edges[1:]*10**edges[:-1])**0.5
            else:
                histo[:,0]= (edges[1:]+edges[:-1])/2
    
            if self.volume=='False' or self.volume==False:
                self.volume=(self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_box_size']/0.6777)**3
                print '--> set volume:', self.volume
                
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
       

                