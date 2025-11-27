from __future__ import print_function
import arangeData as aD
myData = aD.ArangeData()

import myFuncs as mF
myFuncs = mF.MyFunctions()

import myLib as mL
import dataPipeline as dP
import loadObs as lO
import outputData as oD

myOutput = oD.OutputData(config=False)

# system
import numpy as np
import os
import time
system_info = os.getcwd()
start = system_info.find('anaconda')
mycomp = system_info[:start]
import pandas as pd
ts = time.time()
import datetime 
date_time_stamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H:%M:%S')

class DoAnalysis:

    def __init__(
                self,
                myconfig_array,
                snap_array,
                plot_map_array,
                physics_specs,
                histo_configs,
                mycond_configs,
                SAM_scale_factor_map,
                name_conv_map):
        
        self.myconfig_array   = myconfig_array
        self.snap_array   = snap_array
        self.plot_map_array = plot_map_array
        self.physics_specs  = physics_specs
        self.histo_configs  = histo_configs
        self.mycond_configs = mycond_configs
        self.SAM_scale_factor_map = SAM_scale_factor_map
        self.name_conv_map   = name_conv_map

        self.mycomp         = mycomp
        self.plot_all       = self.physics_specs['plot_all']
        self.single_mix     = self.physics_specs['single_mix']
        self.subplots       = self.physics_specs['subplots']
        self.single         = self.physics_specs['single']
        self.load_from_file = self.physics_specs['load_from_file']
        self.plot_custom_loop = self.physics_specs['plot_custom_loop']
        
        self.myPipe = dP.Pipeline(self.myconfig_array, 
                                  self.snap_array, 
                                  self.mycond_configs,
                                  self.plot_map_array,
                                  self.SAM_scale_factor_map,
                                  self.name_conv_map)

    def cSFRD(self,
            plot_key,
            data,
            filename): 

        my_ssfr_cut=None
        my_sfr_cut_min=0
        my_sfr_cut_max='max'
   
        self.myPipe.cSFRD(data, 
                          data_offset=int(self.plot_map_array['data_offset_'+plot_key]),
                          histo_data_min=self.histo_configs['histo_min_cSFRD_mstar'],
                          histo_data_max=self.histo_configs['histo_max_cSFRD_mstar'],
                          mycond_min=self.histo_configs['mycond_min_cSFRD_mstar'], 
                          mycond_max=self.histo_configs['mycond_max_cSFRD_mstar'],
                          ssfr_cut=my_ssfr_cut,
                          sfr_cut_min=my_sfr_cut_min,
                          sfr_cut_max=my_sfr_cut_max)


        myOutput.writeIntoFile(filename,
                               self.myPipe.histo_data_ZvsSFR[:, self.b*int(self.plot_map_array['data_offset_'+plot_key]) : self.b*int(self.plot_map_array['data_offset_'+plot_key])+int(self.plot_map_array['data_offset_'+plot_key])],
                               myheader= plot_key+' '+self.myconfig_array[self.myconfig_array['catname'+str(self.b)]+'_simulation_name']+' '+self.myconfig_array['catname'+str(self.b)]+\
                               'mstar min/max: '+str(self.histo_configs['mycond_min_cSFRD_mstar'])+'/'+str(self.histo_configs['mycond_max_cSFRD_mstar'])+' [Msun], ssfr cut: '+str(my_ssfr_cut)+' [yr-1], sfr min/max: '+str(my_sfr_cut_min)+'/'+str(my_sfr_cut_max)+'[Msun yr-1]'+\
                               r'\n(1) z (2) cSFRD (sumSFR/volume) [Msun yr-1 Mpc-3] (3) sumSFR [Msun yr-1] (4) cSFRD/cSFRD(z=0) (5) comoving volume [Mpc3]',
                               data_format="%0.3f\t%0.8e\t%0.8e\t%0.8e\t%0.2e")

        return self.myPipe.histo_data_ZvsSFR  

        
    def sfr2z(
            self,
            plot_key,
            data,
            filename,
            tarsel_code_space,
            output_filename_code_space): 

        self.myPipe.filterAndExtractData(self.mycond_configs,
                                         preprocessing_only=True)
        
        mymethod='300'

        if mymethod=='M2':       
            #Galacticus
            #SFH_analysis_keys=['','low', 'high', 'passive', 'active', 'red', 'blue', 'low-zcold', 'high-zcold', 'mhalo12gt', 'mhalo12st', 'zcold9gt', 'zcold9st', 'mstar11gt', 'mstar10st', 'redmstar11', 'redmhalo13', 'lowmstar11', 'lowmhalo13']
            SFH_analysis_keys=['','low-zcold','low-zcold-red','high-zcold','high-zcold-red']
            SFH_analysis_keys=['','mstarst11','mstarge11','mhalost13','mhaloge13','lowmstarst11','lowmhalost13','redmstarst11','redmstarge11','redmhalost13','redmhaloge13','low-zcold-red20','high-zcold-red20']
            
        elif mymethod=='300':
            SFH_analysis_keys=['']
    
            # d=1
            # while d<=324:
            #     SFH_analysis_keys.append('r'+str(d))
            #     d+=1
                
            # print(SFH_analysis_keys

        else:
            SFH_analysis_keys=['','low', 'high', 'passive', 'active', 'red', 'blue', 'low-zcold', 'high-zcold']
            
        #SFH_analysis_keys=['r0001']
        
        #SFH_analysis_keys=['','low-zcold-red-dens','low-zcold-red']
        #LGALAXIES
        #SFH_analysis_keys=['','low', 'high', 'passive', 'active']        
 
        #SFH_analysis_keys=['lowZcold-highMstar']
        #SFH_analysis_keys=[''] 
       
        for count, element in enumerate(SFH_analysis_keys):
            print('\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
            
            for env_code, env_name in zip([-99],['']):
            #for env_code, env_name in zip([-99,3,2],['','k','f']):
                
                element_space=''        
                #create spaces '_' in the filename:
                if element!='':
                    element_space='_' 
                
                env_space='' 
                if env_name!='':
                    env_space='_'
                
                myfilename=filename[0:len(filename)-4]+element_space+element+env_space+env_name+'.txt'
                print('myfilename:',  myfilename)
                error_count=0
                if self.myPipe.i>0:
                    filename_before = self.mycomp+'anaconda/pro/myRun/histos/'+plot_key+'/'+plot_key+'_'\
                                        +self.myconfig_array['catname'+str(self.myPipe.a)]\
                                        +'_z_'+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i-1)]+'_snapid'+str(self.myPipe.i-1)]['z'])))\
                                        +'_tarsel'+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']\
                                        +output_filename_code_space+str(self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_output_filename_code'])\
                                        +element_space+element+env_space+env_name+'.txt'
                    #print('filename before:', filename_before
                    data = myData.readAnyFormat(config=False, mypath=filename_before, data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float64, skiprow=2) 
        
                    
                
                nerrors, myprops, perc_high, perc_low = self.myPipe.SFR2Z(data,
                                                              myfilename,
                                                              count,
                                                              error_count,
                                                              self.mycond_configs,
                                                              data_offset=int(self.plot_map_array['data_offset_'+plot_key]),
                                                              histo_data_min=self.histo_configs['histo_min_sfr2z_mstar'],
                                                              histo_data_max=self.histo_configs['histo_max_sfr2z_mstar'],
                                                              mycond_min=self.histo_configs['mycond_min_sfr2z_mstar'], 
                                                              mycond_max=self.histo_configs['mycond_max_sfr2z_mstar'],
                                                              SFH_key=element,
                                                              method=mymethod,
                                                              env_name=env_name,
                                                              env_code=env_code)
                                          
                print('\n++++++++++++++++++++++++++++++++++++++++++++++++++\nOUTPUT -->', end=' ')
                
                myheader=plot_key+' selection key: '+element+' '+self.myconfig_array[self.myconfig_array['catname'+str(self.b)]+'_simulation_name']+' '+self.myconfig_array['catname'+str(self.b)]+' error main progenitor search: '+str(nerrors)+'\n'
                myformat=''
                
                #Create header and format of the output file (ASCII)
                property_map=mL.get_property_dict2()
                                
                count=1

                for i,item in enumerate(myprops):
                    #print('i:', i, 'count:', count, 'item:', item)
                    new_prop=property_map[item]['avg_calc']+property_map[item]['output_prop_as']+'['+property_map[item]['unit']+']'+'('+str(count)+') '

                    myformat+=property_map[item]['format']+' '
                    if property_map[item]['avg_calc']!='':
                        uncert=perc_high+'th'+property_map[item]['output_prop_as']+'['+property_map[item]['unit']+']'+'('+str(count+1)+') '+\
                               perc_low +'th'+property_map[item]['output_prop_as']+'['+property_map[item]['unit']+']'+'('+str(count+2)+') '
                        new_prop+=uncert
                        
                        myformat+=property_map[item]['format']+' '
                        myformat+=property_map[item]['format']+' '
                        count+=2
                    count+=1
                    
                    myheader+=new_prop

                #print(myheader)
                #print(myformat[:-1])
                myOutput.writeIntoFile(myfilename,
                                       self.myPipe.histo_data_ZvsSFR[:, self.b*int(self.plot_map_array['data_offset_'+plot_key]) : self.b*int(self.plot_map_array['data_offset_'+plot_key])+int(self.plot_map_array['data_offset_'+plot_key])],
                                       myheader=myheader,
                                       data_format=myformat[:-1])                  

        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        return self.myPipe.histo_data_ZvsSFR                             
           
    def initBinAndFrac2D(self,
                         plot_key,
                         data,
                         filename,
                         my_name_x,
                         my_name_y=False,
                         mylog10bin=True,
                         mybinUp2D=False,
                         mydiv_y2x=False,
                         myadd_axis=False,
                         mycumulative=False,
                         mynormalise=False,
                         mystarforming_cut=False):                                                            

        try:
            cond_min_x=self.histo_configs['mycond_min_'+plot_key+'_'+my_name_x]
        except:
            cond_min_x='min'
        try:
            cond_max_x=self.histo_configs['mycond_max_'+plot_key+'_'+my_name_x]
        except:
            cond_max_x='max'
        try:
            cond_min_y=self.histo_configs['mycond_min_'+plot_key+'_'+my_name_y]
        except:
            cond_min_y='min'
        try:
            cond_max_y=self.histo_configs['mycond_max_'+plot_key+'_'+my_name_y]
        except:
            cond_max_y='max'
        try:
            histo_min=self.histo_configs['histo_min_'+plot_key+'_'+my_name_x]
        except:
            histo_min='min'
        try:
            histo_max=self.histo_configs['histo_max_'+plot_key+'_'+my_name_x]
        except:
            histo_max='max'

        if mydiv_y2x!=False or myadd_axis!=False:
            mydata_offset=10
            data=mL.dim_expander(data,mydata_offset)
        else:         
            mydata_offset=int(self.plot_map_array['data_offset_'+plot_key])

        if mystarforming_cut!=False:
            starforming_cut=' starforming cut: '+str(mystarforming_cut)+' '
        else:
            starforming_cut=''

        if self.mycond_configs[my_name_x+'_unit'].find('h-')!=-1:
            myvolume_units='h3 Mpc-3'
        else:
            myvolume_units='Mpc-3'

        if self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_calc_fast_histo']=='True':
            
            mybinning='log'
            for prop in ['sfr']:
                print('prop:', prop)
                if prop=='SHMF':
                    self.myPipe.myData2Process = self.myPipe.myData2Process[np.where(self.myPipe.myData2Process['orphan']<2)[:][0]]
                    self.myPipe.myData2Process['mstar']/=self.myPipe.myData2Process['mhalo']
                    
                    self.mycond_configs[prop+'_unit'] = '-'
                elif prop.find('AB')!=-1:
                    mybinning='lin'
                elif prop.find('ssfr')!=-1:
                    self.myPipe.myData2Process['sfr']=self.myPipe.myData2Process['sfr']/self.myPipe.myData2Process['mstar']
                    
                    self.mycond_configs[prop+'_unit'] = '[yr-1]'                  

                self.myPipe.calcFastHisto(self.myPipe.myData2Process[prop],
                                          filename[0:len(filename)-4]+'_'+prop+'.txt',
                                          prop,
                                          self.mycond_configs[prop+'_unit'],
                                          binning=mybinning)
                
            exit()
        else:

            self.myPipe.binAndFrac2D(data, 
                                    data[:,0].size,
                                    my_name_x,
                                    name_y=my_name_y,
                                    mycond_min_x=cond_min_x,
                                    mycond_max_x=cond_max_x,
                                    mycond_min_y=cond_min_y,
                                    mycond_max_y=cond_max_y,
                                    histo_data_min=histo_min, 
                                    histo_data_max=histo_max, 
                                    data_offset=mydata_offset,
                                    log10bin=mylog10bin,
                                    div_y2x=mydiv_y2x,
                                    binUp2D=mybinUp2D,
                                    add_axis=myadd_axis,
                                    plot_key=plot_key,
                                    cumulative=mycumulative,
                                    normalise=mynormalise,
                                    starforming_cut=mystarforming_cut,
                                    volume_units=myvolume_units)


        if mycumulative==True:
            info_cum='YES'
            Phi='n(>)'
            dex=']'
        else:
            info_cum='NO'
            Phi='Phi'
            dex=' dex-1]'

        if plot_key.find('zgas2')!=-1:
            my_name_y='Zcold'
            self.mycond_configs[my_name_y+'_unit'] = '-'
            if plot_key=='zgas2mcold':
                my_name_x='mcold/mstar'
                self.mycond_configs[my_name_x+'_unit'] = '-'



        if plot_key.find('vsSFR')!=-1:
            header=r'\n(1) '+my_name_x+' ['+self.mycond_configs[my_name_x+'_unit']+'] (2) mstar ['+self.mycond_configs[my_name_y+'_unit']+'] / '+my_name_y+' ['+self.mycond_configs[my_name_y+'_unit']+']\t(3) -dx\t(4) dx\t(5) -dy\t(6) +dy\t(7) N count (8) add axis '+myadd_axis+' ['+self.mycond_configs[myadd_axis+'_unit']+'] (9) dy add axis'               
        elif mydiv_y2x==True or myadd_axis!=False:           
            header=r'\n(1) '+my_name_x+' ['+self.mycond_configs[my_name_x+'_unit']+'] (2) '+my_name_y+' ['+self.mycond_configs[my_name_y+'_unit']+']\t(3) -dx\t(4) dx\t(5) -dy\t(6) +dy\t(7) N count (8) add axis '+myadd_axis+' ['+self.mycond_configs[myadd_axis+'_unit']+'] (9) -d add axis (10) +d add axis'
        elif mybinUp2D==True:
            header=r'\n(1) '+my_name_x+' ['+self.mycond_configs[my_name_x+'_unit']+'] (2) '+my_name_y+' ['+self.mycond_configs[my_name_y+'_unit']+']\t(3) -dx\t(4) dx\t(5) -dy\t(6) +dy\t(7) N count'
        else:
            header=r'\n(1) '+my_name_x+' ['+self.mycond_configs[my_name_x+'_unit']+'] (2) '+Phi+' ['+myvolume_units+dex+r'\t(3) -dx\t(4) dx\t(5) -dy\t(6) +dy\t(7) N count'
      
        #total number of objects after all cuts used in calculation
        ngal_in_histo=sum(self.myPipe.ratio_data_binAndFrac2D[:, mydata_offset*self.myPipe.a+6])
        
        myOutput.writeIntoFile(filename,
                               self.myPipe.ratio_data_binAndFrac2D[:, mydata_offset*self.myPipe.a : mydata_offset*self.myPipe.a + mydata_offset],
                               myheader= plot_key+' '+self.myconfig_array['catname'+str(self.myPipe.a)]+' z='+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))+' cumulative: '+info_cum+starforming_cut+' ngal: '+str(int(ngal_in_histo))+header,
                               data_format="%0.8e",
                               mydelimiter=r'\t')

        return self.myPipe.ratio_data_binAndFrac2D
               
    def MAIN(self):

        def caseSwitcher(plot_key, data, filename):
        
            choose = {
                'SMF':                     MassFunction,
                'HMF':                     HMF,
                'HMF_no':                  HMF,
                'SFRF':                    SFRFunction,
                'sSFRF':                   sSFRFunction,
                'sfr2z':                   SFR2Z,
                'cSFRD':                   cSFRD,
                'ssfr2mstar':              sSFR2Mstar,
                'oh2mstar':                OH2Mstar,                
                'mstar2mhalo':             Mstar2Mhalo,
                'mstar2mhalo_no':          Mstar2Mhalo,                
                'zgas2mstar':              Zgas2Mstar,
                'zgas2mcold':              Zgas2Mcold,
                'mcold2mstar':             Mcold2Mstar,
                'mbh2mstarsph':            Mbh2Mstarsph,
                'twoPCF':                  TwoPCF,
                'HOD':                     HODFunction,
                'filterData':              FilterAndExtractData,
                'analyseTargetSelection':  AnalyseTargetSelection,
                'plotXY':                  PlotXY,
                'mstar2mhalovsSFR':        Mstar2MhalovsSFR,
                'mstar2rhalf':             Mstar2Rhalf                
                }
                
            func = choose.get(plot_key)
            
            return func(plot_key, data, filename)
       
        def MassFunction(plot_key, data, filename):
            return self.initBinAndFrac2D(plot_key, data, filename, 'mstar', mynormalise=True)

        def HMF(plot_key, data, filename):
            return self.initBinAndFrac2D(plot_key, data, filename, 'mhalo_cents_200c', mycumulative=False, mynormalise=True)

        def SFRFunction(plot_key, data, filename):
            return self.initBinAndFrac2D(plot_key, data, filename, 'sfr', mynormalise=True)
            
        def sSFRFunction(plot_key, data, filename):
            return self.initBinAndFrac2D(plot_key, data, filename, 'ssfr', mynormalise=True)
        
        def cSFRD(plot_key, data, filename):
            return self.cSFRD(plot_key, data, filename)
            
        def SFR2Z(plot_key, data, filename):
            return self.sfr2z(plot_key, data, filename, tarsel_code_space, output_filename_code_space)
            
        def sSFR2Mstar(plot_key, data, filename):
            return self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'sfr', mybinUp2D=True, mystarforming_cut=1e-11)
            
        def OH2Mstar(plot_key, data, filename):
            #mystarforming_cut default=1e-11, z=0.1 better to use the formular 0.3/t_hubble(z=0.1)=3.547e-11 (quisent cut)
            starforming_cut=mL.starforming_cut_Henriques20_A1(redshift=float(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z']))           
            
            if self.myconfig_array['catname'+str(self.myPipe.a)].startswith('SAGE'): 
                return self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'zgas_disk', mybinUp2D=True, myadd_axis='mstar_disk', mystarforming_cut=starforming_cut)                
            else:                                
                return self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'zcold', mybinUp2D=True, myadd_axis='mstar_disk', mystarforming_cut=starforming_cut)
 
        def Mstar2Rhalf(plot_key, data, filename):
            return self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'rhalf_mass', mybinUp2D=True, mydiv_y2x=False, myadd_axis='rbulge')
            
        def Mstar2Mhalo(plot_key, data, filename):
            return self.initBinAndFrac2D(plot_key, data, filename, 'mhalo', 'mstar', mybinUp2D=True, mydiv_y2x=True, myadd_axis='mstar')

        def Mstar2MhalovsSFR(plot_key, data, filename):
            return self.initBinAndFrac2D(plot_key, data, filename, 'ssfr', 'mhalo', mybinUp2D=True, mydiv_y2x=False, myadd_axis='mstar')
            
        def Zgas2Mstar(plot_key, data, filename):
            starforming_cut=mL.starforming_cut_Henriques20_A1(redshift=float(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z']))           
            starforming_cut=1e-25       
            if self.myconfig_array['catname'+str(self.myPipe.a)].startswith('SAGE'): 
                return self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'mcold_disk', mybinUp2D=True, myadd_axis='mstar_disk', mystarforming_cut=starforming_cut)
            else:    
                data = self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'zcold', mybinUp2D=True, myadd_axis='ssfr', mystarforming_cut=starforming_cut)

        def Zgas2Mcold(plot_key, data, filename):
            if self.myconfig_array['catname'+str(self.myPipe.a)].startswith('SAGE'): 
                return self.initBinAndFrac2D(plot_key, data, filename, 'mcold_disk', 'Mzgas', mybinUp2D=True, myadd_axis='mstar')
            elif self.myconfig_array['catname'+str(self.myPipe.a)].startswith('Galacticus'):
                return self.initBinAndFrac2D(plot_key, data, filename, 'mcold', 'zgas_spheroid', mybinUp2D=True, myadd_axis='mstar')
            else:    
                return self.initBinAndFrac2D(plot_key, data, filename, 'mcold', 'Mzgas', mybinUp2D=True, myadd_axis='mstar')
                
        def Mcold2Mstar(plot_key, data, filename):
            if self.myconfig_array['catname'+str(self.myPipe.a)].startswith('SAGE'):             
                return self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'mcold_disk', mybinUp2D=True, mydiv_y2x=True, myadd_axis='mcold_disk')
            else:
                return self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'mcold', mybinUp2D=True, mydiv_y2x=True, myadd_axis='mcold')
                
        def Mbh2Mstarsph(plot_key, data, filename):
            return self.initBinAndFrac2D(plot_key, data, filename, 'mstar_spheroid', 'mbh', mybinUp2D=True)
            
            
        def TwoPCF(plot_key, data, filename):          
            #self.myPipe.TwoPCF_new(False, filename, self.mycond_configs)
            self.myPipe.TwoPCF(False, filename, self.mycond_configs)

        def HODFunction(plot_key, data, filename):          
            self.myPipe.HODFunction(False, filename)
                       
        def FilterAndExtractData(plot_key, data, filename):
            self.myPipe.filterAndExtractData(self.mycond_configs)
            
        def AnalyseTargetSelection(plot_key, data, filename):
            return self.myPipe.analyseTargetSelection(data, filename, self.mycond_configs, plot_key)
           
        def PlotXY(plot_key, data, filename):
            self.myPipe.plotXY(plot_key,filename)
          
        def caseSwitcherPlotKey(plot_key, sample, key, prop, mysamples, method, plot_num, mycatalog):

            def myLoadFromFile():
                LoadFromFile(sample, key, prop, mysamples, method, plot_num, mycatalog)
            
            choose = {
                    'loadFromFile':           myLoadFromFile,
                    'mainCalculate':          MainCalculate,
                    'analyseTargetSelection': TargetSelection              
                    }
                
            func = choose.get(plot_key)
            return func()

        def check_filename():
            
            if self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_use_store_register']=='yes':
                register_path = self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_use_store_registerpath']
            else:
                register_path = self.mycomp+'/anaconda/pro/data/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'/'
                
            filename1 = register_path+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_'+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))+'_tarsel'+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']+'.'+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_fileformat']           
            filename2 = register_path+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_'+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']+'.'+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_fileformat']
            
            #filename3 = '/store/erebos/doris/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_'+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))+'_tarsel'+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']+'.'+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_fileformat']                
            #filename4 = '/store/erebos/doris/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_'+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']+'.'+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_fileformat']
            
            if self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_set_register_name'].find('Skies')!=-1:
                filename3=register_path+'/SkiesANDUniverses/MDPL2_'\
                            +self.myconfig_array['catname'+str(self.myPipe.a)][0:len(self.myconfig_array['catname'+str(self.myPipe.a)])-5] \
                            +'_z_'+str(format(float(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z']), '.2f'))+'.hdf5'
            else:
                filename3 = register_path+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_set_register_name']+'/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_'+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']+'.'+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_fileformat']
            
            filename4 = register_path+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_set_register_name']+'/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_'+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))+'_tarsel'+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']+'.'+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_fileformat']
                           
            #print('myfilename:', myfilename1, 'skip reading data? --> ', self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_skip_reading_data']
            
            return filename1, filename2, filename3, filename4#, filename5, filename6

        def MainCalculate():

            print('########################################################################################################\n#')                                                                                                   
            print('#     PROGRESS STATUS: MAIN CALCULATE:  catname:', self.myconfig_array['catname'+str(self.myPipe.a)], 'b:', self.b, 'i:', self.myPipe.i, 'a:', self.myPipe.a, 'plot_key:', plot_key)                                                                                                   
            print('#\n########################################################################################################\n')

            if plot_key=='sfr2z':
                preproc=True
            else:
                preproc=False
                        
            if self.b==0 or plot_key=='sfr2z' or plot_key=='cSFRD' or plot_key=='filterData' or plot_key=='plotXY' or plot_key=='twoPCF':
                if self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_data_format']!='HDF5':# and self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_data_format']!='CROSSMATCH':                     
                    self.myPipe.readData(preprocessing_only=preproc)
                                
                else:
                    myfilename1, myfilename2, myfilename3, myfilename4 = check_filename()
                    #print('myfilename1:', myfilename1)
                    #uncommand for test reason:
                    #--------------------
                    #myfilename=myfilename1
                    #self.myPipe.readData(myfilename=myfilename1)


                    if self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_skip_reading_data']!='yes':
                        if plot_key=='twoPCFd':
                            myfilename1 = self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_twoPCF_path_to_data']+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_'+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']+'.'+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_fileformat']                     
                        for filename in [myfilename1, myfilename2, myfilename3, myfilename4]:                            
                            #try:
                                print('try: myfilename:', filename)
                                self.myPipe.readData(myfilename=filename,
                                                     preprocessing_only=preproc)
                                myfilename=filename
                                break
    
                            # except:
                            #     print(filename, 'not found ... --> check for MEMORY ERRORS!')
                            
            if plot_key=='plotXY':
                filename=myfilename
            else:
                try:
                    filename = self.mycomp+'anaconda/pro/myRun/histos/'+plot_key+'/'+plot_key+'_'\
                                +self.myconfig_array['catname'+str(self.myPipe.a)]\
                                +'_z_'+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))\
                                +'_tarsel'+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']\
                                +output_filename_code_space+str(self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_output_filename_code'])\
                                +'.txt'
                except:
                    print('FILENAME was not found --> set to NONE!')
                    filename=None

            #self.myPipe.showConfig()
           
            self.my_analysed_data = caseSwitcher(plot_key, data, filename)                

        def TargetSelection():
            
            print('########################################################################################################\n#\n')                                                                                                 
            print('#     PROGRESS STATUS: MAIN TargetSelection():  catname:', self.myconfig_array['catname'+str(self.myPipe.a)], 'b:', self.b, 'i:', self.myPipe.i, 'a:', self.myPipe.a, 'plot_key:', plot_key)                                                                                                    
            print('#\n########################################################################################################\n')                        

            self.analysed_data_map['name'+str(self.myPipe.a)] = 'analysed_data_array'+str(self.myPipe.a)
            self.analysed_data_map['catname'+str(self.myPipe.a)] = self.myconfig_array['catname'+str(self.myPipe.a)]
            if self.myPipe.a==0:            
                self.analysed_data_map['mytext'] = ''
            
            filename1, filename2, filename3, filename4 = check_filename()

            print('loading ....', end='') 

            #for filename in [filename1, filename2, filename3, filename4, filename5, filename6]:
                #try:
            self.analysed_data_map['name'+str(self.myPipe.a)+'_data'], \
            self.analysed_data_map['name'+str(self.myPipe.a)+'_legend'], \
            self.analysed_data_map['col_name_x'], \
            self.analysed_data_map['col_name_y'], \
            self.analysed_data_map['col_name_z'], \
            self.analysed_data_map['col_name_weights'], \
            self.analysed_data_map['name_x'], \
            self.analysed_data_map['name_y'], \
            self.analysed_data_map['name_z'], \
            self.analysed_data_map['name_weights'] = caseSwitcher('analyseTargetSelection', \
                                                             self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'], \
                                                             filename1)
                        #print('filename:', filename
                 #   break
                #except:
                  #  print(filename, 'not found ...'
                               
            
            
            self.redshift = self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z']
            print('z:', self.redshift)
            #self.analysed_data_map['mytext']+=self.analysed_data_map['catname'+str(self.myPipe.a)][0:self.analysed_data_map['catname'+str(self.myPipe.a)].find('_')]+' ngal: '+str(len(self.analysed_data_map['name'+str(self.myPipe.a)+'_data']))+r'\n'
            #self.analysed_data_map['mytext']+=self.analysed_data_map['catname'+str(self.myPipe.a)]+' ngal: '+str(len(self.analysed_data_map['name'+str(self.myPipe.a)+'_data']))+r'\n'
            #print(self.analysed_data_map)
            #if self.myPipe.a==self.myconfig_array['nr_cats']-1: 
            plotAnalysedData(self.analysed_data_map)
                      
        def plotAnalysedData(analysed_data_map):

            print('########################################################################################################\n#\n')
            print('#     PROGRESS STATUS: MAIN plotAnalysedData():  catname:', self.myconfig_array['catname'+str(self.myPipe.a)], 'b:', self.b, 'i:', self.myPipe.i, 'a:', self.myPipe.a, 'plot_key:', plot_key)
            print('#\n########################################################################################################\n')
    
            LoadObs()
            #print(analysed_data_map)
            
            
            mycatname=self.myconfig_array['catname'+str(self.myPipe.a)][0:self.myconfig_array['catname'+str(self.myPipe.a)].find('_')]
            myredshift=self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z']
            #print('mycatname:', mycatname, 'myredshift:', myredshift)

            
            mynr_cats=self.myconfig_array['nr_cats']
            myobs={}
            myobs_legend={}
            subplot_loops=self.subplot_loops
            i=0
            a=0           
            while i<self.subplot_loops:
                              
                if self.subplot_is_contour[str(i)]==True:
                    try:
                        obs_data=np.zeros((obs_data_array[0][str(i)][:,0].size,), 
                                          dtype=[(analysed_data_map['col_name_x'], 'f8'), 
                                                 (analysed_data_map['col_name_y'], 'f8'),
                                                 (analysed_data_map['col_name_weights'], 'f8')])
                        
                        obs_data[analysed_data_map['col_name_x']]=obs_data_array[0][str(i)][:,0][:]                    
                        obs_data[analysed_data_map['col_name_y']]=obs_data_array[0][str(i)][:,1][:]
                        obs_data[analysed_data_map['col_name_weights']]=obs_data_array[0][str(i)][:,1][:]
                    except:
                        obs_data=np.zeros((obs_data_array[0][str(i)].size,), 
                                          dtype=[(analysed_data_map['col_name_x'], 'f8'), 
                                                 (analysed_data_map['col_name_y'], 'f8'),
                                                 (analysed_data_map['col_name_weights'], 'f8')])

                        obs_data[analysed_data_map['col_name_x']]=obs_data_array[0][str(i)][analysed_data_map['col_name_x']]               
                        obs_data[analysed_data_map['col_name_y']]=obs_data_array[0][str(i)][analysed_data_map['col_name_y']]
                        obs_data[analysed_data_map['col_name_weights']]=obs_data_array[0][str(i)][analysed_data_map['col_name_weights']]

                    analysed_data_map.update({'catname'+str(a+1): 'obs'+str(i), 
                                              'col_name_y': analysed_data_map['col_name_y'], 
                                              'col_name_x': analysed_data_map['col_name_x'],
                                              'col_name_weights': analysed_data_map['col_name_weights'],                                               
                                              'name'+str(a+1)+'_data': obs_data, 
                                              'name'+str(a+1)+'_legend': obs_data_array[1][str(i)]})
                                                                                 
                subplot_loops+=1

                myobs.update({str(a+1): obs_data_array[0][str(i)]})               
                myobs_legend.update({str(a+1): obs_data_array[1][str(i)]})
                a+=1
                i+=1
   
            myOutput.multiPlot(analysed_data_map,
                                 loops=mynr_cats,
                                 my_add_subplot= myobs,
                                 nr_added_subplots = subplot_loops,
                                 legend_added_subplots = myobs_legend,
                                 config_path=mycomp+'anaconda/pro/myRun/plot_config/', 
                                 myconfig_datafile=self.plot_map_array[plot_key+'_config'],
                                 mydir=mycomp+'anaconda/pro/myRun/plots/analyseTargetSelection/',
                                 myfilename= self.myconfig_array[self.myconfig_array['catname0']+'_simulation_name']+'_'+self.myconfig_array['catname0']+'_'+plot_key+'_z_'+str(myredshift)+'_'+str(date_time_stamp),
                                 data_block_offset=int(self.plot_map_array['data_offset_'+plot_key]),
                                 data_block_subplot_offset=5,
                                 plot_key=plot_key,
                                 print_redshift='z=0.56')#pop (ii), centrals')#, centrals')#\n$SFR<-1$')#CMASS DR12:')#$orphans$')# z='+str(float("{0:.2f}".format(myredshift))))#+' centrals')

        def LoadObs():

            if self.plot_map_array['load_obs_'+plot_key]!='False':
                
                key = mL.multicolTestAlgorithm(self.plot_map_array['load_obs_'+plot_key])
                #print('here:', self.plot_map_array['load_obs_'+plot_key])
                legend={}
                array={}
                self.subplot_is_contour={}
                x=0
                while x<key.size:

                    obs_key = getattr(lO, key[x])                 
                    load_obs = obs_key(load_from_file=False)
                    #print('load_obs:', load_obs                                        
                    array[str(x)]=load_obs[0]
                    legend[str(x)]=load_obs[1]
                    self.subplot_is_contour.update({str(x): load_obs[2]})
                    x+=1
                               
                obs_data_array[0] = array
                obs_data_array[1] = legend
                self.subplot_loops=x
                

           
        def LoadFromFile(mysample, mykey, myprop, mysamples, mymethod, myplot_num, mycatalog):

            self.print_redshift=False
            if self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_use_snapidz_mapping']=='False':
                self.redshift = self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_manual_input_redshift']
            else:
                self.redshift = self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]]['z']
                                           
            if plot_key!='plotXY':
               
                print('########################################################################################################\n')
                print('#     PROGRESS STATUS: MAIN LOAD FROM FILE:  catname:', self.myconfig_array['catname'+str(self.myPipe.a)], 'b:', self.b, 'i:', self.myPipe.i, 'a:', self.myPipe.a, 'plot_key:', plot_key)
                print('#\n########################################################################################################\n')

                if plot_key.find('plotOnly')!=-1:
                    filename=self.mycomp+'anaconda/pro/myRun/histos/analyseTargetSelection/histo_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_'+str(self.redshift)+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']+'_'+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_id_col_array']['name'+str(self.b)]+output_filename_code_space+str(self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_output_filename_code'])+'.txt'                          
                  
                else:                            
                    filename = self.mycomp+'anaconda/pro/myRun/histos/'+plot_key+'/'+plot_key+'_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_'+str(self.redshift)+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']+output_filename_code_space+str(self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_output_filename_code'])+'.txt'
                
                data_array = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim=r'\t', mydtype=np.float64, skiprow=2)                      
                    
                print('plot_key:', plot_key, 'redshift:', self.redshift, 'name:', self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_id_col_array']['name'+str(self.b)], 'filename:', filename, 'data_offset:', int(self.plot_map_array['data_offset_'+plot_key]))


###   ONLY CHANGE HERE   ######################################################################################################################################
###############################################################################################################################################################                
                #self.plot_map_array['data_offset_'+plot_key]=8

                #for LOAD_FROM_FILE create an array which fits the plot_key+Cataloge+redshif.txt file
                if self.myPipe.i==0 and self.myPipe.a==0:
                    self.my_analysed_data = np.zeros((data_array[:,0].size, int(self.plot_map_array['data_offset_'+plot_key])*(self.myconfig_array['nr_cats'])), dtype=np.double)
                    print(self.my_analysed_data.shape)
                try:
                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a : int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a + int(self.plot_map_array['data_offset_'+plot_key])] = data_array
                except:
                    try:
                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a : int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a + int(self.plot_map_array['data_offset_'+plot_key])] = data_array[:,[0,1,2,2,3,4,4,4,4,4]]
                    except:
                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a : int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a + int(self.plot_map_array['data_offset_'+plot_key])] = data_array[:,[0,1,2,2,3,4,4]]
                if self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_plot_cum']=='True':
                    print('plot cummulative histogramm!')
                    data = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a : int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a + int(self.plot_map_array['data_offset_'+plot_key])]
                    
                    i=1
                    while i<data[:,1].size:
                        #print('i:', i, 'bin:', data[:,1].size-i-1, 'data[data[:,1].size-1-i,1]', data[data[:,1].size-1-i,1], 'data[data[:,1].size-i,1]:', data[data[:,1].size-i,1], '=', end='')
                        data[data[:,1].size-1-i,1]+=data[data[:,1].size-i,1]
                        #print(data[data[:,1].size-1-i,1])
                        i+=1
                    data[:,1]*=(np.log10(data[1,0])-np.log10(data[0,0]))
                    #print('binsize:', (np.log10(data[1,0])-np.log10(data[0,0])))
                    #print(sum(data[:,1]))                   
                    
                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a : int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a + int(self.plot_map_array['data_offset_'+plot_key])]=data

                    
#PREPARE LOG-LOG PLOTTING PREPARE LOG-LOG PLOTTING PREPARE LOG-LOG PLOTTING PREPARE LOG-LOG PLOTTING PREPARE LOG-LOG PLOTTING PREPARE LOG-LOG PLOTTING PREPARE LOG-LOG PLOTTING PREPARE LOG-LOG PLOTTING PREPARE LOG-LOG PLOTTING PREPARE LOG-LOG PLOTTING 

                yerr1=4
                yerr2=5
                add_axis=7
    
                if plot_key=='cSFRD':
                    print('cSFRD plot:', self.myPipe.a, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1)
                    #SHF: cSFRD +51, sumSFR +1, sumsSFR +6
                    
                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a]+=1
                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] = np.log10(self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1])                        
                    #self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+15]                       


                    if self.myconfig_array['catname'+str(self.myPipe.a)].find('SAG_1Gpc_v2')!=-1:
                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1]+= np.log10(0.6777)                       
                        
                else:          
                    print('adjust errorbars', plot_key)

                    self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a]= np.log10(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a])                    
                    if plot_key=='oh2mstar' or plot_key.find('zgas')!=-1:
                        #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1]+=0.25 
                        if self.myconfig_array['catname'+str(self.myPipe.a)].find('Galar')!=-1:
                            self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+4] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+yerr1]
                            self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+5] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+yerr2]
                        else:
                            self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+4] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] + self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+yerr1]
                            self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+5] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] - self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+yerr2]

                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+add_axis] = np.log10(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+add_axis]) - self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a]

                        #if self.myconfig_array['catname'+str(self.myPipe.a)].find('SAG_1Gpc')!=-1:
                         #   self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1]=12+np.log10(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+add_axis]/16)
                    else:
                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] = np.log10(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1])                                    
                        #using MAD
                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+4] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] + (self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+yerr1]/10**self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1])*0.434
                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+5] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] - (self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+yerr2]/10**self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1])*0.434
                        #using percentiles
                        #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+4] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] + ((self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+yerr2]-10**self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1])/10**self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1])*0.434
                        #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+5] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] - ((10**self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1]-self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+yerr1])/10**self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1])*0.434

#                        if self.myconfig_array['catname'+str(self.myPipe.a)].find('SAG_1Gpc_v2')!=-1 or self.myconfig_array['catname'+str(self.myPipe.a)].find('Galacticus')!=-1:
#                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1]-= 3*np.log10(0.6777)  

                    if plot_key=='HMFs':
                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] = np.log10(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+5])                
                    if plot_key.find('mstar2mhalo')!=-1 or plot_key=='oh2mstar' or plot_key=='zgas2mcold':
                        #NFW_con
                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1]=np.log10(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+add_axis])
                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+4] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] + (self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+yerr1]/10**self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1])*0.434
                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+5] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] - (self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+yerr2]/10**self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1])*0.434

                        #mstar/mhalo vs. ssfr plot
                        #-------------------------------------
                        #add-x: mstar
                        #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+add_axis] = np.log10(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+add_axis])                     
                        #add-x: sfr
                        #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+add_axis] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a]+np.log10(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+add_axis])
                        #add-y: sfr
                        #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+6] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+add_axis]

                        #Mhalo
                        #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] = np.log10(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+add_axis])
                        
                        #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] = np.log10(10**self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1]*10**self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a])                    
                    #mhalo2sfr
                    #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1]=np.log10(1.0/(10**self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1]/10**self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+6]))
                    #mhalo2ssfr
                    #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a]-=self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+6]
                
                #print(self.my_analysed_data[:,[int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] ]
                #print(self.my_analysed_data[0:3, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a : int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a + 6]
                #print(self.my_analysed_data
###   ONLY CHANGE HERE   ######################################################################################################################################
###############################################################################################################################################################
         
            else:
                print('########################################################################################################\n#')                                                                                                   
                print('#     PROGRESS STATUS: PLOTXY --> User defined method')                                                                                                   
                print('#\n########################################################################################################')

                def SFH():
                    #++++++++++++++++++++++++++++++++++++++++++++
                    def set_my_env(myenv):
                        env={}
                        if len(myenv)<=1:
                            if myenv=='':
                                for l, item in enumerate(mysamples):
                                    env.update({l: myenv})
                            else:
                                for l, item in enumerate(mysamples):
                                        env.update({l: '_'+myenv})
                        else:
                            for l, item in enumerate(myenv):                                
                                if myenv=='':        
                                    env.update({l: item})
                                else:
                                    env.update({l: '_'+item})
                        return env
                    
                    def set_my_samples_to_plot():
                        keys={}
                        for l, items in enumerate(mysamples):
                            if len(items)>0:
                                item='_'+str(mysamples[l])
                            else:
                                item=str(mysamples[l])
                            keys.update({str(l): item})
                        return keys
        
                    #method [m1, m2, SDSS, demo]
                    method=mymethod
                    #plot_num [1,2,3,'demo']
                    plot_num=myplot_num
                    #prop [mhalo, mstar, sfr, ...]
                    prop=myprop
                    #key [SFH, gr, gr_one]
                    key=mykey
                    #sample [all, red, passive, ...]
                    sample=mysample
                    #catalog [Gal-dens, Gal-dens-corr3, 'Gal400-dens, ...]
                    catalog = mycatalog
                      
                    print('mysamples:', mysamples, r'\n  method:\t', method, r'\n  plot_num:\t', plot_num, r'\n  sample:\t', sample,\
                          r'\n  prop:\t\t', prop, r'\n  key:\t\t', key, r'\n  catalog:\t', catalog, end='')
                     
                    mysamples_to_plot = set_my_samples_to_plot() 
                    
                    myenv='f'
                    myenv=''
                    env               = set_my_env(myenv)
                    print('env:', env)
                    
                    cat_props         = mL.get_SFH_catalog_props(method=method)

                    self.myconfig_array['nr_cats']=len(mysamples_to_plot)
                    #TODO
                    if key=='SFH' or key=='gr' or key=='gr_res' or key=='gr_frac':                    
                        i=0
                        while i<self.myconfig_array['nr_cats']:
                            catname           = cat_props[mycatalog[i]]['output_name']                                     
                            folder            = cat_props[mycatalog[i]]['folder_name']+method  
                            orphan            = cat_props[mycatalog[i]]['legend_add_info'] 
                                              
                            sample_name_props = mL.get_SFH_sample_name_dict(sample, orphan)                    
                            sample_code_space = mL.get_sample_code_space(sample)                                      
        
                            filename          = mycomp+'anaconda/pro/myRun/histos/sfr2z/'+folder                            
                            
                            if i==0:
                                self.my_analysed_data = np.zeros((cat_props[mycatalog[i]]['nSnaps'], int(self.plot_map_array['data_offset_'+plot_key])*self.myconfig_array['nr_cats']), dtype=np.double)
  
                            self.print_redshift = mL.set_SFH_print_redshift(key,
                                                                         catname,
                                                                         orphan,
                                                                         prop[1::],
                                                                         method,
                                                                         myenv=env[i])

                            prop_col_map   = mL.SFH_prop_col_map(sample_name_props[mysamples_to_plot[str(i)][1::]]['catalog_name'])
                            col            = prop_col_map[prop]
                            
                            dt = mL.dt_SFH_ASCII(sample_name_props[mysamples_to_plot[str(i)][1::]]['catalog_name'])
                            
                            file=filename+cat_props[mycatalog[i]]['filename_part2']+sample_code_space+sample+mysamples_to_plot[str(i)]+env[i]+'.txt'                               
                            
                            print('i:', i, 'sample:', mysamples_to_plot[str(i)].ljust(15), 'catalog_name:',  sample_name_props[mysamples_to_plot[str(i)][1::]]['catalog_name'].ljust(15),\
                                  'col:', str(col).ljust(3), 'space:', sample_code_space.ljust(2), 'env:', env[i], end='')
                                                       
                            
                            data = pd.read_csv(file,\
                                               names=[dt[c][0] for c in xrange(len(dt))], skiprows=2, sep=r'\t')
                                
                            data = mL.df_to_sarray(data)
                            #print(np.info(data)                             

                            #Get normalisation of y-axis e.g. then sumSFR is calculated
                            norm_y       = mL.get_norm_y(prop, method, mycatalog[i])
                            print('norm_y:', norm_y)
                            #print(data[[prop[1::], '+1sig'+prop[1::], '-1sig'+prop[1::]]])
                            print('filenname:', file)
                            
                            property_map    = mL.get_property_dict()
                            axis_style      = property_map[prop[1::]]['axis_style']
                            plot_error_bars = property_map[prop[1::]]['plot_error_bars']
                            
                            if key.startswith('gr'):
                                
                                if key=='gr_frac':
                                    k=0
                                    while k<self.my_analysed_data[:,0].size-1:
                                        #delta_prop=prop(z)-prop(z-1)
                                        delta_prop=self.my_analysed_data[k,int(self.plot_map_array['data_offset_'+plot_key])*i+col]-self.my_analysed_data[k+1, int(self.plot_map_array['data_offset_'+plot_key])*i+col]
                                        
                                        #100% is prop(z-1) ..... frac difference% is delta_prop
                                        self.my_analysed_data[k, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=100.0/self.my_analysed_data[k,int(self.plot_map_array['data_offset_'+plot_key])*i+col]*delta_prop
                                        k+=1
                                    self.my_analysed_data[k, int(self.plot_map_array['data_offset_'+plot_key])*i+1]='NaN'
                                    
                                    #print(self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]
                                    
                                else:
                                    print(prop[1::], data[prop[1::]], data[prop[1::]][0])
                                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=data[prop[1::]]/data[prop[1::]][0]

                                if key=='gr_res':
                                    #plot residuals of growth M1 and growth M2
                                    method_res='M2'
                                    if folder.find('400')!=-1:
                                        folder_res='Gal400-dens_main_cents_'+method_res
                                        part2_res='/sfr2z_Galacticus_400Mpc_z_4.15_tarsel_SFH_down3_M1_main_cents'     #Here come M1 because it was a mistake in naming the file during calculation                                  

                                    else:
                                        folder_res='Gal-dens_main_cents_'+method_res
                                        part2_res='/sfr2z_Galacticus_1Gpc_z_4.15_tarsel_SFH_down3_'+method_res+'_main_cents'                                        
                                        
                                    filename_res=mycomp+'anaconda/pro/myRun/histos/sfr2z/'+folder_res
                                    
                                    if i==0:
                                        self.my_analysed_data_res = np.zeros((ncols, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)
            
                                    print('filename_res:', filename_res+part2_res+sample_code_space+sample+'.txt')#, 'sample_code_space:', sample_code_space 
                                    
                                    self.my_analysed_data_res[:, int(self.plot_map_array['data_offset_'+plot_key])*i : int(self.plot_map_array['data_offset_'+plot_key])*i + int(self.plot_map_array['data_offset_'+plot_key])] = myData.readAnyFormat(config=False, mypath=filename_res+part2_res+sample_code_space+sample+band[str(i)]+'.txt', data_format='ASCII', data_shape='shaped', delim=r'\t', mydtype=np.float64, skiprow=2)
            
                                    self.my_analysed_data_res[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=self.my_analysed_data_res[:, int(self.plot_map_array['data_offset_'+plot_key])*i+col]/self.my_analysed_data_res[0, int(self.plot_map_array['data_offset_'+plot_key])*i+col]
                                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]/self.my_analysed_data_res[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]-1                                
                         
                            else:
                                try:
                                    self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] = data[prop[1::]]
                                except:
                                    pass
                                                                 
                                if prop[1::]=='SHMR' and mycatalog[i]=='Gal-dens-corr3':
                                    for item in [prop[1::], '+1sig'+prop[1::], '-1sig'+prop[1::]]:
                                        data[item]=1.0/data[item]
                                        
                                elif prop[1::]=='cgf' and mycatalog[i]=='Gal-dens':
                                    data[prop[1::]]=data['cgf']/data['mstar']
                                    data['+1sig'+prop[1::]]=data['+1sig'+prop[1::]]/data['+1sigmstar']
                                    data['-1sig'+prop[1::]]=data['-1sig'+prop[1::]]/data['-1sigmstar']
                                        
                                print(data[prop[1::]])#, '-1sig'+prop[1::]]
                                self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] = data[prop[1::]]
                                        
                                if prop[1::]=='sumSFR':
                                    norm=data['ndensity']*1e-4*(1000.0/0.6777)**3
                                    self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] = data[prop[1::]]/norm/(1000.0/0.6777)**3
                                    #print(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1]
                                if plot_error_bars==True:
                                    if axis_style=='logc':                                
                                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]+data['+1sig'+prop[1::]]/10**self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]*0.434
                                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]-data['-1sig'+prop[1::]]/10**self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]*0.434                                                                      
                                    else:
                                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4] = data['+1sig'+prop[1::]]
                                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5] = data['-1sig'+prop[1::]]
                                
                                        
                                
       
                            self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+0] = data['z']
                            #print(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1]
                            i+=1

                    elif key=='stats_violin':
                        self.print_redshift='z='+sample
                        self.plot_map_array['data_offset_'+plot_key] = 3
                        import pandas as pd
                        i=0
                        for item in mysamples:
                            print('z:', sample, 'sample:', item, 'prop:', prop, end='')

                            data1 = pd.read_csv(mycomp+'anaconda/pro/data/Galacticus_1Gpc/SFH/SFH_Index_'+str(sample)+'_'+item+'_props_M1.txt', skiprows=2, names=['haloid', 'parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar',  'SHMR', 'mcold',  'Mzgas',  'zcold', 'g-i', 'sfr', 'ssfr'], sep=r'\t')
                            data_M1 = mL.df_to_sarray(data1)
                            data2 = pd.read_csv(mycomp+'anaconda/pro/data/Galacticus_1Gpc/SFH/SFH_Index_'+str(sample)+'_'+item+'_props_M2.txt', skiprows=2, names=['haloid', 'parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar',  'SHMR', 'mcold',  'Mzgas',  'zcold', 'g-i', 'sfr', 'ssfr'], sep=r'\t')
                            data_M2 = mL.df_to_sarray(data2)
                            print('size M1:', data_M1.size, 'size M2:', data_M2.size)
                            if i==0:
                                self.my_analysed_data = np.zeros((60000, int(self.plot_map_array['data_offset_'+plot_key])*len(mysamples)), dtype=np.float64)

                            #self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i]=data_M1['z']
                            if prop=='_g-i' or prop=='_zcold':
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1][0:data_M1.size]=data_M1[prop[1::]]
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+2][0:data_M2.size]=data_M2[prop[1::]]
                                
                            elif prop=='_SHMR':

                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1][0:data_M1.size]=np.log10(1.0/data_M1[prop[1::]])
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+2][0:data_M2.size]=np.log10(1.0/data_M2[prop[1::]])
                                
                            else:
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1][0:data_M1.size]=np.log10(data_M1[prop[1::]])
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+2][0:data_M2.size]=np.log10(data_M2[prop[1::]])

                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1][data_M1.size::]=-99.99                                                       
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+2][data_M2.size::]=-99.99                           
                            i+=1

                    elif key.find('stats')!=-1:
                        self.plot_map_array['data_offset_'+plot_key] = 23                        

                        i=0
                        for sample in mysamples:
                            self.print_redshift=prop_unit_map[prop]
                            if i==0:
                                self.my_analysed_data = np.zeros((ncols, int(self.plot_map_array['data_offset_'+plot_key])*len(mysamples)), dtype=np.float64)
                                #print(np.info(self.my_analysed_data)
                            print('i:', i, 'sample:', sample, 'prop:', prop, 'filenname:', filename+'/stats/Galacticus_1Gpc_SFH_z-evolution_stats'+band[str(i)]+prop+'.txt')                              
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i : int(self.plot_map_array['data_offset_'+plot_key])*i + int(self.plot_map_array['data_offset_'+plot_key])] = myData.readAnyFormat(config=False, mypath=filename+'/stats/Galacticus_1Gpc_SFH_z-evolution_stats'+band[str(i)]+prop+'.txt', data_format='ASCII', data_shape='shaped', delim=r'\t', mydtype=np.float64, skiprow=1)                                
                            #np.info(self.my_analysed_data)
                            if key=='stats_xbar_M1found_M1':
                                col=12
                            elif key=='stats_xbar_M1found_M2':
                                col=11
                            elif key=='stats_xbar_M1_M2':
                                col=2
                            elif key=='stats_xbar_M1-M2_M2':
                                col=21
                            elif key=='stats_xbar_MAD':
                                col=22                                
                            else:
                                col=1
                            print('col:', col)
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+col]                                
                            #print(self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]
                            i+=1
                                                    
                    elif key.find('histo')!=-1:

                        band=mL.get_redshift_histos(catname,
                                                    'histo1')                      
                            

                        self.plot_map_array['data_offset_'+plot_key] = 7
                        
                        self.myconfig_array['nr_cats']=len(band)
                        i=0
                        while i<self.myconfig_array['nr_cats']:
                            if i==0:
                                self.my_analysed_data = np.zeros((25, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)
                                
                            print('i:', i, 'prop:', prop, 'filenname:', filename+'/histos'+part2_to_z+str(band[str(i)])+part2_from_z_to_main+sample_code_space+sample+prop+'_histo.txt')       

                            data=myData.readAnyFormat(config=False, mypath=filename+'/histos'+part2_to_z+str(band[str(i)])+part2_from_z_to_main+sample_code_space+sample+prop+'_histo.txt', data_format='ASCII', data_shape='shaped', delim=r'\t', mydtype=np.float64, skiprow=2)
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] = data[:,6]
                            if prop.find('zcold')!=-1 or prop.find('-')!=-1 or prop.find('log')!=-1 or prop.find('age')!=-1:
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i] =data[:,0]
                            elif prop.find('SHMF')!=-1:
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i] =np.log10(data[:,1])
                            else:
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i] =np.log10(data[:,0])

                            
                            #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] + (data[:,4]/data[:,1])*0.43
                            #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] - (data[:,5]/data[:,1])*0.43
                            #print(self.my_analysed_data[:, [ int(self.plot_map_array['data_offset_'+plot_key])*i, int(self.plot_map_array['data_offset_'+plot_key])*i+1]]        

                            i+=1

                    elif key=='one': 

    
                        self.myconfig_array['nr_cats']=len(band)
                        self.plot_map_array['data_offset_'+plot_key]=7 
                        
                        i=0
                        while i<self.myconfig_array['nr_cats']:
                            if i==0:
                                self.my_analysed_data = np.zeros((ncols, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)
   
                            print('i:', i, band[str(i)])
                            if method.find('stats'):
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+0]=data[:,1]    
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=data[:,band[str(i)]]
                                #self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+4]=data[:,band[str(i)]]-data[:,band[str(i)]+2]
                                #self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+5]=data[:,band[str(i)]+3]+data[:,band[str(i)]]
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+4]=data[:,band[str(i)]]+data[:,band[str(i)]+2]
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+5]=data[:,band[str(i)]]-data[:,band[str(i)]+2]                             
                            else:
                                self.my_analysed_data[0:array.size, int(self.plot_map_array['data_offset_'+plot_key])*i+0]=array['z1']    
                                self.my_analysed_data[0:array.size, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=array[band[str(i)]]
                                
                            #print(self.my_analysed_data[0:10, [int(self.plot_map_array['data_offset_'+plot_key])*i, int(self.plot_map_array['data_offset_'+plot_key])*i+1, int(self.plot_map_array['data_offset_'+plot_key])*i+4, int(self.plot_map_array['data_offset_'+plot_key])*i+5]]
    
                            i+=1
                            
                            
                            
                    elif key=='gr_one': 
                        #growth history for one sample but many parameters

                        band=['_zcold']#, '_mstar', '_mcold', '_Mzgas','_mbh','_Tcons', '_zcold', '_cgf', '_sfr', '_ssfr', 'SFHd']                        
                        
                        self.print_redshift=r'\nhigh-$Z_{cold}$\n'
                            
                        self.myconfig_array['nr_cats']=len(band)

                        for i, sample in enumerate(['','low-zcold','high-zcold']):
                            #TODO: make work for various samples and bands! Until now only one sample can be plotted at a time
                            for a, prop in enumerate(band):
                                
                                if a==0:
                                    #print('i:', i, 'filenname:', filename+part2+sample_code_space+sample+myenv+'.txt')
                                    self.my_analysed_data = np.zeros((ncols, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)
                                    #self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i : int(self.plot_map_array['data_offset_'+plot_key])*i + int(self.plot_map_array['data_offset_'+plot_key])] =\
                                    data=myData.readAnyFormat(config=False, 
                                                             mypath=filename+part2+sample_code_space+sample+myenv+'.txt', 
                                                             data_format='ASCII', 
                                                             data_shape='shaped', 
                                                             delim=r'\t', 
                                                             mydtype=np.float64, 
                                                             skiprow=2)
                                
                                print('prop:', prop, prop_col_map[prop], 'a:', a)
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*a+0]=data[:,0]
                                
                                if prop=='_Tcons':
                                    #calculate Tcons=Mcold/SFR
                                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*a+1]=data[:, prop_col_map[prop]]/data[:, prop_col_map['_sfr']]
                                elif prop=='_cgf':
                                   #calculate cgf=Mcold/Mstar
                                   self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*a+1]=data[:, prop_col_map[prop]]/data[:, prop_col_map['_mstar']]

                                else:
                                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*a+1]=data[:, prop_col_map[prop]]
                                k=0
                                while k<self.my_analysed_data[:,0].size-1:
                                    
                                    prop_before = self.my_analysed_data[k+1, int(self.plot_map_array['data_offset_'+plot_key])*a+1]
                                    prop        = self.my_analysed_data[k,   int(self.plot_map_array['data_offset_'+plot_key])*a+1]
                                    
                                    delta_prop=prop-prop_before
        
                                    f_prop=100.0/prop*delta_prop
                                    
                                    self.my_analysed_data[k, int(self.plot_map_array['data_offset_'+plot_key])*a+1]=f_prop
                                    k+=1
                                    
                                self.my_analysed_data[k, int(self.plot_map_array['data_offset_'+plot_key])*a+1]='NaN'
                                    #self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+6] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+24]
                          
                    elif key=='wp_one': 
                        self.print_redshift='z=0.56'#\n$high-Z_{cold}$'
                        self.plot_map_array['data_offset_'+plot_key] = 6

                        if folder.find('SDSS')!=-1:
                            endfix='_centrals_0.1_200_150'                            
                        else:
                            endfix='_centrals_0.5_150_150'

                        self.myconfig_array['nr_cats']=len(prop)
                        i=0
                        while i<self.myconfig_array['nr_cats']:
                            if i==0:
                                self.my_analysed_data = np.zeros((15, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)
                            if i==12:
                                sample='low-zcold'
                            elif i==1 or i==3:
                                sample='high-zcold'
                                
                            if prop[i].find('Pop')!=-1 or prop[i].find('fila')!=-1 or prop[i].find('knot')!=-1 :
                                print('more wp props ...')
                                endfix='_0.5_150_150'
                                myfilename=filename+'/wp/twoPCF_Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3_more_props_test_'+sample+'_centrals_'+prop[i]+endfix+'_wp.txt'
                                #twoPCF_Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3_more_props_test_low-zcold_centrals_PopB+k_0.5_150_150_wp.txt
                            else:                                
                                print('z=', prop[i])
                                myfilename=filename+'/wp/'+part2_to_z+prop[i]+part2_from_z_to_main+sample_code_space+sample+endfix+'_wp.txt'
                            
                            print(myfilename)
                            data = myData.readAnyFormat(config=False, mypath=myfilename, data_format='ASCII', data_shape='shaped', delim=r'\t', mydtype=np.float64, skiprow=2)
                          
                            self.my_analysed_data[:, [int(self.plot_map_array['data_offset_'+plot_key])*i, int(self.plot_map_array['data_offset_'+plot_key])*i+1, int(self.plot_map_array['data_offset_'+plot_key])*i+4, int(self.plot_map_array['data_offset_'+plot_key])*i+5]]=data[:,[2,3,5,5]]

                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]*=data[:,2]
                            
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=np.log10(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1])
                            i+=1

                    elif key.find('wp')!=-1:
                         
                        print('here clustering!')
                        if folder.find('SDSS')!=-1:
                            endfix='_centrals_0.1_200_150'
                            self.print_redshift+='centrals z='+prop
                        else:
                            endfix='_centrals_0.5_150_150'
                            self.print_redshift='z='+prop
#                        plot2PCF(SFH=True,
#                                 my_custom_filename=filename+'/wp/'+part2_to_z+prop+part2_from_z_to_main+sample_code_space+sample+endfix+'_wp.txt')
                        
                        self.plot_map_array['data_offset_'+plot_key] = 6
                        i=0
                        while i<self.myconfig_array['nr_cats']:
                            if i==0:
                                self.my_analysed_data = np.zeros((15, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)

                            print(filename+'/wp/'+part2_to_z+prop+part2_from_z_to_main+sample_code_space+band[str(i)]+endfix+'_wp.txt')
                            data = myData.readAnyFormat(config=False, mypath=filename+'/wp/'+part2_to_z+prop+part2_from_z_to_main+sample_code_space+band[str(i)]+endfix+'_wp.txt', data_format='ASCII', data_shape='shaped', delim=r'\t', mydtype=np.float64, skiprow=2)
                          
                            self.my_analysed_data[:, [int(self.plot_map_array['data_offset_'+plot_key])*i, int(self.plot_map_array['data_offset_'+plot_key])*i+1, int(self.plot_map_array['data_offset_'+plot_key])*i+4, int(self.plot_map_array['data_offset_'+plot_key])*i+5]]=data[:,[2,3,5,5]]

                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]*=data[:,2]
                            
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=np.log10(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1])
                            
                            
                            #print(self.my_analysed_data[:, [int(self.plot_map_array['data_offset_'+plot_key])*i, int(self.plot_map_array['data_offset_'+plot_key])*i+1]]
                            i+=1
                            
                    elif key.find('xi')!=-1:
                      
                        print('here XI clustering! redshift:', prop)
                        
                        self.print_redshift='z='+prop+r'\nknots'
                                                                            
                        self.plot_map_array['data_offset_'+plot_key] = 6
                        i=0
                        while i<self.myconfig_array['nr_cats']:
                            
                            catname           = cat_props[mycatalog[i]]['output_name']                                     
                            folder            = cat_props[mycatalog[i]]['folder_name']+method
                            run_name          = cat_props[mycatalog[i]]['run_name']
                            folder_2          = cat_props[mycatalog[i]]['folder_name_2']  
                            orphan            = cat_props[mycatalog[i]]['legend_add_info'] 
                            clustering_endix  = cat_props[mycatalog[i]]['clustering_endfix_xi']
                                              
                            sample_name_props = mL.get_SFH_sample_name_dict(mysamples[i], orphan)                    
                            sample_code_space = mL.get_sample_code_space(mysamples[i])   
                            
                            filename          = mycomp+'anaconda/pro/myRun/histos/sfr2z/'+folder 
                            
                            endfix=env[i]+clustering_endix
                            if i==0:
                                self.my_analysed_data = np.zeros((25, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)


                            filename_endfix=folder_2+'/xi/sfr2z_Galacticus_1Gpc_z_'+prop+'_tarsel_SFH_down3_M2_main_cents'+run_name+sample_code_space+mysamples[i]+endfix+'.txt'
                            print(filename_endfix)                                             

                                
                            data = myData.readAnyFormat(config=False, mypath=filename+filename_endfix, data_format='ASCII', data_shape='shaped', delim=r'\t', mydtype=np.float64, skiprow=2)                       
                          
                            self.my_analysed_data[:, [int(self.plot_map_array['data_offset_'+plot_key])*i, int(self.plot_map_array['data_offset_'+plot_key])*i+1, int(self.plot_map_array['data_offset_'+plot_key])*i+4, int(self.plot_map_array['data_offset_'+plot_key])*i+5]]=data[:,[2,3,5,5]]

                            if key.find('r2')!=-1:
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=data[:,2]**2*data[:,3]

                            if key.find('ref')!=-1:
                                
                                if i==0:
                                    ref=data[:,2]**2*data[:,3]
                                #if i<self.plot_map_array['data_offset_'+plot_key]:
                                 #   ref=data[:,int(self.plot_map_array['data_offset_'+plot_key])*(i+1)+1]**2*data[:,int(self.plot_map_array['data_offset_'+plot_key])*(i+1)+3]
                                    #print('test:', ref[0:5]

                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=data[:,2]**2*data[:,3]
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]/ref - 1.0
                                #print(self.my_analysed_data[0:5, int(self.plot_map_array['data_offset_'+plot_key])*i+1]
                            else:
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=np.log10(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1])
                                #self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=np.log10(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1])

                            
                            #print(self.my_analysed_data[:, [int(self.plot_map_array['data_offset_'+plot_key])*i, int(self.plot_map_array['data_offset_'+plot_key])*i+1]]
                            i+=1                            
                            
                    elif key=='envr_props':
                        print('key:', key)
                        self.my_analysed_data = np.zeros((15, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)
                        self.print_redshift='low-$Z_{cold}$'


                def caseSwitcher(myplot):
                
                    choose = {
                        'SFH': SFH
                        }
                        
                    func = choose.get(myplot)
                    return func()
                
                caseSwitcher('SFH')

                        
###############################################################################################################################################################
###############################################################################################################################################################  

        def plotOutput(mycustom_plot_key, mycustom_plot_filename):

            ts = time.time()              
            date_time_stamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H:%M:%S')           
            
            if mycustom_plot_filename!=None:
                plot_filename=self.myconfig_array[self.myconfig_array['catname0']+'_simulation_name']+'_'+self.myconfig_array['catname0']+mycustom_plot_filename
            else:
                plot_filename=self.myconfig_array[self.myconfig_array['catname0']+'_simulation_name']+'_'+self.myconfig_array['catname0']+'_'+plot_key+'_z_'+str(self.redshift)+'_'+str(date_time_stamp)
 
            myOutput.multiPlot(self.my_analysed_data,
                                 loops=self.myconfig_array['nr_cats'],
                                 my_add_subplot= obs_data_array[0],
                                 nr_added_subplots = self.subplot_loops,
                                 legend_added_subplots = obs_data_array[1],
                                 config_path=mycomp+'anaconda/pro/myRun/plot_config/', 
                                 myconfig_datafile=self.plot_map_array[plot_key+'_config'],
                                 mydir=mycomp+'anaconda/pro/myRun/plots/'+plot_key+'/',
                                 myfilename=plot_filename,
                                 data_block_offset=int(self.plot_map_array['data_offset_'+plot_key]),
                                 data_block_subplot_offset=5,
                                 print_redshift=self.print_redshift,
                                 plot_key=plot_key,
                                 custom_plot_key=mycustom_plot_key,
                                 multipanel=self.multipanel)
            
######## MAIN           ######################################################################################################################################   
        self.b=0
        self.analysed_data_map = {}
        self.subplot_is_contour=False
        self.multipanel=False

        #mL.print_props_table_format(from_file=True)
        
        while self.b<self.plot_map_array['nr_plot_keys']:
            
            obs_data_array = [[],[]]
            self.subplot_loops=0

            plot_key = self.plot_map_array['plot_map_id'+str(self.b)]
            
            if plot_key!='sfr2z' and plot_key!='cSFRD':
                nbins=int(self.plot_map_array['nbins_'+plot_key])            
            else:
                nbins=self.myconfig_array['nr_zs']
                
            #rint(self.plot_map_array['data_offset_'+plot_key])
           
            data = np.zeros((nbins, int(self.plot_map_array['data_offset_'+plot_key])*(self.myconfig_array['nr_cats'])), dtype=np.double)
      
            self.myPipe.i=0
            while self.myPipe.i<self.myconfig_array['nr_zs']:

                self.myPipe.a=0
                while self.myPipe.a<self.myconfig_array['nr_cats']:
#                    if self.myPipe.a>0:
#                        time.sleep(1)
                        
                    output_filename_code_space='' 
                    tarsel_code_space=''
                    
                    #create spaces '_' in the filename:
                    if self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_output_filename_code']!='':
                         output_filename_code_space='_'
                    if self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']!='':
                         tarsel_code_space='_'                     
                   
                    if self.load_from_file=='False' and plot_key.find('analyseTargetSelection')==-1: caseSwitcherPlotKey('mainCalculate', None, None, None, None, None, None, None)
                                       
                    elif plot_key.find('analyseTargetSelection')!=-1: caseSwitcherPlotKey('analyseTargetSelection', None, None, None, None, None, None, None)
                     
                    else:
                        if self.myPipe.a==0: LoadObs()
                        
                        if self.myPipe.a>0 and plot_key=='plotXY':                               
                            pass
                        else:
                            if self.plot_custom_loop=='True':
                                
                                ######################################################################################################
                                # choose which method should be plotted m1 (method 1, subsample selected only once at z=0.56), m2 (method 2, subsample
                                # selected at every snapshot)
                                #-------------------------------------------
                                mymethod='M2' 
                                #M1: default only [1] for M2: [1,2] because there are much more subsample and the are divided into two plots
                                #myplot_num ['_1','_2','_3,'demo']  1: standard sample 1 (all, low, high, active, passive, red, blue, redmstar11, 'redmhalo13)
                                #                           2: standard sample 2 (all, mhalo12st, mhalo12gt, mstar10st, mstar11gt, zcold9st, zcold9gt, lowmstar11, lowmhalo13)
                                #                           3: SDSS sample (highZcold-highMstar, highZcold-highMstar, highZcold-highMstar, lowZcold-highMstar', lowZcold-highMstar, lowZcold-highMstar)
                                #                           'demo': demo sample (e.g. for PhD thesis or ESO proposal)                                
                                myplot_num=''

                                # The SFH of which properties should be plotten?
                                #------------------------------------------                                 
                                if self.myconfig_array['catname0'].find('run2')!=-1:
                                    myprop = ['_mhalo', '_mstar', '_SHMF', '_sfr', '_ssfr', '_mcold', '_Mzgas', '_zcold',\
                                              '_g-i', '_r-i', '_cgf', '_mbh', '_rhalfmass', '_rbulgevsrdisk', 'SFHd', 'sumSFR', '_vmax', '_vdisp',\
                                              '_mean_age_stars_disk', '_mean_age_stars_spheroid']
                                    
                                else:
                                    myprop = ['_mhalo', '_mstar', '_SHMF', '_sfr', '_ssfr', '_mcold', '_Mzgas', '_zcold',\
                                              '_g-i', '_r-i', '_cgf', '_mbh', '_rhalfmass', '_rbulgevsrdisk', 'SFHd', 'sumSFR', '_Tcons', '_rbulge', '_rdisk']                                
#TODO
                                myprop = ['_mstar', '_mhalo', '_mbh', '_mcold', '_Mzgas', '_zcold', '_g-i', '_r-i', '_sfr', '_ssfr', '_cSFRD', '_SHMR']
                                
                                #only corr3
                                myprop = ['_cgf', '_Tcons', '_mbar', '_bheff', '_jbar', '_vbulge'] 

                                myprop = ['_cgf', '_SHMR', '_mstar']#, '_Tcons', '_bheff', '_vbulge']
                                myprop = ['_sumSFR']
                                # selected which plot should be created
                                #   default: SFH for each selected property in 'myprop'
                                #   gr:      growth histories of all properties 'in myprop'
                                #   gr_one:  growth histories of each property in 'myprops', but for all subsamples in 'mysamples'
                                #   histo:   histrogram of each property in 'myprops' and all subsamples in 'mysamples'
                                #   wp:      2pCF of each subsample in 'mysamples' and at each redshift defind in 'band'
                                #   wp_one:  2pCF of all subsamples in 'mysamples' together in one plot at each redshift defind in 'band'
                                #   r2xi:
                                #   refxi:
                                #   print_stats:  print(stats of galaxy properts of all subsamples in 'mysamples' at each redshift defind in 'band' usind the function mL.test_methods()                               
                                #------------------------------------------                                
                                workflow = ['default', 'one', 'gr', 'gr_one', 'gr_frac', 'histo','wp','wp_one', 'stats_N', 'stats_xbar_M1-M2_M2', 'stats_xbar_M1found_M1', 'stats_xbar_M1found_M2', 'stats_xbar_MAD']
                                workflow='default'
                                if workflow.find('histo')!=-1:
                                    if self.myconfig_array['catname0'].find('run2')!=-1:                                    
                                        myprop = ['_mhalo', '_mstar', '_sfr', '_ssfr', '_zcold', '_g-i',\
                                                  '_mbh', '_rbulgevsrdisk', '_vmax', '_vdisp',\
                                                  '_mean_age_stars_disk', '_mean_age_stars_spheroid']
                                    else:
                                        myprop = ['_mhalo', '_mstar']#, '_sfr', '_ssfr', '_zcold', '_g-i',\
                                                 # '_mbh', '_rbulgevsrdisk', '_cgf']                                        

                                myprop = myprop[::-1]

                                errors=''


                                if self.myconfig_array['catname0'].find('400Mpc')!=-1:
                                    if workflow=='wp':
                                        band = ['0.55','0.59', '0.64', '0.71', '0.78','0.8', '0.84', '0.89','1.0','1.26','1.56','2.02'][::-1]  
                                    else:
                                        #wp_one
                                        band = ['0.55','0.59', '0.64', '0.67', '0.71','0.74','0.78', '0.8','0.82','0.84','0.89','1.0']
                                        band = ['0.55', '0.71', '1.37','1.56','1.78','2.14','3.57']
                                        
                                else:
                                    # choose redshifts for 2pCF
                                    #-------------------------------------------
                                    band = ['0.09','0.12','0.14','0.17','0.19','0.22','0.25','0.28','0.3','0.33','0.36','0.39','0.43','0.46','0.49','0.52','0.56']
                                
#                                    band = ['0.56','0.59', '0.63', '0.66', '0.7','0.74','0.78','0.82','0.86','0.9','0.94','0.99',\
#                                                 '1.03','1.08','1.12','1.17','1.22','1.27', '1.32','1.37','1.43', '1.48',\
#                                                 '1.54','1.59','1.65','1.71','1.77','1.83','1.9','1.96',\
#                                                 '2.03','2.1','2.16','2.24','2.31','2.38','2.46','2.53','2.61','2.7','2.78','2.86','2.95',\
#                                                 '3.04','3.13','3.22','3.31','3.41','3.51','3.61','3.71','3.82','3.93',\
#                                                 '4.04','4.15']


 
                                    band = ['0.56','0.59', '0.63', '0.66', '0.7','0.74','0.78','0.82','0.86','0.9','0.94','0.99']
                                    #band = ['0.56','0.59', '0.63', '0.7', '0.78','0.86', '0.9', '0.94','1.03','1.22','1.54','2.03']
                                     
                                    #band = ['0.56','0.59','0.74','1.17','2.03','3.04','4.15']
                                    #band = ['0.56', '0.7', '1.37','1.54','1.77','2.1','3.51']
                                    band = ['0.56','1.54'] #y-axis, no x
                                    #band = ['0.7', '1.37'] #no y & no x
                                    #band = ['2.7','2.1'] #no y but x 
                                    #band = ['3.51'] #x- and y-axis

                                   # band = ['0.55','1.56'] #y-axis, no x
                                    #band = ['0.71', '1.36'] #no y & no x
                                    #band = ['1.79','2.14'] #no y but x 
                                    #band = ['3.57'] #x- and y-axis

                                    band = ['0.56','0.59', '0.63', '0.66', '0.70','0.74','0.78','0.82','0.86','0.90','0.94','0.99',\
                                                 '1.03','1.08','1.12','1.17','1.22','1.27', '1.32','1.37','1.43', '1.48',\
                                                 '1.54','1.59','1.65','1.71','1.77','1.83','1.90','1.96',\
                                                 '2.03','2.10','2.16','2.24','2.31','2.38','2.46','2.53','2.61','2.70','2.78','2.86','2.95',\
                                                 '3.04','3.13','3.22','3.31','3.41','3.51','3.61','3.71','3.82','3.93',\
                                                 '4.04','4.15']

                                #band = ['1.37', '1.54', '2.7', '3.04']
                                #band = ['0.7']
                                #band = ['1.37', '2.1', '3.51']   
                                #choose a custom filename prefix
                                #custom_plot_filename_prefix='_SFH_300-r0001'
                                band = ['PopB+k', 'PopB+f', 'PopA+k', 'PopA+f']
                                band = ['knots', 'filaments', 'knots', 'filaments']#, 'knots', 'filaments']
                                #band = ['knots', 'knots', 'knots']
                                #band = ['filaments', 'filaments','filaments']
                                band = ['0.7', '1.37','2.1','3.51']
                                band = ['0.56']
                                custom_plot_filename_prefix='_CMASS_SFH_down3_'
                                for_paper='_test'
                                
                                ######################################################################################################


                                if myplot_num=='_demo':
                                    mysample=['red', 'red', 'red', 'blue', 'blue', 'blue']                                
                                elif mymethod=='M1':
                                    #Galacticus
                                    #mysample=['', 'low', 'high', 'passive', 'active', 'red', 'blue']                                    
                                    mysample=['low', 'low', 'low', 'passive', 'low', 'red', 'low', 'low-zcold', 'high-zcold']
                                    mysample=['low', 'passive', 'red', 'low-zcold', 'high-zcold']

                                    #mysample=['low','low']
                                elif mymethod=='M2' and myplot_num=='':       
                                    #Galacticus
                                    #standard plot1 M2 & M2
                                    mysample=['low', 'low', 'high', 'passive', 'active', 'red', 'blue', 'low-zcold', 'high-zcold', 'mstar11gt', 'zcold9st', 'zcold9gt', 'redmhalo13', 'lowmhalo13']
                                    #mysample=['', 'low', 'high', 'passive', 'active', 'red', 'blue', 'zcold9gt', 'zcold9st']                                     
                #                              'mstar11gt','mhalo12gt', 'mstar10st', 'mhalo12st', 'zcold9gt', 'zcold9st', 'redmstar11', 'redmhalo13', 'lowmstar11', 'lowmhalo13']                                   
                                    #mysample=['', 'red', 'blue', 'low-zcold', 'high-zcold','redmhalo13']
                                    #mysample=['', 'passive', 'red', 'low-zcold', 'high-zcold', 'zcold9st', 'zcold9gt', 'redmhalo13' ]
                                    mysample=['low', 'passive', 'red', 'low-zcold-red', 'high-zcold-red']
                                    #mysample=['', 'low-zcold', 'high-zcold', '', 'low-zcold', 'high-zcold']
                                    
                                    
                                    
                                elif mymethod=='M2' and myplot_num=='_2':       
                                    #Galacticus
                                    #standard plot2 only M2
                                    mysample=['', 'mstar11gt', 'mhalo12gt', 'zcold9gt', 'zcold9st', 'redmstar11', 'redmhalo13', 'lowmstar11', 'lowmhalo13']
                                    #reduced sample for plots in paper
                                    mysample=['', 'mstar11gt', 'zcold9gt', 'zcold9st', 'redmhalo13', 'lowmhalo13']
                                    
                                elif mymethod=='_SDSS':
                                    mysample=['highZcold-highMstar','lowZcold-highMstar']
                                    
                                elif mymethod=='Cholla':
                                    
                                    mysample=['MM_treeID19','LH_treeID19','HH_treeID19']#,'treeID3','treeID4','treeID5','treeID6','treeID7']

                                elif mymethod=='Cholla_stats':
                                    
                                    mysample=['']#,'treeID3','treeID4','treeID5','treeID6','treeID7']
                                      
                                                                        
                                else:
                                    mysample=['','low', 'high', 'passive', 'active', 'red', 'blue']

                                #mysample=['high-zcold-red', 'high-zcold-red', 'high-zcold-red']
                                #mysample=['low-zcold-red', 'low-zcold-red', 'low-zcold-red']
                                #mysample=['', 'low-zcold-red', 'high-zcold-red']
                                #mysample=['', '', '']
                                mycatalog=['Gal-dens', 'Gal-dens',  'Gal-dens', 'Gal-dens-corr3', 'Gal-dens-corr3']
                                #mycatalog=['Gal-dens-corr3', 'Gal-dens-corr3', 'Gal-dens-corr3']

                                if workflow=='zevol':
                                    for sample in ['mstar']:#, 'sfr', 'ssfr']:
                                        for prop in ['_mstar']:
                                            print(sample, prop, mysample, mymethod)
                                            self.myPlot=caseSwitcherPlotKey('loadFromFile', sample, 'zevol', prop, mysample, mymethod, myplot_num, mycatalog)
                                            mycustom_plot_filename='HOD-SF_zevolv_'+sample+'_'+prop+for_paper
                                            plotOutput(prop, mycustom_plot_filename)  

                                elif workflow=='gr' or workflow=='gr_res' or workflow=='gr_frac':
                                    #plot property for all samples in one plot
                                    #for prop in ['rvir1', 'mhalo1']:#, '_mstar', '_SHMF', '_mcold', '_Mzgas' , '_cgf', '_mbh', '_rbulgevsrdisk', '_rhalfmass'][::-1]:
                                    #for prop in ['_mhalo', '_mstar', '_mbh', '_SHMF', '_mcold', '_Mzgas', '_zcold', '_cgf', '_r-i', '_Tcons', '_sfr', '_ssfr', 'SFHd']:#
                                    for prop in myprop:    
                                        self.myPlot=caseSwitcherPlotKey('loadFromFile', '', workflow, prop, mysample, mymethod, myplot_num, mycatalog)
                                        mycustom_plot_filename=custom_plot_filename_prefix+mymethod+prop+'_'+workflow+myplot_num+for_paper
                                        plotOutput('_'+workflow, mycustom_plot_filename)                                       

                                elif workflow.find('one')!=-1:
                                    #3 plot growth of one sample
                                    for sample in ['high-zcold']:#mysample:
                                        self.myPlot=caseSwitcherPlotKey('loadFromFile', sample, workflow, '', mysample, mymethod, myplot_num, mycatalog)
                                        mycustom_plot_filename=custom_plot_filename_prefix+sample+'_'+mymethod+'_'+workflow+for_paper
                                        plotOutput('_'+workflow, mycustom_plot_filename)
                                        
                                elif workflow.find('histo')!=-1:
                                    #4 plot histo of all samples                                    
                                    for sample in mysample:
                                        for prop in myprop:
                                            self.myPlot=caseSwitcherPlotKey('loadFromFile', sample, workflow, prop, mysample, mymethod, myplot_num, mycatalog)
                                            mycustom_plot_filename=custom_plot_filename_prefix+mymethod+'_'+sample+prop+'_'+workflow+for_paper
                                            plotOutput(prop+'_'+workflow, mycustom_plot_filename)

                                elif workflow=='stats_calc':
                                    #5 print(statistics of properties in a text-file
                                    myprops=['mhalo', 'SHMR', 'sfr', 'ssfr', 'mcold', 'Mzgas', 'zcold', 'g-i']
                                    for prop in myprops:                                                            
                                        for sample in mysample:
                                            for count, redshift in enumerate(band):                     
                                                print('count:', count, 'redshift:', redshift, 'sample:', sample, 'prop:', prop)
                                                if count==0:
                                                    stats_data = np.zeros((55, 23), dtype=np.float64)  
                                                myfilename=mycomp+'anaconda/pro/myRun/histos/sfr2z/Gal-dens_main_cents_'+mymethod+'/stats/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_SFH_z-evolution_stats_'+sample+'_'+prop+'.txt'

                                                stats_data, myheader=mL.calc_residuals(stats_data, count, redshift, sample, prop)
                                                                 
                                                myOutput.writeIntoFile(myfilename,
                                                                       stats_data,
                                                                       myheader='SF- project statistics of galaxy properties with redshift; cat: '+self.myconfig_array['catname'+str(self.myPipe.a)]+' CMASS-sample: '+sample+myheader,
                                                                       data_format='%0.3f\t%0.5f\t%0.5f\t%0.5f\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%i\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%i\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%i\t%0.5f\t%0.5f\t%0.5f\t%0.5f')

                                elif workflow=='print_stats':
                                    #5 print(statistics of properties in a text-file
                                    myOutput.writeIntoFile(mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_SFH_z-evolution_bad_mstar_count.txt',
                                               ['#(1) sample name (2) N\n'],                                                                                             
                                               myheader='SF-project z-evolution cat: '+self.myconfig_array['catname'+str(self.myPipe.a)]+' counts how many times mstar<1e6 / sample',
                                               append_mytext=False,
                                               data_is_string=False,
                                               data_format='%s')                                    
                                    myprops=['sfr', 'ssfr', 'mcold', 'Mzgas', 'zcold', 'g-i']
                                    for prop in myprops:                                                            
                                        for sample in mysample:
                                            bad_mstar_count=0
                                            for redshift in band:                     
                                                print('sample:', sample, 'prop:', prop, 'redshift:', redshift)
                                                myfilename=mycomp+'anaconda/pro/myRun/histos/sfr2z/Gal-dens_main_cents_'+mymethod+'/stats/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_SFH_z-evolution_stats_'+sample+'_'+prop+'.txt'

                                                string, string_found, ratio, header_prefix, myheader, bad_mstar_count=mL.test_methods(redshift, sample, prop, bad_mstar_count)
        
                                                if redshift==0.56 or redshift=='0.56':                                                          
                                                    myOutput.writeIntoFile(myfilename,
                                                               [string],
                                                               myheader='SF- project statistics of galaxy properties with redshift; cat: '+self.myconfig_array['catname'+str(self.myPipe.a)]+' CMASS-sample: '+sample+header_prefix+'(1) z (2) frac (N_(M1 found)/N_M1) '+myheader,
                                                               append_mytext=False,
                                                               data_is_string=False,
                                                               data_format='%s')
        
        #                                            myOutput.writeIntoFile(mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_SFH_z-evolution_'+mysample+'_found_M1.txt',
        #                                                       [string_found],
        #                                                       myheader='SF- project z-evolution of properties! cat: '+self.myconfig_array['catname'+str(self.myPipe.a)]+' sample:'+mysample,
        #                                                       append_mytext=False,
        #                                                       data_is_string=False,
        #                                                       data_format='%s')
        #
        #                                            myOutput.writeIntoFile(mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_SFH_z-evolution_'+mysample+'_frac_found.txt',
        #                                                       [str(prop)+r'\t'+str(ratio)],
        #                                                       myheader='SF-project z-evolution cat: '+self.myconfig_array['catname'+str(self.myPipe.a)]+' sample:'+mysample,
        #                                                       append_mytext=False,
        #                                                       data_is_string=False,
        #                                                       data_format='%s')
                                                   
                                                else:
                                                    myOutput.writeIntoFile(myfilename,
                                                               string+r'\n',
                                                               append_mytext=True,
                                                               data_is_string=True,
                                                               data_format='%s')
                                                    
        #                                            myOutput.writeIntoFile(mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_SFH_z-evolution_'+mysample+'_found_M1.txt',
        #                                                       string_found+r'\n',
        #                                                       append_mytext=True,
        #                                                       data_is_string=True,
        #                                                       data_format='%s')                                            
        #
#                                                    myOutput.writeIntoFile(mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_SFH_z-evolution_'+mysample+'_frac_found.txt',
#                                                               str(prop)+r'\t'+str(ratio)+r'\n',
#                                                               append_mytext=True,
#                                                               data_is_string=True,
#                                                               data_format='%s')

 

                                            myOutput.writeIntoFile(mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_SFH_z-evolution_bad_mstar_count.txt',
                                                       str(sample)+r'\t'+str(bad_mstar_count)+r'\n',
                                                       append_mytext=True,
                                                       data_is_string=True,
                                                       data_format='%s')                                           

                                elif workflow.find('stats_violin')!=-1:
                                    #1 plot property for all samples in one plot
                                    for z in band: 
                                        for prop in ['_mhalo', '_mstar', '_SHMR']:# '_zcold', '_SHMR', '_mstar', '_mcold', '_Mzgas', '_sfr', '_ssfr', '_g-i'][::-1]:
                                                self.myPlot=caseSwitcherPlotKey('loadFromFile', z , workflow, prop, mysample, mymethod, myplot_num, mycatalog)
                                                mycustom_plot_filename=custom_plot_filename_prefix+mymethod+prop+'_z_'+z[0]+'_'+workflow+myplot_num+for_paper
                                                plotOutput('_'+workflow+prop, mycustom_plot_filename)

                                elif workflow.find('stats')!=-1:
                                    #1 plot property for all samples in one plot
                                    for prop in ['_mhalo']:#, '_mstar', '_SHMR']:# '_zcold', '_SHMR', '_mstar', '_mcold', '_Mzgas', '_sfr', '_ssfr', '_g-i'][::-1]:
                                            self.myPlot=caseSwitcherPlotKey('loadFromFile', '' , workflow, prop, mysample, mymethod, myplot_num, mycatalog)
                                            mycustom_plot_filename=custom_plot_filename_prefix+mymethod+prop+'_'+workflow+myplot_num+for_paper
                                            plotOutput('_'+workflow, mycustom_plot_filename)                                                
                                                
                                elif workflow=='wp':
                                #4 plot wp of all samples of certain redshift as prop
                                    for prop in band:                                             
                                        self.myPlot=caseSwitcherPlotKey('loadFromFile', '', 'wp', prop, mysample, mymethod, myplot_num, mycatalog)
                                        mycustom_plot_filename=custom_plot_filename_prefix+mymethod+'_z_'+prop+'_wp'+myplot_num+for_paper
                                        plotOutput('_wp', mycustom_plot_filename)

                                elif workflow=='xi' or workflow=='r2xi' or workflow=='refxi':
                                #4 plot wp of all samples of certain redshift as prop
                                    for prop in band:                                             
                                        self.myPlot=caseSwitcherPlotKey('loadFromFile', '', workflow, prop, mysample, mymethod, myplot_num, mycatalog)
                                        mycustom_plot_filename=custom_plot_filename_prefix+mymethod+'_z_'+prop+'_'+workflow+myplot_num+for_paper
                                        plotOutput('_'+workflow, mycustom_plot_filename)
                                       
                                elif workflow=='wp_one':
                                    #4 plot 2pCF of all samples
                                    for sample in mysample:
                                        self.myPlot=caseSwitcherPlotKey('loadFromFile', sample, 'wp_one', band, mysample, mymethod, myplot_num, mycatalog)
                                        mycustom_plot_filename=custom_plot_filename_prefix+mymethod+sample+'_wp_one'
                                        plotOutput('_wp_one', mycustom_plot_filename)
                                elif workflow=='envr_props':
                                    #4 plot 2pCF of all samples
                                    self.myPlot=caseSwitcherPlotKey('loadFromFile', '', 'envr_props', '', mysample, mymethod, myplot_num, mycatalog)
                                    mycustom_plot_filename=custom_plot_filename_prefix+mymethod+'_envr_props'
                                    plotOutput('_envr_props', mycustom_plot_filename)                                        
                                else:
                                    #1 plot property for all samples in one plot
                                    #print('here: 2778\n', myprop
                                    for prop in myprop:
                                        self.myPlot=caseSwitcherPlotKey('loadFromFile', '', 'SFH', prop, mysample, mymethod, myplot_num, mycatalog)
                                        mycustom_plot_filename=custom_plot_filename_prefix+mymethod+prop+myplot_num+errors+for_paper
                                        plotOutput(prop, mycustom_plot_filename)                                   
                                exit()
                                    
                                
                            else:
                                self.myPlot=caseSwitcherPlotKey('loadFromFile', None, None, None, None, None, None, None)                                
                                
                                

                        
                     
                    self.myPipe.a+=1
                                   
                self.myPipe.i+=1

                if self.load_from_file=='True' and self.plot_custom_loop=='False': plotOutput(None, None)
                
            self.b+=1
            
######## END MAIN       ##################################################################################################################################### 
