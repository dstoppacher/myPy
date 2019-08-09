# Load packages
import numpy as np
import arangeData as aD
myData = aD.ArangeData()

import myFuncs as mF
myFuncs = mF.MyFunctions()

import myLib as mL
import dataPipeline as dP
import loadObs as lO
import outputData as oD

myOutput = oD.OutputData(config=False)

import os
import time
system_info = os.getcwd()
start = system_info.find('anaconda')
mycomp = system_info[:start]

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
        
    def sfr2z(
            self,
            plot_key,
            data,
            filename,
            tarsel_code_space,
            output_filename_code_space): 

        self.myPipe.filterAndExtractData(self.mycond_configs,preprocessing_only=True)
        perc_low='-1sig'
        perc_high='+1sig'
        
        mymethod='subsamp'

        if mymethod=='subsamples':       
            #Galacticus
            SFH_analysis_keys=['','low', 'high', 'passive', 'active', 'red', 'blue', 'SFR3st', 'SFR3gt']
        else:
            SFH_analysis_keys=['','low', 'high', 'passive', 'active', 'red', 'blue']
        #LGALAXIES
        #SFH_analysis_keys=['','low', 'high', 'passive', 'active']        
        
        for count, element in enumerate(SFH_analysis_keys):
            print '\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'
            element_space=''        
            #create spaces '_' in the filename:
            if element!='':
                element_space='_'      
            print 'element:', element
            myfilename=filename[0:len(filename)-4]+element_space+element+'.txt'
            #print 'myfilename:',  myfilename
            error_count=0
            if self.myPipe.i>0:
                filename_before = self.mycomp+'anaconda/pro/myRun/histos/'+plot_key+'/'+plot_key+'_'\
                                    +self.myconfig_array['catname'+str(self.myPipe.a)]\
                                    +'_z_'+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i-1)]+'_snapid'+str(self.myPipe.i-1)]['z'])))\
                                    +'_tarsel'+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']\
                                    +output_filename_code_space+str(self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_output_filename_code'])\
                                    +element_space+element+'.txt'
                #print 'filename before:', filename_before
                data = myData.readAnyFormat(config=False, mypath=filename_before, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2) 
    
                
            
            nerrors = self.myPipe.SFR2Z(data,
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
                                      method=mymethod)

            print '\n++++++++++++++++++++++++++++++++++++++++++++++++++'    

            print 'OUTPUT -->',
            myOutput.writeIntoFile(myfilename,
                                   self.myPipe.histo_data_ZvsSFR[:, self.b*int(self.plot_map_array['data_offset_'+plot_key]) : self.b*int(self.plot_map_array['data_offset_'+plot_key])+int(self.plot_map_array['data_offset_'+plot_key])],
                                   myheader= plot_key+' selection key: '+element+' '+self.myconfig_array[self.myconfig_array['catname'+str(self.b)]+'_simulation_name']+' '+self.myconfig_array['catname'+str(self.b)]+' error main progenitor search: '+str(nerrors)+'\n'\
                                            +'(1) z (2) sum(SFR [Msun yr-1]) (3) SFR/SFR(z=0.56) (4) 50th SFR [Msun yr-1] (5) '+perc_low+' SFR [Msun yr-1] (6) '+perc_high+' SFR [Msun yr-1] '\
                                            +'(7) sum(sSFR [yr-1]) (8) sSFR/sSFR(z=0.56) (9) 50th sSFR [yr-1] (10) '+perc_low+' sSFR [yr-1] (11) '+perc_high+' sSFR [yr-1] '\
                                            +'(12) 50th u-i (13) '+perc_low+' u-i (14) '+perc_high+' u-i '\
                                            +'(15) sum(Mstar [Msun]) (16) Mstar/Mstar(z=0.56) (17) 50th Mstar [Msun] (18) '+perc_low+' Mstar [Msun] (19) '+perc_high+' Mstar [Msun] '\
                                            +'(20) sum(Mvir [Msun]) (21) Mvir/Mvir(z=0.56) (22) 50th Mvir [Msun] (23) '+perc_low+' Mvir [Msun] (24) '+perc_high+' Mvir [Msun] '\
                                            +'(25) lookback time(z=0.56) [Gyr] (26) frac centrals [-] (27) frac no-sats [-] (28) frac orphans [-] (29) n x 10-4 [Mpc-3] '\
                                            +'(30) 50th r-i (31) '+perc_low+' r-i (32) '+perc_high+' r-i '\
                                            +'(33) 50th spin [-] (34) '+perc_low+' spin [-] (35) '+perc_high+' spin [-] '\
                                            +'(36) 50th Mstar/Mvir [-] (37) '+perc_low+' Mstar/Mvir  [-] (38) '+perc_high+' Mstar/Mvir  [-] '\
                                            +'(39) sum(Mcold [Msun]) (40) Mcold/Mcold(z=0.56) (41) 50th Mcold [Msun] (42) '+perc_low+' Mcold [Msun] (43) '+perc_high+' Mcold [Msun] '\
                                            +'(44) sum(Mzgas [Msun]) (45) Mzgas/Mzgas(z=0.56) (46) 50th Mzgas [Msun] (47) '+perc_low+' Mzgas [Msun] (48) '+perc_high+' Mzgas [Msun] '\
                                            +'(49) 50th zcold [-] (50) '+perc_low+' zcold [-] (51) '+perc_high+' zcold [-] '\
                                            +'(52) cosmic SFRD [Msun yr-1 Mpc-3]',                                            
                                   data_format="%0.4f\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.5f\t%0.5f\t%0.5f\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.5f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.5e\t%0.5e\t%0.5e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.5e\t%0.5e\t%0.5e\t%0.5e")

        print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'
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
            header='\n(1) '+my_name_x+' ['+self.mycond_configs[my_name_x+'_unit']+'] (2) mstar ['+self.mycond_configs[my_name_y+'_unit']+'] / '+my_name_y+' ['+self.mycond_configs[my_name_y+'_unit']+']\t(3) -dx\t(4) dx\t(5) -dy\t(6) +dy\t(7) N count (8) add axis mstar ['+self.mycond_configs[my_name_y+'_unit']+'] (9) dy add axis'               
        elif mydiv_y2x==True or myadd_axis!=False:           
            header='\n(1) '+my_name_x+' ['+self.mycond_configs[my_name_x+'_unit']+'] (2) '+my_name_y+' ['+self.mycond_configs[my_name_y+'_unit']+']\t(3) -dx\t(4) dx\t(5) -dy\t(6) +dy\t(7) N count (8) add axis '+myadd_axis+' ['+self.mycond_configs[myadd_axis+'_unit']+'] (9) -d add axis (10) +d add axis'
        elif mybinUp2D==True:
            header='\n(1) '+my_name_x+' ['+self.mycond_configs[my_name_x+'_unit']+'] (2) '+my_name_y+' ['+self.mycond_configs[my_name_y+'_unit']+']\t(3) -dx\t(4) dx\t(5) -dy\t(6) +dy\t(7) N count'
        else:
            header='\n(1) '+my_name_x+' ['+self.mycond_configs[my_name_x+'_unit']+'] (2) '+Phi+' ['+myvolume_units+dex+'\t(3) -dx\t(4) dx\t(5) -dy\t(6) +dy\t(7) N count'
      
        myOutput.writeIntoFile(filename,
                               self.myPipe.ratio_data_binAndFrac2D[:, mydata_offset*self.myPipe.a : mydata_offset*self.myPipe.a + mydata_offset],
                               myheader= plot_key+' '+self.myconfig_array['catname'+str(self.myPipe.a)]+' z='+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))+' cumulative: '+info_cum+starforming_cut+header,
                               data_format="%0.8e",
                               mydelimiter='\t')

        return self.myPipe.ratio_data_binAndFrac2D
        

    def matchhalocat(self): 
         
        filteredData, myfilename, myheader, mydata_format, mylogfilename, mylog, sel_col_list = self.myPipe.matchHaloCatWithResultCat(self.mycond_configs)
        
        if filteredData[:,0].size<1e6:                   
            myOutput.writeIntoFile(myfilename,
                                   filteredData[:,sel_col_list],
                                   myheader=myheader,
                                   data_format=mydata_format)

        myOutput.writeIntoFile(mylogfilename,
                               mylog,
                               myheader=myheader,
                               data_format="%s",
                               append_mytext=False,
                               data_is_string=True) 

       
    def MAIN(self):

        def caseSwitcher(plot_key, data, filename):
        
            choose = {
                'SMF':                     MassFunction,
                'HMF':                     HMF,
                'HMF_no':                  HMF,
                'SFRF':                    SFRFunction,
                'sSFRF':                   sSFRFunction,
                'sfr2z':                   SFR2Z,
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
                'matchHaloCat':            MatchHaloCat,
                'plotXY':                  PlotXY,
                'mstar2mhalovsSFR':        Mstar2MhalovsSFR,
                'mstar2rhalf':             Mstar2Rhalf                
                }
                
            func = choose.get(plot_key)
            
            return func(plot_key, data, filename)
       
        def MassFunction(plot_key, data, filename):
            return self.initBinAndFrac2D(plot_key, data, filename, 'mstar', mynormalise=True)

        def HMF(plot_key, data, filename):
            return self.initBinAndFrac2D(plot_key, data, filename, 'mhalo', mycumulative=False, mynormalise=True)

        def SFRFunction(plot_key, data, filename):
            return self.initBinAndFrac2D(plot_key, data, filename, 'sfr', mynormalise=True)
            
        def sSFRFunction(plot_key, data, filename):
            return self.initBinAndFrac2D(plot_key, data, filename, 'ssfr', mynormalise=True)
            
        def SFR2Z(plot_key, data, filename):
            return self.sfr2z(plot_key, data, filename, tarsel_code_space, output_filename_code_space)
            
        def sSFR2Mstar(plot_key, data, filename):
            return self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'sfr', mybinUp2D=True, mystarforming_cut=1e-11)
            
        def OH2Mstar(plot_key, data, filename):
            #mystarforming_cut default=1e-11, z=0.1 better to use the formular 0.3/t_hubble(z=0.1)=3.547e-11 (quisent cut)
            if self.myconfig_array['catname'+str(self.myPipe.a)].startswith('Galacticus'):
                if self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_data_format']=='SAMHDF5':
                    data = self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'mcold', mybinUp2D=True, myadd_axis='zgas_disk', mystarforming_cut=1e-11)
                else:
                    data = self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'mcold', mybinUp2D=True, myadd_axis='Mzgas', mystarforming_cut=1e-11)
                    
            elif self.myconfig_array['catname'+str(self.myPipe.a)].startswith('SAGE'):
                print 'SAGE ...'
                data = self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'mcold_disk', mybinUp2D=True, myadd_axis='Mzgas_disk',mystarforming_cut=1e-11)
            elif self.myconfig_array['catname'+str(self.myPipe.a)].startswith('SAG_'):
                try:
                    print 'use OH_gas_disk ... '
                    data = self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'mcold_disk', mybinUp2D=True, myadd_axis='Mzgas_disk',mystarforming_cut=1e-11) 
                except:
                    print 'NOT available --> use OH_gas_disk_bulge ... '
                    data = self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'mcold', mybinUp2D=True, myadd_axis='Mzgas',mystarforming_cut=1e-11)
            else:
                print 'default!'
                data = self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'mcold', mybinUp2D=True, myadd_axis='Mzgas',mystarforming_cut=1e-11)
            return data

        def Mstar2Rhalf(plot_key, data, filename):
            return self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'rhalf_mass', mybinUp2D=True, mydiv_y2x=False, myadd_axis='rbulge')
            
        def Mstar2Mhalo(plot_key, data, filename):
            return self.initBinAndFrac2D(plot_key, data, filename, 'mhalo', 'mstar', mybinUp2D=True, mydiv_y2x=True, myadd_axis='NFW_con')

        def Mstar2MhalovsSFR(plot_key, data, filename):
            return self.initBinAndFrac2D(plot_key, data, filename, 'sfr', 'mhalo', mybinUp2D=True, mydiv_y2x=False, myadd_axis='mstar')
            
        def Zgas2Mstar(plot_key, data, filename):
            if self.myconfig_array['catname'+str(self.myPipe.a)].startswith('SAGE'): 
                return self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'Mzgas', mybinUp2D=True, myadd_axis='mcold_disk')
            elif self.myconfig_array['catname'+str(self.myPipe.a)].startswith('Galacticus'):
                return self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'zgas_spheroid', mybinUp2D=True, myadd_axis='mcold')
            else:    
                return self.initBinAndFrac2D(plot_key, data, filename, 'mstar', 'Mzgas', mybinUp2D=True, myadd_axis='mcold')

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
            self.myPipe.TwoPCF(filename,self.mycond_configs)

        def HODFunction(plot_key, data, filename):          
            self.myPipe.HODFunction(filename)
                       
        def FilterAndExtractData(plot_key, data, filename):
            self.myPipe.filterAndExtractData(self.mycond_configs)
            
        def AnalyseTargetSelection(plot_key, data, filename):
            return self.myPipe.analyseTargetSelection(data, filename, self.mycond_configs, plot_key)
            
        def MatchHaloCat(plot_key, data, filename):
            self.matchhalocat()
 
           
        def PlotXY(plot_key, data, filename):
            self.myPipe.plotXY(plot_key,filename)
          
        def caseSwitcherPlotKey(plot_key, sample, key, prop):

            def myLoadFromFile():
                LoadFromFile(sample, key, prop)
            
            choose = {
                    'loadFromFile':           myLoadFromFile,
                    'mainCalculate':          MainCalculate,
                    'analyseTargetSelection': TargetSelection              
                    }
                
            func = choose.get(plot_key)
            return func()

        def check_filename():
            
            if self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_use_store_register']=='yes':
                register_path = '/store/erebos/doris/'
            else:
                register_path = self.mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'/'
                
            filename1 = register_path+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_'+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))+'_tarsel'+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']+'.'+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_fileformat']           
            filename2 = register_path+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_'+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']+'.'+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_fileformat']
            
            filename3 = '/store/erebos/doris/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_'+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))+'_tarsel'+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']+'.'+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_fileformat']                
            filename4 = '/store/erebos/doris/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_'+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']+'.'+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_fileformat']
            
            if self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_set_register_name'].find('Skies')!=-1:
                filename5=register_path+'/SkiesANDUniverses/MDPL2_'\
                            +self.myconfig_array['catname'+str(self.myPipe.a)][0:len(self.myconfig_array['catname'+str(self.myPipe.a)])-5] \
                            +'_z_'+str(format(float(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z']), '.2f'))+'.hdf5'
            else:
                filename5 = register_path+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_set_register_name']+'/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_'+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']+'.'+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_fileformat']
            
            filename6 = register_path+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_set_register_name']+'/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_'+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))+'_tarsel'+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']+'.'+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_fileformat']
                            
                            

            #print 'myfilename:', myfilename1, 'skip reading data? --> ', self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_skip_reading_data']
            
            return filename1, filename2, filename3, filename4, filename5, filename6

        def MainCalculate():

            print ' '
            print '########################################################################################################'
            print '#                                                                                                      '
            print '#     PROGRESS STATUS: MAIN CALCULATE:  catname:', self.myconfig_array['catname'+str(self.myPipe.a)], 'b:', self.b, 'i:', self.myPipe.i, 'a:', self.myPipe.a, 'plot_key:', plot_key
            print '#                                                                                                      '
            print '########################################################################################################'
            print ' ' 
                        
            if self.b==0 or plot_key=='sfr2z' or plot_key=='filterData' or plot_key=='plotXY' or plot_key=='twoPCF':
                if self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_data_format']!='HDF5':# and self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_data_format']!='CROSSMATCH':                     
                    self.myPipe.readData()
                else:
                    myfilename1, myfilename2, myfilename3, myfilename4, myfilename5, myfilename6 = check_filename()
                    #uncommand for test reason:
                    #--------------------
                    myfilename=myfilename1
                    #self.myPipe.readData(myfilename=myfilename1)


                    if self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_skip_reading_data']!='yes':
                        if plot_key=='twoPCF':
                            myfilename1 = self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_twoPCF_path_to_data']+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_'+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']+'.'+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_fileformat']                     
                        for filename in [myfilename1, myfilename2, myfilename3, myfilename4, myfilename5, myfilename6]:                            
                            try:
                                print 'try: myfilename:', filename
                                self.myPipe.readData(myfilename=filename)
                                myfilename=filename
                                break
    
                            except:
                                print filename, 'not found ...'
                            
            if plot_key=='plotXY':
                filename=myfilename
            else:
                filename = self.mycomp+'anaconda/pro/myRun/histos/'+plot_key+'/'+plot_key+'_'\
                                +self.myconfig_array['catname'+str(self.myPipe.a)]\
                                +'_z_'+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))\
                                +'_tarsel'+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']\
                                +output_filename_code_space+str(self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_output_filename_code'])\
                                +'.txt'                     
            
            self.my_analysed_data = caseSwitcher(plot_key, data, filename)                
                

            #self.myPipe.showConfig() 

        def TargetSelection():
            
            print ' '
            print '########################################################################################################'
            print '#                                                                                                      '
            print '#     PROGRESS STATUS: MAIN TargetSelection():  catname:', self.myconfig_array['catname'+str(self.myPipe.a)], 'b:', self.b, 'i:', self.myPipe.i, 'a:', self.myPipe.a, 'plot_key:', plot_key
            print '#                                                                                                      '
            print '########################################################################################################'
            print ' '                         

            self.analysed_data_map['name'+str(self.myPipe.a)] = 'analysed_data_array'+str(self.myPipe.a)
            self.analysed_data_map['catname'+str(self.myPipe.a)] = self.myconfig_array['catname'+str(self.myPipe.a)]
            if self.myPipe.a==0:            
                self.analysed_data_map['mytext'] = ''
            
            filename1, filename2, filename3, filename4, filename5, filename6 = check_filename()

            print 'loading ....', 

            for filename in [filename1, filename2, filename3, filename4, filename5, filename6]:
                try:
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
                                                                     filename)
                        #print 'filename:', filename
                    break
                except:
                    print filename, 'not found ...'
                               
            print 'z:', self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z']
            #self.analysed_data_map['mytext']+=self.analysed_data_map['catname'+str(self.myPipe.a)][0:self.analysed_data_map['catname'+str(self.myPipe.a)].find('_')]+' ngal: '+str(len(self.analysed_data_map['name'+str(self.myPipe.a)+'_data']))+'\n'
            self.analysed_data_map['mytext']+=self.analysed_data_map['catname'+str(self.myPipe.a)]+' ngal: '+str(len(self.analysed_data_map['name'+str(self.myPipe.a)+'_data']))+'\n'
            #print self.analysed_data_map
            #if self.myPipe.a==self.myconfig_array['nr_cats']-1: 
            plotAnalysedData(self.analysed_data_map)
                      
        def plotAnalysedData(analysed_data_map):

            print ' '
            print '########################################################################################################'
            print '#                                                                                                      '
            print '#     PROGRESS STATUS: MAIN plotAnalysedData():  catname:', self.myconfig_array['catname'+str(self.myPipe.a)], 'b:', self.b, 'i:', self.myPipe.i, 'a:', self.myPipe.a, 'plot_key:', plot_key
            print '#                                                                                                      '
            print '########################################################################################################'
            print ' '  

    
            LoadObs()
            #print analysed_data_map
            
            
            mycatname=self.myconfig_array['catname'+str(self.myPipe.a)][0:self.myconfig_array['catname'+str(self.myPipe.a)].find('_')]
            myredshift=self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z']
            #print 'mycatname:', mycatname, 'myredshift:', myredshift

            
            mynr_cats=self.myconfig_array['nr_cats']
            myobs={}
            myobs_legend={}
            subplot_loops=self.subplot_loops
            i=0
            a=0           
            while i<self.subplot_loops:
                #print 'i:', i, 'a:', a, obs_data_array[1][str(i)]
                #print obs_data_array[0][str(i)]
                              
                if self.subplot_is_contour[str(i)]==True:
                    #print 'here!', self.subplot_is_contour[str(i)]
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
                
                    #print 'catname'+str(a+1), analysed_data_map['catname'+str(a+1)]                                                                       
                subplot_loops+=1

                myobs.update({str(a+1): obs_data_array[0][str(i)]})               
                myobs_legend.update({str(a+1): obs_data_array[1][str(i)]})
                a+=1
                i+=1

            #print analysed_data_map

   
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
                                 print_redshift='')#pop (ii), centrals')#, centrals')#\n$SFR<-1$')#CMASS DR12:')#$orphans$')# z='+str(float("{0:.2f}".format(myredshift))))#+' centrals')

        def LoadObs():
            #print ' '
            #print 'LOADING OBSERVATIONS ....'
            #print ' '
            if self.plot_map_array['load_obs_'+plot_key]!='False':
                
                key = mL.multicolTestAlgorithm(self.plot_map_array['load_obs_'+plot_key])
                #print 'here:', self.plot_map_array['load_obs_'+plot_key]
                legend={}
                array={}
                self.subplot_is_contour={}
                x=0
                while x<key.size:

                    obs_key = getattr(lO, key[x])                 
                    load_obs = obs_key(load_from_file=False)
                    #print 'load_obs:', load_obs                                        
                    array[str(x)]=load_obs[0]
                    legend[str(x)]=load_obs[1]
                    self.subplot_is_contour.update({str(x): load_obs[2]})
                    x+=1
                               
                obs_data_array[0] = array
                obs_data_array[1] = legend
                self.subplot_loops=x
                

           
        def LoadFromFile(mysample, mykey, myprop):

            self.print_redshift=False
            if self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_use_snapidz_mapping']=='False':
                self.redshift = self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_manual_input_redshift']
            else:
                self.redshift = self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]]['z']
                                           
            if plot_key!='plotXY':
               
                print ' '
                print '########################################################################################################'
                print '#                                                                                                      '
                print '#     PROGRESS STATUS: MAIN LOAD FROM FILE:  catname:', self.myconfig_array['catname'+str(self.myPipe.a)], 'b:', self.b, 'i:', self.myPipe.i, 'a:', self.myPipe.a, 'plot_key:', plot_key
                print '#                                                                                                      '
                print '########################################################################################################'
                print ' '

                    
                if plot_key.find('plotOnly')!=-1:
                    filename=self.mycomp+'anaconda/pro/myRun/histos/analyseTargetSelection/histo_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_'+str(self.redshift)+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']+'_'+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_id_col_array']['name'+str(self.b)]+output_filename_code_space+str(self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_output_filename_code'])+'.txt'                          
                  
                else:                            
                    filename = self.mycomp+'anaconda/pro/myRun/histos/'+plot_key+'/'+plot_key+'_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_'+str(self.redshift)+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_tarsel_code']+output_filename_code_space+str(self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_output_filename_code'])+'.txt'
                
                data_array = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)                      
                    
                print 'plot_key:', plot_key, 'redshift:', self.redshift, 'name:', self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_id_col_array']['name'+str(self.b)], 'filename:', filename, 'data_offset:', int(self.plot_map_array['data_offset_'+plot_key])


###   ONLY CHANGE HERE   ######################################################################################################################################
###############################################################################################################################################################                
                #self.plot_map_array['data_offset_'+plot_key]=8

                #for LOAD_FROM_FILE create an array which fits the plot_key+Cataloge+redshif.txt file
                if self.myPipe.i==0 and self.myPipe.a==0:
                    self.my_analysed_data = np.zeros((data_array[:,0].size, int(self.plot_map_array['data_offset_'+plot_key])*(self.myconfig_array['nr_cats'])), dtype=np.double)

                try:
                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a : int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a + int(self.plot_map_array['data_offset_'+plot_key])] = data_array
                except:
                    try:
                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a : int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a + int(self.plot_map_array['data_offset_'+plot_key])] = data_array[:,[0,1,2,2,3,4,4,4,4,4]]
                    except:
                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a : int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a + int(self.plot_map_array['data_offset_'+plot_key])] = data_array[:,[0,1,2,2,3,4,4]]
                if self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_plot_cum']=='True':
                    print 'plot cummulative histogramm!'
                    data = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a : int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a + int(self.plot_map_array['data_offset_'+plot_key])]
                    
                    i=1
                    while i<data[:,1].size:
                        #print 'i:', i, 'bin:', data[:,1].size-i-1, 'data[data[:,1].size-1-i,1]', data[data[:,1].size-1-i,1], 'data[data[:,1].size-i,1]:', data[data[:,1].size-i,1], '=',
                        data[data[:,1].size-1-i,1]+=data[data[:,1].size-i,1]
                        #print data[data[:,1].size-1-i,1]
                        i+=1
                    data[:,1]*=(np.log10(data[1,0])-np.log10(data[0,0]))
                    #print 'binsize:', (np.log10(data[1,0])-np.log10(data[0,0]))
                    #print sum(data[:,1])                   
                    
                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a : int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a + int(self.plot_map_array['data_offset_'+plot_key])]=data

                    
#PREPARE LOG-LOG PLOTTING PREPARE LOG-LOG PLOTTING PREPARE LOG-LOG PLOTTING PREPARE LOG-LOG PLOTTING PREPARE LOG-LOG PLOTTING PREPARE LOG-LOG PLOTTING PREPARE LOG-LOG PLOTTING PREPARE LOG-LOG PLOTTING PREPARE LOG-LOG PLOTTING PREPARE LOG-LOG PLOTTING 

                yerr1=4
                yerr2=5
                add_axis=7
    
                if plot_key=='sfr2z':
                    print 'sfr2z plot:', self.myPipe.a, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1
                    #SHF: cSFRD +51, sumSFR +1, sumsSFR +6
                    
                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a]+=1
                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] = np.log10(self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1])                        
                    #self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+15]                       


                    if self.myconfig_array['catname'+str(self.myPipe.a)].find('SAG_1Gpc_v2')!=-1:
                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1]+= np.log10(0.6777)                       
                        
                else:          
                    print 'adjust errorbars', plot_key

                    self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a]= np.log10(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a])                    
                    if plot_key=='oh2mstar' or plot_key.find('zgas')!=-1:
                                              
                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+4] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] + self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+yerr1]
                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+5] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] - self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+yerr2]
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
                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+7] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+add_axis]
                        #Mhalo
                        #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] = np.log10(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+add_axis])
                        
                        #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] = np.log10(10**self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1]*10**self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a])                    
                    #mhalo2sfr
                    #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1]=np.log10(1.0/(10**self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1]/10**self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+6]))
                    #mhalo2ssfr
                    #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a]-=self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+6]
                
                print self.my_analysed_data[:,[int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] ]
                #print self.my_analysed_data[0:3, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a : int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a + 6]

###   ONLY CHANGE HERE   ######################################################################################################################################
###############################################################################################################################################################
         
            else:
                print ' '
                print '########################################################################################################'
                print '#                                                                                                      '
                print '#     PROGRESS STATUS: PLOTXY --> User defined method'
                print '#                                                                                                      '
                print '########################################################################################################'
                print ' '  


                def plotOH2mstar():                
                    #OH2mstar test of galaxy properties
                    #------------------------------------------------------------------------------------------------------------------------------
                    
                    self.my_analysed_data = myData.selectData2Compute(self.my_analysed_data, 
                                                                selected_col=1, 
                                                                operator='>', 
                                                                condition=0.0)
                    print self.my_analysed_data.shape
                    self.my_analysed_data = myData.selectData2Compute(self.my_analysed_data, 
                                                                selected_col=1, 
                                                                operator='>', 
                                                                condition=0.0)
                    print self.my_analysed_data.shape    
                    if self.myPipe.a==0:
                        correction_Mzgas2Mcold_array = np.zeros((1, 4), dtype=np.float)
                     
                    #print self.my_analysed_data[:, [0,1,3,4]]
                     
                    offset=self.my_analysed_data[:,0] - self.my_analysed_data[:,1]           
                    print 'offset:', offset
                    print 'mean:', np.mean(offset), 'std:', np.std(offset), 'median:', np.median(offset)              
        
                    
                    correction_Mzgas2Mcold_array[:,0]=self.redshift
                    correction_Mzgas2Mcold_array[:,1]=np.mean(offset)
                    correction_Mzgas2Mcold_array[:,2]=np.std(offset)
                    correction_Mzgas2Mcold_array[:,3]=np.median(offset)
                
                    myOutput.writeIntoFile(self.mycomp+'anaconda/pro/myRun/histos/'+plot_key+'/'+plot_key+'_offset_OH_vs_Mzgas2mcold.txt',
                                           correction_Mzgas2Mcold_array,
                                           myheader= plot_key+' '+self.myconfig_array['catname'+str(self.myPipe.a)]+' starforming: log10 ssfr > -11, convertion of Mzgas/Mcold*0.3655 included\n(1) z\t(2) mean\t(3) std\t(4) median',
                                           mydelimiter='\t',
                                           data_format="%0.6f")
    
    
                    filename = mycomp+'anaconda/pro/myRun/histos/SFRF/SFRF_SAGE_1Gpc_z_0.14_tarsel_full_v2.txt'
                    self.plot_map_array['data_offset_'+plot_key] = 6
                    self.my_analysed_data = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
                    
                    self.my_analysed_data[:,1] = self.my_analysed_data[:,1]*0.6777**3
    
                    self.my_analysed_data[:,3] = self.my_analysed_data[:,3]*0.6777**3
                    self.my_analysed_data[:,4] = self.my_analysed_data[:,4]*0.6777**3             
                    
    
                    
                    myOutput.writeIntoFile(mycomp+'anaconda/pro/myRun/histos/SFRF/SFRF_SAGE_1Gpc_z_0.14_tarsel_full_v2_corrected.txt',
                                           self.my_analysed_data,
                                           myheader= plot_key+' '+self.myconfig_array['catname'+str(self.myPipe.a)]+' z=0.09 cumulative: NO\n(1) mstar [Msun] (2) Phi [Mpc-3 dex-1]	(3) dx	(4) -dy	(5) +dy	(6) N count',
                                           data_format="%0.8e",
                                           mydelimiter='\t')
    
                    filename = mycomp+'anaconda/pro/myRun/histos/SMF/SMF_SAGE_1Gpc_z_0.12_tarsel_full_v2.txt'
                    self.plot_map_array['data_offset_'+plot_key] = 6
                    self.my_analysed_data = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
                    
                    self.my_analysed_data[:,1] = self.my_analysed_data[:,1]*0.6777**3
    
                    self.my_analysed_data[:,3] = self.my_analysed_data[:,3]*0.6777**3
                    self.my_analysed_data[:,4] = self.my_analysed_data[:,4]*0.6777**3                
    
                    
                    myOutput.writeIntoFile(mycomp+'anaconda/pro/myRun/histos/SMF/SMF_SAGE_1Gpc_z_0.12_tarsel_full_v2_corrected.txt',
                                           self.my_analysed_data,
                                           myheader= plot_key+' '+self.myconfig_array['catname'+str(self.myPipe.a)]+' z=0.12 cumulative: NO\n(1) mstar [Msun] (2) Phi [Mpc-3 dex-1]	(3) dx	(4) -dy	(5) +dy	(6) N count',
                                           data_format="%0.8e",
                                           mydelimiter='\t')
    

                def plotLF():
                    #plot LF:
                    #------------------------------------------------------------------------------------------------------------------------------
                    LoadObs()
                    key=0
                    preffix={0: '+Kcorr_047_i_', 1: '+Kcorr_047_i_g-i_', 2: '+Kcorr_047_i_g-i_r-i_', 3: '+Kcorr_047_all_bands_', 4: '+Kcorr_047_all_bands_g-i_', \
                             5: '+Kcorr_047_all_bands_g-i_r-i', 6: '_', 7: 'g-i_', 8: 'g-i_r-i_', 9: 'CMASS_dmesa+sliding_only_g-i_gt_2.35_'}
                    annotation={0: '\nall galaxies\n', 1: '\nall galaxieS\n$g-i > 2.35$', 2: '\nall galaxieS\n$g-i > 2.35$, $r-i>1.0$, $g-r > 1.3$', 3: '\nsliding cut only\n', 4: '\n$d_\perp$+sliding cut only\n', \
                                5: '\nall\n$g-i > 2.35$', 6: '\nCMASS\n$g-i > 2.35$', 7: '\n$d_\perp > 0.55$ cut only\n$g-i > 2.35$', 8: '\nsliding cut only\n$g-i > 2.35$', 9: '\n$d_\perp$+sliding cut only\n$g-i > 2.35$'}
                    version= 'new_mags_'+preffix[key]+'M'
                    band='i'
                    self.plot_map_array['data_offset_'+plot_key] = 7
                    self.redshift=0.56
                    self.print_redshift='z=0.55'+annotation[key]#, '+band+'-band'
                    cat=''
                    i=0
                    while i<self.myconfig_array['nr_cats']:
                        if i==0:
                            self.my_analysed_data = np.zeros((40, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)                    
                            
                        if self.myconfig_array['catname'+str(i)]=='SAGE_1Gpcd':
                            #MD-paper
                            #cat='_350Mpc_run_1235'
                            #CMASS-paper
                            cat='_mags_run_1238'                            
                        elif self.myconfig_array['catname'+str(i)]=='Galacticus_1Gpc': 
                            cat='_RS'
                        else:
                            cat=''
    
                        filename = mycomp+'anaconda/pro/myRun/histos/analyseTargetSelection/histo_'+self.myconfig_array['catname'+str(i)]+'_z_'+str(self.redshift)+'_'+version+'AB_dA_total_'+band+cat+'.txt'
                        print 'i:', i, self.myconfig_array['catname'+str(i)], 'filenname:', filename
                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i : int(self.plot_map_array['data_offset_'+plot_key])*i + int(self.plot_map_array['data_offset_'+plot_key])] = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i]+=0  
                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=np.log10(self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1])
                        i+=1
    
                    #print self.my_analysed_data
                    
                def plotAllBands():   
                    #plot all bands
                    band={'0':'u', '1':'g', '2':'r','3':'i','4': 'z','5':'u', '6':'g', '7':'r','8':'i','9': 'z'}
                    version= 'm'                 
                    self.plot_map_array['data_offset_'+plot_key] = 6
                    self.print_redshift='SAGE 1244, z=0.55, k-corrected mags'
                    #self.print_redshift='z=0.09,Con09,slab\n$M_{AB}$(h^2))'
                    #cat='_500Mpc_run_1134'
                    cat='_1244_mags'    
                    self.myconfig_array['nr_cats']=len(band)
                    h_factor=0.0
                    i=0
                    while i<len(band):
                        if i==0:
                            self.my_analysed_data = np.zeros((80, self.plot_map_array['data_offset_'+plot_key]*len(band)), dtype=np.double)
                        elif i>4:
                            cat='_1244_original_test_mags'
                            h_factor=5*np.log10(0.6777)
                        self.myconfig_array.update({'catname'+str(i):'SAGE_1Gpc'})
                        filename = mycomp+'anaconda/pro/myRun/histos/analyseTargetSelection/histo_'+self.myconfig_array['catname'+str(i)]+'_z_0.56_'+version+'AB_dA_total_'+band[str(i)]+cat+'.txt'
                        print 'i:', i, self.myconfig_array['catname'+str(i)], 'filenname:', filename, band[str(i)], 'h_factor:', h_factor
    
                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i : int(self.plot_map_array['data_offset_'+plot_key])*i + int(self.plot_map_array['data_offset_'+plot_key])] = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
                        #self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i]-=h_factor
                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=np.log10(self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1])
                        #print  self.my_analysed_data[0:10, [int(self.plot_map_array['data_offset_'+plot_key])*i,int(self.plot_map_array['data_offset_'+plot_key])*i + 1]]
                        i+=1                    
                        
                    #print self.my_analysed_data

                def SFH():
                    #++++++++++++++++++++++++++++++++++++++++++++
                    #sfr2z                    
                    #band={'0':'mstar_st_1e9', '1': 'mstar_full', '2':'mstar_1e9_1e10', '3':'mstar_1e10_1e11','4':'mstar_gt_1e11', '5':'no_cuts', '6':'ssfr_conTest', '7':'ssfr_st_-11'}
                    #legend={'0':'$\log_{10}(M_*$ $[M_{\odot}])<9$', '1':'$9<\log_{10}(M_*$ $[M_{\odot}])<10$', '2':'$10<\log_{10}(M_*$ $[M_{\odot}])<11$','3':'$11<\log_{10}(M_*$ $[M_{\odot}])<12$','4': '$\log_{10}(M_*$ $[M_{\odot}])>12$'}
                    #band={'0':'', '1': 'v2_OII_disk', '2':'v2_mcold', '3':'v2_mcold_disk','4':'mstar_gt_1e11', '5':'no_cuts', '6':'ssfr_conTest', '7':'ssfr_st_-11'}
#                    from cosmolopy import cd, fidcosmo
#                    lbt = cd.lookback_time([4.15,4.7,4.85], z0=0.5574, **fidcosmo)/3.1536e+16
#                    print lbt
#                    exit()
                
                    
                    
                    self.plot_map_array['data_offset_'+plot_key] = 52
                    self.redshift=''
                    
                    catname='Gal-dens'
                    self.print_redshift=catname+'\n\n'
                    
                    #prop [mhalo, mstar, sfr, ...]
                    prop=myprop
                    #key [SFH, gr, gr_one]
                    key=mykey
                    #sample [all, red, passive, ...]
                    sample=mysample
                    
                    print 'sample:', sample, 'prop:', prop, 'key:', key   
                    
                    #folder='Gal-dens_main_cents'
                    #folder='Gal-dens_massive'
                    folder='Gal-dens_main_method1'  
                    
                    filename=mycomp+'anaconda/pro/myRun/histos/sfr2z/'+folder
                    if folder.find('main_cents')!=-1:
                        #Gal-dens, main progenitor, centrals only
                        part2='/sfr2z_Galacticus_1Gpc_z_4.15_tarsel_SFH_main_cents'
                        orphan='centrals'
                    elif folder.find('massive')!=-1:
                        #Gal-dens, most massive
                        part2='/sfr2z_Galacticus_1Gpc_z_4.15_tarsel_OII_massive'
                        orphan=''
                    elif folder.find('method1')!=-1:
                        #Gal-dens, main
                        part2='/sfr2z_Galacticus_1Gpc_z_4.15_tarsel_SFH_down3_method1_main'
                        part2_to_z='/sfr2z_Galacticus_1Gpc_z_'
                        part2_from_z_to_main='_tarsel_SFH_down3_method1_main'                        
                        orphan='method 1'                        
                    elif folder.find('method2')!=-1:
                        #Gal-dens, main
                        part2='/sfr2z_Galacticus_1Gpc_run2_z_4.15_tarsel_SFH_down3_method2_main'
                        part2_to_z='/sfr2z_Galacticus_1Gpc_run2_z_'
                        part2_from_z_to_main='_tarsel_SFH_down3_method2_main'
                        orphan='method 2'                    

                    prop_col_map={'_mstar': 16, '_mhalo': 21, '_mcold': 40, '_zcold': 48, '_Mzgas': 45, '_sfr': 3, '_ssfr': 8, '_g-i': 11, '_r-i': 29, '_spin': 32, 'SFH': 1, 'SFHd': 51, '_SHMF': 35, '_logSHMF': 35, 'sumSFR': 1}
 
                    prop_unit_map={'_mstar': '$M_*$', '_mhalo': '$M_{vir}$', '_mcold': '$M_{cold}$', '_zcold': '$Z_{cold}$', '_Mzgas': '$M_{z_{cold}}$', '_sfr': 'SFR', '_ssfr': 'sSFR', '_g-i': '$g-i$', '_r-i': '$r-i$', '_spin': '$S_{halo}$', 'SFH': '', 'SFHd': '', '_SHMF': '$M_{*}/$M_{vir}$', '_logSHMF': '$\log_{10}$ $M_{*}/$M_{vir}$', 'sumSFR': '$\sum$ SFR [$M_{\odot}$ $yr^{-1}$]', 'sumMstar': '$\sum$ $M_{*}$ [$M_{\odot}$]','sumMhalo': '$\sum$ $M_{vir}$ [$M_{\odot}$]'}

                    if sample!='':
                        sample_code_space='_'
                    else:
                        sample_code_space=''

                    if key=='SFH':
                        self.print_redshift=catname+'\n'+orphan+'\n'                        
                    elif sample.startswith('low') or sample.startswith('high'):
                        self.print_redshift=catname+'\n'+orphan+'\n'+sample+' SFR'
                    elif sample.startswith('SFR3st'):
                        self.print_redshift=catname+'\n'+orphan+'\nSFR$<2$'
                    elif sample.startswith('SFR3gt'):
                        self.print_redshift='Gal-dens\n'+orphan+'\nSFR$>2$'                        
                    elif key.find('histo')!=-1:
                        if sample.startswith('low') or sample.startswith('high'):
                            self.print_redshift=catname+'\n'+orphan+'s\n'+sample+' SFR'
                        elif sample=='':
                            self.print_redshift=catname+'\n'+orphan+'\nall'                             
                        else:
                            self.print_redshift=catname+'\n'+orphan+'\n'+sample                        
                    elif prop!='':    
                        self.print_redshift=catname+'\n'+orphan+'\n'+prop_unit_map[prop]
                    
                    else:
                        self.print_redshift=catname+'\n'+orphan+'\n'+sample                            
                        if sample=='':
                            sample_code_space=''
                            self.print_redshift=catname+'\n'+orphan+'\nall'                          
                    if folder.find('method1')!=-1:
                        band={'0':'', '1': '_low', '2':'_high', '3':'_passive','4':'_active', '5':'_red', '6':'_blue'}
                    else:
                        band={'0':'', '1': '_low', '2':'_high', '3':'_passive','4':'_active', '5':'_red', '6':'_blue', '7': '_SFR3st', '8': '_SFR3gt'}
                    self.print_redshift='\n'+'\n$M_{*}$'+'\n$M_{vir}$\nmethod 2 (samples selected at z=0.56)'
                    self.print_redshift='\n'+'\n$M_{*}$'+'\n$M_{vir}$\nmethod 1 (samples selected at each z)'
                    band={'0': '_passive' , '1':''}
                    band={'0': '_red', '1': '_red', '2':'_blue', '3':'_blue','4':'', '5':''}#, '6':'_blue', '7': '_SFR3st', '8': '_SFR3gt'}
                    self.myconfig_array['nr_cats']=len(band)
                    
                    if key=='SFH' or key=='gr':                    
                        i=0
                        while i<self.myconfig_array['nr_cats']:
                            if i==0:
                                self.my_analysed_data = np.zeros((55, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)
                            
                            #filename2=filename+part2+sample_code_space+sample+band[str(i)]+'.txt'
                            
                            col=prop_col_map[prop]
                            print 'i:', i, 'col:', col,  'filenname:', filename+part2+sample_code_space+sample+band[str(i)]+'.txt'
                            if i==1 or i==3 or i==5:
                                col=prop_col_map['_mhalo']
                                
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i : int(self.plot_map_array['data_offset_'+plot_key])*i + int(self.plot_map_array['data_offset_'+plot_key])] = myData.readAnyFormat(config=False, mypath=filename+part2+sample_code_space+sample+band[str(i)]+'.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)

                            if prop=='sumSFR':
                                norm_y=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+28]*(1450.0**3)*1e-4
                            else:
                                norm_y=1
                                                     
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+0]+=0
                            #self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+0]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+24]
                            if key.startswith('gr'):
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+col]/self.my_analysed_data[0, int(self.plot_map_array['data_offset_'+plot_key])*i+col]
                            else:
                                if prop=='_logSHMF':
                                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=np.log10(self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+col]/norm_y)
                                    self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]+self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+col+1]/10**self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]*0.434
                                    self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]-self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+col+2]/10**self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]*0.434                                                                      

                                else:
                                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+col]/norm_y
                                    
                            if prop!='SFHd' and prop!='_logSHMF':
                                self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+col+1]
                                self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+col+2]                                                                      
                            self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+6]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+24]
                            i+=1
#                            print self.my_analysed_data[:, [0,1]]
#                            print norm_y
                            
                                                    
                    elif key.find('histo')!=-1:

                        #band_z_mag={'0': 0.56, '1': 0.74, '2': 0.99, '3': 1.48, '4': 2.03, '5': 2.53, '6': 3.04, '7': 3.51, '8': 4.04}
                        if key=='histo1':
                            band_z_mag={'0': 0.59, '1': 0.74, '2': 0.86, '3': 0.9, '4': 1.03}
                        else:
                            band_z_mag={'0': 1.22, '1': 1.54, '2': 2.03, '3': 3.04, '4': 4.15}                            

                        self.plot_map_array['data_offset_'+plot_key] = 7
                        
                        self.myconfig_array['nr_cats']=len(band_z_mag)
                        i=0
                        while i<self.myconfig_array['nr_cats']:
                            if i==0:
                                self.my_analysed_data = np.zeros((25, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)
                                
                            print 'i:', i, 'filenname:', filename+part2_to_z+str(band_z_mag[str(i)])+part2_from_z_to_main+sample_code_space+sample+prop+'_histo.txt'       

                            data=myData.readAnyFormat(config=False, mypath=filename+part2_to_z+str(band_z_mag[str(i)])+part2_from_z_to_main+sample_code_space+sample+prop+'_histo.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] = data[:,6]
                            if prop.find('zcold')!=-1 or prop.find('-')!=-1:
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i] =data[:,0]
                            else:
                                
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i] =np.log10(data[:,0])

                            
                            #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] + (data[:,4]/data[:,1])*0.43
                            #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] - (data[:,5]/data[:,1])*0.43
                            #print self.my_analysed_data[:, [0,1]]        

                            i+=1                           
                    elif key=='gr_one': 
                        #growth history for one sample but many parameters
                        #mstar, mhalo, mcold, Mzgas
                        band={'0': 16, '1': 21, '2': 40, '3': 45}#, '4': 11, '5': 29}
    
                        self.myconfig_array['nr_cats']=len(band)
                        i=0
                        while i<self.myconfig_array['nr_cats']:
                            if i==0:
                                self.my_analysed_data = np.zeros((55, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)
    
                            print 'i:', i, 'filenname:', filename+part2+sample_code_space+sample+'.txt'
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i : int(self.plot_map_array['data_offset_'+plot_key])*i + int(self.plot_map_array['data_offset_'+plot_key])] = myData.readAnyFormat(config=False, mypath=filename+part2+sample_code_space+sample+'.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
        
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+0]+=0
                            #self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+0]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+24]
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+band[str(i)]]/self.my_analysed_data[0, int(self.plot_map_array['data_offset_'+plot_key])*i+band[str(i)]]
                            #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+col+1]
                            #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+col+2]                                                                      
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+6] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+24]
                            i+=1                        





                def plotManyCurves():  

#                    prefix='new_mags'
#                    endfix='i'
#                    lum='M'
##                    band={0:'-Kcorr_047_i_', 1:'-Kcorr_047_i_g-i_2.35_', 2:'-Kcorr_047_i_g-i_2.15_', 3:'-Kcorr_047_i_g-i_2.0_', 4:'-Kcorr_047_i_g-i_1.75_', 5:'-Kcorr_047_i_g-i_1.5_', 6:'-Kcorr_047_i_g-i_1.25_', 7:'-Kcorr_047_i_g-i_1.0_'}
##                    band={0:'+Kcorr_047_i_', 1:'+Kcorr_047_i_g-i_2.35_', 2:'+Kcorr_047_i_g-i_2.15_'}#, 3:'+Kcorr_047_i_g-i_2.0_', 4:'+Kcorr_047_i_g-i_1.75_', 5:'+Kcorr_047_i_g-i_1.5_', 6:'+Kcorr_047_i_g-i_1.25_', 7:'+Kcorr_047_i_g-i_1.0_'}
##                    band={0:'+Kcorr_047_all_bands_', 1:'+Kcorr_047_all_bands_g-i_2.35_', 2:'+Kcorr_047_all_bands_g-i_2.15_', 3:'+Kcorr_047_all_bands_g-i_2.0_', 4:'+Kcorr_047_all_bands_g-i_1.75_', 5:'+Kcorr_047_all_bands_g-i_1.5_', 6:'+Kcorr_047_all_bands_g-i_1.25_', 7:'+Kcorr_047_all_bands_g-i_1.0_'}#, 8:'+Kcorr_047_all_bands_g-i_0.75_'}#, 9:'+Kcorr_047_all_bands_g-i_0.5_'}
##                    band={0:prefix+'', 1:prefix+'g-i_2.35_', 2:prefix+'g-i_2.15_', 3:prefix+'g-i_2.0_'}#, 4:'g-i_1.75_', 5:'g-i_1.5_', 6:'g-i_1.25_', 7:'g-i_1.0_'}
##                    band={0:prefix+'g-i_2.35_', 1:prefix+'density_sample_g-i_2.35_', 2:prefix+'mass_sample_g-i_2.15_'}#, 4:'g-i_1.75_', 5:'g-i_1.5_', 6:'g-i_1.25_', 7:'g-i_1.0_'}
#                    #band={0:prefix+'test_', 1:prefix+'density_sample_test_', 2:prefix+'mass_sample_test_'}#, 4:'g-i_1.75_', 5:'g-i_1.5_', 6:'g-i_1.25_', 7:'g-i_1.0_'}
#                    #band={0:prefix+'CMASS_', 1:prefix+'', 2:prefix+'CMASS_g-i_2.35_', 3:prefix+'CMASS_g-i_st_2.35_'}#,  3:prefix+'CMASS_mass_sample_g-i_2.35_', 1, 4:prefix+'g-i_2.35_'}#, 6:'g-i_1.25_', 7:'g-i_1.0_'}
#                    #band={0:prefix+'_test_CMASS_g-i_gt_2.35_', 1:prefix+'_test_CMASS_', 2:prefix+'_CMASS_density_sample_g-i_gt_2.35_', 3:prefix+'_CMASS_density_sample_', 4:prefix+'_CMASS_mass_sample_g-i_gt_2.35_', 5:prefix+'_CMASS_mass_sample_'}
#                    #band={0:prefix+'_CMASS_g-i_gt_2.35_', 1:prefix+'_test2_', 2:prefix+'_test2_g-i_gt_2.35_'}
#                    #band={0:prefix+'_CMASS_g-i_gt_2.35_', 1:prefix+'_CMASS_density_sample_g-i_gt_2.35_', 2:prefix+'_CMASS_mass_sample_g-i_gt_2.35_'}
#                    band={0:prefix+'_CMASS_', 1:prefix+'_CMASS_down_sample3_', 2:prefix+'_CMASS_mass_sample_'}
#                    #band={0:'+Kcorr_047_i_', 1:'+Kcorr_i_after_g-i_2.35_', 2:'+Kcorr_i_after_g-i_2.15_', 3:'+Kcorr_i_after_g-i_2.0_', 4:'+Kcorr_i_after_g-i_1.75_', 5:'+Kcorr_i_after_g-i_1.5_', 6:'+Kcorr_i_after_g-i_1.25_', 7:'+Kcorr_i_after_g-i_1.0_'}
#                    #band={0:'+Kcorr_047_i_', 1:'CMASS_+Kcorr_047_i_', 2:'+Kcorr_047_i_g-i_2.15_', 3:'CMASS_+Kcorr_047_i_g-i_2.15_'}
#                    #band={0:'mags_+Kcorr_047_all_bands_', 1:'CMASS_-Kcorr_047_all_bands_', 2: 'CMASS_-Kcorr_047_all_bands_g-i_2.35_', 3: 'CMASS_-Kcorr_047_all_bands_g-i_2.15_', 4: 'CMASS_-Kcorr_047_all_bands_g-i_2.0_'}                 
#                    version= 'g-i>2.35'
#                    version= 'Guo+13'
#                    self.redshift=0.56
#                    self.print_redshift='z=0.55, Kcorr: i-band mAB-0.47'#, Galacticus CMASS K-corr +0.47, i-band (after g-i cut)'#+annotation[key]#, '+band+'-band'
#                    #version= endfix+'-band'
#                    self.print_redshift=version
#                    cat='_Guo13'
#                    #cat='_g-i_gt_2.35'
#
#                    self.plot_map_array['data_offset_'+plot_key] = 7                    
#
#                    i=0
#                    while i<len(band):
#                        if i==0:
#                            self.my_analysed_data = np.zeros((30, self.plot_map_array['data_offset_'+plot_key]*len(band)), dtype=np.double)
#                        self.myconfig_array.update({'catname'+str(i):self.myconfig_array['catname'+str(self.myPipe.a)]})
#                        filename = mycomp+'anaconda/pro/myRun/histos/analyseTargetSelection/histo_Galacticus_1Gpc'+'_z_'+str(self.redshift)+'_'+band[i]+lum+'AB_dA_total_'+endfix+cat+'.txt'  
#                        #filename = mycomp+'anaconda/pro/myRun/histos/analyseTargetSelection/histo_'+self.myconfig_array['catname'+str(i)]+'_z_'+str(self.redshift)+'_'+prefix+'MAB_dA_total_'+band[i]+cat+'.txt'  
#
#                        print 'i:', i, 'filenname:', filename, band[i]
#                        data=myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
#                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i] = data[:,0]
#                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] = np.log10(data[:,1])
#                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+2] = data[:,2]
#                        
#                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] + (data[:,4]/data[:,1])*0.43
#                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] - (data[:,5]/data[:,1])*0.43
#    
#                        i+=1                      
#                    self.myconfig_array['nr_cats']=len(band)                    


                    #++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    #oh2mstar
                    #legend={'0':'Galacticus', '1':'SAGv2 OH disk', '2':'SAGv2 mcold','3':'SAGv2 mcold disk','4': 'SAGv3 OH', '5':'SAGv3 mcold', '6':'SAGv3 mcold disk', '7':'SAGE zgas disk','8':'SAGE mcold disk'}
                    #band={'0':'Galacticus_1Gpc_z_0.09', '1': 'SAG_1Gpc_v2_z_0.07_OII_disk_sfr+h', '2':'SAG_1Gpc_v2_z_0.07_mcold_sfr+h', '3':'SAG_1Gpc_v2_z_0.07_mcold_disk_sfr+h','4':'SAG_1Gpc_z_0.09_OII', '5':'SAG_1Gpc_z_0.09_mcold', '6':'SAG_1Gpc_z_0.09_mcold_disk', '7':'SAGE_1Gpc_z_0.09_OII', '8':'SAGE_1Gpc_z_0.09_mcold_disk' }


#                    self.myconfig_array['nr_cats']=9
#                    i=0
#                    while i<self.myconfig_array['nr_cats']:
#                        if i==0:
#                            self.my_analysed_data = np.zeros((60, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)
#                        self.myconfig_array.update({'catname'+str(i):'Galacticus_1Gpc'})
#                        filename = mycomp+'anaconda/pro/myRun/histos/oh2mstar/oh2mstar_'+band[str(i)]+'.txt'   
#
#                        print 'i:', i, 'filenname:', filename, band[str(i)]
#                        
#                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i : int(self.plot_map_array['data_offset_'+plot_key])*i + int(self.plot_map_array['data_offset_'+plot_key])] = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
#                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i]=np.log10(self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i])
##                        if band[str(i)].find('OII')!=-1:
##                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]+=0.25
#                        i+=1

                        
                    #++++++++++++++++++++++++++++++++++++++++++++
                    #SMF CMASS and HOD
                    key='SMF'
                    key2='scd'
                    
                    if key=='SMF' or key.find('SFR')!=-1:
                        prefix='SMF_Galacticus_1Gpc_z_0.56_new_mags'
                        #prefix='SMF_SAG_1Gpc_z_0.56_new'
                        
                        #prefix='SMF_SAGE_1Gpc_z_0.56_tarsel_v3_mags_run_1238'
                        prefix=key+'_LGALAXIES_500Mpc_z_0.56_tarsel'                        
                        prefix=key+'_Galacticus_1Gpc_z_'
                        prefix='SMF_SAG_1Gpc_v2_z_0.07'
                        #prefix=key+'_Galacticus_1Gpc_z_0.56_tarsel_new_mags'
                        #band={'0':'Galacticus_1Gpc_z_0.56_mstar'}#, '1': 'SAG_1Gpc_z_0.56_mstar', '2':'SAGE_1Gpc_z_0.56_1238_mstar'}#, '3':'SAG_1Gpc_v2_z_0.07_mcold_disk_sfr+h','4':'SAG_1Gpc_z_0.09_OII', '5':'SAG_1Gpc_z_0.09_mcold', '6':'SAG_1Gpc_z_0.09_mcold_disk', '7':'SAGE_1Gpc_z_0.09_OII', '8':'SAGE_1Gpc_z_0.09_mcold_disk' }
                        #band={'0': prefix+'_CMASS', '2': prefix+'_CMASS_density_sample', '3': prefix+'_CMASS_mass_sample', '1': prefix}#, '3':'SAG_1Gpc_v2_z_0.07_mcold_disk_sfr+h','4':'SAG_1Gpc_z_0.09_OII', '5':'SAG_1Gpc_z_0.09_mcold', '6':'SAG_1Gpc_z_0.09_mcold_disk', '7':'SAGE_1Gpc_z_0.09_OII', '8':'SAGE_1Gpc_z_0.09_mcold_disk' }
                        #band={'0':'SMF_SAGE_1Gpc_z_0.56_CMASS_sample', '1': 'SAGE_1Gpc_z_0.56_mstar', '2':'SMF_SAGE_1Gpc_z_0.56_densCut_CMASS_mstar'}
                        #band={'0': prefix+'_CMASS', '2': prefix+'_test_no5h_CMASS', '3': prefix+'_test_nozboost_no5h_CMASS', '1': prefix+'_test_nozboost_no5h_noiltdmesa_CMASS'}#, '3':'SAG_1Gpc_v2_z_0.07_mcold_disk_sfr+h','4':'SAG_1Gpc_z_0.09_OII', '5':'SAG_1Gpc_z_0.09_mcold', '6':'SAG_1Gpc_z_0.09_mcold_disk', '7':'SAGE_1Gpc_z_0.09_OII', '8':'SAGE_1Gpc_z_0.09_mcold_disk' }
    
                        #band={'0': prefix+'_centrals', '1': prefix+'_sats', '2': prefix+'_CMASS_centrals', '3': prefix+'_CMASS_sats', '4': prefix+'_CMASS_density_sample_no', '5': prefix+'_CMASS_density_sample_sats'}#, '4': prefix+'_mass_sample_centrals', '5': prefix+'_mass_sample_sats'}
                        band={'0': prefix+'_CMASS', '2': prefix+'_CMASS_down_sample3', '3': prefix+'_CMASS_mass_sample', '1': prefix}
                        #band={'0': prefix+'_CMASS_down_sample3', '1': prefix+'_CMASS_density_sample', '2': prefix+'_CMASS_cross_sample', '3': prefix}                        
                        #band={'0': prefix+'', '1': prefix+'_g-i_gt_2.35', '2': prefix+'_g-i_gt_2.15', '3': prefix+'_g-i_gt_2.0'}
                        
                        #band={'0': prefix+'_CMASS_down_sample', '1': prefix+'_35bins_2e10', '2': prefix+'_CMASS_down_sample_no', '3': prefix+'_35bins_2e10_no'}
                        
                        band={'0': prefix+'_CMASS_down_sample_test', \
                              '1': prefix+'_35bins_2e10', \
                              '2': prefix+'_CMASS_density_mstar', \
                              '3': prefix+'_CMASS_density_vmax'}#, \
                              #'4': prefix+'_density_vmax', \
                              #'5': prefix+'_density_vmax_sat'}
                              
                        band={'0': prefix+'_OII', \
                              '1': prefix+'_test_OII'}#, \
                              #'2': prefix+'_test_OII'}#, \
#                              '4': prefix+'_density_vmax', \
#                              '5': prefix+'_density_vmax_sat'

#                        band={'0': prefix+'_CMASS_down_sample_test', \
#                              '1': prefix+'_CMASS_density_vmax', \
#                              '2': prefix+'_CMASS_density_mstar'}#, \
#                              '4': prefix+'_density_vmax', \
#                              '5': prefix+'_density_vmax_sat'
                        
                        
                    elif key2=='sc':
                        prefix='HOD_Galacticus_1Gpc_z_0.56_tarsel_new_mags'
                        #prefix='HOD_SAG_1Gpc_z_0.56_tarsel_new'
                        band={'0': prefix+'_CMASS', '1': prefix+'_CMASS_down_sample3', '2': prefix+'_CMASS_mass_sample'}#, '2':prefix+'_CMASS_down_sample', '1':prefix+'_CMASS_density_sample'}
                        #band={'0': prefix+'_CMASS_down_sample', '1': prefix+'_CMASS_mass_sample'}

                    else:
                        #HOD for each sample seperatly                   
                        prefix=key+'_Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample_HOD_mhalo_cents_200c'
                        
                        sample='dens $M_{*}$'# f$^{rand}_{sats}$'
                        #sample='dens $V_{max}$'
                        prefix=key+'_LGALAXIES_500Mpc_z_0.56_tarsel_CMASS_density_mstar_no_mhalo'
                        #prefix=key+'_SAG_1Gpc_z_0.56_new_CMASS_down_sample_mhalo_cents'
                        band={'0': prefix+'_all', '1': prefix+'_centrals', '2': prefix+'_sats'}


                    endfix='_mstar'
                    #green butterfly test
                    band={'0': 'SMF_Galacticus_1Gpc_z_0.56_new_mags_CMASS', '1': prefix+'_sample2'+endfix, '2': prefix+endfix, '3': prefix+'_sample2_centrals'+endfix, '4': prefix+'_centrals'+endfix, '5': prefix+'_sample2_orphans'+endfix, '6': prefix+'_orphans'+endfix, '7': prefix+'_sample2_no-sats'+endfix, '8': prefix+'_no-sats'+endfix}
                 
                    if self.myconfig_array['catname'+str(self.myPipe.a)].startswith('SAG_'):                       
                        band={'1':'SMF_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_0.56', '0':'SMF_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_0.56_densCut_CMASS_mstar'}
                    else:
                        band={'0':'SMF_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_0.56_CMASS', '1': 'SMF_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_0.56', '2':'SMF_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_0.56_densCut_CMASS_mstar'}#, '3':'SMF_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_0.56_densCut_CMASS_mstar'}

                    self.print_redshift='Galacticus, z=0.56'                    
                    self.plot_map_array['data_offset_'+plot_key] = 7                    

                    i=0
                    while i<len(band):

                        self.myconfig_array.update({'catname'+str(i):self.myconfig_array['catname'+str(self.myPipe.a)]})
                        if key2=='sc':
                            filename = mycomp+'anaconda/pro/myRun/histos/'+key+'/'+band[str(i)]+'_mbasic_200c_sc_centrals.txt' 
                        else:
                            filename = mycomp+'anaconda/pro/myRun/histos/'+key+'/'+band[str(i)]+'.txt'  

                        print 'i:', i, 'filenname:', filename#, band[str(i)]
                        data=myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
                        if i==0:
                            self.my_analysed_data = np.zeros((data[:,0].size, self.plot_map_array['data_offset_'+plot_key]*len(band)), dtype=np.float32)
                            
                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i] = np.log10(data[:,0])#+0.03925
                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] = np.log10(data[:,1])
                        
                        if key=='HOD' and key2!='sc':
                            #self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] = np.log10((data[:,5]-1)/data[:,2])
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+4] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] + (1.0/np.sqrt(data[:,6]))/data[:,1]*0.434
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+5] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] - (1.0/np.sqrt(data[:,6]))/data[:,1]*0.434                                     
                        else:
                            if key2=='sc':
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] = data[:,1]
                                
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+4] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] + data[:,4]
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+5] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] - data[:,5]              
                            else:
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+4] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] + data[:,4]/data[:,1]*0.434
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+5] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] - data[:,5]/data[:,1]*0.434                  
                                
                                
                        #print self.my_analysed_data[:, [int(self.plot_map_array['data_offset_'+plot_key])*i+0, int(self.plot_map_array['data_offset_'+plot_key])*i+1]]
                        i+=1  

#                    self.print_redshift='z=0.56'#(>)\nnon-orphan sats'
#                    #self.print_redshift='LG-'+sample
#                    self.myconfig_array['nr_cats']=len(band)
##                    #print self.my_analysed_data[:, 0]
#                    #SFR-OII-paper instant starformation analysis
##                    band={'0':'gas', '1':'gas_16', '2': 'gas_16', '3':'apprx_gas', '4':'apprx_gas+star'}
##                    corr={'0': 0, '1': 0, '2': 0.25, '3': 0, '4': 0}
##                    band={'0':'sfr_inst', '1':'sfr_spheroid_inst', '2': 'sfr_quies_inst', '3':'sfr_inst_spheroid+quies'}
#                    
#                    self.plot_map_array['data_offset_'+plot_key] = 8                    
#                    self.myconfig_array['nr_cats']=1#len(band)
#                    i=0
#                    while i<self.myconfig_array['nr_cats']:
#                        if i==0:
#                            self.my_analysed_data = np.zeros((50, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)
#                        self.myconfig_array.update({'catname'+str(i):self.myconfig_array['catname'+str(self.myPipe.a)]})
#                        filename = mycomp+'anaconda/pro/myRun/histos/plotXY/histo_SAG_1Gpc_z_0.94_tarsel_OII_'+band[str(i)]+'.txt'   
#
#                        print 'i:', i, 'filenname:', filename, band[str(i)]
#                        data=myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
#                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i] = np.log10(data[:,0]) #+ corr[str(i)]
#                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] = np.log10(data[:,1])
#                        
#                        i+=1  


                    #++++++++++++++++++++++++++++++++++++++++++++
                    #Mstar vs Rhalf_mass CMASS
#                    prefix='mstar2rhalf_Galacticus_1Gpc_z_0.56_tarsel_new_mags'
#                    band={'0':'Galacticus_1Gpc_z_0.56_mstar'}#, '1': 'SAG_1Gpc_z_0.56_mstar', '2':'SAGE_1Gpc_z_0.56_1238_mstar'}#, '3':'SAG_1Gpc_v2_z_0.07_mcold_disk_sfr+h','4':'SAG_1Gpc_z_0.09_OII', '5':'SAG_1Gpc_z_0.09_mcold', '6':'SAG_1Gpc_z_0.09_mcold_disk', '7':'SAGE_1Gpc_z_0.09_OII', '8':'SAGE_1Gpc_z_0.09_mcold_disk' }
#                    band={'0': prefix+'_CMASS', '1': prefix+'_CMASS_density_sample', '2': prefix+''}#, '3': prefix+''}#, '3':'SAG_1Gpc_v2_z_0.07_mcold_disk_sfr+h','4':'SAG_1Gpc_z_0.09_OII', '5':'SAG_1Gpc_z_0.09_mcold', '6':'SAG_1Gpc_z_0.09_mcold_disk', '7':'SAGE_1Gpc_z_0.09_OII', '8':'SAGE_1Gpc_z_0.09_mcold_disk' }
#                    #band={'0':'SMF_SAGE_1Gpc_z_0.56_CMASS_sample', '1': 'SAGE_1Gpc_z_0.56_mstar', '2':'SMF_SAGE_1Gpc_z_0.56_densCut_CMASS_mstar'}
#                 
##                    if self.myconfig_array['catname'+str(self.myPipe.a)].startswith('SAG_'):                       
##                        band={'1':'SMF_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_0.56', '0':'SMF_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_0.56_densCut_CMASS_mstar'}
##                    else:
##                        band={'0':'SMF_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_0.56_CMASS', '1': 'SMF_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_0.56', '2':'SMF_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_0.56_densCut_CMASS_mstar'}#, '3':'SMF_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_0.56_densCut_CMASS_mstar'}
#
#                    self.print_redshift='Galacticus, z=0.56'                    
#                    self.plot_map_array['data_offset_'+plot_key] = 10                    
#                    self.myconfig_array['nr_cats']=len(band)
#                    i=0
#                    while i<self.myconfig_array['nr_cats']:
#                        if i==0:
#                            self.my_analysed_data = np.zeros((50, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)
#                        self.myconfig_array.update({'catname'+str(i):self.myconfig_array['catname'+str(self.myPipe.a)]})
#                        filename = mycomp+'anaconda/pro/myRun/histos/mstar2rhalf/'+band[str(i)]+'.txt'   
#
#                        print 'i:', i, 'filenname:', filename, band[str(i)]
#                        data=myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
#                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i] = np.log10(data[:,0])
#                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] = np.log10(data[:,1]*1000)
#                        
#                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] - (data[:,4]/data[:,1])*0.434
#                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] + (data[:,5]/data[:,1])*0.434
#
#                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+7] = np.log10(data[:,7]*1000)
#
#    
#                        i+=1  

                def plot2PCF():

                    #plot 2PCF:
                    #------------------------------------------------------------------------------------------------------------------------------
                  
                    LoadObs()
                    cut=''
                    gtype_general='_no'
                    #gtype_general='_down_sample3_no_kroup_mstar_gt_'
                    #gtype_general='_down_sample3_no_'                    
                    #gtype_general= sample+'_sample_test'+gtype_general
                    seltype_general='all'
                    plot='wp'
    
                    mypanels=['wp', 'wp_ref']
                    scale='LSS'
                    calc_prog='CMASS'
                              
                    if seltype_general=='mstar':
                        name='$M_{_*}$'
                    elif seltype_general=='mcold':
                        name='$M_{Cold}$'
                    elif seltype_general=='sfr':
                        name='$SFR$'
                    elif seltype_general=='M_r':
                        name='$M_r$'
                    else:
                        name=''
    
                    self.my_analysed_data={}
                    if len(mypanels)>1:
                        self.multipanel=True
                    else:
                        self.multipanel=False
                    plot_legend=True              
                    add_obs=True

                    #rmin          [0.05, 0.5,  0.1, 1,  10,  50,  75,  100]
                    #rmax          [1.5,  3,    200, 10, 100, 100, 200, 200]
                                        #mstar cuts    [1.5e11, 1.73e11, 2e11, 2.24e11, 3e11, 3.5e11, 4e11, 5e11, min-half, half-max]
                    #mstar cuts new [1.5e11, 1.73e11, 2e11, 2.24e11, 3e11, 3.5e11, 1.05e11, 3.45e11]
                    #mstar bins    [1.5e11-2e11, 1.73e11-3e11, 2e11-3e11, 2.24e11-4e11, 3e11-4e11, 3.5e11-5e11, 1.05e11-3.45e11, 3.45e11-1e13]

                    #new pi_max:   [20, 40, 80, 150]
                    #mstar cuts new [1.5e11, 1.73e11, 2e11, 2.24e11, 3e11 3.5e11, 4e11, 5e11]

                    #new mass bins                          
                    #cut_list=[1.5e11, 1.73e11, 2e11, 2.24e11, 3e11, 3.5e11, 4e11]
                    #cut_list_max_dict={1.5e11: 1.73e11, 1.73e11: 2e11, 2e11: 2.24e11, 2.24e11: 3e11, 3e11: 3.5e11, 3.5e11: 4e11, 4e11: 5e11}                     
                    
                    #small scales 
                    prefix=gtype_general
                    prefix=gtype_general+'_mstar_gt_10.0'
                    pimax='_150'
                    endfix=prefix+'_0.5_150'+pimax
                    #endfix=''

                    mstar_cut=[11.25, 11.35, 11.45, 11.55, 11.65]
                    #mstar_cut=['bin1', 'bin2', 'bin3']
                    band={}
                    #band={0: prefix+'1.5e+11'+endfix, 1:prefix+'1.73e+11'+endfix, 2:prefix+'2e+11'+endfix, 3: prefix+'2.24e+11'+endfix, 4: prefix+'3e+11'+endfix, 5:prefix+'3.5e+11'+endfix, 6:prefix+'4e+11'+endfix}#, 7:prefix+'5e+11'+endfix} 
                    if gtype_general.find('no')==-1 and gtype_general.find('centrals')==-1 and gtype_general.find('all')==-1:
                        mypanels=['wp']
                        endfix='0.5_150'+pimax
                        for i, masses in enumerate(mstar_cut):
                            band.update({i:gtype_general+str(mstar_cut[i])+'_mstar_'+endfix}) 
                    
                    elif (prefix.find('_mstar_st_')!=-1 or prefix.find('mstar_gt')!=-1 or prefix.find('_mstar_sm_')!=-1 or str(mstar_cut[0]).find('bin')!=-1) and prefix.find('_10.0')==-1:                         
                        mypanels=['wp']
                        endfix='_0.5_150'+pimax
                        print gtype_general
                        for i, masses in enumerate(mstar_cut):                                                                                                                             
                            band.update({i:gtype_general+str(mstar_cut[i])+endfix})

                            
                    else:
                        #mypanels=['wp']
                        #band={0: prefix+'11.25_11.36'+endfix, 1:prefix+'11.36_11.48'+endfix, 2: prefix+'11.48_15'+endfix} 
                        #band={0: prefix+'1.5e+11_1.73e+11'+endfix, 1:prefix+'1.73e+11_2e+11'+endfix, 2: prefix+'2e+11_2.24e+11'+endfix, 3: prefix+'2.24e+11_3e+11'+endfix, 4:prefix+'3e+11_3.5e+11'+endfix, 5:prefix+'3.5e+11_4e+11'+endfix}#, 6:prefix+'4e+11_5+e11'+endfix} 
                        #band={1: prefix+'1.5e+11_2e+11'+endfix, 2:prefix+'1.73e+11_3e+11'+endfix, 0:prefix+'1.05e+11_3.45e+11'+endfix}#, 6:prefix+'4e+11_max'+endfix} 
                        #band={0: '_down3_sample'+endfix, 1:'_down_sample3_PS'+endfix, 2: '_mass_sample_all'+endfix}
                        #band={0: '_color_sample'+endfix, 1:'_down_sample3'+endfix, 2: '_mass_sample'+endfix}
                        
                        band={0: '_down_sample'+endfix, 1:'_density_sfr_lowest_no'+endfix, 2: '_density_ssfr_lowest_no'+endfix} 

                    #print band

                    self.myconfig_array['nr_cats']=len(band)

                    count=0
                    for f in mypanels:
                        print 'f:', f
                        self.print_redshift=name+' '+cut+' '+'$'+gtype_general+'$'
                        i=0
                        while i<len(band):
                            
                            gtype=gtype_general
#                            if self.myconfig_array['catname'+str(i)]=='SAGE_1Gpc' and gtype=='all': gtype='non-orphans'
                            seltype=seltype_general
#                            if self.myconfig_array['catname'+str(i)]=='SAGE_1Gpc' and seltype=='mcold': seltype='mcold_disk'                  
                            print 'i:', i, self.myconfig_array['catname'+str(0)], 'z:', self.redshift, 'plot:', plot, 'calc_prog:', calc_prog, gtype
        
                            if calc_prog=='Cute':                    
                                #filename example: all-Galacticus_1Gpc_z_0.09_CUT1_Contreras+13_mcold.2PCF-SSS
                                filename =mycomp+'anaconda/pro/myRun/histos/twoPCF/CUTE/'+cut+'/'+gtype+'-'+self.myconfig_array['catname'+str(i)]+'_z_'+str(self.redshift)+'_'+cut+'_Contreras+13_'+seltype+'.2PCF-'+scale 
                                mydata = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', nr_col=4, nr_rows=54, data_shape='shaped', delim=' ', mydtype=np.float64, skiprow=0)
                            
                            elif calc_prog=='Corrfunc':
                                filename =mycomp+'anaconda/pro/myRun/histos/twoPCF/'+cut+'/twoPCF_'+self.myconfig_array['catname'+str(i)]+'_z_'+str(self.redshift)+'_tarsel_'+cut+'_Contreras+13_'+seltype+'_'+gtype+'_'+plot[0:2]+'.txt'
                                mydata = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', nr_col=6, nr_rows=int(self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_twoPCF_nbins']), data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
                                mydata[:,[0,1]]=mydata[:,[2,3]]
                                #print mydata[:,[0,1]]
                             
                            elif calc_prog=='CMASS':
                                #self.print_redshift='Galacticus non-orphan satellites, z=0.557'
                                
                                self.print_redshift='Galacticusz=0.557\n'
                                self.print_redshift=''

                                if endfix.find('all')!=-1 or gtype_general.find('all')!=-1:
                                     self.print_redshift+='$all$\n'
                                elif endfix.find('centrals')!=-1 or gtype_general.find('centrals')!=-1:
                                     self.print_redshift+='$centrals$\n'
                                elif endfix.find('no-sats')!=-1:
                                     self.print_redshift+='$sats$\n'
                                elif endfix.find('sats')!=-1 or gtype_general.find('sats')!=-1:
                                     self.print_redshift+='$sats$\n'                                     
                                elif endfix.find('no')!=-1 or gtype_general.find('no')!=-1:
                                     self.print_redshift+='$centrals+sats$\n'
                                else:
                                    self.print_redshift+='\n'
                              
                                #pi_max=', $\pi_{max}$: 60.0 Mpc'
                                if endfix.find('gt_-22.0')!=-1:
                                    self.print_redshift+='$M_r>-22$\n'

                                elif endfix.find('gt_-21.0')!=-1:
                                    self.print_redshift+='$M_r>-21$\n'
                                    
                                elif endfix.find('gt_-21.5')!=-1:
                                    self.print_redshift+='$M_r>-21.5$\n'
                                    
                                elif endfix.find('22.5_Mr_-21.5')!=-1:
                                    self.print_redshift+='$-22.5<M_r<-21.5$\n'
                                    
                                elif endfix.find('22.0_Mr_-21.0')!=-1:
                                    self.print_redshift+='$-22<M_r<-21$\n'

                                elif endfix.find('21_Mr_-22')!=-1:
                                    self.print_redshift+='$-22<M_r<-21$\n'                                    
                                    
                                elif seltype_general.find('23.0_Mr_-22.0')!=-1:
                                    self.print_redshift+='$-23<M_r<-22$'
                                      
                                elif seltype_general.find('23.5_Mr_-22.5')!=-1:
                                    self.print_redshift+='$-23.5<M_r<-22.5$'
                                    
                                elif seltype_general.find('Mr_st_-22')!=-1:
                                    self.print_redshift+='$M_r<-22$'
                                     
                                elif seltype_general.find('Mr_st_-23')!=-1:
                                    self.print_redshift+='$M_r<-23$'

                                elif endfix.find('RS')!=-1:
                                     self.print_redshift+='$g-i>2.35$\n'
                                        
                                elif endfix.find('BC')!=-1:
                                     self.print_redshift+='$g-i<2.35$\n'
                                     
                                elif band[0].find('mstar')!=-1 and gtype_general.find('density')!=-1:
                                     #self.print_redshift+='Gal-dens\n'
                                     self.print_redshift+='SAG-dens\n'
                                     
                                elif band[0].find('mstar')!=-1 and gtype_general.find('mass')!=-1:
                                     self.print_redshift+='Gal-mass\n'
                                     
                                elif band[0].find('mstar')!=-1 and gtype_general.find('color')!=-1:
                                     self.print_redshift+='Gal-cols\n'
                                     
                                elif (band[0].find('mstar')!=-1 or band[0].find('bin')!=-1) and gtype_general.find('down')!=-1:
                                    if gtype_general.find('down_sample3')!=-1:
                                        #self.print_redshift+='Guo+13\n'
                                        self.print_redshift+='Gal-dens\n'
                                    else:
                                        #self.print_redshift+='g-i>2.35\n'
                                        self.print_redshift+='SAG-dens PS\n'
                                else:
                                    self.print_redshift+='\n'                                                                          

                                if endfix.find('kroup')!=-1 or gtype_general.find('kroup')!=-1 or band[0].find('kroup')!=-1:
                                    self.print_redshift+='\n'#Kroupa IMF\n'#fixed $n_z$'
                                else:
                                    self.print_redshift+='\n'  
                                    
#                                elif key=='20-19':
#                                    if len(mypanels)==10:
#                                        if self.myconfig_array['catname'+str(i)].find('Galacticus')!=-1:
#                                            self.print_redshift+='z=0.1\n$-20<M_r<-19$\nf$_{cents}$: 53.9 %\nf$_{stats}$: 17.8 %\nf$_{orphs}$: 28.3 %'
#                                        elif self.myconfig_array['catname'+str(i)].find('SAG_')!=-1:
#                                            self.print_redshift+='z=0.1\n$-20<M_r<-19$\nf$_{cents}$: 79.0 %\nf$_{stats}$: 16.2 %\nf$_{orphs}$: 4.8 %'
#                                        elif self.myconfig_array['catname'+str(i)].find('SAGE')!=-1:
#                                            self.print_redshift+='z=0.1\n$-20<M_r<-19$\n$f_{cents}$: 53.9 %\n$f_{stats}$: 17.8 %\n '
#                                else:
#                                    self.print_redshift+=', z=0.557'
                                
#                                if pimax.find('_200')!=-1:
#                                     self.print_redshift+='$pi_{max}=200$'
#                                elif pimax.find('_150')!=-1:
#                                    self.print_redshift+=r'${\mathrm{\pi_{max}}=150}$'                                    
#                                elif pimax.find('_100')!=-1:
#                                    self.print_redshift+='$pi_{max}=100$'
#                                elif pimax.find('_80')!=-1:
#                                    self.print_redshift+=r'${\mathrm{\pi_{max}}=80}$'
#                                elif pimax.find('_60')!=-1:
#                                    self.print_redshift+=r'${\mathrm{\pi_{max}}=60}$'
#                                elif pimax.find('_40')!=-1:
#                                    self.print_redshift+=r'${\mathrm{\pi_{max}}=40}$'
#                                elif pimax.find('_20')!=-1:
#                                    self.print_redshift+=', $pi_{max}=20$'                                    
#                                elif pimax.find('_10')!=-1:
#                                    self.print_redshift+=', $pi_{max}=10$'
#                                elif pimax.find('_5')!=-1:
#                                    self.print_redshift+=', $pi_{max}=5$'
#
#                                
#                                self.print_redshift+=' Mpc'
    
                                #self.print_redshift='z=0.55, $\log_{10}$ $n_z=-3.47$'
                                #self.print_redshift='z=0.55, $\log_{10}$ $M_{*}>10.7$'
                                #self.print_redshift='z=0.55, $\log_{10}$ $M_{*}>11.1$'
                                #MD-paper                                 
                                #filename =mycomp+'anaconda/pro/myRun/histos/plotXY/twoPCF/twoPCF_'+self.myconfig_array['catname'+str(i)]+'_z_'+str(self.redshift)+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(i)]+'_tarsel_code']+'_'+plot[0:2]+'.txt'
                                #Galacticus CMASS Clustering
                                if band[i].find('HAM')!=-1:
                                    filename =mycomp+'anaconda/pro/myRun/histos/twoPCF/Galacticus_CMASS/twoPCF_Z_HAM_z_0.56_BigMD_LC_'+plot[0:2]+'.txt'
                                else:
                                    #filename =mycomp+'anaconda/pro/myRun/histos/twoPCF/Galacticus_CMASS/twoPCF_Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS'+band[i]+'_'+gtype_general+seltype_general+'_'+plot[0:2]+'.txt'
                                    
                                    
                                    #filename =mycomp+'anaconda/pro/myRun/histos/twoPCF/Galacticus_CMASS/twoPCF_Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS'+band[i]+'_'+gtype_general+seltype_general+'_'+plot[0:2]+'.txt'
#                                    if i<2:
                                    model='LGALAXIES'
                                    boxsize='500Mpc'
                                    #filename =mycomp+'anaconda/pro/myRun/histos/twoPCF/Galacticus_CMASS/twoPCF_Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS'+band[i]+'_'+plot[0:2]+'.txt'
                                    filename =mycomp+'anaconda/pro/myRun/histos/twoPCF/'+model+'/twoPCF_'+model+'_'+boxsize+'_z_0.56_tarsel_CMASS'+band[i]+'_'+plot[0:2]+'.txt'

#                                    else:
                                    #filename =mycomp+'anaconda/pro/myRun/histos/twoPCF/SAG_CMASS/twoPCF_SAG_1Gpc_z_0.56_tarsel_new_CMASS'+band[i]+'_'+plot[0:2]+'.txt'


                                data = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', nr_col=6, nr_rows=int(self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_twoPCF_nbins']), data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
                                mydata=np.zeros((data[:,0].size,7), dtype=np.float32)
                                mydata[:,[0,1,4,5]]=data[:,[2,3,5,5]]
                                #print data
                                if plot=='wp':
                                    mydata[:,1]*=data[:,2]
                                    mydata[:,4]=mydata[:,1]+data[:,5]*data[:,2]
                                    mydata[:,5]=mydata[:,1]-data[:,5]*data[:,2]
                                    #mydata[:,4]=mydata[:,1]+mydata[:,4]
                                    #mydata[:,5]=np.log10(mydata[:,1])-(abs(mydata[:,1]-mydata[:,5]))#/mydata[:,1])*0.434

                                if band[i].find('HAM')!=-1 and plot=='wp':
                                    mydata[:,[0,1,4,5]]=data[:,[0,1,4,5]]                                    
                                  
                                #print mydata[:,[0,1,4,5]]
                                #print np.log10(mydata[:,1])
                                #exit()                                
                              
                            elif f=='MD':
                                filename =mycomp+'anaconda/pro/myRun/histos/twoPCF/twoPCF_'+self.myconfig_array['catname'+str(i)]+'_z_'+str(self.redshift)+'_tarsel'+tarsel_code_space+self.myconfig_array[self.myconfig_array['catname'+str(i)]+'_tarsel_code']+output_filename_code_space+str(self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_output_filename_code'])+'_'+plot[0:2]+'.txt'
    
                                mydata = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', nr_col=6, nr_rows=int(self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_twoPCF_nbins']), data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
                                mydata[:,[0,1]]=mydata[:,[2,3]]
                                if self.myconfig_array['catname'+str(i)].find('Galacticus')!=-1:
                                    self.print_redshift='Galacticus                (a)\n'
                                elif self.myconfig_array['catname'+str(i)].find('SAG_')!=-1:
                                    self.print_redshift='SAG                         (b)\n'                                    
                                elif self.myconfig_array['catname'+str(i)].find('SAGE')!=-1:
                                    self.print_redshift='SAGE                       (c)\n'
                                    
                                key='22-21'
                                #pi_max=', $\pi_{max}$: 60.0 Mpc'
                                if key=='22-21':
                                    self.print_redshift='z=0.1, $-22<M_r<-21$                                         (d)'
                                elif key=='21-20':
                                    self.print_redshift='z=0.1, $-21<M_r<-20$                                         (c)'
                                elif key=='20-19':
                                    if len(mypanels)==10:
                                        if self.myconfig_array['catname'+str(i)].find('Galacticus')!=-1:
                                            self.print_redshift+='z=0.1\n$-20<M_r<-19$\nf$_{cents}$: 53.9 %\nf$_{stats}$: 17.8 %\nf$_{orphs}$: 28.3 %'
                                        elif self.myconfig_array['catname'+str(i)].find('SAG_')!=-1:
                                            self.print_redshift+='z=0.1\n$-20<M_r<-19$\nf$_{cents}$: 79.0 %\nf$_{stats}$: 16.2 %\nf$_{orphs}$: 4.8 %'
                                        elif self.myconfig_array['catname'+str(i)].find('SAGE')!=-1:
                                            self.print_redshift+='z=0.1\n$-20<M_r<-19$\n$f_{cents}$: 53.9 %\n$f_{stats}$: 17.8 %\n '
                                    else:
                                        self.print_redshift='z=0.1, $-20<M_r<-19$                                        (b)'  
                                elif key=='19-18':
                                    self.print_redshift='z=0.1, $-19<M_r<-18$                                        (a)'
                                elif key=='18-17':
                                    self.print_redshift='z=0.1, $-18<M_r<-17$'                                    
                                elif key=='-21':
                                    self.print_redshift='z=0.1, $M_r<-21$'
                                
                                #self.print_redshift='$M_r<-21$, no dust, $n$: 0.116x$10^{-2}$'
                                #self.print_redshift='$-21<M_r<-20$, $n$: 0.53x$10^{-2}$'
                                #self.print_redshift='$-21<M_r<-20$, no dust, $n$: 0.53x$10^{-2}$'
                                #self.print_redshift='$-22<M_r<-21$'
                                #self.print_redshift='$\log_{10}$ ($M_{_*}$ $[M_{\odot}]) > 9$, $-22<M_r<-21$'  

                            else:
                                #SMF CMASS
                                prefix='SMF_Galacticus_1Gpc_z_0.56_new_mags'
                               # band={'0':'Galacticus_1Gpc_z_0.56_mstar'}#, '1': 'SAG_1Gpc_z_0.56_mstar', '2':'SAGE_1Gpc_z_0.56_1238_mstar'}#, '3':'SAG_1Gpc_v2_z_0.07_mcold_disk_sfr+h','4':'SAG_1Gpc_z_0.09_OII', '5':'SAG_1Gpc_z_0.09_mcold', '6':'SAG_1Gpc_z_0.09_mcold_disk', '7':'SAGE_1Gpc_z_0.09_OII', '8':'SAGE_1Gpc_z_0.09_mcold_disk' }
                                band={'0': prefix+'_CMASS', '2': prefix+'_CMASS_density_sample', '1': prefix+''}#, '3':'SAG_1Gpc_v2_z_0.07_mcold_disk_sfr+h','4':'SAG_1Gpc_z_0.09_OII', '5':'SAG_1Gpc_z_0.09_mcold', '6':'SAG_1Gpc_z_0.09_mcold_disk', '7':'SAGE_1Gpc_z_0.09_OII', '8':'SAGE_1Gpc_z_0.09_mcold_disk' }
                                #band={'0':'SMF_SAGE_1Gpc_z_0.56_CMASS_sample', '1': 'SAGE_1Gpc_z_0.56_mstar', '2':'SMF_SAGE_1Gpc_z_0.56_densCut_CMASS_mstar'}
                             
                                self.print_redshift='z=0.56'                    
                                self.plot_map_array['data_offset_'+plot_key] = 7                    
                                self.myconfig_array['nr_cats']=len(band)
   
                                self.myconfig_array.update({'catname'+str(i):self.myconfig_array['catname'+str(self.myPipe.a)]})
                                filename = mycomp+'anaconda/pro/myRun/histos/SMF/'+band[str(i)]+'.txt'   

                                mydata=myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
                                mydata[:,4]=np.log10(mydata[:,1])+(abs(mydata[:,1]-mydata[:,4])/mydata[:,1])*0.434 
                                mydata[:,5]=np.log10(mydata[:,1])-(abs(mydata[:,1]-mydata[:,5])/mydata[:,1])*0.434 
                          
                            print 'filenname:', filename                     
        
                            #print mydata[:,0].size, self.plot_map_array['data_offset_'+plot_key], int(self.plot_map_array['data_offset_'+plot_key])*self.myconfig_array['nr_cats']
                            if i==0:
                                my_analysed_data = np.zeros((mydata[:,0].size,int(self.plot_map_array['data_offset_'+plot_key])*self.myconfig_array['nr_cats']), dtype=np.float64)
                            
                            if plot=='wp':
                                #my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i] = np.log10(mydata[:,0])
                                my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i] = mydata[:,0]  
                            
                            #create a plot: r [Mpc] vs xi or wp
                            if f!='xi2':
                                #y-axis=xi(r)
                                my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] = mydata[:,1]
                                                               
                                if f=='wp' and len(mypanels)>1:
                                    add_obs=[0]
    
                            else:
                                #y-axis=r^2 * xi(r)
                                my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+0] = data[:,2]
                                my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] = data[:,2]**2.0 * mydata[:,1]
                                my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4] = data[:,0]**2.0 * mydata[:,1]
                                my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5] = data[:,1]**2.0 * mydata[:,1]
                                
                                if self.myconfig_array['catname'+str(i)].find('HAM')!=-1:
                                    my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+0] = data[:,0]
                                    my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] = data[:,1] 
                                #log error bars
                                #mydata[:,4]=np.log10(my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1])+(abs(my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1]-my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4])/my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1])*0.434 
                                #mydata[:,5]=np.log10(my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1])-(abs(my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1]-my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5])/my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1])*0.434

                                #normal error bars
                                #mydata[:,4]=my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] + (my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1]-my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4])
                                #mydata[:,5]=-my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] + my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1]-my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5]
                                #print my_analysed_data[0:5,[0,1,4,5]]                                                                

                                add_obs=[0]
    
        
                            #create a plot with a reference: r [Mpc] vs xi/xi_ref or wp/wp_ref
                            if f.find('ref')!=-1:
                                mean_ref=np.zeros((my_analysed_data[:,0].size), dtype=np.float32)
                                self.print_redshift=False
                                plot_legend=False 
                                add_obs=False

                                my_new_obs = np.zeros((my_analysed_data[:,0].size,int(self.plot_map_array['data_offset_'+plot_key])), dtype=np.float32)
                                if i==0:
                                    myref_data = np.zeros((my_analysed_data[:,0].size,self.myconfig_array['nr_cats']), dtype=np.float32)
                                        
                                if f.find('mean')!=-1:
                                    #calculate mean of different SAMs and us it as reference:
                                    d=0
                                    while d<my_analysed_data[:,0].size:
                                        e=0
                                        while e<self.myconfig_array['nr_cats']:
                                            #print 'd:', d, 'e:', e, my_analysed_data[d,int(self.plot_map_array['data_offset_'+plot_key])*e+1],
                                            mean_ref[d]+=my_analysed_data[d,int(self.plot_map_array['data_offset_'+plot_key])*e+1]
                                            e+=1
                                        #print 'sum:', mean_ref[d],
                                        mean_ref[d]/=self.myconfig_array['nr_cats']
                                        #print 'mean_ref:', mean_ref[d]
                                        d+=1
                                elif f.find('wp')!=-1:                
                                    print  'calculate residuals using observations as reference!'
                               
#                                    x=10**my_analysed_data[:,0]
#                                    #Calculate Gamma function:
#                                   
#                                    coeff = {'22-21_r_00': 5.98/0.6777, '22-21_r_01': (5.98+0.11)/0.6777, '22-21_r_02': (5.98-0.11)/0.6777, '22-21_gamma0': 1.92, '22-21_gamma1': 1.94, '22-21_gamma2': 1.90,
#                                             '21-20_r_00': 5.46/0.6777, '21-20_r_01': (5.46+0.15)/0.6777, '21-20_r_02': (5.46-0.15)/0.6777, '21-20_gamma0': 1.77, '21-20_gamma1': 1.79, '21-20_gamma2': 1.75,
#                                             '20-19_r_00': 4.82/0.6777, '20-19_r_01': (4.82+0.23)/0.6777, '20-19_r_02': (4.82-0.23)/0.6777, '20-19_gamma0': 1.87, '20-19_gamma1': 1.90, '20-19_gamma2': 1.84,
#                                             '19-18_r_00': 4.14/0.6777, '19-18_r_01': (4.14+0.30)/0.6777, '19-18_r_02': (4.14-0.30)/0.6777, '19-18_gamma0': 1.81, '19-18_gamma1': 1.84, '19-18_gamma2': 1.78,                                             
#                                             '18-17_r_00': 2.09/0.6777, '18-17_r_01': (2.09+0.38)/0.6777, '18-17_r_02': (2.09-0.38)/0.6777, '18-17_gamma0': 1.99, '18-17_gamma1': 2.13, '18-17_gamma2': 1.85,                                             
#                                             '-21_r_00': 5.98/0.6777, '-21_r_01': (5.98+0.12)/0.6777, '-21_r_02': (5.98-0.12)/0.6777, '-21_gamma0': 1.96, '-21_gamma1': 1.98, '-21_gamma2': 1.94}
#    
                                    my_new_obs[:,0]=my_analysed_data[:,0]
                                                                       
                                    a=0                                      
                                    myref_data[:,i]=np.interp(mydata[:,0] ,obs_data_array[0][str(0)][:,0], 10**obs_data_array[0][str(0)][:,1])
                                    my_new_obs[:,1]=np.log10(myref_data[:,i])                                        
                                        
                                    #print my_new_obs                                                                                                          
                                    obs_data_array[0].update({str(1+a): my_new_obs})                                     
                                    obs_data_array[1].update({str(1+a): ''})
                                    
                                else:
                                    print 'here else!'
                                    myref_data[:,i]=obs_data_array[0][str(0)][:,1]
                                    my_new_obs[:,1]=np.log10(myref_data[:,i])
                                    
                                    
                                    #print obs_data_array[0][str(0)][:,1]
                                    
                                  
                            else:
                                if plot=='wp':
                                    my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1]=np.log10(my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1])        
 
                            #error bars
                            my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4] =  np.log10(mydata[:,4]) 
                            my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5] =  np.log10(mydata[:,5])                   
                            i+=1
    
                        print '---------------------\n'

                        if f.find('ref')!=-1:
                            #print myref_data[:,0]
                            i=0
                            while i<self.myconfig_array['nr_cats']:
                                #print my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1]/myref_data[:,0]
                                my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1]= my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1]/myref_data[:,i] - 1.0
                                #print my_analysed_data[:,[int(self.plot_map_array['data_offset_'+plot_key])*i, int(self.plot_map_array['data_offset_'+plot_key])*i+1]]
                                i+=1

                        print 'len(myp):', len(mypanels)    
                        if len(mypanels)>1:                              
                            self.my_analysed_data.update({count: {'name': f, 'data': my_analysed_data, 'print_redshift': self.print_redshift, 'plot_legend': plot_legend, 'add_obs': add_obs}})
                        else:
                            self.my_analysed_data=my_analysed_data
                        count+=1
                    #print mypanels
                    #print self.my_analysed_data
                    #exit()
                    #self.my_analysed_data[0]['data'][:,1]=np.log10(myref_data[:,0])

                def skip():
                    pass

                def caseSwitcher(myplot):
                
                    choose = {
                        'plotLF': plotLF,
                        'plotAllBands': plotAllBands,
                        'plotOH2mstar': plotOH2mstar,
                        'plot2PCF': plot2PCF,
                        'plotManyCurves': plotManyCurves,
                        'skip': skip,
                        'SFH': SFH
                        }
                        
                    func = choose.get(myplot)
                    return func()
                
                caseSwitcher('SFH')

                        
###############################################################################################################################################################
###############################################################################################################################################################  

        def plotOutput(mycustom_plot_key):

            ts = time.time()              
            date_time_stamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H:%M:%S')

            myOutput.multiPlot(self.my_analysed_data,
                                 loops=self.myconfig_array['nr_cats'],
                                 my_add_subplot= obs_data_array[0],
                                 nr_added_subplots = self.subplot_loops,
                                 legend_added_subplots = obs_data_array[1],
                                 config_path=mycomp+'anaconda/pro/myRun/plot_config/', 
                                 myconfig_datafile=self.plot_map_array[plot_key+'_config'],
                                 mydir=mycomp+'anaconda/pro/myRun/plots/'+plot_key+'/',
                                 myfilename= self.myconfig_array[self.myconfig_array['catname0']+'_simulation_name']+'_'+self.myconfig_array['catname0']+'_'+plot_key+'_z_'+str(self.redshift)+'_'+str(date_time_stamp),
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
        
        while self.b<self.plot_map_array['nr_plot_keys']:
            
            obs_data_array = [[],[]]
            self.subplot_loops=0

            plot_key = self.plot_map_array['plot_map_id'+str(self.b)]
            
            if plot_key!='sfr2z':
                nbins=int(self.plot_map_array['nbins_'+plot_key])            
            else:
                nbins=self.myconfig_array['nr_zs']
                
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
                   
                    if self.load_from_file=='False' and plot_key.find('analyseTargetSelection')==-1: caseSwitcherPlotKey('mainCalculate', None, None, None)
                                       
                    elif plot_key.find('analyseTargetSelection')!=-1: caseSwitcherPlotKey('analyseTargetSelection')
                     
                    else:
                        if self.myPipe.a==0: LoadObs()
                        
                        if self.myPipe.a>0 and plot_key=='plotXY':                               
                            pass
                        else:
                            if self.plot_custom_loop=='True':
                                mymethod='method1'
                                if mymethod=='method1':
                                    mysample=['','red','blue','active','passive','low','high']
                                else:
                                    mysample=['','red','blue','passive','active','low','high','SFR3st','SFR3gt']

                                    
#                                for sample in mysample:
#                                    #plot growth of one sample
#                                    self.myPlot=caseSwitcherPlotKey('loadFromFile', sample, 'histo', '_mstar')
#                                    plotOutput('histo')

#                                for sample in mysample:
#                                    #plot growth of one sample
#                                    self.myPlot=caseSwitcherPlotKey('loadFromFile', sample, 'gr_one', '')
#                                    plotOutput('_gr_one')                                    
                                    

                                #plot these without errors
                                #for prop in ['_mhalo', '_mstar', '_SHMF', '_logSHMF', '_sfr', '_ssfr', '_mcold', '_Mzgas', '_zcold', '_g-i', '_r-i', 'SFHd', 'sumSFR']:
                                #plot these with errors
#                                for prop in ['_mhalo', '_mstar']:#, '_SHMF', '_sfr', '_ssfr', '_mcold', '_Mzgas', '_zcold', '_g-i', '_r-i']:
#                                #test plot
                                for prop in ['_mstar']:                                    
#                                    #plot property for all samples in one plot
#                                    #sample, key, prop
#                                    self.myPlot=caseSwitcherPlotKey('loadFromFile', '', 'SFH', prop)
#                                    plotOutput(prop)

                                #for prop in ['_mhalo', '_mstar', '_mcold', '_Mzgas']:
                                    #plot property for all samples in one plot
                                    self.myPlot=caseSwitcherPlotKey('loadFromFile', '', 'gr', prop)
                                    plotOutput('_gr')                                       
                                    
                                exit()
                                    
                                
                            else:
                                self.myPlot=caseSwitcherPlotKey('loadFromFile', None, None, None)                                
                                
                                

                        
                     
                    self.myPipe.a+=1
                                   
                self.myPipe.i+=1

                if self.load_from_file=='True' and self.plot_custom_loop=='False': plotOutput(None)
                
            self.b+=1
            
######## END MAIN       ##################################################################################################################################### 