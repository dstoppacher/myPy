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
                               '\n(1) z (2) cSFRD (sumSFR/volume) [Msun yr-1 Mpc-3] (3) sumSFR [Msun yr-1] (4) cSFRD/cSFRD(z=0) (5) comoving volume [Mpc3]',
                               data_format="%0.3f\t%0.8e\t%0.8e\t%0.8e\t%0.2e")

        return self.myPipe.histo_data_ZvsSFR  

        
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
        
        mymethod='FT'

        if mymethod=='M2':       
            #Galacticus
            SFH_analysis_keys=['','low', 'high', 'passive', 'active', 'red', 'blue', 'low-zcold', 'high-zcold', 'mhalo12gt', 'mhalo12st', 'zcold9gt', 'zcold9st', 'mstar11gt', 'mstar10st', 'redmstar11', 'redmhalo13', 'lowmstar11', 'lowmhalo13']
            SFH_analysis_keys=['','low-zcold','high-zcold']
        else:
            SFH_analysis_keys=['','low', 'high', 'passive', 'active', 'red', 'blue', 'low-zcold', 'high-zcold']
            
        SFH_analysis_keys=['r0001']          
        #LGALAXIES
        #SFH_analysis_keys=['','low', 'high', 'passive', 'active']        
 
        #SFH_analysis_keys=['lowZcold-highMstar']
        #SFH_analysis_keys=['']  
       
        for count, element in enumerate(SFH_analysis_keys):
            print '\n\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'
            element_space=''        
            #create spaces '_' in the filename:
            if element!='':
                element_space='_'      
            #print 'element:', element
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
                                            +'(1) z (2) sum(SFR [Msun yr-1]) (3) SFR/SFR(z=z_start) (4) 50th SFR [Msun yr-1] (5) '+perc_low+' SFR [Msun yr-1] (6) '+perc_high+' SFR [Msun yr-1] '\
                                            +'(7) sum(sSFR [yr-1]) (8) sSFR/sSFR(z=z_start) (9) 50th sSFR [yr-1] (10) '+perc_low+' sSFR [yr-1] (11) '+perc_high+' sSFR [yr-1] '\
                                            +'(12) 50th g-i (13) '+perc_low+' g-i (14) '+perc_high+' g-i '\
                                            +'(15) sum(Mstar [Msun]) (16) Mstar/Mstar(z=z_start) (17) 50th Mstar [Msun] (18) '+perc_low+' Mstar [Msun] (19) '+perc_high+' Mstar [Msun] '\
                                            +'(20) 50th rdisk [Mpc] (21) '+perc_low+' rdisk [Mpc] (22) 50th rbulge/rdisk [-] (23) '+perc_low+' rbulge/rdisk [-] (24) '+perc_high+' rbulge/rdisk [-] '\
                                            +'(25) lookback time(z=z_start) [Gyr] (26) frac centrals [-] (27) frac no-sats [-] (28) frac orphans [-] (29) n x 10-4 [Mpc-3] '\
                                            +'(30) 50th r-i (31) '+perc_low+' r-i (32) '+perc_high+' r-i '\
                                            +'(33) 50th mbh [-] (34) '+perc_low+' mbh [-] (35) '+perc_high+' mbh [-] '\
                                            +'(36) 50th rhalfmass [-] (37) '+perc_low+' rhalfmass  [-] (38) '+perc_high+' rhalfmass  [-] '\
                                            +'(39) sum(Mcold [Msun]) (40) Mcold/Mcold(z=z_start) (41) 50th Mcold [Msun] (42) '+perc_low+' Mcold [Msun] (43) '+perc_high+' Mcold [Msun] '\
                                            +'(44) sum(Mzgas [Msun]) (45) Mzgas/Mzgas(z=z_start (46) 50th Mzgas [Msun] (47) '+perc_low+' Mzgas [Msun] (48) '+perc_high+' Mzgas [Msun] '\
                                            +'(49) 50th zcold [-] (50) '+perc_low+' zcold [-] (51) '+perc_high+' zcold [-] '\
                                            +'(52) cosmic SFRD [Msun yr-1 Mpc-3] '\
                                            +'(53) sum(Mvir [Msun]) (54) Mvir/Mvir(z=z_start) (55) 50th Mvir [Msun] (56) '+perc_low+' Mvir [Msun] (57) '+perc_high+' Mvir [Msun] '\
                                            +'(58) 50th Mstar/Mvir [-] (59) '+perc_low+' Mstar/Mvir  [-] (60) '+perc_high+' Mstar/Mvir  [-] '\
                                            +'(61) 50th mean_age_stars_disk [Gyr] (62) '+perc_low+' mean_age_stars_disk [Gyr] (63) '+perc_high+' mean_age_stars_disk [Gyr] '\
                                            +'(64) 50th mean_age_stars_spheroid [Gyr] (65) '+perc_low+' mean_age_stars_spheroid [Gyr] (66) '+perc_high+' mean_age_stars_spheroid [Gyr] '\
                                            +'(67) 50th vmax [kms-1] (68) '+perc_low+' vmax [kms-1]  (69) '+perc_high+' vmax [kms-1] '\
                                            +'(70) 50th vdisp [kms-1] (71) '+perc_low+' vdisp [kms-1]  (72) '+perc_high+' vdisp [kms-1]',
                                   data_format="%0.4f\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.5f\t%0.5f\t%0.5f\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.5f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f")

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

        if self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_calc_fast_histo']=='True':
            
            mybinning='log'
            for prop in ['sfr']:
                print 'prop:', prop
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
            header='\n(1) '+my_name_x+' ['+self.mycond_configs[my_name_x+'_unit']+'] (2) mstar ['+self.mycond_configs[my_name_y+'_unit']+'] / '+my_name_y+' ['+self.mycond_configs[my_name_y+'_unit']+']\t(3) -dx\t(4) dx\t(5) -dy\t(6) +dy\t(7) N count (8) add axis '+myadd_axis+' ['+self.mycond_configs[myadd_axis+'_unit']+'] (9) dy add axis'               
        elif mydiv_y2x==True or myadd_axis!=False:           
            header='\n(1) '+my_name_x+' ['+self.mycond_configs[my_name_x+'_unit']+'] (2) '+my_name_y+' ['+self.mycond_configs[my_name_y+'_unit']+']\t(3) -dx\t(4) dx\t(5) -dy\t(6) +dy\t(7) N count (8) add axis '+myadd_axis+' ['+self.mycond_configs[myadd_axis+'_unit']+'] (9) -d add axis (10) +d add axis'
        elif mybinUp2D==True:
            header='\n(1) '+my_name_x+' ['+self.mycond_configs[my_name_x+'_unit']+'] (2) '+my_name_y+' ['+self.mycond_configs[my_name_y+'_unit']+']\t(3) -dx\t(4) dx\t(5) -dy\t(6) +dy\t(7) N count'
        else:
            header='\n(1) '+my_name_x+' ['+self.mycond_configs[my_name_x+'_unit']+'] (2) '+Phi+' ['+myvolume_units+dex+'\t(3) -dx\t(4) dx\t(5) -dy\t(6) +dy\t(7) N count'
      
        #total number of objects after all cuts used in calculation
        ngal_in_histo=sum(self.myPipe.ratio_data_binAndFrac2D[:, mydata_offset*self.myPipe.a+6])
        
        myOutput.writeIntoFile(filename,
                               self.myPipe.ratio_data_binAndFrac2D[:, mydata_offset*self.myPipe.a : mydata_offset*self.myPipe.a + mydata_offset],
                               myheader= plot_key+' '+self.myconfig_array['catname'+str(self.myPipe.a)]+' z='+str(float("{0:.2f}".format(self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i)]+'_snapid'+str(self.myPipe.i)]['z'])))+' cumulative: '+info_cum+starforming_cut+' ngal: '+str(int(ngal_in_histo))+header,
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
            self.myPipe.TwoPCF(False, filename, self.mycond_configs)

        def HODFunction(plot_key, data, filename):          
            self.myPipe.HODFunction(False, filename)
                       
        def FilterAndExtractData(plot_key, data, filename):
            self.myPipe.filterAndExtractData(self.mycond_configs)
            
        def AnalyseTargetSelection(plot_key, data, filename):
            return self.myPipe.analyseTargetSelection(data, filename, self.mycond_configs, plot_key)
            
        def MatchHaloCat(plot_key, data, filename):
            self.matchhalocat()
           
        def PlotXY(plot_key, data, filename):
            self.myPipe.plotXY(plot_key,filename)
          
        def caseSwitcherPlotKey(plot_key, sample, key, prop, mysamples, method, plot_num):

            def myLoadFromFile():
                LoadFromFile(sample, key, prop, mysamples, method, plot_num)
            
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
                        
            if self.b==0 or plot_key=='sfr2z' or plot_key=='cSFRD' or plot_key=='filterData' or plot_key=='plotXY' or plot_key=='twoPCF':
                if self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_data_format']!='HDF5':# and self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_data_format']!='CROSSMATCH':                     
                    self.myPipe.readData()
                else:
                    myfilename1, myfilename2, myfilename3, myfilename4, myfilename5, myfilename6 = check_filename()
                    #uncommand for test reason:
                    #--------------------
                    #myfilename=myfilename1
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

            #self.myPipe.showConfig()
           
            self.my_analysed_data = caseSwitcher(plot_key, data, filename)                
                

             

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
                                 print_redshift='z=0.56')#pop (ii), centrals')#, centrals')#\n$SFR<-1$')#CMASS DR12:')#$orphans$')# z='+str(float("{0:.2f}".format(myredshift))))#+' centrals')

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
                

           
        def LoadFromFile(mysample, mykey, myprop, mysamples, mymethod, myplot_num):

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
                    print self.my_analysed_data.shape
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
    
                if plot_key=='cSFRD':
                    print 'cSFRD plot:', self.myPipe.a, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1
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
                
                #print self.my_analysed_data[:,[int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a,int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a+1] ]
                #print self.my_analysed_data[0:3, int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a : int(self.plot_map_array['data_offset_'+plot_key])*self.myPipe.a + 6]
                #print self.my_analysed_data
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
                def plotZevol():
                    
                    prop=myprop
                    sample=mysample
                    
                    prefix1='IllustrisTNG300_HOD-SFgalaxies_z-evolution_'
                    prefix2='Galacticus_400Mpc_HOD-SFgalaxies_z-evolution_'
                    prefix3='Galacticus_1Gpc_HOD-SFgalaxies_z-evolution_'
                    

                    band={'0': prefix3+'CUT1_', '1': prefix3+'CUT2_', '2': prefix3+'CUT3_',\
                          '3': prefix2+'CUT1_', '4': prefix2+'CUT2_', '5': prefix2+'CUT3_',\
                          '6': prefix1+'CUT1_', '7': prefix1+'CUT2_', '8': prefix1+'CUT3_'} 
                    
                    self.myconfig_array['nr_cats']=len(band)
                    print len(band)
                                                         
                    i=0
                    while i<self.myconfig_array['nr_cats']:
                        if i==0:
                            self.my_analysed_data = np.zeros((6, 7*self.myconfig_array['nr_cats']), dtype=np.double)
                        print 'i:', i, 'sample:', sample, 'prop:', prop,
                        if prop.startswith('_mhalo')==True and band[str(i)].find('Ill')!=-1:
                            prop='_mhalo_200c'
                        elif prop.startswith('_mhalo')==True:
                            prop='_mhalo_cents_200c'
                            
                        filename=mycomp+'anaconda/pro/myRun/histos/HOD/HOD-SF/'+band[str(i)]+sample+prop+'.txt'
                        data=myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
                        print 'filenname:', band[str(i)]+sample+prop+'.txt'
                        
                        
                        mysample_name = {'mstar': '$M_*$',\
                                         'sfr': '$SFR$',\
                                         'ssfr': '$sSFR$'}
                                             
                        self.print_redshift=mysample_name[sample]+'-selected'
                        #self.print_redshift=''
                                             
                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i]    = data[:,0]                                                                     
                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1]  = data[:,1] 
                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4]  = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]+data[:,2]
                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5]  = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]-data[:,3]                                                                    
                        i+=1


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
                
                    
                    
                    self.plot_map_array['data_offset_'+plot_key] = 60

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
                    
                    print 'mysamples:', mysamples, '\n  method:\t', method, '\n  plot_num:\t', plot_num, '\n  sample:\t', sample, '\n  prop:\t\t', prop, '\n  key:\t\t', key
                    
                    #folder='Gal-dens_main_cents'
                    #folder='Gal-dens_massive'
                    #folder='Gal-dens_main_'+method
                    folder='Gal-dens_main_cents_'+method  
                    #folder='Galacticus_SDSS' 
                    #folder='Gal-dens_cents_cons_test_'+method
                    #folder='Gal2-dens_main_cents_'+method+'_500'
                    #folder='Gal2-dens_main_cents_'+method+'_Gal-halos'
                    #folder='Gal400-dens_main_cents_'+method
                    #folder='Gal_300_main_cents'
                    filename=mycomp+'anaconda/pro/myRun/histos/sfr2z/'+folder

                    if folder.find('300')!=-1:
                        #Gal-dens, main test suite
                        part2='/sfr2z_Galacticus_1Gpc_z_1.27_tarsel_SFH_300_main_cents'
                        part2_to_z='/sfr2z_Galacticus_1Gpc_z_'
                        part2_from_z_to_main='_tarsel_SFH_300_main_cents'                        
                        orphan=method
                        catname='Gal-r0001'
                        self.plot_map_array['data_offset_'+plot_key] = 72
                    elif folder.find('M1')!=-1 and folder.find('test')!=-1:
                        #Gal-dens, main test suite
                        part2='/sfr2z_Galacticus_1Gpc_z_4.15_tarsel_SFH_down3_method1_cents_test'
                        part2_to_z='/sfr2z_Galacticus_1Gpc_z_'
                        part2_from_z_to_main='_tarsel_SFH_down3_method1_cents_test'                        
                        orphan=method
                        catname='Gal-dens TC'
                    elif folder.find('400')!=-1 and method=='M1':
                        #Gal-dens, main suite
                        part2='/sfr2z_Galacticus_400Mpc_z_4.15_tarsel_SFH_down3_'+method+'_main_cents'
                        part2_to_z='/sfr2z_Galacticus_400Mpc_z_'
                        part2_from_z_to_main='_tarsel_SFH_down3_'+method+'_main_cents'                        
                        orphan=method
                        catname='Gal400-dens'
                        self.plot_map_array['data_offset_'+plot_key] = 72
                    elif folder.find('400')!=-1 and method=='M2':
                        #Gal-dens, main suite
                        part2='/sfr2z_Galacticus_400Mpc_z_4.15_tarsel_SFH_down3_M1_main_cents'
                        part2_to_z='/sfr2z_Galacticus_400Mpc_z_'
                        part2_from_z_to_main='_tarsel_SFH_down3_M1_main_cents'                        
                        orphan=method
                        catname='Gal400-dens'
                        self.plot_map_array['data_offset_'+plot_key] = 72                          
                    elif folder.find('M2')!=-1 and folder.find('test')!=-1:
                        #Gal-dens, main suite
                        part2='/sfr2z_Galacticus_1Gpc_m2_z_4.15_tarsel_SFH_down3_method2_cents_test'
                        part2_to_z='/sfr2z_Galacticus_1Gpc_m2_z_'
                        part2_from_z_to_main='_tarsel_SFH_down3_method2_cents_test'                        
                        orphan=method
                        catname='Gal-dens TC'  
                    elif folder.find('M1')!=-1 and folder.find('Gal2')!=-1:
                        #Gal-dens, main suite
                        part2='/sfr2z_Galacticus_1Gpc_run2_z_4.15_tarsel_SFH_down3_M1_main_cents_Gal-halos'
                        part2_to_z='/sfr2z_Galacticus_1Gpc_run2_z_'
                        part2_from_z_to_main='_tarsel_SFH_down3_M1_main_cents_Gal-halos'                        
                        orphan=method
                        self.plot_map_array['data_offset_'+plot_key] = 72 
                        catname='Gal2-dens H'                     
                    elif folder.find('M1')!=-1 or folder.find('M2')!=-1:
                        #Gal-dens, main
                        part2='/sfr2z_Galacticus_1Gpc_z_4.15_tarsel_SFH_down3_'+method+'_main_cents'
                        part2_to_z='/sfr2z_Galacticus_1Gpc_z_'
                        part2_from_z_to_main='_tarsel_SFH_down3_'+method+'_main_cents'                        
                        orphan=method
                        self.plot_map_array['data_offset_'+plot_key] = 72                        
                        catname='Gal-dens'                        

                    elif folder.find('SDSS')!=-1:
                        #Gal-dens, main
                        part2='/sfr2z_Galacticus_1Gpc_z_0.56_tarsel_CUT3_Contreras+13_mcold_SFH_method2_cents'
                        part2_to_z='/sfr2z_Galacticus_1Gpc_z_'
                        part2_from_z_to_main='_tarsel_CUT3_Contreras+13_mcold_SFH_method2_cents'
                        orphan='centrals'
                        catname='Gal-SDSS'
                    else:
                        method='M1'
                        folder='Gal-dens_main_cents_'+method
                        filename=mycomp+'anaconda/pro/myRun/histos/sfr2z/'+folder
                        part2='/sfr2z_Galacticus_1Gpc_z_4.15_tarsel_SFH_down3_'+method+'_main_cents'
                        part2_to_z='/sfr2z_Galacticus_1Gpc_z_'
                        part2_from_z_to_main='_tarsel_SFH_down3_'+method+'_main_cents'
                        orphan=method
                        self.plot_map_array['data_offset_'+plot_key] = 72                        
                        catname='Gal-dens'                         

                    prop_col_map={'': '',\
                                  'sumSFR': 1,\
                                  'SFH':    2,\
                                  '_sfr':   3,\
                                  '_ssfr':  8,\
                                  '_g-i':   11,\
                                  '_mstar': 16,\
                                  '_rdisk': 19,\
                                  '_rbulgevsrdisk': 21,\
                                  '_rbulge': 21,\
                                  '_r-i':   29,\
                                  '_mbh':   32,\
                                  '_rhalfmass': 35,\
                                  '_mcold': 40,\
                                  '_Mzgas': 45,\
                                  '_cgf':   45,\
                                  '_zcold': 48,\
                                  'SFHd':   51,\
                                  '_mhalo': 54,\
                                  '_SHMF':  57,\
                                  '_mean_age_stars_disk': 60,\
                                  '_mean_age_stars_spheroid': 63,\
                                  '_vmax':  66,\
                                  '_vdisp': 69,\
                                  '_Tcons': 40\
                                  }
 
                    prop_unit_map={'': '',\
                                   '_mean_age_stars_disk': '$age_{disk}$',\
                                   '_mean_age_stars_spheroid': '$age_{bulge}$',\
                                   '_vmax':     '$V_{max}$',\
                                   '_vdisp':    '$V_{disp}',\
                                   '_rbulgevsrdisk': '$r_{bulge}/r_{rdisk}$',\
                                   '_mbh':      '$M_{BH}$',\
                                   '_rdisk':    '$r_{disk}$',\
                                   '_rbulge':    '$r_{bulge}$',\
                                   '_rhalfmass':'$r_{1/2}$',\
                                   '_cgf':      '$M_{cold}/$M_{*}$',\
                                   '_mstar':    '$M_*$',\
                                   '_mhalo':    '$M_{vir}$',\
                                   '_mhalo_200c':'$M_{200c}$',\
                                   '_mcold':    '$M_{cold}$',\
                                   '_zcold':    '$Z_{cold}$',\
                                   '_Mzgas':    '$M_{z_{cold}}$',\
                                   '_sfr':      'SFR',\
                                   '_ssfr':     'sSFR',\
                                   '_g-i':      '$g-i$',\
                                   '_r-i':      '$r-i$',\
                                   'SFH':       '',\
                                   'SFHd':      '',\
                                   '_logSHMF': '$\log_{10}$ $M_{*}/$M_{vir}$',\
                                   'sumSFR':    '$\sum$ SFR [$M_{\odot}$ $yr^{-1}$]',\
                                   'sumMstar': '$\sum$ $M_{*}$ [$M_{\odot}$]',\
                                   'sumMhalo': '$\sum$ $M_{vir}$ [$M_{\odot}$]',\
                                   '_SHMF':     '$M_{*}$/$M_{vir}$',
                                   '_SHMR':     '$M_{*}$/$M_{vir}$',
                                   '_Tcons':    '$M_{cold}$/$SFR$ [$yr$]'}

                    if sample!='':
                        sample_code_space='_'
                    else:
                        sample_code_space=''
                    band={}
                    for l, items in enumerate(mysamples):
                        if l>-1:
                            item='_'+str(mysamples[l])
                        else:
                            item=str(mysamples[l])
                        band.update({str(l): item})

                    if key.find('one')!=-1 or key.find('histo')!=-1:
                        try:
                            mysample_name = {'': catname+'\n'+orphan+'\n$all$',\
                                             'low': catname+'\n'+orphan+'\n$low$-$SFR$',\
                                             'high': catname+'\n'+orphan+'\n$high$-$SFR$',\
                                             'passive': catname+'\n'+orphan+'\n$passive$',\
                                             'active': catname+'\n'+orphan+'\n$active$',\
                                             'red': catname+'\n'+orphan+'\n$red$',\
                                             'blue': catname+'\n'+orphan+'\n$blue$',\
                                             'redmstar11': catname+'\n'+orphan+'\n'+'$red$ + $M_*\sim11$',\
                                             'redmhalo13': catname+'\n'+orphan+'\n'+'$red$ + $M_{vir}\sim13$',\
                                             'lowmstar11': catname+'\n'+orphan+'\n'+'$low$-$SFR$ + $M_*\sim11$',\
                                             'lowmhalo13': catname+'\n'+orphan+'\n'+'$low$-$SFR$ + $M_{vir}\sim13$',\
                                             'mhalo12st': catname+'\n'+orphan+'\n'+'$M_{vir}<12$',\
                                             'mhalo12gt': catname+'\n'+orphan+'\n'+'$M_{vir}>12$',\
                                             'mstar10st': catname+'\n'+orphan+'\n'+'$M_*<11$',\
                                             'mstar11gt': catname+'\n'+orphan+'\n'+'$M_{*}>11$',\
                                             'zcold9st': catname+'\n'+orphan+'\n'+'$Z_{cold}<9$',\
                                             'zcold9gt': catname+'\n'+orphan+'\n'+'$Z_{cold}>9$',\
                                             'low-zcold': catname+'\n'+orphan+'\n'+'$low$-$Z_{cold}$',\
                                             'high-zcold': catname+'\n'+orphan+'\n'+'$high$-$Z_{cold}$',\
                                             'lowZcold-highMstar': catname+'\n'+orphan+'\n'+'$low$-$Z_{cold}$',\
                                             'highZcold-highMstar': catname+'\n'+orphan+'\n'+'$high$-$Z_{cold}$'}                                     
        
                            self.print_redshift=mysample_name[sample]
                        except:
                            print 'mysample_name is not set!'
                                                  
                                     
                    elif key=='gr':
                        self.print_redshift=catname+'\n'+orphan+'\n'+prop_unit_map[prop]                        
                    elif key.find('histo')!=-1:                   
                        if sample.find('zcold')==-1 and (sample.startswith('low') or sample.startswith('high')):
                            self.print_redshift=catname+'\n'+orphan+'\n'+sample+' SFR'
                        elif sample=='':
                            self.print_redshift=catname+'\n'+orphan+'\nall'                             
                        else:
                            self.print_redshift=catname+'\n'+orphan+'\n'+sample

                    else:
                        self.print_redshift=catname+'\n'+method+'\n'

                           
                    band={}

                    for l, items in enumerate(mysamples):
                        if l>-1:
                            item='_'+str(mysamples[l])
                        else:
                            item=str(mysamples[l])
                        band.update({str(l): item})

                    if plot_num=='_demo' and method=='M1':                            
                        self.print_redshift='\n'+'\n$M_{*}$\n$M_{vir}$\n$r_{1/2}$'#+'\nmethod 1 (samples selected at z\sim0.55)'
                    elif plot_num=='_demo' and method=='M2':
                        self.print_redshift=''#\n'+'\n$M_{*}$\n$M_{vir}$\n$M_{BH}$\n$r_{1/2}$'#\nmethod 2 (samples selected at each z)'                    
                    elif folder.find('SDSS')!=-1:
                        self.print_redshift=catname+'\n'+orphan+'\n$M_{cold}$-CUT3'
                        #self.print_redshift=mysample_name[sample] 
       
                    self.myconfig_array['nr_cats']=len(band)
 
                    if folder.find('SDSS')!=-1:
                        ncols=17
                    elif folder.find('300')!=-1:
                        ncols=38                        
                    elif folder.find('400')!=-1:
                        ncols=43
                    else:
                        ncols=55
                   
                    if key=='SFH' or key=='gr' or key=='gr_res':                    
                        i=0
                        while i<self.myconfig_array['nr_cats']:
                            if i==0:
                                self.my_analysed_data = np.zeros((ncols, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)
                            
                            col=prop_col_map[prop]

                            if plot_num=='_demo':
                                print 'choose galaxy props for demo sample!'                            
                                if i==1 or i==4:
                                    col=prop_col_map['_mhalo']
                                if i==2 or i==5:
                                    col=prop_col_map['_rhalfmass']
                                    
                            print 'i:', i, 'col:', col,  'filenname:', filename+part2+sample_code_space+sample+band[str(i)]+'.txt'                               
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i : int(self.plot_map_array['data_offset_'+plot_key])*i + int(self.plot_map_array['data_offset_'+plot_key])] = myData.readAnyFormat(config=False, mypath=filename+part2+sample_code_space+sample+band[str(i)]+'.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
                                
                            if prop=='sumSFR':
                                if folder.find('400')!=-1:
                                    norm_y=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+28]*(590.23**3)*1e-4
                                else:
                                    norm_y=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+28]*(1450.0**3)*1e-4
                                    norm_y=1450.0**3
                            elif prop=='_cgf':
                                norm_y=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+16]
                            elif prop=='_Tcons':
                                norm_y=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+3]
                            elif prop=='_rbulge':
                                norm_y=(1.0/self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+19])/1000.0
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+22]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+22]*self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+20]*1000.0 
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+23]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+23]*self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+20]*1000.0 
                            elif prop=='_rdisk':
                                norm_y=1.0/1000.0
                            elif prop=='_rhalfmass':
                                norm_y=1.0/1000.0                                    
                            else:
                                norm_y=1
                                                     
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+0]+=0
                            #self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+0]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+24]
                            if key.startswith('gr'):
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+col]/self.my_analysed_data[0, int(self.plot_map_array['data_offset_'+plot_key])*i+col]
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
            
                                    print 'filename_res:', filename_res+part2_res+sample_code_space+sample+'.txt'#, 'sample_code_space:', sample_code_space 
                                    
                                    self.my_analysed_data_res[:, int(self.plot_map_array['data_offset_'+plot_key])*i : int(self.plot_map_array['data_offset_'+plot_key])*i + int(self.plot_map_array['data_offset_'+plot_key])] = myData.readAnyFormat(config=False, mypath=filename_res+part2_res+sample_code_space+sample+band[str(i)]+'.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
            
                                    self.my_analysed_data_res[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=self.my_analysed_data_res[:, int(self.plot_map_array['data_offset_'+plot_key])*i+col]/self.my_analysed_data_res[0, int(self.plot_map_array['data_offset_'+plot_key])*i+col]
                                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]/self.my_analysed_data_res[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]-1                                
                         
                            else:
                                if prop=='_logSHMF' or prop=='_logSHMF_200c':
                                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=np.log10(self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+col]/norm_y)
                                    self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]+self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+col+1]/10**self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]*0.434
                                    self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+2]-self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+col+1]/10**self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]*0.434                                                                      
                                else:
                                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+col]/norm_y
                                    
                            if prop!='SFHd' and prop!='_logSHMF' and prop!='_logSHMF_200c':
                                self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+col+1]
                                self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+col+2]                                                                      

#                            share_axis=prop_col_map['_mcold']
#                            self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+6]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+share_axis]/self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+prop_col_map['_mstar']]
#                            share_axis=prop_col_map['_sfr']
#                            self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+6]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+share_axis]
#                            self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+7]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+share_axis+1]
#                            self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+8]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+share_axis+2]

                            i+=1
                            #print self.my_analysed_data[:, [0,1]]
#                            print norm_y

                    elif key=='stats_violin':
                        self.print_redshift='z='+sample
                        self.plot_map_array['data_offset_'+plot_key] = 3
                        i=0
                        for item in mysamples:
                            print 'z:', sample, 'sample:', item, 'prop:', prop,
                            import pandas as pd
                            data1 = pd.read_csv(mycomp+'anaconda/pro/data/Galacticus_1Gpc/SFH/SFH_Index_'+str(sample)+'_'+item+'_props_M1.txt', skiprows=2, names=['haloid', 'parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar',  'SHMR', 'mcold',  'Mzgas',  'zcold', 'g-i', 'sfr', 'ssfr'], sep='\t')
                            data_M1 = mL.df_to_sarray(data1)
                            data2 = pd.read_csv(mycomp+'anaconda/pro/data/Galacticus_1Gpc/SFH/SFH_Index_'+str(sample)+'_'+item+'_props_M2.txt', skiprows=2, names=['haloid', 'parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar',  'SHMR', 'mcold',  'Mzgas',  'zcold', 'g-i', 'sfr', 'ssfr'], sep='\t')
                            data_M2 = mL.df_to_sarray(data2)
                            print 'size M1:', data_M1.size, 'size M2:', data_M2.size
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
                        print 'HERE! 1341 stats ...'
                        self.plot_map_array['data_offset_'+plot_key] = 23                        

                        i=0
                        for sample in mysamples:
                            self.print_redshift=prop_unit_map[prop]
                            if i==0:
                                self.my_analysed_data = np.zeros((ncols, int(self.plot_map_array['data_offset_'+plot_key])*len(mysamples)), dtype=np.float64)
                                #print np.info(self.my_analysed_data)
                            print 'i:', i, 'sample:', sample, 'prop:', prop, 'filenname:', filename+'/stats/Galacticus_1Gpc_SFH_z-evolution_stats'+band[str(i)]+prop+'.txt'                               
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i : int(self.plot_map_array['data_offset_'+plot_key])*i + int(self.plot_map_array['data_offset_'+plot_key])] = myData.readAnyFormat(config=False, mypath=filename+'/stats/Galacticus_1Gpc_SFH_z-evolution_stats'+band[str(i)]+prop+'.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=1)                                
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
                            print 'col:', col
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+col]                                
                            #print self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]
                            i+=1
                                                    
                    elif key.find('histo')!=-1:

                        if folder.find('400')!=-1:
                            if key=='histo1':
                                band={'0': 0.55, '1': 0.74, '2': 0.84, '3': 0.89, '4': 1.0, '5': 1.26}
                            elif key=='histo2':
                                band={'0': 1.36, '1': 1.56, '2': 2.02, '3': 2.56, '4': 3.03, '5': 4.15}
                            else:
                                band={'0': 0.55, '1': 0.74, '2': 0.89, '3': 2.02, '4': 3.03, '5': 4.15} 
                            
                        else:    
                            if key=='histo1':
                                band={'0': 0.56, '1': 0.74, '2': 0.82, '3': 0.9, '4': 0.99, '5': 1.27}
                            elif key=='histo2':
                                band={'0': 1.37, '1': 1.54, '2': 2.03, '3': 2.53, '4': 3.04, '5': 4.15}
                            else:
                                band={'0': 0.56, '1': 0.74, '2': 0.99, '3': 2.03, '4': 3.04, '5': 4.15}                            

                        self.plot_map_array['data_offset_'+plot_key] = 7
                        
                        self.myconfig_array['nr_cats']=len(band)
                        i=0
                        while i<self.myconfig_array['nr_cats']:
                            if i==0:
                                self.my_analysed_data = np.zeros((25, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)
                                
                            print 'i:', i, 'prop:', prop, 'filenname:', filename+'/histos'+part2_to_z+str(band[str(i)])+part2_from_z_to_main+sample_code_space+sample+prop+'_histo.txt'       

                            data=myData.readAnyFormat(config=False, mypath=filename+'/histos'+part2_to_z+str(band[str(i)])+part2_from_z_to_main+sample_code_space+sample+prop+'_histo.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] = data[:,6]
                            if prop.find('zcold')!=-1 or prop.find('-')!=-1 or prop.find('log')!=-1 or prop.find('age')!=-1:
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i] =data[:,0]
                            elif prop.find('SHMF')!=-1:
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i] =np.log10(data[:,1])
                            else:
                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i] =np.log10(data[:,0])

                            
                            #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] + (data[:,4]/data[:,1])*0.43
                            #self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] - (data[:,5]/data[:,1])*0.43
                            #print self.my_analysed_data[:, [ int(self.plot_map_array['data_offset_'+plot_key])*i, int(self.plot_map_array['data_offset_'+plot_key])*i+1]]        

                            i+=1                           
                    elif key=='gr_one': 
                        #growth history for one sample but many parameters
                        #mahlo, mhalo, mbh, rhalfmass, mbh, rhalfmass, rbulgevsrdisk
                        band={'0': 16, '1': 54, '2': 40, '3': 45, '4': 32, '5': 35, '6': 21}
                        #2x mstar, mhalo, mbh, rhalfmass, mbh,
                        band={'0': 16, '1': 54, '2': 32, '3': 35, '4': 16, '5': 54, '6': 32, '7': 35}
                        #mhalo, mstar, mcold, Mzgas, mbh, rhalfmass, zcold
                        band={'0': 54, '1': 16, '2': 40, '3': 45, '4': 32, '5': 35, '6': 48}   
    
                        self.print_redshift='\nlow-$Z_{cold}$\n$true$'   
    
                        self.myconfig_array['nr_cats']=len(band)
                        i=0
                        while i<self.myconfig_array['nr_cats']:
                            if i==0:
                                self.my_analysed_data = np.zeros((ncols, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)
    
                            print 'i:', i, 'filenname:', filename+part2+sample_code_space+sample+'.txt'
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i : int(self.plot_map_array['data_offset_'+plot_key])*i + int(self.plot_map_array['data_offset_'+plot_key])] = myData.readAnyFormat(config=False, mypath=filename+part2+sample_code_space+sample+'.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
        
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+0]+=0
                            #self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+0]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+24]
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+band[str(i)]]/self.my_analysed_data[0, int(self.plot_map_array['data_offset_'+plot_key])*i+band[str(i)]]                                                                    
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+6] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+24]
                            i+=1                        

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
                                print 'more wp props ...'
                                endfix='_0.5_150_150'
                                myfilename=filename+'/wp/twoPCF_Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3_more_props_test_'+sample+'_centrals_'+prop[i]+endfix+'_wp.txt'
                                #twoPCF_Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3_more_props_test_low-zcold_centrals_PopB+k_0.5_150_150_wp.txt
                            else:                                
                                print 'z=', prop[i]
                                myfilename=filename+'/wp/'+part2_to_z+prop[i]+part2_from_z_to_main+sample_code_space+sample+endfix+'_wp.txt'
                            
                            print myfilename
                            data = myData.readAnyFormat(config=False, mypath=myfilename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
                          
                            self.my_analysed_data[:, [int(self.plot_map_array['data_offset_'+plot_key])*i, int(self.plot_map_array['data_offset_'+plot_key])*i+1, int(self.plot_map_array['data_offset_'+plot_key])*i+4, int(self.plot_map_array['data_offset_'+plot_key])*i+5]]=data[:,[2,3,5,5]]

                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]*=data[:,2]
                            
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=np.log10(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1])
                            i+=1

                    elif key.find('wp')!=-1:
                         
                        print 'here clustering!'
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

                            print filename+'/wp/'+part2_to_z+prop+part2_from_z_to_main+sample_code_space+band[str(i)]+endfix+'_wp.txt'
                            data = myData.readAnyFormat(config=False, mypath=filename+'/wp/'+part2_to_z+prop+part2_from_z_to_main+sample_code_space+band[str(i)]+endfix+'_wp.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
                          
                            self.my_analysed_data[:, [int(self.plot_map_array['data_offset_'+plot_key])*i, int(self.plot_map_array['data_offset_'+plot_key])*i+1, int(self.plot_map_array['data_offset_'+plot_key])*i+4, int(self.plot_map_array['data_offset_'+plot_key])*i+5]]=data[:,[2,3,5,5]]

                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]*=data[:,2]
                            
                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]=np.log10(self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1])
                            
                            
                            #print self.my_analysed_data[:, [int(self.plot_map_array['data_offset_'+plot_key])*i, int(self.plot_map_array['data_offset_'+plot_key])*i+1]]
                            i+=1
                    elif key=='envr_props':
                        print 'key:', key
                        self.my_analysed_data = np.zeros((15, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)
                        self.print_redshift='low-$Z_{cold}$'





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
                    band={'0':'0.0', '1': '0.49', '2':'0.99', '3':'2.03','4':'3.04', '5':'4.04', '6':'5.02' }
                    #band={'0':'0.0', '1': '0.49', '2':'1.0', '3':'2.02','4':'3.03', '5':'3.95', '6':'5.7' }

                    self.myconfig_array['nr_cats']=7
                    self.plot_map_array['data_offset_'+plot_key]=10
                    i=0
                    while i<self.myconfig_array['nr_cats']:
                        if i==0:
                            self.my_analysed_data = np.zeros((25, self.plot_map_array['data_offset_'+plot_key]*self.myconfig_array['nr_cats']), dtype=np.double)
                        self.myconfig_array.update({'catname'+str(i):'Galacticus_1Gpc'})
                        #filename = mycomp+'anaconda/pro/myRun/histos/oh2mstar/oh2mstar_'+band[str(i)]+'.txt'
                        filename = mycomp+'anaconda/pro/myRun/histos/oh2mstar/oh2mstar_SAGE_1Gpc_z_'+band[str(i)]+'_tarsel_ssfr_Hen+20.txt' 
                        #filename = mycomp+'anaconda/pro/myRun/histos/zgas2mstar/zgas2mstar_Galacticus_1Gpc_z_'+band[str(i)]+'_tarsel.txt'
                        #filename = mycomp+'anaconda/pro/myRun/histos/zgas2mstar/zgas2mstar_Galacticus_400Mpc_z_'+band[str(i)]+'_tarsel_ssfr_Hen+20.txt'                        
                        print 'i:', i, 'filenname:', filename, band[str(i)]
                        
                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i : int(self.plot_map_array['data_offset_'+plot_key])*i + int(self.plot_map_array['data_offset_'+plot_key])] = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i]=np.log10(self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i])
                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] - self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+4]
                        self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5] = self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+1] + self.my_analysed_data[:,int(self.plot_map_array['data_offset_'+plot_key])*i+5]

#                        if band[str(i)].find('OII')!=-1:
#                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1]+=0.25
                        i+=1

                        
                    #++++++++++++++++++++++++++++++++++++++++++++
                    #SMF CMASS and HOD
#                    key='HOD'
#                    key2='scd'
# 
#                    #sample='Gal400-dens'# f$^{rand}_{sats}$'                        
#                    #sample='Gal2-dens haloid from Gal-dens'# f$^{rand}_{sats}$'
#                    #sample='dens $V_{max}$'
#                    #prefix=key+'_LGALAXIES_500Mpc_z_0.56_tarsel_CMASS_density_mstar_no_mhalo'
#                    #prefix=key+'_SAG_1Gpc_z_0.56_new_CMASS_down_sample_mhalo_cents'
#                    
#                    if key=='SMF'or key=='HMF' or key.find('SFR')!=-1:
#                        prefix='SMF_Galacticus_1Gpc_z_0.56_new_mags'
#                        #prefix='SMF_SAG_1Gpc_z_0.56_new'
#                        
#                        #prefix='SMF_SAGE_1Gpc_z_0.56_tarsel_v3_mags_run_1238'
#                        prefix=key+'_LGALAXIES_500Mpc_z_0.56_tarsel'                        
#                        prefix=key+'_Galacticus_1Gpc_z_'
#                        prefix='SMF_SAG_1Gpc_v2_z_0.07'
#                        #prefix=key+'_Galacticus_1Gpc_z_0.56_tarsel_new_mags'
#                        
#                        
#                        #band={'0':'Galacticus_1Gpc_z_0.56_mstar'}#, '1': 'SAG_1Gpc_z_0.56_mstar', '2':'SAGE_1Gpc_z_0.56_1238_mstar'}#, '3':'SAG_1Gpc_v2_z_0.07_mcold_disk_sfr+h','4':'SAG_1Gpc_z_0.09_OII', '5':'SAG_1Gpc_z_0.09_mcold', '6':'SAG_1Gpc_z_0.09_mcold_disk', '7':'SAGE_1Gpc_z_0.09_OII', '8':'SAGE_1Gpc_z_0.09_mcold_disk' }
#                        #band={'0': prefix+'_CMASS', '2': prefix+'_CMASS_density_sample', '3': prefix+'_CMASS_mass_sample', '1': prefix}#, '3':'SAG_1Gpc_v2_z_0.07_mcold_disk_sfr+h','4':'SAG_1Gpc_z_0.09_OII', '5':'SAG_1Gpc_z_0.09_mcold', '6':'SAG_1Gpc_z_0.09_mcold_disk', '7':'SAGE_1Gpc_z_0.09_OII', '8':'SAGE_1Gpc_z_0.09_mcold_disk' }
#                        #band={'0':'SMF_SAGE_1Gpc_z_0.56_CMASS_sample', '1': 'SAGE_1Gpc_z_0.56_mstar', '2':'SMF_SAGE_1Gpc_z_0.56_densCut_CMASS_mstar'}
#                        #band={'0': prefix+'_CMASS', '2': prefix+'_test_no5h_CMASS', '3': prefix+'_test_nozboost_no5h_CMASS', '1': prefix+'_test_nozboost_no5h_noiltdmesa_CMASS'}#, '3':'SAG_1Gpc_v2_z_0.07_mcold_disk_sfr+h','4':'SAG_1Gpc_z_0.09_OII', '5':'SAG_1Gpc_z_0.09_mcold', '6':'SAG_1Gpc_z_0.09_mcold_disk', '7':'SAGE_1Gpc_z_0.09_OII', '8':'SAGE_1Gpc_z_0.09_mcold_disk' }
#
#                        #band={'0': prefix+'_centrals', '1': prefix+'_sats', '2': prefix+'_CMASS_centrals', '3': prefix+'_CMASS_sats', '4': prefix+'_CMASS_density_sample_no', '5': prefix+'_CMASS_density_sample_sats'}#, '4': prefix+'_mass_sample_centrals', '5': prefix+'_mass_sample_sats'}
#                        band={'0': prefix+'_CMASS', '2': prefix+'_CMASS_down_sample3', '3': prefix+'_CMASS_mass_sample', '1': prefix}
#                        band={'0': prefix+'_CMASS_down_sample3', '1': prefix+'_CMASS_density_sample_mhalo', '2': prefix+'_CMASS_density_sample_mstar', '3': prefix}                        
#                        #band={'0': prefix+'', '1': prefix+'_g-i_gt_2.35', '2': prefix+'_g-i_gt_2.15', '3': prefix+'_g-i_gt_2.0'}
#                        
#                        #band={'0': prefix+'_CMASS_down_sample', '1': prefix+'_35bins_2e10', '2': prefix+'_CMASS_down_sample_no', '3': prefix+'_35bins_2e10_no'}
#                        
#                        band={'0': prefix+'_CMASS_down_sample_test', \
#                              '1': prefix+'_35bins_2e10', \
#                              '2': prefix+'_CMASS_density_mstar', \
#                              '3': prefix+'_CMASS_density_vmax'}#, \
#                              #'4': prefix+'_density_vmax', \
#                              #'5': prefix+'_density_vmax_sat'}
#                              
#                        band={'0': prefix+'_OII', \
#                              '1': prefix+'_test_OII'}#, \
#                              #'2': prefix+'_test_OII'}#, \
##                              '4': prefix+'_density_vmax', \
##                              '5': prefix+'_density_vmax_sat'
#
##                        band={'0': prefix+'_CMASS_down_sample_test', \
##                              '1': prefix+'_CMASS_density_vmax', \
##                              '2': prefix+'_CMASS_density_mstar'}#, \
##                              '4': prefix+'_density_vmax', \
##                              '5': prefix+'_density_vmax_sat'
#                        redshift_MD = '0.0'
#                        redshift_SMD = '0.0'
#                        redshift_TNG = '0.0'
#                        rand_sample='0.008'
#                        band={'0': 'MDPL2_Rockstar_z_'+redshift_MD, '1': 'SMDPL_Rockstar_z_'+redshift_SMD, '2': 'IllustrisTNG300_z_'+redshift_TNG,\
#                              '3': 'Galacticus_1Gpc_z_'+redshift_MD+'_tarsel', '4': 'Galacticus_400Mpc_z_'+redshift_SMD+'_tarsel', '5': 'IllustrisTNG300_z_'+redshift_TNG+'_tarsel_n_0.008'}
#                        
#                    elif key2=='sc':
#                        prefix='HOD_Galacticus_1Gpc_z_0.56_tarsel_new_mags'
#                        #prefix='HOD_SAG_1Gpc_z_0.56_tarsel_new'
#                        
#                        prefix1=key+'_IllustrisTNG300_z_0.0_tarsel_mhalo_200c' 
#                        prefix2=key+'_Galacticus_400Mpc_z_0.0_tarsel_mhalo' 
#                        prefix3=key+'_Galacticus_1Gpc_z_0.0_tarsel_mhalo'
#                    
#                        band={'0': prefix3, '1': prefix2, '2': prefix1}
#                                             
#                        #band={'0': prefix+'_CMASS', '1': prefix+'_CMASS_down_sample3', '2': prefix+'_CMASS_mass_sample'}#, '2':prefix+'_CMASS_down_sample', '1':prefix+'_CMASS_density_sample'}
#                        #band={'0': prefix+'_CMASS_down_sample', '1': prefix+'_CMASS_mass_sample'}
#
#                    else:
#                        # HOD for each sample seperatly
#                        #prefix=key+'_Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3_HOD_mhalo_cents_200c'
#                        #prefix=key+'_Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3_HOD_FOF_halos_mhalo_cents_200c'
#                        #prefix=key+'_Galacticus_1Gpc_run2_z_0.56_tarsel_CMASS_down_sample3_haloid_from_Gal-dens_Rockstar_halos_mhalo_cents_200c'
#                        prefix=key+'_Galacticus_1Gpc_run2_z_0.56_tarsel_CMASS_density_sample_mstar_Rockstar_halos_mhalo_cents_200c'                        
#                        #prefix=key+'_Galacticus_400Mpc_z_0.55_tarsel_CMASS_down_sample3_Rockstar_halos_mhalo_cents_200c' 
#                        #prefix=key+'_Galacticus_1Gpc_run2_z_0.56_tarsel_CMASS_down_sample3_haloid_from_Gal-dens_mhalo_cents_200c'
#                        print 'HERE test 1522', self.redshift
#
#                        #for HOD-SFgalaxies
#                        cut='CUT3'
#                        prop='_ssfr'
#                        prefix1='HOD-SF/'+key+'_IllustrisTNG300_z_0.0_tarsel_HOD-SF_'+cut+prop+'_mhalo_200c' 
#                        prefix2='HOD-SF/'+key+'_Galacticus_400Mpc_z_0.0_tarsel_HOD-SF_'+cut+prop+'_mhalo_cents_200c' 
#                        prefix3='HOD-SF/'+key+'_Galacticus_1Gpc_z_0.0_tarsel_HOD-SF_'+cut+prop+'_mhalo_cents_200c'
#
#                        band={'0': prefix3+'_all', '1': prefix3+'_centrals', '2': prefix3+'_sats',\
#                              '3': prefix2+'_all', '4': prefix2+'_centrals', '5': prefix2+'_sats',\
#                              '6': prefix1+'_all', '7': prefix1+'_centrals', '8': prefix1+'_sats'} 
# 
#
#
#                       
#                        #for CMASS SFH paper
#                        prefix2=key+'_Galacticus_400Mpc_z_0.55_tarsel_CMASS_down_sample3_Rockstar_halos_mhalo_cents_200c' 
#                        prefix1=key+'_Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3_HOD_Rockstar_halos_mhalo_cents_200c'
##
#                        band={'0': prefix1+'_all', '1': prefix1+'_centrals', '2': prefix1+'_sats',\
#                              '3': prefix2+'_all', '4': prefix2+'_centrals', '5': prefix2+'_sats'}   
#                       

                    #endfix='_mstar'
                    #green butterfly test
                    #band={'0': 'SMF_Galacticus_1Gpc_z_0.56_new_mags_CMASS', '1': prefix+'_sample2'+endfix, '2': prefix+endfix, '3': prefix+'_sample2_centrals'+endfix, '4': prefix+'_centrals'+endfix, '5': prefix+'_sample2_orphans'+endfix, '6': prefix+'_orphans'+endfix, '7': prefix+'_sample2_no-sats'+endfix, '8': prefix+'_no-sats'+endfix}
                 
#                    if self.myconfig_array['catname'+str(self.myPipe.a)].startswith('SAG_'):                       
#                        band={'1':'SMF_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_0.56', '0':'SMF_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_0.56_densCut_CMASS_mstar'}
#                    else:
#                        band={'0':'SMF_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_0.56_CMASS', '1': 'SMF_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_0.56', '2':'SMF_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_0.56_densCut_CMASS_mstar'}#, '3':'SMF_'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_z_0.56_densCut_CMASS_mstar'}

                    #self.print_redshift='Gal2-dens mstar'
#                    self.print_redshift='z$\sim0.55$'                   
#                    #self.print_redshift='z$\sim$0.1'
#                    self.plot_map_array['data_offset_'+plot_key] = 7                    
#
#                    i=0
#                    while i<len(band):
#
#                        self.myconfig_array.update({'catname'+str(i):self.myconfig_array['catname'+str(self.myPipe.a)]})
#                        if key2=='sc':
#                            filename = mycomp+'anaconda/pro/myRun/histos/'+key+'/'+band[str(i)]+'_sc_centrals.txt'
#                            
#                        elif key=='HMF':
#                            if i<3:
#                                filename = mycomp+'anaconda/pro/myRun/histos/'+key+'/'+key+'_'+band[str(i)]+'_Nhalos.txt'
#                            else:
#                                filename = mycomp+'anaconda/pro/myRun/histos/'+key+'/'+key+'_'+band[str(i)]+'.txt'
#                        else:
#                            filename = mycomp+'anaconda/pro/myRun/histos/'+key+'/'+band[str(i)]+'.txt'  
#
#                        print 'i:', i, 'filenname:', filename#, band[str(i)]
#                        data=myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
#                        if i==0:
#                            self.my_analysed_data = np.zeros((data[:,0].size, self.plot_map_array['data_offset_'+plot_key]*len(band)), dtype=np.float32)
#                            
#                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i] = np.log10(data[:,0])#+0.03925
#                        self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] = np.log10(data[:,1])
#                        
#                        if key=='HOD' and key2!='sc':
#                            #self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] = np.log10((data[:,5])/data[:,6])
#                            
#                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+4] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] + (1.0/np.sqrt(data[:,6]))/data[:,1]*0.434
#                            self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+5] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] - (1.0/np.sqrt(data[:,6]))/data[:,1]*0.434                                     
#                        else:
#                            if key2=='sc':
#                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] = data[:,1]
#                                
#                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+4] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] + data[:,4]
#                                self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+5] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] - data[:,5]              
#                            else:
#                                try:
#                                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+4] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] + data[:,4]/data[:,1]*0.434
#                                    self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+5] = self.my_analysed_data[:, int(self.plot_map_array['data_offset_'+plot_key])*i+1] - data[:,5]/data[:,1]*0.434                  
#                                except:
#                                    pass
#                                
#                        #print self.my_analysed_data[:, [int(self.plot_map_array['data_offset_'+plot_key])*i+0, int(self.plot_map_array['data_offset_'+plot_key])*i+1]]
#                        i+=1  

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

                def plot2PCF(SFH=False,
                             my_custom_filename=False):

                    #plot 2PCF:
                    #------------------------------------------------------------------------------------------------------------------------------
                  
                    LoadObs()
                    cut=''
                    gtype_general='_centrals'
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
                    #prefix=gtype_general+'_mstar_gt_10.0'
                    pimax='_150'
                    endfix=prefix+'_0.5_150'
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
#                        band={0: 'Galacticus_1Gpc_z_0.56_CMASS_down_sample3'+endfix+pimax, 1: 'Galacticus_1Gpc_run2_z_0.56_CMASS_down_sample3_haloid_from_Gal-dens'+endfix+pimax, 2: 'Galacticus_400Mpc_z_0.55_CMASS_down_sample3'+endfix+pimax}
                        #band={0: 'Galacticus_1Gpc_run2_z_0.56_CMASS_down_sample3'+endfix+pimax, 1: 'Galacticus_1Gpc_run2_z_0.56_CMASS_density_sample_mhalo'+endfix+pimax, 2: 'Galacticus_1Gpc_run2_z_0.56_CMASS_density_sample_mstar'+endfix+pimax}
                        band={0: 'Galacticus_1Gpc_z_0.56_CMASS_down_sample3'+endfix+pimax, 1: 'Galacticus_400Mpc_z_0.55_CMASS_down_sample3'+endfix+pimax}

                        #band={0: '_color_sample'+endfix, 1:'_down_sample3'+endfix, 2: '_mass_sample'+endfix}
                        #band={0: '_down_sample'+endfix, 1:'_density_sfr_lowest_no'+endfix, 2: '_density_ssfr_lowest_no'+endfix}
                        #galaxy_type='no'
#                        band={0: '_new_posCor_-19_Mr_-18_wp', 1: '_new_posCor_-20_Mr_-19_wp', 2: '_new_posCor_-21_Mr_-20_wp', 3: '_new_posCor_-22_Mr_-21_wp', \
#                              4:'_CUT1_Contreras+13_mstar_'+galaxy_type+'_wp', 5: '_CUT2_Contreras+13_mstar_'+galaxy_type+'_wp', 6: '_CUT3_Contreras+13_mstar_'+galaxy_type+'_wp',\
#                              7:'_CUT1_Contreras+13_mcold_'+galaxy_type+'_wp', 8: '_CUT2_Contreras+13_mcold_'+galaxy_type+'_wp', 9: '_CUT3_Contreras+13_mcold_'+galaxy_type+'_wp',\
#                              10:'_CUT1_Contreras+13_sfr_'+galaxy_type+'_wp', 11: '_CUT2_Contreras+13_sfr_'+galaxy_type+'_wp', 12: '_CUT3_Contreras+13_sfr_'+galaxy_type+'_wp'} 

#                        band={0: '_new_posCor_-19_Mr_-18_wp', 1: '_new_posCor_-20_Mr_-19_wp', 2: '_new_posCor_-21_Mr_-20_wp', 3: '_new_posCor_-22_Mr_-21_wp', \
#                              4:'_CUT1_Contreras+13_mstar_'+galaxy_type+'_wp', 5: '_CUT2_Contreras+13_mstar_'+galaxy_type+'_wp', 6: '_CUT3_Contreras+13_mstar_'+galaxy_type+'_wp',\
#                              7:'_CUT3_Contreras+13_mstar_non-orphans_wp', 8: '_CUT3_Contreras+13_mstar_non-orphans_wp', 9: '_CUT3_Contreras+13_mcold_'+galaxy_type+'_wp',\
#                              10:'_CUT3_Contreras+13_mcold_non-orphans_wp', 11: '_CUT3_Contreras+13_sfr_'+galaxy_type+'_wp', 12:'_CUT3_Contreras+13_sfr_non-orphans_wp'}

                        #optimal Galacticus
#                        band={0:'_CUT1_Contreras+13_mstar_'+galaxy_type+'_wp', 4: '_CUT1_Contreras+13_mcold_'+galaxy_type+'_wp', 2: '_CUT3_Contreras+13_mcold_'+galaxy_type+'_wp',\
#                              1:'_CUT1_Contreras+13_sfr_'+galaxy_type+'_wp', 3: '_CUT3_Contreras+13_sfr_'+galaxy_type+'_wp', 5: '_CUT3_Contreras+13_MAB_dA_total_r_all_wp'} 
#
#                        band={0:'_CUT2_mstar_no_wp', 1: '_CUT2_mcold_no_wp', 2: '_CUT3_mcold_no_wp',\
#                              3:'_CUT3_mcold_5e8_no_wp', 4: '_CUT1_sfr_no_wp', 5: '_CUT3_Contreras+13_MAB_dA_total_r_test5_no_0.1_200_150_wp'} 

                        #CUT3 mcold
#                        band={0:'_CUT3_Contreras+13_mcold_highZcold-highMstar_wp', 1: '_CUT3_Contreras+13_mcold_highZcold-lowMstar_wp', 2: '_CUT3_Contreras+13_mcold_lowZcold-highMstar_wp',\
#                              3:'_CUT3_Contreras+13_mcold_lowZcold-lowMstar_wp', 4: '_CUT3_mcold_no_wp', 5: '_CUT3_Contreras+13_mcold_centrals_wp'} 


#                        band={0:'_CUT3_Contreras+13_MAB_dA_total_r_test_'+galaxy_type+'_0.1_200_150_wp', 1: '_CUT3_Contreras+13_MAB_dA_total_r_test2_'+galaxy_type+'_0.1_200_150_wp', 2: '_CUT3_Contreras+13_MAB_dA_total_r_test3_'+galaxy_type+'_0.1_200_150_wp',\
#                              3:'_CUT3_Contreras+13_MAB_dA_total_r_test4_'+galaxy_type+'_0.1_200_150_wp', 4: '_CUT3_Contreras+13_MAB_dA_total_r_test5_'+galaxy_type+'_0.1_200_150_wp', 5: '_CUT3_Contreras+13_MAB_dA_total_r_all_wp'} 
#

                        #optimal SAG
#                        band={0:'_CUT3_Contreras+13_mstar_non-orphans_wp', 1: '_CUT3_Contreras+13_mstar_centrals_wp', 2: '_CUT3_Contreras+13_mcold_all_wp',\
#                              4:'_CUT3_Contreras+13_sfr_all_wp', 3: '_CUT1_Contreras+13_mstar_centrals_wp', 5: '_CUT3_Contreras+13_MAB_dA_total_r_all_wp'}
#                       
                        #optimal SAGE
#                        band={0:'_CUT3_Contreras+13_mstar_non-orphans_wp', 2: '_CUT3_Contreras+13_mcold_disk_non-orphans_wp', 3: '_CUT3_Contreras+13_sfr_centrals_wp',\
#                              1:'_CUT1_Contreras+13_mstar_centrals_wp', 4: '_CUT3_Contreras+13_mstar_centrals_wp', 5: '_CUT3_Contreras+13_mcold_disk_centrals_wp'}
#                    #print band

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
                                self.print_redshift='z$\sim$0.55\n'

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
                                    #filename =mycomp+'anaconda/pro/myRun/histos/twoPCF/'+model+'/twoPCF_'+model+'_'+boxsize+'_z_0.56_tarsel_CMASS'+band[i]+'_'+plot[0:2]+'.txt'
                                    #filename =mycomp+'anaconda/pro/myRun/histos/twoPCF_save/twoPCF_Galacticus_1Gpc_z_0.09_tarsel'+band[i]+'_'+plot[0:2]+'.txt'
                                    filename =mycomp+'anaconda/pro/myRun/histos/twoPCF/twoPCF_'+band[i]+'_'+plot[0:2]+'.txt'
#                                    else:
                                    #filename =mycomp+'anaconda/pro/myRun/histos/twoPCF/SAG_CMASS/twoPCF_SAG_1Gpc_z_0.56_tarsel_new_CMASS'+band[i]+'_'+plot[0:2]+'.txt'

                                if SFH==True:
                                    filename=my_custom_filename

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
                                filename =mycomp+'anaconda/pro/myRun/histos/twoPCF_save/twoPCF_Galacticus_1Gpc_z_0.09_tarsel'+band[i]+'.txt'
                                   
                                mydata = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', nr_col=6, nr_rows=int(self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_twoPCF_nbins']), data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
                                mydata[:,[0,1]]=mydata[:,[2,3]]
                                mydata[:,1]*=mydata[:,0]
                                if self.myconfig_array['catname'+str(i)].find('Galacticus')!=-1:
                                    self.print_redshift='Galacticus                (a)\n'
                                elif self.myconfig_array['catname'+str(i)].find('SAG_')!=-1:
                                    self.print_redshift='SAG                         (b)\n'                                    
                                elif self.myconfig_array['catname'+str(i)].find('SAGE')!=-1:
                                    self.print_redshift='SAGE                       (c)\n'
                                    
                                key='22-21'
                                #pi_max=', $\pi_{max}$: 60.0 Mpc'
                                if key=='22-21':
                                    self.print_redshift='Galacticus\nz=0.1\n$M_{Cold}$-CUT3'
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
                        'SFH': SFH,
                        'zevol': plotZevol
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
        
        while self.b<self.plot_map_array['nr_plot_keys']:
            
            obs_data_array = [[],[]]
            self.subplot_loops=0

            plot_key = self.plot_map_array['plot_map_id'+str(self.b)]
            
            if plot_key!='sfr2z' and plot_key!='cSFRD':
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
                   
                    if self.load_from_file=='False' and plot_key.find('analyseTargetSelection')==-1: caseSwitcherPlotKey('mainCalculate', None, None, None, None, None, None)
                                       
                    elif plot_key.find('analyseTargetSelection')!=-1: caseSwitcherPlotKey('analyseTargetSelection', None, None, None, None, None, None)
                     
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
                                mymethod='M1' 
                                for reds in [0.0,0.1,0.56,0.7,3.0,4.15]:
                                    mL.starforming_cut_Henriques20_A1(redshift=reds)
                                    
                                exit()
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
                                #myprop = ['_SHMF', '_rhalfmass']#, '_mstar']#, '_cgf', '_mcold', '_Mzgas', '_mbh', '_rbulgevsrdisk', '_rhalfmass']
                                #myprop = ['_sfr', '_r-i', '_zcold', 'SFHd']
                                myprop = ['_mhalo']                                 
                                # selected which plot should be created
                                #   default: SFH for each selected property in 'myprop'
                                #   gr:      growth histories of all properties 'in myprop'
                                #   gr_one:  growth histories of each property in 'myprops', but for all subsamples in 'mysamples'
                                #   histo:   histrogram of each property in 'myprops' and all subsamples in 'mysamples'
                                #   wp:      2pCF of each subsample in 'mysamples' and at each redshift defind in 'band'
                                #   wp_one:  2pCF of all subsamples in 'mysamples' together in one plot at each redshift defind in 'band'
                                #   print_stats:  print stats of galaxy properts of all subsamples in 'mysamples' at each redshift defind in 'band' usind the function mL.test_methods()                               
                                #------------------------------------------                                
                                workflow = ['default', 'gr', 'gr_one', 'histo','wp','wp_one', 'stats_N', 'stats_xbar_M1-M2_M2', 'stats_xbar_M1found_M1', 'stats_xbar_M1found_M2', 'stats_xbar_MAD']
                                workflow='stats_xbar_MAD'
                                
                                if workflow.find('histo')!=-1:
                                    if self.myconfig_array['catname0'].find('run2')!=-1:                                    
                                        myprop = ['_mhalo', '_mstar', '_sfr', '_ssfr', '_zcold', '_g-i',\
                                                  '_mbh', '_rbulgevsrdisk', '_vmax', '_vdisp',\
                                                  '_mean_age_stars_disk', '_mean_age_stars_spheroid']
                                    else:
                                        myprop = ['_mhalo', '_mstar', '_sfr', '_ssfr', '_zcold', '_g-i',\
                                                  '_mbh', '_rbulgevsrdisk', '_cgf']                                        

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
                                band = ['0.56', '0.56']                                  
                                custom_plot_filename_prefix='_CMASS_SFH_down3_'
                                for_paper=''
                                
                                ######################################################################################################


                                if myplot_num=='_demo':
                                    mysample=['red', 'red', 'red', 'blue', 'blue', 'blue']                                
                                elif mymethod=='M1':
                                    #Galacticus
                                    #mysample=['', 'low', 'high', 'passive', 'active', 'red', 'blue']                                    
                                    mysample=['low', 'low', 'low', 'passive', 'low', 'red', 'low', 'low-zcold', 'high-zcold']
                                    #mysample=['low', 'passive', 'red', 'low-zcold', 'high-zcold']
                                    #mysample=['low','low']
                                elif mymethod=='M2' and myplot_num=='':       
                                    #Galacticus
                                    #standard plot1 M2 & M2
                                    mysample=['low', 'low', 'high', 'passive', 'active', 'red', 'blue', 'low-zcold', 'high-zcold', 'mstar11gt', 'zcold9st', 'zcold9gt', 'redmhalo13', 'lowmhalo13']
                                    #mysample=['', 'low', 'high', 'passive', 'active', 'red', 'blue', 'zcold9gt', 'zcold9st']                                     
                #                              'mstar11gt','mhalo12gt', 'mstar10st', 'mhalo12st', 'zcold9gt', 'zcold9st', 'redmstar11', 'redmhalo13', 'lowmstar11', 'lowmhalo13']                                   
                                    #mysample=['', 'red', 'blue', 'low-zcold', 'high-zcold','redmhalo13']
                                    #mysample=['', 'passive', 'red', 'low-zcold', 'high-zcold', 'zcold9st', 'zcold9gt', 'redmhalo13' ]
                                elif mymethod=='M2' and myplot_num=='_2':       
                                    #Galacticus
                                    #standard plot2 only M2
                                    mysample=['', 'mstar11gt', 'mhalo12gt', 'zcold9gt', 'zcold9st', 'redmstar11', 'redmhalo13', 'lowmstar11', 'lowmhalo13']
                                    #reduced sample for plots in paper
                                    mysample=['', 'mstar11gt', 'zcold9gt', 'zcold9st', 'redmhalo13', 'lowmhalo13']
                                    
                                elif mymethod=='_SDSS':
                                    mysample=['highZcold-highMstar','lowZcold-highMstar']
                                                                        
                                else:
                                    mysample=['','low', 'high', 'passive', 'active', 'red', 'blue']

                                #mysample=['low-zcold']#, 'low-zcold', 'low-zcold', 'low-zcold']

                                if workflow=='zevol':
                                    for sample in ['mstar']:#, 'sfr', 'ssfr']:
                                        for prop in ['_mstar']:
                                            print sample, prop, mysample, mymethod
                                            self.myPlot=caseSwitcherPlotKey('loadFromFile', sample, 'zevol', prop, mysample, mymethod, myplot_num)
                                            mycustom_plot_filename='HOD-SF_zevolv_'+sample+'_'+prop+for_paper
                                            plotOutput(prop, mycustom_plot_filename)  

                                elif workflow=='gr' or workflow=='gr_res':
                                    #plot property for all samples in one plot
                                    for prop in ['_mhalo', '_mstar', '_SHMF', '_mcold', '_Mzgas' , '_cgf', '_mbh', '_rbulgevsrdisk', '_rhalfmass'][::-1]:
                                        self.myPlot=caseSwitcherPlotKey('loadFromFile', '', workflow, prop, mysample, mymethod, myplot_num)
                                        mycustom_plot_filename=custom_plot_filename_prefix+mymethod+prop+'_'+workflow+myplot_num+for_paper
                                        plotOutput('_'+workflow, mycustom_plot_filename)                                       

                                elif workflow=='gr_one':
                                    #3 plot growth of one sample
                                    for sample in ['low-zcold']:#mysample:
                                        self.myPlot=caseSwitcherPlotKey('loadFromFile', sample, workflow, '', mysample, mymethod, myplot_num)
                                        mycustom_plot_filename=custom_plot_filename_prefix+sample+'_'+mymethod+'_'+workflow+for_paper
                                        plotOutput('_'+workflow, mycustom_plot_filename)
                                        
                                elif workflow.find('histo')!=-1:
                                    #4 plot histo of all samples                                    
                                    for sample in mysample:
                                        for prop in myprop:
                                            self.myPlot=caseSwitcherPlotKey('loadFromFile', sample, workflow, prop, mysample, mymethod, myplot_num)
                                            mycustom_plot_filename=custom_plot_filename_prefix+mymethod+'_'+sample+prop+'_'+workflow+for_paper
                                            plotOutput(prop+'_'+workflow, mycustom_plot_filename)

                                elif workflow=='stats_calc':
                                    #5 print statistics of properties in a text-file
                                    myprops=['mhalo', 'SHMR', 'sfr', 'ssfr', 'mcold', 'Mzgas', 'zcold', 'g-i']
                                    for prop in myprops:                                                            
                                        for sample in mysample:
                                            for count, redshift in enumerate(band):                     
                                                print 'count:', count, 'redshift:', redshift, 'sample:', sample, 'prop:', prop
                                                if count==0:
                                                    stats_data = np.zeros((55, 23), dtype=np.float64)  
                                                myfilename=mycomp+'anaconda/pro/myRun/histos/sfr2z/Gal-dens_main_cents_'+mymethod+'/stats/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_SFH_z-evolution_stats_'+sample+'_'+prop+'.txt'

                                                stats_data, myheader=mL.calc_residuals(stats_data, count, redshift, sample, prop)
                                                                 
                                                myOutput.writeIntoFile(myfilename,
                                                                       stats_data,
                                                                       myheader='SF- project statistics of galaxy properties with redshift; cat: '+self.myconfig_array['catname'+str(self.myPipe.a)]+' CMASS-sample: '+sample+myheader,
                                                                       data_format='%0.3f\t%0.5f\t%0.5f\t%0.5f\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%i\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%i\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%i\t%0.5f\t%0.5f\t%0.5f\t%0.5f')

                                elif workflow=='print_stats':
                                    #5 print statistics of properties in a text-file
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
                                                print 'sample:', sample, 'prop:', prop, 'redshift:', redshift
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
        #                                                       [str(prop)+'\t'+str(ratio)],
        #                                                       myheader='SF-project z-evolution cat: '+self.myconfig_array['catname'+str(self.myPipe.a)]+' sample:'+mysample,
        #                                                       append_mytext=False,
        #                                                       data_is_string=False,
        #                                                       data_format='%s')
                                                   
                                                else:
                                                    myOutput.writeIntoFile(myfilename,
                                                               string+'\n',
                                                               append_mytext=True,
                                                               data_is_string=True,
                                                               data_format='%s')
                                                    
        #                                            myOutput.writeIntoFile(mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_SFH_z-evolution_'+mysample+'_found_M1.txt',
        #                                                       string_found+'\n',
        #                                                       append_mytext=True,
        #                                                       data_is_string=True,
        #                                                       data_format='%s')                                            
        #
#                                                    myOutput.writeIntoFile(mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_SFH_z-evolution_'+mysample+'_frac_found.txt',
#                                                               str(prop)+'\t'+str(ratio)+'\n',
#                                                               append_mytext=True,
#                                                               data_is_string=True,
#                                                               data_format='%s')

 

                                            myOutput.writeIntoFile(mycomp+'anaconda/pro/data/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'/'+self.myconfig_array['catname'+str(self.myPipe.a)]+'_SFH_z-evolution_bad_mstar_count.txt',
                                                       str(sample)+'\t'+str(bad_mstar_count)+'\n',
                                                       append_mytext=True,
                                                       data_is_string=True,
                                                       data_format='%s')                                           

                                elif workflow.find('stats_violin')!=-1:
                                    #1 plot property for all samples in one plot
                                    for z in band: 
                                        for prop in ['_mhalo', '_mstar', '_SHMR']:# '_zcold', '_SHMR', '_mstar', '_mcold', '_Mzgas', '_sfr', '_ssfr', '_g-i'][::-1]:
                                                self.myPlot=caseSwitcherPlotKey('loadFromFile', z , workflow, prop, mysample, mymethod, myplot_num)
                                                mycustom_plot_filename=custom_plot_filename_prefix+mymethod+prop+'_z_'+z[0]+'_'+workflow+myplot_num+for_paper
                                                plotOutput('_'+workflow+prop, mycustom_plot_filename)

                                elif workflow.find('stats')!=-1:
                                    #1 plot property for all samples in one plot
                                    for prop in ['_mhalo']:#, '_mstar', '_SHMR']:# '_zcold', '_SHMR', '_mstar', '_mcold', '_Mzgas', '_sfr', '_ssfr', '_g-i'][::-1]:
                                            self.myPlot=caseSwitcherPlotKey('loadFromFile', '' , workflow, prop, mysample, mymethod, myplot_num)
                                            mycustom_plot_filename=custom_plot_filename_prefix+mymethod+prop+'_'+workflow+myplot_num+for_paper
                                            plotOutput('_'+workflow, mycustom_plot_filename)                                                
                                                
                                elif workflow=='wp':
                                #4 plot wp of all samples of certain redshift as prop
                                    for prop in band:                                             
                                        self.myPlot=caseSwitcherPlotKey('loadFromFile', '', 'wp', prop, mysample, mymethod, myplot_num)
                                        mycustom_plot_filename=custom_plot_filename_prefix+mymethod+'_z_'+prop+'_wp'+myplot_num+for_paper
                                        plotOutput('_wp', mycustom_plot_filename)
                                       
                                elif workflow=='wp_one':
                                    #4 plot 2pCF of all samples
                                    for sample in mysample:
                                        self.myPlot=caseSwitcherPlotKey('loadFromFile', sample, 'wp_one', band, mysample, mymethod, myplot_num)
                                        mycustom_plot_filename=custom_plot_filename_prefix+mymethod+sample+'_wp_one'
                                        plotOutput('_wp_one', mycustom_plot_filename)
                                elif workflow=='envr_props':
                                    #4 plot 2pCF of all samples
                                    self.myPlot=caseSwitcherPlotKey('loadFromFile', '', 'envr_props', '', mysample, mymethod, myplot_num)
                                    mycustom_plot_filename=custom_plot_filename_prefix+mymethod+'_envr_props'
                                    plotOutput('_envr_props', mycustom_plot_filename)                                        
                                else:
                                    #1 plot property for all samples in one plot
                                    for prop in myprop:
                                        #sample, key, prop
                                        self.myPlot=caseSwitcherPlotKey('loadFromFile', '', 'SFH', prop, mysample, mymethod, myplot_num)
                                        mycustom_plot_filename=custom_plot_filename_prefix+mymethod+prop+myplot_num+errors+for_paper
                                        plotOutput(prop, mycustom_plot_filename)                                   
                                exit()
                                    
                                
                            else:
                                self.myPlot=caseSwitcherPlotKey('loadFromFile', None, None, None, None, None, None)                                
                                
                                

                        
                     
                    self.myPipe.a+=1
                                   
                self.myPipe.i+=1

                if self.load_from_file=='True' and self.plot_custom_loop=='False': plotOutput(None, None)
                
            self.b+=1
            
######## END MAIN       ##################################################################################################################################### 