# Load packages
import numpy as np
import config as conf
import outputData as oD
import arangeData as aD
import myFuncs as mF
import myLib as mL

myOutput = oD.OutputData(config=' ')
myConfig = conf.Configuration()
myData = aD.ArangeData()
myFuncs = mF.MyFunctions()

import os
system_info = os.getcwd()
start = system_info.find('anaconda')
mycomp = system_info[:start]


def angM0(load_from_file=False):
    
    obs_data_array = np.zeros((10, 7), dtype=np.float32)

    x_data=np.geomspace(9,13,10)
    obs_data_array[:,0]=x_data  
    obs_data_array[:,1]=2/3.0*x_data-3
    
    #print obs_data_array[:,[0,1]]

    return obs_data_array, 'd=-3, $M^{2/3}$', False   

def angM(load_from_file=False):
    
    obs_data_array = np.zeros((10, 7), dtype=np.float32)

    x_data=np.geomspace(9,13,10)
    obs_data_array[:,0]=x_data  
    obs_data_array[:,1]=2/3.0*x_data-4
    
    #print obs_data_array[:,[0,1]]

    return obs_data_array, 'd=-4, $M^{2/3}$', False     

def angM2(load_from_file=False):
    
    obs_data_array = np.zeros((10, 7), dtype=np.float32)

    x_data=np.geomspace(9,13,10)
    obs_data_array[:,0]=x_data  
    obs_data_array[:,1]=2/3.0*x_data-5
    
    #print obs_data_array[:,[0,1]]

    return obs_data_array, 'd=-5, $M^{2/3}$', False  

def Baldry_08_SMF(load_from_file=False):
    
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Baldry+08_from_adam.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=6)

    obs_data_array = np.zeros((data[:,0].size,6), dtype=np.float)
    
    obs_data_array[:,0] = data[:,0]
    obs_data_array[:,1] = np.log10(data[:,1])
    obs_data_array[:,3] = data[:,4]
    obs_data_array[:,4] = data[:,5]
    
    return obs_data_array, 'Baldry+08', False

def Baldry12(load_from_file=False):
    
    obs_data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Baldry_2012_SMF_z0.dat', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=2)
    
    obs_data_array[:,0] = obs_data_array[:,0]-2*np.log10(0.677)
    obs_data_array[:,1] = np.log10(obs_data_array[:,2]/1000)# -2*np.log10(0.7)
    obs_data_array[:,4] = 10**(np.log10(obs_data_array[:,3]/1000))#-2*np.log10(0.7))
    obs_data_array[:,5] = obs_data_array[:,3]

    return obs_data_array, 'Baldry+12', False

def Behroozi10(load_from_file=False):
   
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Behroozi+10_mstar2mhalo_Table3.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=6)  
    obs_data_array = np.zeros((data[:,0].size,7), dtype=np.float) 

    #z=0.1    
    obs_data_array[:,0] = data[:,0]
    obs_data_array[:,1] = data[:,1]
    obs_data_array[:,4] = data[:,2]
    obs_data_array[:,5] = data[:,3]
 
    #print obs_data_array

    return obs_data_array, 'Behroozi+10', False

def DEEP2_FireFly(load_from_file=True):

    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/DEEP2-FireFly.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)    
    
    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float)
    
    obs_data_array[:,0] = data_array[:,0]           
    obs_data_array[:,1] = data_array[:,1]
    obs_data_array[:,4] = data_array[:,1]-data_array[:,2]
    obs_data_array[:,5] = data_array[:,3]-data_array[:,1]

    print obs_data_array
    
    return obs_data_array, 'DEEP2', False  

def Maudau14_sfr2z(load_from_file=False):
    
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Madau+14_SFH_Fig13_from_Driver+18.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=2)
    
    obs_data_array = np.zeros((data[:,0].size,7), dtype=np.float) 
    x_err = np.zeros((data[:,0].size,), dtype=np.float)
    z = np.zeros((data[:,0].size,), dtype=np.float) 
    
    from cosmolopy import cparam, cd
    import cosmolopy.constants as cc
    
    cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
    cosmo_WMAP = cd.set_omega_k_0(cosmo)
   

    for k in enumerate(data[:,0]):
        #print 'k:', k[0], k[1], ' ',
        z[k[0]] = cd.redshift_d_light(k[1] * cc.c_light_Mpc_Gyr, **cosmo)
        #print obs_data_array[k[0],2]
               
    for k in enumerate(data[:,4]):
        #print 'k:', k[0], k[1], ' ',
        x_err[k[0]] = cd.redshift_d_light(k[1] * cc.c_light_Mpc_Gyr, **cosmo)
        print obs_data_array[k[0],0], x_err[k[0]]

    Vc_WMAP = cd.comoving_volume(z, **cosmo_WMAP)
     
    cosmo_Planck = cparam.Planck(flat=True, extras=True)     
    Vc_Planck = cd.comoving_volume(z, **cosmo_Planck)
     
    dVc=Vc_WMAP/Vc_Planck

    obs_data_array[:,0] = 1.0+z
    obs_data_array[:,1] = np.log10(data[:,1])-np.log10(dVc)
    obs_data_array[:,2] = x_err-z
    obs_data_array[:,3] = obs_data_array[:,2]
    obs_data_array[:,4] = (data[:,1]-data[:,3])/data[:,1]*0.434
    obs_data_array[:,5] = (data[:,2]-data[:,1])/data[:,1]*0.434

    return obs_data_array, 'Madau & Dickinsion (2014)', False

def Driver18_sfr2z(load_from_file=False):
    
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Driver+18_SFH_Fig13_GAMA.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=2)

    obs_data_array = np.zeros((data[:,0].size,7), dtype=np.float) 
    
    from cosmolopy import cparam, cd
    import cosmolopy.constants as cc
    
    cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
    cosmo_WMAP = cd.set_omega_k_0(cosmo)
    #print cosmo
    
    for k in enumerate(data[:,0]):
        #print 'k:', k[0], k[1], ' ',
        obs_data_array[k[0],0] = cd.redshift_d_light(k[1] * cc.c_light_Mpc_Gyr, **cosmo)
        #print obs_data_array[k[0],0]

    Vc_WMAP = cd.comoving_volume(obs_data_array[:,0], **cosmo_WMAP)
     
    cosmo_Planck = cparam.Planck(flat=True, extras=True)     
    Vc_Planck = cd.comoving_volume(obs_data_array[:,0], **cosmo_Planck)
     
    dVc=Vc_WMAP/Vc_Planck  

    obs_data_array[:,0]+=1.0
    obs_data_array[:,1] = np.log10(data[:,1])-np.log10(dVc)
    obs_data_array[:,4] = (data[:,2]-data[:,1])/data[:,1]*0.434
    obs_data_array[:,5] = (data[:,2]-data[:,1])/data[:,1]*0.434

    return obs_data_array, 'Driver et al. (2018)', False

def Behroozi13_sfr2z(load_from_file=False):
    
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Behroozi+13_csfrs_new.dat', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=3)
    
    obs_data_array = np.zeros((data[:,0].size,7), dtype=np.float) 
    
    obs_data_array[:,0] = data[:,0]+1.0
    obs_data_array[:,1] = data[:,1]
    obs_data_array[:,4] = data[:,2]
    obs_data_array[:,5] = data[:,3]
    print obs_data_array
    return obs_data_array, 'Behroozi et al. (2013)', False


def Behroozi13c(load_from_file=False):
    
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/c_smmr_z0.55.boss.dat.converted.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=1)
    #myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Behroozi_2013c_Mstar2Mhalo_z0.1_Fig14.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=7)
    
    obs_data_array = np.zeros((data[:,0].size,7), dtype=np.float) 

    obs_data_array[:,0] = data[:,0]
    obs_data_array[:,1] = data[:,1] + 0.07
    obs_data_array[:,4] = obs_data_array[:,1] + data[:,3]
    obs_data_array[:,5] = obs_data_array[:,1] - data[:,2]    

    #print obs_data_array

    return obs_data_array, 'Behroozi+13', False


def Behroozi13cINRodriguez15(load_from_file=False):
    
    obs_data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Rodriguez_2015_Mstar2Mhalo_Fig13.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=7)
    
    obs_data_array[:,0] = obs_data_array[:,0]
    obs_data_array[:,1] = obs_data_array[:,1]
    obs_data_array[:,4] = obs_data_array[:,1]-obs_data_array[:,4]
    obs_data_array[:,5] = obs_data_array[:,1]+obs_data_array[:,5]

    return obs_data_array, 'Behroozi+13 modified (z~0.1)', False

def Bernadi16_LF_CMASS(load_from_file=False):
        
    filename=mycomp+'anaconda/pro/OBS/Bernadi+16_Fig1_LF_CMASS_all_i-band_observed.txt'        
    data = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)

    filename=mycomp+'anaconda/pro/OBS/Bernadi+16_Fig1_LF_CMASS_all_i-band_observed_+dy.txt'        
    plus_dy = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)
    
    
    obs_data_array = np.zeros((plus_dy[:,0].size,7), dtype=np.float32)
     
    print data[:,0]
    
    obs_data_array[:,0] = data[:,0]+5*np.log10(0.7)-5*np.log10(0.6777)-(-2.5*np.log10(1.0+0.55))
    obs_data_array[:,1] = data[:,1]
    obs_data_array[:,4] = abs(10**data[:,1]-10**plus_dy[:,1])
    obs_data_array[:,5] = obs_data_array[:,4]  
    print obs_data_array[:,[0,1,4,5]]
    
    return obs_data_array, 'Bernadi+16', False

def Boselli14_Mcold2Mstar(load_from_file=False):
    
    #data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Boselli+14_Mgas2Mstar.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Boselli+14_Mgas2Mstar_Fig5a_new.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)
    obs_data_array = np.zeros((data[:,0].size,7), dtype=np.float)     
    
    obs_data_array[:,0] = data[:,0]
    obs_data_array[:,1] = data[:,1] 
    obs_data_array[:,4] = data[:,4]
    obs_data_array[:,5] = data[:,5]

    return obs_data_array, 'Boselli+14', False

def CARNage_cold_gas_fraction_Peeples11(load_from_file=True):
        
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Peeples11_mcold2mstar_Fig2.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=6)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)

    obs_data_array[:,0] = data_array[:,0]
    obs_data_array[:,1] = data_array[:,1] 
    obs_data_array[:,4] = data_array[:,4]
    obs_data_array[:,5] = data_array[:,5]
          
    #print obs_data_array        
         
    return obs_data_array, 'Peeples & Shankar 11', False

def CARNage_sfr_Gruppioni15(load_from_file=False):
                
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/CARNage/CARNage_sfr_Gruppioni15.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=2)

    obs_data_array = np.zeros((data[:,0].size, 6), dtype=np.float)

    obs_data_array[:,0] = np.log10((10**data[:,0]+10**data[:,1])/2.0)      
    obs_data_array[:,1] = data[:,2]

    #shaded region corresponding to error bars
    obs_data_array[:,4] = data[:,2] + data[:,3]
    obs_data_array[:,5] = data[:,2] - data[:,3]

    #error bars
    obs_data_array[:,4] = data[:,3]
    obs_data_array[:,5] = data[:,3]
     
    
    return obs_data_array, 'Gruppioni+15', False

def CARNage_mbh_mbulge_Kormendy13(load_from_file=False):

    obs_data_array = np.zeros((87,5), dtype=np.float)

    if load_from_file ==True:
        obs_data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Kormendy&Ho_mbh_mbulge_corrected.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)

    else:  
        
        data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/CARNage/CARNage_mbh_mbulge_Kormendy13.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=3)
    
        obs_data_array[:,0] = data_array[:,0]          
        obs_data_array[:,1] = data_array[:,2]
        obs_data_array[:,2] = data_array[:,1]
        obs_data_array[:,3] = data_array[:,2] + (data_array[:,3]/10**data_array[:,1])*0.43
        obs_data_array[:,4] = data_array[:,2] - (data_array[:,4]/10**data_array[:,1])*0.43

        #print legend_obs        
        #print obs_data_array        
        
        filename_out = mycomp+'anaconda/pro/OBS/Kormendy&Ho_mbh_mbulge_corrected.txt'
        myOutput.writeIntoFile(filename_out, 
                               obs_data_array,
                               myheader='Kormendy&Ho (2013) arXiv:1304.7762v1 Table 2+3\n(1) log mstarsph [Msun] (2) mbh [Msun] (3) x-err (-) (4) y-err (-)  (4) y-err (+)',
                               data_format="%0.8e",
                               mydelimiter='\t')   
    
    return obs_data_array, 'Kormendy&Ho 13', False

def CARNage_mbh_mbulge_McConnell13(load_from_file=False):

    obs_data_array = np.zeros((36,5), dtype=np.float)

    if load_from_file ==True:
        obs_data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/McConell&Ma_mbh_mbulge_corrected.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)

    else:  
        
        data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/CARNage/CARNage_mbh_mbulge_McConnell13.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=3)
    
        obs_data_array[:,0] = data_array[:,0]          
        obs_data_array[:,1] = data_array[:,2]
        obs_data_array[:,2] = data_array[:,1]
        obs_data_array[:,3] = data_array[:,2] + (data_array[:,3]/10**data_array[:,1])*0.43
        obs_data_array[:,4] = data_array[:,2] - (data_array[:,4]/10**data_array[:,1])*0.43

        #print legend_obs        
        #print obs_data_array        
        
        filename_out = mycomp+'anaconda/pro/OBS/CARNage/CARNage_mbh_mbulge_McConnell13_corrected.txt'
        myOutput.writeIntoFile(filename_out, 
                               obs_data_array,
                               myheader='McConnell & Ma, 2013, ApJ, 764, 184 MASSIVE survey http://blackhole.berkeley.edu/wp-content/uploads/2016/09/current_ascii.txt\n(1) log mstarsph [Msun] (2) mbh [Msun] (3) x-err (-) (4) y-err (-)  (4) y-err (+)',
                               data_format="%0.8e",
                               mydelimiter='\t')  
    
    
    return obs_data_array, 'McConnell&Ma 13', False

def CarnegieOBS_z0(load_from_file=True):

    obs_data_array = np.zeros((17,5), dtype=np.float)

    if load_from_file ==True:
        obs_data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/CARNage_z0.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=1)

    else:  
        
        data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/CalibrationData/DATA/CARNage_z0.dat', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=4)
#        old before log-scaling    
#        obs_data_array[:,0] = 10**(np.log10((10**obs_data_array[:,0]+10**obs_data_array[:,1])/2.0)-2*np.log10(0.7) + np.log10(carnage_small_h))           
#        obs_data_array[:,1] = 10**(np.log10(obs_data_array[:,2])-2*np.log10(0.7)+np.log10(carnage_small_h))
#        obs_data_array[:,3] = abs(obs_data_array[:,1] - np.log10(10**obs_data_array[:,1]-(obs_data_array[:,3]/(2*0.7)))) + np.log10(carnage_small_h)
#        obs_data_array[:,4] = obs_data_array[:,3]

        obs_data_array[:,0] = 10**(np.log10((10**data_array[:,0]+10**data_array[:,1])/2.0)-2*np.log10(0.7))           
        obs_data_array[:,1] = 10**(np.log10(data_array[:,2])+3*np.log10(0.7))#+np.log10(carnage_small_h))
        obs_data_array[:,3] = 10**(np.log10(data_array[:,3])+3*np.log10(0.7))# #*carnage_small_h
        obs_data_array[:,4] = obs_data_array[:,3]        
        #print obs_data_array     
 
        obs_data_array[:,0] = np.log10(((10**data_array[:,0]+10**data_array[:,1])/2.0)/(0.7**2))           
        obs_data_array[:,1] = np.log10(data_array[:,2]*(0.7**3))
        obs_data_array[:,3] = (data_array[:,2]*(0.7**3)/10**obs_data_array[:,1])*0.43
        obs_data_array[:,4] = (data_array[:,2]*(0.7**3)/10**obs_data_array[:,1])*0.43  


        
        filename_out = mycomp+'anaconda/pro/OBS/CARNage_z0.txt'
        myOutput.writeIntoFile(filename_out, 
                               obs_data_array,
                               myheader='CARNage_z0.000.dat corrected 13/11/2015 (1) log Mstar [Msun], (2) log Phi dn/dlog M [Msun Mpc^-3], (3) log Phi error (-),   (4) log Phi error (+)',
                               data_format="%0.5f",
                               mydelimiter='\t')
    
    
    return obs_data_array, 'Baldry08+12, Li&White09', False
    
def CarnegieOBS_z2(load_from_file=True):

    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/CalibrationData/DATA/CARNage_z2.dat', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=4)    
    
    obs_data_array = np.zeros((data_array[:,0].size,5), dtype=np.float)
    
    obs_data_array[:,0] = 10**(np.log10((10**data_array[:,0]+10**data_array[:,1])/2.0)-2*np.log10(0.7))           
    obs_data_array[:,1] = 10**(np.log10(data_array[:,2])+3*np.log10(0.7))#+np.log10(carnage_small_h))
    obs_data_array[:,3] = 10**(np.log10(data_array[:,3])+3*np.log10(0.7))# #*carnage_small_h
    obs_data_array[:,4] = obs_data_array[:,3] 

    #print legend_obs
    print obs_data_array    

    filename_out = mycomp+'anaconda/pro/OBS/CARNage_z2.txt'
    myOutput.writeIntoFile(filename_out, 
                           obs_data_array,
                           myheader='CARNage_z2.000.dat no little h in there (1) Mstar [Msun], (2) Phi [Msun Mpc^-3], (3) -,  (4) Phi error (-),   (5) Phi error (+)',
                           data_format="%0.5f",
                           mydelimiter='\t')
    
    return obs_data_array, 'San+11, Muzz+13, Ilb+13, Tomc+13 z=2.0 (CARNage_z2.dat)', False    


def Gruppioni15_z_010(load_from_file=False):
                
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/gruppioni_test.dat', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=3)

    obs_data_array = np.zeros((data[:,0].size, 7), dtype=np.float)

    obs_data_array[:,0] = np.log10((10**data[:,0]+10**data[:,1])/2.0)+2*np.log10(0.71)-2*np.log10(0.6777)      
    obs_data_array[:,1] = data[:,2]-3*np.log10(0.71)+3*np.log10(0.6777)

    #shaded region corresponding to error bars
    obs_data_array[:,4] = data[:,2]-3*np.log10(0.71)+3*np.log10(0.6777) + data[:,3]
    obs_data_array[:,5] = data[:,2]-3*np.log10(0.71)+3*np.log10(0.6777) - data[:,3]

    #error bars
    obs_data_array[:,4] = data[:,3]
    obs_data_array[:,5] = data[:,3]
    
    return obs_data_array, 'Gruppioni+15', False

def Gruppioni15_z_050(load_from_file=False):
                
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/gruppioni_test.dat', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=3)

    obs_data_array = np.zeros((data[:,0].size, 7), dtype=np.float)
    size_one_set=data[:,0].size
    y_col=6

    obs_data_array[0:size_one_set,0] = np.log10((10**data[:,0]+10**data[:,1])/2.0)+2*np.log10(0.71)-2*np.log10(0.6777)      
    obs_data_array[0:size_one_set,1] = data[:,y_col]-3*np.log10(0.71)+3*np.log10(0.6777)

    #shaded region corresponding to error bars
    obs_data_array[0:size_one_set,4] = data[0:size_one_set,2] + data[0:size_one_set,y_col+1]
    obs_data_array[0:size_one_set,5] = data[0:size_one_set,2] - data[0:size_one_set,y_col+1]

    #error bars
    obs_data_array[0:size_one_set,4] = data[0:size_one_set,y_col+1]
    obs_data_array[0:size_one_set,5] = data[0:size_one_set,y_col+1]
    
    return obs_data_array, 'Gruppioni+15', False


def Gruppioni15_z_100(load_from_file=False):
                
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/gruppioni_test.dat', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=3)

    obs_data_array = np.zeros((data[:,0].size, 7), dtype=np.float)
    size_one_set=data[:,0].size
    y_col=12
    obs_data_array[0:size_one_set,0] = np.log10((10**data[:,0]+10**data[:,1])/2.0)+2*np.log10(0.71)-2*np.log10(0.6777)       
    obs_data_array[0:size_one_set,1] = data[0:size_one_set,y_col]-3*np.log10(0.71)+3*np.log10(0.6777)

    #shaded region corresponding to error bars
    obs_data_array[0:size_one_set,4] = data[:,y_col] + data[:,y_col+1]
    obs_data_array[0:size_one_set,5] = data[:,y_col] - data[:,y_col+1]

    #error bars
    obs_data_array[0:size_one_set,4] = data[0:size_one_set,y_col+1]
    obs_data_array[0:size_one_set,5] = data[0:size_one_set,y_col+1]

    
    return obs_data_array, 'Gruppioni+15', False

def Gruppioni15_z_200(load_from_file=False):
                
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/gruppioni_test.dat', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=3)

    obs_data_array = np.zeros((data[:,0].size*2, 7), dtype=np.float)
    size_one_set=data[:,0].size
    y_col=16

    obs_data_array[0:size_one_set,0] = np.log10((10**data[:,0]+10**data[:,1])/2.0)+2*np.log10(0.71)-2*np.log10(0.6777)       
    obs_data_array[0:size_one_set,1] = data[0:size_one_set,y_col]-3*np.log10(0.71)+3*np.log10(0.6777)

    #shaded region corresponding to error bars
    obs_data_array[0:size_one_set,4] = data[:,y_col] + data[:,y_col+1]
    obs_data_array[0:size_one_set,5] = data[:,y_col] - data[:,y_col+1]

    #error bars
    obs_data_array[0:size_one_set,4] = data[0:size_one_set,y_col+1]
    obs_data_array[0:size_one_set,5] = data[0:size_one_set,y_col+1]

#    y_col=16
#
#    obs_data_array[size_one_set:data[:,0].size,0] = np.log10((10**data[size_one_set:data[:,0].size,0]+10**data[size_one_set:data[:,0].size:,1])/2.0)      
#    obs_data_array[size_one_set:data[:,0].size,1] = data[size_one_set:data[:,0].size,y_col]
#
#    #shaded region corresponding to error bars
#    obs_data_array[size_one_set:data[:,0].size,4] = data[size_one_set:data[:,0].size,2] + data[size_one_set:data[:,0].size,y_col+1]
#    obs_data_array[size_one_set:data[:,0].size,5] = data[size_one_set:data[:,0].size,2] - data[size_one_set:data[:,0].size,y_col+1]
#
#    #error bars
#    obs_data_array[size_one_set:data[:,0].size,4] = data[size_one_set:data[:,0].size,y_col+1]
#    obs_data_array[size_one_set:data[:,0].size,5] = data[size_one_set:data[:,0].size,y_col+1]    
    
    print obs_data_array
    
    return obs_data_array, 'Gruppioni+15', False

def Gruppioni15_cSFRD(load_from_file=False):
                
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Grupponi+15_Table1_cSFRD_UV+IR.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=2)
    print data
    obs_data_array = np.zeros((data[:,0].size, 7), dtype=np.float)

    obs_data_array[:,0] = (data[:,0]+data[:,1])/2.0 +1     
    obs_data_array[:,1] = np.log10(data[:,2])

    #shaded region corresponding to error bars
    obs_data_array[:,4] = data[:,3]/data[:,2]*0.434
    obs_data_array[:,5] = data[:,3]/data[:,2]*0.434
    
    print obs_data_array
    
    return obs_data_array, 'Gruppioni et al. (2015)', False

def load_catalog(path,
                 name_x,
                 name_y,
                 name_add=False,
                 name_add2=False,
                 name_add3=False,
                 name_add4=False,
                 name_add5=False,):

    if name_add==False: name_add=name_x
    if name_add2==False: name_add2=name_x    
    if name_add3==False: name_add3=name_x
    if name_add4==False: name_add4=name_x
    if name_add5==False: name_add5=name_x
    
#    system_info = os.getcwd()
#    start = system_info.find('anaconda')
#    if start==-1:
#        mycomp='/store/erebos/doris/'
#    else:
#        mycomp = system_info[:start]+'/anaconda/pro/data/'
    
    name_add6='weight_tot'
    data = myData.readAnyFormat(config=False, 
                                mypath=path,
                                myfilename=path,
                                nr_rows=70000000, 
                                nr_col=2, 
                                data_format='HDF5', 
                                mydtype=np.float32,
                                id_col_array={'name0': name_x, 
                                              'name1': name_y,
                                              'name2': name_add,
                                              'name3': name_add2,
                                              'name4': name_add3,
                                              'name5': name_add4,
                                              'name6': name_add5,
                                              'name7': name_add6,
                                              name_x+'_col_id': 0, 
                                              name_y+'_col_id': 1,
                                              name_add+'_col_id': 2,
                                              name_add2+'_col_id': 3,
                                              name_add3+'_col_id': 4,
                                              name_add4+'_col_id': 5,
                                              name_add5+'_col_id': 6,
                                              name_add6+'_col_id': 7,                                                 
                                              'data_type0': np.float32, 
                                              'data_type1': np.float32,
                                              'data_type2': np.float32,
                                              'data_type3': np.float32,
                                              'data_type4': np.float32,
                                              'data_type5': np.float32,
                                              'data_type6': np.float32,
                                              'data_type7': np.float32,                                              
                                              'nr_entries': 8})
    print '\n############################################\nLOADING REFERENCES \ OBSERVATIONS!\n############################################\n'
    print 'name_x/y', name_x, '/', name_y, 'ngal:', data.shape, 'min/max x:', min(data[name_x]), '/', max(data[name_x]), 'min/max y:', min(data[name_y]), '/', max(data[name_y])

    return data

def adjust_data_for_plot(data,
                         key):
    
    def r_i_i(data):
        try:
            data['mAB_dA_total_r']-=data['mAB_dA_total_i']            
        except:
            data['MAB_dA_total_r']-=data['MAB_dA_total_i']
            data['MAB_dA_total_i']-=5*np.log10(0.6777)
        return data        
        
    def r_i_mstar(data):

        data['mstar']=np.log10(data['mstar'])       
        
        try:
            data['mAB_dA_total_r']-=data['mAB_dA_total_i']        
            data['mAB_dA_total_i']=data['mstar']
            
            mask=np.where(np.isfinite(data['mAB_dA_total_i']))
        except:
            data['MAB_dA_total_r']-=data['MAB_dA_total_i']        
            data['MAB_dA_total_i']=data['mstar']
            
            
            mask=np.where(np.isfinite(data['MAB_dA_total_i']))
        
        return data[mask[:][0]]        

    def u_r_r(data):

        data['MAB_dA_total_u']-=data['MAB_dA_total_r']
        data['MAB_dA_total_r']-=5*np.log10(0.6777)
        return data

    def g_i_i():

        data['mAB_dA_total_g']-=data['mAB_dA_total_i']
        
    def g_i_mstar(data):

        try:
            data['mAB_dA_total_g']-=data['mAB_dA_total_i']        
            data['mAB_dA_total_i']=np.log10(data['mstar'])
            
            mask=np.where(np.isfinite(data['mAB_dA_total_i']))            
        except:
            data['MAB_dA_total_g']-=data['MAB_dA_total_i']        
            data['MAB_dA_total_i']=np.log10(data['mstar'])
            
            mask=np.where(np.isfinite(data['MAB_dA_total_i']))

        return data[mask[:][0]]            
        

    def u_r_mstar():

        data['mAB_dA_total_u']-=data['mAB_dA_total_r']        
        data['mAB_dA_total_r']=np.log10(data['mstar'])

    def dmesa_i(data):
        return data        

    def dmesa_mstar(data):
        data['mstar']=np.log10(data['mstar'])
        return data
    
    def r_i_g_r(data):
        
        try:
            data['mAB_dA_total_g']-=data['mAB_dA_total_r']
            data['mAB_dA_total_r']-=data['mAB_dA_total_i'] 
        except:
            data['MAB_dA_total_g']-=data['MAB_dA_total_r']
            data['MAB_dA_total_r']-=data['MAB_dA_total_i']
        return data

    def rdisk_rbulge(data):

        data['rbulge']=np.log10(data['rbulge']*1e3)
        data['rdisk']=np.log10(data['rdisk']*1e3)
        mask=np.where(np.isfinite(data['rdisk']))
        
        return data[mask[:][0]]

    def rhalf_mstar(data):


        data['rhalf_mass']=np.log10(data['rhalf_mass']*1e3)
        mask=np.where(np.isfinite(data['rhalf_mass']))
       
        return data[mask[:][0]]

    def r(data):
    
        #data['rbulge']=np.log10(data['rbulge']*1e3)
        data['rdisk']=np.log10(data['rdisk']*1e3)
        #data['rhalf_mass']=np.log10(data['rhalf_mass']*1e3)
        data['mstar']=np.log10(data['mstar'])
        mask=np.where(np.isfinite(data['rdisk']))
       
        return data[mask[:][0]]

    def sfr_mstar(data):

        data['sfr']=np.log10(data['sfr'])        
        data['mstar']=np.log10(data['mstar'])
        mask=np.where(np.isfinite(data['sfr']))
            
        return data[mask[:][0]]

    def ssfr_mstar(data):

        try:
            data['ssfr']=np.log10(data['ssfr'])
            mask=np.where(np.isfinite(data['ssfr']))            
        except:
            data['sfr']=np.log10(data['sfr']/data['mstar'])
            mask=np.where(np.isfinite(data['sfr']))
            
        data['mstar']=np.log10(data['mstar'])
            
        return data[mask[:][0]]

    def ssfr_sfr(data):

        data['ssfr']=np.log10(data['ssfr'])        
        data['sfr']=np.log10(data['sfr'])
        mask=np.where(np.isfinite(data['ssfr']))
                   
        return data[mask[:][0]]

    def zcold_sfr(data):

        data['sfr']=np.log10(data['sfr'])
       
        mask=np.where(np.isfinite(data['sfr']))
        data=data[mask[:][0]]
        mask=np.where(np.isfinite(data['zcold']))                   
        return data[mask[:][0]]

    def ssfr_zcold(data):

        data['ssfr']=np.log10(data['ssfr'])
       
        mask=np.where(np.isfinite(data['ssfr']))
        data=data[mask[:][0]]
        data['zcold']=12+np.log10(0.0365507*data['Mzgas']/data['mcold'])        
        
        mask=np.where(np.isfinite(data['zcold']))

        print 'HERE ssfr2zcold:', min(data['zcold']), max(data['zcold']), min(data['ssfr']), max(data['ssfr'])
                   
        return data[mask[:][0]]

    def mcold_mstar(data):

        data['mcold']/=data['mstar']
        data['mcold']=np.log10(data['mcold'])        
        data['mstar']=np.log10(data['mstar'])
        mask=np.where(np.isfinite(data['mcold']))
            
        return data[mask[:][0]]
    

    def mbh_mbulge(data):

        data['mbh']=np.log10(data['mbh'])        
        data['mstar_spheroid']=np.log10(data['mstar_spheroid'])
        mask=np.where(np.isfinite(data['mstar_spheroid']))
            
        return data[mask[:][0]]

    def oh_mstar(data):

        data['Mzgas']=12+np.log10(0.0365507*data['Mzgas']/data['mcold'])
        data['mstar']=np.log10(data['mstar'])
        mask=np.where(np.isfinite(data['Mzgas']))
            
        return data[mask[:][0]]

    def zgas_mstar(data):

        data['Mzgas']=8.69+np.log10(data['Mzgas']/(data['mcold']*0.0134))
        data['mstar']=np.log10(data['mstar'])
        mask=np.where(np.isfinite(data['mstar']))
        data=data[mask[:][0]]
        mask=np.where(np.isfinite(data['Mzgas']))
        print min(data['Mzgas']), max(data['Mzgas'])
        return data[mask[:][0]]

    def zgas_mcold(data):

        data['Mzgas']=8.69+np.log10(data['Mzgas']/(data['mcold']*0.0134))                    
        data['mcold']=np.log10(data['mcold']/data['mstar'])
        mask=np.where(np.isfinite(data['Mzgas']))
        data=data[mask[:][0]]   
        mask=np.where(np.isfinite(data['mcold']))
         
        return data[mask[:][0]]

    def g_r_ssfr(data):

        try:
            data['ssfr']=data['sfr']/data['mstar']
        except:
            pass
        try:
            data['mAB_dA_total_g']-=data['mAB_dA_total_r']
        except:
            data['MAB_dA_total_g']-=data['MAB_dA_total_r']
            
        data['ssfr']=np.log10(data['ssfr'])             
            
        mask=np.where(np.isfinite(data['ssfr']))
            
        return data[mask[:][0]]

    def dmesa_ssfr(data):
        try:
            data['ssfr']=data['sfr']/data['mstar']
        except:
            pass
        data['ssfr']=np.log10(data['ssfr'])
        mask=np.where(np.isfinite(data['ssfr']))
            
        return data[mask[:][0]]


    def NFW_mhalo(data):

        try:
            data['mhalo_200c']=np.log10(data['mhalo_200c'])
        except:
            data['mbasic_200c']=np.log10(data['mbasic_200c'])
          
        return data
    
    def ssfr_mhalo(data):

        data['ssfr']=np.log10(data['ssfr']) 
        try:
            data['mhalo']=np.log10(data['mhalo_200c'])
        except:
            data['mhalo']=np.log10(data['mbasic_200c'])            
        
        mask=np.where(np.isfinite(data['mhalo']))
        #print min(data['mhalo']), max(data['mhalo']), min(data['ssfr']), max(data['ssfr'])
        
        return data[mask[:][0]]

    def ssfr_cgf(data):

        data['ssfr']=np.log10(data['ssfr']) 
        data['mcold']=np.log10(data['mcold']/data['mstar'])          
        
        mask=np.where(np.isfinite(data['mcold']))
        
        return data[mask[:][0]]

    def ssfr_i(data):

        data['ssfr']=np.log10(data['ssfr'])         
        mask=np.where(np.isfinite(data['ssfr']))
        
        return data[mask[:][0]]

    def mstar_mhalo(data):

        data['mstar']=np.log10(data['mstar']) 
        try:
            data['mhalo']=np.log10(data['mhalo_200c'])
        except:
            data['mhalo']=np.log10(data['mbasic_200c'])            
        
        mask=np.where(np.isfinite(data['mhalo']))
        #print min(data['mhalo']), max(data['mhalo']), min(data['ssfr']), max(data['ssfr'])
        
        return data[mask[:][0]]

    def zcold_mhalo(data):

        try:
            data['mhalo']=np.log10(data['mhalo_200c'])
        except:
            data['mhalo']=np.log10(data['mbasic_200c'])            
        
        mask=np.where(np.isfinite(data['mhalo']))
        #print min(data['mhalo']), max(data['mhalo']), min(data['ssfr']), max(data['ssfr'])
        
        return data[mask[:][0]]

    def NFW_Lr(data):
        try:
            mask=np.where(np.isfinite(data['L_SDSS_dA_total_r']))
        except:
            mask=np.where(np.isfinite(data['MAB_dA_total_r']))

        return data[mask[:][0]]

    def jbarcomp_mbar(data):
 
        try:
            data=data[np.where(data['mcold_disk']>=0)[:][0]]
            data=data[np.where(data['mstar_disk']>=0)[:][0]]
            data['mstar']=np.log10(data['mstar_disk']+data['mcold_disk'])
            data['angM_disk']=np.log10(data['angM_disk']/(data['mstar_disk']+data['mcold_disk']))
            mask=np.where(np.isfinite(data['angM_disk']))
            
        except:
            data=data[np.where(data['mcold_spheroid']>=0)[:][0]]
            data=data[np.where(data['mstar_spheroid']>=0)[:][0]]
            data['mstar']=np.log10(data['mstar_spheroid']+data['mcold_spheroid'])
            data['angM_spheroid']=np.log10(data['angM_spheroid']/(data['mstar_spheroid']+data['mcold_spheroid']))
            mask=np.where(np.isfinite(data['angM_spheroid']))
        data=data[mask[:][0]]            
        
        mask=np.where(np.isfinite(data['mstar']))

        return data[mask[:][0]]
    
    def jbar_mbar(data):
 

        data['angM_disk']=np.log10((data['angM_disk']+data['angM_spheroid'])/(data['mstar']+data['mcold']))
      
        data['mstar']=np.log10(data['mstar']+data['mcold'])
        mask=np.where(np.isfinite(data['mstar']))

        return data[mask[:][0]]

    def bdisk_bbulge(data):
        #data=data[np.where(data['orphan']>0)[:][0]]
        #data=data[np.where(data['mstar']<10**11.2)[:][0]]
        data['angM_disk']=np.log10(data['angM_disk']/(data['mstar_disk']+data['mcold_disk']))-2/3.0*np.log10(data['mstar_disk']+data['mcold_disk'])
        data['angM_spheroid']=np.log10(data['angM_spheroid']/(data['mstar_spheroid']+data['mcold_spheroid']))-2/3.0*np.log10(data['mstar_spheroid']+data['mcold_spheroid'])        
        mask=np.where(np.isfinite(data['angM_disk']))
        data=data[mask[:][0]]

        mask=np.where(np.isfinite(data['angM_spheroid']))                
        return data[mask[:][0]]

    def default(data):
 
        data['mstar']=np.log10(data['mstar'])
        mask=np.where(np.isfinite(data['mstar']))
        return data[mask[:][0]]
    
    def caseSwitcher(plot_type):
    
        choose = {
                'r-i_i': r_i_i,
                'r-i_mstar': r_i_mstar,
                'u-r_r': u_r_r,
                'u-r_mstar': u_r_mstar,
                'g-i_i': g_i_i,
                'g-i_mstar': g_i_mstar,                
                'dmesa_i': dmesa_i,
                'dmesa_mstar': dmesa_mstar,
                'r-i_g-r': r_i_g_r,
                'r': r,
                'sfr_mstar': sfr_mstar,
                'ssfr_mstar': ssfr_mstar,
                'ssfr_sfr': ssfr_sfr,
                'mcold_mstar': mcold_mstar,
                'mbh_mbulge': mbh_mbulge,
                'zgas_mstar': zgas_mstar,
                'zgas_mcold': zgas_mcold,
                'oh_mstar': oh_mstar,                 
                'dmesa_ssfr': dmesa_ssfr,
                'g-r_ssfr': g_r_ssfr,
                'NFW_mhalo': NFW_mhalo,
                'NFW_Lr': NFW_Lr,
                'mstar_mhalo': mstar_mhalo,
                'ssfr_mhalo': ssfr_mhalo,
                'zcold_mhalo': zcold_mhalo,
                'zcold_sfr': zcold_sfr,
                'ssfr_zcold': ssfr_zcold,
                'ssfr_cgf': ssfr_cgf,
                'ssfr_i': ssfr_i,
                'ssfr_r': ssfr_i,
                'jbarcomp_mbar': jbarcomp_mbar,
                'bdisk_bbulge': bdisk_bbulge,
                'jbar_mbar': jbar_mbar,
                'default': default                
                }
            
        func = choose.get(plot_type)
        
        if np.all(data['weight_tot']==-99.0) or np.all(data['weight_tot']==0.0):
            print 'set weights to 1!'
            data['weight_tot']=np.ones((data.size,), dtype=np.int8)
   
            print data['weight_tot']
        
        
        return func(data)

    data=caseSwitcher(key)
#    #data=data[(np.where(data['env_1024']==2))[0][:]]
    #data=data[(np.where(data['orphan']==0))[0][:]]
#    #data=data[(np.where(data['mhalo_200c']<10**13))[0][:]]
#
#    y=np.log10(data['ssfr'])
#    x=np.log10(data['sfr'])
#    prop=(y+11.16)/1.12-x
#    data = data[np.where(prop<=0.0)[:][0]] 
      
    
    print '--> after selection: ngal:', data.shape, '\n'
       
    return data 
    

plot_type = 'ssfr_mstar'

#mycomp='/home/doris/anaconda/pro/data/'


def myCatalog(load_from_file=False): 
    return load_SAM_catalog(mycomp+'anaconda/pro/data/Galacticus_1Gpc_run2/Galacticus_1Gpc_run2_z_0.0_tarsel.hdf5', 'Gal2')

def CMASS_Galacticus_color_catalog(load_from_file=False): 
    try:
        return load_SAM_catalog(mycomp+'anaconda/pro/data/Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS.hdf5', 'Gal-cols')
    except:
        return load_SAM_catalog(mycomp+'Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS.hdf5', 'Gal-cols')
        #return load_SAM_catalog(mycomp+'anaconda/pro/data/Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_test_nozboost_no5h_CMASS.hdf5', 'Galacticus color')


def LG_catalog(load_from_file=False):
    try:
         return load_SAM_catalog(mycomp+'anaconda/pro/data/LGALAXIES_500Mpc/LGALAXIES_500Mpc_z_0.56_tarsel_mags.hdf5', 'LG-all')
    except:
         return load_SAM_catalog('/store/erebos/doris/LGALAXIES_500Mpc_z_0.56_tarsel_mags.hdf5', 'LG-all')

def CMASS_LG_mass_catalog(load_from_file=False):
    try:
         return load_SAM_catalog(mycomp+'anaconda/pro/data/LGALAXIES_500Mpc/LGALAXIES_500Mpc_z_0.56_tarsel_mags_CMASS_mass_sample.hdf5', 'LG-mass')
    except:
         return load_SAM_catalog('/store/erebos/doris/LGALAXIES_500Mpc_z_0.56_tarsel_mags_CMASS_mass_sample.hdf5', 'LG-mass')

def CMASS_LG_down_catalog(load_from_file=False):
    try:
         return load_SAM_catalog(mycomp+'anaconda/pro/data/LGALAXIES_500Mpc/LGALAXIES_500Mpc_z_0.56_tarsel_mags_CMASS_down_sample.hdf5', 'LG-down')
    except:
         return load_SAM_catalog('/store/erebos/doris/LGALAXIES_500Mpc_z_0.56_tarsel_mags_CMASS_down_sample.hdf5', 'LG-down')

def CMASS_LG_density_vmax_catalog(load_from_file=False):
    try:
         return load_SAM_catalog(mycomp+'anaconda/pro/data/LGALAXIES_500Mpc/LGALAXIES_500Mpc_z_0.56_tarsel_mags_CMASS_density_vmax.hdf5', 'LG-dens $V_{max}$')
    except:
         return load_SAM_catalog('/store/erebos/doris/LGALAXIES_500Mpc_z_0.56_tarsel_mags_CMASS_density_vmax.hdf5', 'LG-dens $V_{max}$')

def CMASS_LG_density_mstar_catalog(load_from_file=False):
    try:
         return load_SAM_catalog(mycomp+'anaconda/pro/data/LGALAXIES_500Mpc/LGALAXIES_500Mpc_z_0.56_tarsel_mags_CMASS_density_mstar.hdf5', 'LG-dens $M_*$')
    except:
         return load_SAM_catalog('/store/erebos/doris/LGALAXIES_500Mpc_z_0.56_tarsel_mags_CMASS_density_mstar.hdf5', 'LG-dens $M_*$')



def Galacticus_catalog(load_from_file=False): return load_SAM_catalog('/store/erebos/doris/Galacticus_1Gpc_z_0.56_tarsel_new_mags.hdf5', 'Gal-all')

#def CMASS_Galacticus_density_catalog(load_from_file=False): return load_SAM_catalog(mycomp+'anaconda/pro/data/Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_density_sample.hdf5', 'Gal-dens')
def CMASS_Galacticus_density_catalog(load_from_file=False): return load_SAM_catalog(mycomp+'Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_density_sample.hdf5', '$n_{CMASS}$')
def CMASS_Galacticus_down_catalog(load_from_file=False): return load_SAM_catalog(mycomp+'Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample.hdf5', 'g-i>2.35')
def CMASS_Galacticus_down3_catalog(load_from_file=False): return load_SAM_catalog(mycomp+'Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3.hdf5', 'Gal-dens')
def CMASS_Galacticus_down3_crossed_catalog(load_from_file=False): return load_SAM_catalog(mycomp+'Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3_crossed.hdf5', 'crossed')


def CMASS_Galacticus_mass_catalog(load_from_file=False): return load_SAM_catalog(mycomp+'Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_mass_sample.hdf5', 'Gal-mass')
def CMASS_Galacticus_density_catalog_Mzstar(load_from_file=False): return load_SAM_catalog(mycomp+'Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_density_sample_highest_Mzstar.hdf5', 'density Mzstar')
def CMASS_Galacticus_density_catalog_Mzstar_spin(load_from_file=False): return load_SAM_catalog(mycomp+'/Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_density_sample_spin_Mzstar.hdf5', 'density Mzstar spin')
def CMASS_Galacticus_density_catalog_rdisk_spin(load_from_file=False): return load_SAM_catalog(mycomp+'Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_density_sample_spin_st_0.5_rdisk.hdf5', 'density rdisk spin')
def CMASS_Galacticus_density_catalog_rbulge_PS(load_from_file=False): return load_SAM_catalog(mycomp+'anaconda/pro/data/Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_density_sample_PS_rbulge.hdf5', 'density rbulge PS')
def CMASS_Galacticus_density_catalog_rhalf_PS(load_from_file=False): return load_SAM_catalog(mycomp+'anaconda/pro/data/Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_density_sample_PS_rhalf_mass.hdf5', 'density rhalf PS')

def CMASS_Galacticus_Mr_catalog(load_from_file=False): return load_SAM_catalog(mycomp+'anaconda/pro/data/Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_density_sample_Mr_gt_-21.5.hdf5', 'Gal-dens $M_r>-21.5$')

def CMASS_Galacticus_div_mass_color(load_from_file=False): return load_SAM_catalog(mycomp+'anaconda/pro/data/Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_mass_sample_200c_ngal_div_color.hdf5', 'Gal-dens')

def SAG_catalog(load_from_file=False):
    try:
        return load_SAM_catalog('/store/erebos/doris/SAG_1Gpc_z_0.56_tarsel_new_mags.hdf5', 'SAG-all')
    except:
        return load_SAM_catalog(mycomp+'anaconda/pro/data/SAG_1Gpc/SAG_1Gpc_z_0.56_tarsel_new_mags_CMASS.hdf5', 'SAG-cols')

def SAG_run2_catalog(load_from_file=False):
    try:
        return load_SAM_catalog('/store/erebos/doris/SAG_1Gpc_run2_z_0.56_tarsel_mags.hdf5', 'SAGv4-all')
    except:
        return load_SAM_catalog(mycomp+'anaconda/pro/data/SAG_1Gpc_run2/SAG_1Gpc_run2_z_0.56_tarsel_mags.hdf5', 'SAGv4-all')

def CMASS_SAG_color_catalog(load_from_file=False): return load_SAM_catalog(mycomp+'anaconda/pro/data/SAG_1Gpc/SAG_1Gpc_z_0.56_tarsel_new_mags_CMASS.hdf5', 'SAG-cols')
def CMASS_SAG_density_catalog(load_from_file=False): return load_SAM_catalog(mycomp+'anaconda/pro/data/SAG_1Gpc/SAG_1Gpc_z_0.56_tarsel_new_mags_CMASS_density_sample.hdf5', 'SAG-dens')
def CMASS_SAG_mass_catalog(load_from_file=False): return load_SAM_catalog(mycomp+'anaconda/pro/data/SAG_1Gpc/SAG_1Gpc_z_0.56_tarsel_new_mags_CMASS_mass_sample.hdf5', 'SAG-mass')

def CMASS_SAGv4_color_catalog(load_from_file=False): return load_SAM_catalog(mycomp+'anaconda/pro/data/SAG_1Gpc_run2/SAG_1Gpc_run2_z_0.56_tarsel_mags_CMASS.hdf5', 'SAGv4-cols')



def SAGE_catalog(load_from_file=False): 
    try:
        return load_SAM_catalog('/store/erebos/doris/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS.hdf5', 'Gal-col')
    except:
        return load_SAM_catalog(mycomp+'anaconda/pro/data/Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS.hdf5', 'Gal-col')

def CMASS_SAGE_density_catalog(load_from_file=False): return load_SAM_catalog('/store/erebos/doris/SAGE_1Gpc_z_0.56_CMASS_density_sample_mstar.hdf5', 'SAGE density')

def load_SAM_catalog(filename,
                     legend):
    
    data=load_catalog(filename,
                     'mstar',
                     'sfr',
                     name_add2='mstar_spheroid',
                     name_add3='mcold_spheroid',
                     name_add4='mstar')

    #r-i vs g-r tests
    #-----------------------------------------------------
    #test1 sample1
    #prop=-0.26+(y+0.7)/0.92857-x
    
    #test2 --> sample2 to extract green butterfly
    #prop=-0.3+(y+0.7)/0.92857-x      
    #z=data['mAB_dA_total_g']-data['mAB_dA_total_i']
    #data = data[np.where(z>2.35)[:][0]]
    
    #test3
    #prop=(y+1.05)/1.3-x
    #data = data[np.where(prop<0.0)[:][0]]    
    #------------------------------------------------------    
#    cuts = [0]
#    checknum=0
#    name='orphan'
#    for cut in cuts:
#        print 'select ', name, 'by:', cut, ' -->',
#        data_sorted = myData.selectData2Compute(data, 
#                                         selected_col=name, 
#                                         operator='>', 
#                                         condition=cut)
#        print 'num: ', data_sorted.shape
#        checknum+=data_sorted.size
#        print 'total num:', checknum
#
#
#
#    data=data_sorted
    print 'final shape:', data.shape
#
#    print data_sorted
    
    return adjust_data_for_plot(data, plot_type), legend, True



def CMASS_spall_wis_056_DR12v4(load_from_file=False): return CMASS_DR12_catalog('0.61_spall_wis', '0.51', info='DR12v4_compl_sample', legend='CMASS DR12')
def CMASS_spall_gra_056_DR12v4(load_from_file=False): return CMASS_DR12_catalog('0.61_spall_gra_nodust', '0.51', info='DR12v4_compl_sample', legend='CMASS DR12')
def CMASS_spall_por_mer_050_DR12v4(load_from_file=False): return CMASS_DR12_catalog('0.6_spall_por_merged', '0.5', info='DR12v4_compl_sample', legend='CMASS DR12')
def CMASS_spall_por_mer_056_DR12v4(load_from_file=False): return CMASS_DR12_catalog('0.61_spall_por_merged', '0.51', info='DR12v4_compl_sample', legend='CMASS DR12')
def CMASS_spall_por_SF_056_DR12v4(load_from_file=False): return CMASS_DR12_catalog('0.61_spall_por_SF', '0.51', info='DR12v4_compl_sample', legend='CMASS DR12')

def CMASS_spall_por_PS_all_DR12v4(load_from_file=False): return CMASS_DR12_catalog('0.7_spall_por_PS', '0.43', info='DR12v4_compl_sample', legend='0.43<z<0.70')
def CMASS_spall_por_PS_050_DR12v4(load_from_file=False): return CMASS_DR12_catalog('0.6_spall_por_PS', '0.5', info='DR12v4_compl_sample', legend='CMASS DR12')
def CMASS_spall_por_PS_056_DR12v4(load_from_file=False): return CMASS_DR12_catalog('0.61_spall_por_PS', '0.51', info='DR12v4_compl_sample', legend='0.51<z<0.61')
def CMASS_spall_por_PS_052_DR12v4(load_from_file=False): return CMASS_DR12_catalog('0.59_spall_por_PS', '0.52', info='DR12v4_compl_sample', legend='0.52<z<0.59') 
def CMASS_spall_por_PS_053_DR12v4(load_from_file=False): return CMASS_DR12_catalog('0.58_spall_por_PS', '0.53', info='DR12v4_compl_sample', legend='0.53<z<0.58') 
def CMASS_spall_por_PS_054_DR12v4(load_from_file=False): return CMASS_DR12_catalog('0.58_spall_por_PS', '0.54', info='DR12v4_compl_sample', legend='0.54<z<0.58')
def CMASS_spall_por_PS_055_DR12v4(load_from_file=False): return CMASS_DR12_catalog('0.56_spall_por_PS', '0.54', info='DR12v4_compl_sample', legend='0.54<z<0.56')
def CMASS_spall_por_PS_0554_DR12v4(load_from_file=False): return CMASS_DR12_catalog('0.56_spall_por_PS', '0.55', info='DR12v4_compl_sample', legend='0.554<z<0.56')       

def SMF_CMASS_spall_wis_056_DR12v4(load_from_file=False): return CMASS_DR12('0.61_spall_wis', '0.51', info='DR12v4_compl_sample', legend='WIS CMASS DR12')
def SMF_CMASS_spall_gra_056_DR12v4(load_from_file=False): return CMASS_DR12('0.61_spall_gra', '0.51', info='DR12v4_compl_sample', legend='GRA CMASS DR12')
def SMF_CMASS_spall_por_mer_all_DR12v4(load_from_file=False): return CMASS_DR12('0.7_spall_por_merged', '0.43', info='DR12v4_compl_sample_wg', legend='0.43<z<0.70')
def SMF_CMASS_spall_por_mer_056_DR12v4(load_from_file=False): return CMASS_DR12('0.61_spall_por_merged', '0.51', info='DR12v4_compl_sample', legend='CMASS DR12')
def SMF_CMASS_spall_por_mer_050_DR12v4(load_from_file=False): return CMASS_DR12('0.6_spall_por_merged', '0.5', info='DR12v4_compl_sample_wg_25bins_1e10', legend='CMASS DR12')
def SMF_CMASS_spall_por_mer_050_DR12v4_cross(load_from_file=False): return CMASS_DR12('0.6_spall_por_merged', '0.5', info='DR12v4_compl_sample_cross', legend='CMASS DR12', plot_key='SMF')
def mod_SMF_CMASS_spall_por_mer_050_DR12v4(load_from_file=False): return CMASS_DR12('0.6_spall_por_merged', '0.5', info='DR12v4_compl_sample_wg_fixed_35bins_1e10_-0.2', legend='modified CMASS DR12', plot_key='SMF')
def SFRF_CMASS_spall_por_mer_050_DR12v4_cross(load_from_file=False): return CMASS_DR12('0.6_spall_por_merged', '0.5', info='DR12v4_compl_sample_cross', legend='CMASS DR12', plot_key='SFRF')
def sSFRF_CMASS_spall_por_mer_050_DR12v4_cross(load_from_file=False): return CMASS_DR12('0.6_spall_por_merged', '0.5', info='DR12v4_compl_sample_cross', legend='CMASS DR12', plot_key='sSFRF')


def SMF_CMASS_spall_por_SF_056_DR12v4(load_from_file=False): return CMASS_DR12('0.61_spall_por_SF', '0.51', info='DR12v4_compl_sample', legend='POR SF CMASS DR12')

def SMF_CMASS_spall_por_PS_all_DR12v4(load_from_file=False): return CMASS_DR12('0.7_spall_por_PS', '0.43', info='DR12v4_compl_sample_wg', legend='0.43<z<0.70')
def SMF_CMASS_spall_por_PS_050_DR12v4(load_from_file=False): return CMASS_DR12('0.6_spall_por_PS', '0.5', info='DR12v4_compl_sample_wg', legend='Portsmouth PS')
def SMF_CMASS_spall_por_PS_056_DR12v4(load_from_file=False): return CMASS_DR12('0.61_spall_por_PS', '0.51', info='DR12v4_compl_sample', legend='CMASS DR12 0.51<z<0.61')
def SMF_CMASS_spall_por_PS_052_DR12v4(load_from_file=False): return CMASS_DR12('0.593_spall_por_PS', '0.523', info='DR12v4_compl_sample', legend='0.523<z<0.592')
def SMF_CMASS_spall_por_PS_053_DR12v4(load_from_file=False): return CMASS_DR12('0.58_spall_por_PS', '0.53', info='DR12v4_compl_sample', legend='CMASS DR12 0.53<z<0.58')
def SMF_CMASS_spall_por_PS_054_DR12v4(load_from_file=False): return CMASS_DR12('0.58_spall_por_PS', '0.54', info='DR12v4_compl_sample', legend='CMASS DR12 0.54<z<0.58')
def SMF_CMASS_spall_por_PS_055_DR12v4(load_from_file=False): return CMASS_DR12('0.56_spall_por_PS', '0.54', info='DR12v4_compl_sample', legend='CMASS DR12 0.54<z<0.56')
def SMF_CMASS_spall_por_PS_0554_DR12v4(load_from_file=False): return CMASS_DR12('0.56_spall_por_PS', '0.554', info='DR12v4_compl_sample', legend='CMASS DR12 0.554<z<0.56')

def SMF_CMASS_spall_por_SF_050_DR12v4(load_from_file=False): return CMASS_DR12('0.6_spall_por_PS', '0.5', info='DR12v4_compl_sample_wg', legend='Portsmouth SF')
def SMF_CMASS_spall_gra_050_DR12v4(load_from_file=False): return CMASS_DR12('0.6_spall_gra', '0.5', info='DR12v4_compl_sample_wg', legend='Granada')
def SMF_CMASS_spall_wis_050_DR12v4(load_from_file=False): return CMASS_DR12('0.6_spall_wis', '0.5', info='DR12v4_compl_sample_wg', legend='Wisconsin')



def LF_CMASS_spall_por_PS_050_DR12v4_Mg(load_from_file=False): return load_LF_CMASS_DR12('0.6_spall_por_PS', '0.5', 'g', 'M', info='DR12v4_compl_sample', legend='CMASS DR12')
def LF_CMASS_spall_por_PS_050_DR12v4_Mr(load_from_file=False): return load_LF_CMASS_DR12('0.6_spall_por_PS', '0.5', 'r', 'M', info='DR12v4_compl_sample', legend='CMASS DR12')
def LF_CMASS_spall_por_PS_050_DR12v4_Mi(load_from_file=False): return load_LF_CMASS_DR12('0.6_spall_por_PS', '0.5', 'i', 'M', info='DR12v4_compl_sample', legend='CMASS DR12')

def SMF_wis_055_DR12v4(load_from_file=False): return CMASS_DR12('wis', '0.55', info='DR12v4_compl_sample', legend='CMASS D12 Wis')
def SMF_gra_055_DR12v4(load_from_file=False): return CMASS_DR12('gra_nodust', '0.55', info='DR12v4_compl_sample', legend='Granada CMASS')
def SMF_por_PS_055_DR12v4(load_from_file=False): return CMASS_DR12('por_PS', '0.55', info='DR12v4_compl_sample', legend='Portsmouth PS CMASS')
def SMF_por_SF_055_DR12v4(load_from_file=False): return CMASS_DR12('por_SF', '0.55', info='DR12v4_compl_sample', legend='Portsmouth SF CMASS')
def SMF_por_merged_055_DR12v4(load_from_file=False): return CMASS_DR12('por_merged', '0.55', info='DR12v4_compl_sample', legend='Portsmouth CMASS')

def SMF_wis_043_DR12v4(load_from_file=False): return CMASS_DR12('0.75_wis', '0.43', info='DR12v4_compl_sample', legend='CMASS D12 Wis')
def SMF_gra_043_DR12v4(load_from_file=False): return CMASS_DR12('0.75_gra_nodust', '0.43', info='DR12v4_compl_sample', legend='Granada CMASS')
def SMF_por_PS_043_DR12v4(load_from_file=False): return CMASS_DR12('0.75_por_PS', '0.43', info='DR12v4_compl_sample', legend='Portsmouth PS CMASS')
def SMF_por_SF_043_DR12v4(load_from_file=False): return CMASS_DR12('0.75_por_SF', '0.43', info='DR12v4_compl_sample', legend='Portsmouth SF CMASS')
def SMF_por_merged_043_DR12v4(load_from_file=False): return CMASS_DR12('0.75_por_merged', '0.43', info='DR12v4_compl_sample', legend='Portsmouth CMASS')

def CMASS_DR12(catname,
               redshift,
               info='',
               legend='',
               plot_key='SMF'):
    
    try:        
        filename=mycomp+'anaconda/pro/myRun/histos/'+plot_key+'/'+plot_key+'_CMASS_z_'+redshift+'_tarsel_'+catname+'_'+info+'.txt'        
        #print 'load CMASS DR12', filename, '\n'
        data_array = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)         
    except:
        filename=mycomp+'anaconda/pro/myRun/histos/'+plot_key+'/'+plot_key+'_CMASS_SPALL_z_'+redshift+'_'+catname+'_'+info+'.txt'        
        print filename
        data_array = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)       
        
    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)
    
    #IMF and cosmology already corrected
    obs_data_array[:,0] = np.log10(data_array[:,0])#-0.2#-0.03925
    obs_data_array[:,1] = np.log10((data_array[:,1]))
    obs_data_array[:,4] = (data_array[:,4]/data_array[:,1])*0.434
    obs_data_array[:,5] = (data_array[:,5]/data_array[:,1])*0.434
    
    #print obs_data_array
    
    return obs_data_array, legend, False
    
def CMASS_DR12_catalog(catname,
                       redshift,
                       info='',
                       legend=''):


    filename=mycomp+'anaconda/pro/data/CMASS_SPALL/CMASS_SPALL_z_'+redshift+'_'+catname+'_'+info+'.hdf5'        

    obs_data_array=load_catalog(filename,
                                 'mAB_dA_total_i',
                                 'weight_tot',
                                 name_add='mAB_dA_total_g',
                                 name_add2='mAB_dA_total_r',
                                 name_add3='mstar')

    #obs_data_array = obs_data_array[(np.where(obs_data_array['mAB_dA_total_r']-obs_data_array['mAB_dA_total_i']<2.0))[:][0]]
    #obs_data_array = obs_data_array[(np.where(obs_data_array['mAB_dA_total_cut_dmesa']>0.55))[:][0]]
    #obs_data_array = obs_data_array[(np.where(obs_data_array['mAB_dA_total_i']<19.86 + 1.6*(obs_data_array['mAB_dA_total_cut_dmesa']-0.8)))[:][0]]          
    #obs_data_array = obs_data_array[(np.where(obs_data_array['mAB_dA_total_i']<19.90) and np.where(obs_data_array['mAB_dA_total_i']>17.50))[:][0]]

    #Kroupa-->Chabrier: MCha=0.9125*MKro or log10MCha=log10MKro - 0.03925
    obs_data_array['mstar']*=0.9125 #* (0.70/0.6777)**2
    #mask=np.where(np.isfinite(obs_data_array['mAB_dA_total_cut_dmesa']))
    #obs_data_array[mask[:][0]]
    #print obs_data_array['weight_tot']
    #print '--> after selection: ngal:', obs_data_array[mask[:][0]].shape, '\n'   

    return adjust_data_for_plot(obs_data_array, plot_type), legend, True   


#def wp_CMASS_DR12_N_all(load_from_file=False): return load_wp_CMASS_DR12('bigMD-dr12v4-N-NOfiberCol-errors.wp', legend='CMASS DR12')
#def wp_CMASS_DR12_N_cents(load_from_file=False): return load_wp_CMASS_DR12('bigMD-dr12v4-N-NOfiberCol-OnlyCent-errors.wp', legend='CMASS DR12')
#def wp_CMASS_DR12_N_all_1125(load_from_file=False): return load_wp_CMASS_DR12('cmass-dr12v4-Reid-11.25-errors.wp', legend='$\log_{10}(M_*$ $[M_{\odot}])>11.25$')
#def wp_CMASS_DR12_N_all_1135(load_from_file=False): return load_wp_CMASS_DR12('cmass-dr12v4-Reid-11.35-errors.wp', legend='$\log_{10}(M_*$ $[M_{\odot}])>11.35$')
#def wp_CMASS_DR12_N_all_1145(load_from_file=False): return load_wp_CMASS_DR12('cmass-dr12v4-Reid-11.45-errors.wp', legend='$\log_{10}(M_*$ $[M_{\odot}])>11.45$')
#def wp_CMASS_DR12_N_all_1155(load_from_file=False): return load_wp_CMASS_DR12('cmass-dr12v4-Reid-11.55-errors.wp', legend='$\log_{10}(M_*$ $[M_{\odot}])>11.55$')
#def wp_CMASS_DR12_N_all_1165(load_from_file=False): return load_wp_CMASS_DR12('cmass-dr12v4-Reid-11.65-errors.wp', legend='$\log_{10}(M_*$ $[M_{\odot}])>11.65$')
#
#def wp_CMASS_DR12_all_1125_1136(load_from_file=False): return load_wp_CMASS_DR12('cmass-dr12v4-Reid-bin1.wp', legend='$11.25<\log_{10}(M_*$ $[M_{\odot}])<11.36$')
#def wp_CMASS_DR12_all_1136_1148(load_from_file=False): return load_wp_CMASS_DR12('cmass-dr12v4-Reid-bin2.wp', legend='$11.36<\log_{10}(M_*$ $[M_{\odot}])<11.48$')
#def wp_CMASS_DR12_all_1148(load_from_file=False): return load_wp_CMASS_DR12('cmass-dr12v4-Reid-bin3.wp', legend='$\log_{10}(M_*$ $[M_{\odot}])>11.48$')

def wp_CMASS_DR12_N_all(load_from_file=False): return load_wp_CMASS_DR12('bigMD-dr12v4-N-NOfiberCol-errors.wp', legend='BigMD-LC')
def wp_CMASS_DR12_N_cents(load_from_file=False): return load_wp_CMASS_DR12('bigMD-dr12v4-N-NOfiberCol-OnlyCent-errors.wp', legend='BigMD-LC')
def wp_CMASS_DR12_N_all_1125(load_from_file=False): return load_wp_CMASS_DR12('cmass-dr12v4-Reid-bestmod-11.25-errors.wp', legend='cut1')
def wp_CMASS_DR12_N_all_1135(load_from_file=False): return load_wp_CMASS_DR12('cmass-dr12v4-Reid-bestmod-11.35-errors.wp', legend='cut2')
def wp_CMASS_DR12_N_all_1145(load_from_file=False): return load_wp_CMASS_DR12('cmass-dr12v4-Reid-bestmod-11.45-errors.wp', legend='cut3')
def wp_CMASS_DR12_N_all_1155(load_from_file=False): return load_wp_CMASS_DR12('cmass-dr12v4-Reid-bestmod-11.55-errors.wp', legend='cut4')
def wp_CMASS_DR12_N_all_1165(load_from_file=False): return load_wp_CMASS_DR12('cmass-dr12v4-Reid-bestmod-11.65-errors.wp', legend='cut5')

def wp_CMASS_DR12_all_1125_1136(load_from_file=False): return load_wp_CMASS_DR12('cmass-dr12v4-Reid-bestmod-bin1-errors.wp', legend='bin1')
def wp_CMASS_DR12_all_1136_1148(load_from_file=False): return load_wp_CMASS_DR12('cmass-dr12v4-Reid-bestmod-bin2-errors.wp', legend='bin2')
def wp_CMASS_DR12_all_1148(load_from_file=False): return load_wp_CMASS_DR12('cmass-dr12v4-Reid-bestmod-bin3-errors.wp', legend='bin3')


def load_wp_CMASS_DR12(filename,
                        legend):
    #BigMD-LC SAM binning         
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/'+filename , data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=0)             
        
    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)

    obs_data_array[:,0] = data_array[:,0]/0.6777
    obs_data_array[:,1] = np.log10(data_array[:,1]*data_array[:,0]/0.6777**2)
    obs_data_array[:,4] = data_array[:,2]/data_array[:,1]/0.6777*0.434
    obs_data_array[:,5] = obs_data_array[:,4]
    
    #print obs_data_array[:,[0,1]]

    return obs_data_array, legend, False

def load_LF_CMASS_DR12(catname,
                       redshift,
                       band,
                       mag,
                       info='',
                       legend=''):
    
    filename=mycomp+'anaconda/pro/myRun/histos/analyseTargetSelection/histo_CMASS_SPALL_z_'+redshift+'_'+catname+'_'+info+'_'+mag+'AB_dA_total_'+band+'.txt'        

    data_array = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)        
    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)

    obs_data_array[:,0]=data_array[:,0]
    obs_data_array[:,1]=np.log10(data_array[:,1])
    obs_data_array[:,4]=(data_array[:,4]/data_array[:,1])*0.43
    obs_data_array[:,5]=obs_data_array[:,5]

    return obs_data_array, legend, False


def CMASS_DR12_LF_all_056(load_from_file=False):
    
    #data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/analyseTargetSelection/histo_CMASS_SPALL_z_0.51_0.61_spall_por_merged_DR12v4_compl_sample_MAB_dA_total_i_RS.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/analyseTargetSelection/histo_CMASS_SPALL_z_0.5_0.6_spall_por_PS_DR12v4_compl_sample_MAB_dA_total_i.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)

    obs_data_array = np.zeros((data_array[:,0].size, 7), dtype=np.float32)
    
    obs_data_array[:,0]=data_array[:,0]
    obs_data_array[:,1]=np.log10(data_array[:,1])
    obs_data_array[:,4]=(data_array[:,4]/data_array[:,1])*0.43
    obs_data_array[:,5]=obs_data_array[:,5]
    
    return obs_data_array, 'CMASS DR12', False  


def CMASS_DR12_LF_RS_050(load_from_file=False):
    
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/analyseTargetSelection/histo_CMASS_SPALL_z_0.5_0.6_spall_por_merged_DR12v4_compl_sample_MAB_dA_total_i_Guo13.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)

    obs_data_array = np.zeros((data_array[:,0].size, 7), dtype=np.float32)
    
    obs_data_array[:,0]=data_array[:,0]
    obs_data_array[:,1]=np.log10(data_array[:,1])
    obs_data_array[:,4]=(data_array[:,4]/data_array[:,1])*0.43
    obs_data_array[:,5]=obs_data_array[:,5]
    
    return obs_data_array, 'CMASS DR12', False 
 

def CMASS_DR12_LF_RS_056(load_from_file=False):
    
    #data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/analyseTargetSelection/histo_CMASS_SPALL_z_0.51_0.61_spall_por_merged_DR12v4_compl_sample_g-i_gt_2.35_MAB_dA_total_i_RS.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/analyseTargetSelection/histo_CMASS_SPALL_z_0.5_0.6_spall_por_PS_DR12v4_compl_sample_g-i_gt_2.35_MAB_dA_total_i.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)

    obs_data_array = np.zeros((data_array[:,0].size, 7), dtype=np.float32)
    
    obs_data_array[:,0]=data_array[:,0]
    obs_data_array[:,1]=np.log10(data_array[:,1])
    obs_data_array[:,4]=(data_array[:,4]/data_array[:,1])*0.43
    obs_data_array[:,5]=obs_data_array[:,5]
    
    return obs_data_array, 'CMASS DR12', False   

def CMASS_cut(load_from_file=False):
    
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Maraston+13_dperp_vs_i_CMASS_cut.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)

    obs_data_array = np.zeros((data_array[:,0].size, 7), dtype=np.float32)
    
    obs_data_array[:,0]=data_array[:,0]
    obs_data_array[:,1]=data_array[:,1]

    return obs_data_array, 'CMASS selection', False 

def CMASS_dmesa_cut(load_from_file=False):
    
    obs_data_array = np.zeros((10, 7), dtype=np.float32)

    x_data=np.linspace(0.0,3,10)
    obs_data_array[:,0]=x_data  
    obs_data_array[:,1]= 0.55+(x_data/8.0)
    
    #print obs_data_array

    return obs_data_array, '$d_{\perp} > 0.55$', False

def CMASS_RS_doris_cut1(load_from_file=False):
    
    obs_data_array = np.zeros((10, 7), dtype=np.float32)
    
    x_data=np.linspace(-10,10,10)    
    obs_data_array[:,0]=x_data

    #Test2
    obs_data_array[:,1]=(x_data+0.26)*0.92857-0.7
    
    return obs_data_array, 'exp. RS cut1', False 

def CMASS_RS_doris_cut2(load_from_file=False):
    
    obs_data_array = np.zeros((10, 7), dtype=np.float32)
    
    x_data=np.linspace(-10,10,10)    
    obs_data_array[:,0]=x_data

    #Test1
    obs_data_array[:,1]=(x_data+0.3)*0.92857-0.7
    
    return obs_data_array, 'exp. RS cut2', False 

def CMASS_pop_sep(load_from_file=False):
    
    obs_data_array = np.zeros((10, 7), dtype=np.float32)
    
    x_data=np.linspace(-20,10,10)    
    obs_data_array[:,0]=x_data

    #Test1
    obs_data_array[:,1]=x_data*1.12-11.16
    
    return obs_data_array, 'pop sep', False 

def CMASS_RS_doris_cut3(load_from_file=False):
    
    obs_data_array = np.zeros((10, 7), dtype=np.float32)
    
    x_data=np.linspace(-10,10,10)    
    obs_data_array[:,0]=x_data

    #Test2
    obs_data_array[:,1]=x_data*1.3-1.05
    
    return obs_data_array, 'exp. RS cut3', False 

def CMASS_dmesa_sliding_cut(load_from_file=False):
    
    obs_data_array = np.zeros((10, 7), dtype=np.float32)

    x_data=np.linspace(0.0,3,10)
    obs_data_array[:,0]=x_data    

    M_i_min=19.48
    obs_data_array[:,1]= (M_i_min-19.86)/1.6 + 0.8 + x_data/8.0
    
    print obs_data_array

    return obs_data_array, 'sliding cut, $i$='+str(M_i_min), False

def create_snapzred_list():
    
    filename=mycomp+'anaconda/pro/data/sussing_tree/snapidzred.txt'        

    data_array = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)       
    length=50
        
    obs_data_array = np.zeros((length,5), dtype=np.float32)
         
    i=3
    a=0
    while a<data_array[:,0].size+2:
        try:
            obs_data_array[a,0] = data_array[i,0]
            obs_data_array[a,1] = data_array[i,1]
            obs_data_array[a,2] = data_array[i,2]
            obs_data_array[a,3] = data_array[i,3]
            if obs_data_array[a,1]>5.0:
                i+=4
            elif obs_data_array[a,1]>3.5:
                i+=3                
            else:
                i+=2
        except:
            i+=2                   
        a+=1


    i=0
    while i<obs_data_array[:,0].size-1:
        obs_data_array[i+1,4] = obs_data_array[i,3]-obs_data_array[i+1,3]
        i+=1

    filename_out = mycomp+'anaconda/pro/data/sussing_tree/snapidzred_out.txt'
    myOutput.writeIntoFile(filename_out, 
                           obs_data_array,
                           #myheader='# Elbaz et al. A&A 533, A119 (2011), pbwo5rks.CosmicCarnage2015 file: "/CalibrationData/PROPER_DATA/SpecificStarFormationRatevsStellarMassBlue/Elbaz+11/elbaz_11_masssfr.dat" \n(1) Mstar [Msun], (2) sSFR [Gyr-1], (3) binsize, (4) Phi error (-),  (5) Phi error (+)',
                           data_format='%i\t%0.6f\t%0.6f\t%0.6f\t%0.6f',
                           mydelimiter='\t')        

    import matplotlib.pyplot as plt
    
    fig = plt.figure(figsize=(12,10), dpi=150)
    myax = fig.add_subplot(111) 

    #plt.scatter3D(xs, ys, zs=0, zdir='z', s=20, c=None, depthshade=True, *args, **kwargs)                  
    myax.plot(obs_data_array[:,1],
            obs_data_array[:,5],
            alpha=0.1,
            ls='-',
            marker='o',
            color='r',
            markercol='')
    plt.show()
    plt.savefig(mycomp+'anaconda/pro/myRun/plots/plotXY/scatter_plot.png', format='png', rasterized=True, dpi=100, bbox_inches='tight', pad_inches=0.05) 
    #from matplotlib.backends.backend_pdf import PdfPages
    #pp = PdfPages(mycomp+'anaconda/pro/myRun/plots/plotXY/scatter_plot.pdf')
    #plt.savefig(pp, format='pdf', rasterized=True, dpi=20, pad_inches=0.05, transparent=True)      
    exit()



###################################################################################################################################################################
def DeLucia07_xi_CUT3_mstar(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mstar_xi_CUT3_DeLucia07.txt'   
    
    return Contreras_Fig6_twoPCF('DeLucia07', filename)

def Bower06_xi_CUT3_mstar(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mstar_xi_CUT3_Bower06.txt'   
    
    return Contreras_Fig6_twoPCF('Bower06', filename)

def DeLucia07_xi_CUT3_mcold(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mcold_xi_CUT3_DeLucia07.txt'   
    
    return Contreras_Fig6_twoPCF('DeLucia07', filename)

def Bower06_xi_CUT3_mcold(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mcold_xi_CUT3_Bower06.txt'   
    
    return Contreras_Fig6_twoPCF('Bower06', filename)
    
def DeLucia07_xi_CUT3_sfr(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_sfr_xi_CUT3_DeLucia07.txt'   
    
    return Contreras_Fig6_twoPCF('DeLucia07', filename)

def Bower06_xi_CUT3_sfr(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_sfr_xi_CUT3_Bower06.txt'   
    
    return Contreras_Fig6_twoPCF('Bower06', filename)
    


def DeLucia07_xi_CUT2_mstar(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mstar_xi_CUT2_DeLucia07.txt'   
    
    return Contreras_Fig6_twoPCF('DeLucia07', filename)

def Bower06_xi_CUT2_mstar(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mstar_xi_CUT2_Bower06.txt'   
    
    return Contreras_Fig6_twoPCF('Bower06', filename)

def DeLucia07_xi_CUT2_mcold(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mcold_xi_CUT2_DeLucia07.txt'   
    
    return Contreras_Fig6_twoPCF('DeLucia07', filename)

def Bower06_xi_CUT2_mcold(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mcold_xi_CUT2_Bower06.txt'   
    
    return Contreras_Fig6_twoPCF('Bower06', filename)
    
def DeLucia07_xi_CUT2_sfr(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_sfr_xi_CUT2_DeLucia07.txt'   
    
    return Contreras_Fig6_twoPCF('DeLucia07', filename)

def Bower06_xi_CUT2_sfr(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_sfr_xi_CUT2_Bower06.txt'   
    
    return Contreras_Fig6_twoPCF('Bower06', filename) 



def DeLucia07_xir2_CUT3_mstar(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mstar_xir2_CUT3_DeLucia07.txt'   
    
    return Contreras_Fig6_twoPCF('DeLucia07', filename)

def Bower06_xir2_CUT3_mstar(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mstar_xir2_CUT3_Bower06.txt'   
    
    return Contreras_Fig6_twoPCF('Bower06', filename)

def DeLucia07_xir2_CUT3_mcold(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mcold_xir2_CUT3_DeLucia07.txt'   
    
    return Contreras_Fig6_twoPCF('DeLucia07', filename)

def Bower06_xir2_CUT3_mcold(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mcold_xir2_CUT3_Bower06.txt'   
    
    return Contreras_Fig6_twoPCF('Bower06', filename)
    
def DeLucia07_xir2_CUT3_sfr(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_sfr_xir2_CUT3_DeLucia07.txt'   
    
    return Contreras_Fig6_twoPCF('DeLucia07', filename)

def Bower06_xir2_CUT3_sfr(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_sfr_xir2_CUT3_Bower06.txt'   
    
    return Contreras_Fig6_twoPCF('Bower06', filename)
    

def DeLucia07_xir2_CUT2_mstar(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mstar_xir2_CUT2_DeLucia07.txt'   
    
    return Contreras_Fig6_twoPCF('DeLucia07', filename)

def Bower06_xir2_CUT2_mstar(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mstar_xir2_CUT2_Bower06.txt'   
    
    return Contreras_Fig6_twoPCF('Bower06', filename)

def DeLucia07_xir2_CUT2_mcold(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mcold_xir2_CUT2_DeLucia07.txt'   
    
    return Contreras_Fig6_twoPCF('DeLucia07', filename)

def Bower06_xir2_CUT2_mcold(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mcold_xir2_CUT2_Bower06.txt'   
    
    return Contreras_Fig6_twoPCF('Bower06', filename)
    
def DeLucia07_xir2_CUT2_sfr(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_sfr_xir2_CUT2_DeLucia07.txt'   
    
    return Contreras_Fig6_twoPCF('DeLucia07', filename)

def Bower06_xir2_CUT2_sfr(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_sfr_xir2_CUT2_Bower06.txt'   
    
    return Contreras_Fig6_twoPCF('Bower06', filename)    

def DeLucia07_xir2_CUT1_mstar(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mstar_xir2_CUT1_DeLucia07.txt'   
    
    return Contreras_Fig6_twoPCF('DeLucia07', filename)

def Bower06_xir2_CUT1_mstar(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mstar_xir2_CUT1_Bower06.txt'   
    
    return Contreras_Fig6_twoPCF('Bower06', filename)

def DeLucia07_xir2_CUT1_mcold(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mcold_xir2_CUT1_DeLucia07.txt'   
    
    return Contreras_Fig6_twoPCF('DeLucia07', filename)

def Bower06_xir2_CUT1_mcold(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_mcold_xir2_CUT1_Bower06.txt'   
    
    return Contreras_Fig6_twoPCF('Bower06', filename)
    
def DeLucia07_xir2_CUT1_sfr(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_sfr_xir2_CUT1_DeLucia07.txt'   
    
    return Contreras_Fig6_twoPCF('DeLucia07', filename)

def Bower06_xir2_CUT1_sfr(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Contreras+13_sfr_xir2_CUT1_Bower06.txt'   
    
    return Contreras_Fig6_twoPCF('Bower06', filename)   

def Sridhar_xi(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/Sridhar+17_mstar_xi.txt'   
    
    return Contreras_Fig6_twoPCF('Sridhar+17', filename) 

def vanDaalen_wrp(load_from_file=False):

    filename=mycomp+'anaconda/pro/OBS/vanDaalen+15_pro2PCF_SDSS_Fig2_10.77-11.27.txt'   
    
    return Contreras_Fig6_twoPCF('vanDaalen+15 SDSS DR7', filename) 

def Contreras_Fig6_twoPCF(reference, filename):
    
    data = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)
    obs_data_array = np.zeros((data[:,0].size,7), dtype=np.float32)

    
    obs_data_array[:,0] = data[:,0]/0.6777
    
    if filename.find('xir2')!=-1:
        if filename.find('mstar')!=-1:
            obs_data_array[:,1] = data[:,1]/(0.6777**2)
        else:          
           obs_data_array[:,1] = (data[:,1]/(data[:,0]**1.3))*obs_data_array[:,0]**2

    elif filename.find('pro')!=-1:
        obs_data_array[:,0] = data[:,0]
        obs_data_array[:,1] = data[:,1]
    else:    
        obs_data_array[:,1] = data[:,1]

    obs_data_array[:,[0,1]]=np.log10(obs_data_array[:,[0,1]])
    obs_data_array[:,4] = (data[:,4]/data[:,1])*0.434
    obs_data_array[:,5] = (data[:,5]/data[:,1])*0.434

    return obs_data_array, reference, False  
###################################################################################################################################################################
    

###################################################################################################################################################################
def Bower06_mstar(load_from_file=False):

    return Contreras_Fig3_cumhistos('Bower06', 'mstar')

def DeLucia07_mstar(load_from_file=False):

    return Contreras_Fig3_cumhistos('DeLucia07', 'mstar')

def Bower06_mcold(load_from_file=False):

    return Contreras_Fig3_cumhistos('Bower06', 'mcold')

def DeLucia07_mcold(load_from_file=False):

    return Contreras_Fig3_cumhistos('DeLucia07', 'mcold')

def Bower06_sfr(load_from_file=False):

    return Contreras_Fig3_cumhistos('Bower06', 'sfr')

def DeLucia07_sfr(load_from_file=False):

    return Contreras_Fig3_cumhistos('DeLucia07', 'sfr')

def Contreras_Fig3_cumhistos(reference,
                             name):

    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Contreras+13_CSMF_'+name+'_'+reference+'_Fig3.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)
    obs_data_array = np.zeros((data[:,0].size,7), dtype=np.float)
    
    if name=='sfr':
        obs_data_array[:,0] = np.log10(data[:,0])
    else:    
        obs_data_array[:,0] = np.log10(data[:,0]/0.6777)
        
    obs_data_array[:,1] = np.log10(data[:,1]*0.6777**3)
    obs_data_array[:,4] = (data[:,4]/data[:,1])*0.434
    obs_data_array[:,5] = (data[:,5]/data[:,1])*0.434
 
    return obs_data_array, reference, False  
###################################################################################################################################################################

def Elbaz11_function(load_from_file=False):

    reference={}
    reference  = {'name': 'Elbaz+11', 
                  'filename': 'elbaz_11_masssfr.dat', 
                  'legend': 'Elbaz+11', 
                  'paper': 'Elbaz et al. A&A 533, A119 (2011)', 
                  'little_h': 0.6777,
                  'IMF_corr': -0.24}
    print mycomp+'anaconda/pro/OBS/CalibrationData/DATA/'+reference['filename']              
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/CalibrationData/DATA/'+reference['filename'], data_format='ASCII', nr_col=2, data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=2)                   

    data[:,0]=10**(np.log10(data[:,0]/reference['little_h'])+reference['IMF_corr'])   
    data[:,1] = data[:,1]/data[:,0]

    return ssfr2mstar_ref(reference,data)
    
def Elbaz11_dots(load_from_file=False):

    reference={}
    reference  = {'name': 'Elbaz+11', 
                  'filename': 'elbaz_11_masssfr.dat', 
                  'legend': 'Elbaz+11 compilation', 
                  'paper': 'Elbaz et al. A&A 533, A119 (2011)', 
                  'little_h': 0.6777,
                  'IMF_corr': -0.24}
                   
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/CalibrationData/DATA/'+reference['filename'], data_format='ASCII', nr_col=2, data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=2)                   

    data[:,0]=10**(np.log10(data[:,0]/reference['little_h'])+reference['IMF_corr'])   
    data[:,1]=data[:,1]/data[:,0]
    print data
    
    return data, reference['legend'], True

def Elbaz112(load_from_file=False):
    
    if load_from_file==True:

        obs_data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/GOOD-Herschel_Elbaz+11_sSFR2Mstar_z_0.1_2.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)


    else:
    
        obs_data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/CalibrationData/PROPER_DATA/SpecificStarFormationRatevsStellarMassBlue/Elbaz+11/elbaz_11_masssfr.dat', data_format='ASCII', nr_col=2, nr_rows=28610, data_shape='shaped', delim=' ', mydtype=np.float, skiprow=1)
        #obs_data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/elbaz+11_ex1.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)
    
        myFuncs = mF.MyFunctions()
    
        obs_data_array[:,0] = 10**obs_data_array[:,0]
        
        obs_data_array[:,1] = 10**obs_data_array[:,1]*1e9
        
        data = myData.selectData2Compute(obs_data_array, 
                                      selected_col=1, 
                                      operator='>', 
                                      condition=0.01)


 
        obs_data_array = myFuncs.binUp(data, 
                                        20,
                                        {'0': 'mstar', '1': 'sfr'},
                                        histo_min=min(obs_data_array[:,0]),
                                        histo_max=max(obs_data_array[:,0]),
                                        log10=True)

    
        obs_data_array[:,0] = obs_data_array[:,0]/0.6777
    
                           
        filename_out = mycomp+'anaconda/pro/OBS/GOOD-Herschel_Elbaz+11_sSFR2Mstar_z_0.1_2.txt'
        myOutput.writeIntoFile(filename_out, 
                               obs_data_array,
                               myheader='# Elbaz et al. A&A 533, A119 (2011), pbwo5rks.CosmicCarnage2015 file: "/CalibrationData/PROPER_DATA/SpecificStarFormationRatevsStellarMassBlue/Elbaz+11/elbaz_11_masssfr.dat" \n(1) Mstar [Msun], (2) sSFR [Gyr-1], (3) binsize, (4) Phi error (-),  (5) Phi error (+)',
                               data_format='%0.1f\t%0.5f\t%0.1f\t%0.5f\t%0.5f',
                               mydelimiter='\t')                              

    return obs_data_array, 'GOOD-Herschel: Elbaz+11 (z~0.1)', False


def Guo18_SMF_CMASS(load_from_file=True):
        
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Guo+18_SMF_BOSS_z_0.5_0.6.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=2)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)

    obs_data_array[:,0] = data_array[:,0]
    obs_data_array[:,1] = np.log10(data_array[:,1]*1e-5)
    obs_data_array[:,4] = obs_data_array[:,1] - (data_array[:,2]/data_array[:,1])*0.434
    obs_data_array[:,5] = obs_data_array[:,1] + (data_array[:,2]/data_array[:,1])*0.434
          
    print obs_data_array        
         
    return obs_data_array, 'Guo+18', False


def Henriques15_SMF_1_0(load_from_file=True):
        
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Henriques+15_SMF_z_1.0_Fig2.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=6)

    obs_data_array = np.zeros((data_array[:,0].size,5), dtype=np.float32)

    obs_data_array[:,0] = data_array[:,0] - 2*np.log10(0.6777)
    obs_data_array[:,1] = np.log10((10**data_array[:,1])*(0.6777**3))
    obs_data_array[:,3] = data_array[:,4]
    obs_data_array[:,4] = data_array[:,5]
          
    print obs_data_array        
         
    return obs_data_array, 'Compilation Henriques+15 z=1.0', False

def LieAndWhite09(load_from_file=True):
   
    masses = np.arange(8.0, 9.33, 0.33)   
    data_schechterMFit = mL.schechterFit(masses, 0.01465, -1.1309, 9.6124)
    
    mass = masses
    fit= data_schechterMFit
    
    masses = np.arange(9.67, 10.67, 0.33)   
    data_schechterMFit = mL.schechterFit(masses, 0.0132, -0.9004, 10.3702)

    mass = np.concatenate((mass, masses))
    fit = np.concatenate((fit, data_schechterMFit))

    masses = np.arange(11.0, 12.0, 0.33)   
    data_schechterMFit = mL.schechterFit(masses, 0.00446, -1.9918, 10.7104)

    mass = np.concatenate((mass, masses))
    fit = np.concatenate((fit, data_schechterMFit))
    
    obs_data_array = np.zeros((13,5), dtype=np.float)
    
    mass=mass.astype(np.float)
    fit=fit.astype(np.float)

    obs_data_array[:,0] = mass
    obs_data_array[:,1] = np.log10(fit)
  
    return obs_data_array, 'Lie&White09', False

def Maraston13_BOSS_CMASS_model(load_from_file=False):

    legend_obs='Maraston+13 BOSS-CMASS'

    if load_from_file ==True:

        obs_data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/'+str(legend_obs)+'.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=1)
    else:

        data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Maraston+13_SMF_BOSS_CMASS_model_blue_squars_z_0.5_0.6.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)

        obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float) 

        #without little-h
        obs_data_array[:,0] = 10**(data_array[:,0])
        obs_data_array[:,1] = data_array[:,1]*(0.719**3)

        #with little-h
        obs_data_array[:,0] = data_array[:,0]+np.log10(0.73/0.6777)+np.log10(0.6777)-0.03925 #Corrected from Kroupa to Charbier Lacey+16 from Violeta
        obs_data_array[:,1] = np.log10(data_array[:,1])-np.log10(0.7/0.6777)*3
        
        #print 'mydata:', obs_data_array[:,[0,1]]
   
    return obs_data_array, legend_obs, False

def Maraston13_BOSS_model(load_from_file=False):
    #model

    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Maraston+13_SMF_model_blue_dashed_line_z_0.5_0.6.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float) 

    obs_data_array[:,0] = data_array[:,0]+np.log10(0.73/0.6777)+np.log10(0.6777)-0.03925 #Corrected from Kroupa to Charbier Lacey+16 from Violeta
    obs_data_array[:,1] = np.log10(data_array[:,1])-np.log10(0.73/0.6777)*3
        
    
    return obs_data_array, 'Maraston+13', False

def Maraston13_M09_BOSS_CMASS_model(load_from_file=False):

    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Maraston+13_SMF_BOSS_CMASS_M09_red_dots_z_0.5_0.6', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=1)
  
    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float) 


    obs_data_array[:,0] = data_array[:,0]+np.log10(0.73/0.6777)+np.log10(0.6777)-0.03925
    obs_data_array[:,1] = np.log10(data_array[:,1])-np.log10(0.73/0.6777)*3
    obs_data_array[:,4] = (data_array[:,4]/data_array[:,1])*0.434
    obs_data_array[:,5] = (data_array[:,5]/data_array[:,1])*0.434
          
    return obs_data_array, 'Maraston+09 CMASS', False

def Marastaon13_dperp_vs_i_contour0(load_from_file=False):
    
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Maraston+13_dperp_vs_i_contour_level0.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)

    obs_data_array = np.zeros((data_array[:,0].size, 7), dtype=np.float32)
    
    obs_data_array[:,0]=data_array[:,0]
    obs_data_array[:,1]=data_array[:,1]
    
    return obs_data_array, 'Maraston+13 0.54<z<0.56', False

def Marastaon13_dperp_vs_i_contour1(load_from_file=False):
    
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Maraston+13_dperp_vs_i_contour_level1.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)

    obs_data_array = np.zeros((data_array[:,0].size, 7), dtype=np.float32)
    
    obs_data_array[:,0]=data_array[:,0]
    obs_data_array[:,1]=data_array[:,1]
    
    return obs_data_array, 'Maraston+13 0.54<z<0.56', False

def Marastaon13_dperp_vs_i_contour2(load_from_file=False):
    
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Maraston+13_dperp_vs_i_contour_level2.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)

    obs_data_array = np.zeros((data_array[:,0].size, 7), dtype=np.float32)
    
    obs_data_array[:,0]=data_array[:,0]
    obs_data_array[:,1]=data_array[:,1]
    
    return obs_data_array, 'Maraston+13 0.54<z<0.56', False

def Marastaon13_dperp_vs_i_contour3(load_from_file=False):
    
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Maraston+13_dperp_vs_i_contour_level3.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)

    obs_data_array = np.zeros((data_array[:,0].size, 7), dtype=np.float32)
    
    obs_data_array[:,0]=data_array[:,0]
    obs_data_array[:,1]=data_array[:,1]
    
    return obs_data_array, 'Maraston+13 0.54<z<0.56', False

def Marastaon13_dperp_vs_i_contour4(load_from_file=False):
    
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Maraston+13_dperp_vs_i_contour_level4.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)

    obs_data_array = np.zeros((data_array[:,0].size, 7), dtype=np.float32)
    
    obs_data_array[:,0]=data_array[:,0]
    obs_data_array[:,1]=data_array[:,1]
    
    return obs_data_array, 'Maraston+13 0.54<z<0.56', False

def Marastaon13_rmini_vs_mstar_model_contour(load_from_file=False):
    
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Maraston+13_model_CMASS_0.50_z_0.60_r-r_vs_i_contour_levels.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)

    obs_data_array = np.zeros((data_array[:,0].size, 7), dtype=np.float32)
    
    obs_data_array[:,0]=data_array[:,0]
    obs_data_array[:,1]=data_array[:,1]
    
    return obs_data_array, 'Maraston+13 Model 0.50<z<0.60', False

def Marastaon13_rmini_vs_mstar_BOSS_CMASS_contour(load_from_file=False):
    
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Maraston+13_BOSS_CMASS_0.50_z_0.60_r-r_vs_i_contour_levels.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)

    obs_data_array = np.zeros((data_array[:,0].size, 7), dtype=np.float32)
    
    obs_data_array[:,0]=data_array[:,0]
    obs_data_array[:,1]=data_array[:,1]
    
    return obs_data_array, 'BOSS CMASS 0.50<z<0.60 (Maraston+13)', False

def Montero09_blue_red_sep(load_from_file=False):
   
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Montero_SDSS_DR6_blue-red_separation.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)
    obs_data_array = np.zeros((data[:,0].size,7), dtype=np.float)     
    obs_data_array[:,0] = data[:,0]-5.0*np.log10(0.6777)
    obs_data_array[:,1] = data[:,1]

    return obs_data_array, 'Strateva+01', False

def Montero16_dperp_cut(load_from_file=False):
   
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Montero_CMASS_dperp_cut.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)
    obs_data_array = np.zeros((data[:,0].size,7), dtype=np.float)     
    obs_data_array[:,0] = data[:,0]
    obs_data_array[:,1] = data[:,1]

    return obs_data_array, 'Montero-Dorta+16', False

def Montero16_LF_RS_CMASS_DR10(load_from_file=False):
   
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Montero-Dorta+16_Fig4_LF_RS_CMASS_DR19_z_0.55.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)
    obs_data_array = np.zeros((data[:,0].size,7), dtype=np.float)     
    obs_data_array[:,0] = data[:,0]+5*np.log10(0.7)-5*np.log10(0.6777)
    obs_data_array[:,1] = data[:,1]

    return obs_data_array, 'CMASS DR10', False

def Montero16_sliding_dperp_cut(load_from_file=False):
   
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Montero_CMASS_sliding_dperp_cut.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)
    obs_data_array = np.zeros((data[:,0].size,7), dtype=np.float)     
    obs_data_array[:,0] = data[:,0]
    obs_data_array[:,1] = data[:,1]

    return obs_data_array, 'Montero-Dorta+16', False



def Montero16_LF_i(load_from_file=False): return Montero16_LF_intrinic('i', 0.7857e-3*(0.6777/0.7)**3, 0.0247e-3*(0.6777/0.7)**3, -21.7151+5*np.log10(0.7)-5*np.log10(0.6777), 0.0123-5*np.log10(0.7)+5*np.log10(0.6777), -1.0, 0.555)

def Montero16_LF_intrinic(band, 
                          phiStar, 
                          phiStar_error, 
                          LStar,
                          LStar_error,
                          alpha, 
                          redshift):
   
    obs_data_array = np.zeros((100,7), dtype=np.float)  

    obs_data_array[:,0] = np.linspace(-25,-20,100)

    from cosmolopy import luminosityfunction as LF
    
    obs_data_array[:,1] = np.log10(LF.schechterM(obs_data_array[:,0], phiStar, alpha, LStar)) 
    obs_data_array[:,4] = obs_data_array[:,1] - LF.schechterM(obs_data_array[:,0], phiStar_error, alpha, LStar+LStar_error)/10**obs_data_array[:,1]*0.434
    obs_data_array[:,5] = obs_data_array[:,1] + LF.schechterM(obs_data_array[:,0], phiStar_error, alpha, LStar+LStar_error)/10**obs_data_array[:,1]*0.434
    
    #print obs_data_array[:, [0,1,4,5]]

    return obs_data_array, 'Montero-Dorta+16', False

def Montero16_RS(load_from_file=False):

    obs_data_array = np.zeros((5,7), dtype=np.float)     

    obs_data_array[:,0]= [1.6831398562741853, 1.6716417910447763,1.667910447761194,1.6641791044776117,1.6603920695032306]
    obs_data_array[:,1]= [1.0407407407407407, 1.0261194029850749,1.0106965174129354,0.9995946194951169,0.9925373134328361]    
    obs_data_array[:,3]= [18.9875, 19.2125,19.4375,19.6625,19.8875]  
    
    return obs_data_array, 'Montero-Dorta+16', False

def Moustakas13_PRIMUS_z025(load_from_file=False):

    legend_obs='PRIMUS (z~0.25)'

    if load_from_file ==True:


        obs_data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/'+str(legend_obs)+'.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=1)


    else:

        obs_data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Moustakas+13_PRIMUS_Table4.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)
        
        print 10**obs_data_array[:,0]
        print obs_data_array[:,1]
        print obs_data_array[:,2]

        obs_data_array[:,0] = 10**obs_data_array[:,0]
        obs_data_array[:,1] = 10**obs_data_array[:,1]
        obs_data_array[:,3] = obs_data_array[:,1]/(obs_data_array[:,2]**0.5)
        obs_data_array[:,4] = obs_data_array[:,3]
       
        print 'mydata:', obs_data_array      
       
        filename_out = mycomp+'anaconda/pro/OBS/'+str(legend_obs)+'.txt'
        myOutput.writeIntoFile(filename_out, 
                               obs_data_array,
                               myheader='Moustakas et al. 2013 (AJ 767:50 (34pp)) PRIMUS Table 4: 0.2 < z 0.3 /n #(1) Mstar [Msun], (2) Phi [Msun Mpc^-3], (3) Phi error (-),   (4) Phi error (+)',
                               data_format="%0.5f",
                               mydelimiter='\t')

    return obs_data_array, legend_obs, False
    
def Moustakas13_SDSS_GALEX_z010(load_from_file=False):

    if load_from_file ==True:

        data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Moustakas+13_SMF_all_Fig4_corrected.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=1)

    else:

        #obs_data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/CalibrationData/DATA/MF_Moustakas13_z0001_all.dat', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=1)
        data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Moustakas+13_SMF_all_Fig4.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)
        
        obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float) 

        #log-scaled axis
        obs_data_array[:,0] = 10**data_array[:,0]
        obs_data_array[:,1] = 10**data_array[:,1]
        #obs_data_array[:,3] = obs_data_array[:,1]/(obs_data_array[:,5]**0.5)

        #log values ploted
        obs_data_array[:,0] = data_array[:,0]
        obs_data_array[:,1] = data_array[:,1]

        obs_data_array[:,4] = data_array[:,1] + data_array[:,4]
        obs_data_array[:,5] = data_array[:,1] - data_array[:,5]        

        obs_data_array[:,4] = data_array[:,5]*0.43
        obs_data_array[:,5] = data_array[:,5]*0.43 
        
        #print 'mydata:', obs_data_array
       
        filename_out = mycomp+'anaconda/pro/OBS/Moustakas+13_SMF_all_Fig4_corrected.txt'
        myOutput.writeIntoFile(filename_out, 
                               obs_data_array,
                               myheader='Moustakas et al. 2013 (AJ 767:50 (34pp)) SDSS-GALEX Table 3: z~0.1 /n #(1) Mstar [Msun], (2) Phi [Msun Mpc^-3], (3) X-error, (4) Phi error (-), (5) Phi error (+), (6) N-count',
                               data_format="%0.5f",
                               mydelimiter='\t')
    #'SDSS-GALEX (z~0.1)'
    return obs_data_array, 'SDSS-GALEX', False

def Moustakas13_PRIMUS_z050_065(load_from_file=False):

    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Moustakas+13_PRIMUS_z_0.50_0.65_Table4_all.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)
    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float)
    
    obs_data_array[:,0] = data_array[:,0]#-2*np.log10(0.6777)
    obs_data_array[:,1] = data_array[:,1]+3*np.log10(0.6777)
    obs_data_array[:,4] = (10**obs_data_array[:,1]/(data_array[:,2]**0.5))/10**obs_data_array[:,1]*0.434
    obs_data_array[:,5] = obs_data_array[:,4]
   
    return obs_data_array, 'PRIMUS 0.50<z<0.65', False

def NYU_VAGC_SDSS_LOW_Z(load_from_file=False):

   
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/NYU_VAGC_LOWZ_DR4_Blanton+05.txt', data_format='ASCII', data_shape='shaped', delim='  ', mydtype=np.float, skiprow=5)
    obs_data_array = np.zeros((data[:,0].size,7), dtype=np.float)     
    #u=2
    #g=3
    #r=4
    #i=5
    #z=6
        
    #i-band=5
    obs_data_array[:,0] = data[:,5]
    #r-i
    obs_data_array[:,1] = data[:,4]-data[:,5]

    #r-band
    obs_data_array[:,0] = data[:,4]
    #g-r
    obs_data_array[:,1] = data[:,3]-data[:,4]
    
    #u-r
    obs_data_array[:,1] = data[:,2]-data[:,4]


    return obs_data_array, 'NYU-VAGC Blanton+05', True

def Ports_BOSS_DR12(load_from_file=False):

    #SDSS BOSS DR12 - Portsmouth Catalog passive Salpeter
    #orignial name: portsmouth_stellarmass_starforming_salp-DR12.fits
    #downloaded: 11/28/2016 from https://data.sdss.org/datamodel/files/BOSS_GALAXY_REDUX/GALAXY_VERSION/portsmouth_stellarmass.html
    #comment: magnitude in the SDSS bands in observed frame see webpage for definition of extinction
    # (1) RA	(2) DEC	(3) redshift   (4) u	  (5) g   (6) r	  (7) i	(8) z	   (9) log mstar [Msun] (10) sfr [Msunyr-1]
   
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/portsmouth_stellarmass_starforming_salp-DR12.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=5)
    

    data = myData.selectData2Compute(data, 
                                  selected_col=2, 
                                  operator='>', 
                                  condition=0.50) 

    data = myData.selectData2Compute(data, 
                                  selected_col=2, 
                                  operator='<', 
                                  condition=0.60)
                                  
    print 'data after redshift selection!', data.shape
    obs_data_array = np.zeros((data[:,0].size,7), dtype=np.float)      
    #mstar
    obs_data_array[:,0] = data[:,8]-0.24
    #r-i
    obs_data_array[:,1] = data[:,5]-data[:,6]
    
    #mass
    obs_data_array[:,1] = data[:,8]

    return obs_data_array, 'SDSS BOSS DR12', False

def red_blue_cut(load_from_file=False):
    
    obs_data_array = np.zeros((10, 7), dtype=np.float32)

    x_data=np.linspace(0.0,3,10)
    
    obs_data_array[:,0]=x_data   

    obs_data_array[:,1]=-x_data+2.35 
    
    #print obs_data_array

    return obs_data_array, '$g-i>2.35$', False
  

def Guo13_cut(load_from_file=False):
    
    obs_data_array = np.zeros((10, 7), dtype=np.float32)

    x_data=np.linspace(0.0,3,10)

       
    obs_data_array[:,0]=x_data   
  
    obs_data_array[:,1]=0.679-0.082*(-24.69+20.0)  
    
    #print obs_data_array

    return obs_data_array, 'Guo+13', False

def Guo13_cut2(load_from_file=False):
    
    obs_data_array = np.zeros((10, 7), dtype=np.float32)

    x_data=np.linspace(0.0,3,10)

       
    obs_data_array[:,0]=x_data   
  
    obs_data_array[:,1]=0.679-0.082*(-20.37+20.0)  
    
    #print obs_data_array

    return obs_data_array, 'Guo+13', False

def Rodriguez15BMstar2Mhalo(load_from_file=False):
    #model Mstar/Mhalo
    
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/Rodriguez+15_Fig13_BigMD-BOSS_LC_asymmetric_errorbars.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float)  
  
    obs_data_array[:,0] = np.log10(data_array[:,0])#+np.log10(0.70/0.6777)+np.log10(0.6777)-0.03925
    obs_data_array[:,1] = np.log10(data_array[:,1])#-3*np.log10(0.70/0.6777) - 3*np.log10(0.6777)
    obs_data_array[:,4] = (data_array[:,4]/data_array[:,1])*0.434
    obs_data_array[:,5] = (data_array[:,5]/data_array[:,1])*0.434
    
    #print obs_data_array
    
    return obs_data_array, 'Rodriguez-Torres+16', False

##############################################################################################################################################################
def SMF_CMASS_Granada_055(load_from_file=False): return Rodriguez15SMF_CMASS('gra', '0.55', info='tot', legend='Granada CMASS S')
def SMF_CMASS_Wisconsin_055(load_from_file=False): return Rodriguez15SMF_CMASS('wis', '0.55', info='tot', legend='Wisconsin CMASS S')
def SMF_CMASS_Portsmouth_055(load_from_file=False): return Rodriguez15SMF_CMASS('por', '0.55', info='tot', legend='Portsmouth CMASS S')

def Rodriguez15SMF_CMASS(catname,
                         redshift,
                         info='',
                         legend=''):
        
    filename=mycomp+'anaconda/pro/data/CMASS_SAMS/SM_DISTRIBUTIONS/MF_'+catname+'_cmassport_Planck_'+redshift+'_'+info+'.dat'        
    print 'load CMASS from Rodirguez-Torres et al. (2015)', filename

    data_array = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=0)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)
     
    obs_data_array[:,0] = data_array[:,0]+np.log10(0.6777)-0.03925  #Corrected from Kroupa to Charbier Lacey+16 from Violeta
    obs_data_array[:,1] = np.log10(data_array[:,1]) - 3*np.log10(0.6777)
    obs_data_array[:,4] = (data_array[:,2]/data_array[:,1])*0.434
    obs_data_array[:,5] = (data_array[:,2]/data_array[:,1])*0.434
    
    return obs_data_array, legend, False
##############################################################################################################################################################
    
def Rodriguez15SMF(load_from_file=False):
    #model   
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/Rodriguez+15_model_black_line.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)
  
    obs_data_array[:,0] = data_array[:,0]+np.log10(0.6777)+np.log10(0.6777)-0.03925  #Corrected from Kroupa to Charbier Lacey+16 from Violeta
    obs_data_array[:,1] = np.log10(data_array[:,1])-3*np.log10(0.6777) - 3*np.log10(0.6777)
    
    return obs_data_array, 'Rodriguez+16 BigMD', False

def Rodriguez16_SMHM_obs_plus_dy(load_from_file=False):
        
    filename=mycomp+'anaconda/pro/OBS/Rodriguez-Torres+16_Fig13_Shan+15_-dy.txt'        
    minus_dy = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)

    filename=mycomp+'anaconda/pro/OBS/Rodriguez-Torres+16_Fig13_Shan+15_+dy.txt'        
    plus_dy = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)
    
    
    obs_data_array = np.zeros((plus_dy[:,0].size,7), dtype=np.float32)
     
    obs_data_array[:,0] = np.log10(plus_dy[:,0])
    obs_data_array[:,1] = (plus_dy[:,1]-minus_dy[:,1][-1::])
    obs_data_array[:,4] = np.log10(plus_dy[:,1])
    obs_data_array[:,5] = np.log10(minus_dy[:,1][::-1]) 
    print obs_data_array[:,[0,1,4,5]]
    
    return obs_data_array, 'Shan+17', False

def Rodriguez16_SMF_z055_plus_dy(load_from_file=False):
        
    filename=mycomp+'anaconda/pro/OBS/Rodriguez-Torres+16_Fig5_0.55_-dy.txt'        
    minus_dy = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=0)

    filename=mycomp+'anaconda/pro/OBS/Rodriguez-Torres+16_Fig5_0.55_+dy.txt'        
    plus_dy = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=0)
    
    
    obs_data_array = np.zeros((plus_dy[:,0].size,7), dtype=np.float32)
     
    obs_data_array[:,0] = plus_dy[:,0]
    obs_data_array[:,1] = np.log10((plus_dy[:,1]-minus_dy[:,1][-1::]))
    obs_data_array[:,4] = np.log10(plus_dy[:,1])
    obs_data_array[:,5] = np.log10(minus_dy[:,1][::-1]) 
    print obs_data_array[:,[0,1,4,5]]
    
    return obs_data_array, 'BigMD LC z=0.55', False

def Rodriguez15_CMASS_wp(load_from_file=False):
    #observational data
    #data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/data/CMASS_SAMS/CMASS_CLUSTERING/cmass-dr12v4-N-Reid-errors.wp', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=0)

    #unit test
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Rodriguez-Torres+16_Fig10_CMASS_DR12_wp.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)
  
    obs_data_array[:,0] = np.log10(data_array[:,0]/0.6777)
    obs_data_array[:,1] = np.log10(data_array[:,1]*data_array[:,0]/0.6777**2)
    obs_data_array[:,4] = data_array[:,2]/data_array[:,1]/0.6777*0.434
    obs_data_array[:,5] = obs_data_array[:,4]
    
    return obs_data_array, 'CMASS DR12', False

def Rodriguez15_CMASS_xi(load_from_file=False):
    #data
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/data/CMASS_SAMS/CMASS_CLUSTERING/cmass-dr12v4-N-errors.2Dbin5.corr.mono', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=0)    
    #data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/data/CMASS_SAMS/CMASS_CLUSTERING/dr12v4-CMASS-NGC-log-errors.0.01-80Mpc.corr.mono', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=0)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)
  
    obs_data_array[:,0] = data_array[:,0]*0.6777
    obs_data_array[:,1] = data_array[:,1]*(data_array[:,0]*0.6777)**2
    obs_data_array[:,4] = data_array[:,2]*(data_array[:,0]*0.6777)**2
    obs_data_array[:,5] = obs_data_array[:,4]
    #print obs_data_array
        
    return obs_data_array, 'CMASS DR12', False

def Rodriguez15_BigMD_LC_wp(load_from_file=False):

    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/data/CMASS_SAMS/HAM_CLUSTERING/bigMD-cmass-dr12v4-RST-standHAM-Vpeak-errors.wp', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=0)
    print data_array
    obs_data_array = np.zeros((19,7), dtype=np.float32)

    data_size=data_array[:,0].size
    obs_data_array[0:data_size,0] = data_array[:,0]/0.6777
    obs_data_array[0:data_size,1] = data_array[:,1]*data_array[:,0]/0.6777**2
    obs_data_array[0:data_size,4] = np.log10(obs_data_array[0:data_size,1]) + (data_array[:,2]/data_array[:,1]/0.6777*0.434)
    obs_data_array[0:data_size,5] = np.log10(obs_data_array[0:data_size,1]) - (data_array[:,2]/data_array[:,1]/0.6777*0.434)
    obs_data_array[data_size::,:] = 'nan'  
    filename_out = mycomp+'anaconda/pro/myRun/histos/twoPCF/Galacticus_CMASS/twoPCF_Z_HAM_z_0.56_BigMD_LC_wp.txt'

    myOutput.writeIntoFile(filename_out, 
                           obs_data_array,
                           myheader='Rodriguez+16 BigMD z=0.55 /CMASS_SAMS/HAM_CLUSTERING/bigMD-cmass-dr12v4-RST-standHAM-Vpeak-errors.wp\n(1) rmin [Mpc]  (2) rmax [Mpc]  (3) r_avg [Mpc] (4) r*wp [Mpc2] (5) xi error min (6) xi error max (7) -',
                           data_format="%0.8f",
                           mydelimiter='\t')
    exit()    
    return obs_data_array, 'Rodriguez+16 HAM BigMD', False

def Rodriguez15_BigMD_LC_xi(load_from_file=False):
    #data 
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/data/CMASS_SAMS/HAM_CLUSTERING/bigMD-cmass-dr12v4-RST-standHAM-Vpeak-errors.2Dbin5.corr.mono', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float, skiprow=0)

    obs_data_array = np.zeros((72,7), dtype=np.float32)
    data_size=data_array[:,0].size  
    obs_data_array[0:data_size,0] = data_array[:,0]
    obs_data_array[0:data_size,1] = data_array[:,1]*data_array[:,0]**2
    obs_data_array[0:data_size,4] = obs_data_array[0:data_size,1] + data_array[:,2]*data_array[:,0]**2
    obs_data_array[0:data_size,5] = obs_data_array[0:data_size,1] - data_array[:,2]*data_array[:,0]**2
    obs_data_array[data_size::,:] = 'nan'

    filename_out = mycomp+'anaconda/pro/myRun/histos/plotXY/twoPCF/twoPCF_Z_HAM_z_0.56_BigMD_LC-cmass_xi.txt'
    myOutput.writeIntoFile(filename_out, 
                           obs_data_array,
                           myheader='Rodriguez+16 BigMD z=0.55 /CMASS_SAMS/HAM_CLUSTERING/bigMD-cmass-dr12v4-RST-standHAM-Vpeak-errors.2Dbin5.corr.mono\n(1) rmin [h-1Mpc]  (2) rmax [h-1Mpc]  (3) r_avg [h-1Mpc] (4) r^2*xi (5) xi error min (6) xi error max (7) -',
                           data_format="%0.8f",
                           mydelimiter='\t')
    
    return obs_data_array, 'Rodriguez+16 HAM BigMD', False

def Rodriguez15PRIMUS(load_from_file=False):
    #observations
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/Rodriguez+15_PRIMUS_z_0.5_0.6.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float)
  
    obs_data_array[:,0] = data_array[:,0]-0.03925 #Corrected from Kroupa to Charbier Lacey+16 from Violeta
    obs_data_array[:,1] = np.log10(data_array[:,1])
    
    return obs_data_array, 'PRIMUS 0.5<z<0.6', False


def Salim14(load_from_file=False):

    reference = {'name': 'Salim+14', 
                'filename': 'salim_14', 
                'legend': 'SDSS DR7 z=0.0', 
                'paper': 'Salim et al. ApJ 797:126 (19pp), 2014', 
                'little_h': 0.7,
                'IMF_corr': 0.0}

    data = myData.readAnyFormat(config=False, 
                                mypath=mycomp+'anaconda/pro/OBS/CalibrationData/DATA/'+reference['filename'], 
                                data_format='ASCII', 
                                nr_col=2, 
                                data_shape='shaped', 
                                delim=' ', 
                                mydtype=np.float32, 
                                skiprow=2)

    
    #Correct from Salpeter+55 to Charbier+03 --> Refernece from Violeta!    
    data[:,0]=10**(np.log10(data[:,1]/reference['little_h'])+reference['IMF_corr'])
    data[:,1] = data[:,3]/data[:,0] 

    return ssfr2mstar_ref(reference,data), False
                
def Santini09(load_from_file=False):
                
    reference = {'name': 'Santini+09', 
                'filename': 'santini_09', 
                'legend': 'GOODS-Music z=0.0', 
                'paper': 'Santini et al. AA 504, 751-767 (2009)', 
                'little_h': 0.7,
                'IMF_corr': -0.24}

    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/CalibrationData/DATA/'+reference['filename'], data_format='ASCII', nr_col=2, data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=2)

                                           
    return ssfr2mstar_ref(reference,data), False
 
def SDSS_DR7_MHU_wp(load_from_file=False):

    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/twoPCF_SDSS_DR7_z_0.1_wp_mocks.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=2)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float)
  
    obs_data_array[:,0] = np.log10(data_array[:,0])
    obs_data_array[:,1] = np.log10(data_array[:,1])
    
    return obs_data_array, 'SDSS DR7 z~0.1', False


def SDSS_ELG_z01(load_from_file=False):

    obs_data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/SDSS_ELG_z_0.1_from_Ginevra.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)
    
    return obs_data_array, 'SDSS ELG', False    

def ssfr2mstar_ref(reference,
                   data):                            
    
    print reference['name'], 'data.shape:', data.shape, 'little-h:', reference['little_h'], 'IMF correction:', reference['IMF_corr']

    myFuncs = mF.MyFunctions()
    
    data=data[np.where(data[:,1] > 1.0e-11)[:][0]]
 
    data = data[np.where(data[:,0]>=np.percentile(data[:,0],0.01))[:][0]]
    data = data[np.where(data[:,0]<np.percentile(data[:,0],99.9))[:][0]]                          

    print data
                                        
    mybinned_array, binsize = myFuncs.binUp(data,
                                            20,
                                            histo_min='min',
                                            histo_max='max',
                                            log10bin=True,
                                            binup_2D=True,
                                            use_MAD=True,
                                            norm_by_binsize=False,
                                            equal_bins=False)
    
    #print mybinned_array

    obs_data_array = np.zeros((mybinned_array[:,0].size,7), dtype=np.float32)  
    #Errorbars in log(y) vs log(x) plot
    obs_data_array[:,0] = np.log10(mybinned_array[:,0])
    obs_data_array[:,1] = np.log10(mybinned_array[:,1])
    obs_data_array[:,2] = mybinned_array[:,2]
    obs_data_array[:,3] = mybinned_array[:,2]    
    
    #This is the only way to do log-error bars --> with this calculation the log-error bars in a plot
    #where the quantities are log-values, the error bars are always SYMETRIC
    
    #shaded region corresponding to error bars
    obs_data_array[:,4] = obs_data_array[:,1] + ((mybinned_array[:,4]/mybinned_array[:,1]))*0.43
    obs_data_array[:,5] = obs_data_array[:,1] - ((mybinned_array[:,5]/mybinned_array[:,1]))*0.43

    #error bars                
    obs_data_array[:,4] = ((mybinned_array[:,4]/mybinned_array[:,1]))*0.43
    obs_data_array[:,5] = ((mybinned_array[:,4]/mybinned_array[:,1]))*0.43
               
#    filename_out = mycomp+'anaconda/pro/OBS/'+reference['legend0']+'_'+reference['filename0']+'_sSFR2Mstar.txt'
#    myOutput.writeIntoFile(filename_out, 
#                           mybinned_array,
#                           myheader='#'+reference['paper0']+', pbwo5rks.CosmicCarnage2015 file: "/CalibrationData/PROPER_DATA/'reference['filneame0']+'_masssfr.dat" \n(1) Mstar [Msun], (2) sSFR [yr-1], (3) binsize, (4) Phi error (-),  (5) Phi error (+)',
#                           data_format="%0.8e",
#                           mydelimiter='\t')

    #print obs_data_array[:, [0,1]]

    return obs_data_array, reference['legend'], False

def ssfr2mstar_SDSS_DR8(load_from_file=False):


    redshift='0.0'
    data = myData.readAnyFormat(config=False, 
                                myfilename=mycomp+'anaconda/pro/data/SDSS_DR8/SDSS_DR8_z_'+redshift+'_Portsmouth_starforming_Salpeter.hdf5',
                                nr_rows=1000000, 
                                nr_col=2, 
                                data_format='HDF5', 
                                mydtype=np.float32,
                                id_col_array={'name0': 'mstar', 'name1': 'sfr', 'mstar_col_id': 0, 'sfr_col_id': 1})

    data[:,0]-=0.24
    
    print data[0:10,:]
    data[:,1]=data[:,1]/data[:,0]
          
    data = myData.selectData2Compute(data, 
                                  selected_col=1, 
                                  operator='>', 
                                  condition=0.0)

    data = myData.selectData2Compute(data, 
                                  selected_col=0, 
                                  operator='>', 
                                  condition=1e7)
                                        
    mybinned_array, binsize = myFuncs.binUp(data, 
                                    20,
                                    histo_min=min(data[:,0]),
                                    histo_max=max(data[:,0]),
                                    log10bin=True,
                                    binup_2D=True,
                                    use_median=True)

    obs_data_array = np.zeros((mybinned_array[:,0].size,7), dtype=np.float32)  
    #Errorbars in log(y) vs log(x) plot
    obs_data_array[:,0] = np.log10(mybinned_array[:,0])
    obs_data_array[:,1] = np.log10(mybinned_array[:,1])
    obs_data_array[:,2] = mybinned_array[:,2]
    obs_data_array[:,3] = mybinned_array[:,2]    
    
    #This is the only way to do log-error bars --> with this calculation the log-error bars in a plot
    #where the quantities are log-values, the error bars are always SYMETRIC
    
    #shaded region corresponding to error bars
    obs_data_array[:,4] = obs_data_array[:,1] + ((mybinned_array[:,3]/mybinned_array[:,1]))*0.43
    obs_data_array[:,5] = obs_data_array[:,1] - ((mybinned_array[:,3]/mybinned_array[:,1]))*0.43

    #error bars                
    obs_data_array[:,4] = ((mybinned_array[:,3]/mybinned_array[:,1]))*0.43
    obs_data_array[:,5] = ((mybinned_array[:,3]/mybinned_array[:,1]))*0.43
               
    filename_out = mycomp+'anaconda/pro/OBS/SDSS_DR8_'+redshift+'_sSFR2Mstar.txt'
    myOutput.writeIntoFile(filename_out, 
                           mybinned_array,
                           myheader='#SDSS DR8 Porthsmouth starforming salpeter 26 download 11/28/2016 "http://www.sdss.org/dr12/spectro/galaxy_portsmouth/" \n(1) Mstar(Charbier) [Msun] (2) sSFR [yr-1] (3) dx (4) -dy  (5) +dy (6) count',
                           data_format="%0.8e",
                           mydelimiter='\t')


    return obs_data_array, 'SDSS DR8 z='+redshift, False

def Tremonti_04_OH_mstar(load_from_file=False):
  
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/Tremonti+04_OH_Mstar_relation.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=6)

    obs_data_array = np.zeros((data_array[:,0].size,10), dtype=np.float)
    
    obs_data_array[:,0] = data_array[:,0]-0.03925 #Corrected from Kroupa to Charbier Lacey+16 from Violeta
    obs_data_array[:,1] = data_array[:,3]
    obs_data_array[:,4] = obs_data_array[:,1]-data_array[:,5]
    obs_data_array[:,5] = obs_data_array[:,1]-data_array[:,5]
    
    #print obs_data_array[:,0:5]

    return obs_data_array, 'SDSS z~0.1', False

##############################################################################################################################################################
def SAGE_load_add(load_from_file=False): return load_the_shit2('SAGE_1Gpc', 'SMF', '0.09', info='_v3_mstar+IC', legend='SAGE mstar+IC')

def SAGE_load_add2(load_from_file=False): return load_the_shit2('SAGE_1Gpc', 'SMF', '0.09', info='_v3_IC', legend='SAGE IC')

def Gal_ssfr2mstar(load_from_file=False): return load_the_shit2('Galacticus_1Gpc', 'ssfr2mstar', '0.0', info='_stf_cut_1e-11', legend='Gal')

def Gal2_ssfr2mstar(load_from_file=False): return load_the_shit2('Galacticus_1Gpc_run2', 'ssfr2mstar', '0.0', info='_tarsel', legend='Gal2')
    
def SAG_ssfr2mstar(load_from_file=False): return load_the_shit2('SAG_1Gpc',  'ssfr2mstar', '0.0', info='_stf_cut_1e-11_v3')

def SAGE_ssfr2mstar(load_from_file=False): return load_the_shit2('SAGE_1Gpc', 'ssfr2mstar', '0.0', info='_stf_cut_1e-11_v3')
 
def Gal_mbh2mstarsph(load_from_file=False): return load_the_shit2('Galacticus_1Gpc', 'mbh2mstarsph', '0.0') 

def SAG_mbh2mstarsph(load_from_file=False): return load_the_shit2('SAG_1Gpc', 'mbh2mstarsph', '0.0', info='_v3') 

def SAGE_mbh2mstarsph(load_from_file=False): return load_the_shit2('SAGE_1Gpc', 'mbh2mstarsph', '0.0', info='_v3')

def Gal_mcold2mstar(load_from_file=False): return load_the_shit2('Galacticus_1Gpc', 'mcold2mstar', '0.0')

def SAG_mcold2mstar(load_from_file=False): return load_the_shit2('SAG_1Gpc', 'mcold2mstar', '0.0', info='_v3') 

def SAGE_mcold2mstar(load_from_file=False): return load_the_shit2('SAGE_1Gpc', 'mcold2mstar', '0.0', info='_v3')

def Gal_zgas2mstar(load_from_file=False): return load_the_shit2('Galacticus_1Gpc', 'zgas2mstar', '0.09')

def SAG_zgas2mstar(load_from_file=False): return load_the_shit2('SAG_1Gpc', 'zgas2mstar', '0.09', info='_v3') 

def SAGE_zgas2mstar(load_from_file=False): return load_the_shit2('SAGE_1Gpc', 'zgas2mstar', '0.09', info='_v3')    

#2PCF
##############################################################################################################################################################
band='_-20_Mr_-19'

catname='Galacticus_1Gpc'
key_info='_new_posCor_MAB'
catname='SAG_1Gpc_v2'
key_info='_v3_new'

def twoPCF_cents(load_from_file=False): return load_the_shit2(catname, 'twoPCF', '0.09', info=key_info+band+'_Pi60_centrals_wp', legend='centrals', key='Mr/')   
def twoPCF_satsOrphs(load_from_file=False): return load_the_shit2(catname, 'twoPCF', '0.09', info=key_info+band+'_Pi60_sats+orphans_wp', legend='satellites+orphans', key='Mr/') 
def twoPCF_sats(load_from_file=False): return load_the_shit2(catname, 'twoPCF', '0.09', info=key_info+band+'_Pi60_sats_wp', legend='satellites', key='Mr/')
def twoPCF_orphs(load_from_file=False): return load_the_shit2(catname, 'twoPCF', '0.09', info=key_info+band+'_Pi60_orphans_wp', legend='orphans', key='Mr/')
def oh2mstar(load_from_file=False): return load_the_shit2(catname, 'oh2mstar', '0.07', info='_OII_disk_sfr+h', legend='SAGv2 OII + 0.25', key='')

def load_the_shit(load_from_file=False):

    #filename=mycomp+'/anaconda/pro/myRun/histos/HOD/HOD_Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3_HOD_test_HOD_mhalo_cents_200c'
    #filename=mycomp+'/anaconda/pro/myRun/histos/HOD/HOD_Galacticus_1Gpc_run2_z_0.56_tarsel_CMASS_down_sample3_mhalo_cents_200c'
    redshift='3.04'
    
    for item in ['CUT1_mstar', 'CUT2_mstar', 'CUT3_mstar', 'CUT1_sfr', 'CUT2_sfr', 'CUT3_sfr', 'CUT1_ssfr', 'CUT2_ssfr', 'CUT3_ssfr']:
        #filename=mycomp+'/anaconda/pro/myRun/histos/HOD/HOD_Galacticus_400Mpc_z_'+redshift+'_tarsel_rand-sample_n_0.008_mhalo_cents_200c'
        #filename=mycomp+'/anaconda/pro/myRun/histos/HOD/HOD_Galacticus_1Gpc_z_'+redshift+'_tarsel_rand-sample_n_0.016_mhalo_cents_200c'
        filename=mycomp+'/anaconda/pro/myRun/histos/HOD/HOD-SF/HOD_Galacticus_1Gpc_z_'+redshift+'_tarsel_HOD-SF_'+item+'_mhalo_cents_200c'    
        #filename=mycomp+'/anaconda/pro/myRun/histos/HOD/HOD-SF/HOD_Galacticus_400Mpc_z_'+redshift+'_tarsel_HOD-SF_'+item+'_mhalo_cents_200c'
        #filename=mycomp+'/anaconda/pro/myRun/histos/HOD/HOD-SF/HOD_IllustrisTNG300_z_'+redshift+'_tarsel_HOD-SF_'+item+'_mhalo_200c'
        
        data_all = myData.readAnyFormat(config=False, mypath=filename+'_all.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=2)
        
        data_centrals = myData.readAnyFormat(config=False, mypath=filename+'_centrals.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=2)    
        
        data_all[:,1]=data_all[:,1]-data_centrals[:,1]
        
        header='\n(1) mhalo_200c [Msun] (2) <ngal> [-] (3) -dx \t(4) dx\t(5) N unique\t(6) N sats (7) Nhalos M200c from cosmosim.org' # from IllustrisTNG300 FULLHalos cats (Celeste)'
        #header='\n(1) mhalo_200c [Msun] (2) <ngal> [-] (3) -dx \t(4) dx\t(5) N unique\t(6) N sats (7) from IllustrisTNG300 FULLHalos cats (Celeste)'
    
        myOutput.writeIntoFile(filename+'_sats.txt',
                               data_all,
                               myheader='HODFunction MDPL2 Galacticus 1h-1Gpc z='+redshift+' cumulative: no'+header,
                               #myheader='HODFunction SMDPL Galacticus 400h-1Gpc z='+redshift+' cumulative: no'+header,
                               #myheader='HODFunction Illustris TNG300-1 205h-1Gpc z='+redshift+' cumulative: no'+header,
                               data_format="%0.8e",
                               mydelimiter='\t')
    exit()    

def load_the_shit2(catname,
                   plot_key,
                   redshift,
                   info='',
                   legend='binned function',
                   key=''):

    filename=mycomp+'/anaconda/pro/myRun/histos/'+plot_key+'/'+key+plot_key+'_'+catname+'_z_'+redshift+info+'.txt'        
    print 'load binned-function:', filename
    data_array = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=2)
    #print data_array
    #if plot_key.find('zgas')!=-1:
    #only take every second datopoint
    i=0
    while i<data_array[:,0].size:
        data_array[i,0]=-1      
        i+=2

    if plot_key.find('zgas')==-1: 
        obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float)
    else:
        obs_data_array = np.zeros((data_array[:,0].size,10), dtype=np.float)
        obs_data_array[:,6] = np.log10(data_array[:,6])  

    if plot_key.find('twoPCF')!=-1:
        obs_data_array[:,0] = np.log10(data_array[:,2])
        obs_data_array[:,1] = np.log10(data_array[:,3])
    else:
                  
        obs_data_array[:,0] = np.log10(data_array[:,0])
        if plot_key.find('zgas')==-1 and plot_key.find('oh')==-1:
            obs_data_array[:,1] = np.log10(data_array[:,1])
            obs_data_array[:,4] = (data_array[:,4]/data_array[:,1])*0.43    
            obs_data_array[:,5] = obs_data_array[:,4]#(data_array[:,5]/data_array[:,1])*0.43
        else:
            obs_data_array[:,1] = data_array[:,1] #+0.25     
            obs_data_array[:,4] = data_array[:,4]   
            obs_data_array[:,5] = data_array[:,5]

            obs_data_array[:,4] = obs_data_array[:,1] + data_array[:,4]    
            obs_data_array[:,5] = obs_data_array[:,1] - data_array[:,5]
         
    #print obs_data_array
   
    return obs_data_array, legend, False      

def load_the_shit3(load_from_file=False):

    filename=mycomp+'/anaconda/pro/data/Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS'
        
    data_mass=load_catalog(filename+'_mass_sample_200c_mstar_gt_4e11.hdf5',
                     'haloid',
                     'hostid',
                     name_add='orphan',
                     name_add2='mAB_dA_total_i',
                     name_add3='mstar',
                     name_add4='MAB_dA_total_r')

    data_cmass=load_catalog(filename+'_200c_mstar_gt_4e11.hdf5',
                     'haloid',
                     'hostid',
                     name_add='orphan',
                     name_add2='mAB_dA_total_i',
                     name_add3='mstar',
                     name_add4='MAB_dA_total_r')

    test = np.in1d(data_mass['hostid'],data_cmass['hostid'])
    data=data_mass[np.where(test==False)[:][0]]
    print data
    exit()

def format_MTs_Rockstar(load_from_file=False):

    path=mycomp+'/anaconda/pro/data/ROCKSTAR/MDPL2_merger_tree_info_300-r0002.csv'
    path=mycomp+'/anaconda/pro/data/ROCKSTAR/rtest2.csv'
    IDs = myData.readAnyFormat(config=False, mypath=path, data_format='ASCII', data_shape='shaped', delim=',', mydtype=np.str_, skiprow=1)
    cluster_name='rtest2'#path[path.find('-r')+2: -4]

    for i, item in enumerate(IDs[:,2]):
        #print item, i
        IDs[i,2]=format(mL.expfactor_to_redshift(float(IDs[i,2])), '0.2f')

    myOutput.writeIntoFile(mycomp+'/anaconda/pro/data/ROCKSTAR/MDPL2_merger_tree_info_300-r'+cluster_name+'.txt',
                           IDs[:,[1,2,3]],
                           myheader='MDPL2.Rockstar main progenitor of central halo from 300-cluster region: '+cluster_name+' from cosmosim.org\n(1) snapnum (2) z (3) haloid/rockstarId',
                           data_format="%s",
                           mydelimiter='\t')



    exit()        

##############################################################################################################################################################  

def Reid14_do_fit():

    #MedRes    
    Mcut        = (10**13.159)/0.6777
    log10Mmin   = 13.011-np.log10(0.6777)
    M1          = (10**14.068)/0.6777
    sigmaM      = 0.358
    alpha       = 0.89

    #COMV   - lowest n_CMASS  
    #Mcut        = (10**11.89)/0.6777
    #log10Mmin   = 12.983-np.log10(0.6777)
    #M1          = (10**14.23)/0.6777
    #sigmaM      = 0.31
    #alpha       = 1.15

    #data = np.linspace(11,15, 30) - np.log10(0.6777)
    data = np.arange(12,15.3458,0.1034)  
    #print 'Reid+14:', data
    
    from scipy.special import erf as erf
    
    f=0.5*(1.0+erf((data-log10Mmin)/sigmaM))

    g=(((10**data - Mcut)/M1)**alpha)
    
    g[np.where(np.isfinite(g)==False)[:][0]]=0.0

    e=(f+g*(data>np.log10(Mcut)))
    
    return data,f,g,e

#MD Reid resc=1.31   
resc=1.31

def Reid14_HOD_fit_all(load_from_file=False):
           
    data_x,f,g,e = Reid14_do_fit()
    obs_data_array = np.zeros((data_x.size, 7), dtype=np.float32)

    obs_data_array[:,0] = data_x
    #/1.31 is the correction to number density of CMASS at z=0.56 (Reid+14 uses maximum n with 4.23e-4 h3Mp3)
    obs_data_array[:,1] = np.log10((f+g)/resc)

    return obs_data_array, 'Reid+14', False    
    
def Reid14_HOD_fit_cents(load_from_file=False):
    
    data_x,f,g,e = Reid14_do_fit()   
    obs_data_array = np.zeros((data_x.size, 7), dtype=np.float32)

    obs_data_array[:,0] = data_x
        #/1.31 is the correction to number density of CMASS at z=0.56 (Reid+14 uses maximum n with 4.23e-4 h3Mp3)
    obs_data_array[:,1] = np.log10(f/resc)

    return obs_data_array, 'Reid+14 centrals', False  

def Reid14_HOD_fit_sats(load_from_file=False):
    
    data_x,f,g,e = Reid14_do_fit()  
    obs_data_array = np.zeros((data_x.size, 7), dtype=np.float32)

    obs_data_array[:,0] = data_x
        #/1.31 is the correction to number density of CMASS at z=0.56 (Reid+14 uses maximum n with 4.23e-4 h3Mp3)
    obs_data_array[:,1] = np.log10(g/resc)

    return obs_data_array, 'Reid+14 sats', False  

def R_I_mstar_zoomOut(load_from_file=False):
    
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/r-i_vs_mstar_zoomed_region.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=1)

    obs_data_array = np.zeros((data_array[:,0].size, 7), dtype=np.float32)
    
    obs_data_array[:,0]=data_array[:,0]
    obs_data_array[:,1]=data_array[:,1]

    return obs_data_array, ' ', False

def BigMDPL_LC_HOD_all(load_from_file=False):
    #from Sergio
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/BigMDPL-LC_CMASS-hod.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=1)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)
  
    obs_data_array[:,0] = data_array[:,0]
    obs_data_array[:,1] = np.log10(data_array[:,3])
    obs_data_array[:,4] = (1.0/np.sqrt(data_array[:,6]))/(data_array[:,3])*0.434
    obs_data_array[:,5] = obs_data_array[:,4]
    #print 'BigMD-LC:', obs_data_array[:,0]
    return obs_data_array, 'BigMD-LC', False 

def BigMDPL_LC_HOD_cents(load_from_file=False):
    #from Sergio
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/BigMDPL-LC_CMASS-hod.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=1)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)
  
    obs_data_array[:,0] = data_array[:,0]
    obs_data_array[:,1] = np.log10(data_array[:,2])
    obs_data_array[:,4] = (1.0/np.sqrt(data_array[:,5]))/(data_array[:,2])*0.434
    obs_data_array[:,5] = obs_data_array[:,4]
    
    return obs_data_array, 'BigMD-LC', False

def BigMDPL_LC_HOD_sats(load_from_file=False):
    #from Sergio
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/BigMDPL-LC_CMASS-hod.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=1)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)
  
    obs_data_array[:,0] = data_array[:,0]
    obs_data_array[:,1] = np.log10(data_array[:,1])
    obs_data_array[:,4] = (1.0/np.sqrt(data_array[:,4]))/(data_array[:,1])*0.434
    obs_data_array[:,5] = obs_data_array[:,4]
    
    return obs_data_array, 'BigMD-LC', False

def MD_LC_HOD_all(load_from_file=False):
    #from Sergio
    #first try with scatter = 0.2
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/MD-LC_CMASS-hod_0.56_3.30e-4.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=1)
    #second try with scatter = 0.4
    #data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/MD-LC_CMASS-hod_0.56_3.30e-4_second.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=1)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)
  
    obs_data_array[:,0] = data_array[:,0]
    obs_data_array[:,1] = np.log10(data_array[:,3])
    obs_data_array[:,4] = (1.0/np.sqrt(data_array[:,6]))/(data_array[:,3])*0.434
    obs_data_array[:,5] = obs_data_array[:,4]
    #print 'BigMD-LC:', obs_data_array[:,0]
    return obs_data_array, 'MD-LC', False 

def MD_LC_HOD_cents(load_from_file=False):
    #from Sergio
    #first try with scatter = 0.2
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/MD-LC_CMASS-hod_0.56_3.30e-4.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=1)
    #second try with scatter = 0.4
    #data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/MD-LC_CMASS-hod_0.56_3.30e-4_second.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=1)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)
  
    obs_data_array[:,0] = data_array[:,0]
    obs_data_array[:,1] = np.log10(data_array[:,2])
    obs_data_array[:,4] = (1.0/np.sqrt(data_array[:,5]))/(data_array[:,2])*0.434
    obs_data_array[:,5] = obs_data_array[:,4]
    
    return obs_data_array, 'MD-LC', False

def MD_LC_HOD_sats(load_from_file=False):
    #from Sergio
    #first try with scatter = 0.2
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/MD-LC_CMASS-hod_0.56_3.30e-4.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=1)
    #second try with scatter = 0.4
    #data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/MD-LC_CMASS-hod_0.56_3.30e-4_second.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=1)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)
  
    obs_data_array[:,0] = data_array[:,0]
    obs_data_array[:,1] = np.log10(data_array[:,1])
    obs_data_array[:,4] = (1.0/np.sqrt(data_array[:,4]))/(data_array[:,1])*0.434
    obs_data_array[:,5] = obs_data_array[:,4]

    return obs_data_array, 'MD-LC', False

def BigMDPL_LC_HOD_all2(load_from_file=False):
    #from Sergio
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/BigMDPL-LC_CMASS-hod-correct.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=2)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)
  
    obs_data_array[:,0] = data_array[:,0]-np.log10(0.6777)
    obs_data_array[:,1] = np.log10(data_array[:,2]+data_array[:,3])
    #obs_data_array[:,4] = np.sqrt(data_array[:,2]+data_array[:,3])/(data_array[:,2]+data_array[:,3])*0.434
    #obs_data_array[:,5] = obs_data_array[:,4]
    #print 'BigMD-LC:', obs_data_array[:,0]
    return obs_data_array, 'BigMD-LC', False 

def BigMDPL_LC_HOD_cents2(load_from_file=False):
    #from Sergio
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/BigMDPL-LC_CMASS-hod-correct.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=2)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)
  
    obs_data_array[:,0] = data_array[:,0]-np.log10(0.6777)
    obs_data_array[:,1] = np.log10(data_array[:,2])
    
    return obs_data_array, 'BigMD-LC', False

def BigMDPL_LC_HOD_sats2(load_from_file=False):
    #from Sergio
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/BigMDPL-LC_CMASS-hod-correct.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=2)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)
  
    obs_data_array[:,0] = data_array[:,0]-np.log10(0.6777)
    obs_data_array[:,1] = np.log10(data_array[:,3])

    return obs_data_array, 'BigMD-LC', False

def HOD_all_dummy(load_from_file=False):
    #from Sergio
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/BigMDPL-LC_CMASS-hod.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=2)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)
  
    return obs_data_array, 'centrals+sats', False   

def HOD_cents_dummy(load_from_file=False):
    #from Sergio
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/BigMDPL-LC_CMASS-hod.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=2)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)

    return obs_data_array, 'centrals', False

def HOD_sats_dummy(load_from_file=False):
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/BigMDPL-LC_CMASS-hod.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=2)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)
    return obs_data_array, 'sats', False   

def CUT1_dummy(load_from_file=False):
  
    return np.zeros((2,7), dtype=np.float32), 'CUT1', False   

def CUT2_dummy(load_from_file=False):

    return np.zeros((2,7), dtype=np.float32), 'CUT2', False

def CUT3_dummy(load_from_file=False):

    return np.zeros((2,7), dtype=np.float32), 'CUT3', False   

def White11_HOD_fit(load_from_file=False):

    Mcut        = 10**13.159
    log10Mmin   = 13.011
    M1          = 10**14.068
    sigmaM      = 0.358
    alpha       = 0.89
    
    #h(x)=0.5*erfc(log(10**13.08/10**x)/(sqrt(2)*0.98))
    #k(x)=h(x)*((10**x-1.13*10**13.08)/10**14.06)**0.9
    
    #j(x)=h(x)+k(x)*(x>log10(1.13*10**(13.08)))


def White11_wp(load_from_file=False):
    
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/White+11_BOSS_wp_table1.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=2)

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)
  
    obs_data_array[:,0] = np.log10(data_array[:,0]/0.70*0.6777)
    obs_data_array[:,1] = np.log10(data_array[:,1]/0.70*0.6777)

    obs_data_array[:,4] = obs_data_array[:,1] + (data_array[:,2]/data_array[:,1])*0.434
    obs_data_array[:,5] = obs_data_array[:,1] - (data_array[:,2]/data_array[:,1])*0.434
   
    #print obs_data_array
    return obs_data_array, 'White+11', False


def Zehavi11_18_17_wp(load_from_file=False):
    
    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/Zehavi+11_Table7_18_M_r_17.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=2)
    
    return load_Zehavi11_wp(data_array, '$-18<M_r<-17$')

def Zehavi11_19_18_wp(load_from_file=False):

    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/Zehavi+11_Table7_19_M_r_18.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=2)

    return load_Zehavi11_wp(data_array, '$-19<M_r<-18$')

def Zehavi11_20_19_wp(load_from_file=False):

    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/Zehavi+11_Table7_20_M_r_19.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=2)

    return load_Zehavi11_wp(data_array, '$-20<M_r<-19$')

def Zehavi11_21_20_wp(load_from_file=False):

    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/Zehavi+11_Table7_21_M_r_20.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=2)

    return load_Zehavi11_wp(data_array, '$-21<M_r<-20$')

def Zehavi11_22_21_wp(load_from_file=False): 

    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/Zehavi+11_Table7_22_M_r_21.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=2)
    
    return load_Zehavi11_wp(data_array, '$-22<M_r<-21$')

def Zehavi11_21_wp(load_from_file=False): 

    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/Zehavi+11_Table8_M_r_21.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=2)   

    return load_Zehavi11_wp(data_array, '$M_r>-22$')

def Zehavi11_22_wp(load_from_file=False): 

    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'/anaconda/pro/OBS/Zehavi+11_Table8_M_r_22.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=2)
    
    return load_Zehavi11_wp(data_array, '$M_r>-22$')

def load_Zehavi11_wp(data_array, legend):    

    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float32)
  
    obs_data_array[:,0] = data_array[:,0]/0.6777
    obs_data_array[:,1] = np.log10(data_array[:,1]/0.6777*obs_data_array[:,0])

    obs_data_array[:,4] = (data_array[:,2]/data_array[:,1])*0.434
    obs_data_array[:,5] = obs_data_array[:,4]
   
    #print obs_data_array
    return obs_data_array, legend, False

def ZFOURGE_150_z_250_SF(load_from_file=True):

    data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Tomazak+14_1.5_z_2.5_SF_best_to_fit_OII-SAMs.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)    
    
    obs_data_array = np.zeros((data_array[:,0].size,7), dtype=np.float)
    
    obs_data_array[:,0] = data_array[:,0]+2*np.log10(0.7)-2*np.log10(0.6777)          
    obs_data_array[:,1] = data_array[:,1]-3*np.log10(0.7)+3*np.log10(0.6777)
    obs_data_array[:,4] = -data_array[:,3]*0.434
    obs_data_array[:,5] = data_array[:,2]*0.434

    print obs_data_array
    
    return obs_data_array, 'ZFOURGE/CANDLES', False  

def loadDiverse(load_from_file=False):
    
    redshift='0.09'            
    filename = mycomp+'anaconda/pro/OBS/HMF_MDPL2_Rockstar_z_'+redshift+'_M200c_Msun_Nhalos.txt'
    print filename
    data = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=11)
    obs_data_array = np.zeros((data[:,0].size,7), dtype=np.float)
    obs_data_array[:,0] = 10**data[:,0]
    
    volume=(1000/0.6778)**3
    obs_data_array[:,1] = data[:,1]/0.1/volume
    
    obs_data_array[:,4] = obs_data_array[:,1] + obs_data_array[:,1]/(data[:,1]**0.5)
    obs_data_array[:,5] = obs_data_array[:,1] - obs_data_array[:,1]/(data[:,1]**0.5)

    obs_data_array[:,6] = data[:,1]
    
    filename_out = mycomp+'anaconda/pro/myRun/histos/HMF/HMF_MDPL2_Rockstar_z_'+redshift+'_Nhalos.txt'
    myOutput.writeIntoFile(filename_out, 
                               obs_data_array,
                               myheader='HMF MDPL2 Rockstar 1h-1Gpc z='+redshift+' cumulative: NO\n(1) mhalo_200c [Msun] (2) Phi [Mpc-3 dex-1]	(3) -dx	(4) dx	(5) -dy	(6) +dy	(7) N count',
                               data_format="%0.8e",
                               mydelimiter='\t')
    
    exit()
    return obs_data_array, legend, False


#Luminosity functions
##########################################################################################################################################################    
def loadFileLF(reference, band):
                 
    filename = mycomp+'anaconda/pro/OBS/'+reference+'_best_fit_LF_'+band+'.txt'
    filename = mycomp+'anaconda/pro/OBS/'+reference+'_data_LF_'+band+'.txt'
    print 'band:', band, filename   
    
    data=myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=6)    

    obs_data_array = np.zeros((data[:,0].size,6), dtype=np.float)
    
    obs_data_array[:,0] = data[:,0]+5*np.log10(0.6777)
    obs_data_array[:,1] = np.log10(10**data[:,1]*(0.6777**3))
    obs_data_array[:,3] = data[:,5]
    obs_data_array[:,4] = data[:,5]
     
    #print obs_data_array    
    
    return obs_data_array, reference, False

def Blanton05_LF_i(load_from_file=False):
    
   obs_data_array = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/Blanton+05_best_fit_LF_i.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=6)    
    
   obs_data_array[:,1]=np.log10(obs_data_array[:,1])
   
   return obs_data_array, 'Blanton+05', False

def Blanton03_LF_u(load_from_file=False):
    
   return loadFileLF('Blanton+03', 'u')

def Blanton03_LF_g(load_from_file=False):
    
    return loadFileLF('Blanton+03', 'g')   
    
def Blanton03_LF_r(load_from_file=False):
    
    return loadFileLF('Blanton+03', 'r')

def Blanton03_LF_i(load_from_file=False):
    
    return loadFileLF('Blanton+03', 'i')

def Blanton03_LF_z(load_from_file=False):
    
    return loadFileLF('Blanton+03', 'z')

def Montero09_LF_u(load_from_file=False):
    
   return loadFileLF('Montero+09', 'u')

def Montero09_LF_g(load_from_file=False):
    
    return loadFileLF('Montero+09', 'g')   
    
def Montero09_LF_r(load_from_file=False):
    
    return loadFileLF('Montero+09', 'r')

def Montero09_LF_i(load_from_file=False):
    
    return loadFileLF('Montero+09', 'i')

def Montero09_LF_z(load_from_file=False):
    
    return loadFileLF('Montero+09', 'z')
########################################################################################################################################################


#Luminosty functions SAGE
########################################################################################################################################################
def SAGE_u(load_from_file=False):
    
    return load_LF_SAGE('SAGE_1Gpc','u')

def SAGE_g(load_from_file=False):
    
    return load_LF_SAGE('SAGE_1Gpc','g')

def SAGE_r(load_from_file=False):
    
    return load_LF_SAGE('SAGE_1Gpc','r')
    
def SAGE_i(load_from_file=False):
    
    return load_LF_SAGE('SAGE_1Gpc','i')    

def SAGE_z(load_from_file=False):
    
    return load_LF_SAGE('SAGE_1Gpc','z')
    
def OBS_SAGE_u(load_from_file=False):
    
    return load_LF_SAGE('OBS','u')

def OBS_SAGE_g(load_from_file=False):
    
    return load_LF_SAGE('OBS','g')

def OBS_SAGE_r(load_from_file=False):
    
    return load_LF_SAGE('OBS','r')

def OBS_SAGE_i(load_from_file=False):
    
    return load_LF_SAGE('OBS','i')

def OBS_SAGE_z(load_from_file=False):
    
    return load_LF_SAGE('OBS','z')

def load_LF_SAGE(name,
                band):
    
    data = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/OBS/'+name+'_'+band+'-band_Adam.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float32, skiprow=6)
    obs_data_array = np.zeros((data[:,0].size,6), dtype=np.float)
   
    obs_data_array[:,0] = data[:,0]
    obs_data_array[:,1] = np.log10(data[:,1])
    obs_data_array[:,3] = data[:,4]
    obs_data_array[:,4] = data[:,5]
    
    return obs_data_array, name+' '+band+'-band', False
########################################################################################################################################################


#Z-Evolution
########################################################################################################################################################    
#Paper I
plot_key='zgas2mcold'
catname='Galacticus_1Gpc'
version={'0': 'tarsel_full_mstar_5e8_1e9',
          '1': 'tarsel_full_mstar_1e9_5e9',
          '2': 'tarsel_full_mstar_5e9_1e10',
          '3': 'tarsel_full_mstar_1e10_5e10',
          '4': 'tarsel_full_mstar_5e10_1e11',
          '5': 'tarsel_full_mstar_1e11'       
          }

legend={'0': '$8.7 < \log$ $M_{_*} < 9.0$',
          '1': '$9.0 < \log$ $M_{_*} < 9.7$',
          '2': '$9.7 < \log$ $M_{_*} < 10.0$',
          '3': '$10.0 < \log$ $M_{_*} < 10.7$',
          '4': '$10.7 < \log$ $M_{_*} < 11.0$',
          '5': '$\log$ $M_{_*} > 11.0$'       
          }

redshift={'0': 0.0,
          '1': 0.0,
          '2': 0.0,
          '3': 0.0,
          '4': 0.0,
          '5': 0.0
          }    

#OII-SFR Paper
#z=(0,0.1) 0.6,0.75,0.9,1.2

#plot_key='ssfr2mstar'
#Galacticus and SAGE
#redshift={'0': 0.09,
#          '1': 0.59,
#          '2': 0.74,
#          '3': 0.9,
#          '4': 1.22       
#          }
          
#SAG:
#redshift={'0': 0.07,
#          '1': 0.59,
#          '2': 0.7,
#          '3': 0.94
#          #'4': 1.22       
#          }     

def loadFile(plot_key,
             redshift,
             version,
             legend):
                 
    filename = mycomp+'anaconda/pro/myRun/histos/'+plot_key+'/'+plot_key+'_'+catname+'_z_'+redshift+'_'+version+'.txt'
    print filename
    obs_data_array = myData.readAnyFormat(config=False, mypath=filename, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2)
    
    obs_data_array[:,0] = np.log10(obs_data_array[:,0])
    
    if plot_key=='zgas2mcold':
        obs_data_array[:,1] = obs_data_array[:,1]
        
        obs_data_array[:,3] = obs_data_array[:,1] + obs_data_array[:,3]
        obs_data_array[:,3] = obs_data_array[:,1] - obs_data_array[:,3]
    else:
        obs_data_array[:,1] = np.log10(obs_data_array[:,1])
    
        obs_data_array[:,3] = obs_data_array[:,1] + (obs_data_array[:,3]/10**obs_data_array[:,1])*0.43
        obs_data_array[:,4] = obs_data_array[:,1] - (obs_data_array[:,4]/10**obs_data_array[:,1])*0.43

    return obs_data_array, legend, False

def zevol1(load_from_file=False):
    
    return loadFile(plot_key, str(redshift['0']), version['0'], legend['0'])
    
def zevol2(load_from_file=False):
    
    return loadFile(plot_key, str(redshift['1']), version['1'], legend['1'])

def zevol3(load_from_file=False):
    
    return loadFile(plot_key, str(redshift['2']), version['2'], legend['2'])
    
def zevol4(load_from_file=False):
    
    return loadFile(plot_key, str(redshift['3']), version['3'], legend['3'])

def zevol5(load_from_file=False):
    
    return loadFile(plot_key, str(redshift['4']), version['4'], legend['4'])

def zevol6(load_from_file=False):
    
    return loadFile(plot_key, str(redshift['5']), version['5'], legend['5'])