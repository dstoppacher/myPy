# Load packages
import arangeData as aD
import numpy as np
import os
system_info = os.getcwd()
start = system_info.find('anaconda')
mycomp = system_info[:start]

import time
ts = time.time()
import datetime 
date_time_stamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H:%M:%S')

def calc_physical_distance(velocity,z,H_0):
    
    #use Hubble law!
    #   c*z = H_0/D
    
    D = H_0 / (velocity * z)
      
    return D



def check_datatype(value):
    #print 'here check_datatype:', value
    try:
        return int(value)
    except ValueError:
        try: 
            return float(value)
        except ValueError:
            return value

def check_filter_band(filter_band):
    try:
        return float(filter_band)
    except ValueError:
            return get_filter_wavelenght(filter_band)

def check_magAB_band(magAB):
    try:
        return float(magAB)
    except ValueError:
            return get_magAB_flux0(magAB)

def conv_filtermag_to_ABmag(data,
                         band):
    
    """
    Purpose: Convert apparent Magnitude to AB-System
    last update: 17.02.2016
    
    input:
    --------
    mAB_array: numpy array(n,)
    band: string common name combination of Telescope and bandpass (e.g. SDSS_r)
    
    output:
    ---------
    mAB_array: numpy array(n,) apparent Magnitudes in AB-System
    
    
    sources used:
    --------------
    http://www.sdss.org/dr12/algorithms/fluxcal/#SDSStoAB
    """
    #print 'Convert filter magnitude to real ABmag! band -->', band
    if band.startswith('SDSSu'):
        data-=0.04
    elif band.startswith('SDSSg'):
        pass
    elif band.startswith('SDSSr'):
        pass
    elif band.startswith('SDSSi'):
        pass
    elif band.startswith('SDSSz'):
        data+=0.02
    else:
        print 'no filter found!'
        
    return data
       
def conv_and_find_Msun_filter(
        w_band='R',
        conv_to=False):
        
    #Read in solar abs magintudes in various filter bands
    myData = aD.ArangeData()
    data_array = myData.readAnyFormat(config=False, data_format='ASCII', data_shape='shaped', mypath=mycomp+'anaconda/pro/data/wavebands.txt', delim=',')

    print 'chosen w_band:', w_band, ' --> converted to:', conv_to   

    r = float(data_array[np.where(data_array[:,0]=='r'),1][0][0])
    g = float(data_array[np.where(data_array[:,0]=='g'),1][0][0])
    
    if conv_to == False:
        #find chosen filter band
        M_sun_w_band = float(data_array[np.where(data_array[:,0]==w_band),1][0][0])
           
    else:
        #convert wavebands:
        if w_band=='R':
            M_sun_w_band = r + 0.51 + 0.15 * (g-r)
        elif w_band=='B':
            M_sun_w_band = g + 0.51 + 0.60*(g-r) 
        elif w_band=='V': 	
            M_sun_w_band = g - 0.03 - 0.42*(g-r)
        elif w_band=='I': 	
            M_sun_w_band = float(data_array[np.where(data_array[:,0]=='i'),1][0][0]) - 0.75
        elif w_band=='J': 	
            M_sun_w_band = g + 0.39 + 0.37*(g-r)
        elif w_band=='H': 	
            M_sun_w_band = float(data_array[np.where(data_array[:,0]=='H'),1][0][0])	
        elif w_band=='F':    
            M_sun_w_band = r - 0.25 + 0.17 * (g-r)                                            
        elif w_band=='u':  	
            M_sun_w_band = float(data_array[np.where(data_array[:,0]=='u'),1][0][0])
        elif w_band=='g':	
            M_sun_w_band = (float(data_array[np.where(data_array[:,0]=='V'),1][0][0])	+ 0.03 + 0.42*r) / 0.42
        elif w_band=='r': 	
            M_sun_w_band = (-(float(data_array[np.where(data_array[:,0]=='V'),1][0][0])) - 0.03 + 0.58*g ) / 0.42
        elif w_band=='i': 	
            M_sun_w_band = float(data_array[np.where(data_array[:,0]=='I'),1][0][0]) + 0.75
        elif w_band=='z': 	
            M_sun_w_band = float(data_array[np.where(data_array[:,0]=='z'),1][0][0])
        
    return M_sun_w_band

def conv_absMag_to_appMag(absMag,
                          lum_dist):
    """Converts absolute magnitudes to apparent magnitudes of a certain luminosity distance [pc]
        
        absMag:     numpy array (n,) [-]
        lum_disk:   numpy array (n,) [Mpc]
        """    
                          
    return absMag + 5*np.log10(lum_dist*1e6/10.0)                       

def conv_appMag_to_absMag(appMag,
                          lum_dist):
    """Converts apparent magnitudes to absolut magnitudes of a certain luminosity distance [pc]
        
        appMag:     numpy array (n,) [-]
        lum_disk:   numpy array (n,) [Mpc]
        """   
                          
    return appMag - 5 *np.log10(lum_dist*1e6/10.0)

def conv_Watt_to_erg(data):
  
    """Converts from SI sytem [Watt] to cgs system [erg]"""
    
    return data*1e7





def conv_fluxdensity_to_appMag(fluxdensity,
                               filter_info=False):     
    
    """Convert flux density to apparent magnitudes"""    

    flux_ref_band = get_filter_flux_m0()
    
    #appMag = appMag_ref - 2.5 log10 (fluxdensity / fluxdensity_ref)
    #appMag_ref = 0 for the fluxdensity in a certain band
    print 'filter:', filter_info, 'f_m0:', flux_ref_band[filter_info+'_flux_m0']
    print 'flux[0:25]'
    print fluxdensity[0:25]
    
    return - 2.5*np.log10(fluxdensity / flux_ref_band[filter_info+'_flux_m0'])

     
def conv_fluxdensity_to_ABmag(flux_density,
                              band):
    
    """Converts flux density [erg s-1 cm-2 Hz-1] to apparent magnitudes in the AB System"""
                 
    #mAB = -2.5 log10(fluxdensity) - 48.60
    #or mAB = -2.5 log10(fluxdensity / 3631 Jy)
    #log10 ( f[erg s-1 cm-2 Hz-1] / 3631 Jy[10-23 erg s-1 cm-2 Hz-1] ) --> log10(f/3631Jy)
                       
    return - 2.5*np.log10(flux_density/3631e-23)
                              

def conv_fluxdensity_to_flux(flux_density,
                             lambda_filter):

    print 'lambda_filter:', lambda_filter
    wavelenght = check_filter_band(lambda_filter)            
                                 
    return flux_density*(wavelenght[lambda_filter]/2.998e8)                             

def conv_fluxdensity_to_maggies(flux_density,
                                flux_0):

    """Converts the flux [erg s-1 cm-2 Hz-1] to maggies applying mu = flux/flux_0
    
    flux:     numpy array (n,) [erg s-1 cm-2 Hz-1]
    flux_0:   zero-point of the flux in the assumed system (e.g. AB --> 3161 Jy [10-23 erg s-1 cm-2 Hz-1])
    """
                                                           
    return flux_density/flux_0 

def conv_lum_to_appMag(lum,
                       lum_dist,
                       filter_info):     
    
    """Convert luminosty to absolute magnitudes"""

    absMag_sun_band = get_filter_absMag_Sun()

    #calculate absolute Magnitudes
    solar_lum_ergs = conv_Watt_to_erg(3.845e26)    #Solar luminosity [erg s-1]
    print 'sol lum:', solar_lum_ergs, 'M_sun:', absMag_sun_band[filter_info]
    
#    i=0
#    while i<lum.size:
#        if i<27:  print lum[i], solar_lum_ergs
#        try:
#            np.log10(lum[i] / solar_lum_ergs)
#        except:
#            print 'Error:', lum[i], solar_lum_ergs
#        i+=1
        
      
    return absMag_sun_band[filter_info] - 2.5*np.log10( (lum / solar_lum_ergs) * (10.0 / lum_dist*1e6)**2 )

def conv_lum_to_flux(lum,
                     distance):

    """Converts the luminosity [W s-1] or [erg s-1] to flux [W s-1 cm-2] or [erg s-1 cm-2]
    
    lum:        numpy array (n,) [erg s-1]
    distance:   numpy array (n,) [Mpc]
    """
                                                           
    return lum/(4.0*np.pi*(conv_Mpc_to_cm(distance))**2) 

def conv_maggies_to_mag(mu):

    """Converts the maggies mu to mag
    
    mu:     numpy array (n,) [-])
    """
                                                           
    return -2.5*np.log10(mu)

def conv_Mpc_to_cm(value):

    """Convert Mpc to cm

    distance [Mpc]  --> distance [cm]
    
    1 Mpc = 1e6 pc
    1 pc = 3.0857e16 m
    1 m = 100 cm
    
    distance * 1e6 * 3.0857e16 * 100
    """
    
    return value*1e6*3.0857e16*100 

def conv_lumdensity_to_fluxdensity(lum,
                                   filter_info):

    """Converts the luminosity density [W Hz-1] or [erg s-1 Hz-1] to flux density [W s-1 cm-2] or [erg s-1 cm-2]

    lum:             numpy array (n,) [erg s-1]
    lambda_filter:   float [m]
    """

    lambda_filter = get_filter_wavelenght(filter_info)     
    
    print 'lambda filter [m]:', lambda_filter[filter_info]
    print 'lum[0:25]'
    print lum[0:25]
    print 'lum[0:25]/(lambda_filter*100)**2'
    print lum[0:25]/(lambda_filter[filter_info]*100)**2
                                                    
    return lum/(lambda_filter[filter_info]*100)**2

def conv_lum_to_solarflux(lum,
                          lambda_filter):

    """converts luminosities [Watt]==[Joule/s] to flux density [erg s-1 cm-2] using the wavelenght in the corrisponding filter band"""
            
    #calculate absolute Magnitudes
    solar_lum_erg = conv_Watt_to_erg(3.845e26)              #[Watt s-1]-->[erg s-1]
    d_sun = 1.49e8*1e5                                       #[cm]
    flux_density_sun = solar_lum_erg / (4.*np.pi*(d_sun**2)) #[erg s-1 cm-2]
    waveband = lambda_filter/2.998e8                         #[lambda(r-band) m / light speed ms-1]--> [s]
    
    return lum * flux_density_sun * waveband  

def conv_mag_to_lum(mAB,                      
                    lum_dist,
                    filter_info,                      
                    lambda_filter=0.0,
                    system_is_AB=False):

    lum = []
                        
    return lum
                            
def conv_flux_to_asinhMag(flux,
                         band):
    
    """Converts Luminosities to apparent Magintudes
        
    according to SDSS Early Data Relase, Stoughton+02, Astronomical Journal 123:485-548, Table 21 'asinh Magnitude Softening Parameters'
    band    b    0-Flux     mag [m(f/f0=0)]	 m(f/f0=10b)
    --------------------------------------------------------
    u	  1.4e-10	  24.63	                     22.12
    g       0.9e-10	  25.11	                     22.60
    r	  1.2e-10	  24.80                        22.29
    i	  1.8e-10	  24.36	                     21.85
    z	  7.4e-10	  22.83   	                20.32
    
    equation used to calculate mAB --> Eq(6) in Paper Stoughton+02
    """
#    print 'here: conv_to_asinhMag'
#    print '++++++++++++++++++++++++++++++++'
    
    myData = aD.ArangeData()
    data_array = myData.readAnyFormat(config=False, data_format='ASCII', data_shape='shaped', mypath=mycomp+'anaconda/pro/data/reference/asinhMag_parameters.txt', delim='\t')

    para_map = {}

    i=0
    while i<data_array[:,0].size:
        para_map[data_array[i,0]+'_b'] = float(data_array[i,1])
        para_map[data_array[i,0]+'_m0_1b'] = float(data_array[i,2])
        para_map[data_array[i,0]+'_m0_10b'] = float(data_array[i,3])
        para_map[data_array[i,0]+'_lambda_filter'] = float(data_array[i,4])
        para_map[data_array[i,0]+'_limiting_mag'] = float(data_array[i,5])
        i+=1  

    flux0 = get_filter_flux_m0()

    print 'band:', band, 'm0:', para_map[band+'_m0_1b'], flux0[band+'_flux_m0']
    print 'flux[0:25]'
    print flux[0:25]

    
    return - 2.5/np.log(10) * (np.arcsinh( (flux/flux0[band+'_flux_m0']) / (2*para_map[band+'_b']) ) + np.log(para_map[band+'_b']))
    
def choose_random_sample(data,
                         number_of_objects):

    """Chooses a random sample of objects from a coluum in an array by masking randomly selected row-ids (e.g. row nr=1,4,5,9,11 etc.) to the data"""                          
    import random 
    mask = np.asarray([random.sample(range(0,data.size,1),number_of_objects)])
   
    return data[mask[:][0]]

def convert_units_Galacticus(catname,
                              data,
                              id_col_array,
                              redshift,
                              telescope_name=' ',
                              apply_k_corr_app=False,
                              apply_k_corr_abs=False,
                              apply_z_boost=False,
                              hubble_par=1,
                              bands='',
                              unit_code='DS',
                              use_kcorrect=None,
                              which_kcorrect='approx',
                              cosmology='Planck'):

    print 'CONVERT UNITS GALACTICUS\n###################################################################################\n'
      
    i=0
    conv_name=''
    while i<int(id_col_array['nr_entries']):
        try:
            name=id_col_array['name'+str(i)]
            if np.sum(data[name])==0.0 and name.find('cut')==-1 and name!='Z':     
                def caseSwitcher(conv_name):
    
                    def L2MAB():
                        #print 'name:', name, 'index:', index
                        data[name] = -2.5*np.log10(data[index])
                        #print name, '--> Lum to abolute mag + SDSS/GALEX filter band correction!\n--------------------\n'
    
                    def L2mAB():
                        data[name] = -2.5*np.log10(data[index]) + 5*np.log10(lum_dist*1e6/10.0)+(-2.5*np.log10(1+redshift))    
                        #print name, '--> Lum to apparent mag + SDSS/GALEX filter band correction!\n--------------------\n'
                    
                    choose = {
                        'L_SDSS_2MAB': L2MAB,
                        'L_SDSS_2mAB': L2mAB,
                        'L_Galex_2MAB': L2MAB,
                        'L_Galex_2mAB': L2mAB 
                        }
                        
                    func = choose.get(conv_name)
                    return func()
                
                filters={'u': 'SDSS', 'g': 'SDSS', 'r': 'SDSS', 'i': 'SDSS', 'z': 'SDSS', 'B': 'RGO', 'NUV': 'Galex', 'FUV': 'Galex'}
                
                filter_name = id_col_array['name'+str(i)][len(id_col_array['name'+str(i)])-id_col_array['name'+str(i)][::-1].find('_')::]  
                if (id_col_array['name'+str(i)].startswith('mA') or id_col_array['name'+str(i)].startswith('MA') or id_col_array['name'+str(i)].startswith('mag') or id_col_array['name'+str(i)].startswith('Mag')) and catname.startswith('Gal'):
                    conv_name = 'L_'+filters[filter_name]+'_'
                    #print 'filters[]:', filters[filter_name], id_col_array['name'+str(i)], conv_name
                    lum_col_name=conv_name+id_col_array['name'+str(i)][4::]
                    #print lum_col_name, conv_name, filter_name
    
                else:
                    conv_name = 'MAB_'
                    lum_col_name=conv_name+id_col_array['name'+str(i)][4::]
                    
                c=0
                while c<int(id_col_array['nr_entries']):
                   
                    if id_col_array['name'+str(c)]==lum_col_name:                 
                        index=id_col_array['name'+str(c)]
                        #print 'c:', c, id_col_array['name'+str(c)], 'lum_col_name:', lum_col_name, 'INDEX:', index
                    c+=1
                  
                conv_name+='2'+id_col_array['name'+str(i)][0:3]
                #print 'conv_name:', conv_name
                #use the calculated redshift z_box of every galaxy to calculate the luminosity distance --> use luminosity distance to calculate magnitude  

                from cosmolopy import cparam, cd
                fidcosmo = cparam.Planck(flat=True, extras=True)
                #print 'fidcosmo:', fidcosmo
          
                if redshift==0.0:
                    lum_dist=cd.luminosity_distance(redshift, **fidcosmo)
                    lum_dist=1e-5 #[Mpc]             
                else:
                    lum_dist=cd.luminosity_distance(redshift, **fidcosmo)
                    #print 'lum_dist:', lum_dist, '[Mpc] z:', redshift
                        
                if apply_z_boost=='True':
                    z_boost=1+redshift
                else:
                    z_boost=1
                
                print 'i:', i, 'name:', id_col_array['name'+str(i)], 'mag_col_id:', id_col_array['col_id'+str(i)], 'filtername:', filter_name,\
                      'lum_col_name:', lum_col_name, 'lum_dist:', "{0:.2f}".format(lum_dist), 'z_boost:', "{0:.2f}".format(z_boost), 'unit_code:', unit_code
                               
                caseSwitcher(conv_name)
            
        except:
            pass
            #print name, 'is not able to convert ...\n'  
    
        i+=1
    
    return data
    
def correct_units(catname,
                  data,
                  id_col_array,
                  hubble_par,
                  scale_factor):

    """Corrects all catalogues data columns to one unified unit set and converts magnitudes to the AB system
    
    unified set: distance [Mpc], SFR [yr], no little-h in e.g. masses, position etc.
    corrects magnitudes to the AB system
    
    function name: correct_units
       
    catname:      name of the input catalogue, needed to identify the units used in the catalogue
    data:         numpy array where the catalogues data is storred
    id_col_array: maps which address to each column in 'data' the unit and the correction (if correction=1 --> no correction applied)
    hubble_par:   little-h applied to the whole 'data' array (if no little-h in the data hubble_par=1)
    """
    print ' '
    print '###################################################################################'
    print 'CORRECT UNITS'
    print ' '
    #print catname, 'correct general "little-h"!    hubble parameter:', hubble_par

    i=0
    while i<len(data[0]): 
        #print 'i:', i, id_col_array['corr_type'+str(i)], id_col_array['name'+str(i)], hubble_par
        if id_col_array['corr_type'+str(i)]!='id_num' and id_col_array['corr_type'+str(i)]!='int_num' and id_col_array['corr_type'+str(i)]!='no_corr':
            data[id_col_array['name'+str(i)]] = data[id_col_array['name'+str(i)]]/float(hubble_par)
        i+=1
                       
    i=0

    while i<len(id_col_array):

        try:
            #print 'i:', i, 'name:', id_col_array['name'+str(i)], 'corr_type:', id_col_array['corr_type'+str(i)], 'unit_corr:', id_col_array['unit_corr'+str(i)] 
            if id_col_array['corr_type'+str(i)]!='no_corr' and id_col_array['corr_type'+str(i)]!='int_num' and id_col_array['corr_type'+str(i)]!='id_num':        
                print 'i:', i, 'corr_type:', id_col_array['corr_type'+str(i)], 'name:', id_col_array['name'+str(i)], 'unit_corr:', id_col_array['unit_corr'+str(i)], 'col_id:', id_col_array['col_id'+str(i)]

                if id_col_array['name'+str(i)].find('pos')!=-1:

                    if id_col_array['corr_type'+str(i)]=='COMV' and id_col_array['unit_corr'+str(i)]=='True':
                        #print 'HERE comv corr: --> True, scale factor:', scale_factor
                        data[id_col_array['name'+str(i)]]/=scale_factor
                        
                    elif id_col_array['unit_corr'+str(i)]!='False':
                        data[id_col_array['name'+str(i)]]/=id_col_array['unit_corr'+str(i)]
                  
                    #print 'min:', min(data[:,id_col_array['col_id'+str(i)]]), 'max:', max(data[:,id_col_array['col_id'+str(i)]])
                elif id_col_array['unit_corr'+str(i)]=='power10':
                    #print 'power10'
                    data[id_col_array['name'+str(i)]]=10**data[id_col_array['name'+str(i)]]
                elif id_col_array['corr_type'+str(i)]=='lum':
                    #print 'here! lum:', id_col_array['unit_corr'+str(i)] 
                    data[id_col_array['name'+str(i)]]*=id_col_array['unit_corr'+str(i)] 
                elif id_col_array['corr_type'+str(i)]=='MAB':
                    pass
                    #print 'MAB-corr ... nothing to do!'         
                else:
                    #print 'default: /'
                    data[id_col_array['name'+str(i)]]/=id_col_array['unit_corr'+str(i)]
        except:
                #print 'no more entries to correct ...'
                break
        i+=1
        
    return data


    

def convert_units(catname,
                  data,
                  id_col_array,
                  redshift,
                  telescope_name=' ',
                  apply_k_corr_app=False,
                  apply_k_corr_abs=False,
                  apply_z_boost=False,
                  hubble_par=1,
                  bands='',
                  unit_code='DS',
                  use_kcorrect=None,
                  cosmology='Planck'):

    def caseSwitcher_kcorr(kcorrect):
    
        choose = {
            'approx': kcorrect_approx,
            'Blanton': kcorrect_Blanton,
            'kcorr_mAB': kcorrect_ABmag
            }
            
        func = choose.get(kcorrect)
        return func()
    
    def kcorrect_bands_mAB(k_corr_array):

        print 'mAB_obs=MAB_rest+Kcorr'
        for band in bands:
            print 'band:', band,
            #print 'before:', data[band][0:10]
            print 'kcorrection: median =', "{0:.2f}".format(np.median(k_corr_array[band])), '16th/84th', "{0:.2f}".format(np.percentile(k_corr_array[band], 16)), '/', "{0:.2f}".format(np.percentile(k_corr_array[band], 84))
  
            data[band]+=k_corr_array[band]

    def kcorrect_bands_MAB(k_corr_array):
 
        print 'MAB_rest=mAB_app-Kcorr'
        for band in bands:
            print 'band:', band,
            #print 'before:', data[band][0:10]
            print 'kcorrection: median =', "{0:.2f}".format(np.median(k_corr_array[band])), '16th/84th', "{0:.2f}".format(np.percentile(k_corr_array[band], 16)), '/', "{0:.2f}".format(np.percentile(k_corr_array[band], 84))

            data[band]-=k_corr_array[band]
        
    
    print 'CONVERT UNITS\n###################################################################################\n'
  
    from cosmolopy import cparam, cd
    fidcosmo = cparam.Planck(flat=True, extras=True)
    print fidcosmo   
    i=0
    bands=[]
    conv_name=''
    while i<int(id_col_array['nr_entries']):
        try:
            name=id_col_array['name'+str(i)]
            if np.sum(data[name])==0.0 and name.find('cut')==-1 and name!='Z':           
            #if name.find('MAB_dA_total')!=-1 and name.find('cut')==-1 and name!='Z':
                print '\n--------------------\n'
                
                def L2mag():
                    print 'Lum to apparent Mag!\n--------------------\n'
                    
                    flux_density = conv_lum_to_flux(data[index], lum_dist)
                    data[name] = conv_fluxdensity_to_appMag(flux_density*z_boost, filter_name)
                    
                        
                def MAB2mAB():
                    print 'MAB to mAB! Via DM + SDSS filter band correction!'
                    
                    print 'name:', name, 'hubble_par=', hubble_par, 'lum_dist=', lum_dist, '[Mpc] z=', redshift #, 'Mpc3 --> correct units to h-1Mpc3: lum_dist/h=', lum_dist/hubble_par
                    data[name]= data[index] + 5*np.log10(lum_dist*1e6/10.0) - 5*np.log10(hubble_par) +(-2.5*np.log10(1+redshift))
                    
                    if redshift!=0.0 and apply_k_corr_app!='False':
                        #calculate k-corrections
                        k_corr_array=caseSwitcher_kcorr(apply_k_corr_app)
                        
                        #apply k-correction on magnitudes
                        kcorrect_bands_mAB(k_corr_array)
    
                def L2mAsinh():
                    print 'Lum to asinhMag!\n--------------------\n'
                    
                    flux_density = conv_lum_to_flux(data[index], lum_dist)
                    data[name] = conv_flux_to_asinhMag(flux_density*z_boost, filter_name)
                    data[name] = conv_filtermag_to_ABmag(data[name], telescope_name+filter_name)
                                        
    
                def L2ABmag():
    
                    print 'Lum to apparent mag + SDSS filter band correction!\n--------------------\n'
    
                    flux_density = conv_lum_to_flux(data[index], lum_dist)
    
                    if unit_code=='MD':
                        data[name] = conv_fluxdensity_to_ABmag(flux_density*4.4659e20*z_boost, filter_name)
                    else:
                        data[name] = conv_fluxdensity_to_ABmag(flux_density*z_boost, filter_name)                        
                    data[name] = conv_filtermag_to_ABmag(data[name], telescope_name+filter_name)
                  
                    
                    if redshift!=0.0 and apply_k_corr_app!='False':
                        #calculate k-corrections
                        k_corr_array=caseSwitcher_kcorr(use_kcorrect)
                        
                        #apply k-correction on magnitudes
                        kcorrect_bands_mAB(k_corr_array)

                def app2abs():
    
                    print 'Calculate abs magnitude via DM! name:', name, 'lum_dist', lum_dist, '[Mpc], hubble_par:', hubble_par
                    data[name] = conv_appMag_to_absMag(data[name], lum_dist) -5*np.log10(hubble_par)

                    if redshift!=0.0 and apply_k_corr_abs!='False':
                        #calculate k-corrections
                        k_corr_array=caseSwitcher_kcorr(use_kcorrect)
                        
                        #apply k-correction on magnitudes
                        kcorrect_bands_MAB(k_corr_array)
  
    
                def caseSwitcher(conv_name):
                
                    choose = {
                        'L_SDSS_2mag': L2mag,
                        'L_SDSS_2mAB': L2ABmag,
                        'L_SDSS_2mAs': L2mAsinh,
                        'MAB_2mAB': MAB2mAB,
                        'MAB_2mag': MAB2mAB,
                        'mAB_2Mag': app2abs,
                        'mAB_2MAB': app2abs,
                        'L_SDSS_2MAB': L2ABmag,
                        'L_RGO_2Mag': L2ABmag,
                        'MAB_2kco': MAB2mAB
                        }
                        
                    func = choose.get(conv_name)
                    return func()
                
                filters={'u': 'SDSS', 'g': 'SDSS', 'r': 'SDSS', 'i': 'SDSS', 'z': 'SDSS', 'B': 'RGO'}
                
                filter_name = id_col_array['name'+str(i)][len(id_col_array['name'+str(i)])-1]    
                if (id_col_array['name'+str(i)].startswith('mA') or id_col_array['name'+str(i)].startswith('MA') or id_col_array['name'+str(i)].startswith('mag') or id_col_array['name'+str(i)].startswith('Mag')) and catname.startswith('Gal'):
                    conv_name = 'L_'+filters[filter_name]+'_'     
                    lum_col_name=conv_name+id_col_array['name'+str(i)][4::]                  
                    #print lum_col_name, conv_name
    
                else:
                    conv_name = 'MAB_'
                    lum_col_name=conv_name+id_col_array['name'+str(i)][4::]
                    
                c=0
                while c<int(id_col_array['nr_entries']):
                      
                    if id_col_array['name'+str(c)]==lum_col_name:                 
                        index=id_col_array['name'+str(c)]
                        #print 'c:', c, id_col_array['name'+str(c)], 'lum_col_name:', lum_col_name, 'INDEX:', index
                    c+=1
                  
                conv_name+='2'+id_col_array['name'+str(i)][0:3]
                #print 'conv_name:', conv_name
                #use the calculated redshift z_box of every galaxy to calculate the luminosity distance --> use luminosity distance to calculate magnitude
                                
                if redshift==0.0:
                    lum_dist=cd.luminosity_distance(redshift, **fidcosmo)
                    lum_dist=1e-5 #[Mpc]             
                else:
                    lum_dist=cd.luminosity_distance(redshift, **fidcosmo)
                    print 'lum_dist:', lum_dist, '[Mpc] z:', redshift
                        
                if apply_z_boost=='True':
                    z_boost=1+redshift
                else:
                    z_boost=1
                
                print 'i:', i, 'name:', id_col_array['name'+str(i)], 'mag_col_id:', id_col_array['col_id'+str(i)], 'index:', index, 'filtername:', filter_name,\
                      'lum_col_name:', lum_col_name, 'lum_dist:', "{0:.2f}".format(lum_dist), 'z_boost:', "{0:.2f}".format(z_boost), 'unit_code:', unit_code
                               
                caseSwitcher(conv_name)
    
                bands.extend([name])
                
        except:
            print name, 'is not able to convert ...'


        i+=1
                
    return data

def kcorrect_ABmag(data,
                  redshift):
    
    """Function applies the K-correction to all apparent magnitudes in the data array. The formular for the K-correction is taken from
    Blanton&Roweis (2007) eq. (8), see also Oke&Sandage (1968)

    magnitude relation:
    
        m_B = M_Q + DM[h-1pc] + Kcorr - 5*log10(h)
    
    K-correction formular to convert from on bandpass at a certain redshift to another at a certain redshift can be assumed with
    
        Kcorr = -2.5*log10(1+z)

    data:       magnitudes where the DM is already applied
    redshift:   desired bandpass to which the data should be K-correted
    bands:      list of bands where the K-correction should be applied
    hubble_par: little h of the used cosmology

    output: apparent AB magnitude in the bandpass of z, K-corrected with a simple formular from Blanton&Roweis (2007)
              
    """

    print '###################################################################################'
    print 'K-correction for AB magnitudes --> to desired bassband, NOT to REST FRAME!'
         
    return -2.5*np.log10(1.0+redshift)
  
    
def kcorrect_approx(data,
                  redshift,
                  bands):


    """Function applies the K-correction to all apparent magnitudes in the data array. The formular for the K-correction is taken from
    Blanton&Roweis (2007) eq. (1), see also Oke&Sandage (1968)
    
    magnitude relation:
    
        m_B = M_Q + DM[h-1pc] + Kcorr - 5*log10(h)
        
    m_B:        apparent magnitude in the desired band pass
    M_Q:        absolute magnitude in the rest frame
    DM:         distance modulus
    Kcorr:     pyhton package from Chiligarian+2010 MRAS 405. 1409-1420 (NOTE ONLY FOR REDSHIFT z<0.5)
    hubble_par: little h of the used cosmology
        
    This formular is applied to the data in order to get the apparent/absolute magnitudes!       
    """
    print ' '
    print '###################################################################################'
    print 'K-CORRECTION for MAGNITUDES approximation '
    print ' '

    k_corr_array = data[bands].copy()
    for band in bands:
        k_corr_array[band]=0.0
    
    i=0
    for name in bands:

        filter_name = name[len(name)-1::]
        print 'name:', name, 'filter:', filter_name,
                                                           
        #find second filter to do k-corr!
        colour_value = filtermap_kcorr(str(filter_name))
        filter_name1 = str(colour_value[0:0+colour_value.find(' - ')])
        filter_name2 = str(colour_value[colour_value.find(' - ')+3::])

        if filter_name1==filter_name:
            filter1 = name
            filter2 = name[0:len(name)-1]+str(filter_name2)
        else:
            filter2 = name             
            filter1 = name[0:len(name)-1]+str(filter_name1)
        
        #if redshift>0.50: redshift=0.5               
        print 'i:', i, 'name:', name, 'FILTERNAME to correct:', filter_name, 'redshift:', redshift
        print 'colour_value:', colour_value,'filter_name1:', filter_name1, 'filter_name2:', filter_name2, 'filter1:', filter1, 'index2:', filter2         
    
        k_corr_array[name] = calc_kcor(filter_name, redshift, colour_value, data[filter1]-data[filter2])
           
        print 'k_corr_array:', k_corr_array[name][0:6], '\n-----------\n'

    return k_corr_array
  

def kcorrect_Blanton(redshift):

    """ python_kcorrect from https://pypi.python.org/pypi/kcorrect_python/2013.10.12 via pip install (04/12/2017)
    
        available functions
        -------------------
        
        The following functions are currently available for this version:
        
        o :func:`load_templates`
        o :func:`load_filters`
        o :func:`fit_coeffs_from_file`
        o :func:`fit_coeffs`
        o :func:`reconstruct_maggies`
        o :func:`reconstruct_maggies_from_files`
        o :func:`fit_photoz`
        o :func:`fit_photoz_from_file`
        
        examples
        --------
        
        The example below uses the data shipped with kcorrect.v4_2.
        You can use the module as follow::
        
        >>> import kcorrect, numpy
        >>> kcorrect.load_templates()
        >>> kcorrect.load_filters()
        >>> a=[0.03077382, 1.144068e-08, 5.262234e-08, 8.210213e-08, 8.744532e-08, 1.017738e-07, 6.216309e+16, 3.454767e+17, 1.827409e+17, 1.080889e+16, 3163927000000000.0]
        >>> c = kcorrect.fit_coeffs(a)
        >>> c
        array([ 3.07738204e-02, 2.02254747e-14, 1.49129165e-35,
        2.15513887e-06, 6.94462278e-06, 1.78061924e-13], dtype=float32)
        >>> m = kcorrect.reconstruct_maggies(c)
        >>> m
        array([ 3.07738204e-02, 1.44426586e-08, 5.28384980e-08,
        8.09117182e-08, 9.51680121e-08, 1.10408600e-07], dtype=float32)
        
        The example above successively loads the module,
        loads the default templates, *vmatrix.default.dat*
        and *lambda.default.dat*, loads the default filter,
        *sdss_filters.dat*, then computes the coeffs and
        reconstructs maggies.
        
        To compute the reconstructed maggies at rest-frame with bandpasses
        shifted by 0.1, you need first reload the filters with the given
        band_shift, then compute the coeffs and the maggies::
        
        >>> kcorrect.load_filters(band_shift=0.1)
        >>> m0 = kcorrect.reconstruct_maggies(c, redshift=0.)
        
        If the redshifs, maggies and maggies_invvar are stored
        in a file like *sample.dat* found in the *test* directory
        of kcorrect package, you can use :func:`fit_coeffs_from_file`
        and :func:`reconstruct_maggies_from_files` to perform the
        computation::
        
        >>> kcorrect.fit_coeffs_from_file('some_file.dat', outfile='output_coeffs.dat')
        >>> kcorrect.reconstruct_maggies_from_files('output_coeffs.dat', outfile='computed_maggies.dat')
        
        these produce 2 files *output_coeffs.dat* and *computed_maggies.dat*
        
        To use different templates, you load them as follow::
        
        >>> kcorrect.load_templates(v='vmatrix.goods.dat',l='lambda.goods.dat')
        
        If templates and filters are not loaded before calling the other
        functions, error is raised::
        
        >>> kcorrect.fit_coeffs(range(11))
        Traceback (most recent call last):
        File "<stdin>", line 1, in <module>
        File "kcorrect.py", line 37, in fit_coeffs
        return _kcorrect.fit_coeffs(c)
        _kcorrect.error: no filters loaded.
        
        :func:`fit_photoz` and :func:`fit_photoz_from_file` can be used
        as follow, after loading the appropriate templates and filter::
        
        >>> p = kcorrect.fit_photoz(a[1:])
        >>> p
        array([ 1.41886109e-02, 5.18920551e-09, 6.65258128e-36,
        2.18073205e-06, 5.97664302e-06, 4.88666385e-14], dtype=float32)
        
        if the data are from a file *photoz.dat*::
        
        >>> fit_photoz('photoz.dat', outfile='photoz.out')
        
        which produces the result to the output file *photoz.out*
        """
    
    import kcorrect
    kcorrect.load_templates()
    kcorrect.load_filters()
    a=[0.03077382, 1.144068e-08, 5.262234e-08, 8.210213e-08, 8.744532e-08, 1.017738e-07, 6.216309e+16, 3.454767e+17, 1.827409e+17, 1.080889e+16, 3163927000000000.0]
    c = kcorrect.fit_coeffs(a)
    print 'output c:', c

    m = kcorrect.reconstruct_maggies(c)
    print 'output reconstr maggies: ', m

    print 'load filters, band shift to :', redshift
    kcorrect.load_filters(band_shift=redshift)
 
    return kcorrect.reconstruct_maggies(c, redshift=0.0)

def filtermap_kcorr(filtername):
    
    filtermap = {}
    
    filtermap['u'] = {'u'+str(0): 'u - r', 'u'+str(1): 'u - i', 'u'+str(2): 'u - z'}
    filtermap['g'] = {'g'+str(0): 'g - r', 'g'+str(1): 'g - i', 'g'+str(2): 'g - z'}
    filtermap['r'] = {'r'+str(0): 'g - r', 'r'+str(1): 'u - r'}
    filtermap['i'] = {'i'+str(0): 'g - i', 'i'+str(1): 'u - i'}
    filtermap['z'] = {'z'+str(0): 'g - z', 'z'+str(1): 'u - z', 'z'+str(2): 'i - z'}
    
    return filtermap[filtername][filtername+str(0)] 

def convert_coordinates(data):
    
    dist            =cd.comoving_distance(data['Z'], **fidcosmo)
    conv_deg_to_rad =np.pi/180.0
    
    data['x_pos']=dist*np.sin(data['DEC']*conv_deg_to_rad)*np.sin(data['RA']*conv_deg_to_rad)
    data['y_pos']=dist*np.sin(data['DEC']*conv_deg_to_rad)*np.cos(data['RA']*conv_deg_to_rad)
    data['z_pos']=dist*np.cos(data['DEC']*conv_deg_to_rad)
    
    return data

def processInput(a,data,centrals,pos):
    data[pos][np.where((data['haloid']==centrals['haloid'][a]) & (data['orphan']==2))]=centrals[pos][a]
    #return np.where((data['haloid']==centrals['haloid'][a]) & (data['orphan']==2)), centrals[pos][a]
    
def correctOrphan2CentralPos_not_working(data,
                                       prop=False, 
                                       name=False):
                                      
    print 'here! Galacticuts set orphan positions ...'

    data_centrals=data[np.where(data['haloid']==data['hostid'])]
    data_type0=data[np.where(data['orphan']==0)]
    data_sats=data[np.where(data['orphan']==1)]
    data_orphans=data[np.where(data['orphan']==2)]
    
    test = np.in1d(data_centrals['haloid'], data_orphans['haloid'])
    
    print 'ngal haloid==hostid', len(data_centrals), 'ngal type 0:', len(data_type0), 'sats:', len(data_sats), 'orphans:', len(data_orphans), 'orphans in centrals:', len(np.where(test==True))           
    centrals=data_centrals[np.where(test==True)]
    print 'centrals with orphan satellite:', len(centrals)                           

    from joblib import Parallel, delayed
    import multiprocessing
    print data.flags
    num_cores = multiprocessing.cpu_count()   
    Parallel(n_jobs=8)(delayed(processInput)(i,data,centrals,'x_pos') for i in range(5))
    #print result
    #data['x_pos'][[0]]=Parallel(n_jobs=num_cores)(delayed(processInput)(i,data,centrals) for i in range(len(centrals)))[1]
    print data[0:5]
    
    return data

def correctOrphan2CentralPos_slow(data,
                                   prop=False, 
                                   name=False):
    print 'here! Galacticuts set orphan positions ...'

    data_centrals=data[np.where(data['haloid']==data['hostid'])]
    data_type0=data[np.where(data['orphan']==0)]
    data_sats=data[np.where(data['orphan']==1)]
    data_orphans=data[np.where(data['orphan']==2)]
    
    test = np.in1d(data_centrals['haloid'], data_orphans['haloid'])
    
    print 'ngal haloid==hostid', len(data_centrals), 'ngal type 0:', len(data_type0), 'sats:', len(data_sats), 'orphans:', len(data_orphans), 'orphans in centrals:', len(np.where(test==True))           
    centrals=data_centrals[np.where(test==True)]
    print 'centrals with orphan satellite:', len(centrals)                           
    a=0
    while a<len(centrals):
        orphans=np.where((data['haloid']==centrals['haloid'][a]) & (data['orphan']==2))
        data['x_pos'][orphans]=centrals['x_pos'][a]        
        data['y_pos'][orphans]=centrals['y_pos'][a]       
        data['z_pos'][orphans]=centrals['z_pos'][a]                            
        a+=1

    return data    


def calc_colour_cut_parameter(data,
                              id_col_array,
                              input_name):
   

    def r_i():
        print 'r-i!'
        data[name+'cut_r_i'] = data[name+'r'] - data[name+'i']

    def g_r():
        print 'g_r!'
        data[name+'cut_g_r'] = data[name+'g'] - data[name+'r']

    def g_i():
        print 'g_i!'
        data[name+'cut_g_i'] = data[name+'g'] - data[name+'i']

    def dmesa():
        print 'dmesa!'
        data[name+'cut_dmesa'] = data[name+'r'] - data[name+'i'] - (data[name+'g'] - data[name+'r'])/8.0
      
    def cpar():
        print 'cpar!'
        data[name+'cut_cpar'] = 0.7*(data[name+'g'] - data[name+'r']) + 1.2*((data[name+'r'] - data[name+'i']) - 0.18)

    def cmesa():
        print 'cperp!'
        data[name+'cut_cperp'] =  abs((data[name+'r'] - data[name+'i']) - (data[name+'g'] - data[name+'r'])/4.0 - 0.18)
        
    def i_lt_dmesa():
        print 'i_lt_dmesa!'
        i_lt_dmesa = 19.86 + 1.6*(data[name+'cut_dmesa'] - 0.8) - data[name+'i']
    
        data[name+'cut_i_lt_dmesa'][np.where(i_lt_dmesa >= 0.0)[0][:]]=1.01      
        #print 'positive:', np.where(i_lt_dmesa > 0.0)[0]
        data[name+'cut_i_lt_dmesa'][np.where(i_lt_dmesa < 0.0)[0][:]]=0
        #print 'negative:', np.where(i_lt_dmesa < 0.0)[0]

    def i_lt_dmesa_sparse():
        print 'i_lt_dmesa!'
        i_lt_dmesa_sparse_min = 19.86 + 1.6*(data[name+'cut_dmesa'] - 0.8)
        i_lt_dmesa_sparse_max = 20.14 + 1.6*(data[name+'cut_dmesa'] - 0.8)

            
    def r_lt_cpar():
        print 'r_lt_cpar!'
        r_lt_cpar = 13.5 + data[:, id_col_array[name+'cut_cpar_col_id']]/0.3

    
    def caseSwitcher(key):
    
        choose = {
            'r_i': r_i,
            'g_r': g_r,
            'g_i': g_i,
            'dmesa': dmesa,
            'cmesa': cmesa,
            'cpar': cpar,           
            'i_lt_dmesa': i_lt_dmesa,
            'i_lt_dmesa_sparse': i_lt_dmesa_sparse,
            'r_lt_cpar': r_lt_cpar
            }
            
        func = choose.get(key)
        return func()

    print 'CALC COLOUR CUT PARAMETERS\n###################################################################################\n'

    name    = input_name[0:input_name.find('cut_')]
    key     = input_name[input_name.find('cut_')+4::]
    print 'name:', name, 'key:', key, '-->',
            
    caseSwitcher(key)
    '\n-----------\n'
                                  
    return data                              
    
def dim_expander(data,
                col,
                id_of_col_to_expand=0):

    """Expands the dimions of an (n,n)-dim numpy array to one with (n, n+my_expansion)-dim
          
    data:                n-dim numpy array
    nr_cols_exapnd_to:   totalnumber of coluumns the array should have after expanding. The number of coluumns I am expanding my array to!
    
    keywords:
    ---------
    id_of_col_to_expand: id of the coluumn which schould be taken to expand the array (default: first column==0)
    """

    if np.prod(data.shape)==np.sum(data.shape):
        data = np.expand_dims(data, axis=1)

    while np.prod(data.shape)<=col*data[:,0].size:
        data = np.concatenate((data, np.expand_dims(data[:,id_of_col_to_expand], axis=1)), axis=1)
        if np.prod(data.shape)>=col*data[:,0].size:
            break
            
    return data

def dim_expander_struct(data,
                        name,
                        name_new,
                        mydtype=None):

    """Expands the dimions of an (n,n)-dim numpy array to one with (n, n+my_expansion)-dim
 
    data:                n-dim numpy array
    nr_cols_exapnd_to:   totalnumber of coluumns the array should have after expanding. The number of coluumns I am expanding my array to!
    
    keywords:
    ---------
    id_of_col_to_expand: id of the coluumn which schould be taken to expand the array (default: first column==0)
    """

    print 'name:', name, 'name_new:', name_new    
    
    import numpy.lib.recfunctions as rcfuncs
    data = rcfuncs.append_fields([data], name_new, data[name], usemask=False, dtypes=mydtype)
    
    #print data
    #print np.info(data)
            
    return data

def redshift_to_expfactor(redshift):
    
    """Converts the redshift to the expansion factor"""
    
    return 1.0/(redshift + 1.0) 


def expfactor_to_redshift(exp_factor):
    
    """Converts the expansion factor to the redshift""" 
   
    return 1.0 / exp_factor - 1


def filter_data_before_analysis(data,
                                mycond_min,
                                mycond_max,
                                myselected_col,
                                myplot_key='',
                                mycol_orphan=False):

    myData = aD.ArangeData()
    #print data.shape
    #find max/min-values for calculate histo/binned array       
    cond_min, cond_max = find_min_and_max_values(data[myselected_col], mycond_min, mycond_max)
     
    #print 'check: min/max --> data:', min(data[myselected_col]), '/', max(data[myselected_col]), data.shape
    
    #print 'cond_min:', cond_min, 'cond_max', cond_max, 'sel_col:', myselected_col, 'myplot_key:', myplot_key, 'col_orphan_status:', mycol_orphan
    data = myData.selectData2Compute(data, 
                                      selected_col=myselected_col, 
                                      operator='>=', 
                                      condition=cond_min)
                                                                                                
    data = myData.selectData2Compute(data, 
                                      selected_col=myselected_col, 
                                      operator='<=', 
                                      condition=cond_max)     

    #print 'check: min/max --> data:', min(data[myselected_col]), '/', max(data[myselected_col]), data.shape
    return data
    

def filter_data_before_write2file(data,
                                  myconds_array,
                                  name_col_array,
                                  myheader,
                                  mylog):
    
    sel_col_list = [] 
    mydataformat = ''
    check_size   = True

    print 'start filter_data_before_write2file():'
    print '----------------------------------------'
    print ' '

    i=0
    count=0
    while i<myconds_array['nr_entries']:
                        
        name=myconds_array[str(i)+'_name']           
        sel_col=myconds_array[name+'_col_id']
        mycond_min = myconds_array[name+'_min']
        mycond_max = myconds_array[name+'_max']         
           
        print 'name:', name, 'sel_col:', sel_col, 'i:', i, 'real col id:', myconds_array[name+'_col_id'], 'data.shape:', data.shape, 'mycond_min/max:', mycond_min, '/',mycond_max

        if myconds_array[name+'_output_col_name']!='False':
            col_name=myconds_array[name+'_output_col_name']
        else:
            col_name=name
           
        if myconds_array[name+'_exclude']!='yes': 
            sel_col_list+=[name]
          
            myheader+='('+str(count+1)+') '+col_name+' ['+str(myconds_array[name+'_unit'])+'] '         
            mydataformat+=str(myconds_array[name+'_format'])
            log_name_des = '             name: '
            if myconds_array['nr_entries']>1 and i!=myconds_array['nr_entries']-1:  mydataformat+='  '
            count+=1
        else:
            log_name_des = 'Not in file! name: '

        if mycond_min=='min' and mycond_max=='max':
            print 'nothing to filter -->  mycond_min:', mycond_min, 'mycond_max:', mycond_max
            print '-------------------'
            print ' '
        elif mycond_min=='inf' or mycond_max=='inf':
            print 'data.shape before test isfinite:', data.shape
            mask=np.where(np.isfinite(data[name]))
            print 'mask:', mask
            data = data[mask[:][0]]
            print 'data.shape after test isfinite:', data.shape
        else:
            data = filter_data_before_analysis(data,
                                               mycond_min,
                                               mycond_max,
                                               name)

            print '-------------------'
            print ' '
            
        if str(myconds_array[name+'_format'])=='%i':                                                            
            max_data = int(max(data[name]))
            min_data = int(min(data[name]))
        else:
            if name.startswith('L_') or name.find('age')!=-1 or name.find('ang')!=-1 or name.startswith('satelliteMergeTim') or name.startswith('mbasi') or name.startswith('sfr') or name.startswith('ssfr') or name.startswith('r') or name.startswith('mstar') or name.startswith('mhalo') or name.startswith('mcold') or name.startswith('mhot') or name.startswith('mbh') or name.find('bh_acc')!=-1 or name.startswith('zgas') or name.startswith('zstar') or name.startswith('zhot') or name.startswith('Mz'):
                
                max_data = "{0:.2e}".format(max(data[name]))
                min_data = "{0:.2e}".format(min(data[name]))
            else:    
                max_data = "{0:.2f}".format(max(data[name]))
                min_data = "{0:.2f}".format(min(data[name]))
                                              
        mylog+=log_name_des+col_name.ljust(25)+'['+(str(myconds_array[name+'_unit'])+']').ljust(25)+'ngal: '+str(data[name].size).ljust(20)+'min/max: '+(str(min_data)+'/'+str(max_data)).ljust(35)+'cut: '+str(mycond_min)+'/'+str(mycond_max)+'\n'
                                        
        if data[name].size==0:
            check_size=False                
            break
 
        i+=1
        
    #print 'sel_col_list:', sel_col_list
    #print 'mydataformat:', mydataformat
    return data, myheader, mydataformat, mylog, check_size, sel_col_list


def find_min_and_max_values(data,
                        value_min,
                        value_max):
    
    #find max/min-values for calculate histo/binned array       
    if value_min=='min' and value_max=='max':
        #print 'min,max'
        value_min = min(data)
        value_max = max(data)         
    elif value_min=='min' and value_max!='max':
        #print 'min,float'
        value_min = min(data)
        value_max = float(value_max)
    elif value_min!='min' and value_max=='max':
        #print 'float,max'
        value_min = float(value_min)
        value_max = max(data)
    else:
        #print 'float,float'
        value_min = float(value_min)
        value_max = float(value_max)
    
    return value_min, value_max

def get_filter_absMag_Sun():
 
    """Returns the a map with the absolute Magnitudes of the sun in a certain bass band filter
    
    from: https://www.astro.umd.edu/~ssm/ASTR620/mags.html (03/01/2016)
    
    Filter	Msun 	Source
    U 	5.61 	B&M
    B 	5.48 	B&M
    V 	4.83 	B&M
    R 	4.42 	B&M
    I 	4.08 	B&M
    J 	3.64 	B&M
    H 	3.32 	B&M
    K 	3.28 	B&M
    K' 	3.27 	*
    Spitzer 		
    3.6mu 	3.24 	Oh
    4.5mu 	3.27 	Oh
    SDSS 		
    u 	6.55 	S&G
    g 	5.12 	S&G
    r 	4.68 	S&G
    i 	4.57 	S&G
    z 	4.60 	S&G

    B&M = Binney & Merrifield
    S&G = Sparke & Gallagher
    Oh = Oh et al (2008) AJ, 136, 2761
    *My (SSM) estimate for the K' filter used by 2MASS after long and painful hunting through the calibration literature.        
    """    
    return {'u': 6.55, 'g': 5.12, 'i': 4.57, 'r': 4.68, 'z': 4.60, 'return_unit': '-'} 

def get_filter_flux_m0():
 
    """Returns the a map with the absolute Magnitudes of the sun in a certain bass band filter

    from: https://www.astro.umd.edu/~ssm/ASTR620/mags.html (03/01/2016)
    
    Band 	lambda_c 	dlambda/lambda 	Flux at m=0 	Reference
    	um 		                          Jy 	
    U 	0.36 	0.15 	         1810 	          Bessel (1979)
    B 	0.44 	0.22 	         4260 	          Bessel (1979)
    V 	0.55 	0.16 	         3640 	          Bessel (1979)
    R 	0.64 	0.23 	         3080 	          Bessel (1979)
    I 	0.79 	0.19 	         2550 	          Bessel (1979)
    J 	1.26 	0.16 	         1600 	          Campins, Reike, & Lebovsky (1985)
    H 	1.60 	0.23 	         1080 	          Campins, Reike, & Lebovsky (1985)
    K 	2.22 	0.23 	         670 	          Campins, Reike, & Lebovsky (1985)
    g 	0.52 	0.14 	         3730 	          Schneider, Gunn, & Hoessel (1983)
    r 	0.67 	0.14 	         4490 	          Schneider, Gunn, & Hoessel (1983)
    i 	0.79 	0.16 	         4760 	          Schneider, Gunn, & Hoessel (1983)
    z 	0.91 	0.13 	         4810 	          Schneider, Gunn, & Hoessel (1983)
    
    Also useful are these identities:
    
    1 Jy = 10^-23 erg sec^-1 cm^-2 Hz^-1
    1 Jy = 1.51e7 photons sec^-1 m^-2 (dlambda/lambda)^-1
    
    Example: How many V-band photons are incident per second on an area of 1 m^2 at the top of the atmosphere from a V=23.90 star? From the table, the flux at V=0 is 3640 Jy; hence, at V=23.90 the flux is diminished by a factor 10^(-0.4*V)=2.75e-10, yielding a flux of 1.e-6 Jy. Since dlambda/lambda=0.16 in V, the flux per second on a 1 m^2 aperture is
    
    f=1.e-6 Jy * 1.51e7 * 0.16 = 2.42 photons sec^-1       
    """   
    
    return {'u_flux_m0': 3631e-23, 'g_flux_m0': 3730e-23, 'r_flux_m0': 4490e-23, 'i_flux_m0': 4760e-23, 'z_flux_m0': 4810e-23, 'description': 'flux at m=0 for filter band', 'return_unit': 'erg s-1 cm-2 Hz-1'}

def get_filter_wavelenght(filter_band):
 
    """#SDSS Early Data Relase, Stoughton+02, Astronomical Journal 123:485-548, Table 21 "asinh Magnitude Softening Parameters"
    #band	b	0-Flux mag [m(f/f0=0)]	m(f/f0=10b)	filter wavelenght [m]	limiting Magnitude (95% completeness)
    u	1.4e-10	24.63	22.12	3551e-10	22.0
    g	0.9e-10	25.11	22.60	4686e-10	22.2
    r	1.2e-10	24.80	22.29	6166e-10	22.2
    i	1.8e-10	24.36	21.85	7480e-10	21.3
    z	7.4e-10	22.83	20.32	8932e-10	20.5
    """
   
    return {'u': 3551e-10, 'g': 4686e-10, 'r': 6166e-10, 'i': 7480e-10, 'z': 8932e-10, 'description': 'wavelenght', 'return_unit': 'm'}

def get_magAB_flux0(mAB_flux0):

    """#SDSS Early Data Relase, Stoughton+02, Astronomical Journal 123:485-548, Table 21 "asinh Magnitude Softening Parameters"
    #band	b	0-Flux mag [m(f/f0=0)]	m(f/f0=10b)	filter wavelenght [m]	limiting Magnitude (95% completeness)
    u	1.4e-10	24.63	22.12	3551e-10	22.0
    g	0.9e-10	25.11	22.60	4686e-10	22.2
    r	1.2e-10	24.80	22.29	6166e-10	22.2
    i	1.8e-10	24.36	21.85	7480e-10	21.3
    z	7.4e-10	22.83	20.32	8932e-10	20.5
    """
   
    return {'u': 24.63, 'g': 25.11, 'r': 24.80,  'i': 24.36, 'z': 22.83, 'description': 'AB magntude in filter band with flux=0','return_unit': '-'}

def myMultiColDetector(myconfig_datafile):

        myData = aD.ArangeData()
        #print 'len myconfig datafile:', myconfig_datafile, len(myconfig_datafile)
        data_array = myData.readAnyFormat(config=False, mypath=myconfig_datafile, mydtype=np.str_, data_format='ASCII', data_shape='shaped', delim='= ')     
        myconfig_array = np.empty([data_array[:,0].size + 1000, 3], dtype='S'+str(len(myconfig_datafile)+2000))        
        
        
        i=0
        count=0
        while i<data_array[:,1].size:         
            
            length_str = len(data_array[i][1])
            #print 'length_str:', length_str
            count_my_char = data_array[i][1].count(',')
            if count_my_char!=0:
                #print 'count_my_cahr:', count_my_char#, 'bucket', bucket.shape
                
                start=0
                j=0
                while j<count_my_char:             
                    
                    #print 'j:', j, data_array[i,1][start:length_str], 'start:', start
                    multicol_test = data_array[i,1][start:length_str].find(',')
    
                    if multicol_test!=-1:
                        
                        #print 'i:', i, 'count:', count
                        #print 'mydata:', data_array[i,1][start:start+multicol_test], 'len:', length_str, 'count_my_char', count_my_char, 'test:', multicol_test
                        
                        myconfig_array[count,0] = str(data_array[i,1][start:start+multicol_test])
                        myconfig_array[count,1] = i
                        myconfig_array[count,2] = data_array[i,0]+str(j) 
                        
                        
                        #print 'bucket[j]:', myconfig_array[count,0]
                        
                        if j==count_my_char-1:
                            myconfig_array[count+1,0] = data_array[i,1][start+multicol_test+1:]
                            myconfig_array[count+1,1] = i
                            myconfig_array[count+1,2] = data_array[i,0]+str(j+1)
                            #print 'last bucket:', myconfig_array[count+1,0]
                            count+=1
                            
                        start = start + multicol_test+1
                        count+=1
                        

                    j+=1
            else:
                #print 'i:', i, 'count:', count
                #print 'mydata:', data_array[i,1]
                myconfig_array[count,0] = data_array[i,1]
                myconfig_array[count,1] = i
                myconfig_array[count,2] = data_array[i,0]
                count+=1
            
                       
            i+=1

        return myconfig_array[0:count,:]

def multicolTestAlgorithm(data,
                          delimiter=False):

    if delimiter==False:
        count_my_char = data.count(',')
        delimiter=','
    else:
        count_my_char = data.count(delimiter)

    multicol_result = np.empty([count_my_char+1], dtype='S'+str(len(data)+200))

    i=0
    count=0
    while i<data.size:         
        
        length_str = len(data)
        
        if count_my_char!=0:
            
            start=0
            j=0
            while j<count_my_char:             
                
                multicol_test = data[start:length_str].find(delimiter)

                if multicol_test!=-1:                  
                   
                    multicol_result[count] = str(data[start:start+multicol_test])
                    
                    if j==count_my_char-1:
                        multicol_result[count+1] = data[start+multicol_test+1:]

                        count+=1
                        
                    start = start + multicol_test+1
                    count+=1
                    

                j+=1
        else:

            multicol_result[count] = data
            count+=1
        i+=1


    return multicol_result




def calculate_NFW_con_with_fit(mhalo,
                               redshift,
                               cosmology='Planck',
                               overdens='vir'):

    def approx_NFW_con_z000_Planck():
        """approximates NFW concentration
        from Klypin+16 Eq. (24) for PLANCK cosmology, Mhalo[200c], Table 2/A1, all halos selected by mass, z=0.0
        Parameter: c0=7.40, gamma=0.120, M0=5.5e5/1e12h-1
        """
        M12=1e12/0.6778
        M0=5.5e5
        
        return  7.40*(mhalo/M12)**-0.120 * (1+(mhalo/M0/M12)**0.4)
    
    def approx_NFW_con_z050_Planck():
        """approximates NFW concentration
        from Klypin+16 Eq. (24) for PLANCK cosmology, Mhalo[200c], Table 2/A1, all halos selected by mass, z=0.5
        Parameter: c0=6.5, gamma=0.105, M0=1e4/1e12h-1
        """
        M12=1e12/0.6778
        M0=1e4
        
        return  6.5*(mhalo/M12)**-0.105 * (1+(mhalo/M0/M12)**0.4)
    
    def approx_NFW_con_z050_WMAP7():
        """approximates NFW concentration
        from Klypin+16 Eq. (24) for WMAP7 cosmology, Mhalo[200c], Table A5, all halos selected by mass, z=0.5
        Parameter: c0=5.25, gamma=0.105, M0=6e4/1e12h-1
        """
    
        M12=1e12/0.6778    
        M0=6e4
        
        return  5.25*(mhalo/M12)**-0.105 * (1+(mhalo/M0/M12)**0.4)
    
    def approx_NFW_con_z000_Planck_vir():
        """approximates NFW concentration
        from Klypin+16 Eq. (24) for WMAP7 cosmology, Mhalo[vir], Table A3, all halos selected by mass, z=0.0
        Parameter:
        """
        c0=9.75
        gamma=0.110
        M0=5e5
        M12=1e12/0.6778    
        
        return  c0*(mhalo/M12)**-gamma * (1+(mhalo/M0/M12)**0.4)
    
    def approx_NFW_con_z035_Planck_vir():
        """approximates NFW concentration
        from Klypin+16 Eq. (24) for WMAP7 cosmology, Mhalo[vir], Table A3, all halos selected by mass, z=0.35
        Parameter:
        """
        c0=7.25
        gamma=0.107
        M0=2.2e4
        M12=1e12/0.6778    
        
        return  c0*(mhalo/M12)**-gamma * (1+(mhalo/M0/M12)**0.4)
    
    def approx_NFW_con_z050_Planck_vir():
        """approximates NFW concentration
        from Klypin+16 Eq. (24) for WMAP7 cosmology, Mhalo[vir], Table A3, all halos selected by mass, z=0.5
        Parameter:
        """
        c0=6.5
        gamma=0.105
        M0=1e4
        M12=1e12/0.6778    
        
        return  c0*(mhalo/M12)**-gamma * (1+(mhalo/M0/M12)**0.4)
    
    def approx_NFW_con_z100_Planck_vir():
        """approximates NFW concentration
        from Klypin+16 Eq. (24) for WMAP7 cosmology, Mhalo[vir], Table A3, all halos selected by mass, z=1.0
        Parameter:
        """

        c0=4.75
        gamma=0.100
        M0=1e3
        M12=1e12/0.6778    
        
        return  c0*(mhalo/M12)**-gamma * (1+(mhalo/M0/M12)**0.4)
    
    def approx_NFW_con_z144_Planck_vir():
        """approximates NFW concentration
        from Klypin+16 Eq. (24) for WMAP7 cosmology, Mhalo[vir], Table A3, all halos selected by mass, z=1.44
        Parameter:
        """
        c0=3.8
        gamma=0.095
        M0=210
        M12=1e12/0.6778    
        
        return  c0*(mhalo/M12)**-gamma * (1+(mhalo/M0/M12)**0.4)
    
    def approx_NFW_con_z215_Planck_vir():
        """approximates NFW concentration
        from Klypin+16 Eq. (24) for WMAP7 cosmology, Mhalo[vir], Table A3, all halos selected by mass, z=2.15
        Parameter:
        """
        c0=3.0
        gamma=0.085
        M0=43
        M12=1e12/0.6778    
        
        return  c0*(mhalo/M12)**-gamma * (1+(mhalo/M0/M12)**0.4)
    
    def approx_NFW_con_z250_Planck_vir():
        """approximates NFW concentration
        from Klypin+16 Eq. (24) for WMAP7 cosmology, Mhalo[vir], Table A3, all halos selected by mass, z=2.5
        Parameter:
        """
        c0=2.65
        gamma=0.080
        M0=18
        M12=1e12/0.6778    
        
        return  c0*(mhalo/M12)**-gamma * (1+(mhalo/M0/M12)**0.4)
    
    def approx_NFW_con_z290_Planck_vir():
        """approximates NFW concentration
        from Klypin+16 Eq. (24) for WMAP7 cosmology, Mhalo[vir], Table A3, all halos selected by mass, z=2.9
        Parameter:
        """
        c0=2.42
        gamma=0.080
        M0=9
        M12=1e12/0.6777    
        
        return  c0*(mhalo/M12)**-gamma * (1+(mhalo/M0/M12)**0.4)
    
    def approx_NFW_con_z410_Planck_vir():
        """approximates NFW concentration
        from Klypin+16 Eq. (24) for WMAP7 cosmology, Mhalo[vir], Table A3, all halos selected by mass, z=4.1
        Parameter:
        """
        c0=2.10
        gamma=0.080
        M0=1.9
        M12=1e12/0.6778    
        
        return  c0*(mhalo/M12)**-gamma * (1+(mhalo/M0/M12)**0.4)
    
    def approx_NFW_con_z540_Planck_vir():
        """approximates NFW concentration
        from Klypin+16 Eq. (24) for WMAP7 cosmology, Mhalo[vir], Table A3, all halos selected by mass, z=5.4
        Parameter:
        """
        c0=1.86
        gamma=0.080
        M0=0.42
        M12=1e12/0.6778    
        
        return  c0*(mhalo/M12)**-gamma * (1+(mhalo/M0/M12)**0.4)
    
    
    def caseSwitcher(key):
    
        choose = {

            'WMAP7200c050': approx_NFW_con_z050_WMAP7,
            'Planck200c000': approx_NFW_con_z000_Planck,
            'Planck200c050': approx_NFW_con_z050_Planck,
            'Planckvir000': approx_NFW_con_z000_Planck_vir,
            'Planckvir035': approx_NFW_con_z035_Planck_vir,
            'Planckvir050': approx_NFW_con_z050_Planck_vir,
            'Planckvir100': approx_NFW_con_z100_Planck_vir,
            'Planckvir144': approx_NFW_con_z144_Planck_vir,
            'Planckvir215': approx_NFW_con_z215_Planck_vir,
            'Planckvir250': approx_NFW_con_z250_Planck_vir,
            'Planckvir290': approx_NFW_con_z290_Planck_vir,
            'Planckvir410': approx_NFW_con_z410_Planck_vir,
            'Planckvir540': approx_NFW_con_z540_Planck_vir
            }
            
        func = choose.get(key)
        return func()

    if redshift>=0.0 and redshift<=0.175:
        z_key='000'
    elif redshift>0.175 and redshift<=0.42:
        z_key='035'
    elif redshift>0.42 and redshift<=0.6:
        z_key='050'        
    elif redshift>0.6 and redshift<=1.2:
        z_key='100'
    elif redshift>1.2 and redshift<=1.8:
        z_key='144'         
    elif redshift>1.8 and redshift<=2.25:
        z_key='215' 
    elif redshift>2.25 and redshift<=2.75:
        z_key='250'
    elif redshift>2.75 and redshift<=3.1:
        z_key='290' 
    elif redshift>3.10 and redshift<=4.60:
        z_key='410'
    elif redshift>4.6:
        z_key='540'

    print 'key:', z_key          
        
    return caseSwitcher(cosmology+overdens+z_key)

def convert_halo_mass(mhalo,
                      NFW_con,
                      orphan,
                      Dv_old=139,
                      Dv_new=200,
                      redshift=0.0):
    print '--> convert Mvir to M200c!',
    #Formular to get the halo concentration from Klypin+16 (1411.4001v2), Eq.: 24
    #Convertion of halo mass from 139 to 200 x critical overdensity follows Lokas&Mamon01 MNRAS 321, 155
    
    NFW_con[np.where(np.isfinite(NFW_con)==False)[:][0]] = 99.99


    g    = 1.0/(np.log(1.0+NFW_con)-NFW_con/(1.0+NFW_con))    
    s    = np.linspace(0.05,2,10000)
    
    s_new=np.zeros((mhalo.size,), np.float32)

    import matplotlib.pyplot as plt
    #formulae
    for i,c in enumerate(NFW_con):
        #print 'i:', i, 'con', c,
        M_in_s   = g[i]*(np.log(1.0+c*s)-c*s/(1.0+c*s))
        #rho_in_s = M_in_s/(4.0*np.pi/3.0*s**3)*(Dv_old*4.0*np.pi/3.0)
        rho_in_s = M_in_s*Dv_old/s**3
        rho_inter = np.interp(s,s,rho_in_s)      
        #try:
        mask1=np.where(rho_inter>(Dv_new-0.1))[:][0]

        #print mask1, rho_inter[mask1]
        mask2=np.where(rho_inter[mask1]<(Dv_new+0.1))[:][0]
        #print mask2, rho_inter[mask2]
    
        #print 's_new:', s[mask2[0]]
        s_new[i] = s[mask2[0]]
#        except:
#            print 'i', i, 'gal type:', orphan[i], 'with mhalo:', mhalo[i], 'not able to convert ...'
        #plt.plot(s, rho_in_s, '.')
        #plt.plot(s_new[i], rho_inter[mask2[0]], 'x')

    plt.plot(s_new[i], Dv_new, '^')
    plt.ylim((0,1000))
    plt.axhline(y=Dv_new, xmin=-100, xmax=100, color='k', ls='--', lw=2.0)
    plt.axhline(y=Dv_old, xmin=-100, xmax=100, color='k', ls='--', lw=2.0)
    plt.savefig(mycomp+'anaconda/pro/data/Galacticus_1Gpc/HMF_c_test.png', format='png', rasterized=True, dpi=100, bbox_inches='tight', pad_inches=0.1)
             
    ratio  = g*(np.log(1.0+NFW_con*s_new)-NFW_con*s_new/(1.0+NFW_con*s_new))
    print 'median ratio:', np.median(ratio[np.where(np.isfinite(ratio)==True)[:][0]]), '+/-', np.percentile(ratio[np.where(np.isfinite(ratio)==True)[:][0]], 75)-np.median(ratio[np.where(np.isfinite(ratio)==True)[:][0]]), '/', np.median(ratio[np.where(np.isfinite(ratio)==True)[:][0]])-np.percentile(ratio[np.where(np.isfinite(ratio)==True)[:][0]],25) 
    print 'CHECK!/n'
    return ratio*mhalo
   
   
def calc_survey_volume(skycoverage,
                       zmin,
                       zmax,
                       little_h_out=False):

    """calculates the volume of shell between two redshifts for a certain skycoverage.
        
        skycoverage:     float [deg2]
        zmin:            float lower redshift bound
        zmax:            float higher redshift bound
        
        keywords:
        ---------
        little_h_out:    bool factor the little-h out to get the following volume units: [Mpc3 h-3]
        
        
        returns the volume [Mpc3 h-3] or [Mpc3]
        """
    
    from astropy.cosmology import Planck as cos
    if little_h_out==True:
        unit_volume='Mpc3h-3'
        h=cos.H0.value/100.00
    else:
        unit_volume='Mpc3'
        h=1.0
    print 'h:', h, 'zmin/max', zmin,'/', zmax, cos.comoving_volume(zmax).value, '/', cos.comoving_volume(zmin).value, 'unit:', cos.comoving_volume(zmin).unit
#    print 'reshift | Volume [x1e9 Mpc3]'
#    print '0.43-0.7', float("{0:.3f}".format(skycoverage/41253.0*(cos.comoving_volume(0.7).value-cos.comoving_volume(0.43).value)*h**3/1e9))
#    print '0.43-0.5', float("{0:.3f}".format(skycoverage/41253.0*(cos.comoving_volume(0.5).value-cos.comoving_volume(0.43).value)*h**3/1e9))
#    print '0.44-0.54', float("{0:.3f}".format(skycoverage/41253.0*(cos.comoving_volume(0.54).value-cos.comoving_volume(0.44).value)*h**3/1e9)) 
#    print '0.46-0.53', float("{0:.3f}".format(skycoverage/41253.0*(cos.comoving_volume(0.53).value-cos.comoving_volume(0.46).value)*h**3/1e9))
#    print '0.56-0.63', float("{0:.3f}".format(skycoverage/41253.0*(cos.comoving_volume(0.63).value-cos.comoving_volume(0.56).value)*h**3/1e9))
#    print '0.50-0.60', float("{0:.3f}".format(skycoverage/41253.0*(cos.comoving_volume(0.60).value-cos.comoving_volume(0.50).value)*h**3/1e9))  
#    print '0.51-0.61', float("{0:.3f}".format(skycoverage/41253.0*(cos.comoving_volume(0.61).value-cos.comoving_volume(0.51).value)*h**3/1e9))
#    print '0.54-0.64', float("{0:.3f}".format(skycoverage/41253.0*(cos.comoving_volume(0.64).value-cos.comoving_volume(0.54).value)*h**3/1e9))    
#    print '0.05-0.7', float("{0:.3f}".format(skycoverage/41253.0*(cos.comoving_volume(0.7).value-cos.comoving_volume(0.05).value)*h**3/1e9))    
#    exit()
    return skycoverage/41253.0*(cos.comoving_volume(zmax).value-cos.comoving_volume(zmin).value)*h**3, unit_volume
    
def survey_VolumeSqDeg_to_Mpc(
        box_size=62.5,
        redshift=0.0,
        h=1.0,
        box_unit='Mpc',
        scale_factor=None):
    
    """Calculates the volume of the survey or box at a given redshift/comoving distance"""    
    
    if scale_factor==None:
        scale_factor = redshift_to_expfactor(redshift)
    
    if box_unit == 'deg':
        from astropy.cosmology import Planck as cos        
        comov_dist = cos.comoving_distance(redshift)    #Comoving Distance [Mpc] --> USE ALWAYS FOR VOLUMES
        V_sphere = 4*np.pi/3.0 * ((comov_dist.value)**3)
        #Volume of sphere 4Pi*r^3/3    
        V_survey = V_sphere / (129600.0 / np.pi) * box_size/h     #[Mpc+3 h-3]
    else:
        V_survey = (box_size/h)**3                   #[Mpc+3h-3]

    print   
    print 'GEOMTRY DETAILS:'
    print '--------------------------------'
    print 'redshift:', redshift, 'scale_factor:', scale_factor, 'Volume:', V_survey, '[Mpc3] hubble_par:', h

    return V_survey
    

def calc_kcor(filter_name, redshift, colour_name, colour_value):
    """
    K-corrections calculator in Python. See http://kcor.sai.msu.ru for the 
    reference (Chiligarian+2010 MRAS 405. 1409-1420). Available filter-colour combinations must be present in the 
    `coeff` dictionary keys.
    

    @type   filter_name: string    
    @param  filter_name: Name of the filter to calculate K-correction for, e.g. 
                         'u', 'g', 'r' for some of the SDSS filters, or 'J2', 
                         'H2', 'Ks2' for 2MASS filters (must be present in 
                         `coeff` dictionary)
    @type      redshift: float    
    @param     redshift: Redshift of a galaxy, should be between 0.0 and 0.5 (no
                         check is made, however)
    @type   colour_name: string    
    @param  colour_name: Human name of the colour, e.g. 'u - g', 'g - r', 
                         'V - Rc', 'J2 - Ks2' (must be present in `coeff` dictionary)
    @type  colour_value: float    
    @param colour_value: Value of the galaxy's colour, specified in colour_name    
    @rtype:              float
    @return:             K-correction in specified filter for given redshift and 
                         colour
    @version:            2012
    @author:             Chilingarian, I., Melchior. A.-L., and Zolotukhin, I.
    @license:            Simplified BSD license, see http://kcor.sai.msu.ru/license.txt

    Usage example:
    
        >>> calc_kcor('g', 0.2, 'g - r', 1.1)
        0.5209713975999992
        >>> calc_kcor('Ic', 0.4, 'V - Ic', 2.0)
        0.310069919999993
        >>> calc_kcor('H', 0.5, 'H - K', 0.1)
        -0.14983142499999502
        
    """
    coeff = {

        'B_BRc': [
            [0,0,0,0],
            [-1.99412,3.45377,0.818214,-0.630543],
            [15.9592,-3.99873,6.44175,0.828667],
            [-101.876,-44.4243,-12.6224,0],
            [299.29,86.789,0,0],
            [-304.526,0,0,0],
        ],
        
        'B_BIc': [
            [0,0,0,0],
            [2.11655,-5.28948,4.5095,-0.8891],
            [24.0499,-4.76477,-1.55617,1.85361],
            [-121.96,7.73146,-17.1605,0],
            [236.222,76.5863,0,0],
            [-281.824,0,0,0],
        ],

        'H2_H2Ks2': [
            [0,0,0,0],
            [-1.88351,1.19742,10.0062,-18.0133],
            [11.1068,20.6816,-16.6483,139.907],
            [-79.1256,-406.065,-48.6619,-430.432],
            [551.385,1453.82,354.176,473.859],
            [-1728.49,-1785.33,-705.044,0],
            [2027.48,950.465,0,0],
            [-741.198,0,0,0],
        ],

        'H2_J2H2': [
            [0,0,0,0],
            [-4.99539,5.79815,4.19097,-7.36237],
            [70.4664,-202.698,244.798,-65.7179],
            [-142.831,553.379,-1247.8,574.124],
            [-414.164,1206.23,467.602,-799.626],
            [763.857,-2270.69,1845.38,0],
            [-563.812,-1227.82,0,0],
            [1392.67,0,0,0],
        ],

        'Ic_VIc': [
            [0,0,0,0],
            [-7.92467,17.6389,-15.2414,5.12562],
            [15.7555,-1.99263,10.663,-10.8329],
            [-88.0145,-42.9575,46.7401,0],
            [266.377,-67.5785,0,0],
            [-164.217,0,0,0],
        ],

        'J2_J2Ks2': [
            [0,0,0,0],
            [-2.85079,1.7402,0.754404,-0.41967],
            [24.1679,-34.9114,11.6095,0.691538],
            [-32.3501,59.9733,-29.6886,0],
            [-30.2249,43.3261,0,0],
            [-36.8587,0,0,0],
        ],

        'J2_J2H2': [
            [0,0,0,0],
            [-0.905709,-4.17058,11.5452,-7.7345],
            [5.38206,-6.73039,-5.94359,20.5753],
            [-5.99575,32.9624,-72.08,0],
            [-19.9099,92.1681,0,0],
            [-45.7148,0,0,0],
        ],

        'Ks2_J2Ks2': [
            [0,0,0,0],
            [-5.08065,-0.15919,4.15442,-0.794224],
            [62.8862,-61.9293,-2.11406,1.56637],
            [-191.117,212.626,-15.1137,0],
            [116.797,-151.833,0,0],
            [41.4071,0,0,0],
        ],

        'Ks2_H2Ks2': [
            [0,0,0,0],
            [-3.90879,5.05938,10.5434,-10.9614],
            [23.6036,-97.0952,14.0686,28.994],
            [-44.4514,266.242,-108.639,0],
            [-15.8337,-117.61,0,0],
            [28.3737,0,0,0],
        ],

        'Rc_BRc': [
            [0,0,0,0],
            [-2.83216,4.64989,-2.86494,0.90422],
            [4.97464,5.34587,0.408024,-2.47204],
            [-57.3361,-30.3302,18.4741,0],
            [224.219,-19.3575,0,0],
            [-194.829,0,0,0],
        ],

        'Rc_VRc': [
            [0,0,0,0],
            [-3.39312,16.7423,-29.0396,25.7662],
            [5.88415,6.02901,-5.07557,-66.1624],
            [-50.654,-13.1229,188.091,0],
            [131.682,-191.427,0,0],
            [-36.9821,0,0,0],
        ],

        'U_URc': [
            [0,0,0,0],
            [2.84791,2.31564,-0.411492,-0.0362256],
            [-18.8238,13.2852,6.74212,-2.16222],
            [-307.885,-124.303,-9.92117,12.7453],
            [3040.57,428.811,-124.492,-14.3232],
            [-10677.7,-39.2842,197.445,0],
            [16022.4,-641.309,0,0],
            [-8586.18,0,0,0],
        ],

        'V_VIc': [
            [0,0,0,0],
            [-1.37734,-1.3982,4.76093,-1.59598],
            [19.0533,-17.9194,8.32856,0.622176],
            [-86.9899,-13.6809,-9.25747,0],
            [305.09,39.4246,0,0],
            [-324.357,0,0,0],
        ],

        'V_VRc': [
            [0,0,0,0],
            [-2.21628,8.32648,-7.8023,9.53426],
            [13.136,-1.18745,3.66083,-41.3694],
            [-117.152,-28.1502,116.992,0],
            [365.049,-93.68,0,0],
            [-298.582,0,0,0],
        ],

        'FUV_FUVNUV': [
            [0,0,0,0],
            [-0.866758,0.2405,0.155007,0.0807314],
            [-1.17598,6.90712,3.72288,-4.25468],
            [135.006,-56.4344,-1.19312,25.8617],
            [-1294.67,245.759,-84.6163,-40.8712],
            [4992.29,-477.139,174.281,0],
            [-8606.6,316.571,0,0],
            [5504.2,0,0,0],
        ],

        'FUV_FUVu': [
            [0,0,0,0],
            [-1.67589,0.447786,0.369919,-0.0954247],
            [2.10419,6.49129,-2.54751,0.177888],
            [15.6521,-32.2339,4.4459,0],
            [-48.3912,37.1325,0,0],
            [37.0269,0,0,0],
        ],

        'g_gi': [
            [0,0,0,0],
            [1.59269,-2.97991,7.31089,-3.46913],
            [-27.5631,-9.89034,15.4693,6.53131],
            [161.969,-76.171,-56.1923,0],
            [-204.457,217.977,0,0],
            [-50.6269,0,0,0],
        ],

        'g_gz': [
            [0,0,0,0],
            [2.37454,-4.39943,7.29383,-2.90691],
            [-28.7217,-20.7783,18.3055,5.04468],
            [220.097,-81.883,-55.8349,0],
            [-290.86,253.677,0,0],
            [-73.5316,0,0,0],
        ],

        'g_gr': [
            [0,0,0,0],
            [-2.45204,4.10188,10.5258,-13.5889],
            [56.7969,-140.913,144.572,57.2155],
            [-466.949,222.789,-917.46,-78.0591],
            [2906.77,1500.8,1689.97,30.889],
            [-10453.7,-4419.56,-1011.01,0],
            [17568,3236.68,0,0],
            [-10820.7,0,0,0],
        ],

        'H_JH': [
            [0,0,0,0],
            [-1.6196,3.55254,1.01414,-1.88023],
            [38.4753,-8.9772,-139.021,15.4588],
            [-417.861,89.1454,808.928,-18.9682],
            [2127.81,-405.755,-1710.95,-14.4226],
            [-5719,731.135,1284.35,0],
            [7813.57,-500.95,0,0],
            [-4248.19,0,0,0],
        ],

        'H_HK': [
            [0,0,0,0],
            [0.812404,7.74956,1.43107,-10.3853],
            [-23.6812,-235.584,-147.582,188.064],
            [283.702,2065.89,721.859,-713.536],
            [-1697.78,-7454.39,-1100.02,753.04],
            [5076.66,11997.5,460.328,0],
            [-7352.86,-7166.83,0,0],
            [4125.88,0,0,0],
        ],

        'i_gi': [
            [0,0,0,0],
            [-2.21853,3.94007,0.678402,-1.24751],
            [-15.7929,-19.3587,15.0137,2.27779],
            [118.791,-40.0709,-30.6727,0],
            [-134.571,125.799,0,0],
            [-55.4483,0,0,0],
        ],

        'i_ui': [
            [0,0,0,0],
            [-3.91949,3.20431,-0.431124,-0.000912813],
            [-14.776,-6.56405,1.15975,0.0429679],
            [135.273,-1.30583,-1.81687,0],
            [-264.69,15.2846,0,0],
            [142.624,0,0,0],
        ],

        'J_JH': [
            [0,0,0,0],
            [0.129195,1.57243,-2.79362,-0.177462],
            [-15.9071,-2.22557,-12.3799,-2.14159],
            [89.1236,65.4377,36.9197,0],
            [-209.27,-123.252,0,0],
            [180.138,0,0,0],
        ],

        'J_JK': [
            [0,0,0,0],
            [0.0772766,2.17962,-4.23473,-0.175053],
            [-13.9606,-19.998,22.5939,-3.99985],
            [97.1195,90.4465,-21.6729,0],
            [-283.153,-106.138,0,0],
            [272.291,0,0,0],
        ],

        'K_HK': [
            [0,0,0,0],
            [-2.83918,-2.60467,-8.80285,-1.62272],
            [14.0271,17.5133,42.3171,4.8453],
            [-77.5591,-28.7242,-54.0153,0],
            [186.489,10.6493,0,0],
            [-146.186,0,0,0],
        ],

        'K_JK': [
            [0,0,0,0],
            [-2.58706,1.27843,-5.17966,2.08137],
            [9.63191,-4.8383,19.1588,-5.97411],
            [-55.0642,13.0179,-14.3262,0],
            [131.866,-13.6557,0,0],
            [-101.445,0,0,0],
        ],

        'NUV_NUVr': [
            [0,0,0,0],
            [2.2112,-1.2776,0.219084,0.0181984],
            [-25.0673,5.02341,-0.759049,-0.0652431],
            [115.613,-5.18613,1.78492,0],
            [-278.442,-5.48893,0,0],
            [261.478,0,0,0],
        ],

        'NUV_NUVg': [
            [0,0,0,0],
            [2.60443,-2.04106,0.52215,0.00028771],
            [-24.6891,5.70907,-0.552946,-0.131456],
            [95.908,-0.524918,1.28406,0],
            [-208.296,-10.2545,0,0],
            [186.442,0,0,0],
        ],

        'r_gr': [
            [0,0,0,0],
            [1.83285,-2.71446,4.97336,-3.66864],
            [-19.7595,10.5033,18.8196,6.07785],
            [33.6059,-120.713,-49.299,0],
            [144.371,216.453,0,0],
            [-295.39,0,0,0],
        ],

        'r_ur': [
            [0,0,0,0],
            [3.03458,-1.50775,0.576228,-0.0754155],
            [-47.8362,19.0053,-3.15116,0.286009],
            [154.986,-35.6633,1.09562,0],
            [-188.094,28.1876,0,0],
            [68.9867,0,0,0],
        ],

        'u_ur': [
            [0,0,0,0],
            [10.3686,-6.12658,2.58748,-0.299322],
            [-138.069,45.0511,-10.8074,0.95854],
            [540.494,-43.7644,3.84259,0],
            [-1005.28,10.9763,0,0],
            [710.482,0,0,0],
        ],

        'u_ui': [
            [0,0,0,0],
            [11.0679,-6.43368,2.4874,-0.276358],
            [-134.36,36.0764,-8.06881,0.788515],
            [528.447,-26.7358,0.324884,0],
            [-1023.1,13.8118,0,0],
            [721.096,0,0,0],
        ],

        'u_uz': [
            [0,0,0,0],
            [11.9853,-6.71644,2.31366,-0.234388],
            [-137.024,35.7475,-7.48653,0.655665],
            [519.365,-20.9797,0.670477,0],
            [-1028.36,2.79717,0,0],
            [767.552,0,0,0],
        ],

        'Y_YH': [
            [0,0,0,0],
            [-2.81404,10.7397,-0.869515,-11.7591],
            [10.0424,-58.4924,49.2106,23.6013],
            [-0.311944,84.2151,-100.625,0],
            [-45.306,3.77161,0,0],
            [41.1134,0,0,0],
        ],

        'Y_YK': [
            [0,0,0,0],
            [-0.516651,6.86141,-9.80894,-0.410825],
            [-3.90566,-4.42593,51.4649,-2.86695],
            [-5.38413,-68.218,-50.5315,0],
            [57.4445,97.2834,0,0],
            [-64.6172,0,0,0],
        ],

        'z_gz': [
            [0,0,0,0],
            [0.30146,-0.623614,1.40008,-0.534053],
            [-10.9584,-4.515,2.17456,0.913877],
            [66.0541,4.18323,-8.42098,0],
            [-169.494,14.5628,0,0],
            [144.021,0,0,0],
        ],

        'z_rz': [
            [0,0,0,0],
            [0.669031,-3.08016,9.87081,-7.07135],
            [-18.6165,8.24314,-14.2716,13.8663],
            [94.1113,11.2971,-11.9588,0],
            [-225.428,-17.8509,0,0],
            [197.505,0,0,0],
        ],

        'z_uz': [
            [0,0,0,0],
            [0.623441,-0.293199,0.16293,-0.0134639],
            [-21.567,5.93194,-1.41235,0.0714143],
            [82.8481,-0.245694,0.849976,0],
            [-185.812,-7.9729,0,0],
            [168.691,0,0,0],
        ],

    }

    c = coeff[filter_name + '_' + colour_name.replace(' - ', '')]

    kcor_array=np.zeros((colour_value[:].size), dtype=np.float32)

    i=0
    while i<colour_value[:].size:

        #print 'i:', i, 'kcor:', kcor
        for x, a in enumerate(c):
            for y, b in enumerate(c[x]):
                #print 'x:', x, 'y:', y, 'a:', a, 'b:', b
                kcor_array[i] += c[x][y] * redshift**x * colour_value[i]**y

        i+=1
    
    return kcor_array
    
def return_lists_by_galaxy_type(data):
    central_ind   = (np.where(data['orphan'] == 0))[0] 
    data_centrals = data[central_ind]
    sat_ind       = (np.where(data['orphan'] == 1))[0] 
    data_sats     = data[sat_ind]        
    orphan_ind    = (np.where(data['orphan'] == 2))[0]
    data_orphans  = data[orphan_ind]
    
    return data_centrals, data_sats, data_orphans

def quick_plot(data_x,
               data_y,
               data,
               x_title='X',
               y_title='Y'):

    import matplotlib.pyplot as plt          
    fig = plt.figure(figsize=(12,10), dpi=150)
    myax = fig.add_subplot(111) 
       
    myax.plot(np.log10(data_x),
            np.log10(data_y),
            lw=2.0,
            marker='o',
            ls='-',
            color='r')

    fig.text(0.5, 0.03, x_title, ha='center', fontsize=20)
    fig.text(0.01, 0.5, y_title, va='center', rotation='vertical', fontsize=20)                           

    myax.tick_params(axis='x', which='major', top='on', bottom='on', pad=10, labelsize=16, length=12, width=2.0, direction='in', zorder=20) 
    myax.tick_params(axis='x', which='minor', top='on', bottom='on', pad=10, labelsize=16, length=8, width=2.0, direction='in', zorder=20) 

    myax.tick_params(axis='y', which='minor', left='on', right='on', pad=10, labelsize=16, length=8, width=2.0, direction='in', zorder=20) 
    myax.tick_params(axis='y', which='major', left='on', right='on', pad=10, labelsize=16, length=12, width=2.0, direction='in', zorder=20) 
          
    #plt.savefig(mycomp+'anaconda/pro/myRun/plots/plotXY/plot_only.png', format='png', rasterized=True, dpi=100, bbox_inches='tight', pad_inches=0.05) 
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(mycomp+'anaconda/pro/myRun/plots/plotXY/plot_only.pdf')
    plt.savefig(pp, format='pdf', rasterized=True, dpi=50, pad_inches=0.05, transparent=True)
    pp.close()
    print 'QUCIK PLOT ... done!'
  
    rand_mstar= choose_random_sample(data, 5000)                
    data_int_test=np.interp(rand_mstar, data_x, data_y)
    
    fig = plt.figure(figsize=(12,10), dpi=10)
    myax = fig.add_subplot(111) 

    myax.plot(np.log10(data_x),
                np.log10(data_y),
                lw=2.0,
                marker='o',
                ls='-',
                color='k') 

    myax.plot(np.log10(rand_mstar),
                np.log10(data_int_test),
                lw=2.0,
                marker='.',
                ls='',
                color='r')         

    pp = PdfPages(mycomp+'anaconda/pro/myRun/plots/plotXY/plot_only2.pdf')
    plt.savefig(pp, format='pdf', rasterized=True, dpi=10, pad_inches=0.05, transparent=True)
    pp.close() 
 
    
    
def downsample_SMF(data):
    
    myData = aD.ArangeData()
 
    data_mask=data[np.where(data['mstar']>1e10)[0][:]]
    print data.shape    
    #data_mask=data_mask[np.where(data_mask['mstar']<1.1e12)[0][:]]
        
#
#        #Guo+13 cut
#        #(r-i)>0.679 -0.082(Mi+20)
    #cut=0.679-0.082*(data_mask['MAB_dA_total_i']-5*np.log10(0.6777)+20.0)      
    #data_mask=data_mask[np.where(data_mask['mAB_dA_total_r']-data_mask['mAB_dA_total_i']>cut)[0][:]]
    #data_mask=data_mask[np.where(data_mask['sfr']/data_mask['mstar']<1e-11)[0][:]]
    #print data_mask.shape    
    #standard g-i>2.35
#    try:
#        data_mask=data_mask[np.where(data_mask['mAB_dA_total_g']-data_mask['mAB_dA_total_i']>2.35)[0][:]]        
#    except:
#        data_mask=data_mask[np.where(data_mask['mAB_total_g']-data_mask['mAB_total_i']>2.35)[0][:]]

    cut=0.679-0.082*(data_mask['MAB_dA_total_i']-5*np.log10(0.6777)+20.0)
    data_mask=data_mask[np.where(data_mask['mAB_dA_total_r']-data_mask['mAB_dA_total_i']>cut)[0][:]]         
    print 'selection: --> r-i > ', data_mask.shape
    
    histo_data=myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/SMF/SMF_CMASS_SPALL_z_0.5_0.6_spall_por_merged_DR12v4_compl_sample_wg_fixed_25bins_1e10.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)        
    #histo_data=myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/SMF/SMF_CMASS_SPALL_z_0.5_0.6_spall_por_merged_DR12v4_compl_sample_wg_fixed_35bins_1e10_-0.2.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)        

    #histo_model=myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/SMF/SMF_Galacticus_1Gpc_z_0.56_tarsel_new_mags_g-i_gt_2.35_fixed_40bins.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)        
    #histo_model=myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/SMF/SMF_Galacticus_1Gpc_z_0.56_tarsel_new_mags_red_fixed_40bins.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)        

    #histo_model=myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/SMF/SMF_SAG_1Gpc_z_0.56_tarsel_new_red_fixed_40bins.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)        
    #histo_model=myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/SMF/SMF_SAGE_1Gpc_z_0.56_tarsel_v3_mags_run_1238_red_fixed_40bins.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)        

    #histo_model=myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/SMF/SMF_LGALAXIES_500Mpc_z_0.56_fixed_35bins_1e10.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)        
    #histo_model=myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/SMF/SMF_LGALAXIES_500Mpc_z_0.56_tarsel_35bins_2e10_no.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)        

    #histo_model=myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/SMF/SMF_Galacticus_1Gpc_run2_z_0.56_tarsel_fixed_25bins_1e10_test2.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)        
    histo_model=myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/SMF/SMF_Galacticus_1Gpc_run2_z_0.56_tarsel_Guo13_cut_fixed_25bins_1e10.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)        

    #histo_model=myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/SMF/SMF_Galacticus_400Mpc_z_0.55_tarsel_g-i_gt_2.35_fixed_25bins_1e10.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)        
    #histo_model=myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/SMF/SMF_Galacticus_400Mpc_z_0.55_tarsel_fixed_25bins_1e10_test2.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)        
    #histo_model=myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/myRun/histos/SMF/SMF_Galacticus_400Mpc_z_0.55_tarsel_Guo13_cut_fixed_25bins_1e10.txt', data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float, skiprow=2)        


    frac_data=histo_data[:,1]/histo_model[:,1]

    print 'frac_data min/max:', min(frac_data), '/', max(frac_data)

    quick_plot(histo_model[:,0], frac_data, data['mstar'], x_title='log10 (Mstar [Msun])', y_title='log10 P(Mstar)')
          
    rand_sample = np.random.random(data_mask.size)
    #print 'rand_sample min/max:', min(rand_sample), '/', max(rand_sample)
   
    data_int    = np.interp(data_mask['mstar'], histo_model[:,0], frac_data)
    print 'data_int min/max:', min(data_int), '/', max(data_int)
   
    data_mask = data_mask[np.where(rand_sample<data_int)[0][:]] 
    print data_mask.shape             
    test      = np.in1d(data['hostid'], data_mask['hostid'])             
    data      = data[np.where(test==True)[0][:]]
           
    print 'new sample:', data.size, 'data log10 mstar min/max:' , np.log10(min(data['mstar'])), '/', np.log10(max(data['mstar']))

    #data2cross = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/data/Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_haloid.txt', data_format='ASCII', data_shape='shaped', delim='  ', mydtype=np.uint64, skiprow=2) 
    #print data2cross[0:10,:]
#
    #test = np.in1d(data['hostid'],data2cross[:,1], invert=True)
    #data=data[test]
    #cut=0.679-0.082*(data['MAB_dA_total_i']+20.0)
    #data=data[np.where(data['mAB_dA_total_r']-data['mAB_dA_total_i']>cut)[0][:]]        
    
    #data=data[np.where(data['orphan']<2)[:][0]]
    #data=data[np.where(data['mstar']>3e10)[:][0]]
    
    #print data.size
    
    return data

def give_galaxy_type_info(data):
    
       # myData = aD.ArangeData()
       
        data_cents=data[np.where(data==0)[0][:]]    
        data_sats=data[np.where(data==1)[0][:]]
        data_os=data[np.where(data==2)[0][:]]

        print 'all [N]:\t', data.size, 'n: ', format(data.size/(400/0.6778)**3, '.4f')
        print 'counts [N]\tcents:', data_cents.size, '\ttot sats:', data_sats.size+data_os.size, '\tno-sats:', data_sats.size, '\torphans:', data_os.size
        print 'fractions \tcents:', format(data_cents.size/float(data.size), '.3f'), \
                            '\ttot sats:',  format((data_sats.size+data_os.size)/float(data.size), '.3f'), \
                            '\tno-sats:',  format(data_sats.size/float(data.size), '.3f'), \
                            '\torphans:',  format(data_os.size/float(data.size), '.3f'), \
                            '\tof tot sats:',  format(data_os.size/float(data_os.size+data_sats.size), '.3f')
    
    
def give_env_info(data):
    
        myData = aD.ArangeData()
       
        #data['rbulge']/=data['rdisk']
        data_cents=data[np.where(data['orphan']==0)[0][:]]
        data_cents=data_cents[np.where(data_cents['pop']==2)[0][:]]        
        data_sats=data[np.where(data['orphan']==1)[0][:]]
        data_os=data[np.where(data['orphan']==2)[0][:]]
        #print 'all:', data.size, 'centrals:', data_cents.size, 'sats:', data_sats.size+data_os.size, 'orphans:', data_os.size
        #print 'fractions cent:', data_cents.size/float(data.size), 'sats:', (data_sats.size+data_os.size)/float(data.size), 'orphans in sats:', data_os.size/float(data_sats.size)
        
        data_test = data_cents
        data_test_node=data_test[np.where(data_test['env_1024']==3)[0][:]]
        data_test_fila=data_test[np.where(data_test['env_1024']==2)[0][:]]
        data_test_sheet=data_test[np.where(data_test['env_1024']==1)[0][:]]
        
        print 'all [N]:\t', data_test.size
        print 'counts [N]\tn:', data_test_node.size, '\tf:', data_test_fila.size, '\tsh:', data_test_sheet.size
        print 'fractions [%]\tn:', format(100*data_test_node.size/float(data_test.size), '.0f'), '\t\tf:',  format(100*data_test_fila.size/float(data_test.size), '.0f'), '\t\tsh:',  format(100*data_test_sheet.size/float(data_test.size), '.0f')




        data_highZ=data_cents[np.where(data_cents['zcold']>9.5)[0][:]]

        data_test_node=data_highZ[np.where(data_highZ['env_1024']==3)[0][:]]
        data_test_fila=data_highZ[np.where(data_highZ['env_1024']==2)[0][:]]
        data_test_sheet=data_highZ[np.where(data_highZ['env_1024']==1)[0][:]]        

        print 'HIGH-Zcold>9.5 all [N]:\t', data_highZ.size, '\t', format(100*data_highZ.size/float(data_cents.size), '.0f'), '%'
        print 'counts [N]\tn:', data_test_node.size, '\tf:', data_test_fila.size, '\tsh:', data_test_sheet.size
        print 'fractions [%]\tn:', format(100*data_test_node.size/float(data_test.size), '.0f'), '\t\tf:',  format(100*data_test_fila.size/float(data_test.size), '.0f'), '\t\tsh:',  format(100*data_test_sheet.size/float(data_test.size), '.0f')
        
        
        data_lowZ=data_cents[np.where(data_cents['zcold']<=9.5)[0][:]]

        data_test_node=data_lowZ[np.where(data_lowZ['env_1024']==3)[0][:]]
        data_test_fila=data_lowZ[np.where(data_lowZ['env_1024']==2)[0][:]]
        data_test_sheet=data_lowZ[np.where(data_lowZ['env_1024']==1)[0][:]]  
        
        print 'LOW-Zcold<=9.5 all [N]:\t', data_lowZ.size, '\t', format(100*data_lowZ.size/float(data_cents.size), '.0f'), '%'
        print 'counts [N]\tn:', data_test_node.size, '\tf:', data_test_fila.size, '\tsh:', data_test_sheet.size
        print 'fractions [%]\tn:', format(100*data_test_node.size/float(data_test.size), '.0f'), '\t\tf:',  format(100*data_test_fila.size/float(data_test.size), '.0f'), '\t\tsh:',  format(100*data_test_sheet.size/float(data_test.size), '.0f')

        exit()   
    
    
        data_test_node=data_test[np.where(data_test['env_1024']==3)[0][:]]
        data_test_fila=data_test[np.where(data_test['env_1024']==2)[0][:]]
        data_test_sheet=data_test[np.where(data_test['env_1024']==1)[0][:]]
        print 'all [N]:\t', data_test.size
        print 'counts [N]\tn:', data_test_node.size, '\tf:', data_test_fila.size, '\tsh:', data_test_sheet.size
        print 'fractions [%]\tn:', format(100*data_test_node.size/float(data_test.size), '.0f'), '\t\tf:',  format(100*data_test_fila.size/float(data_test.size), '.0f'), '\t\tsh:',  format(100*data_test_sheet.size/float(data_test.size), '.0f')

        pop_dict=          {'all': 10,    'pop2': 2, 'pop1': 1}       
        pop_dict_operators={'all': '<', 'pop2': '==', 'pop1': '=='}
        
        from statsmodels import robust
        
        for pop in ['all', 'pop1','pop2']:
            mydata = myData.selectData2Compute(data_cents, 
                                            selected_col='pop', 
                                            operator=pop_dict_operators[pop], 
                                            condition=pop_dict[pop])
            
            data_test_node=mydata[np.where(mydata['env_1024']==3)[0][:]]
            data_test_fila=mydata[np.where(mydata['env_1024']==2)[0][:]]
            data_test_sheet=mydata[np.where(mydata['env_1024']==1)[0][:]]
            
            mydata_dict={'k': data_test_node, 'f': data_test_fila, 's':data_test_sheet} 
            
            print '\nsample:', pop, 'size:', mydata.size, '\n'
            
            
            #for prop in ['mhalo_200c', 'sfr', 'ssfr', 'zcold', 'mstar', 'cgf', 'rhalf_mass', 'rbulge', 'mbh']:
            for prop in ['rbulge']:                
                for envir in ['k', 'f', 's']:
                    data2use=mydata_dict[envir]
                    if prop=='zcold':
                        print prop.ljust(10), 'size:', str(data2use.size).ljust(6), envir, ' ', 'median:', float("{0:.2f}".format(np.nanmedian(data2use[prop]))),'\t', 'MAD:', float("{0:.2f}".format(robust.mad(data2use[prop]))), '\t16/84',  float("{0:.2f}".format(np.percentile(data2use[prop], 16))), '/', float("{0:.2f}".format(np.percentile(data2use[prop], 84))), '\t25/75',  float("{0:.2f}".format(np.percentile(data2use[prop], 25))), '/', float("{0:.2f}".format(np.percentile(data2use[prop], 75))), '\t10/90',  float("{0:.2f}".format(np.percentile(data2use[prop], 10))), '/', float("{0:.2f}".format(np.percentile(data2use[prop], 90))) 
                    elif prop.startswith('lr'):
                        print prop.ljust(10)+'/rhalf', 'size:', str(data2use.size).ljust(6), envir, ' ', 'median:', float("{0:.2f}".format(np.nanmedian(data2use[prop])*1000.0)),'\t', 'MAD:', float("{0:.2f}".format(robust.mad(data2use[prop])*1000.0)), '\t16/84',  float("{0:.2f}".format(np.percentile(data2use[prop], 16)*1000.0)), '/', float("{0:.2f}".format(np.percentile(data2use[prop], 84)*1000.0)), '\t25/75',  float("{0:.2f}".format(np.percentile(data2use[prop], 25)*1000.0)), '/', float("{0:.2f}".format(np.percentile(data2use[prop], 75)*1000.0)), '\t10/90',  float("{0:.2f}".format(np.percentile(data2use[prop], 10)*1000.0)), '/', float("{0:.2f}".format(np.percentile(data2use[prop], 90)*1000.0))            
                    else:
                        print prop.ljust(10), 'size:', str(data2use.size).ljust(6), envir, ' ', 'median:', float("{0:.2f}".format(np.log10(np.nanmedian(data2use[prop])))),'\t', 'MAD: ', float("{0:.2f}".format(np.log10(robust.mad(data2use[prop])))), '\t16/84',  float("{0:.2f}".format(np.log10(np.percentile(data2use[prop], 16)))), '/', float("{0:.2f}".format(np.log10(np.percentile(data2use[prop], 84)))), '\t25/75',  float("{0:.2f}".format(np.log10(np.percentile(data2use[prop], 25)))), '/', float("{0:.2f}".format(np.log10(np.percentile(data2use[prop], 75)))), '\t10/90',  float("{0:.2f}".format(np.log10(np.percentile(data2use[prop], 10)))), '/', float("{0:.2f}".format(np.log10(np.percentile(data2use[prop], 90)))) 
                print '\n'
            print '---------------------------'

        exit()

def manual_selection(data):
    print 'MANUAL SELECTION!'
#    try:
#        data=data[np.where((data['sfr_spheroid']+data['sfr_disk'])>1e-4)[0][:]]
#    except:
#        data=data[np.where(data['sfr']>1e-4)[0][:]]
#    print 'SFR CUT ', data.shape
    
    data=data[np.where(data['orphan']<2)[0][:]]
    print 'data.shape after orphan:', data.shape
 
    return data

def reduce_rand_sats(data):
    print 'Redruce randomly galaxies from the catalog!'
    data_sats=data[np.where(data['orphan']>=1)[0][:]]
    
    
    sat_frac=data_sats.size/float(data.size)
    
    ngal2reduce=data_sats.size-np.int(data_sats.size/sat_frac*0.10)
    print ngal2reduce
    
    data_sats=np.sort(data_sats, order=['mhalo'])
    
    data2reduce=choose_random_sample(data_sats,ngal2reduce)
    
    data2reduce=data_sats[0:ngal2reduce]

    test, indices_basic, indices_check = np.intersect1d(data['hostid'], data2reduce['hostid'], return_indices=True)       

    data=np.delete(data, np.s_[indices_basic], axis=0)                

                   
    return data

def reduce_rand(data):
    print 'Redruce randomly galaxies from the catalog!'
    
    
    data2reduce=choose_random_sample(data,int(data.size*0.9))
    
    test, indices_basic, indices_check = np.intersect1d(data['hostid'], data2reduce['hostid'], return_indices=True)       

    data=np.delete(data, np.s_[indices_basic], axis=0)                

                   
    return data

def manual_manipulation(data):
    print 'MANUAL MANIPULATION!'
#    key='mAB_dA_total_'
#    for band in ['u','g','r','i']:
#        print key+band
#        data[key+band]+=5*np.log10(0.704)
        
    #data['mstar']*=1.16
    data=data[np.where(data['mstar']>3e10)[0][:]]        
        
    return data

def df_to_sarray(df):
    """
    from: https://stackoverflow.com/questions/13187778/convert-pandas-dataframe-to-numpy-array
    Convert a pandas DataFrame object to a numpy structured array.
    This is functionally equivalent to but more efficient than
    np.array(df.to_array())

    :param df: the data frame to convert
    :return: a numpy structured array representation of df
    """

    #df=df.drop(['row_id','index'], axis=1)
    
    v = df.values
    cols = df.columns
    
    types = [(cols[i].encode(), df[k].dtype.type) for (i, k) in enumerate(cols)]

    dtype = np.dtype(types)
    z = np.zeros(v.shape[0], dtype)
    for (i, k) in enumerate(z.dtype.names):
        z[k] = v[:, i]
    return z    

def show_mag_info(data):
    
    print 'MAB'
    from statsmodels import robust     
    for band in ['u','g','r','i']:
        diff = data['MAB_dA_total_'+band]-data['MAB_total_'+band]
        print 'median:', np.nanmedian(diff), 'std:', robust.mad(diff)

    print 'mAB'
    for band in ['u','g','r','i']:
        diff = data['mAB_dA_total_'+band]-data['mAB_total_'+band]
        print 'median:', np.nanmedian(diff), 'std:', robust.mad(diff)            
        
    exit()
    
def test_haloids(data):

    import pandas as pd
    data1 = pd.read_csv(mycomp+'anaconda/pro/data/Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3_haloid.txt', skiprows=2, names=['haloid', 'hostid', 'orphan', 'env_512', 'env_1024'], sep='  ')
    data_test_against = df_to_sarray(data1)
    
    #print data_test_against[0:1]
    print np.info(data_test_against)
    
#    data2 = pd.read_csv(mycomp+'anaconda/pro/myRun/histos/sfr2z/Galacticus_SDSS/haloids/sfr2z_Galacticus_1Gpc_z_0.56_tarsel_CUT3_Contreras+13_mcold_SFH_method2_cents_lowZcold-highMstar_haloids.txt', skiprows=2, names=['haloid', 'hostid', 'orphan', 'mhalo', 'mstar'], sep='\t')
#    data_test = df_to_sarray(data2)

#    print 'Gal-dens:\n'
#    for prop in ['mhalo', 'mstar', 'sfr', 'ssfr', 'zcold', 'mbh']:
#        print prop, '\t\t', "{0:.2f}".format(np.log10(np.nanmedian(data_test_against[prop]))), '\t', "{0:.2f}".format(np.log10(np.nanmedian(data_test_against[prop]))-np.log10(np.nanpercentile(data_test_against[prop], 32))), '/', "{0:.2f}".format(np.log10(np.nanpercentile(data_test_against[prop], 68))-np.log10(np.nanmedian(data_test_against[prop])))

#    print 'Gal2-dens:\n'
#    for prop in ['mhalo', 'mstar', 'sfr', 'ssfr', 'zcold', 'mbh']:
#        if prop=='zcold': 10**data[prop]
#        print prop, '\t\t', "{0:.2f}".format(np.log10(np.nanmedian(data[prop]))), '\t', "{0:.2f}".format(np.log10(np.nanmedian(data[prop]))-np.log10(np.nanpercentile(data[prop], 32))), '/', "{0:.2f}".format(np.log10(np.nanpercentile(data[prop], 68))-np.log10(np.nanmedian(data[prop])))



#    import numpy.lib.recfunctions as rcfuncs        
#    data = rcfuncs.append_fields([data], ['haloid_Gal_dens','hostid_Ga_-dens','orphan_Gal_dens', 'env_1024', 'env_512'],[data['haloid'],data['haloid'],data['orphan'],data['haloid'],data['haloid']], usemask=False)

    #test = np.in1d(data_test_against['hostid'], data['hostid'])
#    #print test
#    
    #data_test_against=data_test_against[np.where(test==True)[:][0]]

    for ID in data['hostid']:
        print 'ID:', ID
        try:
            data[['haloid_Gal_dens','hostid_Gal_dens', 'orphan_Gal_dens', 'env_512', 'env_1024']][np.where(data['hostid']==ID)[:][0]]=data_test_against[['haloid','hostid','orphan','env_512','env_1024']][np.where(data_test_against['hostid']==ID)[:][0]]
        except:
            print 'no ID found!'
    #data[::-1].sort(order=['hostid'], axis=0)
    #data_test_against[::-1].sort(order=['hostid'], axis=0)    
    
    #data[['haloid_Gal_dens','hostid_Gal_dens', 'orphan_Gal_dens', 'env_512', 'env_1024']]=data_test_against[['haloid','hostid','orphan','env_512','env_1024']]
    
    print data[['hostid','hostid_Gal_dens', 'orphan', 'orphan_Gal_dens', 'env_512', 'env_1024']][0:3]
    #exit()

    return data
    #return data_SDSS    

def calc_residuals(data, count, redshift, sample, prop):
    #print 'count:', count, 'redshift:', redshift, 'sample:', sample, 'prop:', prop 
    
    import pandas as pd
    data1 = pd.read_csv(mycomp+'anaconda/pro/data/Galacticus_1Gpc/SFH/SFH_Index_'+str(redshift)+'_'+sample+'_props_M1.txt', skiprows=2, names=['haloid', 'parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar',  'SHMR', 'mcold',  'Mzgas',  'zcold', 'g-i', 'sfr', 'ssfr'], sep='\t')
    data_M1 = df_to_sarray(data1)
        
    #print data_test_against[0:1]
    #print np.info(data_M1)
    
    data2 = pd.read_csv(mycomp+'anaconda/pro/data/Galacticus_1Gpc/SFH/SFH_Index_'+str(redshift)+'_'+sample+'_props_M2.txt', skiprows=2, names=['haloid', 'parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar',  'SHMR', 'mcold',  'Mzgas',  'zcold', 'g-i', 'sfr', 'ssfr'], sep='\t')
    data_M2 = df_to_sarray(data2)   
    #print np.info(data_M2)

    #data_M1['mhalo']=data_M1['mstar']/data_M1['mhalo']
    #data_M2['mhalo']=data_M2['mstar']/data_M2['mhalo']
    
    test_frac_in_M1 = np.in1d(data_M1['hostid'], data_M2['hostid'])
    test = test_frac_in_M1[np.where(test_frac_in_M1==True)[:][0]]
    
    print 'data_M1 / test', data_M1.size , '/', test.size

    data_M1_found=data_M1[np.where(test_frac_in_M1==True)[:][0]]

    myprop_unit_map={'mstar':  'log10(mstar [Msun])',\
                   'mhalo':    'log10(mhalo [Msun])',\
                   'mcold':    'log10(mcold [Msun])',\
                   'zcold':    'zcold [-]',\
                   'Mzgas':    'log10(Mzcold [Msun])',\
                   'sfr':      'log10(SFR [Msunyr-1])',\
                   'ssfr':     'log10(sSFR [yr-1])',\
                   'g-i':      'g-i [-]',\
                   'SHMR':     'log10(mstar/mhalo) [-])'}  
    
    i=1
    if prop=='SHMR':
        prop='mhalo'
        myprop_unit_map.update({'mhalo': 'log10(mstar/mhalo) [-])'})
        
    myheader='('+str(i)+') z ('+str(i+1)+') N_M1found/N_M1 ('+str(i+2)+') N_M1found/N_M2 ('+str(i+3)+') N_M1/N_M1'+\
             '('+str(i+4)+') '+myprop_unit_map[prop]+' M1found ('+str(i+5)+') +dy ('+str(i+6)+') -dy ('+str(i+7)+') N_M1found [-] '+\
             '('+str(i+8)+') '+myprop_unit_map[prop]+' M1 ('+str(i+9)+') +dy ('+str(i+10)+') -dy ('+str(i+11)+') N_M1 [-] '+\
             '('+str(i+12)+') '+myprop_unit_map[prop]+' M2 ('+str(i+13)+') +dy ('+str(i+14)+') -dy ('+str(i+15)+') N_M2 [-] '+\
             '('+str(i+16)+') '+myprop_unit_map[prop]+' median(M1_found)/median(M1)-1 [-] '+\
             '('+str(i+17)+') '+myprop_unit_map[prop]+' median(M1_found)/median(M2)-1 [-] '+\
             '('+str(i+18)+') '+myprop_unit_map[prop]+' median(M1)/median(M2)-1 [-] '+\
             '('+str(i+19)+') '+myprop_unit_map[prop]+' (median(M1found)-median(M1))/median(M2) [-]'


    data[count, 0] = redshift
    data[count, 1] = data_M1_found.size/float(data_M1.size)
    data[count, 2] = data_M1_found.size/float(data_M2.size)
    data[count, 3] = data_M1.size/float(data_M2.size)
    data[count, 4] = np.nanmedian(data_M1_found[prop])
    data[count, 5] = np.nanmedian(data_M1_found[prop])-np.nanpercentile(data_M1_found[prop], 32)
    data[count, 6] = np.nanpercentile(data_M1_found[prop], 68)-np.nanmedian(data_M1_found[prop])
    data[count, 7] = data_M1_found.size
    data[count, 8] = np.nanmedian(data_M1[prop])
    data[count, 9] = np.nanmedian(data_M1[prop])-np.nanpercentile(data_M1[prop], 32)
    data[count, 10] = np.nanpercentile(data_M1[prop], 68)-np.nanmedian(data_M1[prop])
    data[count, 11] = data_M1.size
    data[count, 12] = np.nanmedian(data_M2[prop])
    data[count, 13] = np.nanmedian(data_M2[prop])-np.nanpercentile(data_M2[prop], 32)
    data[count, 14] = np.nanpercentile(data_M2[prop], 68)-np.nanmedian(data_M2[prop])
    data[count, 15] = data_M2.size
    data[count, 16] = np.nanmedian(data_M1_found[prop])/np.nanmedian(data_M1[prop])-1.0    
    data[count, 17] = np.nanmedian(data_M1_found[prop])/np.nanmedian(data_M2[prop])-1.0
    data[count, 18] = np.nanmedian(data_M1[prop])/np.nanmedian(data_M2[prop])-1.0
    data[count, 19] = (np.nanmedian(data_M1_found[prop])-np.nanmedian(data_M1[prop]))/np.nanmedian(data_M2[prop]) 

#    if prop=='zcold' or prop=='g-i':
#        min_prop_M1=str("{0:.2f}".format(min(data_M1[prop])))
#        max_prop_M1=str("{0:.2f}".format(max(data_M1[prop])))
#        min_prop_M2=str("{0:.2f}".format(min(data_M2[prop])))
#        max_prop_M2=str("{0:.2f}".format(max(data_M2[prop])))
#        min_prop_M1_found=str("{0:.2f}".format(min(data_M1_found[prop])))
#        max_prop_M1_found=str("{0:.2f}".format(max(data_M1_found[prop])))
#    else:
#        min_prop_M1=str("{0:.2f}".format(min(np.log10(data_M1[prop]))))
#        max_prop_M1=str("{0:.2f}".format(max(np.log10(data_M1[prop]))))
#        min_prop_M2=str("{0:.2f}".format(min(np.log10(data_M2[prop]))))
#        max_prop_M2=str("{0:.2f}".format(max(np.log10(data_M2[prop]))))
#        min_prop_M1_found=str("{0:.2f}".format(min(np.log10(data_M1_found[prop]))))
#        max_prop_M1_found=str("{0:.2f}".format(max(np.log10(data_M1_found[prop]))))        

    #header_prefix=', '+myprop_unit_map[prop]+' M1found min/max: '+min_prop_M1_found+'/'+max_prop_M1_found+', M1 min/max: '+min_prop_M1+'/'+max_prop_M1+', M2 min/max: '+min_prop_M2+'/'+max_prop_M2+'\n'
        
    return data,  myheader

def test_methods(redshift, sample, props, bad_mstar_count):

    print 'z:', redshift,
        
    import pandas as pd
    data1 = pd.read_csv(mycomp+'anaconda/pro/data/Galacticus_1Gpc/SFH/SFH_Index_'+str(redshift)+'_'+sample+'_props_M1.txt', skiprows=2, names=['haloid', 'parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar',  'SHMR', 'mcold',  'Mzgas',  'zcold', 'g-i', 'sfr', 'ssfr'], sep='\t')
    data_M1 = df_to_sarray(data1)
        
    #print data_test_against[0:1]
    #print np.info(data_M1)
    
    data2 = pd.read_csv(mycomp+'anaconda/pro/data/Galacticus_1Gpc/SFH/SFH_Index_'+str(redshift)+'_'+sample+'_props_M2.txt', skiprows=2, names=['haloid', 'parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar',  'SHMR', 'mcold',  'Mzgas',  'zcold', 'g-i', 'sfr', 'ssfr'], sep='\t')
    data_M2 = df_to_sarray(data2)   
    #print np.info(data_M2)

    data_M1['mhalo']=data_M1['mstar']/data_M1['mhalo']
    data_M2['mhalo']=data_M2['mstar']/data_M2['mhalo']
    
    test_frac_in_M1 = np.in1d(data_M1['hostid'], data_M2['hostid'])
    test = test_frac_in_M1[np.where(test_frac_in_M1==True)[:][0]]
    
    print 'data_M1 / test', data_M1.size , '/', test.size
    
#    stats_props_all='\multicolumn{1}{c|}{\multirow{2}{*}{\parbox{0.04\linewidth}{\centering'+str(redshift)+'}}}\t&M1&'
#    for prop in ['mhalo', 'mstar', 'SHMR', 'sfr', 'ssfr', 'mcold', 'Mzgas', 'zcold', 'g-i']:
#        if prop=='zcold' or prop=='g-i':
#            stats_props_all+=str("{0:.2f}".format(np.nanmedian(data_M1[prop])))+'$_{-'+str("{0:.2f}".format(np.nanmedian(data_M1[prop])-np.nanpercentile(data_M1[prop], 32)))+'}^{+'+str("{0:.2f}".format(np.nanpercentile(data_M1[prop], 68)-np.nanmedian(data_M1[prop])))+'}$\t& ' 
#        else:
#            stats_props_all+=str("{0:.2f}".format(np.log10(np.nanmedian(data_M1[prop]))))+'$_{-'+str("{0:.2f}".format(np.log10(np.nanmedian(data_M1[prop]))-np.log10(np.nanpercentile(data_M1[prop], 32))))+'}^{+'+str("{0:.2f}".format(np.log10(np.nanpercentile(data_M1[prop], 68))-np.log10(np.nanmedian(data_M1[prop]))))+'}$\t& ' 
#
#    stats_props_all+='\multicolumn{1}{c}{\multirow{2}{*}{\parbox{0.04\linewidth}{\centering'+str("{0:.2f}".format(test.size/float(data_M1.size)))+'}}}///'
#
#    stats_props_all+='\n\cline{2-11}\n&M2\t&'
#    for prop in ['mhalo', 'mstar', 'SHMR', 'sfr', 'ssfr', 'mcold', 'Mzgas', 'zcold', 'g-i']:
#        if prop=='zcold' or prop=='g-i':
#            stats_props_all+=str("{0:.2f}".format(np.nanmedian(data_M2[prop])))+'$_{-'+str("{0:.2f}".format(np.nanmedian(data_M2[prop])-np.nanpercentile(data_M2[prop], 32)))+'}^{+'+str("{0:.2f}".format(np.nanpercentile(data_M2[prop], 68)-np.nanmedian(data_M1[prop])))+'}$\t& ' 
#        else:
#            stats_props_all+=str("{0:.2f}".format(np.log10(np.nanmedian(data_M2[prop]))))+'$_{-'+str("{0:.2f}".format(np.log10(np.nanmedian(data_M2[prop]))-np.log10(np.nanpercentile(data_M2[prop], 32))))+'}^{+'+str("{0:.2f}".format(np.log10(np.nanpercentile(data_M2[prop], 68))-np.log10(np.nanmedian(data_M2[prop]))))+'}$\t& ' 
#
#    stats_props_all+='\n///\n\hline'
#    #print "{0:.2f}".format(test.size/float(data_M1.size))
#
#    #print data_M1[0:10]
#
    stats_props_found_M1='\multicolumn{1}{c|}{\multirow{4}{*}{\parbox{0.03\linewidth}{\centering'+str(redshift)+'}}}\t&M1&' 
#    for prop in ['mhalo', 'mstar', 'SHMR', 'sfr', 'ssfr', 'mcold', 'Mzgas', 'zcold', 'g-i']:
#        if prop=='zcold' or prop=='g-i':
#            stats_props_found_M1+=str("{0:.2f}".format(np.nanmedian(data_M1[prop])))+'$_{-'+str("{0:.2f}".format(np.nanmedian(data_M1[prop])-np.nanpercentile(data_M1[prop], 32)))+'}^{+'+str("{0:.2f}".format(np.nanpercentile(data_M1[prop], 68)-np.nanmedian(data_M1[prop])))+'}$\t& ' 
#        else:
#            stats_props_found_M1+=str("{0:.2f}".format(np.log10(np.nanmedian(data_M1[prop]))))+'$_{-'+str("{0:.2f}".format(np.log10(np.nanmedian(data_M1[prop]))-np.log10(np.nanpercentile(data_M1[prop], 32))))+'}^{+'+str("{0:.2f}".format(np.log10(np.nanpercentile(data_M1[prop], 68))-np.log10(np.nanmedian(data_M1[prop]))))+'}$\t& ' 
#
#
#    stats_props_found_M1+='\multicolumn{1}{c}{\multirow{4}{*}{\parbox{0.03\linewidth}{\centering'+str("{0:.2f}".format(test.size/float(data_M1.size)))+'}}}///\n'

    data_M1_found=data_M1[np.where(test_frac_in_M1==True)[:][0]]
#    stats_props_found_M1+='\cline{2-11}\n&M2\t&'
#    for prop in ['mhalo', 'mstar', 'SHMR', 'sfr', 'ssfr', 'mcold', 'Mzgas', 'zcold', 'g-i']:
#        if prop=='zcold' or prop=='g-i':
#            stats_props_found_M1+=str("{0:.2f}".format(np.nanmedian(data_M2[prop])))+'$_{-'+str("{0:.2f}".format(np.nanmedian(data_M2[prop])-np.nanpercentile(data_M2[prop], 32)))+'}^{+'+str("{0:.2f}".format(np.nanpercentile(data_M2[prop], 68)-np.nanmedian(data_M2[prop])))+'}$\t& ' 
#        else:
#            stats_props_found_M1+=str("{0:.2f}".format(np.log10(np.nanmedian(data_M2[prop]))))+'$_{-'+str("{0:.2f}".format(np.log10(np.nanmedian(data_M2[prop]))-np.log10(np.nanpercentile(data_M2[prop], 32))))+'}^{+'+str("{0:.2f}".format(np.log10(np.nanpercentile(data_M2[prop], 68))-np.log10(np.nanmedian(data_M2[prop]))))+'}$\t& ' 
#
#    stats_props_found_M1+='///\n'
#
#    stats_props_found_M1+='\cline{2-11}\n&M1$_{/rm{true}}$\t&'
#    for prop in ['mhalo', 'mstar', 'SHMR', 'sfr', 'ssfr', 'mcold', 'Mzgas', 'zcold', 'g-i']:
#        if prop=='zcold' or prop=='g-i':
#            stats_props_found_M1+=str("{0:.2f}".format(np.nanmedian(data_M1_found[prop])))+'$_{-'+str("{0:.2f}".format(np.nanmedian(data_M1_found[prop])-np.nanpercentile(data_M1_found[prop], 32)))+'}^{+'+str("{0:.2f}".format(np.nanpercentile(data_M1_found[prop], 68)-np.nanmedian(data_M1_found[prop])))+'}$\t& ' 
#        else:
#            stats_props_found_M1+=str("{0:.2f}".format(np.log10(np.nanmedian(data_M1_found[prop]))))+'$_{-'+str("{0:.2f}".format(np.log10(np.nanmedian(data_M1_found[prop]))-np.log10(np.nanpercentile(data_M1_found[prop], 32))))+'}^{+'+str("{0:.2f}".format(np.log10(np.nanpercentile(data_M1_found[prop], 68))-np.log10(np.nanmedian(data_M1_found[prop]))))+'}$\t& ' 
#
#    stats_props_found_M1+='///\n'
#    
#    stats_props_found_M1+='\cline{2-11}\n&f$_{/rm{M1_{/rm{true}}}/rm{M2}}$\t&'
#    for prop in ['mhalo', 'mstar', 'SHMR', 'sfr', 'ssfr', 'mcold', 'Mzgas', 'zcold', 'g-i']:
#        data_M1_found_median=np.nanmedian(data_M1_found[prop])/np.nanmedian(data_M2[prop])
#        if prop=='zcold' or prop=='g-i':
#            stats_props_found_M1+=str("{0:.2f}".format(data_M1_found_median))+'\t& ' 
#        else:
#            stats_props_found_M1+=str("{0:.2f}".format(np.log10(data_M1_found_median)))+'\t& ' 
#
#    stats_props_found_M1+='///\n\hline'
    
    stats_props_tab=str(redshift)+'  '+str("{0:.2f}".format(data_M1_found.size/float(data_M1.size)))+'  '

    i=3

    myprop_unit_map={'mstar':  'log10(mstar [Msun])',\
                   'mhalo':    'log10(mhalo [Msun])',\
                   'mcold':    'log10(mcold [Msun])',\
                   'zcold':    'zcold [-]',\
                   'Mzgas':    'log10(Mzcold [Msun])',\
                   'sfr':      'log10(SFR [Msunyr-1])',\
                   'ssfr':     'log10(sSFR [yr-1])',\
                   'g-i':      'g-i [-]',\
                   'SHMR':     'log10(mstar/mhalo) [-])'}  
    
    myheader=''
    for prop in [props]:

        myheader+='('+str(i)+') '+myprop_unit_map[prop] +'M1 ('+str(i+1)+') +dy ('+str(i+2)+') -dy ('+str(i+3)+') '+myprop_unit_map[prop] +' M2 ('+str(i+4)+') +dy ('+str(i+5)+') -dy ('+str(i+6)+') '+myprop_unit_map[prop] +' M1_found ('+str(i+7)+') +dy ('+str(i+8)+') -dy ('+str(i+9)+') '+myprop_unit_map[prop]+'  median(M1_found)/median(M2) [-] ('+str(i+10)+') '+myprop_unit_map[prop]+' median(M1_found)/median(M1) [-] '
        if prop=='SHMR':
            prop='mhalo'
            myprop_unit_map.update({'mhalo': 'log10(mstar/mhalo) [-])'})

        #print data_M1_found[prop].size
        
        if min(data_M1_found['mstar'])<1e6:
            print '------------------------------> bad mstar!!!!', min(data_M1_found['mstar'])
            bad_mstar_count+=1
            #print 'total count: ', bad_mstar_count

        data_M1_found_median=np.nanmedian(data_M1_found[prop])       
        data_M1_found_32=np.nanpercentile(data_M1_found[prop], 32)
        data_M1_found_68=np.nanpercentile(data_M1_found[prop], 68)
            
        median_data_M1_found_data_M2=data_M1_found_median/np.nanmedian(data_M2[prop])
        median_data_M1_found_data_M1=data_M1_found_median/np.nanmedian(data_M1[prop])
        
        if prop=='zcold' or prop=='g-i':
            stats_props_tab+=str("{0:.2f}".format(np.nanmedian(data_M1[prop])))+'  '+str("{0:.2f}".format(np.nanmedian(data_M1[prop])-np.nanpercentile(data_M1[prop], 32)))+'  '+str("{0:.2f}".format(np.nanpercentile(data_M1[prop], 68)-np.nanmedian(data_M1[prop])))+'  '+\
                             str("{0:.2f}".format(np.nanmedian(data_M2[prop])))+'  '+str("{0:.2f}".format(np.nanmedian(data_M2[prop])-np.nanpercentile(data_M2[prop], 32)))+'  '+str("{0:.2f}".format(np.nanpercentile(data_M2[prop], 68)-np.nanmedian(data_M2[prop])))+'  '+\
                             str("{0:.2f}".format(data_M1_found_median))+'  '+str("{0:.2f}".format(data_M1_found_median-data_M1_found_32))+'  '+str("{0:.2f}".format(data_M1_found_68-data_M1_found_median))+'  '+\
                             str("{0:.2f}".format(median_data_M1_found_data_M2))+'  '+str("{0:.2f}".format(median_data_M1_found_data_M1))+'  '
            min_prop_M1=str("{0:.2f}".format(min(data_M1[prop])))
            max_prop_M1=str("{0:.2f}".format(max(data_M1[prop])))
            min_prop_M2=str("{0:.2f}".format(min(data_M2[prop])))
            max_prop_M2=str("{0:.2f}".format(max(data_M2[prop])))
            min_prop_M1_found=str("{0:.2f}".format(min(data_M1_found[prop])))
            max_prop_M1_found=str("{0:.2f}".format(max(data_M1_found[prop])))            
        else:
            stats_props_tab+=str("{0:.2f}".format(np.log10(np.nanmedian(data_M1[prop]))))+'  '+str("{0:.2f}".format(np.log10(np.nanmedian(data_M1[prop]))-np.log10(np.nanpercentile(data_M1[prop], 32))))+'  '+str("{0:.2f}".format(np.log10(np.nanpercentile(data_M1[prop], 68))-np.log10(np.nanmedian(data_M1[prop]))))+'  '+\
                             str("{0:.2f}".format(np.log10(np.nanmedian(data_M2[prop]))))+'  '+str("{0:.2f}".format(np.log10(np.nanmedian(data_M2[prop]))-np.log10(np.nanpercentile(data_M2[prop], 32))))+'  '+str("{0:.2f}".format(np.log10(np.nanpercentile(data_M2[prop], 68))-np.log10(np.nanmedian(data_M2[prop]))))+'  '+\
                             str("{0:.2f}".format(np.log10(data_M1_found_median)))+'  '+str("{0:.2f}".format(np.log10(data_M1_found_median)-np.log10(data_M1_found_32)))+'  '+str("{0:.2f}".format(np.log10(data_M1_found_68)-np.log10(data_M1_found_median)))+'  '+\
                             str("{0:.2f}".format(median_data_M1_found_data_M2))+'  '+str("{0:.2f}".format(median_data_M1_found_data_M1))+'  '
            min_prop_M1=str("{0:.2f}".format(np.log10(min(data_M1[prop]))))
            max_prop_M1=str("{0:.2f}".format(np.log10(max(data_M1[prop]))))
            min_prop_M2=str("{0:.2f}".format(np.log10(min(data_M2[prop]))))
            max_prop_M2=str("{0:.2f}".format(np.log10(max(data_M2[prop]))))
            min_prop_M1_found=str("{0:.2f}".format(np.log10(min(data_M1_found[prop]))))
            max_prop_M1_found=str("{0:.2f}".format(np.log10(max(data_M1_found[prop]))))
        i+=11

    header_prefix=', '+myprop_unit_map[prop]+' M1 min/max: '+min_prop_M1+'/'+max_prop_M1+', M2 min/max: '+min_prop_M2+'/'+max_prop_M2+', M1_found min/max: '+min_prop_M1_found+'/'+max_prop_M1_found+', bad mstar count: '+str(bad_mstar_count)+'\n'
        
    return stats_props_tab[:-2], stats_props_found_M1[:-1], "{0:.2f}".format(test.size/float(data_M1.size)), header_prefix, myheader, bad_mstar_count

def property_stats(data):
    
    #data=data[np.where(data['orphan']==0)[:][0]]
    size=data.size
      
    pop=2
    if pop==1 or pop==2:
        data=data[np.where(data['pop']==pop)[:][0]]

      
    print 'pop:', pop, 'ngal:', data.size, '\b', "{0:.2f}".format(100.0/size*data.size), '% of all galaxies', size
#
    for env, env_name in zip([-99, 0,1,2,3],['all', 'void','sheet','filament','knot']):
    #for env, env_name in zip([-99],['all']):            
        print '\n'
        if env_name!='all':
            sample=data[np.where(data['env_1024']==env)[:][0]]
        else:
            sample=data

        print 'env', env, 'name:', env_name, 'ngal:', sample.size, '\t', "{0:.0f}".format(100.0/data.size*sample.size),'% of all galaxies'
        print '////////////////////////////////////////////////'
    
        print 'prop\t\tmedian\t32th / 68th\n-------------------------------\n'
        for prop in ['mhalo', 'mhalo_200c', 'NFW_con', 'mstar', 'zcold', 'sfr', 'ssfr', 'mcold', 'Mzgas', 'cgf', 'mbh','MAB_dA_total_g', 'MAB_dA_total_r', 'MAB_dA_total_i', 'mAB_dA_total_cut_r_i','mAB_dA_total_cut_g_r','mAB_dA_total_cut_g_i', 'mean_age_stars_disk', 'mean_age_stars_spheroid']:
        #for prop in ['mcold', 'mstar', 'zcold', 'mbh']:
            print prop,
            if prop.find('AB')==-1 and prop!='zcold' and prop!='NFW_con' and prop.find('age')==-1:
                print '\t\t', "{0:.2f}".format(np.log10(np.nanmedian(sample[prop]))), '\t', "{0:.2f}".format(np.log10(np.nanmedian(sample[prop]))-np.log10(np.nanpercentile(sample[prop], 32))), '/', "{0:.2f}".format(np.log10(np.nanpercentile(sample[prop], 68))-np.log10(np.nanmedian(sample[prop])))
            else:
                print '\t', "{0:.2f}".format(np.nanmedian(sample[prop])), '\t', "{0:.2f}".format(np.nanmedian(sample[prop])-np.nanpercentile(sample[prop], 32)), '/', "{0:.2f}".format(np.nanpercentile(sample[prop], 68)-np.nanmedian(sample[prop]))


    print '--> DONE!\n'
    exit()
   