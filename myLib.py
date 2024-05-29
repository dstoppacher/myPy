import myFuncs as mF
import arangeData as aD
import outputData as oD
import numpy as np
import os
system_info = os.getcwd()
start = system_info.find('anaconda')
mycomp = system_info[:start]

import time
ts = time.time()
import datetime 
date_time_stamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H:%M:%S')

#Template of function documentation
"""Description of function

    input
    =========
        
    input1:   input1
               
    output
    =========
    
    output1:    output1
    
    usage:

"""

def calc_norm_vector(data):
    
    """ calculated the norm of the vector data input

        data[x,y,z]: either vector or vector array, three column
        
    """
    #print data.dtype.names
    
    product=0
    for key in data.dtype.names:
        #print key
        product+= data[key]**2       
    
    return (product)**0.5

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
    #myData = aD.ArangeData()
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


def conv_radius_pc_to_arcsec(radius): 
        
    """converts the radius of an object to arcsec at a distance of 10pc see http://arcsec2parsec.joseonorbe.com/about.html
    input: radius in pc
    output: arcsec of radius at standard distance of 10pc
    """
    return 206265.0*(radius/10.0)


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

     
def conv_fluxdensity_to_ABmag(flux_density):
    
    """Converts flux density [erg s-1 cm-2 Hz-1] to apparent magnitudes in the AB System"""
                 
    #mAB = -2.5 log10(fluxdensity) - 48.60
    #or mAB = -2.5 log10(fluxdensity / 3631 Jy)
    #log10 ( f[erg s-1 cm-2 Hz-1] / 3631 Jy[10-23 erg s-1 cm-2 Hz-1] ) --> log10(f/3631Jy)
                       
    return -2.5*np.log10(flux_density/3631e-23)
 
def conv_fluxJy_to_ABmag(fluxJy):
    
    """Converts flux in Jy [10-23 erg s-1 cm-2 Hz-1] to apparent magnitudes in the AB System"""
                 
    #or mAB = -2.5 log10(flux Jy / 3631 Jy)
    #log10 ( f[10-23 erg s-1 cm-2 Hz-1] / 3631 Jy[10-23 erg s-1 cm-2 Hz-1] ) --> log10(f/3631Jy)
                       
    return -2.5*np.log10(fluxJy/3631.0)                             

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
 
def convert_SDSS_to_Johnson_Lupton2005(filter_band,g,r):
    """Acording to Lupton+05 https://www.sdss3.org/dr10/algorithms/sdssUBVRITransform.php
    
           B = u - 0.8116*(u - g) + 0.1313;  sigma = 0.0095
           B = g + 0.3130*(g - r) + 0.2271;  sigma = 0.0107
        
           V = g - 0.2906*(u - g) + 0.0885;  sigma = 0.0129
           V = g - 0.5784*(g - r) - 0.0038;  sigma = 0.0054
        
           R = r - 0.1837*(g - r) - 0.0971;  sigma = 0.0106
           R = r - 0.2936*(r - i) - 0.1439;  sigma = 0.0072
        
           I = r - 1.2444*(r - i) - 0.3820;  sigma = 0.0078
           I = i - 0.3780*(i - z)  -0.3974;  sigma = 0.0063
    
    """
    
    if filter_band=='V':
        return g - 0.5784*(g - r) - 0.0038
    elif filter_band=='B':
        return g + 0.3130*(g - r) + 0.2271
    elif filter_band=='R':
        return r - 0.1837*(g - r) - 0.0971
    elif filter_band=='I':
        return r - 1.2444*(r - i) - 0.3820
    else:
        return False
   
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
            if np.sum(data[name])==0.0 and name.find('cut')==-1 and name!='Z' and name.find('AB')!=-1:     
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
                    #lum_dist=cd.luminosity_distance(redshift, **fidcosmo)
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
            #pass
            print name, 'is not able to convert ...\n'  
    
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
            print 'i:', i, 'name:', id_col_array['name'+str(i)], 'corr_type:', id_col_array['corr_type'+str(i)], 'unit_corr:', id_col_array['unit_corr'+str(i)] 
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
  
    i=0
    bands=[]
    conv_name=''
    while i<int(id_col_array['nr_entries']):
        try:
            name=id_col_array['name'+str(i)]
            if np.sum(data[name])==0.0 and name.find('cut')==-1 and name.find('AB')!=-1 or name.find('ag')!=-1 :           
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
    
                from cosmolopy import cparam, cd
                fidcosmo = cparam.Planck(flat=True, extras=True)
                                
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

def find_nearest_value(array,value):
    idx,val = min(enumerate(array), key=lambda x: abs(x[1]-value))
    return idx

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
                                  mylog,
                                  set_header_Topcat_format='yes'):
    
    sel_col_list = [] 
    mydataformat = ''
    check_size   = True

    print 'start filter_data_before_write2file():'
    print '----------------------------------------'
    print ' '

    i=0
    count=0
    
    #print myconds_array
    #exit()
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
          
            if set_header_Topcat_format=='yes':
                myheader+=col_name+'['+str(myconds_array[name+'_unit'])+']('+str(count+1)+') '
                
            else:
                myheader+='('+str(count+1)+') '+col_name+' ['+str(myconds_array[name+'_unit'])+'] '
                
            mydataformat+=str(myconds_array[name+'_format'])
            log_name_des = '             name'+str(count+1)+': '
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
            
        elif str(myconds_array[name+'_format'])=='%s':
            max_data = '-'
            min_data = '-'
            
        else:
            if name.startswith('delta_mhal')  or name.startswith('Mdot') or name.startswith('E')\
                or name.startswith('mPart_') or name.startswith('Mgas') or name.startswith('mgas') or name.startswith('DeltaSHMR')\
                or name.startswith('Mpseu') or name.startswith('L_')\
                or name.startswith('ang') or name.startswith('satelliteMergeTim') or name.startswith('mbasi')\
                or name.startswith('sfr') or name.startswith('ssfr') or name.startswith('r')\
                or name.startswith('mstar') or name.startswith('mhalo') or name.startswith('mcold')\
                or name.startswith('mhot') or name.startswith('mbh') or name.startswith('bh_acc')\
                or name.startswith('zgas') or name.startswith('zstar') or name.startswith('zhot')\
                or name.startswith('Mz') or name=='SHMR' or name=='bheff' or name.startswith('Sigma') or name.startswith('sfe') or name=='vpec_norm':         
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

def give_galaxy_type_info(data,
                          box_side_lenght,
                          hubble_par):
    
       # myData = aD.ArangeData()
       
        data_cents=data[np.where(data==0)[0][:]]    
        data_sats=data[np.where(data==1)[0][:]]
        data_os=data[np.where(data==2)[0][:]]

        print 'all [N]:\t', data.size, 'n: ', format(data.size/(box_side_lenght/hubble_par)**3, '.4f')
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
    Function converts a pandas data frame to a numpy stuctured array
    """    

    v = df.to_numpy()
    cols = df.columns
    #print(cols)
    
    types = [(k, df[k].dtype.type) for (i, k) in enumerate(cols)]
    #print(types)
    
    z = np.zeros((v[:,0].size,), dtype=types)
    
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
  

def binaryToDecimal(binary): 
      
    decimal, i = 0, 0, 0
    while(binary != 0): 
        dec = binary % 10
        decimal = decimal + dec * pow(2, i) 
        binary = binary//10
        i += 1
    
    #print 'decimal:', decimal   
    return decimal
  
def simple_hdf52struct_array(path,
                             cols=[]):
    """Reads a set of columns dedicated by col from and hdf5-file to a structured array
    thereby it uses the dtype of each column to generate the array:
        
        input
        =========
            
        path:   the path to the hdf5-file
        col:    keyword, exact name of the column to be accessed in the hdf5-file,
                if not set, the column name from the hdf5-file is used and all columns
                in the hdf5-file are read to the structured array
                
        output
        =========
        
        returns a structured array with column names and dtypes and size as in the hdf5-file.
        The array dimension is taken from the hdf5-file. The shape of the structured array
        needs to be identical for all column which should be read in.
    
    """
    import h5py as hdf5
    
    f=hdf5.File(path, "r")
    print '\nSimle read-in hdf5 data to structured array! -->',
    print 'path:', path
    
    if cols==[]:
        cols=f.keys()

    dt = np.dtype([(k, f[k].dtype) for k in cols])

    data = np.zeros(f[cols[0]].shape, dtype=dt)
   
    #print np.info(data)
    #print 'properties read into structured array!'
    for name in cols:
        #print 'property:', name,
        for dim in range(len(f[cols[0]].shape)):
            #print 'dim:', dim
            try:
                data[name][dim,:]=f[name][dim]
            except:
               data[name]=f[name][:] 
    
    f.close()
        
    return data
    
def get_prop_unit_map_latex():
    
    myprop_unit_map={'mstar':           ('$\log_{10}(M_{*}$ [$M_{\odot}$])', 'log'),\
                   'mstar_30kpc':     ('$\log_{10}(M_{\mathrm{*,30kpc}}$ [$M_{\odot}$])', 'log'),\
                   'mstar_1.5ropt':   ('$\log_{10}(M_{\mathrm{*,1.5}R_{\mathrm{opt}}}$ [$M_{\odot}$])', 'log'),\
                   'mstar_half_reff_gas_disk':      ('$\log_{10}(M_{\mathrm{*,1/2,eff,disk}}$ [$M_{\odot}$])', 'log'),\
                   'mPart_stars':     ('$\log_{10}(M_{\mathrm{*,part}}$ [$M_{\odot}$])', 'log'),\
                   'mPart_gasSF':     ('$\log_{10}(M_{\mathrm{SF-gas,part}}$ [$M_{\odot}$])', 'log'),\
                'mhalo':             ('$\log_{10}(M_{\mathrm{vir}}$ [$M_{\odot}$])', 'log'),\
                   'mhalo_200c':        ('$\log_{10}(M_{\mathrm{200c}}$ [$M_{\odot}$])', 'log'),\
                   'mcold':             ('$\log_{10}(M_{\mathrm{cold}}$ [$M_{\odot}$])', 'log'),\
                   'mhot':              ('$\log_{10}(M_{\mathrm{hot}}$ [$M_{\odot}$])', 'log'),\
                   'mbh':               ('$\log_{10}(M_{\mathrm{BH}}$ [$M_{\odot}$])', 'log'),\
                   'zcold':             ('$Z_{\mathrm{cold}}$', 'lin'),\
                   'zcold_gasSF':       ('$12+(\mathrm{O/H})$', 'lin'),\
                    'spinParameter':       ('spin_{\mathrm{DM}}', 'lin'),\
                       'OH_gas_30kpc':       ('$12+\log_{10}(\mathrm{O/H})$', 'lin'),\
                   'zstar':             ('$Z_{*}$', 'lin'),\
                   'zcold_zstar':       ('$Z_{\mathrm{cold}}$-$Z_{*}$', 'lin'),\
                   'zgasSF_H':          ('$Z_{\mathrm{SF-gas,H}}$', 'lin'),\
                   'zgasSF_O':          ('$Z_{\mathrm{SF-gas,O}}$', 'lin'),\
                'Mzgas':             ('$\log_{10}(M_{Z_{\mathrm{cold}}}$ [$M_{\odot}$])', 'log'),\
                   'Mzstar':            ('$\log_{10}(M_{Z_{*}}$ [$M_{\odot}$])', 'log'),\
                   'Mzhot_halo':        ('$\log_{10}(M_{Z_{\mathrm{hot,halo}}}$ [$M_{\odot}$])', 'log'),\
                   'Mgas_HI':           ('$\log_{10}(M_{\mathrm{H}_\mathrm{I}}$ [$M_{\odot}$])', 'log'),\
                   'Mgas_H2':          ('$\log_{10}(M_{\mathrm{H}_2}$ [$M_{\odot}$])', 'log'),\
                   'Mgas_HI_30kpc_GK11': ('$\log_{10}(M_{\mathrm{H}_\mathrm{I}\mathrm{,30kpc,GK11}}$ [$M_{\odot}$])', 'log'),\
                   'Mgas_H2_30kpc_GK11': ('$\log_{10}(M_{\mathrm{H}_2\mathrm{,30kpc,GK11}}$ [$M_{\odot}$])', 'log'),\
                   'Mgas_HIH2_30kpc_GK11': ('$\log_{10}(M_{\mathrm{H}_\mathrm{I}+\mathrm{H}_2\mathrm{,30kpc,GK11}}$ [$M_{\odot}$])', 'log'),\
                   'mgas_30kpc':        ('$\log_{10}(M_{\mathrm{gas,30kpc}}$ [$M_{\odot}$])', 'log'),\
                   'mgas_1.5ropt':      ('$\log_{10}(M_{\mathrm{gas,1.5}R_{\mathrm{opt}}}$ [$M_{\odot}$])', 'log'),\
                   'mgas_half_reff_gas_disk':      ('$\log_{10}(M_{\mathrm{gas,1/2,eff,disk}}$ [$M_{\odot}$])', 'log'),\
                   'sfr':               ('$\log_{10}(SFR$ [$M_{\odot}$ yr$^{-1}$])', 'log'),\
                   'sfr_1.5ropt':       ('$\log_{10}(SFR_{1.5R_{\mathrm{opt}}}$ [$M_{\odot}$ yr$^{-1}$])', 'log'),\
                   'sfr_30kpc':         ('$\log_{10}(SFR_{\mathrm{30kpc}}$ [$M_{\odot}$ yr$^{-1}$])', 'log'),\
                   'sfr_gasSF':         ('$\log_{10}(SFR_{\mathrm{SF-gas}}$ [$M_{\odot}$ yr$^{-1}$])', 'log'),\
                   'ssfr':              ('$\log_{10}(sSFR$ [yr$^{-1}$])', 'log'),\
                   'ssfr_1.5ropt':      ('$\log_{10}(sSFR_{1.5R_{\mathrm{opt}}}$ [yr$^{-1}$])', 'log'),\
                   'ssfr_30kpc':        ('$\log_{10}(sSFR_{30kpc}}$ [yr$^{-1}$])', 'log'),\
                   'ssfr_gasSF':        ('$\log_{10}(sSFR_{\mathrm{SF-gas}}$ [yr$^{-1}$])', 'log'),\
                'g-i':               ('($g$-$i$)', 'lin'),\
                   'r-i':               ('($r$-$i$)', 'lin'),\
                   'Sersic_n':          ('$n_{\mathrm{Sersic}}$', 'lin'),\
                    'BvT':               ('$M_{\mathrm{bulge}}/M_*$', 'lin'),\
                'SHMR':              ('$\log_{10}$(SHMR)', 'log'),\
                   'SHMR_1.5ropt':      ('$\log_{10}(M_{*,1.5R_{\mathrm{opt}}}/M_{\mathrm{200c}}$)', 'log'),\
                   'SHMR_30kpc':        ('$\log_{10}(M_{\mathrm{*,30kpc}}/M_{\mathrm{200c}}$)', 'log'),\
               'fbar':               ('$\log_{10}(f_{\mathrm{bar}})$', 'log'),\
                'cgf':               ('$\log_{10}$(CGF)', 'log'),\
                    'fatom_30kpc_GK11':               ('$f_{\mathrm{atom,30kpc,GK11}}$', 'lin'),\
                        'fmol_30kpc_GK11':               ('$f_{\mathrm{mol,30kpc,GK11}}$', 'lin'),\
                   'cgf_30kpc':         ('$\log_{10}(M_{\mathrm{cold}}/M_*)_{\mathrm{30kpc}}$', 'log'),\
                   'cgf_HI_30kpc':         ('$\log_{10}(M_{\mathrm{cold}}/M_*)_{\mathrm{H}_\mathrm{I}\mathrm{,30kpc}}$', 'log'),\
                   'cgf_HIH2_30kpc':         ('$\log_{10}(M_{\mathrm{cold}}/M_*)_{\mathrm{H}_\mathrm{I}+\mathrm{H}_2\mathrm{,30kpc}}$', 'log'),\
                   'cgf_1.5ropt':         ('$\log_{10}(M_{\mathrm{cold}}/M_*)_{1.5R_{\mathrm{opt}}}$', 'log'),\
                'cgf_gasSF':         ('$\log_{10}(M_{\mathrm{cold}}/M_*)_{\mathrm{SF-gas}}$', 'log'),\
                'bheff':             ('$\log_{10}(M_{\mathrm{BH}}/M_{\mathrm{200c}}$)', 'log'),\
                   'Tcons':             ('$\log_{10}(T_{\mathrm{cons}}$ [Gyr])', 'log'),\
                   'Tcons_1.5ropt':             ('$T_{\mathrm{cons,1.5}R_{\mathrm{opt}}}$ [Gyr]', 'lin'),\
                'Tcons_30kpc':             ('$T_{\mathrm{cons,30kpc}}$ [Gyr]', 'lin'),\
                'Tcons_HI_30kpc':             ('$T_{\mathrm{cons,H}_\mathrm{I}\mathrm{,30kpc}}$ [Gyr]', 'lin'),\
                'Tcons_HIH2_30kpc':             ('$T_{\mathrm{cons,H}_\mathrm{I}+\mathrm{H}_2\mathrm{,30kpc}}$ [Gyr]', 'lin'),\
                   'Tcons_gasSF':             ('$T_{\mathrm{cons,SF-gas}}$ [Gyr]', 'lin'),\
                'rhalf_mass':        ('$r_{\mathrm{1/2,DM halo}}$ [kpc]', 'lin'),\
                   'NFW_con':           ('c$_{\mathrm{NFW}}$', 'lin'),\
                   'r200c':             ('$r_{\mathrm{200c}}$ [kpc]', 'lin'),\
                   'rhalf_stars':       ('$R_{\mathrm{1/2,*}}$ [kpc]', 'lin'),\
                   'rhalf_stars_30kpc': ('$R_{\mathrm{1/2,*,30kpc}}$ [kpc]', 'lin'),\
                   'rhalf_stars_1.5ropt':       ('$R_{\mathrm{1/2,*,1.5}R_{\mathrm{opt}}}$ [kpc]', 'lin'),\
                   'rhalf_gas':         ('$R_{\mathrm{1/2,gas}}$ [kpc]', 'lin'),\
                   'reff_gas_1.5ropt':         ('$R_{\mathrm{eff,gas,1.5}R_{\mathrm{opt}}}$ [kpc]', 'lin'),\
                   'reff_gas_disk_1.5ropt':         ('$R_{\mathrm{eff,gas,disk,1.5}R_{\mathrm{opt}}}$ [kpc]', 'lin'),\
                   'ropt':              ('$R_{\mathrm{opt,80\% }M_*}$ [kpc]', 'lin'),\
                   'rVmax':              ('$R_{\mathrm{V}_{\mathrm{max}}}$ [kpc]', 'lin'),\
                       'r2voidCenter':              ('$r_{\mathrm{dist2centre}}$ [kpc]', 'lin'),\
                   'vbulge':            ('$V_{\mathrm{bulge}}$ [kms$^{-1}$]', 'lin'),\
                    'vdisk':            ('$V_{\mathrm{disk}}$ [kms$^{-1}$]', 'lin'),\
                   'vmax':              ('$V_{\mathrm{max}}$ [kms$^{-1}$]', 'lin'),\
                    'v200c':              ('$V_{\mathrm{200c}}$ [kms$^{-1}$]', 'lin'),\
                    'vdisp_30kpc_T19':   ('$V_{\mathrm{disp,30kpc,T19}}$ [kms$^{-1}$]', 'lin'),\
                   'Vrot2Vdisp_30kpc_T19':   ('$(V_{\mathrm{rot}}/V_{\mathrm{disp}})_{\mathrm{30kpc,T19}}$', 'lin'),\
                'jbar':              ('$\log_{10}(j_{\mathrm{bar}}$ [kpc kms$^{-1}$])', 'log'),\
                    'jbulge':              ('$\log_{10}(j_{\mathrm{bulge}}$ [kpc kms$^{-1}$])', 'log'),\
                        'jdisk':              ('$\log_{10}(j_{\mathrm{disk}}$ [kpc kms$^{-1}$])', 'log'),\
                   'angM_norm_1.5ropt_stars': ('$\log_{10}(j_{*,1.5R_{\mathrm{opt}}}$ [kpc kms$^{-1}$])', 'log'),\
                   'angM_norm_1.5ropt_gas': ('$\log_{10}(j_{\mathrm{gas,1.5}R_{\mathrm{opt}}}$ [kpc kms$^{-1}$])', 'log'),\
                   'angM_stars': ('$\log_{10}(j_{\mathrm{*}}$ [kpc kms$^{-1}$])', 'log'),\
                   'angM_SFgas': ('$\log_{10}(j_{\mathrm{SF-gas}}$ [kpc kms$^{-1}$])', 'log'),\
                       'angM_NSFgas': ('$\log_{10}(j_{\mathrm{NSF-gas}}$ [kpc kms$^{-1}$])', 'log'),\
                           'angM_disk': ('$\log_{10}(j_{*,\mathrm{disk}}$ [kpc kms$^{-1}$])', 'log'),\
                               'angM_spheroid': ('$\log_{10}(j_{*,\mathrm{bulge}}$ [kpc kms$^{-1}$])', 'log'),\
                       'DvT_stars_c0.4_1.5ropt':         ('$(M_{\mathrm{disk}}/M_*)_{\mathrm{c>0.4}}$', 'lin'),\
                   'DvT_stars_c0.5_1.5ropt':     ('$(M_{disk}/M_*)_{\mathrm{c>0.5}}$', 'lin'),\
                   'DvT_gas_c0.5_1.5ropt':           ('$(M_{\mathrm{disk,gas}}/M_{\mathrm{gas}})_{\mathrm{c>0.5}}$', 'lin'),\
                   'mAB_total_cut_g_i': ('($g$-$i$)', 'lin'),\
                   'mAB_total_cut_r_i': ('($r$-$i$)', 'lin'),\
                   'mAB_dA_total_cut_g_i': ('($g$-$i$)', 'lin'),\
                   'mAB_dA_total_cut_r_i': ('($r$-$i$)', 'lin'),\
                   'mAB_dA_total_cut_g_r': ('($g$-$r$)', 'lin'),\
                   'mAB_dA_total_cut_B_V': ('($B$-$V$)', 'lin'),\
                'MAB_dA_SDSS_g':     ('$M_{\mathrm{AB}_g}$', 'lin'),\
                   'MAB_dA_SDSS_r':     ('$M_{\mathrm{AB}_r}$', 'lin'),\
                   'MAB_dA_total_i':     ('$M_{\mathrm{AB}_i}$', 'lin'),\
                   'MAB_dA_Johnson_B':     ('$M_{\mathrm{AB}_B}$', 'lin'),\
                   'MAB_dA_Johnson_V':     ('$M_{\mathrm{AB}_V}$', 'lin'),\
                'lastMinorM':        ('$t_{\mathrm{LB,last mm}}$ [Gyr]', 'lin'),\
                   'lastMajorM':        ('$t_{\mathrm{LB,last Mm}}$ [Gyr]', 'lin'),\
                   'lastMerger':        ('$t_{\mathrm{LB,last M}}$ [Gyr]', 'lin'),\
                   'ratio_lastM':        ('ratio$_{\mathrm{last M}}$', 'lin'),\
                'timeLastIsolated':  ('$t_{\mathrm{last isolated}}$ [Gyr]', 'lin'),\
                   'Sigma_sfr_1.5ropt':  ('$\log_{10}(\Sigma_{\mathrm{SFR,1.5}R_{\mathrm{opt}}}$ [$M_{\odot}$ kpc$^{-2}$])', 'log'),\
                   'Sigma_sfr_30kpc':  ('$\log_{10}(\Sigma_{\mathrm{SFR,30kpc}}$ [$M_{\odot}$ kpc$^{-2}$])', 'log'),\
                'Sigma_stars_1.5ropt':  ('$\log_{10}(\Sigma_{\mathrm{*,1.5}R_{\mathrm{opt}}}$ [$M_{\odot}$ kpc$^{-2}$])', 'log'),\
                   'Sigma_stars_30kpc':  ('$\log_{10}(\Sigma_{\mathrm{*,30kpc}}$ [$M_{\odot}$ kpc$^{-2}$])', 'log'),\
                'Sigma_gas_1.5ropt':  ('$\log_{10}(\Sigma_{\mathrm{gas,1.5}R_{\mathrm{opt}}}$ [$M_{\odot}$ pc$^{-2}$])', 'log'),\
                   'Sigma_gas_reff_1.5ropt':  ('$\log_{10}(\Sigma_{\mathrm{gas,eff,1.5}R_{\mathrm{opt}}}$ [$M_{\odot}$ pc$^{-2}$])', 'log'),\
                'Sigma_gas_reff_disk_1.5ropt':  ('$\log_{10}(\Sigma_{\mathrm{gas,disk,eff,1.5}R_{\mathrm{opt}}}$ [$M_{\odot}$ pc$^{-2}$])', 'log'),\
                   'Sigma_HI_30kpc':  ('$\log_{10}(\Sigma_{\mathrm{H}_\mathrm{I}\mathrm{,30kpc}}$ [$M_{\odot}$ pc$^{-2}$])', 'log'),\
                   'Sigma_H2_30kpc':  ('$\log_{10}(\Sigma_{\mathrm{H}_2\mathrm{,30kpc}}$ [$M_{\odot}$ pc$^{-2}$])', 'log'),\
                   'Sigma_HIH2_30kpc':  ('$\log_{10}(\Sigma_{\mathrm{H}_\mathrm{I}+\mathrm{H}_2\mathrm{,30kpc}}$ [$M_{\odot}$ pc$^{-2}$])', 'log'),\
                   'sfe_gas_reff_disk_1.5ropt':  ('$\log_{10}$(SFE$_{\mathrm{gas,eff,disk,1.5}R_{\mathrm{opt}}}$ [yr$^{-1}$])', 'log'),\
                   'sfe_gas_reff_1.5ropt':  ('$\log_{10}$(SFE$_{\mathrm{gas,eff,1.5}R_{\mathrm{opt}}}$ [yr$^{-1}$])', 'log'),\
                'sfe_gas_1.5ropt':  ('$\log_{10}$(SFE$_{\mathrm{gas,1.5}R_{\mathrm{opt}}}$ [yr$^{-1}$])', 'log'),\
                   'sfe_HIH2_30kpc':  ('$\log_{10}$(SFE$_{\mathrm{H}_\mathrm{I}+\mathrm{H}_2\mathrm{,30kpc}}$ [yr$^{-1}$])', 'log'),\
                'SB_mu_eff_gas_disk_B': ('$\mu_{\mathrm{eff,gas,disk,1.5}R_{\mathrm{opt}},B}$', 'lin'),\
                   'SB_mu_opt_B': ('$\mu_{\mathrm{opt,}B}$', 'lin'),\
                   'SB_mu_1.5opt_B': ('$\mu_{1.5R_{\mathrm{opt}},B}$', 'lin'),\
                'SB_mu_eff_stars_1.5ropt_r': ('$\mu_{\mathrm{eff,*,1.5}R_{\mathrm{opt}},r}$', 'lin'),\
                   'SB_mu_eff_stars_30kpc_r': ('$\mu_{\mathrm{eff,*,30kpc,}r}$', 'lin'),\
                   'age_stars_rband_r502D': ('age$_{*,R_{\mathrm{1/2}},r}$ [Gyr]', 'lin'),\
                   'age_stars_rband_2r502D': ('age$_{*,2R_{\mathrm{1/2}},r}$ [Gyr]', 'lin'),\
                   'age_mean_stars': ('age$_{\mathrm{*,mean}}$ [Gyr]', 'lin'),\
                'ell_30kpc_T19': ('ell$_{\mathrm{30kpc,T19}}$', 'lin'),\
                   'z_mean_birth_stars':  ('$z_{\mathrm{mean,birth,*}}$', 'lin'),\
                'nr_sample_keys': ('$N_{\mathrm{samples}}$', 'lin'),\
                    't50_stars': ('$t_{50\%,*}$ [Gyr]', 'lin'),\
                        't70_stars': ('$t_{70\%,*}$ [Gyr]', 'lin'),\
                            't50_halo': ('$t_{50\%,halo}$ [Gyr]', 'lin'),\
                                't70_halo': ('$t_{70\%,*}$ [Gyr]', 'lin'),\
                                    't50_bh': ('$t_{50\%,bh}$ [Gyr]', 'lin'),\
                                        't70_bh': ('$t_{70\%,bh}$ [Gyr]', 'lin'),\
                    'delta_tform_stars': ('$\Delta_{t_{50\%}-t_{70\%},*}$ [Gyr]', 'lin'),\
                        'delta_tform_halo': ('$\Delta_{t_{50\%}-t_{70\%},halo}$ [Gyr]', 'lin'),\
                            'delta_tform_bh': ('$\Delta_{t_{50\%}-t_{70\%},bh}$ [Gyr]', 'lin'),\
                                'delta_t50_st2ha': ('$\Delta_{t_{50\%,halo}-t_{50\%},*}$ [Gyr]', 'lin'),\
                                    'delta_t50_bh2ha': ('$\Delta_{t_{50\%,halo}-t_{50\%},bh}$ [Gyr]', 'lin'),\
                                        'delta_age_stars_rband': ('$\Delta_{age_{R_{1/2}}-age_{2R_{1/2}},*,r}$ [Gyr]', 'lin'),\
                    }
                   
    return myprop_unit_map

def calc_residuals(data, count, redshift, sample, prop):
    #print 'count:', count, 'redshift:', redshift, 'sample:', sample, 'prop:', prop 
    from statsmodels import robust 
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
        
    myheader='\n('+str(i)+') z ('+str(i+1)+') N_M1found/N_M1 ('+str(i+2)+') N_M1found/N_M2 ('+str(i+3)+') N_M1/N_M1'+\
             '('+str(i+4)+') '+myprop_unit_map[prop]+' M1found ('+str(i+5)+') +dy ('+str(i+6)+') -dy ('+str(i+7)+') MAD ('+str(i+8)+') N_M1found [-] '+\
             '('+str(i+9)+') '+myprop_unit_map[prop]+' M1 ('+str(i+10)+') +dy ('+str(i+11)+') -dy ('+str(i+12)+') MAD ('+str(i+13)+') N_M1 [-] '+\
             '('+str(i+14)+') '+myprop_unit_map[prop]+' M2 ('+str(i+15)+') +dy ('+str(i+16)+') -dy ('+str(i+17)+') MAD ('+str(i+18)+') N_M2 [-] '+\
             '('+str(i+19)+') '+myprop_unit_map[prop]+' median(M1_found)/median(M1)-1 [-] '+\
             '('+str(i+20)+') '+myprop_unit_map[prop]+' median(M1_found)/median(M2)-1 [-] '+\
             '('+str(i+21)+') '+myprop_unit_map[prop]+' median(M1)/median(M2)-1 [-] ('+str(i+22)+') '+myprop_unit_map[prop]+' MAD(M1)/MAD(M2)-1'


    data[count, 0] = redshift
    data[count, 1] = data_M1_found.size/float(data_M1.size)
    data[count, 2] = data_M1_found.size/float(data_M2.size)
    data[count, 3] = data_M1.size/float(data_M2.size)
    data[count, 4] = np.nanmedian(data_M1_found[prop])
    data[count, 5] = np.nanmedian(data_M1_found[prop])-np.nanpercentile(data_M1_found[prop], 32)
    data[count, 6] = np.nanpercentile(data_M1_found[prop], 68)-np.nanmedian(data_M1_found[prop])
    data[count, 7] = robust.mad(data_M1_found[prop])
    data[count, 8] = data_M1_found.size
    data[count, 9] = np.nanmedian(data_M1[prop])
    data[count, 10] = np.nanmedian(data_M1[prop])-np.nanpercentile(data_M1[prop], 32)
    data[count, 11] = np.nanpercentile(data_M1[prop], 68)-np.nanmedian(data_M1[prop])
    data[count, 12] = robust.mad(data_M1[prop])    
    data[count, 13] = data_M1.size
    data[count, 14] = np.nanmedian(data_M2[prop])
    data[count, 15] = np.nanmedian(data_M2[prop])-np.nanpercentile(data_M2[prop], 32)
    data[count, 16] = np.nanpercentile(data_M2[prop], 68)-np.nanmedian(data_M2[prop])
    data[count, 17] = robust.mad(data_M2[prop])    
    data[count, 18] = data_M2.size
    data[count, 19] = np.nanmedian(data_M1_found[prop])/np.nanmedian(data_M1[prop])-1.0 
    data[count, 20] = np.nanmedian(data_M1_found[prop])/np.nanmedian(data_M2[prop])-1.0
    data[count, 21] = np.nanmedian(data_M1[prop])/np.nanmedian(data_M2[prop])-1.0
    data[count, 22] = robust.mad(data_M1[prop])/robust.mad(data_M2[prop])-1.0

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

def property_stats(data,props=[]):
    
    #data=data[np.where(data['orphan']==0)[:][0]]
    size=data.size
      
    pop=0
    if pop==1 or pop==2:
        data=data[np.where(data['pop']==pop)[:][0]]

      
    print 'pop:', pop, 'ngal:', data.size, '\b', "{0:.2f}".format(100.0/size*data.size), '% of all galaxies', size
#
    #for env, env_name in zip([-99,1,2,3],['all', 'sheet','filament','knot']):
    for env, env_name in zip([-99],['all']):            
        print '\n'
        if env_name!='all':
            sample=data[np.where(data['env_1024']==env)[:][0]]
        else:
            sample=data

        print 'env', env, 'name:', env_name, 'ngal:', sample.size, '\t', "{0:.0f}".format(100.0/data.size*sample.size),'% of all galaxies'
        print '////////////////////////////////////////////////'
    
        print 'prop\t\tmedian\t32th / 68th\n-------------------------------\n'
        #for prop in ['mhalo', 'mhalo_200c', 'NFW_con', 'mstar', 'zcold', 'sfr', 'ssfr', 'mcold', 'Mzgas', 'cgf', 'mbh','MAB_dA_total_g', 'MAB_dA_total_r', 'MAB_dA_total_i', 'mAB_dA_total_cut_r_i','mAB_dA_total_cut_g_r','mAB_dA_total_cut_g_i', 'mean_age_stars_disk', 'mean_age_stars_spheroid']:
        for prop in props:#, 'mstar', 'mcold', 'Mzgas', 'mbh', 'zcold', 'rhalf_mass', 'BvsT', 'sfr', 'ssfr', 'Tcons', 'cgf', 'fbar', 'zcold_zstar']:
            print prop,
            if prop.find('AB')==-1 and prop.find('zcold')==-1 and prop!='NFW_con' and prop.find('age')==-1 and prop!='fbar' and prop!='Tcons' and prop!='BvsT' and prop!='ecc' and prop!='a' and prop!='b' and prop!='rdisk' and prop!='rbulge':
                print '(log10)\t', "{0:.2f}".format(np.log10(np.nanmedian(sample[prop]))), '\t', "{0:.2f}".format(np.log10(np.nanmedian(sample[prop]))-np.log10(np.nanpercentile(sample[prop], 16))), '/', "{0:.2f}".format(np.log10(np.nanpercentile(sample[prop], 84))-np.log10(np.nanmedian(sample[prop])))
            else:
                print '\t', "{0:.2f}".format(np.nanmedian(sample[prop])), '\t', "{0:.2f}".format(np.nanmedian(sample[prop])-np.nanpercentile(sample[prop], 32)), '/', "{0:.2f}".format(np.nanpercentile(sample[prop], 68)-np.nanmedian(sample[prop]))


    print '--> DONE!\n'


def starforming_cut_Henriques20_A1(redshift):
    """calculate starforming cut after LGALAXIES SAM paper from Henriques+20 Appendix 1"""
    
    from cosmolopy import cparam, cd,cc
    fidcosmo = cparam.Planck(flat=True, extras=False)

    t_H0=(1.0/cd.hubble_z(redshift, **fidcosmo))/cc.yr_s    
    
    print 'z=', redshift, 't_H0:' , t_H0, '[yr]'
    print 'Henr+20:', format(10**(np.log10(2 * (1+redshift)**2/t_H0) - 1),'0.3e')
    print 'Franx+08:', format(0.3/t_H0,'0.3e')
    
    return 10**(np.log10(2 * (1+redshift)**2/t_H0) - 1)

def print_dt_to_screen(config_array,
                       catname):

    myhelper_dict={}
    for i in range(config_array[catname+'_id_col_array']['nr_entries']):           
        col_id = config_array[catname+'_id_col_array']['col_id'+str(i)]
        name = config_array[catname+'_id_col_array']['name'+str(i)]
        data_dtype = config_array[catname+'_id_col_array']['data_type'+str(i)]
        #print 'i:', 'col_id:', col_id, 'name:', name
        myhelper_dict.update({col_id: {'name': name, 'dtype': data_dtype}})

    print myhelper_dict
    
    for item in myhelper_dict.items():
        print '("'+item[1]['name']+'",'+str(item[1]['dtype'])+'),'
     
def print_props_table_format(data,
                             myprops=[],
                             envr_prop_name=None,
                             envr_props={},
                             percentiles=[25,75],
                             order_data_dict=False,
                             order_envr_props_dict=False,
                             order_myprops=False,
                             plot_median=False,
                             print_boot_err=False,
                             output_filename=None,
                             output_filename_key=''):

    import copy 
    import collections
    myFuncs = mF.MyFunctions() 
    
    myprop_unit_map  = get_prop_unit_map_latex()
    sample_name_dict = get_SFH_sample_name_dict()
    header_sample_names='\t&\t'
    header_stats_props=''
    stats_props=''
 
    if order_data_dict==True:
        data = collections.OrderedDict(sorted(data.items()))
    else:
        ordered_envr_props = envr_props
        
    if order_envr_props_dict==True:
        ordered_envr_props = collections.OrderedDict(sorted(envr_props.items()))
    else:
        ordered_envr_props = envr_props        

    if order_myprops==True:
        myprops = sorted(myprops, key=str.lower)
               
    #print sample_name_dict

    plot_stats_dict={}
    
    boot_mean= -99
    boot_std = -99

    print ordered_envr_props.items()
    print 'envr_prop_name:', envr_prop_name

    #loop over all properties
    for i,prop in enumerate(myprops):
        
        stats_props+=myprop_unit_map[prop][0]+'&\t'
        
        #loop for all data arrays in dict data={}
        for j, sample_name in enumerate(data.keys()):
            
            mydata=copy.deepcopy(data[sample_name])
               
            if i==0:
                header_sample_names+=sample_name_dict[sample_name]['name_Latex']+'\t'
                if j==len(data.keys())-1:
                    header_sample_names+='break'
                else:
                    header_sample_names+='&\t'
                    
            k=0
            #loop over target selections flag in envr_props={}
            for env, env_name in ordered_envr_props.items():              
  
                if env==-99:
                    data2analyse=mydata
                else:
                    data2analyse=mydata[np.where(mydata[envr_prop_name]==env)[:][0]]                  

                if i==0:
                    header_stats_props+=sample_name_dict[env_name]['name_Latex']+' $N_{\mathrm{gal}}$: '+str(data2analyse.size)+', '+str("{0:.1f}".format(100.0/mydata.size*data2analyse.size))+'\%\t'
 
                    if k==len(envr_props.keys())-1 and j==len(data.keys())-1 :
                        header_stats_props+='break'
                    else:
                        header_stats_props+='&\t'                       
                        
                data2analyse=data2analyse[np.where(data2analyse[prop]>-99)[:][0]]                

                print 'j:', j, 'i:', i, 'k:', k, 'envr:', env, env_name, data2analyse.size, 'sample:', sample_name, 'prop:', prop,
                if myprop_unit_map[prop][1]=='lin':

                    median  = np.nanmedian(data2analyse[prop])
                    lower   = np.nanmedian(data2analyse[prop])-np.nanpercentile(data2analyse[prop], percentiles[0])
                    upper   = np.nanpercentile(data2analyse[prop], percentiles[1])-np.nanmedian(data2analyse[prop])
                    
                    if print_boot_err!=False:
                        boot_mean, boot_std = myFuncs.bootstrap_error(data2analyse[prop])
                    
                else:
                    median  = np.log10(np.nanmedian(data2analyse[prop]))
                    lower   = np.log10(np.nanmedian(data2analyse[prop]))-np.log10(np.nanpercentile(data2analyse[prop], percentiles[0]))
                    upper   = np.log10(np.nanpercentile(data2analyse[prop], percentiles[1]))-np.log10(np.nanmedian(data2analyse[prop]))
                    
                    if print_boot_err!=False:
                        boot_mean, boot_std = myFuncs.bootstrap_error(np.log10(data2analyse[prop]))
                
                print 'median:', "{0:.2f}".format(median), '-/+', "{0:.2f}".format(lower), '/', "{0:.2f}".format(upper)
                if print_boot_err!=False:
                    stats_props+=str("{0:.2f}".format(median))+'$_{-'+str("{0:.2f}".format(lower))+'}^{+'+str("{0:.2f}".format(upper))+'}$\t & \t'+str("{0:.3f}".format(boot_std))+'\t'
                else:
                    stats_props+=str("{0:.2f}".format(median))+'$_{-'+str("{0:.2f}".format(lower))+'}^{+'+str("{0:.2f}".format(upper))+'}$\t'                    
                    
                if k==len(envr_props.keys())-1 and j==len(data.keys())-1:
                    stats_props+='break'
                else:
                    stats_props+='&\t'                  
                    
                    

                plot_stats_dict.update({prop+'_'+sample_name+env_name: {'sample_name_latex': sample_name_dict[sample_name]['name_Latex'],
                                                                        'env_name_latex': sample_name_dict[env_name]['name_Latex'],
                                                                       'median': median,
                                                                       'lower': lower,
                                                                       'upper': upper,
                                                                       'boot_err': boot_mean,
                                                                       'boot_lower': boot_std,
                                                                       'boot_upper': boot_std}})

                                                              
                k+=1
        
        if plot_median==True:
            plot_table_stats(plot_stats_dict,
                             prop,
                             sample_list=data.keys(),
                             envr_list=ordered_envr_props.values(),                         
                             output_filename=output_filename,
                             output_filename_key=output_filename_key)
            
    
    print '\n\n\\newpage\n\n\\begin{table*}\n\t\\begin{center}\n\t\setlength{\\tabcolsep}{3pt}\n\t\t\\begin{tabular}{l||c|c|c|c|c|c|c}'    
    print header_sample_names
    print header_stats_props
    print stats_props
    print' \t\t\t\hline\n\t\t\t (i)	\t & (ii)\t & (iii)\t & (iv)\t  & (v) \t & \t  (vi) \t  &  \t (vii) \t  & \t  (viii) \t break \n\t\t\end{tabular}\n\t\t\caption{\DS{Table just for checking, will not be published!}.}\label{tab:props_subsamples_stats_boxplots_test}\n\t\end{center}\n\end{table*}\n\n'

    print '\\begin{figure*}\n\t\centering\n\t\includegraphics[width=0.75\\textwidth,angle=0]{plots/boxplots/SFH_box_plot_stats_test.pdf}\n\t\caption{\DS{Table just for checking, will not be published!}.}\label{fig:box_plot_mhalo}\n\end{figure*}\n\n'
     
            
    exit()

def plot_table_stats(data_dict,
                     prop,
                     sample_list=[],
                     plot_median=False,
                     envr_list=[],
                     xaxis_length='large',
                     output_filename=None,
                     output_filename_key=''):
    
    #Load plot related packages
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.ticker import MultipleLocator
    from matplotlib.ticker import NullFormatter
    nullfmt = NullFormatter()
    mpl.style.use('classic')            
    mpl.mathtext.use_cm = False
    mpl.rcParams['mathtext.fontset'] = 'custom'
    mpl.rcParams['mathtext.tt'] = 'Typewriter'
    mpl.mathtext.fallback_to_cm = True
    
    nr_samples=len(sample_list)
    
    if nr_samples<=2:
        colours = ['k', 'r']
    elif nr_samples==3:
        colours = ['#AA4499', '#117733', '#696969' ]        
    elif nr_samples>3 and nr_samples<=6:        
        colours = ['#332288', '#CC6677', '#ffab00', '#117733', '#88CCEE', '#AA4499']
    else:
        colours = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', '#ffab00', '#CC6677', '#AA4499', '#7439A6', '#661100', '#333333', '#d2691e']
    frame_lw = 4.0

    myprop_unit_map  = get_prop_unit_map_latex()
    sample_name_dict = get_SFH_sample_name_dict()
    
    mycol=colours[0:nr_samples]
    
    print envr_list
    print sample_list
    print 'property:', prop, 'latex:',  myprop_unit_map[prop][0],


    if xaxis_length=='large':
          xsize=8.5
          ysize=3.0
          left=0.05
          right=0.66
          bottom=0.4
    else:
          xsize=5
          ysize=3
          left=0.05
          right=0.90
          bottom=0.25
          


    for i, env in enumerate(envr_list):
        
        fig = plt.figure(figsize=(xsize,ysize), dpi=150)
        ax0 = fig.add_subplot(111) 

        fig.subplots_adjust(hspace=0.01, wspace=0.01, left=left, bottom=bottom, right=right, top=0.98)    

        legend=[]
        axes=[]
        
        for k, sample in enumerate(sample_list):
            
            print 'i:', i, 'k:', k, env, sample,
            key1 = ax0            
            key2 = 'axis' + str(i)
    
            key2, = key1.plot([], color=mycol[k], linewidth=3.0, alpha=1, ls='-')
            
            if plot_median==True:
                lower = data_dict[prop+'_'+sample+env]['median']-data_dict[prop+'_'+sample+env]['lower']
                upper = data_dict[prop+'_'+sample+env]['median']+data_dict[prop+'_'+sample+env]['upper']
            else:
                #use bootstrap error
                lower = data_dict[prop+'_'+sample+env]['boot_err']-data_dict[prop+'_'+sample+env]['boot_lower']
                upper= data_dict[prop+'_'+sample+env]['boot_err']+data_dict[prop+'_'+sample+env]['boot_upper']
            
            if k==0:
                check_max=upper
                check_min=lower
            
            #print mycol[k], 'median:', "{0:.2f}".format(data_dict[prop+'_'+sample]['median']), '-/+', "{0:.2f}".format(lower), '/', "{0:.2f}".format(upper)    
    
            plt.axvline(data_dict[prop+'_'+sample+env]['median'], ymin=-15, ymax=15, color=mycol[k], ls='-', lw=6, zorder=0)
            plt.axvline(lower, ymin=-10, ymax=10, color=mycol[k], ls='--', lw=3, zorder=0)
            plt.axvline(upper, ymin=-10, ymax=10, color=mycol[k], ls='--', lw=4, zorder=0)
            
            if upper>check_max:   check_max=upper
            if lower<check_min:   check_min=lower        
                        
            axes.extend([key2])
            legend.extend([sample_name_dict[sample]['name_Latex']])
        
        fig.text(0.35, 0.075, myprop_unit_map[prop][0], ha='center', fontsize=30)
        plt.setp(ax0.get_yticklabels(), visible=False)
        
        fig.text(0.7, 0.85, sample_name_dict[env]['name_Latex'], ha='left', fontsize=32)   
    
        #print 'max-min:', check_max-check_min
        if check_max>=10.0 and (check_max-check_min)>5.0:
            float_format= 0
            x_minor     = 0.5
            ax0.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{0:.0f}".format(x)))
    
        if check_max>=10.0 and (check_max-check_min)>1.0:
            float_format= 1
            x_minor     = 0.2
            ax0.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{0:.1f}".format(x)))
        elif check_max>=10.0 and (check_max-check_min)>0.5:
            float_format= 1
            x_minor     = 0.1
            ax0.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{0:.1f}".format(x)))       
        else:
            float_format= 2
            x_minor     = 0.05
            ax0.get_xaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: "{0:.2f}".format(x))) 
    
        add_range_xaxis     = abs(check_max-check_min)/100.0*20.0
        xmin                = round(check_min-add_range_xaxis, float_format)
        xmax                = round(check_max+add_range_xaxis, float_format)
        
        if prop=='Tcons':
            xmin=0
            x_minor=5
        #zcolds mhalo props
        # xmin=12.75
        # xmax=14.25        
        # x_minor = 0.1
        # xticks=[13.0, 13.5, 14.0] 
        
        minor_x = MultipleLocator(x_minor)
        ax0.xaxis.set_minor_locator(minor_x)
        plt.xlim(xmin, xmax)
        plt.ylim(0,10)                         
   
        print 'xlim:', xmin, '/', xmax, 'minor_x:', x_minor, 'xticks:',   
        xticks=[]
        l=0
        while l<len(ax0.get_xticks()[:]):
            #print 'i:', i, ax0.get_xticks()[k]
            xticks+=[round(ax0.get_xticks()[:][l],float_format)]
            l+=2       
            
        print xticks,     
        ax0.set_xticks(xticks)
        ax0.set_xticklabels(xticks)
        
        plt.setp(ax0.get_xticklabels(), visible=True)
    
        ax0.tick_params(axis='x', which='major', top='on', bottom='on', pad=10, labelsize=28, length=12, width=frame_lw, direction='in', zorder=20) 
        ax0.tick_params(axis='x', which='minor', top='on', bottom='on', pad=10, labelsize=28, length=8, width=frame_lw, direction='in', zorder=20) 
               
        ax0.tick_params(axis='y', which='both', left='on', right='on', pad=10, labelsize=16, length=8, width=0.0, direction='in', zorder=20) 
     
        for axis in ['top','bottom','left','right']:
            ax0.spines[axis].set_linewidth(frame_lw)   
     
        fig.legend(axes,
                    legend,
                    loc=[0.72,0.58-nr_samples*0.1],
                    ncol=1,
                    fontsize=25,
                    borderaxespad=3,
                    facecolor='w',
                    numpoints=1, 
                    frameon=False,
                    labelspacing=0.3,
                    handlelength=1.0,
                    columnspacing=0.95,
                    handletextpad=0.5)
          
        #plt.savefig(mycomp+'anaconda/pro/myRun/plots/plotXY/plot_only.png', format='png', rasterized=True, dpi=100, bbox_inches='tight', pad_inches=0.05) 
        from matplotlib.backends.backend_pdf import PdfPages
        pp = PdfPages(output_filename+'stats_'+prop+output_filename_key+'_'+env+'.pdf')
        plt.savefig(pp, format='pdf', rasterized=True, dpi=50, pad_inches=0.05, transparent=True)
        pp.close()
        print ' --> TABLE STATS PLOT', prop, ' ... done!\n'

        plt.close(fig)

def generate_snapidzred_file(start_snapid,end_snapid=0):
    print 'generate sanpidzred file!'
    print 'start:', start_snapid, 'end:', end_snapid
    import h5py as hdf5
    myOutput = oD.OutputData(config=False)

#    for k in range(0,len(f.attrs.keys()),1):
#        print '\t', f.attrs.keys()[k], ':\t\t', f.attrs.values()[k][0]
    
    path='/data/256_hydro_Cholla50Mpc_halo_tests/'
    path='/data/groups/comp-astro/bruno/cosmo_sims/halo_tests/256_hydro_50Mpc/output_files/'    
    for snapid in range(start_snapid,end_snapid-1,-1):
        #print 'snapid:', snapid,
        filename=path+str(snapid)+'.h5.0'
        #print 'filename', filename
        f=hdf5.File(filename, "r")
    
        output_filename=mycomp+'anaconda/pro/data/Cholla256_50Mpc/snapidzred.txt'
        
        if snapid==start_snapid:
            myOutput.writeIntoFile(
                       output_filename,
                       [str(snapid)+' '+str(f.attrs['Current_z'][0])+' '+str(f.attrs['Current_a'][0])+' '+str(f.attrs['t'][0])],
                       myheader='256 hydro Cholla 50 h-1Mpc from "halo tests" series from Bruno\n(1) snapid (2) z (3) a (4) t',
                       append_mytext=False,
                       data_is_string=False,
                       data_format='%s')
        else:
            myOutput.writeIntoFile(
                       output_filename,
                       str(snapid)+' '+str(f.attrs['Current_z'][0])+' '+str(f.attrs['Current_a'][0])+' '+str(f.attrs['t'][0])+'\n',
                       append_mytext=True,
                       data_is_string=True,
                       data_format='%s')         
    
def write_snapidzred_as_bash(start_snapid,end_snapid=0):
    """Function writes the redshift and related properties as if statements in bash used in configuratoin files cat_config.sh
    
            input
        =========
            
        start_snapid:   integer number of the start snapshot id number
        end_snapid:     integer number of the end snapshot id number (usually id=0 for z=0 and a=1) --> end_snapid must be smaller than start_snapid
          
        output
        =========
        
        print on screen:    		else if ($item  =~ xx*) then
                            			set SAM_scale_factor = xx
                                    set SAM_redshift = xx
                                else if ($item  == xx) then
                            			set SAM_scale_factor = xx
                            			set SAM_redshift = xx
    
    """
 
    import h5py as hdf5
    
    print 'write redshift props as bash!\nSnapshot start:', start_snapid, 'end:', end_snapid
    
    path='/data/256_hydro_Cholla50Mpc_halo_tests/'
    path='/data/groups/comp-astro/bruno/cosmo_sims/halo_tests/256_hydro_50Mpc/output_files/' 
     
    for snapid in range(start_snapid,end_snapid-1,-1):
        #print 'snapid:', snapid,
        filename=path+str(snapid)+'.h5.0'
        #print 'filename', filename
        f=hdf5.File(filename, "r")
 
        print 'else if ($item  =~ '+str(snapid)+'*) then\n' \
            '\tset SAM_scale_factor = '+str(f.attrs['Current_a'][0])+'\n' \
			'\tset SAM_redshift = '+str(format(f.attrs['Current_z'][0],'0.2f'))+'\n' \
		'else if ($item  == '+str(snapid)+') then\n'\
			'\tset SAM_scale_factor = '+str(f.attrs['Current_a'][0])+'\n' \
			'\tset SAM_redshift = '+str(format(f.attrs['Current_z'][0],'0.2f'))+'\n'
  
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

def assign_progenitor_indices_old(data_RS_now, redshift_before):
    """Function assigns predestant indices to the descdent indices and the haloids from a list of halo properties from ROCKSTAR out_xxx.list halo catalog. xxx stands for the 
    the snapshot number. The function goes through every snapshot and connect the haloid and descdent index of the main progenitor (also most massive progenitor) of a halo with
    the predestant index a snapshot before. It further assigns a unique tree index (treeID -- an ascending number of the first occurence of the halo)
    
    input:
    ==========
        data_RS_now:     (array) structured array of the ROCKSTAR halo properties catalog of the current snapshot
        redshift_before: (float) redshift of one snapshot before, if the first snapshot is processed it will take this redshift
        
    output:
    ==========        
        data:   (array) catalog of main progenitor halos with assigned predestant (predIndex) at the current snapshot, unique first progenitor ID (firstProgenitorID), number of progenitors where the main progenitor
                is chosen from, difference in halo mass (delta_mhalo) and virial radius (delta_rvir) to the halo on the tree
                one snapshot before
    
    """
    
    print 'redshift_before:', redshift_before 
    
 
    #Read the data you want to compare with one snapshot before!
    import read_MDGal as rMD

    try:
        #If you want to already processed data meaning that the predIndex and firstProgenitorIDs are already assigned that that this
        data_RS_before= rMD.readMDGal(myfilename=mycomp+'anaconda/pro/data/ROCKSTAR_50Mpc/ROCKSTAR_50Mpc_z_'+str(redshift_before)+'_tarsel_tree.hdf5', mycol=[1,13,38,40,41,42,43,44,45,48,49,50,51,52,53,54,0,56])
        print 'data_RS_before.size:', data_RS_before['npros'].size
    except:
        #If you want to assign them again, take this one
        data_RS_before= rMD.readMDGal(myfilename=mycomp+'anaconda/pro/data/ROCKSTAR_50Mpc/ROCKSTAR_50Mpc_z_'+str(redshift_before)+'_tarsel.hdf5', mycol=[1,13,38,40,41,42,43,44,45,48,49,50,51,52,53])

        #Add three more columns to your structured array
        import numpy.lib.recfunctions as rcfuncs
        data_RS_before = rcfuncs.append_fields([data_RS_before], ['npros','firstProgenitorID'] , [np.zeros(data_RS_before.size,),np.zeros(data_RS_before.size,),np.zeros(data_RS_before.size,)], dtypes=['i4','i4','i8'], usemask=False)
        data_RS_before['npros']= -99
        data_RS_before['firstProgenitorID']= -99       
        print '-- set npros and orphan values at start redshift of tree! --> data_RS_before.size:', data_RS_before.size
#    for k in range(0,100,1):
 #       print 'k:', k, data_RS_now[np.where(data_RS_now['haloid']==k)[:][0]][['haloid','descIndex', 'mhalo','x_pos', 'y_pos', 'z_pos']], data_RS_before[np.where(data_RS_before['descIndex']==k)[:][0]][['haloid','descIndex', 'mhalo','x_pos', 'y_pos', 'z_pos']]

    data_RS_now['delta_mhalo']=0.0
    data_RS_now['delta_rvir']=0.0

    #Sort the data you want to assign a predIndex, the halo on top is the main progenitor (also the most massive one)    
    data_RS_before[::-1].sort(order=['descIndex','mhalo','haloid'], axis=0)
    print 'start values --> data_RS_before:\n', data_RS_before[['haloid','descIndex','firstProgenitorID', 'npros', 'mhalo', 'x_pos', 'y_pos', 'z_pos']]
    
    #Look for the unique IDs, the functions selects the ID which appreas at first which are also our main progenitors. Store in count how often a certain descIndex is occuring
    #that is the number of progenitors of the current halos. Assign the satellite status as orphan==0 which means it is a central halo (note that this does not matter right now for
    #this analysis)
    non_unique, index, count = np.unique(data_RS_before['descIndex'], return_index=True, return_counts=True)
    data_RS_before['npros'][index]=count
  
    print 'after finding uniques and setting npros, orphans and haloid --> data_RS_before:\n', data_RS_before[['haloid','descIndex','firstProgenitorID','orphan', 'npros', 'mhalo', 'x_pos', 'y_pos', 'z_pos']]     
    print 'start values --> data_RS_now:\n', data_RS_now[['haloid','descIndex','predIndex','firstProgenitorID','orphan', 'npros', 'mhalo', 'delta_mhalo', 'x_pos', 'y_pos', 'z_pos']]
    
    
    #Intersect now the haloids from this snapshot with the descdent idexes (descIndex) a snapshot before. Store in index_now and index_before where those halos
    #can be found in our structured array
    parentIndices, index_now, index_before = np.intersect1d(data_RS_now['haloid'], data_RS_before['descIndex'], return_indices=True)
    
    index_start=0
    index_end=len(index_now)
    print index_now[index_start:index_end], '\n', data_RS_now[['haloid','descIndex', 'firstProgenitorID','mhalo','x_pos', 'y_pos', 'z_pos']][index_now[index_start:index_end]],'\n', index_before[index_start:index_end], '\n', data_RS_before[['haloid','descIndex', 'mhalo','x_pos', 'y_pos', 'z_pos']][index_before[index_start:index_end]]


    #The predestant index (predIndex) is the haloid of the halo on the main progenitor one snapshot before. We set this haloid as predIndex in the current snapshot       
    data_RS_now['predIndex'][index_now]=data_RS_before['haloid'][index_before]
    
    #We calculate the difference in halo mass/virial radius of the same halo found now and one snaphsot before
    data_RS_now['delta_mhalo'][index_now]=data_RS_now['mhalo'][index_now]-data_RS_before['mhalo'][index_before]
    data_RS_now['delta_rvir'][index_now]=data_RS_now['rvir'][index_now]-data_RS_before['rvir'][index_before]
    
    #We transfer further data from the halo before to the current halos e.g. the firstProgenitorID which also serves as unique identification of the merger tree
    data_RS_now['npros'][index_now]=data_RS_before['npros'][index_before]
    data_RS_now['firstProgenitorID'][index_now]=data_RS_before['firstProgenitorID'][index_before]

    #and sort them    
    data_RS_now[::-1].sort(order=['predIndex','mhalo','haloid'], axis=0)
    data_RS_now[::-1].sort(order=['descIndex','mhalo','haloid'], axis=0)
    
    
    print 'after linking snapshot --> data_RS_before:\n', data_RS_before[['haloid','descIndex','firstProgenitorID','orphan', 'npros', 'mhalo', 'x_pos', 'y_pos', 'z_pos']]#[np.where(data_RS_before['descIndex']==0)[:][0]]    
    print 'data_RS_now:\n', data_RS_now[['haloid','descIndex','predIndex', 'firstProgenitorID', 'orphan', 'npros', 'mhalo', 'delta_mhalo', 'x_pos', 'y_pos', 'z_pos']]#[np.where(data_RS_now['predIndex']==0)[:][0]]


    #We set firstProgenitorID which is also the root of the merger tree (treeID) and a unique identification of the whole tree!
    #This ID is an ascendent count of all main progenitor trees in the data, if a new halo was find, then the firstProgenitorID is assigned at first -99
    #those IDs need to be assigned now. We find them as data_wo_root (halos in data without a root ID)      
    data_wo_root=data_RS_now[np.where(data_RS_now['firstProgenitorID']==-99)[:][0]]

    #We sort them by the haloid which is also ascending      
    data_wo_root.sort(order=['haloid'], axis=0)
    
    #We assign the firstProgenitorID of the newly found halos an ascending number count. That means that newly found halos with lower haloids get lower firstProgenitorIDs
    data_wo_root['firstProgenitorID']=range(data_RS_now.size-data_wo_root.size,data_RS_now.size,1)

    #The last step is to assign now the right halos their new firstProgenitorID (unique identification of the main progrenitor merger tree) by comparing the haloids       
    parentIndices, index_root, index_now = np.intersect1d(data_wo_root['haloid'], data_RS_now['haloid'], return_indices=True)
    data_RS_now['firstProgenitorID'][index_now]=data_wo_root['firstProgenitorID'][index_root]
    
    print 'data_RS_now:\n', data_RS_now[['haloid','descIndex','predIndex', 'firstProgenitorID', 'npros', 'mhalo', 'delta_mhalo', 'x_pos', 'y_pos', 'z_pos']]    

    return data_RS_now

def get_property_list(project_name, 
                      list_type='default'):
    """returns the list of properties processed and caculated within a certain project
        input:      
            project_name   e.g. 'SFH_corr3' the list of propertie in the star formation history project and the correction 3 run)
            list_type:     some properties are calcuated with error estimations e.g. within percetiles, then the the subsequent
                           2 column correspond to +dy and -dy in the percentiles defined at the caculation run see                     
                    
        output:
            list with property names --> []
            
    """
    
    if project_name=='SFH_corr3':
        """note that the entry 'cSFRD' is acutally only 'sum_sfr' while the normalisation with the volume needs to be done when plotted!"""
        if list_type!='default':          
            
            return ['BvT', 'g-i', 'r-i', 'bheff', 'Tcons',\
                       'mbar', 'jbar', 'jbulge', 'jdisk', 'cgf', 'fbar', 'rbulgevsrdisk', 'sfr', 'ssfr', 'mstar', 'mhalo', 'mbh', 'SHMR',\
                       'mcold', 'Mzgas', 'zcold', 'rdisk', 'rbulge', 'rhalf_mass', 'mean_age_stars_disk', 'mean_age_stars_spheroid',\
                       'vmax', 'vdisp', 'vdisk', 'vbulge']           
        else:
            return ['z', 't', 'n_cents', 'n_sats', 'n_orphans', 'density', 'cSFRD', 'sum_mstar']
        
    elif project_name=='SFH_old':
        """note that the entry 'cSFRD' is acutally the 'sum_sfr' and viceversa, therefore the normalisation with the volume needs to be done when plotted!"""
        return ['z',\
                'cSFRD', 'sfr2sfr_zstart', 'sfr', '-1sigsfr', '+1sigsfr',\
                'sum_ssfr', 'ssfr2ssfr_zstart', 'ssfr', '-1sigssfr', '+1sigssfr',\
                'g-i', '-1sigg-i', '+1sigg-i',\
                'sum_mstar', 'mstar2mstar_zstart', 'mstar', '-1sigmstar', '+1sigmstar',\
                'rdisk', '-1sigrdisk',\
                'rbulgevsrdisk', '-1sigrbulgevsrdisk', '+1sigrbulgevsrdisk',\
                't', 'n_cents', 'n_sats', 'n_orphans', 'density',\
                'r-i', '-1sigr-i', '+1sigr-i',\
                'mbh', '-1sigmbh', '+1sigmbh',\
                'rhalf_mass', '-1sigrhalf_mass', '+1sigrhalf_mass',\
                'sum_mcold', 'mcold2mcold_zstart', 'mcold', '-1sigmcold', '+1sigmcold',\
                'sum_Mzgas', 'Mzgas2Mzgas_zstart', 'Mzgas', '-1sigMzgas', '+1sigMzgas',\
                'zcold', '-1sigzcold', '+1sigzcold',\
                'sum_sfr',\
                'sum_mhalo', 'mhalo2mhalo_zstart', 'mhalo', '-1sigmhalo', '+1sigmhalo',\
                'SHMR', '-1sigSHMR', '+1sigSHMR',\
                'mean_age_stars_disk', '-1sigmean_age_stars_disk', '+1sigmean_age_stars_disk',\
                'mean_age_stars_spheroid', '-1sigmean_age_stars_spheroid', '+1sigmean_age_stars_spheroid',\
                'vmax', '-1sigvmax', '+1sigvmax',\
                'vdisp', '-1sigvdisp', '+1sigvdisp']
    else:
        print 'Project name NOT FOUND!'
        return False

def get_SFH_catalog_props(method=False):

    return {'Gal-dens':         {'ncols': 72, 'nSnaps': 55, 'output_name': 'Gal-dens', 'volume': (1000.0/0.6778)**3, 'legend_add_info': 'centrals',\
                                 'folder_name':                     'Gal-dens_main_cents_',\
                                 'folder_name_2':                   '',\
                                 'run_name':                        '',\
                                 'filename_part2':                  '/sfr2z_Galacticus_1Gpc_z_4.15_tarsel_SFH_down3_'+method+'_main_cents',\
                                 'filename_part2_to_z':             '/sfr2z_Galacticus_1Gpc_z_',\
                                 'filename_part2_from_z_to_main':   '_tarsel_SFH_down3_'+method+'_main_cents',\
                                 'clustering_endfix_xi':            '_test_centrals_0.5_200.0_60.0_xi_test'},\
            'Gal-dens-corr3':   {'ncols': 98, 'nSnaps': 55, 'output_name': 'Gal-dens', 'volume': (1000.0/0.6778)**3, 'legend_add_info': 'centrals',\
                                 'folder_name':                     'Gal-dens_main_cents_',\
                                 'folder_name_2':                   '/corr3',\
                                 'run_name':                        '_corr3',\
                                 'filename_part2':                  '/sfr2z_Galacticus_1Gpc_z_4.15_tarsel_SFH_down3_'+method+'_main_cents',\
                                 'filename_part2_to_z':             '/sfr2z_Galacticus_1Gpc_z_',\
                                 'filename_part2_from_z_to_main':   '_tarsel_SFH_down3_'+method+'_main_cents',\
                                 'clustering_endfix_xi':            '_centrals_0.5_200.0_60.0_xi_test'},\
            'Gal2-dens':       {'ncols': 72, 'nSnaps': 55, 'output_name': 'Gal2-dens', 'volume': (1000.0/0.6778)**3,\
                                 'folder_name':                     'Gal2-dens_main_cents_',\
                                 'filename_part2':                  '/sfr2z_Galacticus_1Gpc_run2_z_4.15_tarsel_SFH_down3_'+method+'_main_cents',\
                                 'filename_part2_to_z':             '/sfr2z_Galacticus_1Gpc_run2_z_',\
                                 'filename_part2_from_z_to_main':   '_tarsel_SFH_down3_'+method+'_main_cents'},\
            'Gal-SDSS':         {'ncols': 17, 'nSnaps': 55, 'output_name': 'Gal-SDSS', 'volume': (1000.0/0.6778)**3,\
                                 'folder_name':                     'Galacticus_SDSS',\
                                 'filename_part2':                  '/sfr2z_Galacticus_1Gpc_z_0.56_tarsel_CUT3_Contreras+13_mcold_SFH_method2_cents',\
                                 'filename_part2_to_z':             '/sfr2z_Galacticus_1Gpc_z_',\
                                 'filename_part2_from_z_to_main':   '_tarsel_CUT3_Contreras+13_mcold_SFH_method2_cents'},\
            'Gal400-dens':      {'ncols': 40, 'nSnaps': 55, 'output_name': 'Gal400-dens', 'volume': (400.0/0.6778)**3,\
                                 'folder_name':                     'Gal400-dens_main_cents_',\
                                 'filename_part2':                  '/sfr2z_Galacticus_400Mpc_z_4.15_tarsel_SFH_down3_'+method+'_main_cents',\
                                 'filename_part2_to_z':             '/sfr2z_Galacticus_400Mpc_z_',\
                                 'filename_part2_from_z_to_main':   '_tarsel_SFH_down3_'+method+'_main_cents'},\
            'Gal-dens-300':     {'ncols': 38, 'nSnaps': 55, 'output_name': 'Gal-r0001', 'volume': (1000.0/0.6778)**3,\
                                 'folder_name':                     'Gal_300_main_cents_',\
                                 'filename_part2':                  '/sfr2z_Galacticus_1Gpc_z_1.27_tarsel_SFH_300_main_cents',\
                                 'filename_part2_to_z':             '/sfr2z_Galacticus_1Gpc_z_',\
                                 'filename_part2_from_z_to_main':   '_tarsel_SFH_300_main_cents'}
          }
                
def get_SFH_sample_name_dict(catname='', orphan=''):
    
    return  {'EAGLE':       {'name_Latex': 'EAGLE', 'name_print_redshift': catname+'\n'+orphan+'\n$all$', 'catalog_name': 'EAGLE_100Mpc'},\
             'S':          {'name_Latex': 'skeleton (S)', 'name_print_redshift': 'EAGLE-S', 'catalog_name': 'EAGLE_100Mpc'},\
             'W':          {'name_Latex': 'walls (W)', 'name_print_redshift': 'EAGLE-W', 'catalog_name': 'EAGLE_100Mpc'},\
             'IV':          {'name_Latex': 'IV', 'name_print_redshift': 'EAGLE-IV', 'catalog_name': 'EAGLE_100Mpc'},\
             'OV':          {'name_Latex': 'OV', 'name_print_redshift': 'EAGLE-OV', 'catalog_name': 'EAGLE_100Mpc'},\
             'V':          {'name_Latex': 'voids (V)', 'name_print_redshift': 'EAGLE-V', 'catalog_name': 'EAGLE_100Mpc'},\
             'B':          {'name_Latex': 'LSB', 'name_print_redshift': 'EAGLE-B', 'catalog_name': 'EAGLE_100Mpc'},\
             'C':          {'name_Latex': 'HSB', 'name_print_redshift': 'EAGLE-C', 'catalog_name': 'EAGLE_100Mpc'},\
             'D':          {'name_Latex': 'LSFE', 'name_print_redshift': 'EAGLE-D', 'catalog_name': 'EAGLE_100Mpc'},\
             'E':          {'name_Latex': 'HSFE', 'name_print_redshift': 'EAGLE-E', 'catalog_name': 'EAGLE_100Mpc'},\
             'B-S':          {'name_Latex': 'LSB-S', 'name_print_redshift': 'EAGLE-B-S', 'catalog_name': 'EAGLE_100Mpc'},\
             'C-S':          {'name_Latex': 'HSB-S', 'name_print_redshift': 'EAGLE-C-S', 'catalog_name': 'EAGLE_100Mpc'},\
             'D-S':          {'name_Latex': 'LSFE-S', 'name_print_redshift': 'EAGLE-D-S', 'catalog_name': 'EAGLE_100Mpc'},\
             'E-S':          {'name_Latex': 'HSFE-S', 'name_print_redshift': 'EAGLE-E-S', 'catalog_name': 'EAGLE_100Mpc'},\
             'B-W':          {'name_Latex': 'LSB-W', 'name_print_redshift': 'EAGLE-B-W', 'catalog_name': 'EAGLE_100Mpc'},\
             'C-W':          {'name_Latex': 'HSB-W', 'name_print_redshift': 'EAGLE-C-W', 'catalog_name': 'EAGLE_100Mpc'},\
             'D-W':          {'name_Latex': 'LSFE-W', 'name_print_redshift': 'EAGLE-D-W', 'catalog_name': 'EAGLE_100Mpc'},\
             'E-W':          {'name_Latex': 'HSFE-W', 'name_print_redshift': 'EAGLE-E-W', 'catalog_name': 'EAGLE_100Mpc'},\
             'B-V':          {'name_Latex': 'LSB-V', 'name_print_redshift': 'EAGLE-B-V', 'catalog_name': 'EAGLE_100Mpc'},\
             'C-V':          {'name_Latex': 'HSB-V', 'name_print_redshift': 'EAGLE-C-V', 'catalog_name': 'EAGLE_100Mpc'},\
             'D-V':          {'name_Latex': 'LSFE-V', 'name_print_redshift': 'EAGLE-D-V', 'catalog_name': 'EAGLE_100Mpc'},\
             'E-V':          {'name_Latex': 'HSFE-V', 'name_print_redshift': 'EAGLE-E-V', 'catalog_name': 'EAGLE_100Mpc'},\
             'BD':          {'name_Latex': 'LSB+LSFE', 'name_print_redshift': 'EAGLE-B+D', 'catalog_name': 'EAGLE_100Mpc'},\
             'BE':          {'name_Latex': 'LSB+HSFE', 'name_print_redshift': 'EAGLE-B+E', 'catalog_name': 'EAGLE_100Mpc'},\
             'CD':          {'name_Latex': 'HSB+LSFE', 'name_print_redshift': 'EAGLE-C+D', 'catalog_name': 'EAGLE_100Mpc'},\
             'CE':          {'name_Latex': 'HSB+HSFE', 'name_print_redshift': 'EAGLE-C+E', 'catalog_name': 'EAGLE_100Mpc'},\
            'BD-S':          {'name_Latex': 'LSB+LSFE-S', 'name_print_redshift': 'EAGLE-B+D-S', 'catalog_name': 'EAGLE_100Mpc'},\
            'BD-W':          {'name_Latex': 'LSB+LSFE-W', 'name_print_redshift': 'EAGLE-B+D-W', 'catalog_name': 'EAGLE_100Mpc'},\
            'BD-V':          {'name_Latex': 'LSB+LSFE-V', 'name_print_redshift': 'EAGLE-B+D-V', 'catalog_name': 'EAGLE_100Mpc'},\
            'BE-S':          {'name_Latex': 'LSB+HSFE-S', 'name_print_redshift': 'EAGLE-B+E-S', 'catalog_name': 'EAGLE_100Mpc'},\
            'BE-W':          {'name_Latex': 'LSB+HSFE-W', 'name_print_redshift': 'EAGLE-B+E-W', 'catalog_name': 'EAGLE_100Mpc'},\
            'BE-V':          {'name_Latex': 'LSB+HSFE-V', 'name_print_redshift': 'EAGLE-B+E-V', 'catalog_name': 'EAGLE_100Mpc'},\
            'CD-S':          {'name_Latex': 'HSB+LSFE-S', 'name_print_redshift': 'EAGLE-C+D-S', 'catalog_name': 'EAGLE_100Mpc'},\
            'CD-W':          {'name_Latex': 'HSB+LSFE-W', 'name_print_redshift': 'EAGLE-C+D-W', 'catalog_name': 'EAGLE_100Mpc'},\
            'CD-V':          {'name_Latex': 'HSB+LSFE-V', 'name_print_redshift': 'EAGLE-C+D-V', 'catalog_name': 'EAGLE_100Mpc'},\
            'CE-S':          {'name_Latex': 'HSB+HSFE-S', 'name_print_redshift': 'EAGLE-C+E-S', 'catalog_name': 'EAGLE_100Mpc'},\
            'CE-W':          {'name_Latex': 'HSB+HSFE-W', 'name_print_redshift': 'EAGLE-C+E-W', 'catalog_name': 'EAGLE_100Mpc'},\
            'CE-V':          {'name_Latex': 'HSB+HSFE-V', 'name_print_redshift': 'EAGLE-C+E-V', 'catalog_name': 'EAGLE_100Mpc'},\
            'LSB':          {'name_Latex': 'LSB', 'name_print_redshift': 'LSB', 'catalog_name': 'EAGLE_100Mpc'},\
            'HSB':          {'name_Latex': 'HSB', 'name_print_redshift': 'HSB', 'catalog_name': 'EAGLE_100Mpc'},\
            'low':          {'name_Latex': 'low', 'name_print_redshift': 'low', 'catalog_name': 'EAGLE_100Mpc'},\
            'inter':          {'name_Latex': 'inter', 'name_print_redshift': 'inter', 'catalog_name': 'EAGLE_100Mpc'},\
            'high':          {'name_Latex': 'high', 'name_print_redshift': 'high', 'catalog_name': 'EAGLE_100Mpc'},\
            'LSB-low':          {'name_Latex': 'LSB-low', 'name_print_redshift': 'LSB-low', 'catalog_name': 'EAGLE_100Mpc'},\
            'LSB-inter':          {'name_Latex': 'LSB-inter', 'name_print_redshift': 'LSB-inter', 'catalog_name': 'EAGLE_100Mpc'},\
            'LSB-high':          {'name_Latex': 'LSB-high', 'name_print_redshift': 'LSB-high', 'catalog_name': 'EAGLE_100Mpc'},\
            'HSB-low':          {'name_Latex': 'HSB-low', 'name_print_redshift': 'HSB-low', 'catalog_name': 'EAGLE_100Mpc'},\
            'HSB-inter':          {'name_Latex': 'HSB-inter', 'name_print_redshift': 'HSB-inter', 'catalog_name': 'EAGLE_100Mpc'},\
            'HSB-high':          {'name_Latex': 'HSB-high', 'name_print_redshift': 'HSB-high', 'catalog_name': 'EAGLE_100Mpc'},\
            'all':          {'name_Latex': 'Gal-dens', 'name_print_redshift': catname+'\n'+orphan+'\n$all$', 'catalog_name': 'Gal-dens-corr3'},\
            '':             {'name_Latex': 'Gal-dens', 'name_print_redshift': catname+'\n'+orphan+'\n$all$', 'catalog_name': 'Gal-dens-corr3'},\
             'Gal-dens':    {'name_Latex': 'Gal-dens', 'name_print_redshift': catname+'\n'+orphan+'\nGal-dens', 'catalog_name': 'Gal-dens-corr3'},\
            'low':          {'name_Latex': 'low-$SFR$', 'name_print_redshift': catname+'\n'+orphan+'\n$low$-$SFR$','catalog_name': 'Gal-dens'},\
            'high':         {'name_Latex': 'high-$SFR$', 'name_print_redshift': catname+'\n'+orphan+'\n$high$-$SFR$','catalog_name': 'Gal-dens'},\
            'passive':      {'name_Latex': '$passive$', 'name_print_redshift': catname+'\n'+orphan+'\n$passive$','catalog_name': 'Gal-dens'},\
            'active':       {'name_Latex': '$active$', 'name_print_redshift': catname+'\n'+orphan+'\n$active$','catalog_name': 'Gal-dens'},\
            'red':          {'name_Latex': '$red$', 'name_print_redshift': catname+'\n'+orphan+'\n$red$','catalog_name': 'Gal-dens'},\
            'blue':         {'name_Latex': '$blue$', 'name_print_redshift': catname+'\n'+orphan+'\n$blue$','catalog_name': 'Gal-dens'},\
            'redmstar11':   {'name_print_redshift': catname+'\n'+orphan+'\n'+'$red$ + $M_*\sim11$','catalog_name': 'Gal-dens'},\
            'redmhalo13':   {'name_print_redshift': catname+'\n'+orphan+'\n'+'$red$ + $M_{vir}\sim13$','catalog_name': 'Gal-dens'},\
            'lowmstar11':   {'name_print_redshift': catname+'\n'+orphan+'\n'+'$low$-$SFR$ + $M_*\sim11$','catalog_name': 'Gal-dens'},\
            'lowmhalo13':   {'name_print_redshift': catname+'\n'+orphan+'\n'+'$low$-$SFR$ + $M_{vir}\sim13$','catalog_name': 'Gal-dens'},\
            'mhalo12st':    {'name_print_redshift': catname+'\n'+orphan+'\n'+'$M_{vir}<12$','catalog_name': 'Gal-dens'},\
            'mhalo12gt':    {'name_print_redshift': catname+'\n'+orphan+'\n'+'$M_{vir}>12$','catalog_name': 'Gal-dens'},\
            'mstar10st':    {'name_print_redshift': catname+'\n'+orphan+'\n'+'$M_*<11$','catalog_name': 'Gal-dens'},\
            'mstar11gt':    {'name_print_redshift': catname+'\n'+orphan+'\n'+'$M_{*}>11$','catalog_name': 'Gal-dens'},\
            'zcold9st':     {'name_print_redshift': catname+'\n'+orphan+'\n'+'$Z_{\mathrm{cold}}<9$','catalog_name': 'Gal-dens'},\
            'zcold9gt':     {'name_print_redshift': catname+'\n'+orphan+'\n'+'$Z_{\mathrm{cold}}>9$','catalog_name': 'Gal-dens'},\
            'low-zcold':    {'name_Latex': 'low-$Z_{\mathrm{cold}}$', 'name_print_redshift': catname+'\n'+orphan+'\n'+'$low$-$Z_{\mathrm{cold}}$','catalog_name': 'Gal-dens-corr3'},\
            'high-zcold':   {'name_Latex': 'high-$Z_{\mathrm{cold}}$', 'name_print_redshift': catname+'\n'+orphan+'\n'+'$high$-$Z_{\mathrm{cold}}$','catalog_name': 'Gal-dens-corr3'},\
            'low-zcold-red': {'name_Latex': 'low-$Z_{\mathrm{cold}}^{red}$', 'name_print_redshift': catname+'\n'+orphan+'\n'+'$low$-$Z_{\mathrm{cold}}$','catalog_name': 'Gal-dens-corr3'},\
            'high-zcold-red': {'name_Latex': 'high-$Z_{\mathrm{cold}}^{red}$', 'name_print_redshift': catname+'\n'+orphan+'\n'+'$high$-$Z_{\mathrm{cold}}$','catalog_name': 'Gal-dens-corr3'},\
            'lowZcold-highMstar': {'name_print_redshift': catname+'\n'+orphan+'\n'+'$low$-$Z_{\mathrm{cold}}$','catalog_name': 'Gal-dens'},\
            'highZcold-highMstar': {'name_print_redshift': catname+'\n'+orphan+'\n'+'$high$-$Z_{\mathrm{cold}}$', 'catalog_name': 'Gal-dens'},\
            'knots':             {'name_Latex': 'knots', 'name_print_redshift': catname+'\n'+orphan+'\nknots', 'catalog_name': 'Gal-dens-corr3'},\
            'filaments':         {'name_Latex': 'filaments', 'name_print_redshift': catname+'\n'+orphan+'\nfilaments', 'catalog_name': 'Gal-dens-corr3'},\
            }  

def SFH_prop_col_map(catname):
    
    if catname=='Gal-dens-corr3':
        
        return {'': '',\
                'sumSFR': 6,\
                '_BvT':      8,\
                '_sfr':   44,\
                '_ssfr':  47,\
                '_g-i':   11,\
                '_mstar': 50,\
                '_r-i':   14,\
                '_mbh':   56,\
                '_mbar':    23,\
                '_jbar':    26,\
                '_jbulge':    29,\
                '_jdisk':    32,\
                '_rhalfmass': 77,\
                '_mcold': 62,\
                '_Mzgas': 65,\
                '_cgf':   35,\
                '_fbar':  38,\
                '_zcold': 68,\
                '_cSFRD': 6,\
                '_mhalo': 53,\
                '_SHMR':  59,\
                '_mean_age_stars_disk': 80,\
                '_mean_age_stars_spheroid': 83,\
                '_vmax':  68,\
                '_vdisp': 89,\
                '_vdisk': 92,\
                '_vbulge': 95,\
                '_Tcons': 20,\
                '_bheff': 17\
                }
        
    else:
        return {'': '',\
                      'sumSFR': 51,\
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
                      '_cgf':   40,\
                      '_zcold': 48,\
                      '_cSFRD': 2,\
                      '_mhalo': 54,\
                      '_SHMR':  57,\
                      '_mean_age_stars_disk': 60,\
                      '_mean_age_stars_spheroid': 63,\
                      '_vmax':  66,\
                      '_vdisp': 69,\
                      '_Tcons': 40\
                      }

def get_sample_code_space(sample):
    
    if sample!='':
        return '_'
    else:
        return ''

def get_norm_y(prop,
               method,
               catalog):

    if prop=='_cSFRD':
         cat_props = get_SFH_catalog_props(method)
         return cat_props[catalog]['volume']                                 
    else:
         return 1.0

def get_redshift_histos(catname,
                        histo_key):
       
    if catname.find('400')!=-1:
        if histo_key=='histo1':
            redshifts={'0': 0.55, '1': 0.74, '2': 0.84, '3': 0.89, '4': 1.0, '5': 1.26}
        elif histo_key=='histo2':
            redshifts={'0': 1.36, '1': 1.56, '2': 2.02, '3': 2.56, '4': 3.03, '5': 4.15}
        else:
            redshifts={'0': 0.55, '1': 0.74, '2': 0.89, '3': 2.02, '4': 3.03, '5': 4.15} 
        
    else:    
        if histo_key=='histo1':
            redshifts={'0': 0.56, '1': 0.74, '2': 0.82, '3': 0.9, '4': 0.99, '5': 1.27}
        elif histo_key=='histo2':
            redshifts={'0': 1.37, '1': 1.54, '2': 2.03, '3': 2.53, '4': 3.04, '5': 4.15}
        else:
            redshifts={'0': 0.56, '1': 0.74, '2': 0.99, '3': 2.03, '4': 3.04, '5': 4.15}
            
    return redshifts

def set_SFH_print_redshift(key,
                           catname,
                           orphan,
                           prop,
                           method,
                           myenv=''):

    if myenv.find('k')!=-1:
        return '\nknots\n'
    elif myenv.find('f')!=-1:
        return '\nfilaments\n'
    
    elif key.find('one')!=-1 or key.find('histo')!=-1:
        mysample_name = get_SFH_sample_name_dict(catname=catname, orphan=orpham)                           

        return mysample_name[sample]                                                                                      
    elif key=='gr':
        prop_unit_map=get_prop_unit_map_latex()
        return catname+'\n'+prop_unit_map[prop]+'\n' 
               
    elif key.find('histo')!=-1:                   
        if sample.find('zcold')==-1 and (sample.startswith('low') or sample.startswith('high')):
            return catname+'\n'+orphan+'\n'+sample+' SFR'
        elif sample=='':
            return catname+'\n'+orphan+'\nall'                             
        else:
            return catname+'\n'+orphan+'\n'+sample
    else:
        return catname+'\n'+method+'\n' 
    
def dt_SFH_ASCII(catalog):
    """get dt for SFH sample"""
    
    property_dict       = get_property_dict()
    
    #print property_dict['mstar']
    #print property_dict['mstar']['dtype']
    
    if catalog=='Gal-dens-corr3':
    
        property_list       = get_property_list('SFH_corr3')
        
        property_list_perc  = get_property_list('SFH_corr3', list_type='percentiles')

        dt = [(item, property_dict[item]['dtype']) for item in property_list]
        
        for item in property_list_perc:
            dt.extend([(item, property_dict[item]['dtype']), ('+1sig'+item, property_dict[item]['dtype']), ('-1sig'+item, property_dict[item]['dtype'])])
        
    else:

        property_list       = get_property_list('SFH_old')
        
        dt=[]
        for item in property_list:
            #print 'item:', item
            try:
                dt.extend([(item, property_dict[item]['dtype'])])
            except:
                dt.extend([(item, property_dict[item[5::]]['dtype'])])
    
    return dt  

def dt_MDPL_Galacticus_CMASS_down_sample3_original():
    """get dt for Gal-dens file more properties
    usual filename: Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3.txt
    
    # (1) haloid [idnumber] (2) hostid [idnumber] (3) mstar [Msun] (4) mhalo [Msun] (5) orphan [0=cent,1=sat,2=orph] (6) mcold [Msun] (7) Mzgas [Msun] (8) Mzstar [Msun] (9) mbh [Msun] (10) mstar_spheroid [Msun] (11) mstar_disk [Msun] (12) mcold_disk [Msun] (13) mcold_spheroid [Msun] (14) mhot [Msun] (15) sfr_spheroid [Msunyr-1] (16) sfr_disk [Msunyr-1] (17) sfr [Msunyr-1] (18) Mzgas_spheroid [Msun] (19) Mzgas_disk [Msun] (20) Mzstar_spheroid [Msun] (21) Mzstar_disk [Msun] (22) Mzhot_halo [Msun] (23) x_pos [comvMpc] (24) y_pos [comvMpc] (25) z_pos [comvMpc] (26) x_vel [kms-1] (27) y_vel [kms-1] (28) z_vel [kms-1] (29) L_SDSS_dA_total_u [4.4659e13WHz-1] (30) L_SDSS_dA_total_g [4.4659e13WHz-1] (31) L_SDSS_dA_total_r [4.4659e13WHz-1] (32) L_SDSS_dA_total_i [4.4659e13WHz-1] (33) L_SDSS_dA_total_z [4.4659e13WHz-1] (34) rdisk [Mpc] (35) rbulge [Mpc] (36) rhalf_mass [Mpc] (37) spinParameter [-] (38) MAB_dA_total_u [-] (39) MAB_dA_total_g [-] (40) MAB_dA_total_r [-] (41) MAB_dA_total_i [-] (42) MAB_dA_total_z [-] (43) mAB_dA_total_u [-] (44) mAB_dA_total_g [-] (45) mAB_dA_total_r [-] (46) mAB_dA_total_i [-] (47) mAB_dA_total_z [-] (48) mAB_dA_total_cut_r_i [-] (49) mAB_dA_total_cut_dmesa [-] (50) mAB_dA_total_cut_i_lt_dmesa [0=False,1=True] (51) mhalo_sat [Msun] (52) nodeIndex [idnumber] (53) parentIndex [idnumber] (54) satelliteNodeIndex [idnumber] (55) satelliteIndex [idnumber] (56) siblingIndex [idnumber] (57) satelliteMergeTime [Gyr] (58) isolated [0=non-iso,1=iso] (59) timeLastIsolated [Gyr] (60) mhalo_200c [Msun] (61) zcold [-] (62) cgf [-] (63) mbasic [Msun] (64) mbasic_200c [Msun] (65) NFW_con [-] (66) ssfr [yr-1] (67) weight_tot [-] (68) env_512 [0=vo,1=sh,2=fi,3=no] (69) env_1024 [0=vo,1=sh,2=fi,3=no] (70) pop [-]     
    """
    property_list = ['haloid', 'hostid', 'mstar', 'mhalo', 'orphan', 'mcold', 'Mzgas', 'Mzstar', 'mbh', 'mstar_spheroid', 'mstar_disk',\
                     'mcold_disk', 'mcold_spheroid', 'mhot', 'sfr_spheroid', 'sfr_disk', 'sfr', 'Mzgas_spheroid', 'Mzgas_disk',\
                     'Mzstar_spheroid', 'Mzstar_disk', 'Mzhot_halo', 'x_pos', 'y_pos', 'z_pos', 'x_vel', 'y_vel', 'z_vel',\
                     'L_SDSS_dA_total_u', 'L_SDSS_dA_total_g', 'L_SDSS_dA_total_r', 'L_SDSS_dA_total_i', 'L_SDSS_dA_total_z',\
                     'rdisk', 'rbulge', 'rhalf_mass', 'spinParameter',\
                     'MAB_dA_total_u', 'MAB_dA_total_g', 'MAB_dA_total_r', 'MAB_dA_total_i', 'MAB_dA_total_z',\
                     'mAB_dA_total_u', 'mAB_dA_total_g', 'mAB_dA_total_r', 'mAB_dA_total_i', 'mAB_dA_total_z',\
                     'mAB_dA_total_cut_r_i', 'mAB_dA_total_cut_dmesa', 'mAB_dA_total_cut_i_lt_dmesa',\
                     'mhalo_sat', 'nodeIndex', 'parentIndex', 'satelliteNodeIndex', 'satelliteIndex', 'siblingIndex', 'satelliteMergeTime',\
                     'isolated', 'timeLastIsolated',\
                     'mhalo_200c', 'zcold', 'cgf', 'mbasic', 'mbasic_200c', 'NFW_con', 'ssfr', 'weight_tot', 'env_512', 'env_1024', 'pop']
        
    property_dict     = get_property_dict()
    
    dt=[]
    for item in property_list:
        dt.extend([(item, property_dict[item]['dtype'])])
                            
    return dt   


def dt_MDPL_Galacticus_CMASS():
    """get dt for Gal-dens file more properties
    usual filename: Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3_more_props_test.txt
    
    # (1) haloid [idnumber] (2) hostid [idnumber] (3) mstar [Msun] (4) mhalo [Msun] (5) orphan [0=cent,1=sat,2=orph] (6) mcold [Msun] (7) Mzgas [Msun] (8) Mzstar [Msun] (9) mbh [Msun] (10) mstar_spheroid [Msun] (11) mstar_disk [Msun] (12) mcold_disk [Msun] (13) mcold_spheroid [Msun] (14) mhot [Msun] (15) sfr_spheroid [Msunyr-1] (16) sfr_disk [Msunyr-1] (17) sfr [Msunyr-1] (18) Mzgas_spheroid [Msun] (19) Mzgas_disk [Msun] (20) Mzstar_spheroid [Msun] (21) Mzstar_disk [Msun] (22) Mzhot_halo [Msun] (23) x_pos [comvMpc] (24) y_pos [comvMpc] (25) z_pos [comvMpc] (26) x_vel [kms-1] (27) y_vel [kms-1] (28) z_vel [kms-1] (29) L_SDSS_dA_total_u [4.4659e13WHz-1] (30) L_SDSS_dA_total_g [4.4659e13WHz-1] (31) L_SDSS_dA_total_r [4.4659e13WHz-1] (32) L_SDSS_dA_total_i [4.4659e13WHz-1] (33) L_SDSS_dA_total_z [4.4659e13WHz-1] (34) rdisk [Mpc] (35) rbulge [Mpc] (36) rhalf_mass [Mpc] (37) spinParameter [-] (38) MAB_dA_total_u [-] (39) MAB_dA_total_g [-] (40) MAB_dA_total_r [-] (41) MAB_dA_total_i [-] (42) MAB_dA_total_z [-] (43) mAB_dA_total_u [-] (44) mAB_dA_total_g [-] (45) mAB_dA_total_r [-] (46) mAB_dA_total_i [-] (47) mAB_dA_total_z [-] (48) mAB_dA_total_cut_r_i [-] (49) mAB_dA_total_cut_dmesa [-] (50) mAB_dA_total_cut_i_lt_dmesa [0=False,1=True] (51) mhalo_sat [Msun] (52) nodeIndex [idnumber] (53) parentIndex [idnumber] (54) satelliteNodeIndex [idnumber] (55) satelliteIndex [idnumber] (56) siblingIndex [idnumber] (57) satelliteMergeTime [Gyr] (58) isolated [0=non-iso,1=iso] (59) timeLastIsolated [Gyr] (60) mhalo_200c [Msun] (61) zcold [-] (62) cgf [-] (63) mbasic [Msun] (64) mbasic_200c [Msun] (65) NFW_con [-] (66) ssfr [yr-1] (67) weight_tot [-] (68) env_512 [0=vo,1=sh,2=fi,3=no] (69) env_1024 [0=vo,1=sh,2=fi,3=no] (70) pop [-] (71) mAB_dA_total_cut_g_r [-] (72) mAB_dA_total_cut_g_i [-] (73) Tcons [Gyr] (74) fbar [-] (75) BvsT [-] (76) rbulgevsrdisk [-] (77) MzgasvsMzstar [-] (78) zstar [-] (79) zcold_zstar [-]     
    
    """
    property_list = ['haloid', 'hostid', 'mstar', 'mhalo', 'orphan', 'mcold', 'Mzgas', 'Mzstar', 'mbh', 'mstar_spheroid', 'mstar_disk',\
                     'mcold_disk', 'mcold_spheroid', 'mhot', 'sfr_spheroid', 'sfr_disk', 'sfr', 'Mzgas_spheroid',\
                     'Mzgas_disk', 'Mzstar_spheroid', 'Mzstar_disk', 'Mzhot_halo',\
                     'x_pos', 'y_pos', 'z_pos', 'x_vel', 'y_vel', 'z_vel',\
                     'L_SDSS_dA_total_u', 'L_SDSS_dA_total_g', 'L_SDSS_dA_total_r', 'L_SDSS_dA_total_i', 'L_SDSS_dA_total_z',\
                     'rdisk', 'rbulge', 'rhalf_mass', 'spinParameter',\
                     'MAB_dA_total_u', 'MAB_dA_total_g', 'MAB_dA_total_r', 'MAB_dA_total_i', 'MAB_dA_total_z',\
                     'mAB_dA_total_u', 'mAB_dA_total_g', 'mAB_dA_total_r', 'mAB_dA_total_i', 'mAB_dA_total_z',\
                     'mAB_dA_total_cut_r_i', 'mAB_dA_total_cut_dmesa', 'mAB_dA_total_cut_i_lt_dmesa',\
                     'mhalo_sat', 'nodeIndex', 'parentIndex', 'satelliteNodeIndex', 'satelliteIndex', 'siblingIndex', 'satelliteMergeTime', 'isolated', 'timeLastIsolated',
                     'mhalo_200c', 'zcold', 'cgf', 'mbasic', 'mbasic_200c', 'NFW_con', 'ssfr', 'weight_tot', 'env_512', 'env_1024', 'pop',\
                     'mAB_dA_total_cut_g_r', 'mAB_dA_total_cut_g_i', 'Tcons', 'fbar', 'BvT', 'rbulgevsrdisk', 'MzgasvsMzstar', 'zstar', 'zcold_zstar']
        
    property_dict       = get_property_dict()
    
    dt=[]
    for item in property_list:
        dt.extend([(item, property_dict[item]['dtype'])])
                            
    return dt    

def dt_MDPL_Galacticus_CMASS_v2():
    """get dt for Gal-dens file more properties version 2
    usual filename: Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3_more_props_sample_keys.txt    
    
    # haloid[idnumber](1) hostid[idnumber](2) env_512[0=vo,1=sh,2=fi,3=ko](3) orphan[0=cent,1=sat,2=orph](4) parentIndex[idnumber](5) mhalo[Msun](6) pop[flag](7) CMASS_sample_key[string_flag](8) nr_sample_keys[count](9) env_1024[0=vo,1=sh,2=fi,3=ko](10) x_pos[comvMpc](11) y_pos[comvMpc](12) z_pos[comvMpc](13) mbh[Msun](14) bheff[-](15) mstar[Msun](16) BvT[frac](17) SHMR[-](18) mcold[Msun](19) mhot[Msun](20) sfr[Msunyr-1](21) ssfr[yr-1](22) Mzgas[Msun](23) Mzstar[Msun](24) zcold[-](25) Tcons[Gyr](26) cgf[-](27) fbar[-](28) MAB_dA_total_g[-](29) MAB_dA_total_r[-](30) MAB_dA_total_i[-](31) mAB_dA_total_cut_r_i[-](32) mAB_dA_total_cut_g_r[-](33) mAB_dA_total_cut_g_i[-](34) Mzhot_halo[Msun](35) NFW_con[-](36) flag_CMASS_l[flag](37) flag_CMASS_h[flag](38) flag_CMASS_p[flag](39) flag_CMASS_a[flag](40) flag_CMASS_r[flag](41) flag_CMASS_b[flag](42) flag_CMASS_lz[flag](43) flag_CMASS_hz[flag](44) flag_CMASS_lzr[flag](45) flag_CMASS_hzr[flag](46) zstar[-](47) zcold_zstar[-](48) timeLastIsolated[Gyr](49) 
    
    """
    property_list = ['haloid', 'hostid', 'env_512', 'orphan', 'parentIndex', 'mhalo', 'pop', 'CMASS_sample_key', 'nr_sample_keys', 'env_1024',\
                     'x_pos', 'y_pos', 'z_pos', 'mbh', 'bheff', 'mstar', 'BvT', 'SHMR', 'mcold', 'mhot', 'sfr', 'ssfr', 'Mzgas', 'Mzstar',\
                     'zcold', 'Tcons', 'cgf', 'fbar',\
                     'MAB_dA_total_g', 'MAB_dA_total_r', 'MAB_dA_total_i', 'mAB_dA_total_cut_r_i', 'mAB_dA_total_cut_g_r', 'mAB_dA_total_cut_g_i',\
                     'Mzhot_halo', 'NFW_con',\
                     'flag_CMASS_l', 'flag_CMASS_h', 'flag_CMASS_p', 'flag_CMASS_a', 'flag_CMASS_r', 'flag_CMASS_b', 'flag_CMASS_lz',\
                     'flag_CMASS_hz', 'flag_CMASS_lzr', 'flag_CMASS_hzr',\
                     'zstar', 'zcold_zstar', 'timeLastIsolated', 'flag_CMASS_lsm', 'flag_CMASS_hsm', 'flag_CMASS_lhm', 'flag_CMASS_hhm']
        
    property_dict       = get_property_dict()
    
    dt=[]
    for item in property_list:
        dt.extend([(item, property_dict[item]['dtype'])])
                            
    return dt 

def dt_merger_trees_ASCII():

    dt =    [
            ('haloid1'          , np.int64),  
            ('haloid2'          , np.int64),                     
            ('descIndex1'       , np.int64), 
            ('descIndex2'       , np.int64),                  
            ('rootIndex'           , np.int64),
            ('snapid1'            , np.int32), 
            ('snapid2'            , np.int32),  
            ('z1'            , np.float32), 
            ('z2'            , np.float32),                   
            ('mhalo1'            , np.float32),
            ('mhalo2'            , np.float32),
            ('delta_mhalo'      , np.float32),
            ('delta_mhalo_perc'   , np.float32),                    
            ('rvir1'            , np.float32),
            ('rvir2'            , np.float32), 
            ('delta_rvir',        np.float32),
            ('delta_rvir_perc'   , np.float32),                     
            ('n_particles1'      , np.int64),
            ('n_particles2'      , np.int64),
            ('n_particles_shared', np.float32),
            ('n_particles_shared_perc1', np.float32),
            ('n_particles_shared_perc2', np.float32),                      
            ('x_pos1'             , np.float32),
            ('y_pos1'             , np.float32), 
            ('z_pos1'             , np.float32),
            ('x_pos2'             , np.float32),
            ('y_pos2'             , np.float32), 
            ('z_pos2'             , np.float32),                      
            ('delta_x_pos_perc'   , np.float32),                              
            ('delta_y_pos_perc'   , np.float32),                     
            ('delta_z_pos_perc'   , np.float32)  
            ]
    
    return dt


def fill_bucket_sfr2z(histo_bucket,
                      data,
                      data_offset,
                      redshift,
                      redshift_count,
                      column_count,
                      volume,
                      percent_low=31.7,
                      percent_high=68.3):


    header = '(1) z (2) sum(SFR [Msun yr-1]) (3) SFR/SFR(z=z_start) (4) 50th SFR [Msun yr-1] (5) '+perc_low+' SFR [Msun yr-1] (6) '+perc_high+' SFR [Msun yr-1] '\
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
             +'(70) 50th vdisp [kms-1] (71) '+perc_low+' vdisp [kms-1]  (72) '+perc_high+' vdisp [kms-1]'
    data_format='%0.4f\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.5f\t%0.5f\t%0.5f\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.5f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.5e\t%0.5e\t%0.5e\t%0.5e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.8e\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f\t%0.5f'



    #calculate properties        
    histo_bucket[redshift_count, data_offset*(column_count)]    = redshift
    histo_bucket[redshift_count, data_offset*(column_count) +1] = np.nansum(myfiltered_data['sfr'])  

    #normalise sfr to sfr(z=0)
    histo_bucket[redshift_count, data_offset*(column_count) +2] =  histo_bucket[redshift_count, data_offset*(column_count) +1] / histo_bucket[0, data_offset*(column_count) +1]
    histo_bucket[redshift_count, data_offset*(column_count) +3] =  np.nanmedian(myfiltered_data['sfr']) 
    histo_bucket[redshift_count, data_offset*(column_count) +4] =  np.nanpercentile(myfiltered_data['sfr'],percent_low)
    histo_bucket[redshift_count, data_offset*(column_count) +5] =  np.nanpercentile(myfiltered_data['sfr'],percent_high)
 
 
    histo_bucket[redshift_count, data_offset*(column_count) +6] =  np.nansum(myfiltered_data['ssfr'])
    histo_bucket[redshift_count, data_offset*(column_count) +7] =  histo_bucket[redshift_count, data_offset*(column_count) +6]/histo_bucket[0, data_offset*(column_count) +6]
    histo_bucket[redshift_count, data_offset*(column_count) +8] =  np.nanmedian(myfiltered_data['ssfr']) 
    histo_bucket[redshift_count, data_offset*(column_count) +9] =  np.nanpercentile(myfiltered_data['ssfr'],percent_low)
    histo_bucket[redshift_count, data_offset*(column_count) +10] =  np.nanpercentile(myfiltered_data['ssfr'],percent_high)        
    
    try:
        histo_bucket[redshift_count, data_offset*(column_count) +11]= np.nanmedian(myfiltered_data[mag_prefix+'g']-myfiltered_data[mag_prefix+'i'])
        histo_bucket[redshift_count, data_offset*(column_count) +12]= np.nanpercentile((myfiltered_data[mag_prefix+'g']-myfiltered_data[mag_prefix+'i']),percent_low)
        histo_bucket[redshift_count, data_offset*(column_count) +13]= np.nanpercentile((myfiltered_data[mag_prefix+'g']-myfiltered_data[mag_prefix+'i']),percent_high)
    except:
        print 'g-i color calculation faild!'

    try:
        histo_bucket[redshift_count, data_offset*(column_count) +29]= np.nanmedian(myfiltered_data[mag_prefix+'r']-myfiltered_data[mag_prefix+'i']) 
        histo_bucket[redshift_count, data_offset*(column_count) +30]= np.nanpercentile((myfiltered_data[mag_prefix+'r']-myfiltered_data[mag_prefix+'i']),percent_low)
        histo_bucket[redshift_count, data_offset*(column_count) +31]= np.nanpercentile((myfiltered_data[mag_prefix+'r']-myfiltered_data[mag_prefix+'i']),percent_high)
    except:
        print 'r-i color calculation faild!' 

       
    histo_bucket[redshift_count, data_offset*(column_count) +14]= np.nansum(myfiltered_data['mstar'])
    histo_bucket[redshift_count, data_offset*(column_count) +15]= histo_bucket[redshift_count, data_offset*(column_count) +14]/histo_bucket[0, data_offset*(column_count) +14]
    histo_bucket[redshift_count, data_offset*(column_count) +16]= np.nanmedian(myfiltered_data['mstar'])
    histo_bucket[redshift_count, data_offset*(column_count) +17]= np.nanpercentile(myfiltered_data['mstar'],percent_low)
    histo_bucket[redshift_count, data_offset*(column_count) +18]= np.nanpercentile(myfiltered_data['mstar'],percent_high)

    histo_bucket[redshift_count, data_offset*(column_count) +19]= np.nanmedian(myfiltered_data['rdisk'])
    histo_bucket[redshift_count, data_offset*(column_count) +20]= np.nanpercentile(myfiltered_data['rdisk'],percent_low)
    histo_bucket[redshift_count, data_offset*(column_count) +21]= np.nanmedian(myfiltered_data['rbulge']/myfiltered_data['rdisk'])
    histo_bucket[redshift_count, data_offset*(column_count) +22]= np.nanpercentile((myfiltered_data['rbulge']/myfiltered_data['rdisk']),percent_low)
    histo_bucket[redshift_count, data_offset*(column_count) +23]= np.nanpercentile((myfiltered_data['rbulge']/myfiltered_data['rdisk']),percent_high)        

    histo_bucket[redshift_count, data_offset*(column_count) +24] = cd.lookback_time(histo_bucket[redshift_count, data_offset*(column_count)], z0=histo_bucket[0, data_offset*(column_count)], **fidcosmo)/3.1536e+16 #seconds to Gyrs

    data_cents=myfiltered_data[np.where(myfiltered_data['orphan']==0)[0][:]]    
    data_sats=myfiltered_data[np.where(myfiltered_data['orphan']==1)[0][:]]
    data_os=myfiltered_data[np.where(myfiltered_data['orphan']==2)[0][:]]

    histo_bucket[redshift_count, data_offset*(column_count) +25] = format(data_cents.size/float(myfiltered_data.size), '0.3f')
    histo_bucket[redshift_count, data_offset*(column_count) +26] = format(data_sats.size/float(myfiltered_data.size), '0.3f')
    histo_bucket[redshift_count, data_offset*(column_count) +27] = format(data_os.size/float(myfiltered_data.size), '0.3f')

    histo_bucket[redshift_count, data_offset*(column_count) +28] = format(myfiltered_data.size/float(volume)/1e-4,'0.3f')

    try:        
        histo_bucket[redshift_count, data_offset*(column_count) +32]= np.nanmedian(myfiltered_data['mbh'])
        histo_bucket[redshift_count, data_offset*(column_count) +33]= np.nanpercentile(myfiltered_data['mbh'],percent_low)
        histo_bucket[redshift_count, data_offset*(column_count) +34]= np.nanpercentile(myfiltered_data['mbh'],percent_high)
    except:
        pass
    try:             
        histo_bucket[redshift_count, data_offset*(column_count) +35]= np.nanmedian(myfiltered_data['rhalf_mass'])
        histo_bucket[redshift_count, data_offset*(column_count) +36]= np.nanpercentile(myfiltered_data['rhalf_mass'],percent_low)
        histo_bucket[redshift_count, data_offset*(column_count) +37]= np.nanpercentile(myfiltered_data['rhalf_mass'],percent_high)          
    except:
        pass
    try: 
        histo_bucket[redshift_count, data_offset*(column_count) +38]= np.nansum(myfiltered_data['mcold'])
        histo_bucket[redshift_count, data_offset*(column_count) +39]= histo_bucket[redshift_count, data_offset*(column_count) +38]/histo_bucket[0, data_offset*(column_count) +38]
        histo_bucket[redshift_count, data_offset*(column_count) +40]= np.nanmedian(myfiltered_data['mcold'])
        histo_bucket[redshift_count, data_offset*(column_count) +41]= np.nanpercentile(myfiltered_data['mcold'],percent_low)
        histo_bucket[redshift_count, data_offset*(column_count) +42]= np.nanpercentile(myfiltered_data['mcold'],percent_high)        
    except:
        pass    
    try: 
        histo_bucket[redshift_count, data_offset*(column_count) +43]= np.nansum(myfiltered_data['Mzgas'])
        histo_bucket[redshift_count, data_offset*(column_count) +44]= histo_bucket[redshift_count, data_offset*(column_count) +43]/histo_bucket[0, data_offset*(column_count) +43]
        histo_bucket[redshift_count, data_offset*(column_count) +45]= np.nanmedian(myfiltered_data['Mzgas'])
        histo_bucket[redshift_count, data_offset*(column_count) +46]= np.nanpercentile(myfiltered_data['Mzgas'],percent_low)
        histo_bucket[redshift_count, data_offset*(column_count) +47]= np.nanpercentile(myfiltered_data['Mzgas'],percent_high)
    except:
        pass    
    try: 
        histo_bucket[redshift_count, data_offset*(column_count) +48]= np.nanmedian(myfiltered_data['zcold'])
        histo_bucket[redshift_count, data_offset*(column_count) +49]= np.nanpercentile(myfiltered_data['zcold'],percent_low)
        histo_bucket[redshift_count, data_offset*(column_count) +50]= np.nanpercentile(myfiltered_data['zcold'],percent_high)
    except:
        pass
    try: 
        histo_bucket[redshift_count, data_offset*(column_count) +52]= np.nansum(myfiltered_data['mhalo'])
        histo_bucket[redshift_count, data_offset*(column_count) +53]= histo_bucket[redshift_count, data_offset*(column_count) +52]/histo_bucket[0, data_offset*(column_count) +52]
        histo_bucket[redshift_count, data_offset*(column_count) +54]= np.nanmedian(myfiltered_data['mhalo'])
        histo_bucket[redshift_count, data_offset*(column_count) +55]= np.nanpercentile(myfiltered_data['mhalo'],percent_low)
        histo_bucket[redshift_count, data_offset*(column_count) +56]= np.nanpercentile(myfiltered_data['mhalo'],percent_high)
    except:
        pass
    try: 
        histo_bucket[redshift_count, data_offset*(column_count) +57]= np.nanmedian(myfiltered_data['mstar']/myfiltered_data['mhalo'])
        histo_bucket[redshift_count, data_offset*(column_count) +58]= np.nanpercentile(myfiltered_data['mstar']/myfiltered_data['mhalo'],percent_low)
        histo_bucket[redshift_count, data_offset*(column_count) +59]= np.nanpercentile(myfiltered_data['mstar']/myfiltered_data['mhalo'],percent_high)            
    except:
        pass

    try: 
        histo_bucket[redshift_count, data_offset*(column_count) +60]= np.nanmedian(myfiltered_data['mean_age_stars_disk'])
        histo_bucket[redshift_count, data_offset*(column_count) +61]= np.nanpercentile(myfiltered_data['mean_age_stars_disk'],percent_low)
        histo_bucket[redshift_count, data_offset*(column_count) +62]= np.nanpercentile(myfiltered_data['mean_age_stars_disk'],percent_high)            
    except:
        pass
    try: 
        histo_bucket[redshift_count, data_offset*(column_count) +63]= np.nanmedian(myfiltered_data['mean_age_stars_spheroid'])
        histo_bucket[redshift_count, data_offset*(column_count) +64]= np.nanpercentile(myfiltered_data['mean_age_stars_spheroid'],percent_low)
        histo_bucket[redshift_count, data_offset*(column_count) +65]= np.nanpercentile(myfiltered_data['mean_age_stars_spheroid'],percent_high)            
    except:
        pass
    try: 
        histo_bucket[redshift_count, data_offset*(column_count) +66]= np.nanmedian(myfiltered_data['vmax'])
        histo_bucket[redshift_count, data_offset*(column_count) +67]= np.nanpercentile(myfiltered_data['vmax'],percent_low)
        histo_bucket[redshift_count, data_offset*(column_count) +68]= np.nanpercentile(myfiltered_data['vmax'],percent_high)            
    except:
        pass
    try: 
        histo_bucket[redshift_count, data_offset*(column_count) +69]= np.nanmedian(myfiltered_data['vdisp'])
        histo_bucket[redshift_count, data_offset*(column_count) +70]= np.nanpercentile(myfiltered_data['vdisp'],percent_low)
        histo_bucket[redshift_count, data_offset*(column_count) +71]= np.nanpercentile(myfiltered_data['vdisp'],percent_high)            
    except:
        pass
    
    histo_bucket[redshift_count, data_offset*(column_count) +51]= histo_bucket[redshift_count, data_offset*(column_count) +1]/volume
    
    return histo_bucket




def get_property_dict():
    """method returns a dictionary of all formats, units, and dtyes of all available properties. The
    property names serve as keys.
    """
    
    h=''
    format_float='%0.6f'
    format_exp='%0.7e'
    format_int='%i'
    dtype_prop=np.float32
    dtype_exp=np.float64
    dtype_ID_large=np.int64
    dtype_ID_small=np.int32
    dtype_flag=np.int8
    avg_calc='50th'
    
    property_dict={
                #standard properties
                'mhalo':            {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mvir', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'mhalo_200c':       {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mhalo_200c', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'mhalo_sat':        {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mvir_sat', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'mbasic':           {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mbasic', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'mbasic_200c':      {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mbasic_200c', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},


                'mstar':            {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mstar', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'mbh':              {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mbh', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'mcold':            {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mcold', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'mhot':             {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mhot', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'Mzgas':            {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mzgas', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'Mzstar':           {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mzstar', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'Mzhot_halo':       {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mzhot_halo', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},

                'mstar_spheroid':   {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mstar_spheroid', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'mstar_disk':       {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mstar_disk', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'mcold_spheroid':   {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mcold_spheroid', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'mcold_disk':       {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mcold_disk', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'Mzgas_spheroid':   {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mzgas_spheroid', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'Mzgas_disk':       {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mzgas_disk', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'Mzstar_spheroid':  {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mzstar_spheroid', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'Mzstar_disk':      {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mzstar_disk', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'Mzhot':            {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mzhot', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},

                'NFW_con':  {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'NFW_con', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},               
                
                'sfr':              {'format': format_exp, 'unit': h+'Msunyr-1', 'dtype': dtype_prop, 'output_prop_as': 'SFR', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},              
                'sfr_spheroid':     {'format': format_exp, 'unit': h+'Msunyr-1', 'dtype': dtype_prop, 'output_prop_as': 'SFR_spheroid', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},              
                'sfr_disk':         {'format': format_exp, 'unit': h+'Msunyr-1', 'dtype': dtype_prop, 'output_prop_as': 'SFR_disk', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},              
               
                'rvir':         {'format': format_float, 'unit': h+'kpc', 'dtype': dtype_prop, 'output_prop_as': 'Rvir', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'rdisk':        {'format': format_float, 'unit': h+'kpc', 'dtype': dtype_prop, 'output_prop_as': 'Rdisk', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'rbulge':       {'format': format_float, 'unit': h+'kpc', 'dtype': dtype_prop, 'output_prop_as': 'Rbulge', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'rhalf_mass':   {'format': format_float, 'unit': h+'kpc', 'dtype': dtype_prop, 'output_prop_as': 'R1/2_halo', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},

                'vmax':         {'format': format_float, 'unit': h+'kms-1', 'dtype': dtype_prop, 'output_prop_as': 'Vmax', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},
                'vdisp':        {'format': format_float, 'unit': h+'kms-1', 'dtype': dtype_prop, 'output_prop_as': 'Vdisp', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},
                'vdisp_30kpc_T19':        {'format': format_float, 'unit': h+'kms-1', 'dtype': dtype_prop, 'output_prop_as': 'Vdisp', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},
                'Vrot2Vdisp_30kpc_T19':   {'format': format_float, 'unit': '', 'dtype': dtype_prop, 'output_prop_as': 'Vdisp', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},


                'vdisk':        {'format': format_float, 'unit': h+'kms-1', 'dtype': dtype_prop, 'output_prop_as': 'Vdisk', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'vbulge':       {'format': format_float, 'unit': h+'kms-1', 'dtype': dtype_prop, 'output_prop_as': 'Vbulge', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},

                'x_vel':       {'format': format_float, 'unit': h+'kms-1', 'dtype': dtype_prop, 'output_prop_as': 'Vx', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},
                'y_vel':       {'format': format_float, 'unit': h+'kms-1', 'dtype': dtype_prop, 'output_prop_as': 'Vy', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},
                'z_vel':       {'format': format_float, 'unit': h+'kms-1', 'dtype': dtype_prop, 'output_prop_as': 'Vz', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},

                      
                'T_U':          {'format': format_float, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': '', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},

                #calculated properties
                'mbar':         {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Mbar', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'sum_mstar':    {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'sum(Mstar)', 'avg_calc': '', 'axis_style': 'log', 'plot_error_bars': True},
                'sum_mhalo':    {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'sum(Mvir)', 'avg_calc': '', 'axis_style': 'log', 'plot_error_bars': True},
                'sum_mcold':    {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'sum(Mcold)', 'avg_calc': '', 'axis_style': 'log', 'plot_error_bars': True},
                'sum_Mzgas':    {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'sum(Mzgas)', 'avg_calc': '', 'axis_style': 'log', 'plot_error_bars': True},

                'cSFRD':    {'format': format_exp, 'unit': h+'Msunyr-1Mpc-3', 'dtype': dtype_prop, 'output_prop_as': 'cSFRD', 'avg_calc': '', 'axis_style': 'log', 'plot_error_bars': False},
                
                'ssfr':         {'format': format_exp, 'unit': 'yr-1', 'dtype': dtype_prop, 'output_prop_as': 'sSFR', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'sum_sfr':      {'format': format_exp, 'unit': 'Msunyr-1', 'dtype': dtype_prop, 'output_prop_as': 'sum(SFR)', 'avg_calc': '', 'axis_style': 'log', 'plot_error_bars': False},
                'Tcons':        {'format': format_float, 'unit': 'Gyr', 'dtype': dtype_prop, 'output_prop_as': 'Mcold/SFR', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'mean_age_stars_spheroid': {'format': format_float, 'unit': 'Gyr', 'dtype': dtype_prop, 'output_prop_as': 'mean_age_stars_spheroid', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},
                'mean_age_stars_disk':     {'format': format_float, 'unit': 'Gyr', 'dtype': dtype_prop, 'output_prop_as': 'mean_age_stars_disk', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},                

                'spinParameter':{'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'spinParameter', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': False},


                'jbar':        {'format': format_exp, 'unit': 'Msunkpckms-1', 'dtype': dtype_prop, 'output_prop_as': 'jbar', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'jbulge':      {'format': format_exp, 'unit': 'Msunkpckms-1', 'dtype': dtype_prop, 'output_prop_as': 'jbulge', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'jdisk':       {'format': format_exp, 'unit': 'Msunkpckms-1', 'dtype': dtype_prop, 'output_prop_as': 'jdisk', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},                 
                
                'BvT':         {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'Mbulge/Mstar', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},
                'bheff':       {'format': format_exp, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'Mbh/Mvir', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'SHMR':        {'format': format_exp, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'Mstar/Mhalo', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},                
                'cgf':         {'format': format_exp, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'Mcold/Mstar', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'fbar':        {'format': format_exp, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'Mcold/Mbar', 'avg_calc': avg_calc, 'axis_style': 'log', 'plot_error_bars': True},
                'zcold':       {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'Zcold', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},   

                'g-i':            {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'g-i', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},
                'r-i':            {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'r-i', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},
                'zstar':       {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'zstar', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},   
                'zcold_zstar': {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'zcold-zstar', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},   

                'rbulgevsrdisk':  {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'Rbulge/Rdisk', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},
                'MzgasvsMzstar':  {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'Mzgas/Mzstar', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},

                #growth/variation depending properties
                'delta_mhalo':      {'format': format_exp, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': '', 'avg_calc': avg_calc},
                'delta_mhalo_perc': {'format': format_float, 'unit': '%', 'dtype': dtype_prop, 'output_prop_as': '', 'avg_calc': avg_calc},
                'delta_rvir':       {'format': format_float, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': '', 'avg_calc': avg_calc},
                'delta_rvir_perc':  {'format': format_float, 'unit': '%', 'dtype': dtype_prop, 'output_prop_as': '', 'avg_calc': avg_calc},                 

                'mhalo2mhalo_zstart':  {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'Mvir/Mvir_zstart', 'avg_calc': ''},                 
                'mstar2mstar_zstart':  {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'Mstar/star_zstart', 'avg_calc': ''},
                'mcold2mcold_zstart':  {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'Mcold/Mcold_zstart', 'avg_calc': ''}, 
                'Mzgas2Mzgas_zstart':  {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'Mvir/Mvir_zstart', 'avg_calc': ''}, 
                'sfr2sfr_zstart':  {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'sfr/sfr_zstart', 'avg_calc': ''}, 
                'ssfr2ssfr_zstart':  {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'ssfr/ssfr_zstart', 'avg_calc': ''}, 
                
                #temporal axis, positions & demographics
                'count':        {'format': format_int, 'unit': '-', 'dtype': dtype_ID_small, 'output_prop_as': '', 'avg_calc': ''},
                'n_count':      {'format': format_int, 'unit': '-', 'dtype': dtype_ID_small, 'output_prop_as': '', 'avg_calc': ''},
                'n_cents':      {'format': format_int, 'unit': 'count', 'dtype': dtype_ID_small, 'output_prop_as': 'N(centrals)', 'avg_calc': ''},                
                'n_sats':       {'format': format_int, 'unit': 'count', 'dtype': dtype_ID_small, 'output_prop_as': 'N(sats)', 'avg_calc': ''},               
                'n_orphans':    {'format': format_int, 'unit': 'count', 'dtype': dtype_ID_small, 'output_prop_as': 'N(orphans)', 'avg_calc': ''},
                'density':      {'format': format_exp, 'unit': h+'x10-4Mpc-3', 'dtype': dtype_prop, 'output_prop_as': 'n', 'avg_calc': ''},
                 
                'snapid':       {'format': format_int, 'unit': '-', 'dtype': dtype_ID_small, 'output_prop_as': '', 'avg_calc': ''},
                'z':            {'format': '%0.4f', 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'z', 'avg_calc': ''},
                'a':            {'format': '%0.4f', 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': '', 'avg_calc': ''},
                'a_desc':       {'format': '%0.4f', 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': '', 'avg_calc': ''},
                'a_lastMM':     {'format': '%0.4f', 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': '', 'avg_calc': ''},
                't':            {'format': '%0.4f', 'unit': 'Gyr', 'dtype': dtype_prop, 'output_prop_as': 'lookbacktime(z=z_start)', 'avg_calc': ''},
                
                'x_pos':        {'format': format_float, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'X', 'avg_calc': ''},
                'y_pos':        {'format': format_float, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Y', 'avg_calc': ''},
                'z_pos':        {'format': format_float, 'unit': h+'Msun', 'dtype': dtype_prop, 'output_prop_as': 'Z', 'avg_calc': ''},

                'satelliteMergeTime':  {'format': format_float, 'unit': 'Gyr', 'dtype': dtype_prop, 'output_prop_as': 'satelliteMergeTime', 'avg_calc': ''},
                'timeLastIsolated':  {'format': format_float, 'unit': 'Gyr', 'dtype': dtype_prop, 'output_prop_as': 'timeLastIsolated', 'avg_calc': ''},
                
                #IDs                     
                'haloid':       {'format': format_int, 'unit': 'ID', 'dtype': dtype_ID_large, 'output_prop_as': '', 'avg_calc': ''},
                'hostid':       {'format': format_int, 'unit': 'ID', 'dtype': dtype_ID_large, 'output_prop_as': '', 'avg_calc': ''},
                'haloid_CT':         {'format': format_int, 'unit': 'ID', 'dtype': dtype_ID_large, 'output_prop_as': '', 'avg_calc': ''},
                'orignal_haloid_RS': {'format': format_int, 'unit': 'ID', 'dtype': dtype_ID_large, 'output_prop_as': '', 'avg_calc': ''},                     
                'descID':       {'format': format_int, 'unit': 'ID', 'dtype': dtype_ID_large, 'output_prop_as': '', 'avg_calc': ''},
                'DFirstID':     {'format': format_int, 'unit': 'ID', 'dtype': dtype_ID_large, 'output_prop_as': '', 'avg_calc': ''}, 
                'LastDFirstID': {'format': format_int, 'unit': 'ID', 'dtype': dtype_ID_large, 'output_prop_as': '', 'avg_calc': ''}, 
                'LastMLDFirstID':{'format': format_int, 'unit': 'ID', 'dtype': dtype_ID_large, 'output_prop_as': '', 'avg_calc': ''},                    
                'rootIndex':    {'format': format_int, 'unit': 'ID', 'dtype': dtype_ID_large, 'output_prop_as': '', 'avg_calc': ''},
                'subTreeID':    {'format': format_int, 'unit': 'ID', 'dtype': dtype_ID_large, 'output_prop_as': '', 'avg_calc': ''},
                'nodeIndex':    {'format': format_int, 'unit': 'ID', 'dtype': dtype_ID_large, 'output_prop_as': 'nodeIndex', 'avg_calc': ''},
                'parentIndex':    {'format': format_int, 'unit': 'ID', 'dtype': dtype_ID_large, 'output_prop_as': 'parentIndex', 'avg_calc': ''},
                'satelliteNodeIndex':    {'format': format_int, 'unit': 'ID', 'dtype': dtype_ID_large, 'output_prop_as': 'satelliteNodeIndex', 'avg_calc': ''},
                'satelliteIndex':    {'format': format_int, 'unit': 'ID', 'dtype': dtype_ID_large, 'output_prop_as': 'satelliteIndex', 'avg_calc': ''},
                'siblingIndex':    {'format': format_int, 'unit': 'ID', 'dtype': dtype_ID_large, 'output_prop_as': 'siblingIndex', 'avg_calc': ''},

                #luminosities and magnitudes
                'L_SDSS_dA_total_u':    {'format': format_exp, 'unit': '4.4659e13WHz-1', 'dtype': dtype_prop, 'output_prop_as': 'L_SDSS_dA_total_u', 'avg_calc': ''},
                'L_SDSS_dA_total_g':    {'format': format_exp, 'unit': '4.4659e13WHz-1', 'dtype': dtype_prop, 'output_prop_as': 'L_SDSS_dA_total_g', 'avg_calc': ''},
                'L_SDSS_dA_total_r':    {'format': format_exp, 'unit': '4.4659e13WHz-1', 'dtype': dtype_prop, 'output_prop_as': 'L_SDSS_dA_total_r', 'avg_calc': ''},
                'L_SDSS_dA_total_i':    {'format': format_exp, 'unit': '4.4659e13WHz-1', 'dtype': dtype_prop, 'output_prop_as': 'L_SDSS_dA_total_i', 'avg_calc': ''},
                'L_SDSS_dA_total_z':    {'format': format_exp, 'unit': '4.4659e13WHz-1', 'dtype': dtype_prop, 'output_prop_as': 'L_SDSS_dA_total_z', 'avg_calc': ''},

                'MAB_dA_total_u':    {'format': format_float, 'unit': 'mag', 'dtype': dtype_prop, 'output_prop_as': 'MAB_dA_total_u', 'avg_calc': ''},
                'MAB_dA_total_g':    {'format': format_float, 'unit': 'mag', 'dtype': dtype_prop, 'output_prop_as': 'MAB_dA_total_g', 'avg_calc': ''},
                'MAB_dA_total_r':    {'format': format_float, 'unit': 'mag', 'dtype': dtype_prop, 'output_prop_as': 'MAB_dA_total_r', 'avg_calc': ''},
                'MAB_dA_total_i':    {'format': format_float, 'unit': 'mag', 'dtype': dtype_prop, 'output_prop_as': 'MAB_dA_total_i', 'avg_calc': ''},
                'MAB_dA_total_z':    {'format': format_float, 'unit': 'mag', 'dtype': dtype_prop, 'output_prop_as': 'MAB_dA_total_z', 'avg_calc': ''},

                'mAB_dA_total_u':    {'format': format_float, 'unit': 'mag', 'dtype': dtype_prop, 'output_prop_as': 'mAB_dA_total_u', 'avg_calc': ''},
                'mAB_dA_total_g':    {'format': format_float, 'unit': 'mag', 'dtype': dtype_prop, 'output_prop_as': 'mAB_dA_total_g', 'avg_calc': ''},
                'mAB_dA_total_r':    {'format': format_float, 'unit': 'mag', 'dtype': dtype_prop, 'output_prop_as': 'mAB_dA_total_r', 'avg_calc': ''},
                'mAB_dA_total_i':    {'format': format_float, 'unit': 'mag', 'dtype': dtype_prop, 'output_prop_as': 'mAB_dA_total_i', 'avg_calc': ''},
                'mAB_dA_total_z':    {'format': format_float, 'unit': 'mag', 'dtype': dtype_prop, 'output_prop_as': 'mAB_dA_total_z', 'avg_calc': ''},

                'mAB_dA_total_cut_g_i':         {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'g-i', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},
                'mAB_dA_total_cut_g_r':         {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'g-r', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},
                'mAB_dA_total_cut_r_i':         {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'r-i', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},
 
                #CMASS colour-cuts
                'mAB_dA_total_cut_dmesa':       {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'dmesa', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},
                'mAB_dA_total_cut_i_lt_dmesa':  {'format': format_float, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'i_lt_dmesa', 'avg_calc': avg_calc, 'axis_style': 'lin', 'plot_error_bars': True},

                #flags
                'orphan':        {'format': format_int, 'unit': '0=cent,1=sat,2=orph', 'dtype': dtype_flag, 'output_prop_as': 'orphan', 'avg_calc': ''},
                'env_512':        {'format': format_int, 'unit': '0=vo,1=sh,2=fi,3=no', 'dtype': dtype_flag, 'output_prop_as': 'env_512', 'avg_calc': ''},                
                'env_1024':        {'format': format_int, 'unit': '0=vo,1=sh,2=fi,3=no', 'dtype': dtype_flag, 'output_prop_as': 'env_1024', 'avg_calc': ''},
                'pop':        {'format': format_int, 'unit': 'flag', 'dtype': dtype_flag, 'output_prop_as': 'pop', 'avg_calc': ''},
                'isolated':        {'format': format_int, 'unit': 'isolated[0=non-iso,1=iso]', 'dtype': dtype_flag, 'output_prop_as': 'pop', 'avg_calc': ''},                 
                'CMASS_sample_key':  {'format': format_int, 'unit': 'flag', 'dtype': dtype_flag, 'output_prop_as': 'CMASS_sample_key', 'avg_calc': ''},
                'flag_CMASS_l':  {'format': format_int, 'unit': 'flag', 'dtype': dtype_flag, 'output_prop_as': 'flag_CMASS_l', 'avg_calc': ''},
                'flag_CMASS_h':  {'format': format_int, 'unit': 'flag', 'dtype': dtype_flag, 'output_prop_as': 'flag_CMASS_h', 'avg_calc': ''},
                'flag_CMASS_p':  {'format': format_int, 'unit': 'flag', 'dtype': dtype_flag, 'output_prop_as': 'flag_CMASS_p', 'avg_calc': ''},
                'flag_CMASS_a':  {'format': format_int, 'unit': 'flag', 'dtype': dtype_flag, 'output_prop_as': 'flag_CMASS_a', 'avg_calc': ''},
                'flag_CMASS_r':  {'format': format_int, 'unit': 'flag', 'dtype': dtype_flag, 'output_prop_as': 'flag_CMASS_r', 'avg_calc': ''},
                'flag_CMASS_b':  {'format': format_int, 'unit': 'flag', 'dtype': dtype_flag, 'output_prop_as': 'flag_CMASS_b', 'avg_calc': ''},
                'flag_CMASS_lz':  {'format': format_int, 'unit': 'flag', 'dtype': dtype_flag, 'output_prop_as': 'flag_CMASS_lz', 'avg_calc': ''},
                'flag_CMASS_hz':  {'format': format_int, 'unit': 'flag', 'dtype': dtype_flag, 'output_prop_as': 'flag_CMASS_hz', 'avg_calc': ''},
                'flag_CMASS_lzr':  {'format': format_int, 'unit': 'flag', 'dtype': dtype_flag, 'output_prop_as': 'flag_CMASS_lzr', 'avg_calc': ''},
                'flag_CMASS_hzr':  {'format': format_int, 'unit': 'flag', 'dtype': dtype_flag, 'output_prop_as': 'flag_CMASS_hzr', 'avg_calc': ''},
                'flag_CMASS_lsm':  {'format': format_int, 'unit': 'flag', 'dtype': dtype_flag, 'output_prop_as': 'flag_CMASS_lsm', 'avg_calc': ''},
                'flag_CMASS_hsm':  {'format': format_int, 'unit': 'flag', 'dtype': dtype_flag, 'output_prop_as': 'flag_CMASS_hsm', 'avg_calc': ''},
                'flag_CMASS_lhm':  {'format': format_int, 'unit': 'flag', 'dtype': dtype_flag, 'output_prop_as': 'flag_CMASS_lhm', 'avg_calc': ''},
                'flag_CMASS_hhm':  {'format': format_int, 'unit': 'flag', 'dtype': dtype_flag, 'output_prop_as': 'flag_CMASS_hhm', 'avg_calc': ''},
                #other
                'weight_tot':        {'format': format_int, 'unit': '-', 'dtype': dtype_prop, 'output_prop_as': 'weight_tot', 'avg_calc': ''},\
                'nr_sample_keys':    {'format': format_int, 'unit': 'count', 'dtype': dtype_prop, 'output_prop_as': 'nr_sample_keys', 'avg_calc': ''}
                }

    return property_dict

def crossmatch_catalogs(data,
                        prop_to_cross=False,
                        key=False,
                        redshift=None):

    import pandas as pd

    if key=='vweb':    
        print '\n+++++++++++++++++++\nset environment',
        #data2cross = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/data/Galacticus_1Gpc/Galac_with_Environments_new.txt', data_format='ASCII', data_shape='shaped', delim=' ', mydtype=np.float32, skiprow=2) 
                        
        try:
            data['env_512']=data2cross[:,3]
            print 'env_512 set!'
        except:
            print 'env_512 not excisting ...'                        
        try:
            data['env_1024']=data2cross[:,4]
            print 'env_1024 set!'
        except:
            print 'env_1024 not excisting ...\nalternative properties'
            data2cross = myData.readAnyFormat(config=False, mypath=mycomp+'anaconda/pro/data/Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3_haloid.txt', data_format='ASCII', data_shape='shaped', delim='  ', mydtype=np.uint64, skiprow=2) 
            print data2cross[0:10,0]

            test = np.in1d(data2cross[:,1],data['hostid'])
            data2cross=data2cross[test]
            print data2cross[:,0].size
            for i, hostid in enumerate(data2cross[:,1]):
                #print 'hostid', hostid, 'i:', i, 'check hostid:', data2cross[i,1], 'env_512:',  data2cross[i,3], 'env_1024:', data2cross[i,4]
                data['env_512'][np.where(data['hostid']==hostid)[:][0]]=data2cross[i,3]
            data['env_1024'][np.where(data['hostid']==hostid)[:][0]]=data2cross[i,4]
            print '--> DONE!\n'
 
    elif key=='Yetli':
        #Crossmatch catalog from Patricia and Yetli for environment affiliation in extended galaxies
        
        names=['haloid', 'hostid', 'envr', 'lastMajorM', 'lastMinorM', 'mstar_30kpc',\
               'Mgas_HI_60kpc_GK11', 'Mgas_H2_60kpc_GK11', 'Mgas_HI_60kpc_K13', 'Mgas_H2_60kpc_K13',\
               'Mgas_HI_30kpc_GK11', 'Mgas_H2_30kpc_GK11', 'Mgas_HI_30kpc_K13', 'Mgas_H2_30kpc_K13',\
               'MAB_SDSS_u', 'MAB_SDSS_g', 'MAB_SDSS_r', 'MAB_SDSS_i', 'MAB_SDSS_z', 'mAB_total_cut_r_i', 'mAB_total_cut_g_r', 'mAB_total_cut_g_i']
            
        myfilename2cross=mycomp+'anaconda/pro/data/EAGLE_100Mpc/EAGLE_100Mpc_z_0.0_tarsel_Mgas_Yetli_January2023.txt'
               
        data2cross = pd.read_csv(myfilename2cross,\
                                 names=names, skiprows=2, sep='  ')
        data2cross = df_to_sarray(data2cross)
                
        print 'CROSSMATCH key: "', key, '"\n-----------------\nparent catalog:', data.size, 'data2cross:', data2cross.size,

        test = np.in1d(data['haloid'],data2cross['haloid'])

        print '-->', test[np.where(test==True)].size, 'found to crossmatch!'
        
        #set all columns of the parent data to -99, so only the properties of the galaxies which are exsist in the sample
        #which should be crossmatched will be filled! But do not set the 'fofID' to -99 or the properties which are in both catalogs!!!
        
        new_props=['envr']#, 'Mgas_HI_30kpc_GK11', 'Mgas_H2_30kpc_GK11']#,\
                    #'Mgas_HI_60kpc_GK11', 'Mgas_H2_60kpc_GK11',\
                    #'Mgas_HI_30kpc_K13', 'Mgas_H2_30kpc_K13',  'Mgas_HI_60kpc_K13', 'Mgas_H2_60kpc_K13',\
                    #'MAB_SDSS_u', 'MAB_SDSS_g', 'MAB_SDSS_r', 'MAB_SDSS_i', 'MAB_SDSS_z']#,\
                    # #'mAB_total_cut_r_i', 'mAB_total_cut_g_r']
     
        data[new_props]=-99
               
        parentIndices, index1, index2 = np.intersect1d(data['haloid'], data2cross['haloid'], return_indices=True) 
        
        for prop in new_props:
            #print 'property:', prop
            data[prop][index1]=data2cross[prop][index2]
            
        # for env in [-99, 0, 1, 2, 3]:            
        #     print 'envr:', env, data[np.where(data['envr']==env)[:][0]].size
               
        
        #print data[['haloid','fofID', 'envr', 'mstar_30kpc', 'Mgas_HI_30kpc_GK11', 'MAB_total_i']][0:20]
        print '--> DONE!\n-----------------\n\n'

    elif key=='Yetli_ev':
        #Crossmatch catalog from Patricia and Yetli for environment affiliation in extended galaxies
        
        names=['haloid', 'hostid', 'envr', 'mhalo_200c', 'mstar_30kpc', 'sfr_30kpc', 'mgas_30kpc', 'DvT_stars',\
               'mbh', 'bh_acc_rate', 't50_stars', 't70_stars', 'rVoid', 'DtoVoidCent', 'OH_gas_30kpc', 
               'z_grad', 'L_bolo', 'lastMinorM', 'lastMajorM']
            
        myfilename2cross=mycomp+'anaconda/pro/data/EAGLE_ev_100Mpc/EAGLE_ev_100Mpc_z_0.0_tarsel_7Nov2023.txt'
               
        data2cross = pd.read_csv(myfilename2cross,\
                                 names=names, skiprows=2, sep='  ')
        data2cross = df_to_sarray(data2cross)
        
        
        print 'CROSSMATCH key: "', key, '"\n-----------------\nparent catalog:', data.size, 'data2cross:', data2cross.size,
        
        test = np.in1d(data['haloid'],data2cross['haloid'])

        print '-->', test[np.where(test==True)].size, 'found to crossmatch!'
       
        #set all columns of the parent data to -99, so only the properties of the galaxies which are exsist in the sample
        #which should be crossmatched will be filled! But do not set the 'fofID' to -99 or the properties which are in both catalogs!!!
        
        new_props=['DtoVoidCent', 'rVoid']
     
        data[new_props]=-99   
        
        parentIndices, index1, index2 = np.intersect1d(data['haloid'], data2cross['haloid'], return_indices=True) 
        
        for prop in new_props:
            #print 'property:', prop
            data[prop][index1]=data2cross[prop][index2]            
            
        # for env in [-99, 0, 1, 2, 3]:            
        #     print 'env:', env, data[np.where(data['envr']==env)[:][0]].size
               
        
        #print data[['haloid','fofID', 'envr', 'mstar_30kpc', 'Mgas_HI_30kpc_GK11', 'MAB_total_i']][0:20]
        print '--> DONE!\n-----------------\n\n'

    elif key=='EAGLE_MKC':
   
        names=['fofID', 'haloid', 'hostid', 'galaxyID', 'nsubhalos', 'mPart_stars', 'mstar_30kpc', 'sfr_30kpc', 'vdisp_30kpc',\
               'mhalo_200c', 'r200c', 'rhalf_stars', 'rhalf_stars_2D', 'rhalf_gas', 'rhalf_gas_2D', 'vmax', 'rVmax',\
               'xpos', 'ypos', 'zpos', 'DvT_30kpc_T19', 'vdisp_30kpc_T19', 'ell_30kpc_T19', 'Vrot2Vdisp_30kpc_T19', 'ellDM_30kpc_T19',\
               'rhalf_stars_30kpc', 'rhalf_stars_30kpc_2D',\
               'F_dA_Johnson_U', 'F_dA_Johnson_B', 'F_dA_Johnson_V', 'F_dA_Johnson_R',\
               'F_dA_SDSS_g', 'F_dA_SDSS_r', 'F_dA_SDSS_i',\
               'MAB_dA_Johnson_U', 'MAB_dA_Johnson_B', 'MAB_dA_Johnson_V', 'MAB_dA_Johnson_R',\
               'MAB_dA_SDSS_g', 'MAB_dA_SDSS_r', 'MAB_dA_SDSS_i',\
               'MAB_SDSS_g', 'MAB_SDSS_r', 'MAB_SDSS_i',\
               'z_mean_birth_stars', 'age_mean_stars']

            
        myfilename2cross=mycomp+'anaconda/pro/data/EAGLE_100Mpc/original_catalogs/EAGLE_100Mpc_full_catalog_MKC.dat'
        
        data2cross = pd.read_csv(myfilename2cross,\
                                 names=names, skiprows=2, sep=' ')
        data2cross = df_to_sarray(data2cross)
        #print np.info(data2cross)
        
        #data2cross['mAB_dA_total_U'] =conv_fluxJy_to_ABmag(data2cross['mAB_dA_total_U'])
        # for filter_band in ['U', 'B', 'V', 'g', 'r', 'i']:
        #     data2cross['mAB_dA_total_'+filter_band] = -2.5*np.log10(data2cross['mAB_dA_total_'+filter_band])  + 8.926
        
        #print data2cross[['fofID', 'haloid', 'mstar_30kpc']]
        print 'CROSSMATCH key: "', key, '"\n-----------------\nparent catalog:', data.size, 'data2cross:', data2cross.size,

        test = np.in1d(data['fofID'],data2cross['fofID'])

        #data2cross=data2cross[test]
        #data2cross[props[2::]]
        print '-->', test[np.where(test==True)].size, 'found to crossmatch!' 
        #exit()
        #set all columns of the parent data to -99, so only the properties of the galaxies which are exsist in the sample
        #which should be crossmatched will be filled! But do not set the 'fofID' to -99 or the properties which are in both catalogs!!!
        

        new_props=['rhalf_stars_30kpc',\
                   'MAB_dA_Johnson_B', 'MAB_dA_Johnson_V', 'MAB_dA_SDSS_g', 'MAB_dA_SDSS_r', 'MAB_dA_SDSS_i',\
                   'vmax', 'rVmax',\
                   'vdisp_30kpc_T19', 'Vrot2Vdisp_30kpc_T19', 'ellDM_30kpc_T19']
            
        data[new_props]=-99       
        
        parentIndices, index1, index2 = np.intersect1d(data['fofID'], data2cross['fofID'], return_indices=True) 
        
        
        for prop in new_props:
            #print 'property:', prop
            data[prop][index1]=data2cross[prop][index2]
                       
        print '--> DONE!\n-----------------\n\n'

    elif key=='EAGLE_ev_full':
   
        myfilename2cross=mycomp+'anaconda/pro/data/EAGLE_ev_100Mpc/EAGLE_ev_100Mpc_full_merger_trees.txt'
        
        data2cross= pd.read_csv(myfilename2cross, skiprows=1,\
                          names=['fofID', 'descfofID', 'progRedshift', 'progHaloid', 'descHaloid', 'progHostid', 'GalaxyID', 'lastProgID', 'topLeafID',\
                                 'mhalo_200c', 'r200c', 'nsats', 'xpos_subhalo', 'ypos_subhalo', 'zpos_subhalo',\
                                 'mPart_DM', 'mPart_BH', 'mPart_gas', 'mPart_stars', 'sfr', 'vmax', 'rVmax',\
                                 'mstar_30kpc', 'sfr_30kpc', 'vdisp_30kpc', 'mgas_30kpc', 'mhalo_30kpc',
                                 'mgas_SF', 'mgas_NSF', 'spinNSFGas_x', 'spinNSFGas_y', 'spinNSFGas_z',\
                                 'spinSFGas_x', 'spinSFGas_y', 'spinSFGas_z', 'spinStars_x', 'spinStars_y', 'spinStars_z'],\
                                 sep=',') 
        data2cross = df_to_sarray(data2cross)
                
        #print data2cross[['fofID', 'haloid', 'mstar_30kpc']]
        print 'CROSSMATCH key: "', key, '"\n-----------------\nparent catalog:', data.size, 'data2cross:', data2cross.size,

        from cosmolopy import cparam, cd, cc

        fidcosmo = cparam.PlanckEAGLE(flat=True, extras=False)
        #print fidcosmo
        

        from scipy import interpolate
        
        count_haloid0_missing=0
        array_lost_fofIdz0=[]
        
        for fofID in data['fofID']:
            print 'fofID:', fofID,
            #print  data2cross[['redshift', 'fofID', 'mstar_30kpc', 'vmax', 'rVmax', 'rhalf_stars_30kpc', 'mbh']][np.where(data2cross['fofID']==fofID)[:][0]]
            #print  data[['fofID', 'mstar_30kpc', 'vmax', 'rVmax', 'rhalf_stars_30kpc', 'mbh']][np.where(data['fofID']==fofID)[:][0]]
            topLeafIDz0=data2cross['topLeafID'][(np.where((data2cross['fofID']==fofID) & (data2cross['redshift']<=5e-16) & (data2cross['progHostid']==0)))[:][0]]                       
            
            print 'topLeafIDz0:', topLeafIDz0
            
            if len(topLeafIDz0)==0:
                print 'central not found, moving to lowest progenitor hostid available! Choosen topLeafID -->',
                try:
                    checkDataTopLeafID=data2cross[['redshift','topLeafID', 'progHostid']][(np.where((data2cross['fofID']==fofID) & (data2cross['redshift']<=5e-16)))[:][0]]
                  
                    checkDataTopLeafID.sort(order=['redshift','progHostid'], axis=0)        
                    #print checkDataTopLeafID,
                    topLeafIDz0=checkDataTopLeafID['topLeafID'][0]
                    #print ' new topLeafID:', topLeafIDz0                   
    
                    count_haloid0_missing+=1
                
                except:
                    print 'fofID not present at z=0!\n'
                    array_lost_fofIdz0.append(fofID)
            try:              
                for prop, key in zip(['mhalo_200c', 'mstar_30kpc', 'mbh'], ['halo', 'stars', 'bh']):
                    
                    #print '\tcalculate formation times --> fofid:', fofID, 'prop:', prop, 'key:', key
                    
                    data_ev_fofID=data2cross[np.in1d(data2cross['topLeafID'], topLeafIDz0)]
                                        
                    data_ev_fofID.sort(order=['redshift'], axis=0)
    
                    #print data_ev_fofID[['fofID', 'redshift', 'topLeafID', prop]]
                    
                    data_ev_fofID[prop] = data_ev_fofID[prop]/data_ev_fofID[prop][0]
                    #print data_ev_fofID[['fofID', 'redshift', prop]]
                    f = interpolate.interp1d(data_ev_fofID[prop], data_ev_fofID['redshift'])
                    #print f(0.5)                                 
                    
                    data['t50_'+key][np.where(data['fofID']==fofID)[:][0]] = cd.lookback_time(f(0.5), z0=0.0, **fidcosmo)/3.1536e+16 #seconds to Gyrs                           
                    data['t70_'+key][np.where(data['fofID']==fofID)[:][0]] = cd.lookback_time(f(0.7), z0=0.0, **fidcosmo)/3.1536e+16 #seconds to Gyrs
                    
                    
                    data['delta_tform_'+key][np.where(data['fofID']==fofID)[:][0]] = data['t50_'+key][np.where(data['fofID']==fofID)[:][0]] - data['t70_'+key][np.where(data['fofID']==fofID)[:][0]]
                    
                    #print  data_ev_fofID[['fofID', prop]][np.where(data_ev_fofID['fofID']==fofID)[:][0]]
                    #print  data[['fofID', 'haloid', 'hostid', 'galaxyID', prop, 't50_'+key, 't70_'+key, 'delta_tform_'+key]][np.where(data['fofID']==fofID)[:][0]]
    

                for prop, key in zip(['bh_acc_rate','rVmax', 'vmax', 'vdisp', 'sfr', 'ssfr', 'angM_stars', 'angM_SFgas', 'angM_NSFgas', 'mgas_SF', 'mgas_NSF'], ['bh_acc_rate','rVmax', 'vmax', 'vdisp', 'sfr', 'ssfr', 'angM_stars', 'angM_SFgas', 'angM_NSFgas', 'mgas_SF', 'mgas_NSF']):
                    
                    #when did the properties hold their maximal values:
                    #vmax, angular momenta, vdisp, sfr, ssfr, gas masses, velocities
                    
                    data_ev_fofID=data2cross[np.in1d(data2cross['topLeafID'], topLeafIDz0)]
    
    
                    data[prop][np.where(data['fofID']==fofID)[:][0]] = data_ev_fofID[prop][0]                   
    
                    #print data_ev_fofID[['fofID', 'redshift', 'topLeafID', prop]]
                                       
                    data_ev_fofID.sort(order=[prop], axis=0)
    
                    #print data_ev_fofID[['fofID', 'redshift', 'topLeafID', prop]]                             
    
    
                    data['tmin_'+key][np.where(data['fofID']==fofID)[:][0]] = cd.lookback_time(data_ev_fofID['redshift'][0], z0=0.0, **fidcosmo)/3.1536e+16 #seconds to Gyrs                           
                    data['tmax_'+key][np.where(data['fofID']==fofID)[:][0]] = cd.lookback_time(data_ev_fofID['redshift'][len(data_ev_fofID['redshift'])-1], z0=0.0, **fidcosmo)/3.1536e+16 #seconds to Gyrs
                                        
                    data['delta_tvar_'+key][np.where(data['fofID']==fofID)[:][0]] = data['tmax_'+key][np.where(data['fofID']==fofID)[:][0]] - data['tmin_'+key][np.where(data['fofID']==fofID)[:][0]]

            except:
                pass
        print 'Nr missing haloid=0 at z=0:', count_haloid0_missing, 'Nr missing fofID at z=0:', len(array_lost_fofIdz0), '\nIDs:', array_lost_fofIdz0
                
        data['delta_t50_st2ha'] = data['t50_halo'] - data['t50_stars']
        data['flag_t50_st2ha'] = 1
        data['flag_t50_st2ha'][np.where(data['delta_t50_st2ha']<0.0)[:][0]]=0
        
        data['delta_t50_bh2ha'] = data['t50_halo'] - data['t50_bh']
        data['flag_t50_bh2ha'] = 1
        data['flag_t50_bh2ha'][np.where(data['delta_t50_bh2ha']<0.0)[:][0]]=0
        
        print data[['haloid', 'envr', 'mstar_30kpc', 't50_halo', 't50_stars', 'delta_t50_st2ha', 'delta_tvar_rVmax', 'tmin_angM_stars']][0:100]
        
        print '--> DONE!\n---------------------------------'

                
    elif key=='EAGLE_full_centrals':
   
        names=['fofID', 'haloid', 'hostid', 'galaxyID', 'nSubhalos',  'mhalo_fof', 'mhalo_200c', 'r200c',\
               'x_pos_subhalo', 'y_pos_subhalo', 'z_pos_subhalo', 'mPart_DM', 'rhalf_DM', 'rhalf_DM_2D', 'mPart_stars', 'mstar_30kpc',\
               'rhalf_stars',  'rhalf_stars_2D', 'mPart_gas', 'rhalf_gas', 'rhalf_gas_2D', 'mPart_bh', 'rhalf_bh', 'rhalf_bh_2D',\
               'vmax', 'rVmax', 'x_pos', 'y_pos', 'z_pos', 'x_vel', 'y_vel', 'z_vel', 'spinGas_x', 'spinGas_y', 'spinGas_z',\
               'mbh', 'vdisp', 'mstar_birth', 'sfr_total', 'Etot',  'Ekin', 'Etherm']
            
        myfilename2cross=mycomp+'anaconda/pro/data/EAGLE_100Mpc/EAGLE_100Mpc_z_0.0_tarsel_centrals.txt'
        
        data2cross = pd.read_csv(myfilename2cross,\
                                 names=names, skiprows=2, sep='  ')
        data2cross = df_to_sarray(data2cross)
        
        #print np.info(data2cross)
        
        #print data2cross[['fofID', 'haloid', 'mstar_30kpc']]
        print 'CROSSMATCH key: "', key, '"\n-----------------\nparent catalog:', data.size, 'data2cross:', data2cross.size,

        test = np.in1d(data['fofID'],data2cross['fofID'])

        #data2cross=data2cross[test]
        #data2cross[props[2::]]
        print '-->', test[np.where(test==True)].size, 'found to crossmatch!'
        #exit()
        #set all columns of the parent data to -99, so only the properties of the galaxies which are exsist in the sample
        #which should be crossmatched will be filled! But do not set the 'fofID' to -99 or the properties which are in both catalogs!!!
        

        new_props=['haloid', 'hostid', 'galaxyID', 'nSubhalos',  'mhalo_fof',\
                   'x_pos_subhalo', 'y_pos_subhalo', 'z_pos_subhalo', 'mPart_DM', 'rhalf_DM', 'rhalf_DM_2D', 'mPart_stars', 'mstar_30kpc',\
                   'rhalf_stars',  'rhalf_stars_2D', 'mPart_gas', 'rhalf_gas', 'rhalf_gas_2D', 'mPart_bh', 'rhalf_bh', 'rhalf_bh_2D',\
                   'vmax', 'rVmax', 'x_pos', 'y_pos', 'z_pos', 'x_vel', 'y_vel', 'z_vel',\
                   'vdisp', 'mstar_birth', 'sfr_total']
        data[new_props]=-99       
       
        parentIndices, index1, index2 = np.intersect1d(data['fofID'], data2cross['fofID'], return_indices=True) 
        
        
        for prop in new_props:
            print 'property:', prop
            data[prop][index1]=data2cross[prop][index2]
                       
        #print data[['jsub','fofID', 'haloid', 'galaxyID', 'mstar_1.5ropt', 'mhalo_200c', 'DvT_gas_c0.4','DvT_gas_c0.5', 'mstar_30kpc', 'x_pos_subhalo', 'y_pos_subhalo', 'z_pos_subhalo']][0:50]
           
        print '--> DONE!\n---------------------------------'
        
    elif key.startswith('EAGLE_AP'):

        ap=key[key.find('AP')+2::]
        print 'aperture is:', ap, 'pkpc'
  
        names=['fofID','haloid', 'hostid', 'galaxyID',
                'x_pos_subhalo', 'y_pos_subhalo', 'z_pos_subhalo',\
                'mhalo_'+ap+'kpc','mstar_'+ap+'kpc','mgas_'+ap+'kpc','sfr_'+ap+'kpc','vdisp_'+ap+'kpc']
            
        myfilename2cross=mycomp+'anaconda/pro/data/EAGLE_100Mpc/original_catalogs/EAGLE_100Mpc_full_catalog_ap'+ap+'pkpc.dat'
        
        data2cross = pd.read_csv(myfilename2cross,\
                                 names=names, skiprows=2, sep=' ')
        data2cross = df_to_sarray(data2cross)
     
        #print data2cross[['fofID', 'haloid', 'mstar_30kpc']]
        print 'CROSSMATCH key: "', key, '"\n-----------------\nparent catalog:', data.size, 'data2cross:', data2cross.size,

        test = np.in1d(data['fofID'],data2cross['fofID'])

        #data2cross=data2cross[test]
        #data2cross[props[2::]]
        print '-->', test[np.where(test==True)].size, 'found to crossmatch!'
        #exit()
        #set all columns of the parent data to -99, so only the properties of the galaxies which are exsist in the sample
        #which should be crossmatched will be filled! But do not set the 'fofID' to -99 or the properties which are in both catalogs!!!       

        new_props=['mstar_'+ap+'kpc','mgas_'+ap+'kpc']
        data[new_props]=-99       
        
        parentIndices, index1, index2 = np.intersect1d(data['fofID'], data2cross['fofID'], return_indices=True) 
        
        
        for prop in new_props:
            #print 'property:', prop
            data[prop][index1]=data2cross[prop][index2]
                       
        print data[['fofID','mstar_30kpc','mstar_50kpc']][0:50]
           
        print '--> DONE!\n-----------------\n\n'        

    elif key.startswith('EAGLE_SFgas'):
  
        names=['fofID','haloid', 'hostid', 'galaxyID',\
               'mPart_gas', 'mPart_stars', 'mstar_30kpc', 'mgas_30kpc', 'mPart_gasSF', 'mPart_gasNSF',\
               'Etherm_gasSF', 'Etherm_gasNSF', 'spinGasSF_x', 'spinGasSF_y', 'spinGasSF_z',\
               'metalfrac_SF', 'metalfrac_NSF', 'zgasSF_H', 'zgasSF_O']
            
        myfilename2cross=mycomp+'anaconda/pro/data/EAGLE_100Mpc/original_catalogs/EAGLE_100Mpc_full_catalog_gas+star-masses.dat'
        
        data2cross = pd.read_csv(myfilename2cross,\
                                 names=names, skiprows=2, sep=' ')
        data2cross = df_to_sarray(data2cross)
     
        #print data2cross[['fofID', 'haloid', 'mstar_30kpc']]
        print 'CROSSMATCH key: "', key, '"\n-----------------\nparent catalog:', data.size, 'data2cross:', data2cross.size,

        test = np.in1d(data['fofID'],data2cross['fofID'])

        #data2cross=data2cross[test]
        #data2cross[props[2::]]
        print '-->', test[np.where(test==True)].size, 'found to crossmatch!'
        #exit()
        #set all columns of the parent data to -99, so only the properties of the galaxies which are exsist in the sample
        #which should be crossmatched will be filled! But do not set the 'fofID' to -99 or the properties which are in both catalogs!!!       

        #new_props=['mgas_30kpc', 'mPart_gasSF', 'zgasSF_H', 'zgasSF_O']
        new_props=['mPart_gasSF', 'zgasSF_H', 'zgasSF_O']
        data[new_props]=-99       
        
        parentIndices, index1, index2 = np.intersect1d(data['fofID'], data2cross['fofID'], return_indices=True) 
        
        
        for prop in new_props:
            #print 'property:', prop
            data[prop][index1]=data2cross[prop][index2]
                                  
        print '--> DONE!\n-----------------\n\n'   


    elif key=='Silvio':
   
        names=['haloid','hostid','mPart_stars','rhalf_stars_2D',\
               'vdisp_r502D_edgeOn', 'VvS_r502D_edgeOn', 'vdisp_r502D_random', 'VvS_r502D_random',\
               'ell_r502D_edgeOn','ell_r502D_random', 'BvT_Lagos17b', 'Sersic_n',\
               'lambda_r502D_edgeOn','lambda_r502D_random',\
               'Mgas_HI_60kpc_GK11', 'Mgas_HI_60kpc_K13', 'Mgas_H2_60kpc_GK11', 'Mgas_H2_60kpc_K13',\
               'vdisp_2r502D_edgeOn', 'VvS_2r502D_edgeOn', 'vdisp_2r502D_random', 'VvS_2r502D_random', \
               'ell_2r502D_edgeOn', 'ell_2r502D_random',\
               'lambda_2r502D_edgeOn', 'lambda_2r502D_random', 'age_stars_rband_r502D', 'age_stars_rband_2r502D',\
               'sfr_30kpc', 'mhalo_unknown','lambda_r502D_random_ltc', 'lambda_2r502D_random_ltc', 'sfr_ltc', 'ssfr_ltc', 'mstar_ltc',\
               'VvS_r502D_edgeOn_ltc', 'VvS_r502D_random_ltc', 'VvS_2r502D_edgeOn_ltc', 'VvS_2r502D_random_ltc', 'lastTimeCentral',\
               'x_pos_subhalo', 'y_pos_subhalo', 'z_pos_subhalo', 'x_vel', 'y_vel', 'z_vel',\
               'rcenter2r200c', 'sDensity_N10', 'sDensity_N7', 'sDensity_N5', 'lastMerger', 'ratio_lastM', 'rhalf_stars', 'rhalf_sfr']
            
        myfilename2cross=mycomp+'anaconda/pro/data/EAGLE_100Mpc/original_catalogs/EAGLE_L100N1504_ForSilvio_format_Topcat.dat'
        
        data2cross = pd.read_csv(myfilename2cross,\
                                 names=names, skiprows=2, sep=' ')
            
        data2cross = df_to_sarray(data2cross)

        print 'CROSSMATCH key: "', key, '"\n-----------------\nparent catalog:', data.size, 'data2cross:', data2cross.size,
        
        data2cross['haloid'] = data2cross['haloid'].astype(np.int64)
        
        #data2cross=data2cross[np.where(data2cross['hostid']==0)[0][:]]
        #data=data[np.where(data['hostid']==0)[0][:]]
        test = np.in1d(data['haloid'],data2cross['haloid'])
        #data2cross[props[2::]]
        print '-->', test[np.where(test==True)].size, 'found to crossmatch!'
        
        #set all columns of the parent data to -99, so only the properties of the galaxies which are exsist in the sample
        #which should be crossmatched will be filled! But do not set the 'fofID' to -99 or the properties which are in both catalogs!!!
       

        # new_props=['lambda_r502D_edgeOn', 'lambda_2r502D_edgeOn',\
        #            #'age_stars_rband_r502D', 'age_stars_rband_2r502D',\
        #            'lastMerger', 'ratio_lastM']
            
            
        new_props=['lastMerger', 'ratio_lastM']
        data[new_props]=-99       
        
        parentIndices, index1, index2 = np.intersect1d(data['haloid'], data2cross['haloid'], return_indices=True) 
        
        
        for prop in new_props:
            #print 'property:', prop
            data[prop][index1]=data2cross[prop][index2]
            
            
        # data['delta_age_stars_rband'] = data['age_stars_rband_r502D']-data['age_stars_rband_2r502D']
        # data['delta_lambda_edgeOn'] = data['lambda_2r502D_edgeOn']-data['lambda_r502D_edgeOn']
        
        # for item in np.arange(data.size):
        #     if data[item]['age_stars_rband_r502D']==-99 and data[item]['age_stars_rband_2r502D']==-99:
        #         data[['delta_age_stars_rband']][item]=-99

        # for item in np.arange(data.size):
        #     if data[item]['lambda_2r502D_edgeOn']==-99 and data[item]['lambda_r502D_edgeOn']==-99:
        #         data[['delta_lambda_edgeOn']][item]=-99                
                
        
        for item in np.arange(data.size):
            #print item, data[['lastMerger','ratio_lastM','flag_majorM']][item],
            if data[item]['ratio_lastM']>0.0 and data[item]['ratio_lastM']<=0.2:
                data[item]['flag_majorM']=0
            elif data[item]['ratio_lastM']>0.2:
                data[item]['flag_majorM']=1
            else:
                data[['lastMerger','ratio_lastM','flag_majorM']][item]=-99
                
            #print '-->', data[['lastMerger','ratio_lastM','flag_majorM']][item]
                
        #exit()
        
        
        print '--> DONE!\n-----------------\n\n'

    elif key=='LSBGs':
 
        myfilename2cross=mycomp+'anaconda/pro/data/EAGLE_100Mpc/EAGLE_100Mpc_z_0.0_tarsel_centrals_crossmatched_bins_regplot.txt'
        names=['fofID','haloid','hostid','galaxyID','jsub',\
               'tmin_angM_stars','tmax_angM_stars','delta_tvar_angM_stars','tmin_rVmax',\
               'lastMinorM','lastMajorM','tmax_rVmax','mhalo_200c','vmax','NFW_con',\
               'ellDM_30kpc_T19','vdisp_30kpc_T19','Vrot2Vdisp_30kpc_T19','v200c','envr','r200c','rVmax','r2voidCenter',\
               'ropt','reff_gas_1.5ropt','reff_gas_disk_1.5ropt','rhalf_stars_1.5ropt',\
               'age_stars_rband_r502D','age_stars_rband_2r502D','delta_age_stars_rband','angM_stars','angM_SFgas','angM_NSFgas',\
               't50_stars','t70_stars','t50_halo','t70_halo','t50_bh','t70_bh',\
               'delta_tform_stars','delta_tform_halo','delta_tform_bh','delta_t50_st2ha','delta_t50_bh2ha',\
               'mbh','bheff','rhalf_stars_30kpc','DvT_stars_c0.5_1.5ropt','DvT_gas_c0.5_1.5ropt',\
               'SHMR_1.5ropt','Sigma_gas_1.5ropt','Sigma_gas_reff_disk_1.5ropt','Sigma_stars_1.5ropt','Sigma_sfr_1.5ropt',\
               'sfe_gas_reff_disk_1.5ropt','sfe_gas_1.5ropt','SB_mu_1.5opt_B','SB_mu_eff_stars_1.5ropt_r',\
               'mgas_1.5ropt','sfr_1.5ropt','ssfr_1.5ropt','OH_gas_30kpc','z_grad','Tcons_1.5ropt','cgf_1.5ropt','mPart_stars',\
               'MAB_dA_SDSS_g','MAB_dA_SDSS_r','MAB_dA_SDSS_i','MAB_dA_Johnson_B','MAB_dA_Johnson_V',\
               'fmol_30kpc_GK11','flag_SHMR_bin','mAB_dA_total_cut_g_i','mAB_dA_total_cut_B_V',\
               'flag_SB_B','flag_SB_r','flag_SFE','flag_mstar_bin','flag_mhalo_bin','flag_sfe_bin','flag_age_bin',\
               'angM_norm_1.5ropt_stars','angM_norm_1.5ropt_gas','lambda_r502D_edgeOn','lambda_2r502D_edgeOn','delta_lambda_edgeOn',\
               'flag_t50_st2ha','flag_t50_bh2ha','Mgas_HI_30kpc_GK11','Mgas_H2_30kpc_GK11','fatom_30kpc_GK11',\
               'mstar_30kpc','mstar_1.5ropt','mstar_half_reff_gas_disk','mgas_half_reff_gas_disk',\
               'sfe_gas_reff_1.5ropt','Sigma_gas_reff_1.5ropt','delta_tvar_rVmax']
        data2cross = pd.read_csv(myfilename2cross,\
                                 names=names, skiprows=2, sep='  ')
            

            
        data2cross = df_to_sarray(data2cross)
        
        #print np.info(data2cross)
        
        print 'CROSSMATCH key: "', key, '"\n-----------------\nparent catalog:', data.size, 'data2cross:', data2cross.size,



        data2cross['mAB_dA_total_cut_B_V']=data2cross['MAB_dA_Johnson_B']-data2cross['MAB_dA_Johnson_V']
        data2cross['mAB_dA_total_cut_g_i']=data2cross['MAB_dA_SDSS_g']-data2cross['MAB_dA_SDSS_i']
        
        new_props=['jsub', 'SB_mu_1.5opt_B','rhalf_stars_30kpc', 'MAB_dA_SDSS_r','MAB_dA_Johnson_B', 'ropt', 'rhalf_stars_1.5ropt',\
                   'NFW_con', 'Mgas_HI_30kpc_GK11','Mgas_H2_30kpc_GK11','fatom_30kpc_GK11','angM_norm_1.5ropt_stars','angM_norm_1.5ropt_gas',\
                   'z_grad','lambda_r502D_edgeOn','lambda_2r502D_edgeOn','delta_lambda_edgeOn','mAB_dA_total_cut_g_i','mAB_dA_total_cut_B_V', 'v200c',\
                   'DvT_stars_c0.5_1.5ropt','DvT_gas_c0.5_1.5ropt', 'Sigma_gas_1.5ropt','Sigma_gas_reff_disk_1.5ropt','Sigma_stars_1.5ropt','Sigma_sfr_1.5ropt',\
                   'sfe_gas_reff_disk_1.5ropt','sfe_gas_1.5ropt']
            
        new_props=['jsub', 'SB_mu_1.5opt_B','rhalf_stars_30kpc', 'MAB_dA_SDSS_r','MAB_dA_Johnson_B', 'ropt', 'rhalf_stars_1.5ropt',\
                   'NFW_con', 'Mgas_HI_30kpc_GK11','Mgas_H2_30kpc_GK11','fatom_30kpc_GK11','angM_norm_1.5ropt_stars','angM_norm_1.5ropt_gas',\
                   'z_grad','lambda_r502D_edgeOn','lambda_2r502D_edgeOn','delta_lambda_edgeOn','mAB_dA_total_cut_g_i','mAB_dA_total_cut_B_V', 'v200c',\
                   'DvT_stars_c0.5_1.5ropt','DvT_gas_c0.5_1.5ropt', 'Sigma_gas_1.5ropt','Sigma_gas_reff_disk_1.5ropt','Sigma_stars_1.5ropt','Sigma_sfr_1.5ropt',\
                   'sfe_gas_reff_disk_1.5ropt','sfe_gas_1.5ropt']
            
        data[new_props]=-99       
        
        parentIndices, index1, index2 = np.intersect1d(data['haloid'], data2cross['haloid'], return_indices=True) 
        #print data2cross[new_props][0:25]
        
        for prop in new_props:
            #print 'property:', prop
            data[prop][index1]=data2cross[prop][index2]        
        
        print '--> DONE!\n-----------------\n\n'
        
        #exit()

    elif key=='t50':
        mypath=mycomp+'anaconda/pro/data/EAGLE_100Mpc/EAGLE_100Mpc_z_0.0_tarsel_tree_roots_history_fit_assessment_good_trees_consistent_HF-PaperII.txt'

        print mypath
        data2cross = pd.read_csv(mypath, skiprows=2,\
                          names=['fofID', 'haloid', 'galaxyID', 'jsub', 'topLeafID', 'envr', 'mhalo_200c',\
                                 'r200c', 'ropt', 'rhalf_stars_1.5ropt',\
                                 't50_stars', 't70_stars', 't50_halo', 't70_halo', 't50_bh', 't70_bh',\
                                 'delta_tform_stars', 'delta_tform_halo', 'delta_tform_bh',\
                                 'delta_t50_st2ha', 'delta_t50_bh2ha', 't50_stars_std',\
                                 't70_stars_std', 't50_halo_std', 't70_halo_std', 't50_bh_std', 't70_bh_std',\
                                 'flag_std_fit_t70_stars', 'flag_std_fit_t70_halo', 'flag_std_fit_t70_bh',\
                                 'flag_std_fit_t50_stars', 'flag_std_fit_t50_halo', 'flag_std_fit_t50_bh',\
                                 'mbh', 'mstar_30kpc', 'SB_mu_1.5opt_B', 'SB_mu_eff_stars_1.5ropt_r', 'MAB_dA_SDSS_r', 'MAB_dA_Johnson_B',\
                                 'lastMerger', 'ratio_lastM', 'rVoid', 'DtoVoidCent', 'flag_majorM',\
                                 'flag_t50_st2ha', 'flag_t50_bh2ha', 'flag_SB_B', 'flag_SB_r',\
                                 'flag_use_inter_stars', 'flag_use_inter_halo', 'flag_use_inter_bh',\
                                 'flag_use_fit_stars', 'flag_use_fit_halo', 'flag_use_fit_bh',\
                                 'nr_strikes_stars', 'nr_strikes_halo', 'nr_strikes_bh',\
                                 'delta_t50_bh2st', 'flag_t50_bh2st', 'delta_t70_st2ha', 'delta_t70_bh2ha', 'delta_t70_bh2st',\
                                 'flag_t70_st2ha', 'flag_t70_bh2ha', 'flag_t70_bh2st','rhalf_stars_30kpc','NFW_con',\
                                 'Mgas_HI_30kpc_GK11','Mgas_H2_30kpc_GK11', 'z_grad','fatom_30kpc_GK11', 'fmol_30kpc_GK11',\
                                 'angM_norm_1.5ropt_stars','angM_norm_1.5ropt_gas','lambda_r502D_edgeOn','lambda_2r502D_edgeOn','delta_lambda_edgeOn',\
                                 'mAB_dA_total_cut_g_i','mAB_dA_total_cut_B_V', 'v200c',\
                                 'DvT_stars_c0.5_1.5ropt','DvT_gas_c0.5_1.5ropt',\
                                  'Sigma_gas_1.5ropt','Sigma_gas_reff_disk_1.5ropt','Sigma_stars_1.5ropt','Sigma_sfr_1.5ropt',\
                                 'sfe_gas_reff_disk_1.5ropt','sfe_gas_1.5ropt'
                                 ],\
                                 sep='  ') 
        data2cross = df_to_sarray(data2cross)
        

        print '\nLOAD HALF-MASS ASSEMBLY TIMES & CROSSMATCH WITH CATALOG AT z=0!\n-----------------\nparent catalog:',\
        data.size, 'data2cross:', data2cross.size
        #print np.info(data2cross)
        
        # new_props=['t50_stars', 't70_stars', 't50_halo', 't70_halo', 't50_bh', 't70_bh', 'delta_tform_stars', 'delta_tform_halo', 'delta_tform_bh',\
        #            'delta_t50_st2ha', 'delta_t50_bh2ha', 'delta_t50_bh2st', 'flag_t50_st2ha', 'flag_t50_bh2ha', 'flag_t50_bh2st',\
        #            'delta_t70_st2ha', 'delta_t70_bh2ha', 'delta_t70_bh2st',\
        #            'flag_t70_st2ha', 'flag_t70_bh2ha', 'flag_t70_bh2st',\
        #            'flag_SB_B', 'flag_SB_r', 'envr', 'lastMerger', 'ratio_lastM', 'rVoid', 'DtoVoidCent', 'flag_majorM','SB_mu_1.5opt_B', 'SB_mu_eff_stars_1.5ropt_r']     

        new_props=['t50_stars', 't50_halo', 't50_bh', 'delta_tform_stars', 'delta_tform_halo', 'delta_tform_bh',\
                   'delta_t50_st2ha', 'delta_t50_bh2ha', 'delta_t50_bh2st', 'flag_t50_st2ha', 'flag_t50_bh2ha', 'flag_t50_bh2st',\
                   'flag_SB_B', 'envr', 'lastMerger', 'ratio_lastM', 'rVoid', 'DtoVoidCent', 'flag_majorM',\
                   'z_grad', 'angM_norm_1.5ropt_stars', 'angM_norm_1.5ropt_gas', 'lambda_r502D_edgeOn', 'lambda_2r502D_edgeOn', 'delta_lambda_edgeOn',\
                    'DvT_stars_c0.5_1.5ropt', 'DvT_gas_c0.5_1.5ropt', 'Mgas_HI_30kpc_GK11', 'Mgas_H2_30kpc_GK11', 'fatom_30kpc_GK11', 'fmol_30kpc_GK11',\
                   'mAB_dA_total_cut_B_V', 'mAB_dA_total_cut_g_i', 'NFW_con', 'rhalf_stars_30kpc']     


        data[new_props]=-99       
        
        parentIndices, index1, index2 = np.intersect1d(data['haloid'], data2cross['haloid'], return_indices=True) 
        
        
        for prop in new_props:
            #print 'property:', prop
            data[prop][index1]=data2cross[prop][index2]
        
        print '--> DONE!\n-----------------\n\n'  
     

    elif key=='Patricia_reff_gas_disk':
 
        names=['jsub','reff_gas_disk']     
 
        myfilename2cross=mycomp+'anaconda/pro/data/EAGLE_100Mpc/original_catalogs/Doris_tabla_complemento.txt'
        
        data2cross = pd.read_csv(myfilename2cross,\
                                 names=names, skiprows=2, sep=' ')
            
        data2cross = df_to_sarray(data2cross)
        #print 'here: 3860', np.info(data2cross)
        data2cross['jsub'] = data2cross['jsub'].astype(np.int64)
                
        print 'CROSSMATCH key: "', key, '"\n-----------------\nparent catalog:', data.size, 'data2cross:', data2cross.size,

        test = np.in1d(data['jsub'],data2cross['jsub'])
        #data2cross[props[2::]]
        print '-->', test[np.where(test==True)].size, 'found to crossmatch!'
        
        #set all columns of the parent data to -99, so only the properties of the galaxies which are exsist in the sample
        #which should be crossmatched will be filled! But do not set the 'fofID' to -99 or the properties which are in both catalogs!!!    

        data['reff_gas_disk']=-99       
        
        parentIndices, index1, index2 = np.intersect1d(data['jsub'], data2cross['jsub'], return_indices=True) 
        
        data['reff_gas_disk'][index1]=data2cross['reff_gas_disk'][index2]
        
        
        print '--> DONE!\n-----------------\n\n'  


    elif key=='Patricia2':
 
        names=['jsub','np_disk','np_stars_1.5ropt','np_gas_1.5ropt','sfr_1.5ropt','mstar_1.5ropt','DvT_stars_c0.4','DvT_gas_c0.5','ropt','mhalo_200c','r200c',\
               'ssfr_1.5ropt','rhalf_stars_1.5ropt', 'mbh', 'bh_acc_rate', 'haloid', 'DvT_stars_c0.5', 'sfr_int_gas',\
               'reff_gas', 'angM_norm_1.5ropt_stars', 'angM_norm_1.5ropt_gas', 'x_pos', 'y_pos', 'z_pos']     
 
        myfilename2cross=mycomp+'anaconda/pro/data/EAGLE_100Mpc/original_catalogs/REF25_generaltable_z0_new.dat'
        
        data2cross = pd.read_csv(myfilename2cross,\
                                 names=names, skiprows=2, sep=' ')
            
        data2cross = df_to_sarray(data2cross)
        #print 'here: 3860', np.info(data2cross)
        data2cross['haloid'] = data2cross['haloid'].astype(np.int64)
        
        #Correct Mpc to kpc
        data2cross['rhalf_stars_1.5ropt']*=1000.0
        data2cross['ropt']*=1000.0
        data2cross['reff_gas']*=1000.0
        
        print 'CROSSMATCH key: "', key, '"\n-----------------\nparent catalog:', data.size, 'data2cross:', data2cross.size,

        test = np.in1d(data['haloid'],data2cross['haloid'])
        #data2cross[props[2::]]
        print '-->', test[np.where(test==True)].size, 'found to crossmatch!'
        
        #set all columns of the parent data to -99, so only the properties of the galaxies which are exsist in the sample
        #which should be crossmatched will be filled! But do not set the 'fofID' to -99 or the properties which are in both catalogs!!!

        new_props=['reff_gas']     

        data[new_props]=-99       
        
        parentIndices, index1, index2 = np.intersect1d(data['haloid'], data2cross['haloid'], return_indices=True) 
        
        
        for prop in new_props:
            #print 'property:', prop
            data[prop][index1]=data2cross[prop][index2]
        
        print '--> DONE!\n-----------------\n\n'       

    elif key=='Patricia3':
        """
        #Data from 27/06/2023
            properties:
                isub (internal indexx patricia), ndisks (nstars in the disc),nopts(nstars at roptical-83%), noptg(ngas particles at roptical), 
                sfr(stas <2Gyr within ropt),
                total(mstaropt), total(mgassopt), DTstars(at roptical, elim=0.4), DTgas(at roptical, elim=0.5),
                roptical, mvirial,rvirial,sSFR (sfr/total(mstaropt)),
                refft(half-stellar mass ratio within roptical),
                MblackH(isub),MblackHAccretion(isub),IDhalo(isub),
                DTstars1(elim=0.5),sfrgastotal(estimated with all SF gas),reffg(half gas mass radius within roptical),
                jmodgas,jmodstars,xcm,ycm,zcm,
                mhalft(half of the stellar mass when the total stellar mass is defined within 1.5 x optical radius),
                reffg_disc (effective radius (half-mass radius) of the disk component when the corresponding gas mass was defined within 1.5 x optical radius,
                mhalft_disc (half of the gas mass within reffg_disc) ,
                mstar_reffg_disc (is the stellar mass within reffg_disc)
                
           notes:   when it says roptical it is actual 1.5 roptical)
                    radius are in Mpc, h is already included.
        """
        names=['jsub','np_disk','np_stars_1.5ropt','np_gas_1.5ropt','sfr_1.5ropt','mstar_1.5ropt', 'mgas_1.5ropt',\
               'DvT_stars_c0.4_1.5ropt','DvT_gas_c0.5_1.5ropt','ropt','mhalo_200c','r200c',\
               'ssfr_1.5ropt','rhalf_stars_1.5ropt', 'mbh', 'bh_acc_rate', 'haloid', 'DvT_stars_c0.5_1.5ropt', 'sfr_gasSF',\
               'reff_gas_1.5ropt', 'angM_norm_1.5ropt_stars', 'angM_norm_1.5ropt_gas', 'x_pos', 'y_pos', 'z_pos',\
               'mgas_half_1.5ropt', 'reff_gas_disk_1.5ropt','mgas_half_reff_gas_disk','mstar_half_reff_gas_disk']     
 
        myfilename2cross=mycomp+'anaconda/pro/data/EAGLE_100Mpc/original_catalogs/REF25_generaltable_z0_27062023.dat'
        
        data2cross = pd.read_csv(myfilename2cross,\
                                 names=names, skiprows=2, sep=' ')
            
        data2cross = df_to_sarray(data2cross)
        #print 'here: 3860', np.info(data2cross)
        data2cross['haloid'] = data2cross['haloid'].astype(np.int64)
        
        #Correct Mpc to kpc
        data2cross['rhalf_stars_1.5ropt']*=1000.0
        data2cross['ropt']*=1000.0
        data2cross['reff_gas_1.5ropt']*=1000.0
        data2cross['reff_gas_disk_1.5ropt']*=1000.0
        
        print 'CROSSMATCH key: "', key, '"\n-----------------\nparent catalog:', data.size, 'data2cross:', data2cross.size,

        test = np.in1d(data['haloid'],data2cross['haloid'])
        #data2cross[props[2::]]
        print '-->', test[np.where(test==True)].size, 'found to crossmatch!'
        
        #set all columns of the parent data to -99, so only the properties of the galaxies which are exsist in the sample
        #which should be crossmatched will be filled! But do not set the 'fofID' to -99 or the properties which are in both catalogs!!!

        # new_props=['jsub','np_disk','np_stars_1.5ropt','np_gas_1.5ropt',\
        #            'mstar_1.5ropt', 'mstar_half_reff_gas_disk', 'mgas_1.5ropt', 'mgas_half_1.5ropt', 'mgas_half_reff_gas_disk',\
        #            'DvT_stars_c0.4_1.5ropt', 'DvT_stars_c0.5_1.5ropt', 'DvT_gas_c0.5_1.5ropt',\
        #            'ropt', 'rhalf_stars_1.5ropt','reff_gas_1.5ropt', 'reff_gas_disk_1.5ropt',\
        #            'sfr_gasSF', 'sfr_1.5ropt', 'ssfr_1.5ropt',\
        #            'angM_norm_1.5ropt_stars', 'angM_norm_1.5ropt_gas', 'mbh', 'bh_acc_rate']
            
        # new_props=['jsub', 'ropt', 'mstar_1.5ropt', 'mgas_1.5ropt', 'rhalf_stars_1.5ropt','reff_gas_1.5ropt',\
        #            'sfr_1.5ropt']
        #            # 'angM_norm_1.5ropt_stars', 'angM_norm_1.5ropt_gas']              

        new_props=['jsub']#, 'ropt','rhalf_stars_1.5ropt']
        data[new_props]=-99       
        
        parentIndices, index1, index2 = np.intersect1d(data['haloid'], data2cross['haloid'], return_indices=True) 
        
        
        for prop in new_props:
            #print 'property:', prop
            data[prop][index1]=data2cross[prop][index2]
        
        print '--> DONE!\n-----------------\n\n'       

    elif key=='CMASS_new':
        #Crossmatch total catalog of Galacticus with original CMASS_down3 to get more propeties!
        
        import pandas as pd
        dt = dt_MDPL_Galacticus_CMASS_down_sample3_original()
                
        data2cross = pd.read_csv(mycomp+'anaconda/pro/data/Galacticus_1Gpc/SFH/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3.txt',\
                           names=[dt[c][0] for c in xrange(len(dt))], skiprows=2, sep='  ')
            
        data2cross = df_to_sarray(data2cross)
        indices = np.in1d(data['hostid'], data2cross['hostid'])
        
        print indices
        data=data[np.where(indices==True)[:][0]]
        
        print data.size
    
    elif key=='CMASS_new_envr':
        #Crossmatch new catalog with environment from original!

        filename_env=mycomp+'anaconda/pro/data/Galacticus_1Gpc/SFH/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3.txt'
        mydata_names_props= ['haloid','hostid','mstar','mhalo','orphan','mcold','Mzgas','Mzstar','mbh','mstar_spheroid','mstar_disk','mcold_disk',\
                              'mcold_spheroid','mhot','sfr_spheroid','sfr_disk','sfr','Mzgas_spheroid','Mzgas_disk','Mzstar_spheroid','Mzstar_disk',\
                              'Mzhot_halo','x_pos','y_pos','z_pos','x_vel','y_vel','z_vel','L_SDSS_dA_total_u','L_SDSS_dA_total_g','L_SDSS_dA_total_r',\
                              'L_SDSS_dA_total_i','L_SDSS_dA_total_z','rdisk','rbulge','rhalf_mass','spinParameter','MAB_dA_total_u','MAB_dA_total_g',\
                              'MAB_dA_total_r','MAB_dA_total_i','MAB_dA_total_z','mAB_dA_total_u','mAB_dA_total_g','mAB_dA_total_r','mAB_dA_total_i',\
                              'mAB_dA_total_z','mAB_dA_total_cut_r_i','mAB_dA_total_cut_dmesa','mAB_dA_total_cut_i_lt_dmesa','mhalo_sat','nodeIndex',\
                              'parentIndex','satelliteNodeIndex','satelliteIndex','siblingIndex','satelliteMergeTime','isolated','timeLastIsolated',\
                              'mhalo_200c','zcold','cgf','mbasic','mbasic_200c','NFW_con','ssfr','weight_tot','env_512','env_1024','pop']
        
        data2cross = df_to_sarray(pd.read_csv(filename_env, skiprows=2, names=mydata_names_props, sep='  '))       
        
        parentIndices, indices, index_before = np.intersect1d(data['hostid'], data2cross['hostid'], return_indices=True)
        
        data[['env_1024', 'env_512', 'pop','MAB_dA_total_g','MAB_dA_total_r','MAB_dA_total_i','mhot','Mzhot_halo','timeLastIsolated', 'NFW_con','mstar_spheroid']][indices]\
            =data2cross[['env_1024', 'env_512', 'pop','MAB_dA_total_g','MAB_dA_total_r','MAB_dA_total_i','mhot','Mzhot_halo','timeLastIsolated','NFW_con','mstar_spheroid']][index_before]
        

    elif key=='300MCG':
        #Crossmatch 300 Cluster Catalog "Properties" with Galacticus output        
        
        import pandas as pd
        props=['haloid', 'parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar',\
                'SHMR', 'mcold', 'Mzgas', 'zcold', 'g-i',\
                'sfr', 'ssfr', 'X', 'Y', 'Z']  
        
        myfile2cross=mycomp+'anaconda/pro/data/300/SFH_Properties/SFH_Properties_300_z_'+str("{0:.2f}".format(redshift))+'__props.txt'
        

        print 'myfile2cross:', myfile2cross
       
        data2cross = pd.read_csv(myfile2cross,\
                                 names=props, skiprows=2, sep='\t')
        data2cross = df_to_sarray(data2cross)

        #data2cross = data[0:324] #Uncomment when testing only!!!
        data2cross[::-1].sort(order=['mhalo'], axis=0)

        data2cross['SHMR']=data2cross['mstar']/data2cross['mhalo']
        
        SHMR_50=np.nanmedian(data2cross['SHMR'])
        SHMR_25=np.nanpercentile(data2cross['SHMR'], 25)
        SHMR_75=np.nanpercentile(data2cross['SHMR'], 75)      
        
        import numpy.lib.recfunctions as rcfuncs       
        data2cross = rcfuncs.append_fields([data2cross], ['DeltaSHMR','cgf','Tcons','clusterIndex','z50_mstar','z50_mhalo','DeltaSHMR_Seb','Delta_z50',\
                                                          'mbh','bheff','haloidz0','g-i_z50Mstar','g-i_z50Mhalo','Delta_g-iz50'], [np.zeros(data2cross.size,),np.zeros(data2cross.size,),np.zeros(data2cross.size,),data2cross['parentIndex'],np.zeros(data2cross.size,),np.zeros(data2cross.size,),np.zeros(data2cross.size,),np.zeros(data2cross.size,),np.zeros(data2cross.size,),np.zeros(data2cross.size,),data2cross['parentIndex'],np.zeros(data2cross.size,),np.zeros(data2cross.size,),np.zeros(data2cross.size,)], usemask=False)
        
        data2cross['DeltaSHMR']=data2cross['SHMR']-SHMR_50
        data2cross['cgf']=data2cross['mcold']/data2cross['mstar']
        data2cross['Tcons']=data2cross['mcold']/data2cross['sfr']/1e9
         
        #assign clusterIndices
        if redshift==0.0:
            #Set for the first time if z=0
            i=0        
            while i<data2cross.size:
                data2cross['clusterIndex'][i]=i+1                
                i+=1
                
            data2cross['haloidz0']=data2cross['haloid']
            
        else:
            #load from file for z>0
            clusterIds = pd.read_csv(mycomp+'/anaconda/pro/data/300/SFH_Indices/SFH_clusterIndices.txt',\
                                     names=['haloid','clusterIndex','haloidz0'], skiprows=2, sep='\t')
                
            parentIndices, index1, index2 = np.intersect1d(data2cross['parentIndex'], clusterIds['haloid'], return_indices=True) 
        
            data2cross[['clusterIndex','haloidz0']][index1]=clusterIds[['clusterIndex','haloidz0']][index2]

        #write formation times and DeltaSHMR calculated by Sebastian to data set which should be crossmatched!
        data_z50 = pd.read_csv(mycomp+'anaconda/pro/data/300/haloid_z50.txt',\
                                 names=['haloid','z50_mstar','z50_mhalo','DeltaSHMR_Seb'], skiprows=1, sep='  ')
        data_z50 = df_to_sarray(data_z50)  
        
        parentIndices, index1, index2 = np.intersect1d(data2cross['haloidz0'], data_z50['haloid'], return_indices=True) 
        
        data2cross[['z50_mstar','z50_mhalo','DeltaSHMR_Seb']][index1]=data_z50[['z50_mstar','z50_mhalo','DeltaSHMR_Seb']][index2]       
          
        data2cross['Delta_z50']=data2cross['z50_mhalo']-data2cross['z50_mstar']
        #print data2cross[['haloid','z50_mstar','z50_mhalo','DeltaSHMR_Seb','Delta_z50']][0:10]        
                                                         
                   
        #print data2cross[['haloid','hostid','parentIndex','clusterIndex']][0:10]
        
        print 'data to crossmatch:\n-----------------\nparent catalog:', data.size, 'data2cross:', data2cross.size,
        #look up all galaxies in from the set with should be crossmatched with the main data set        
        test = np.in1d(data['hostid'], data2cross['hostid'])
        #if that is true we include those indices in 'index1' too. In other words we have found all satellites to a certain halo                    
        data=data[np.where(test==True)[:][0]]     
 
        print '-->', data.size, 'found --> DONE!\n'
        #Assign properties to the galaxies which should be crossmatched!
        parentIndices, index1, index2 = np.intersect1d(data['hostid'], data2cross['hostid'], return_indices=True) 
        
        #Write formation times and DeltaSHMR values from Sebastian in into the main data set
        data[['parentIndex','clusterIndex','DeltaSHMR','z50_mstar','z50_mhalo','DeltaSHMR_Seb','Delta_z50','haloidz0']][index1]=data2cross[['haloid','clusterIndex','DeltaSHMR','z50_mstar','z50_mhalo','DeltaSHMR_Seb','Delta_z50','haloidz0']][index2]     
        
        #Write the properties which are missing to the 300-cluster data set
        data2cross[['mbh','bheff']][index2]=data[['mbh','bheff']][index1] 

        data_colorz50 = pd.read_csv(mycomp+'anaconda/pro/data/300/Colors_at_z50_and_other_properties.txt',\
                                 names=['haloid','z50_mstar','z50_mhalo','g-i_z50Mstar', 'g-i_z50Mhalo', 'dummy1', 'dummy2', 'dummy3', 'dummy4', 'dummy5' ,'dummy6','dummy7'], skiprows=1, sep='  ')
        data_colorz50 = df_to_sarray(data_colorz50)  

        parentIndices, index1, index2 = np.intersect1d(data2cross['haloid'], data_colorz50['haloid'], return_indices=True) 
        
        data2cross[['g-i_z50Mstar','g-i_z50Mhalo']][index1]=data_colorz50[['g-i_z50Mstar','g-i_z50Mhalo']][index2]
        
        data2cross['Delta_g-iz50']=data2cross['g-i_z50Mhalo']-data2cross['g-i_z50Mstar']

        myOutput = oD.OutputData(config=False)
        keyID='writeIndicesd'
        if keyID=='writeIndices':
                     
            #print data[['haloid','hostid','parentIndex','clusterIndex']][0:20] 
                                 
            myOutput.writeIntoFile(mycomp+'/anaconda/pro/data/300/SFH_Indices/SFH_clusterIndices_'+str(format(redshift, '0.2f'))+'.txt',
                       data[['haloid','parentIndex','clusterIndex','haloidz0','x_pos','y_pos','z_pos']],
                       myheader= 'The ThreeHundred (300) MD-Galacticus 1Gpc merger trees with clusterIDs z='+str(format(redshift, '0.2f'))+'\n(1) haloid (2) parentIndex (3) clusterIndex (4) haloid at z=0 (5) X [Mpc] (6) Y [Mpc] (7) Z [Mpc]',
                       data_format="%i\t%i\t%i\t%i\t%0.5f\t%0.5f\t%0.5f",
                       mydelimiter='\t')
            
            myOutput.writeIntoFile(mycomp+'/anaconda/pro/data/300/SFH_Indices/SFH_clusterIndices.txt',
                       data[['haloid','clusterIndex','haloidz0']],
                       myheader= 'last z='+str(format(redshift, '0.2f'))+'\n(1) haloid (2) clusterIndex (3) haloid at z=0',
                       data_format="%i\t%i\t%i",
                       mydelimiter='\t')

            
        #Copy data to produce a output with some columns as log10 (easier to handle in Topcat!)
        import copy
        data2cross_log_output = copy.deepcopy(data2cross)        
        
        for item in ['mhalo', 'mstar','SHMR', 'mcold', 'Mzgas', 'sfr','ssfr', 'mbh', 'bheff']:            
            data2cross_log_output[item]=np.log10(data2cross_log_output[item])
        
        data2cross_log_output = rcfuncs.append_fields([data2cross_log_output], ['z'], [np.zeros(data2cross.size,)], usemask=False)
        data2cross_log_output['z']=redshift
        
        myOutput.writeIntoFile(myfile2cross[:-4]+'_new_log10_test.txt',
                               data2cross_log_output[['haloid','parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar', 'SHMR', 'mcold', 'Mzgas', 'zcold', 'g-i', 'sfr', 'ssfr','X','Y','Z','DeltaSHMR','cgf','Tcons','clusterIndex','mbh','bheff','z50_mstar','z50_mhalo','DeltaSHMR_Seb','Delta_z50','g-i_z50Mstar','g-i_z50Mhalo', 'Delta_g-iz50']],
                               myheader='z='+str(format(redshift, '0.2f'))+' SHMR_50='+str("{0:.6f}".format(SHMR_50))+' 25th/75th='+str("{0:.6f}".format(SHMR_25))+'/'+str("{0:.6f}".format(SHMR_75))+'\n(1) haloid  (2) parentIndex (3) hostid (4) nodeIndex (5) satelliteNodeIndex (6) orphan '\
                                                                           '(7) log10(Mvir [Msun]) (8) log10(Mstar [Msun]) (9) log10(SHMR [-]) (10) log10(mcold [Msun]) (11) log10(Mzgas [Msun]) '\
                                                                           '(12) zcold [-] (13) g-i [-] (14) log10(sfr [Msun yr-1]) (15) log10(ssfr [yr-1]) (16) X [Mpc] (17) Y [Mpc] (18) Z [Mpc] '\
                                                                           '(19) Delta_SHMR (20) Mcold/Mstar [-] (21) Tcons (Mcold/sfr) [Gyr-1] (22) clusterIndex (23) log10(mbh [Msun]) (24) log10(bheff (mbh/mhalo) [-]) '\
                                                                           '(25) z50_mstar [-] (26) z50_mhalo [-] (27) DeltaSHMR_Seb [-] (28) Delta_z50 [-] (29) g-i_z50-mstar (30) g-i_z50-mhalo (31) Delta_g-iz50',
                               data_format="%i\t%i\t%i\t%i\t%i\t%i\t%0.5e\t%0.5e\t%0.5f\t%0.8e\t%0.8e\t%0.5f\t%0.5f\t%0.8e\t%0.8e\t%0.5f\t%0.5f\t%0.5f\t%0.6e\t%0.6e\t%0.6e\t%i\t%0.6e\t%0.6e\t%0.4f\t%0.4f\t%0.6e\t%0.4f\t%0.5f\t%0.5f\t%0.5f",
                               mydelimiter='\t')


        exit()
        #write all data into one file!
        if redshift==0.0:
            data2output=data2cross_log_output[['haloid','parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar', 'SHMR', 'mcold', 'Mzgas', 'zcold', 'g-i', 'sfr', 'ssfr','X','Y','Z','DeltaSHMR','cgf','Tcons','clusterIndex','mbh','bheff','z50_mstar','z50_mhalo','DeltaSHMR_Seb','Delta_z50']]
           
        else:
            data2output=np.append(data2output,data2cross_log_output[['haloid','parentIndex', 'hostid', 'nodeIndex', 'satelliteNodeIndex', 'orphan', 'mhalo', 'mstar', 'SHMR', 'mcold', 'Mzgas', 'zcold', 'g-i', 'sfr', 'ssfr','X','Y','Z','DeltaSHMR','cgf','Tcons','clusterIndex','mbh','bheff','z50_mstar','z50_mhalo','DeltaSHMR_Seb','Delta_z50']], axis=0)
            

        myOutput.writeIntoFile(mycomp+'/anaconda/pro/data/300/SFH_300-MDPL2-Galacticus_main_log10.txt',
                               data2output,
                               myheader='MDPL2-Galacticus "TheThreeHundred" central cluster galaxies, main progenitors: 0.00 < z < '+str(format(redshift, '0.2f'))+'\n(1) haloid  (2) parentIndex (3) hostid (4) nodeIndex (5) satelliteNodeIndex (6) orphan '\
                                                                           '(7) log10(Mvir [Msun]) (8) log10(Mstar [Msun]) (9) log10(SHMR [-]) (10) log10(mcold [Msun]) (11) log10(Mzgas [Msun]) '\
                                                                           '(12) zcold [-] (13) g-i [-] (14) log10(sfr [Msun yr-1]) (15) log10(ssfr [yr-1]) (16) X [Mpc] (17) Y [Mpc] (18) Z [Mpc] '\
                                                                           '(19) Delta_SHMR (20) Mcold/Mstar [-] (21) Tcons (Mcold/sfr) [Gyr-1] (22) clusterIndex (23) log10(mbh [Msun]) (24) log10(bheff (mbh/mhalo) [-]) '\
                                                                           '(25) z50_mstar [-] (26) z50_mhalo [-] (27) DeltaSHMR_Seb [-] (28) Delta_z50 [-]',
                               data_format="%i\t%i\t%i\t%i\t%i\t%i\t%0.5e\t%0.5e\t%0.5f\t%0.8e\t%0.8e\t%0.5f\t%0.5f\t%0.8e\t%0.8e\t%0.5f\t%0.5f\t%0.5f\t%0.6e\t%0.6e\t%0.6e\t%i\t%0.6e\t%0.6e\t%0.4f\t%0.4f\t%0.6e\t%0.4f",
                               mydelimiter='\t')

    elif key=='filter300':
        #Crossmatch 300 Cluster Catalog "Properties" with Galacticus output        
        
        cat='SAGE'
        
        import pandas as pd
        props=['regionName', 'regionX', 'regionY', 'regionZ', 'regionRadius', 'ngals', 'ncents', 'nsats', 'norphans',\
               'mhalo', 'X', 'Y', 'Z', 'haloid', 'offsetX', 'offsetY', 'offsetZ'] 
        
        myfile2cross=mycomp+'anaconda/pro/data/workshop_300/00_MDPL2_'+cat+'_z0.00_subsamples_log.txt'
        myfile2cross=mycomp+'anaconda/pro/data/workshop_300/00_MDPL2_'+cat+'_z0.00_subsamples_r_1h-1Mpc_log.txt'
        

        print 'myfile2cross:', myfile2cross
       
        data2cross = df_to_sarray(pd.read_csv(myfile2cross,\
                                 names=props, skiprows=8, sep='\t'))
        
        # print data2cross
        print np.info(data2cross)
        print data2cross[['regionName','haloid']]
                       
        # data2cross=choose_random_sample(data,324)
        # data2cross['regionName']=data2cross_test['regionName']
        # print data2cross.size
 
        print data[['haloid', 'hostid', 'CentralGalaxyIndex','GalaxyIndex', 'orphan', 'mhalo']][0:100]
        print data[np.where(data['haloid']==12695175512)[:][0]][['haloid', 'hostid', 'CentralGalaxyIndex','GalaxyIndex', 'orphan', 'mhalo']]
  
 
        test = np.in1d(data['hostid'], data2cross['haloid'])
        #if that is true we include those indices in 'index1' too. In other words we have found all satellites to a certain halo                    
        data=data[np.where(test==True)[:][0]]
        

        
 
        print '-->', data.size, 'found --> DONE!\n'
        #Assign properties to the galaxies which should be crossmatched!
        parentIndices, index1, index2 = np.intersect1d(data['hostid'], data2cross['haloid'], return_indices=True) 
        data['regionName'][index1]=data2cross['regionName'][index2]
        
        for regName in data2cross['regionName']:
            print regName, data[np.where(data['regionName']==regName)[:][0]][['haloid', 'hostid', 'CentralGalaxyIndex','GalaxyIndex', 'orphan', 'mhalo']]
        #print data2cross['regionName']
        
        #exit()
 
    else:

        import numpy.lib.recfunctions as rcfuncs
        data = rcfuncs.append_fields([data], ['env_1024'], [np.zeros(data.size,)], usemask=False)


        filename_env=mycomp+'anaconda/pro/data/Galacticus_1Gpc/SFH/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3.txt'
        mydata_names_props= ['haloid','hostid','mstar','mhalo','orphan','mcold','Mzgas','Mzstar','mbh','mstar_spheroid','mstar_disk','mcold_disk',\
                              'mcold_spheroid','mhot','sfr_spheroid','sfr_disk','sfr','Mzgas_spheroid','Mzgas_disk','Mzstar_spheroid','Mzstar_disk',\
                              'Mzhot_halo','x_pos','y_pos','z_pos','x_vel','y_vel','z_vel','L_SDSS_dA_total_u','L_SDSS_dA_total_g','L_SDSS_dA_total_r',\
                              'L_SDSS_dA_total_i','L_SDSS_dA_total_z','rdisk','rbulge','rhalf_mass','spinParameter','MAB_dA_total_u','MAB_dA_total_g',\
                              'MAB_dA_total_r','MAB_dA_total_i','MAB_dA_total_z','mAB_dA_total_u','mAB_dA_total_g','mAB_dA_total_r','mAB_dA_total_i',\
                              'mAB_dA_total_z','mAB_dA_total_cut_r_i','mAB_dA_total_cut_dmesa','mAB_dA_total_cut_i_lt_dmesa','mhalo_sat','nodeIndex',\
                              'parentIndex','satelliteNodeIndex','satelliteIndex','siblingIndex','satelliteMergeTime','isolated','timeLastIsolated',\
                              'mhalo_200c','zcold','cgf','mbasic','mbasic_200c','NFW_con','ssfr','weight_tot','env_512','env_1024','pop']
        
        mydata_props = df_to_sarray(pd.read_csv(filename_env, skiprows=2, names=mydata_names_props, sep='  '))       
        
        parentIndices, indices, index_before = np.intersect1d(data['hostid'], mydata_props['hostid'], return_indices=True)
        
        data['env_1024'][indices]=mydata_props['env_1024'][index_before]

        
        if redshift_str==0.56:            
            data_histo = np.zeros((55,72), dtype=np.float32)
        else:
            redshift_before=self.snap_array[self.myconfig_array[self.myconfig_array['catname'+str(self.myPipe.a)]+'_filename'+str(self.myPipe.i-1)]+'_snapid'+str(self.myPipe.i-1)]['z']
        
            filename_before=mycomp+'anaconda/pro/myRun/histos/sfr2z/Gal-dens_main_cents_'+method+'/xi/sfr2z_Galacticus_1Gpc_z_'+str(format(redshift_before, '0.2f'))+'_tarsel_SFH_down3_'+method+'_main_cents_'+sel+'_'+env_name+'.txt'
            #print 'filename before:', filename_before
            data_histo = myData.readAnyFormat(config=False, mypath=filename_before, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.float64, skiprow=2) 
        
        for env, env_name in zip([2,3],['f','k']):
            data_sel=data[np.where(data['env_1024']==env)[0][:]]  
        
        
            for sel in ['','low','passive']:#,'red','blue','low-zcold','high-zcold']:
                
                print 'selection', sel,
                method='M2'
                    
                filename=mycomp+'anaconda/pro/data/Galacticus_1Gpc/SFH/SFH_Index_'+str(format(redshift, '0.2f'))+'_'+sel+'_props_'+method+'_'+env_name+'.txt'
                print 'z:', str(format(redshift, '0.2f')), filename
                
                # (0) haloid  (2) parentIndex (1) hostid (3) nodeIndex (4) satelliteNodeIndex (5) orphan (6) Mvir [Msun] (7) Mstar [Msun] (8) SHMR [-] (9) mcold [Msun] (10) Mzgas [Msun] (11) zcold [-] (12) g-i [-] (13) sfr [Msun yr-1] (14) ssfr [yr-1]
                mydata_names=['haloid','parentIndex','hostid','nodeIndex','satelliteNodeIndex','orphan','Mvir','Mstar','SHMR','mcold','Mzgas','zcold',\
                              'g-i','sfr','ssfr']
                
                mydata_sel = df_to_sarray(pd.read_csv(filename, skiprows=2, names=mydata_names, sep='\t'))
                
                parentIndices, indices, index_before = np.intersect1d(data_sel['hostid'], mydata_sel['hostid'], return_indices=True)
        
                mydata=data_sel[indices]
                
        
        
                percent_low = 31.7
                percent_high = 68.3
        
                data_histo, header_bucket, data_format = fill_bucket_sfr2z(data_histo,
                                                                myfiltered_data,
                                                                data_offset,
                                                                redshift,
                                                                self.i,
                                                                self.a,
                                                                self.volume,
                                                                percent_low=percent_low,
                                                                percent_high=percent_hight)              
        
                
                filename=mycomp+'anaconda/pro/myRun/histos/sfr2z/Gal-dens_main_cents_'+method+'/xi/sfr2z_Galacticus_1Gpc_z_'+str(format(redshift, '0.2f'))+'_tarsel_SFH_down3_'+method+'_main_cents_'+sel+'_'+env_name+'.txt'
        
                #self.TwoPCF(mydata, filename, self.mycond_array)
                
                myOutput.writeIntoFile(filename,
                                        data_histo,
                                        myheader=sel+' selection key: '+element+' '\
                                                +self.myconfig_array[self.myconfig_array['catname'+str(self.a)]+'_simulation_name']+' '\
                                                +self.myconfig_array['catname'+str(self.a)]\
                                                +' error main progenitor search: not applied\n'\
                                                +header_bucket,
                                        data_format=data_format)
             
    

    return data
    
def split_sample(data):
    """
        Divide sample in two populations as published in paper Stoppacher+19 (SAM-CMASS, arXiv:1902.05496)
    """
    print 'divide sample!'
    y=np.log10(data['ssfr'])
    x=np.log10(data['sfr'])
    prop=(y+11.16)/1.12-x
    print prop
    data['pop'][np.where(prop<=0.0)]=2
    data['pop'][np.where(prop>0.0)]=1

    print 'min/max', min(data['pop']), max(data['pop'])      

    return data            

def organise_300Clusters(data):

    import pandas as pd
    
    for n_cl in list(range(1,324)):
        print 'cluster number:', n_cl,
            
        filename=mycomp+'anaconda/pro/data/300/Galacticus/MDPL2_Galacticus_z0.00_region'+str(n_cl).zfill(4)+'.txt'
        print filename
        mydata_names=['z','hostid','haloid','orphan','x_pos','y_pos','z_pos','x_vel','y_vel','z_vel','mstar_spheroid','mstar_disk',\
                      'mcold_spheroid','mcold_disk','mhot','mbh','sfr','sfr_spheroid','sfr_disk','mean_age_stars','mhalo','rest',\
                      'Vmax','Vpeak','NFW_con','spin','MZstar_spheroid','Mzstar_disk', 'Mzgas_spheroid', 'Mzgas_disk','mzhot_halo',\
                      'L_SDSS_dA_total_u','L_SDSS_dA_total_g','L_SDSS_dA_total_r','L_SDSS_dA_total_i','L_SDSS_dA_total_z',\
                      'MAB_dA_total_u','MAB_dA_total_g','MAB_dA_total_r','MAB_dA_total_i','MAB_dA_total_z',\
                      'rdisk','rbulge','rhalf_mass','mhot_outflow']
        
        mydata = df_to_sarray(pd.read_csv(filename, skiprows=2, names=mydata_names, sep='  '))                
        
        if n_cl==1:
            data=mydata
        else:
            data = np.append(data, mydata, axis=0)
            
        print 'size:', data.size
    
    try:
        if np.all(data['zcold']==-99.0):
            print '\n+++++++++++++++++++\ncalculate gas-phase metallicity',
            data['zcold']=8.69+np.log10(data['Mzgas']/(data['mcold']*0.0134))
            #data['OH_gas_disk_bulge']=12.0+np.log10(data['OH_gas_disk_bulge']/16.0)
            print '--> DONE!\n'
    except:
        print 'SAGE:',
        #data['zcold'] = 8.69+np.log10(data['Mzgas_disk']/(data['mcold_disk']*0.0134))
        data['zcold'] = 8.69+np.log10(data['zgas_disk']/0.0134)
        print '--> DONE\n'

    try:
        if np.all(data['cgf']==-99.0):
            print '\n+++++++++++++++++++\ncalculate gas-phase metallicity',
            data['cgf']=data['mcold']/data['mstar']
            print '--> DONE!\n'
    except:
        print 'SAGE:',
        data['cgf'] = (data['mcold_disk']/data['mstar'])
        print '--> DONE'

    try:
        mags='mAB_dA_total_'
        bands=[mags+'u', mags+'g', mags+'r', mags+'i', mags+'z']
        
        if np.all(data['kcorr_u']==-99.0):

            #data=data[np.where(data['mstar']>1e11)[0][:]]
            #print 'MSTAR CUT ', data.shape            
            
            print '\n+++++++++++++++++++\ncalculate k-correction with approximate formula!',
            k_corr_array = kcorrect_approx(data,
                                      redshift,
                                      bands)


            print np.info(k_corr_array)
            for k, band in enumerate(bands):
                print 'band:', band, 'num:', k
                data['kcorr_'+band[-1::]] = k_corr_array[band]
                print 'kcorrection: median =', "{0:.2f}".format(np.median(data['kcorr_'+band[-1::]])), '16th/84th', "{0:.2f}".format(np.percentile(data['kcorr_'+band[-1::]], 16)), '/', "{0:.2f}".format(np.percentile(data['kcorr_'+band[-1::]], 84))

            print '--> DONE!\n'
    except:
      pass


    return data

def count_nr_sats(data,
                  name_ID_central='fofID'):
    
    #sort so that the centrals are always first
    data.sort(order=[name_ID_central,'orphan'], axis=0)
    #print '\n', data[['haloid','hostid','orphan','nsats']][0:50]
    
    non_unique, index, count = np.unique(data[name_ID_central], return_index=True, return_counts=True)
    
    data['nsats'][index]=count-1
    
    print 'name of the ID to identfy centrals:', name_ID_central, '\n', data[[name_ID_central,'hostid','orphan','nsats']][0:50]
    
    
    return data

def assign_sats(data,
                name_ID_central='fofID'):
    
    #sort so that the centrals are always first
    data.sort(order=[name_ID_central,'hostid'], axis=0)

    
    #Set all types to default being satellites
    data['orphan']=1
    print '\n', data[[name_ID_central,'hostid','orphan','nsats']][0:300]
    non_unique, index, count = np.unique(data[name_ID_central], return_index=True, return_counts=True)
    
    #Set the central now to being centrals
    data['orphan'][index]=0
    
    data['nsats'][index]=count-1
    
    print 'name of the ID to identfy centrals:', name_ID_central, '\n', data[[name_ID_central,'hostid','orphan','nsats']][0:50]
    exit()
    
    return data

def sample_selection(data,
                 selection_key=None,
                 sample_key=None,
                 prop='sfr',
                 filename_indices='',
                 redshift=False,
                 ext_redshift=False,
                 parent_data=None,
                 color_index='color'): 


    print 'SFH SAMPLE SELECTION:', 'selection_key:', selection_key, 'sample_key:', sample_key
    
    if sample_key.startswith('r0'):
        #print self.parent_data
        print '300-region selection: ', sample_key,
        try:
            haloid_this_region=parent_data['haloid'][np.int(sample_key[1::])-1]
            data=data[np.where(data['hostid']==haloid_this_region)[:][0]]
            #print sample_key[1::], data['haloid'], haloid_this_region
        except:
            print 'No parent_data halos found!'
            return False
    
    elif sample_key=='passive2' or sample_key=='active2':
        Hubble_z=cd.hubble_z(redshift, **fidcosmo)
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
        print 'select galaxies by color!\n index:', color_index
        
        frac_to_select=int(np.floor(data.size/100.0*20.0)) 
        data = data[np.argsort(data[color_index])]
      
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

    elif sample_key=='mstarst11':

        print sample_key, 'select --> log10 (Mstar [Msun]) < 11.3',
        data = data[np.where(data['mstar']<10**11.3)[0][:]]
        print 'size:', data.size
            
    elif sample_key=='mstarge11':

        print sample_key, 'select --> log10 (Mstar [Msun]) >= 11.3',
        data = data[np.where(data['mstar']>=10**11.3)[0][:]]
        print 'size:', data.size          
                      
    elif sample_key=='mhalost13':
        
        print sample_key, 'select --> >= 12.3 log10 (Mvir [Msun]) < 13.6',
        data = data[np.where(data['mhalo']>=10**12.3)[0][:]]
        data = data[np.where(data['mhalo']<10**13.6)[0][:]]
        print 'size:', data.size
            
    elif sample_key=='mhaloge13':
        
        print sample_key, 'select --> log10 (Mvir [Msun]) >= 13.6',
        data = data[np.where(data['mhalo']>=10**13.6)[0][:]]
        print 'size:', data.size
    
    elif sample_key=='high-zcold':
        frac_to_select=int(np.floor(data.size/100.0*20.0)) 
        data = data[np.argsort(data['zcold'])]               
    
        print '--> select 20% highest zcold metallicities >>'
        start_size=data.size-frac_to_select                       
        data=data[start_size:data.size]
    
    elif sample_key=='high-zcold-red':
        frac_to_select=int(np.floor(data.size/100.0*20.0)) 
        data = data[np.argsort(data['zcold'])]               
    
        print '--> select 20% highest zcold metallicities >> '
        start_size=data.size-frac_to_select                       
        data=data[start_size:data.size]
    
        data=data[np.where(data[color_index]>=2.35)[0][:]]
        
    elif sample_key=='high-zcold-red20':
        
        print '1st step --> colour >= 2.35',
        data=data[np.where(data[color_index]>=2.35)[0][:]]
        
        frac_to_select=int(np.floor(data.size/100.0*20.0)) 
        data = data[np.argsort(data['zcold'])]               
    
        print '2nd --> select 20% highest zcold metallicities >> '
        start_size=data.size-frac_to_select                       
        data=data[start_size:data.size]   

    elif sample_key=='low-zcold':
        frac_to_select=int(np.floor(data.size/100.0*20.0)) 
        data = data[np.argsort(data['zcold'])]
        
        print '--> select 20% lowest zcold metallicities <<'
        data=data[0:int(frac_to_select)]
    
    elif sample_key=='low-zcold-red':
        frac_to_select=int(np.floor(data.size/100.0*20.0)) 
        data = data[np.argsort(data['zcold'])]
        
        print '--> select 20% lowest zcold metallicities <<'
        data=data[0:int(frac_to_select)]
        
        data=data[np.where(data[color_index]>=2.35)[0][:]]

    elif sample_key=='low-zcold-red20':
        
        print '1st step --> colour >= 2.35',
        data=data[np.where(data[color_index]>=2.35)[0][:]] 
        
        frac_to_select=int(np.floor(data.size/100.0*20.0)) 
        data = data[np.argsort(data['zcold'])]
        
        print '2nd step --> select 20% lowest zcold metallicities <<'
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
    
    elif sample_key=='redmstarst11':
       
        print  sample_key, '1st step --> red sample log10 (Mstar [Msun]) < 11.3',
        data = data[np.where(data['mstar']<10**11.3)[0][:]]                    
        print 'size:', data.size        

        frac_to_select=int(np.floor(data.size/100.0*20.0)) 
        data = data[np.argsort(data[color_index])]
    
        print '2nd step --> select 20% reddest (g-i)>>',
        start_size=data.size-frac_to_select                       
        data=data[start_size:data.size]
        print 'size:', data.size
        
    elif sample_key=='redmstarge11':
       
        print  sample_key, '1st step --> red sample log10 (Mstar [Msun]) >= 11.3',
        data = data[np.where(data['mstar']>=10**11.3)[0][:]]                    
        print 'size:', data.size        

        frac_to_select=int(np.floor(data.size/100.0*20.0)) 
        data = data[np.argsort(data[color_index])]
    
        print '2nd step --> select 20% reddest (g-i)>>',
        start_size=data.size-frac_to_select                       
        data=data[start_size:data.size]
        print 'size:', data.size        
        
    elif sample_key=='redmhalost13':
        
        print  sample_key, '1st step --> red sample >= 12.3  log10 (Mvir [Msun]) < 13.6',
        data = data[np.where(data['mhalo']>=10**12.3)[0][:]]  
        data = data[np.where(data['mhalo']<10**13.6)[0][:]]                    
        print 'size:', data.size      
        
        frac_to_select=int(np.floor(data.size/100.0*20.0)) 
        data = data[np.argsort(data[color_index])]
    
        print '2nd step --> select 20% reddest (g-i)>>',
        start_size=data.size-frac_to_select                       
        data=data[start_size:data.size]
        print 'size:', data.size       
 
    elif sample_key=='redmhaloge13':
        
        print sample_key, '1st step--> red sample  log10 (Mvir [Msun]) >= 13.6',
        data = data[np.where(data['mhalo']>=10**13.6)[0][:]]                    
        print 'size:', data.size      
        
        frac_to_select=int(np.floor(data.size/100.0*20.0)) 
        data = data[np.argsort(data[color_index])]
    
        print '2nd step --> select 20% reddest (g-i)>>',
        start_size=data.size-frac_to_select                       
        data=data[start_size:data.size]
        print 'size:', data.size                     
    
    elif sample_key=='lowmstarst11':

        print  sample_key, '1st step --> low sfr sample  log10 (Mstar [Msun]) < 11.3',
        data = data[np.where(data['mstar']<10**11.3)[0][:]]                    
        print 'size:', data.size           

        frac_to_select=int(np.floor(data.size/100.0*20.0)) 
        data = data[np.argsort(data['sfr'])]
    
        print '2nd --> select 20% lowest sfr',
        start_size=data.size-frac_to_select                       
        data=data[0:int(frac_to_select)]
        print 'size:', data.size        
             

    elif sample_key=='lowmhalost13':
        test = str(format(redshift, '0.2f'))
        print sample_key, '1st step',
        
        print  ' --> low sample >= 12.3  log10 (Mvir [Msun]) < 13.6',
        data = data[np.where(data['mhalo']>=10**12.3)[0][:]]  
        data = data[np.where(data['mhalo']<10**13.6)[0][:]]                    
        print 'size:', data.size

        frac_to_select=int(np.floor(data.size/100.0*20.0)) 
        data = data[np.argsort(data['sfr'])]

        print '2nd step--> select 20% lowest sfr',
        start_size=data.size-frac_to_select                       
        data=data[0:int(frac_to_select)]
        print 'size:', data.size    

         
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
    
        print 'filename_indices:', filename_indices
        myDataIndi = aD.ArangeData()
        indices = myDataIndi.readAnyFormat(config=False, mypath=filename_indices, data_format='ASCII', data_shape='shaped', delim='\t', mydtype=np.uint64, skiprow=2) 
        #indices[:,0]-> parentIndex one snapshot before
        #indices[:,1]-> nodeIndex one snapshot before
        
        #expand dimension of the numpy-array if it has only 1 row!
        if np.sum(indices.size)==2:
            indices = np.expand_dims(indices, axis=0)
        elif np.sum(indices.size)==0:
            #If there are not more progenitors found at this redshift, creat a dummy array to not interrupt the calculations
            indices=np.zeros((1,2), np.int8)                    
    
       # print 'nr indices:', np.sum(indices.size)
        data = data[np.where(data['orphan']==0)[:][0]] 
        #general selection
        #-------------------------
        #we intersect the current parentIndices with the nodeIndices a snapshot before, a unique array of parentIndices
        #which are present in both arrays are returned. Index1 gives the first occurance of the common IDs in parentIndex,
        #index2 of the same in the nodeIndex a snapshot before
        parentIndices, index1, index2 = np.intersect1d(data['parentIndex'], indices[:,1], return_indices=True)             
    
        #we test wheater the nodeIndex of this snapshot, in particular, those we already found in the first step 
        #(stored in index1) are linked to the a parentIndex. Thereby we check if those nodeIndices can be found in the
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
        
        if sample_key=='main':
            print 'SELECT MAIN PROGENITOR!'                               
            #the most massive is always the one with the lowest ID!
            
            data.sort(order=['parentIndex','satelliteNodeIndex'], axis=0)
       
          
        else:
            print 'SELECT MOST MASSIVE!' 
     
            data[::-1].sort(order=['parentIndex','mhalo'], axis=0)
       
        unique_parent_id, index1, count_parentIndex = np.unique(data['parentIndex'], return_index=True, return_counts=True)
        data=data[index1][::-1]
        
        myOutputIndi = oD.OutputData(config=False)
        myOutputIndi.writeIntoFile(filename_indices,
                   data[['parentIndex','nodeIndex']],
                   myheader= 'z='+str(format(redshift, '0.2f'))+'\nparentIndex(1) nodeIndex(2)',
                   data_format="%i\t%i",
                   mydelimiter='\t')
    
    elif sample_key=='valid':
        print 'check if hostids are valid!'
        data = data[np.where(data['hostid']>-2)[:][0]]
       
     
    else:
        print '--> NO SAMPLE SELECTION DONE!'
        
    return data

def assign_SFH_samples(data,
                       col_sample_key,
                       color_index='mAB_dA_total_cut_g_i'):
      
    import copy
    
    data[col_sample_key] = 'all:'
    data2analyse = copy.deepcopy(data)

    indices2track = data2analyse['hostid']
    
    list_samples        = ['low', 'high', 'passive', 'active', 'red', 'blue', 'low-zcold', 'high-zcold', 'low-zcold-red', 'high-zcold-red', 'mstarst11','mstarge11','mhalost13','mhaloge13']
    list_sample_keys    = ['l', 'h', 'p', 'a', 'r', 'b', 'lz','hz', 'lzr', 'hzr', 'lsm', 'hsm', 'lhm', 'hhm']

    import collections
    from collections import OrderedDict as od
        
    data2print={'Gal-dens': data2analyse[np.where(data2analyse['orphan']==0)[:][0]]}
    odict = od()
    
    odict[0]='Gal-dens'
    
    for i, sample, key in zip(np.arange(len(list_samples)),list_samples,list_sample_keys):
         
        print 'i', i, 'sample:', sample, 'key:', key,
        data_sample = sample_selection(data2analyse,
                                       sample_key=sample,
                                       color_index=color_index)
        

        #check assignment of samples
        #for item in indices2track:
            #print 'data:', data[['hostid', 'orphan', col_sample_key]][np.where(data['hostid']==item)[:][0]],
            #print 'data_sample:', data_sample[['hostid', 'orphan', col_sample_key]][np.where(data_sample['haloid']==item)[:][0]]
                   
        data_sample[col_sample_key]=key+':'
                    
        #print np.info(data_sample)
        data2print.update({sample: data_sample[np.where(data_sample['orphan']==0)[:][0]]})

        test = np.in1d(data_sample['hostid'],data['hostid'])

        print '-->', test[np.where(test==True)].size, 'found to crossmatch!\n\n'

        parentIndices, index1, index2 = np.intersect1d(data['hostid'], data_sample['hostid'], return_indices=True)

        data[col_sample_key][index1] = np.char.add(data[col_sample_key][index1],data_sample[col_sample_key][index2])
        
        #put flag
        data['flag_CMASS_'+key][index1] = 1
        odict[i+1]=sample
        
    #print odict
       
    return data, data2print, odict

def get_list_all_properties(data):
    
    myprops=[]
    for prop in data.dtype.names:
        prop_dtype = data[prop].dtype
        if prop_dtype!='int64' and prop_dtype!='int8' and prop.find('pos')==-1 and prop.startswith('zgasSF')==False and prop.startswith('env')==False and prop.startswith('flag')==False:
            #print prop, prop_dtype
            myprops.append(prop)
            
    return myprops

def calc_full_set_stats_boxplot(data,
                                log10=False):
    
    """functions calculates the tipical uncertainty estimations used to be presented in boxplots"""
    
    # try: 
    #     'min/max:', min(data), '/', max(data), np.nanmedian(data), data.size
    # except:
    #     print 'NO DATA -->', data.size
    
    if log10==True:
        data = data[np.where(data>0.0)[0][:]]
        data = np.log10(data)
        data = data[np.where(np.isfinite(data)==True)[0][:]]
        #print np.isfinite(data)
    else:
        data = data[np.where(np.isfinite(data)==True)[0][:]]
        
    # try:                
    #     print 'min/max:', min(data), '/', max(data), np.nanmedian(data), data.size 
    # except:
    #     print 'NO DATA -->', data.size 
        
    median  = np.nanmedian(data)
    Q1   = np.nanpercentile(data, 25)
    Q3   = np.nanpercentile(data, 75)
    IQR     = abs(Q3 - Q1)
    lower = Q1 - 1.5 * IQR
    upper = Q3 + 1.5 * IQR
    
    return median, Q1, Q3, IQR, lower, upper, data.size

def handle_HF_EAGLE_data(data,
                         key=None):

    def CrossMatch(data):
        
        redshift=0.0

        #data = crossmatch_catalogs(data, key='Yetli', redshift=redshift)
        #data = crossmatch_catalogs(data, key='Yetli_ev', redshift=redshift)
        #data = crossmatch_catalogs(data, key='Patricia3', redshift=redshift)
        #data = crossmatch_catalogs(data, key='Silvio', redshift=redshift)
        #data = crossmatch_catalogs(data, key='LSBGs', redshift=redshift)
        data = crossmatch_catalogs(data, key='t50', redshift=redshift)
        #data = crossmatch_catalogs(data, key='EAGLE_MKC', redshift=redshift)
        #data = crossmatch_catalogs(data, key='EAGLE_SFgas', redshift=redshift)       

        #data=data[np.where(data['mPart_stars']>=1e9)[:][0]]
        #data=data[np.where(data['hostid']==0)[:][0]]
        #data=data[np.where(data['redshift']==0)[:][0]]
       
        #data = crossmatch_catalogs(data, key='EAGLE_ev', redshift=redshift)
        
        return data


    def CalcProfile(data,
                    prop,
                    log_prop,
                    aperture_sizes=['1','10','100']):
        
        """Returns statistics of a property as a function of different aperture sizes
        
            input: a property measured in various apertures (property name must be the same for all apertures in the same structured array)
            output: histogram first column ascending the aperture sizes and second column the median value of
                    the properties (uses the same histo calculation as for SFHs!)
            
        """
        #print '--> calcualte aperture profile!'
        histo = np.zeros((len(aperture_sizes), 8), dtype=np.float32)
        for i, AP in enumerate(aperture_sizes):
            
            myprop=prop+'_AP'+AP    
            
            #print '\t\tz:', z, 'size:', data_snap.size,
            #print 'i:', i, 'AP:', AP, 'myprop:', myprop
                                          
            histo[i,0]=int(AP)
            #median, Q1, Q3, IQR, lower, upper, N[count]
            histo[i,1], histo[i,2], histo[i,3], histo[i,4], histo[i,5], histo[i,6], histo[i,7] = calc_full_set_stats_boxplot(data[myprop], log10=log_prop)
                
        return histo                   
                  
    
    def CalcSFH(data,
                flag_name='',
                sample_names=[''],                
                sample_flags=[-99],
                output_key='',
                props='',
                calc_profile=False,
                aperture_sizes=['1','10','100']):
        
        import copy
        
        if output_key!='':
            output_key_space='_'
        else:
            output_key_space=''
            
        z_array = np.unique(data['redshift'])
        z_array.sort()

        mylist=[]
        for item in z_array: 
            mylist.append(format(item, '0.2f'))
        
        #print mylist
        #exit()

        if props=='':
            props = ['r200c', 'mbh','mstar_30kpc', 'mhalo_200c', 'mhalo_30kpc','mPart_BH', 'mgas_SF', 'SHMR', 'bheff', 'cgf', 'Tcons', 'bh_acc_rate',\
                  'rVmax', 'angM_norm_stars', 'angM_norm_SFgas', 'angM_norm_NSFgas', 'sfr_30kpc', 'ssfr_30kpc', 'SHMR_30kpc',\
                  'vmax', 'vdisp_30kpc', 'mgas_30kpc', 'age_mean_stars','rhalf_stars','mstar_birth','z_mean_birth_stars',\
                  'Sigma_stars', 'Sigma_sfr','vpec', 'vpec_norm']
          

        for item, item_flag in zip(sample_names,sample_flags):
            
            mydata=copy.deepcopy(data)
            
            print 'flag_name:', flag_name, 'item:', item, 'flag:', item_flag
            if item!='':
                sample_space='_'
                name_space='_'
            else:
                sample_space=''
                name_space=''
                
            if item_flag!=-99:
                data_item = mydata[np.where(mydata[flag_name]==item_flag)[0][:]]
            else:
                data_item = mydata

            for prop in props:
                if prop.startswith('r') or prop.startswith('MAB') or prop.startswith('n') or prop.startswith('12') or prop.startswith('Fe') or prop.startswith('vm') or prop.find('-')!=-1 or prop.startswith('SB') or prop.startswith('vdis') or prop=='vpec' or prop=='age_mean_stars' or prop.startswith('z_mean'):
                    log_prop=False
                else:
                    log_prop=True
                    
                for key, flag in zip(['', 'LSB', 'HSB'],[-99,1,0]):
                    print 'key:', key, 'flag:', flag,
                    
                    if key!='':
                        key_space='_'
                    else:
                        key_space=''
                    
                    mydata_key=copy.deepcopy(data_item)
                    #print 'min/max:', min(mydata_key[prop]), '/', max(mydata_key[prop])
                    
                    header='EAGLE 100Mpc HF-II '+key+' sample: '+str(item)+' flag_name: '+flag_name+' flag: '+str(item_flag)+'\n'
                    
                    if calc_profile==True:
                        
                        for i,z in enumerate(z_array):
                            if flag!=-99:
                                data_profile = mydata_key[(np.where((mydata_key['redshift']==z) & (mydata_key['flag_SB_B']==flag)))[0][:]]
                            else:
                                data_profile = mydata_key[np.where(mydata_key['redshift']==z)[0][:]]
                                                            
                            array_profile = CalcProfile(data_profile, prop, log_prop, aperture_sizes=aperture_sizes)
                            
                            if i==0: 
                                z=0.0
                            #print array_profile
                            myOutput = oD.OutputData(config=False)
                            myOutput.writeIntoFile(mycomp+'Documents/Drafts/Hidden-Figures/Paper-II/profile/EAGLE_100Mpc_profile_z_'+str(format(z, '0.2f'))+'_'+output_key+output_key_space+prop+key_space+key+sample_space+flag_name+name_space+item+'.txt',
                                        array_profile,
                                        myheader= header+'aperture[pkpc]\tmedian\tQ1\tQ3\tIQR\tlower\tupper\tN[count]',
                                        data_format="%0.5e",
                                        mydelimiter='\t') 
                        
                 
                    else:
                        array_SFH = np.zeros((len(z_array), 8), dtype=np.float32)
                        for i,z in enumerate(z_array):                     
    
                            if flag!=-99:
                                data_snap = mydata_key[prop][(np.where((mydata_key['redshift']==z) & (mydata_key['flag_SB_B']==flag)))[0][:]]
                            else:
                                data_snap = mydata_key[prop][np.where(mydata_key['redshift']==z)[0][:]]
    
                                                            
                            #print '\t\tz:', z, 'size:', data_snap.size,
                            #print 'min/max:', min(data_snap), '/', max(data_snap), np.nanmedian(data_snap)
                            #data_snap = data_snap[np.where(data_snap>0.0)[0][:]]
                            #print 'min/max:', min(data_snap), '/', max(data_snap), np.nanmedian(data_snap)                                            
                            array_SFH[i,0]=z
                            #median, Q1, Q3, IQR, lower, upper
                            array_SFH[i,1], array_SFH[i,2], array_SFH[i,3], array_SFH[i,4], array_SFH[i,5], array_SFH[i,6], array_SFH[i,7] = calc_full_set_stats_boxplot(data_snap, log10=log_prop)
                         
                                
                                               
                        array_SFH[0,0]= 0.05
                        
                        myOutput = oD.OutputData(config=False)
                        myOutput.writeIntoFile(mycomp+'Documents/Drafts/Hidden-Figures/Paper-II/SFH/EAGLE_100Mpc_SFH_'+output_key+output_key_space+prop+key_space+key+sample_space+flag_name+name_space+item+'.txt',
                                    array_SFH,
                                    myheader= header+'z\tmedian\tQ1\tQ3\tIQR\tlower\tupper\tN[count]',
                                    data_format="%0.5e",
                                    mydelimiter='\t')
                       
                    
        #exit()
        
    def CalcAssemblyHistory(data,
                            key=None):
        
        def load_full_centrals():
            mypath=mycomp+'anaconda/pro/data/EAGLE_ev_100Mpc/EAGLE_ev_100Mpc_full_centrals_main_progenitor_trees.txt'
            print mypath
            data2cross = pd.read_csv(mypath, skiprows=1,\
                              names=['fofID',  'progFofID',  'redshift',  'haloid',  'hostid',  'galaxyID',  'lastProgID',  'topLeafID',\
                                     'mhalo_200c',  'r200c',  'nSubhalos', 
                                     'xpos',  'ypos',  'zpos',  'x_pos_subhalo',  'y_pos_subhalo',  'z_pos_subhalo',\
                                     'xpos_cof',  'ypos_cof',  'zpos_cof',
                                     'mPart_DM',  'mPart_BH',  'mPart_gas',  'mPart_stars',\
                                     'sfr',  'vmax',  'rVmax',  'vdisp',  'mstar_birth',  'mbh',  'bh_acc_rate',\
                                     'mstar_30kpc',  'sfr_30kpc',  'vdisp_30kpc',  'mgas_30kpc',  'mhalo_30kpc',\
                                     'rhalf_stars',  'rhalf_gas',  'mgas_SF',  'mgas_NSF',\
                                     'spinNSFGas_x',  'spinNSFGas_y',  'spinNSFGas_z',\
                                     'spinSFGas_x',  'spinSFGas_y',  'spinSFGas_z',  'spinStars_x',  'spinStars_y',  'spinStars_z',\
                                     'z_mean_birth_stars',  'age_mean_stars',\
                                     'x_vel',  'y_vel',  'z_vel'],\
                                     sep=',')
                
            data2cross = df_to_sarray(data2cross)
        
            data2cross = rcfuncs.append_fields([data2cross],
                                               ['envr','flag_t50_st2ha', 'flag_t50_bh2ha', 'flag_t50_bh2st', 'flag_tform_halo','vpec_norm', 'vpec','Sigma_sfr','Sigma_stars', 'flag_dens','Tcons', 'cgf', 'bheff', 'SHMR','flag_SB_B','flag_tform_stars', 'flag_tform_bh', 'angM_norm_NSFgas','angM_norm_SFgas','angM_norm_stars', 'ssfr_30kpc', 'SHMR_30kpc'],
                                               [data2cross['nSubhalos'],data2cross['nSubhalos'],data2cross['nSubhalos'],data2cross['nSubhalos'],data2cross['nSubhalos'],data2cross['sfr_30kpc'],data2cross['sfr_30kpc'],data2cross['sfr_30kpc'],data2cross['sfr_30kpc'],data2cross['nSubhalos'],data2cross['sfr_30kpc'],data2cross['sfr_30kpc'],data2cross['sfr_30kpc'],data2cross['mhalo_200c'],data2cross['nSubhalos'],data2cross['nSubhalos'],data2cross['nSubhalos'],data2cross['spinSFGas_x'],data2cross['spinSFGas_x'],data2cross['spinStars_x'],data2cross['sfr_30kpc'], data2cross['mhalo_30kpc']], 
                                               usemask=False)
                             
            #calculate angular momenta
            data2cross['angM_norm_NSFgas'] = calc_norm_vector(data2cross[['spinNSFGas_x','spinNSFGas_y','spinNSFGas_z']])
            data2cross['angM_norm_SFgas'] = calc_norm_vector(data2cross[['spinSFGas_x','spinSFGas_y','spinSFGas_z']])
            data2cross['angM_norm_stars'] = calc_norm_vector(data2cross[['spinStars_x','spinStars_y','spinStars_z']])

            data2cross['vpec'] = calc_norm_vector(data2cross[['x_vel','y_vel','z_vel']])
            data2cross['vpec_norm'] = calc_norm_vector(data2cross[['x_vel','y_vel','z_vel']])/data2cross['mstar_30kpc']
            
            data2cross['ssfr_30kpc']=data2cross['sfr_30kpc']/data2cross['mstar_30kpc']
            data2cross['SHMR_30kpc']=data2cross['mstar_30kpc']/data2cross['mhalo_30kpc']
            data2cross['SHMR']      =data2cross['mstar_30kpc']/data2cross['mhalo_200c']
            data2cross['cgf']       =data2cross['mgas_SF']/data2cross['mstar_30kpc']
            data2cross['bheff']     =data2cross['mbh']/data2cross['mhalo_200c']
            data2cross['Tcons']     =data2cross['mgas_SF']/data2cross['sfr_30kpc']/1e9
            data2cross['Sigma_stars']=(data2cross['mstar_30kpc']/2.0)/(2*np.pi*(data2cross['rhalf_stars'])**2)
            data2cross['Sigma_sfr']=data2cross['sfr']/(2*np.pi*(data2cross['rhalf_stars'])**2)
            data2cross['z_mean_birth_stars'] = cd.lookback_time(data2cross['z_mean_birth_stars'], z0=0.0, **fidcosmo)/3.1536e+16
            
            return data2cross
            
        def load_full_centrals_metals():
            mypath=mycomp+'anaconda/pro/data/EAGLE_ev_100Mpc/EAGLE_ev_100Mpc_full_centrals_main_progenitor_trees_dens3_consistence_test.txt'
            print mypath
            data2cross = pd.read_csv(mypath, skiprows=0,\
                              names=['fofID', 'propFofID', 'redshift', 'haloid', 'hostid', 'GalaxyID', 'lastProgID', 'topLeafID',
                                     'mhalo_200c', 'mhalo_500c', 'mhalo_2500c', 'r200c', 'r500c', 'r2500c', 'nsats',
                                     'sfr', 'rhalf_stars',
                                     'zgasSF_H', 'zgasSF_O', 'zgasSF_N', 'zgasSF_Fe',
                                     'zgasNSF_H', 'zgasNSF_O','zgasNSF_N', 'zgasNSF_Fe',
                                     'zcold', 'zhot', 'zstars',
                                     'Etot', 'Ekin', 'Etherm', 'temp_SF', 'temp_NSF'],\
                                     sep=' ')
                
            data2cross = df_to_sarray(data2cross)          
 
            data2cross = rcfuncs.append_fields([data2cross],
                                               ['envr','flag_SB_B','flag_dens', '12+log10OH', 'FeH','flag_tform_stars', 'flag_tform_bh', 'flag_tform_halo'],
                                               [data2cross['nsats'], data2cross['nsats'], data2cross['nsats'], data2cross['zcold'],data2cross['zcold'],\
                                                np.zeros(data2cross.size,),np.zeros(data2cross.size,),np.zeros(data2cross.size,)],
                                                dtypes=['i4','i4','i4','f8','f8','i4','i4','i4'],
                                               usemask=False)
            
            data2cross['12+log10OH'] = 12+np.log10(data2cross['zgasSF_O']/data2cross['zgasSF_H'])
            data2cross['FeH']      = np.log10(data2cross['zgasSF_Fe']/data2cross['zgasSF_H'])                  
            
            myprops=['mhalo_500c', 'mhalo_2500c', 'r500c', 'r2500c', 'nsats',
                    'FeH', '12+log10OH', 'zcold', 'zhot', 'zstars', 'Etot', 'Ekin', 'Etherm', 'temp_SF', 'temp_NSF']         
            
            return data2cross, myprops

        
        def load_full_centrals_mags():
           mypath=mycomp+'anaconda/pro/data/EAGLE_ev_100Mpc/EAGLE_ev_100Mpc_full_centrals_main_progenitor_trees_dens3_mags.txt'
           print mypath
           data2cross = pd.read_csv(mypath, skiprows=0,\
                             names=['fofID', 'propFofID', 'redshift', 'haloid', 'hostid', 'GalaxyID', 'lastProgID', 'topLeafID',
                                    'mhalo_200c','sfr', 'rhalf_stars',
                                    'MAB_dA_Johnson_B', 'MAB_dA_Johnson_V', 'MAB_dA_SDSS_g', 'MAB_dA_SDSS_r', 'MAB_dA_SDSS_i'],\
                                    sep=' ')
               
           data2cross = df_to_sarray(data2cross)
           
           data2cross = rcfuncs.append_fields([data2cross],
                                              ['envr', 'flag_SB_B','flag_dens', 'SB_mu_eff_r', 'r-i', 'g-i', 'g-r', 'B-V','flag_tform_stars', 'flag_tform_bh', 'flag_tform_halo'],
                                              [np.zeros(data2cross.size,),np.zeros(data2cross.size,),np.zeros(data2cross.size,),np.zeros(data2cross.size,),np.zeros(data2cross.size,),\
                                               np.zeros(data2cross.size,),np.zeros(data2cross.size,),np.zeros(data2cross.size,),\
                                               np.zeros(data2cross.size,),np.zeros(data2cross.size,),np.zeros(data2cross.size,)],\
                                              dtypes=['i4','i4','i4','f8','f8','f8','f8','f8','i4','i4','i4'],
                                              usemask=False)
       
           radius = conv_radius_pc_to_arcsec(data2cross['rhalf_stars']*1000.0)
           data2cross['SB_mu_eff_r'] = data2cross['MAB_dA_SDSS_r']    + 2.5*np.log10(2.0*np.pi*(radius)**2.0)
           
           #print data2cross
           
           data2cross['r-i']=data2cross['MAB_dA_SDSS_r']-data2cross['MAB_dA_SDSS_i']
           data2cross['g-i']=data2cross['MAB_dA_SDSS_g']-data2cross['MAB_dA_SDSS_i']
           data2cross['g-r']=data2cross['MAB_dA_SDSS_g']-data2cross['MAB_dA_SDSS_r']
           data2cross['B-V']=data2cross['MAB_dA_Johnson_V']-data2cross['MAB_dA_Johnson_B']
           
           myprops=['MAB_dA_Johnson_B', 'MAB_dA_Johnson_V', 'MAB_dA_SDSS_g', 'MAB_dA_SDSS_r', 'MAB_dA_SDSS_i',\
                    'SB_mu_eff_r', 'r-i', 'g-i', 'g-r', 'B-V']
               
           return data2cross, myprops

            
        def load_full_centrals_AP():
            
            import copy
            from itertools import chain
            
            def load_AP_data(AP_list):
                #print 'AP_list:', AP_list, AP_list[0], AP_list[-1]
                mypath=mycomp+'anaconda/pro/data/EAGLE_ev_100Mpc/EAGLE_ev_100Mpc_full_centrals_main_progenitor_trees_dens3_AP_'+AP_list[0]+'-'+AP_list[-1]+'.txt'
                print mypath              
                names=['fofID', 'propFofID', 'redshift', 'haloid', 'hostid', 'GalaxyID', 'lastProgID', 'topLeafID', 'mhalo_200c']
                for item in AP_list:
                    #print 'item:', item
                    names+=['vdisp_AP'+item, 'sfr_AP'+item,'mPart_bh_AP'+item,'mPart_DM_AP'+item,'mPart_gas_AP'+item,'mPart_stars_AP'+item]
                    
                #print 'names:', names
                mydata = pd.read_csv(mypath, skiprows=0,\
                                         names=names,\
                                         sep=' ')
                
                return df_to_sarray(mydata)


            
            unique_AP_list=[]

            for k, AP_list in enumerate([['1','3'],['5','10','20'],['30','40','50'],['70','100']]):
                print 'k:', k, 'AP_list', AP_list,
                
                mydata = load_AP_data(AP_list)
                
                unique_AP_list.extend(AP_list)
                
                print 'size:', mydata.size
                #print np.info(mydata) #list(mydata.dtype.names[9::])
                                 
                if k==0:
                    data2cross = copy.deepcopy(mydata)
                else:               
                    col_names = list(mydata.dtype.names[9::])
                    dtype_list=[]
                    for name in col_names:
                        dtype_list.append(np.array((mydata.size,), dtype=np.float32))
                    
                    #print 'col_names:', col_names, '\n dtype_list', dtype_list

                    data2cross = rcfuncs.append_fields([data2cross],
                                                        col_names,
                                                        dtype_list, 
                                                        usemask=False)                    
                    
                    parentIndices, index1, index2 = np.intersect1d(data2cross['GalaxyID'], mydata['GalaxyID'], return_indices=True)
                    
                    data2cross[list(mydata.dtype.names[9::])][index1]=mydata[list(mydata.dtype.names[9::])][index2]
                 
                    
            data2cross = rcfuncs.append_fields([data2cross],
                                               ['envr', 'flag_SB_B','flag_dens','flag_tform_stars', 'flag_tform_bh', 'flag_tform_halo'],
                                               [np.zeros(data2cross.size,),np.zeros(data2cross.size,),np.zeros(data2cross.size,),\
                                                np.zeros(data2cross.size,),np.zeros(data2cross.size,),np.zeros(data2cross.size,)], 
                                               usemask=False)
            
            myprops=[]    
            for ap in aperture_sizes:
                data2cross = rcfuncs.append_fields([data2cross],
                                                    ['ssfr_AP'+ap,'cgf_AP'+ap, 'bheff_AP'+ap, 'SHMR_AP'+ap, 'Tcons_AP'+ap],
                                                    [data2cross['sfr_AP'+ap],data2cross['sfr_AP'+ap],data2cross['sfr_AP'+ap],data2cross['sfr_AP'+ap],data2cross['sfr_AP'+ap]], 
                                                    usemask=False)
            
                data2cross['ssfr_AP'+ap]=data2cross['sfr_AP'+ap]/data2cross['mPart_stars_AP'+ap]
                data2cross['cgf_AP'+ap]=data2cross['mPart_gas_AP'+ap]/data2cross['mPart_stars_AP'+ap]
                data2cross['bheff_AP'+ap]=data2cross['mPart_bh_AP'+ap]/data2cross['mPart_DM_AP'+ap]
                data2cross['SHMR_AP'+ap]=data2cross['mPart_stars_AP'+ap]/data2cross['mPart_DM_AP'+ap]
                data2cross['Tcons_AP'+ap]=data2cross['mPart_gas_AP'+ap]/data2cross['sfr_AP'+ap]/1e9
            
                myprops+=['vdisp_AP'+ap, 'sfr_AP'+ap,'mPart_bh_AP'+ap,'mPart_DM_AP'+ap,'mPart_gas_AP'+ap,'mPart_stars_AP'+ap]
                
            myprops_profile=['vdisp', 'sfr', 'mPart_stars', 'mPart_bh', 'mPart_DM', 'mPart_gas','ssfr', 'bheff', 'SHMR', 'cgf', 'Tcons']
            #print np.info(data2cross)
                
            return data2cross, unique_AP_list, myprops, myprops_profile
        
        def find_and_assign_flag(data1, data2, flag_name, default_fill=-99):
            
            data2[flag_name]=default_fill
            
            parentIndices, index1, index2 = np.intersect1d(data1['fofID'], data2['fofID'], return_indices=True)
            for index, id1 in zip(parentIndices, index1):
                #print 'fofID:', index, 'id1:', id1,
                flag = data1[flag_name][id1]
                #print 'flag:', flag
                data2[flag_name][np.where(data2['fofID']==index)[:][0]] = flag
                
            return data2
        
        def calc_fractions(data):
            
            """calc procentages of number in each sub-sample and print to screen"""
            
            print 'ngals:', data.size,  format(100.0/data.size*data.size, '0.0f'), '%'
            for mass in ['halo', 'stars', 'bh']:
                
                print '\nmass:', mass, '\n##################\n' 
                for t50, t50_key in zip(['all', 'early', 'late'], [-99,0,1]):   
                    print 't50:', t50, t50_key,
                    if t50=='all':
                        basis = data
                    else:
                        basis = data[np.where(data['flag_tform_'+mass]==t50_key)[:][0]]
                    
                    print 'ngals:', basis.size,  format(100.0/data.size*basis.size, '0.0f'), '%' 

                    for envr, envr_key in zip(['all','S', 'W', 'V'], [-99,3,2,1]):
                        print '\tenvr:', envr,
                        if envr=='all':
                            basis_envr = basis
                        else:
                            basis_envr = basis[np.where(basis['envr']==envr_key)[:][0]]                    
                        
                        print 'ngals:', basis_envr.size,  format(100.0/basis.size*basis_envr.size, '0.0f'), '%'                    
                        
                        for SB, SB_key in zip(['all', 'HSB', 'LSB'], [-99,0,1]):
                            print '\t\tSB:', SB, SB_key,
                            if SB=='all':
                                basis_envr_SB = basis_envr
                            else:
                                basis_envr_SB = basis_envr[np.where(basis_envr['flag_SB_B']==SB_key)[:][0]]
                                                           
                            print 'ngals:', basis_envr_SB.size,  format(100.0/basis_envr.size*basis_envr_SB.size, '0.0f'), '%'  
                            
            exit()
        
        def func(x, a, b, c):

            return a * np.exp(-b * x) + c
        
        import pandas as pd        
        import numpy.lib.recfunctions as rcfuncs        
        from cosmolopy import cparam, cd, cc
        fidcosmo = cparam.PlanckEAGLE(flat=True, extras=False)
        #print fidcosmo

        #data2cross = load_full_centrals()
        #data2cross, myprops = load_full_centrals_metals()
        data2cross, myprops = load_full_centrals_mags()
        #data2cross, aperture_sizes,myprops, myprops_profile = load_full_centrals_AP()


        print '\nFOLLOW MERGER TREES AND CALCULATE HALF-MASS ASSEMLBY TIMES & CROSSMATCH WITH CATALOG AT z=0!\n-----------------\nparent catalog:',\
        data.size, 'data2cross:', data2cross.size
        #print np.info(data2cross)
                                             
        for flag in ['envr','flag_dens','flag_SB_B', 'flag_tform_stars', 'flag_tform_bh', 'flag_tform_halo']:
            print 'assign flag -->', flag, 
            data2cross = find_and_assign_flag(data, data2cross, flag)
            print '\t ASSIGNED!\t',
            print data2cross[np.where(data2cross[flag]==0)[:][0]].size, data2cross[np.where(data2cross[flag]==1)[:][0]].size, data2cross[np.where(data2cross[flag]-99)[:][0]].size

        data2cross = data2cross[np.where(data2cross['flag_dens']==1)[:][0]]
            
        
        #data2cross = data2cross[np.where(data2cross['envr']==1)[:][0]]
        #print data2cross.size, data2cross['envr']
        
        #data2cross = data2cross[np.where(data2cross['flag_tform_halo']==0)[:][0]]
        #print data2cross.size, data2cross['flag_tform_halo']
              
        # CalcSFH(data2cross,
        #         props=myprops,
        #         output_key='dens3_t50_halo_early',
        #         calc_profile=True,
        #         aperture_sizes=aperture_sizes)

        CalcSFH(data2cross,
                flag_name='flag_tform_bh',
                props=myprops,
                output_key='dens3',
                sample_names=['early','late'],                
                sample_flags=[0,1])
        
        # CalcSFH(data2cross,
        #         flag_name='flag_tform_bh',
        #         sample_names=['early','late'],                
        #         sample_flags=[0,1])
        
        # CalcSFH(data2cross,
        #         flag_name='flag_tform_halo',
        #         output_key='dens3',
        #         props=myprops,
        #         sample_names=['early', 'late'],
        #         sample_flags=[0, 1])

        # CalcSFH(data2cross,
        #         flag_name='envr',
        #         sample_names=['S','W', 'V'],                
        #         sample_flags=[3,2,1])
        
        exit()
        
        CalcSFH(data2cross,
                flag_name='flag_t50_st2ha',
                sample_names=['mstar1st','mhalo1st'],                
                sample_flags=[0,1])
                
        CalcSFH(data2cross,
                flag_name='flag_t50_bh2st',
                sample_names=['mbh1st','mstar1st'],                
                sample_flags=[0,1])
        
        CalcSFH(data2cross,
                flag_name='flag_t50_bh2ha',
                sample_names=['mbh1st','mhalo1st'],                
                sample_flags=[0,1])
         
        exit()
        #Follow trees!
        from scipy import interpolate
        from scipy.optimize import curve_fit
        from scipy.interpolate import CubicSpline
        array_short_trees_fofIDz0 = []
        array_error_inter         = []
        array_error_curve_fit     = []
        array_error_cs            = []
        array_excluded_fofIDz0    = []
        array_variation_excluded  = []
        array_unrealistic_fit     = []
        array_low_delta_tform     = []
        
        array_fit_prefered    = []
        array_inter_prefered    = []
        
        test_fofIDs = [28000000000000,28000000000002,28000000000010,28000000000099,28000000000249,\
                       28000000001726,28000000002846,28000000004498,28000000005171,28000000004698]
            
        threshold_data_fit_z0 = 2.0
        print 'RUNNING THROUGH TREES!\n---------------------------\n'
        for fofID in data['fofID']:
        #for fofID in test_fofIDs:
            #print 'fofID:', fofID, 
            topLeafIDz0=data['topLeafID'][np.where(data['fofID']==fofID)[:][0]]                      
            
            #print 'topLeafIDz0:', topLeafIDz0
            
            #Check is the merger tree is sufficiently long to get reliable results!                        
            if len(data2cross[np.where(data2cross['fofID']==fofID)[:][0]])>=18:

                #Detect tree
                data_ev_fofID=data2cross[np.in1d(data2cross['topLeafID'], topLeafIDz0)]
                #print data_ev_fofID[['redshift','haloid','fofID', 'flag_SB_B']]
                
                data_ev_fofID.sort(order=['redshift'], axis=0)
                data_ev_fofID['redshift'][0]=0.0
                #print data_ev_fofID[['redshift','haloid','topLeafID', 'galaxyID']]
                
                bucket = np.zeros((30,13), np.float32)
                #print data_ev_fofID[['redshift','fofID', 'mhalo_200c', 'mstar_30kpc', 'mbh']]
                header=''
                
                key_calc='other'
                if key_calc=='masses':
                    for prop, key in zip(['mhalo_200c', 'mstar_30kpc', 'mbh'], ['halo', 'stars', 'bh']):
                        
                        #print '\t half-assembly redshifts --> prop:', prop, 'key:', key
                        
                        check_fit   = True
                        check_cs    = True
                        use_inter   = False
                        data['nr_strikes_'+key][np.where(data['fofID']==fofID)[:][0]] = 0
                        strikes     = 0
                        f_fit       = False
                        f_cs        = False
                        f_inter_50  = False
                        f_inter_70  = False
                        
                        t50_std     = -99
                        t70_std     = -99
                        
    
                        if key=='bh':
                            header+=prop+'\t\t\t--> '
                        else:
                            header+=prop+'\t--> '
    
                        if key=='halo':
                            col_fit = 1
                            col_int = 2
                            col_cs = 3
                            col_data = 4
                        elif key=='stars':
                            col_fit = 5
                            col_int = 6
                            col_cs = 7
                            col_data = 8
                        else:
                            col_fit = 9
                            col_int = 10
                            col_cs = 11
                            col_data = 12                  
                       
                        data2process = data_ev_fofID[np.where(data_ev_fofID[prop]>0.0)[:][0]]                               
                        
                        data2process.sort(order=['redshift'], axis=0)                   
                        try:       
                            data2process[prop] = data2process[prop]/data2process[prop][0]
                            
                            bucket[0:data2process[prop].size,  0]       = data2process['redshift']
                            bucket[0:data2process[prop].size, col_data] = data2process[prop]
                                              
                            try:                        
                                x_new = np.logspace(np.log10(0.005),np.log10(20),num=500)
                                data_inter = np.interp(x_new, data2process['redshift'], data2process[prop])
            
                                f_inter_50 = x_new[np.where(data_inter==data_inter[np.abs(data_inter-0.5).argmin()])[:][0]][0]
                                f_inter_70 = x_new[np.where(data_inter==data_inter[np.abs(data_inter-0.7).argmin()])[:][0]][0]
                                                        
                                data['t50_'+key][np.where(data['fofID']==fofID)[:][0]] = cd.lookback_time(f_inter_50, z0=0.0, **fidcosmo)/3.1536e+16 #seconds to Gyrs                           
                                data['t70_'+key][np.where(data['fofID']==fofID)[:][0]] = cd.lookback_time(f_inter_70, z0=0.0, **fidcosmo)/3.1536e+16 #seconds to Gyrs
                                                           
                                bucket[0:data2process[prop].size, col_int]  = np.interp(data2process['redshift'], data2process['redshift'], data2process[prop]) 
                                
                                header+= '\tz 50%/70%\t\t inter: '+str(format(f_inter_50, '0.3f'))+'/'+format(f_inter_70, '0.3f')
                            except:
                                header+= '\tz 50%/70%\t\t inter: --/--'
                                array_error_inter.append(fofID)
        
            
                            try:
                                #Curvefit
                                popt, _ = curve_fit(func, data2process['redshift'], data2process[prop])                    
                                # summarize the parameter values
                                a, b, c = popt
                                #print 'y = %.5f * x + %.5f * x^2 + %.5f' % (a, b, c)  
                                
                                f_fit = interpolate.interp1d(func(data2process['redshift'], a, b, c) , data2process['redshift'])
                                data_fit = func(data2process['redshift'], a, b, c)
                                bucket[0:data2process[prop].size, col_fit]  = data_fit
                                
                                header+='\tcurve fit: '+str(format(f_fit(0.5), '0.3f'))+'/'+format(f_fit(0.7), '0.3f')
                            except:
                                header+='\tcurve fit: --/--'
                                array_error_curve_fit.append(fofID)
                                check_fit=False
                             
                            try:                        
                                #CubisSpline (cs) fit                    
                                cs = CubicSpline(data2process['redshift'], data2process[prop])
                                
                                f_cs = interpolate.interp1d(cs(data2process['redshift']), data2process['redshift'])
                                bucket[0:data2process[prop].size, col_cs]   = cs(data2process['redshift'])
                                
                                header+='\tcs: '+format(f_cs(0.5), '0.3f')+'/'+format(f_cs(0.7), '0.3f')
                            except:
                                header+='\tcs: --/--'
                                array_error_cs.append(fofID)
                                check_cs=False
                                
                            header+='\t std: '    
                            try:   
                                data['t50_'+key+'_std'][np.where(data['fofID']==fofID)[:][0]] = np.std([f_inter_50,f_fit(0.5).flat[0],f_cs(0.5).flat[0]])                        
                                #t50_std = np.std([f_inter_50,f_fit(0.5).flat[0],f_cs(0.5).flat[0]])                        
    
                                header+=format(data['t50_'+key+'_std'][np.where(data['fofID']==fofID)[:][0]].flat[0], '0.4f')+'/'
                            except:
                                header+='--/'                               
                                    
                                
                            try:
                                data['t70_'+key+'_std'][np.where(data['fofID']==fofID)[:][0]] = np.std([f_inter_70,f_fit(0.7).flat[0],f_cs(0.7).flat[0]])
                                #t70_std = np.std([f_inter_70,f_fit(0.7).flat[0],f_cs(0.7).flat[0]])
                               
                                header+=format(data['t70_'+key+'_std'][np.where(data['fofID']==fofID)[:][0]].flat[0], '0.4f')     
                            except:
                                header+='--'
                                                         
                            header+='\t Delta(inter-cs): '
                            try:                           
                                header+=format(abs(f_inter_50-f_cs(0.5)), '0.3f')+'/'+format(abs(f_inter_70-f_cs(0.7)), '0.3f')                            
                            except:
                                header+='--'
                                
                            if check_fit==False and check_cs==True:
                                use_inter=True
                                array_inter_prefered.append(fofID)
                                
                            #Check tree consistency
                            if abs(data2process[prop][0]-data2process[prop][1])>0.5:
                                strikes+=1
                            if abs(data2process[prop][0]-data2process[prop][2])>0.5:
                                strikes+=1
                            if abs(data2process[prop][1]-data2process[prop][2])>0.5:
                                strikes+=1
        
                            if strikes>1: 
                                array_variation_excluded.append(fofID)
                                use_inter=False
                                                                                           
                            #print 'check_fit:', check_fit, 'check_cs:', check_cs, 'use_inter: ', use_inter, 'srikes:', strikes, 'data_fit_z0:', data_fit[0].flat[0],
                            
                            if use_inter==True:
                                header+=' use inter! '
                                data['flag_use_inter_'+key][np.where(data['fofID']==fofID)[:][0]] = 1   
                                
                            data['flag_std_fit_t50_'+key][np.where(data['fofID']==fofID)[:][0]] = -99
                            if (data['t50_'+key+'_std'][np.where(data['fofID']==fofID)[:][0]].flat[0]!=-99 and strikes<2) or use_inter==True:
                            #if (t50_std!=-99 and strikes<2) or use_inter==True:
                                try:
                                    if data['t50_'+key+'_std'][np.where(data['fofID']==fofID)[:][0]].flat[0]>0.05 and abs(f_inter_50-f_cs(0.5))>0.05 and data_fit[0]<threshold_data_fit_z0:
                                        data['t50_'+key][np.where(data['fofID']==fofID)[:][0]] = cd.lookback_time(f_fit(0.5), z0=0.0, **fidcosmo)/3.1536e+16 #seconds to Gyrs
                                        data['flag_std_fit_t50_'+key][np.where(data['fofID']==fofID)[:][0]] = 1
                                    else:
                                        
                                        data['flag_std_fit_t50_'+key][np.where(data['fofID']==fofID)[:][0]] = 0
                                except:
                                    pass
                            else:
                                data['t50_'+key][np.where(data['fofID']==fofID)[:][0]] = -99                              
        
        
                            data['flag_std_fit_t70_'+key][np.where(data['fofID']==fofID)[:][0]] = -99
                            if (data['t70_'+key+'_std'][np.where(data['fofID']==fofID)[:][0]].flat[0]!=-99 and strikes<2) or use_inter==True:
                                try:
                                    if data['t70_'+key+'_std'][np.where(data['fofID']==fofID)[:][0]].flat[0]>0.05 and abs(f_inter_70-f_cs(0.7))>0.05:
                                        data['t70_'+key][np.where(data['fofID']==fofID)[:][0]] = cd.lookback_time(f_fit(0.7), z0=0.0, **fidcosmo)/3.1536e+16 #seconds to Gyrs   
                                        data['flag_std_fit_t70_'+key][np.where(data['fofID']==fofID)[:][0]] = 1
                                    else:
                                        data['flag_std_fit_t70_'+key][np.where(data['fofID']==fofID)[:][0]] = 0
                                except:
                                    pass
                            else:
                                data['t70_'+key][np.where(data['fofID']==fofID)[:][0]] = -99                                
    
                            #header+='\t std flag: '+str(data['flag_std_fit_t50_'+key][np.where(data['fofID']==fofID)[:][0]].flat[0])+'/'+str(data['flag_std_fit_t70_'+key][np.where(data['fofID']==fofID)[:][0]].flat[0])
                            if data_fit[0]>threshold_data_fit_z0:
                                array_unrealistic_fit.append(fofID)
                                
                            if data['t50_'+key][np.where(data['fofID']==fofID)[:][0]] < data['t70_'+key][np.where(data['fofID']==fofID)[:][0]] or abs(data['t50_'+key][np.where(data['fofID']==fofID)[:][0]] - data['t70_'+key][np.where(data['fofID']==fofID)[:][0]]) < 0.05:
                                array_low_delta_tform.append(fofID)
                                
                            if data['t50_'+key][np.where(data['fofID']==fofID)[:][0]]==-99 and data['t70_'+key][np.where(data['fofID']==fofID)[:][0]]==-99:
                                data['delta_tform_'+key][np.where(data['fofID']==fofID)[:][0]] = -99
                            else:
                                if data['t50_'+key][np.where(data['fofID']==fofID)[:][0]] < data['t70_'+key][np.where(data['fofID']==fofID)[:][0]] or abs(data['t50_'+key][np.where(data['fofID']==fofID)[:][0]] - data['t70_'+key][np.where(data['fofID']==fofID)[:][0]]) < 0.05:
                                    if strikes<2 and data_fit[0]<threshold_data_fit_z0:
                                        #print 'FIT IS BETTER OPTION!',format(data['t50_'+key][np.where(data['fofID']==fofID)[:][0]].flat[0], '0.3f'), format(data['t70_'+key][np.where(data['fofID']==fofID)[:][0]].flat[0], '0.3f')
                                        data['t50_'+key][np.where(data['fofID']==fofID)[:][0]] = cd.lookback_time(f_fit(0.5), z0=0.0, **fidcosmo)/3.1536e+16 #seconds to Gyrs 
                                        data['t70_'+key][np.where(data['fofID']==fofID)[:][0]] = cd.lookback_time(f_fit(0.7), z0=0.0, **fidcosmo)/3.1536e+16 #seconds to Gyrs 
                                        data['delta_tform_'+key][np.where(data['fofID']==fofID)[:][0]] = data['t50_'+key][np.where(data['fofID']==fofID)[:][0]] - data['t70_'+key][np.where(data['fofID']==fofID)[:][0]]
                                        header+=' use fit! '
                                        data['flag_use_fit_'+key][np.where(data['fofID']==fofID)[:][0]] = 1
                                        array_fit_prefered.append(fofID)
                                    else:
                                        
                                        data['delta_tform_'+key][np.where(data['fofID']==fofID)[:][0]] = -99
        
                                else:
                                    data['delta_tform_'+key][np.where(data['fofID']==fofID)[:][0]] = data['t50_'+key][np.where(data['fofID']==fofID)[:][0]] - data['t70_'+key][np.where(data['fofID']==fofID)[:][0]]
    
                            data['nr_strikes_'+key][np.where(data['fofID']==fofID)[:][0]] = strikes
                            header+='\n'
                            
                            if data['delta_tform_'+key][np.where(data['fofID']==fofID)[:][0]]==-99:
                                array_excluded_fofIDz0.append(fofID)
                            
                        except:
                            array_excluded_fofIDz0.append(fofID)
                            
                    # myOutput = oD.OutputData(config=False)
                    # myOutput.writeIntoFile(mycomp+'/anaconda/pro/data/EAGLE_ev_100Mpc/fit_inter_test/Interpolation_test_FoFID_'+str(fofID)+'.txt',
                    #             bucket,
                    #             myheader= header+'z\tfit_mhalo\tint_mhalo\tcs_mhalo\tmhalo_z/z0\tfit_mstar\tint_mstar\tcs_mstar\tmstar_z/z0\tfit_mb\tint_mbh\tcs_mbh\tmbh_z/z0',
                    #             data_format="%0.3f",
                    #             mydelimiter='\t')
                
                else:
                    for prop, key in zip(['bh_acc_rate', 'rVmax', 'angM_norm_stars', 'angM_norm_SFgas', 'angM_norm_NSFgas', 'sfr_30kpc', 'ssfr_30kpc', 'SHMR_30kpc', 'vmax', 'vdisp_30kpc', 'mgas_30kpc', 'age_mean_stars'],['bh_acc_rate', 'rVmax', 'angM_norm_stars', 'angM_norm_SFgas', 'angM_norm_NSFgas', 'sfr', 'ssfr', 'SHMR', 'vmax', 'vdisp', 'mgas', 'age']):
                        
                        #print '\t min/max --> prop:', prop, 'key:', key
                        data2process = data_ev_fofID[np.where(data_ev_fofID[prop]>0.0)[:][0]] 
                                                                            
                        if prop.startswith('angM') or prop.startswith('vdis'):
                            data2process = data2process[np.where(data2process[prop]>10.0)[:][0]]
                        elif prop=='rVmax':
                                data2process = data2process[np.where(data2process[prop]>2.0)[:][0]]                              
                        else:
                            #make sure that an non-exsiting vale (aka 0.0 at a certain redshift) is set as the minimum
                            data2process = data2process[np.where(data2process[prop]>0.0)[:][0]]
                            
                        data2process.sort(order=[prop], axis=0)
                        #print data2process[['fofID', 'redshift', prop]]
                                           
                        try:
                            data['tmin_'+key][np.where(data['fofID']==fofID)[:][0]] = cd.lookback_time(data2process['redshift'][0], z0=0.0, **fidcosmo)/3.1536e+16 #seconds to Gyrs                           
                            data['tmax_'+key][np.where(data['fofID']==fofID)[:][0]] = cd.lookback_time(data2process['redshift'][len(data2process['redshift'])-1], z0=0.0, **fidcosmo)/3.1536e+16 #seconds to Gyrs
                        
                            data['delta_tvar_'+key][np.where(data['fofID']==fofID)[:][0]] = data['tmax_'+key][np.where(data['fofID']==fofID)[:][0]] - data['tmin_'+key][np.where(data['fofID']==fofID)[:][0]]
                        except:
                            pass
            else:
                print ' short tree detected!\n',                   
                array_short_trees_fofIDz0.append(fofID)
                array_excluded_fofIDz0.append(fofID)
        
        if key_calc=='masses':                       
            print '\n-----------------------\nNr of short merger trees:\t', len(array_short_trees_fofIDz0), '\nIDs:', array_short_trees_fofIDz0
            print '\n-----------------------\nNr of np.interp errors:\t', len(np.unique(array_error_inter)), '\nIDs:', np.unique(array_error_inter)
            print '\n-----------------------\nNr of curve fit errors:\t', len(np.unique(array_error_curve_fit)), '\nIDs:', np.unique(array_error_curve_fit)[0:1000], np.unique(array_error_curve_fit)[1000::]   
            print '\n-----------------------\nNr of cs errors:\t\t', len(np.unique(array_error_cs)), '\nIDs:', np.unique(array_error_cs)
            print '\n-----------------------\nNr where inter was prefered over fit', len(np.unique(array_inter_prefered)), '\nIDs:', np.unique(array_inter_prefered)       
            print '\n-----------------------\nNr where fit was prefered over inter', len(np.unique(array_fit_prefered)), '\nIDs:', np.unique(array_fit_prefered)       
            print '\n-----------------------\nNr to high discuntinuity for first threee snapshots', len(np.unique(array_variation_excluded)), '\nIDs:', np.unique(array_variation_excluded)       
            print '\n-----------------------\nNr fit at z=0 exceed true values more than 50%', len(np.unique(array_unrealistic_fit)), '\nIDs:', np.unique(array_unrealistic_fit)
            print '\n-----------------------\nNr t50<t70 or very close together', len(np.unique(array_low_delta_tform)), '\nIDs:', np.unique(array_low_delta_tform)
           
            
            unique_fit_array    = np.unique(array_fit_prefered)
            test_both_corr      = np.in1d(unique_fit_array,np.unique(array_inter_prefered),assume_unique=True)
            
            print '\n-----------------------\nNr where inter was set first and then corrected to fit: ', len(unique_fit_array[np.where(test_both_corr==True)[:][0]]), '\nIDs:', unique_fit_array[np.where(test_both_corr==True)[:][0]]
           
            for mass in ['halo', 'stars', 'bh']:
                for key in ['50', '70']:
                    print 't'+key+'_'+mass, data[np.where(data['t'+key+'_'+mass]==-99)[:][0]].size
                    
                print 'delta_tform_'+mass+':', data[np.where(data['delta_tform_'+mass]==-99)[:][0]].size
                
            print 'Ngals where one of the delta_tforms is -99:', data[(np.where((data['delta_tform_halo']==-99) | (data['delta_tform_stars']==-99) |(data['delta_tform_bh']==-99)))[:][0]].size


            array_thalf_exluded = []
            
            for key in ['50', '70']:
                #print 'key:', key
                data['delta_t'+key+'_st2ha'] = -99.0
                data['flag_t'+key+'_st2ha']   = -99
                for l,value_halo, value_stars in zip(np.arange(data['t'+key+'_halo'].size),data['t'+key+'_halo'],data['t'+key+'_stars']):
                    #print l, value_halo, value_stars
                    if value_halo > 0.0 and value_stars>0.0:
                        data['delta_t'+key+'_st2ha'][l] = value_halo - value_stars
                        if data['delta_t'+key+'_st2ha'][l]>0.0:
                            data['flag_t'+key+'_st2ha'][l] = 1
                        else:
                            data['flag_t'+key+'_st2ha'][l] = 0
                    else:
                        array_thalf_exluded.append(fofID)
                        array_excluded_fofIDz0.append(fofID)
                
                data['delta_t'+key+'_bh2ha']  = -99.0
                data['flag_t'+key+'_bh2ha']   = -99
                for l, fofID, value_halo, value_bh in zip(np.arange(data['t'+key+'_halo'].size), data['fofID'], data['t'+key+'_halo'], data['t'+key+'_bh']):
                    #print l, value_halo, value_stars            
                    if value_halo>0.0 and value_bh>0.0:
                        data['delta_t'+key+'_bh2ha'][l] = value_halo - value_bh
                        if data['delta_t'+key+'_bh2ha'][l]>0.0:
                            data['flag_t'+key+'_bh2ha'][l] = 1
                        else:
                            data['flag_t'+key+'_bh2ha'][l] = 0
                    else:
                        array_thalf_exluded.append(fofID)
                        array_excluded_fofIDz0.append(fofID)
                        
                data['delta_t'+key+'_bh2st']  = -99.0
                data['flag_t'+key+'_bh2st']   = -99
                for l, fofID, value_stars, value_bh in zip(np.arange(data['t'+key+'_stars'].size), data['fofID'], data['t'+key+'_stars'], data['t'+key+'_bh']):
                    #print l, value_halo, value_stars            
                    if value_stars>0.0 and value_bh>0.0:
                        data['delta_t'+key+'_bh2st'][l] = value_stars - value_bh
                        if data['delta_t'+key+'_bh2st'][l]>0.0:
                            data['flag_t'+key+'_bh2st'][l] = 1
                        else:
                            data['flag_t'+key+'_bh2st'][l] = 0
                    else:
                        array_thalf_exluded.append(fofID)
                        array_excluded_fofIDz0.append(fofID)
                        
            print '\n-----------------------\nNr of not-exsiting t50/t70 --> finally excluded from merger tree study!', len(np.unique(array_thalf_exluded)), '\nIDs:', np.unique(array_thalf_exluded)[0:1000], np.unique(array_thalf_exluded)[1000::]
            print '\n++++++++++++++++++++++++++++++\nTotal Nr of excluded IDs at z=0:', len(np.unique(array_excluded_fofIDz0)), '\nIDs:', np.unique(array_excluded_fofIDz0)       
    
    
            print '\n++++++++++++++++++++++++++++++\n'
        
            data = data[np.in1d(data['fofID'], np.unique(array_excluded_fofIDz0), invert=True)]
            print 'final lenght data:', data.size 
            
     
            print '\n--> DONE!\n---------------------------------'
            
        else:
            
            new_props=['angM_norm_NSFgas','angM_norm_SFgas','angM_norm_stars', 'ssfr_30kpc', 'SHMR_30kpc', 'Tcons', 'cgf', 'bheff', 'SHMR']
            data[new_props]=-99       
            
            parentIndices, index1, index2 = np.intersect1d(data['haloid'], data2cross['haloid'], return_indices=True) 
            #print data2cross[new_props][0:25]
            
            for prop in new_props:
                #print 'property:', prop
                data[prop][index1]=data2cross[prop][index2]

        return data

    def FilterData(data):            
  
        #data=data[np.where(data['mPart_stars']>=1e9)[:][0]]
        #data=data[np.where(data['rhalf_stars_30kpc']>=0)[:][0]]
                
        #data=data[np.where(data['SB_mu_eff_stars_1.5ropt_r']>0.0)[:][0]]

        # data['flag_SB_B'] = -99
        # data['flag_SB_B'][(np.where((data['SB_mu_1.5opt_B']<25.6) & (data['SB_mu_1.5opt_B']>0.0)))[:][0]]=0
        # data['flag_SB_B'][np.where(data['SB_mu_1.5opt_B']>=25.6)[:][0]]=1

        # data['flag_SB_r'] = -99        
        # data['flag_SB_r'][(np.where((data['SB_mu_eff_stars_30kpc_r']<22.0) &  (data['SB_mu_eff_stars_30kpc_r']>0.0)))[:][0]]=0
        # data['flag_SB_r'][np.where(data['SB_mu_eff_stars_30kpc_r']>=22.0)[:][0]]=1

        data = data[np.where(data['flag_SB_B']!=-99)[:][0]]
        # #put inner and outer void galaxie into one category
        #data = data[np.where(data['envr']>=0)[:][0]]
        #data['envr'][np.where(data['envr']==0)[:][0]] = 1
        
        data['flag_tform_stars'] = -99
        data['flag_tform_stars'][(np.where((data['t50_stars']>=4) & (data['t50_stars']<7)))[:][0]]= 1
        data['flag_tform_stars'][(np.where((data['t50_stars']>=7) & (data['t50_stars']<10)))[:][0]]= 0
        
        data['flag_tform_bh'] = -99
        data['flag_tform_bh'][(np.where((data['t50_bh']>=2) & (data['t50_bh']<7)))[:][0]]= 1
        data['flag_tform_bh'][(np.where((data['t50_bh']>=7) & (data['t50_bh']<12)))[:][0]]= 0
        
        data['flag_tform_halo'] = -99
        data['flag_tform_halo'][(np.where((data['t50_halo']>=6) & (data['t50_halo']<9)))[:][0]]= 1
        data['flag_tform_halo'][(np.where((data['t50_halo']>=9) & (data['t50_halo']<12)))[:][0]]= 0
                
        #data['flag_mhalo'] = -99
        #data['flag_mhalo'][(np.where((data['mhalo_200c']>=10**11.5) & (data['mhalo_200c']<10**11.7)))[:][0]]= 0
        #data['flag_mhalo'][(np.where((data['mhalo_200c']>=10**11.7) & (data['mhalo_200c']<10**11.9)))[:][0]]= 1

        #for flag in ['flag_mhalo', 'flag_tform_stars', 'flag_tform_bh', 'flag_SB_B']:
         #   print 'flag:', flag, data[np.where(data[flag]==0)[:][0]].size, data[np.where(data[flag]==1)[:][0]].size, data[np.where(data[flag]-99)[:][0]].size
        #xit()

        # import numpy.lib.recfunctions as rcfuncs
        # data = rcfuncs.append_fields([data], 'environment', data['envr'], usemask=False, dtypes='S1')
        # data = rcfuncs.append_fields([data], 'flag_mass_bin', data['envr'], usemask=False, dtypes='S5')
        # print data.size
        # seperate_by_bins = True
        # data['flag_dens'] = 0
        # #selection den4 with SHMR (galaxy halo connection fixed)
        # bins = np.linspace(-2.6, -1.4, 49)

        # data['SHMR'] = data['mstar_30kpc']/data['mhalo_200c']

        # data = select_sample_dens(data, bins, 'SHMR')



        print data.size
        data['flag_dens'] = 0
        
        #selection dens3 with M200c
        #bins = np.linspace(11, 13.5, 51)
        
        #selection dens3 with mstar
        bins = np.linspace(9, 11.3, 47)
        #print min(data['mstar_30kpc']), max(data['mstar_30kpc'])
        print bins

        data = select_sample_dens(data, bins, 'mstar_30kpc')
        
        data = data[np.where(data['flag_dens'] == 1)[:][0]]


        # separation_dict={'flag_mstar_bin': {'prop': 'mstar_30kpc',   'bins_low': [9.4, 9.7, 10.0],   'bins_high': [9.7, 10.0, 10.3]},
        #                   'flag_mhalo_bin': {'prop': 'mhalo_200c',      'bins_low': [11.1, 11.4, 11.7], 'bins_high': [11.4, 11.7, 12.1]},
        #                   'flag_sfe_bin':   {'prop': 'sfe_gas_1.5ropt', 'bins_low': [-10.3, -10.0, -9.7], 'bins_high': [-10.0, -9.7, -9.3]},
        #                   'flag_age_bin':   {'prop': 't50_stars',       'bins_low': [4.0, 6.0, 8.0],   'bins_high': [6.0, 8.0, 10.0]},
        #                   'flag_SHMR_bin':  {'prop': 'SHMR_30kpc',    'bins_low': [-2.2, -2.0, -1.8],   'bins_high': [-2.0, -1.8, -1.6]}}
 
    
        # dict_data_print2table={}
        # for myflag in separation_dict.keys():
        #     myflag_prop=separation_dict[myflag]['prop']
        #     print 'myflag:', myflag, 'myflag_prop:', myflag_prop
        #     data[myflag]='None'
        #     list_samples=[]
        #     for SB, flag_SB in zip(['HSB', 'LSB'], [0,1]):
        #         for mybin, bin_low, bin_high in zip(['low', 'inter', 'high'], separation_dict[myflag]['bins_low'], separation_dict[myflag]['bins_high']): #sfe_HIH2_30kpc
        #             #sample_name = SB+'-'+mybin
        #             sample_name = mybin
        #             print SB, 'SB_flag:', flag_SB, 'sample_name:', sample_name, 'low/high bin:', bin_low, ' / ', bin_high,
                    
        #             if myflag.find('age')!=-1:
        #                 data[myflag][(np.where((data['flag_SB_B']==flag_SB) & (data[myflag_prop]>bin_low) & (data[myflag_prop]<=bin_high)))[:][0]] = sample_name                         
        #             else:
        #                 #print 'HERE:', myflag_prop, min(np.log10(data[myflag_prop])), max(np.log10(data[myflag_prop]))
        #                 data[myflag][(np.where(((data['flag_SB_B'])==flag_SB) & (np.log10(data[myflag_prop])>bin_low) & (np.log10(data[myflag_prop])<=bin_high)))[:][0]] = sample_name 

        #             print 'ngals: ', data[np.where(data[myflag]==sample_name)[:][0]].size
                    
        #             list_samples.append(sample_name)
                    
        return data

    def select_sample_dens(data,
                            bins,
                            prop,
                            sample_flag='flag_SB_B'):
        
        print 'bin data and select same number density for two different samples!\nprop:', prop, 'sample_flag:', sample_flag, '\nbins:', bins
        
        j=0
        while j < len(bins)-1:
            print 'j:', j, 'bin: [', bins[j], ',', bins[j+1], ']',
            data_in_bin = data[(np.where((np.log10(data[prop])>=bins[j]) & (np.log10(data[prop])<bins[j+1]) & (data[prop]>=0)))[:][0]]
            nHSB = data_in_bin[np.where(data_in_bin[sample_flag]==0)[:][0]].size
            nLSB = data_in_bin[np.where(data_in_bin[sample_flag]==1)[:][0]].size 

            delta_ngals = nLSB-nHSB
            
            
            print 'ngals:', data_in_bin.size, 'HSBs:', nHSB, 'LSBs:', nLSB, 'delta:', delta_ngals
            
            if delta_ngals!=0:
                flag=1
                if delta_ngals<0:
                    flag=0
                    
                random_ids = choose_random_sample(data_in_bin['haloid'][np.where(data_in_bin[sample_flag]==flag)[:][0]], abs(delta_ngals))
    
                #Equal numbers of galaxies in both samples!
                data_in_bin = data_in_bin[np.in1d(data_in_bin['haloid'], random_ids, invert=True)]
                
                #print '--> after LSB:', data_in_bin[np.where(data_in_bin['flag_SB_B']==1)[:][0]].size, 'HSB:', data_in_bin[np.where(data_in_bin['flag_SB_B']==0)[:][0]].size
                
                parentIndices, index1, index2 = np.intersect1d(data['haloid'], data_in_bin['haloid'], return_indices=True) 
    
                data['flag_dens'][index1]=1
                
            j+=1
        return data

    def AssignSubsamples(data):
        
        data['flag_SB+SFE'] = -99
        data['flag_SB+SFE'][(np.where((data['flag_SB_B']==1) & (data['flag_SFE']==1)))[:][0]]=1 #'BD'
        data['flag_SB+SFE'][(np.where((data['flag_SB_B']==1) & (data['flag_SFE']==0)))[:][0]]=2 #'BE'        
        data['flag_SB+SFE'][(np.where((data['flag_SB_B']==0) & (data['flag_SFE']==1)))[:][0]]=3 #'CD'         
        data['flag_SB+SFE'][(np.where((data['flag_SB_B']==0) & (data['flag_SFE']==0)))[:][0]]=4 #'CE'

        data['flag_sample'] = -99       
        data['flag_sample'][(np.where((data['flag_SB_B']==1) & (data['flag_SFE']==1) & (data['envr']==3)))[:][0]]=0 #'BD-S
        data['flag_sample'][(np.where((data['flag_SB_B']==1) & (data['flag_SFE']==1) & (data['envr']==2)))[:][0]]=1 #'BD-W
        data['flag_sample'][(np.where((data['flag_SB_B']==1) & (data['flag_SFE']==1) & (data['envr']==1)))[:][0]]=2 #'BD-V        
    
        data['flag_sample'][(np.where((data['flag_SB_B']==1) & (data['flag_SFE']==0) & (data['envr']==3)))[:][0]]=3 #'BE-S'
        data['flag_sample'][(np.where((data['flag_SB_B']==1) & (data['flag_SFE']==0) & (data['envr']==2)))[:][0]]=4 #'BE-W'
        data['flag_sample'][(np.where((data['flag_SB_B']==1) & (data['flag_SFE']==0) & (data['envr']==1)))[:][0]]=5 #'BE-V'

        
        data['flag_sample'][(np.where((data['flag_SB_B']==0) & (data['flag_SFE']==1) & (data['envr']==3)))[:][0]]=6 #'CD-S'
        data['flag_sample'][(np.where((data['flag_SB_B']==0) & (data['flag_SFE']==1) & (data['envr']==2)))[:][0]]=7 #'CD-W'
        data['flag_sample'][(np.where((data['flag_SB_B']==0) & (data['flag_SFE']==1) & (data['envr']==1)))[:][0]]=8 #'CD-V'
         
        data['flag_sample'][(np.where((data['flag_SB_B']==0) & (data['flag_SFE']==0) & (data['envr']==3)))[:][0]]=9 #'CE-S'
        data['flag_sample'][(np.where((data['flag_SB_B']==0) & (data['flag_SFE']==0) & (data['envr']==2)))[:][0]]=10 #'CE-W'
        data['flag_sample'][(np.where((data['flag_SB_B']==0) & (data['flag_SFE']==0) & (data['envr']==1)))[:][0]]=11 #'CE-V'
        
        data['flag_SB_envr'] = -99                      
        data['flag_SB_envr'][(np.where((data['flag_SB_B']==1) & (data['envr']==3)))[:][0]]=0 #'B-S
        data['flag_SB_envr'][(np.where((data['flag_SB_B']==1) & (data['envr']==2)))[:][0]]=1 #'B-W
        data['flag_SB_envr'][(np.where((data['flag_SB_B']==1) & (data['envr']==1)))[:][0]]=2 #'B-V
        
        data['flag_SB_envr'][(np.where((data['flag_SB_B']==0) & (data['envr']==3)))[:][0]]=3 #'C-S
        data['flag_SB_envr'][(np.where((data['flag_SB_B']==0) & (data['envr']==2)))[:][0]]=4 #'C-W
        data['flag_SB_envr'][(np.where((data['flag_SB_B']==0) & (data['envr']==1)))[:][0]]=5 #'C-V

        return data

    def GenerateBoxPlot(data):
                
        list_samples        = ['B-S', 'B-W', 'B-V', 'C-S', 'C-W', 'C-V']
        list_sample_keys    = [0, 1, 2, 3, 4, 5]
 
        #boxplot  environment with hue LSB/HSB
        list_samples        = ['S', 'W', 'V']
        list_sample_keys    = [3,2,1]
        
        
        #boxplot  environment with hue mass bins
        #list_samples        = ['low', 'inter', 'high']
        #list_sample_keys    = ['low', 'inter', 'high']
        
        list_samples        = ['LSB', 'HSB']
        list_sample_keys    = [1,0]
 
        #boxplot  environment with hue mass bins
        #list_samples        = ['HSB-low', 'LSB-low', 'HSB-inter', 'LSB-inter', 'HSB-high', 'LSB-high'] 
        #list_sample_keys    = list_samples
               
 
        # for item in list_samples:
        #     print item, data['flag_mass_bin'][np.where(data['flag_mass_bin']==item)].size, format(min(data['mstar_1.5ropt'][np.where(data['flag_mass_bin']==item)]), '0.2e'), format(max(data['mstar_1.5ropt'][np.where(data['flag_mass_bin']==item)]), '0.2e')
              
        
 
        myprops=get_list_all_properties(data)       
 
        myprops=['mhalo_200c','mstar_1.5ropt', 'SHMR_1.5ropt', 'mbh', 'bheff', 'cgf_1.5ropt', 'sfr_1.5ropt', 'ssfr_1.5ropt', 'Tcons_1.5ropt',\
                  'reff_gas_1.5ropt', 'reff_gas_disk_1.5ropt', 'rhalf_stars_1.5ropt', 'sfe_gas_reff_disk_1.5ropt', 'sfe_gas_1.5ropt',\
                  'Mgas_HI_30kpc_GK11','Mgas_H2_30kpc_GK11','fatom_30kpc_GK11', 'fmol_30kpc_GK11', 'OH_gas_30kpc',\
                  'Sigma_sfr_1.5ropt', 'Sigma_stars_1.5ropt', 'Sigma_gas_reff_disk_1.5ropt',\
                  'lastMinorM', 'lastMajorM',\
                  'age_stars_rband_r502D', 'age_stars_rband_2r502D', 'delta_age_stars_rband', 't50_stars', 'delta_tform_stars', 't50_halo', 'delta_tform_halo', 't50_bh', 'delta_tform_bh',\
                  'delta_t50_st2ha', 'delta_t50_bh2ha',\
                  'r2voidCenter',\
                  'r200c', 'v200c', 'vmax', 'rVmax', 'vdisp_30kpc_T19', 'Vrot2Vdisp_30kpc_T19', 'NFW_con',\
                  'angM_norm_1.5ropt_stars', 'angM_norm_1.5ropt_gas', 'angM_stars', 'angM_SFgas', 'angM_NSFgas', 'mAB_dA_total_cut_g_i', 'mAB_dA_total_cut_B_V']        
        myprops=myprops[::-1]
            
        fix_ndensity=False
        myflag_prop=['flag_age_bin']
        if fix_ndensity==True:
            for myflag in myflag_prop:
                data = data[np.where(data[myflag]!='0')[:][0]]
            
                for bin_flag in ['low','inter', 'high']:
    
                    delta_ngal = data[np.where(data[myflag]=='HSB-'+bin_flag)[:][0]].size - data[np.where(data[myflag]=='LSB-'+bin_flag)[:][0]].size
                                                                                              
                    print   'prop:', myflag, 'flag:', bin_flag,\
                            'LSB:', data[np.where(data[myflag]=='LSB-'+bin_flag)[:][0]].size,\
                            'HSB:', data[np.where(data[myflag]=='HSB-'+bin_flag)[:][0]].size,\
                              delta_ngal,
                    if  delta_ngal< 0:
                        prefix_through_out='LSB'
                    else:
                        prefix_through_out='HSB'
                    
                    random_ids = choose_random_sample(data['haloid'][np.where(data[myflag]==prefix_through_out+'-'+bin_flag)[:][0]], abs(delta_ngal))
        
                    #Equal numbers of galaxies in both samples!
                    data = data[np.in1d(data['haloid'], random_ids, invert=True)]
                    
                    print '--> after LSB:', data[np.where(data[myflag]=='LSB-'+bin_flag)[:][0]].size,\
                          'HSB:', data[np.where(data[myflag]=='HSB-'+bin_flag)[:][0]].size   
        
        stats_test_samples(data,
                           myprops,
                           list_samples,
                           list_sample_keys,
                           print_stats_table=False,
                           prefix_flag_name='',
                           is_flag_bool=False,
                           custom_bool_flag='flag_SB_B',
                           my_hue='flag_SB_B',
                           my_flag_name='flag_SB_envr',
                           box_plot_output_key='LSB_vs_HSB_envr_samples')
                           #box_plot_output_key='LSB_vs_HSB_envr_all_age_bin_samples_fixed_ndensity',
                           #mybin_flag='flag_age_bin')
        
        
        
        exit()
        
        
        

        
    def Print2TableBins(data):


        myprops=get_list_all_properties(data)
        
        myprops=myprops[0:3]

        envr_props={-99: 'EAGLE'}
        envr_prop_name=False
        
        data2print={}
        
        list_samples        = ['HSB-low', 'LSB-low', 'HSB-inter', 'LSB-inter', 'HSB-high', 'LSB-high'] 

        mybin_flag = 'flag_mhalo_bin'
        
        for item in ['HSB', 'LSB']:
            for mybin in ['low', 'inter', 'high']:
                bin_flag = item+'-'+mybin                 
                mydata = data[np.where(data[mybin_flag]==bin_flag)[:][0]]
                data2print.update({bin_flag: mydata})
                
        print "{0:.2e}".format(min(data2print['HSB-low']['mhalo_200c'])), "{0:.2e}".format(max(data2print['HSB-low']['mhalo_200c']))
        print "{0:.2e}".format(min(data2print['HSB-low']['mstar_1.5ropt'])), "{0:.2e}".format(max(data2print['HSB-low']['mstar_1.5ropt']))
                
        H, xedges, yedges = np.histogram2d(data2print['HSB-low']['mhalo_200c'], data2print['HSB-low']['mstar_1.5ropt'], bins=20)      
        
        myFuncs = mF.MyFunctions()   
        nbins=20
        histo, binsize = myFuncs.binUp(data2print['HSB-low'][['mhalo_200c','mstar_1.5ropt']],
                                            nbins,
                                            histo_min=min(data2print['HSB-low']['mhalo_200c']),
                                            histo_max=max(data2print['HSB-low']['mhalo_200c']),
                                            log10bin=True,
                                            binup_2D=True,
                                            div_y2x=False,
                                            use_MAD=True,
                                            equal_bins=True)
        
        histo[:,2] = xedges[:-1]
        histo[:,3] = yedges[:-1]
        
                                                        
        print histo[:, [0,1,2,3]]
        print binsize                                                 
        
        exit()
                
                                                    
        print_props_table_format(data2print,
                                myprops=myprops,
                                envr_prop_name=envr_prop_name,
                                envr_props=envr_props,
                                percentiles=[25.0,75.0],
                                order_myprops=True,
                                plot_median=False,
                                order_data_dict=False,
                                output_filename=mycomp+'anaconda/pro/myRun/plots/plots_table_stats/Hidden_Figures_I/',
                                output_filename_key='boot_err')
        


        exit()

    def Print2Table(data):


       myprops=get_list_all_properties(data)

                   
       myprops=['mhalo_200c', 'mstar_1.5ropt', 'SHMR_1.5ropt', 'mbh', 'bheff', 'cgf_1.5ropt', 'ssfr_1.5ropt', 'Tcons_1.5ropt',\
                'reff_gas_disk_1.5ropt', 'rhalf_stars_1.5ropt', 'sfe_gas_reff_disk_1.5ropt',\
                'lastMerger', 'lastMinorM', 'lastMajorM', 'ratio_lastM', 'age_stars_rband_r502D',\
                'vmax', 'rVmax', 'vdisp_30kpc_T19', 'Vrot2Vdisp_30kpc_T19', 'DvT_stars_c0.5_1.5ropt', 'DvT_gas_c0.5_1.5ropt', 'angM_norm_1.5ropt_gas']

           
           
       # #2 sub-samples L/HSB in B-band
       # #evnr_prop_name='flag_SB_B'
       # #envr_props={-99: 'EAGLE', 1: 'B', 0: 'C'}
       
       # #2 sub-samples L/HSFE in star formation efficiency
       # #evnr_prop_name='flag_SFE'
       # #envr_props={-99: 'EAGLE', 1: 'D', 0: 'E'}
       
       # #4 sub-samples L/HSB & L/HSFE
       # evnr_prop_name='flag_sample'
       # envr_props={-99: 'EAGLE', 0: 'BD', 1: 'BE', 2: 'CD', 3: 'CE'}
       
       # #4 environments: 0=IV,1=OV,2=W,3=S
       evnr_prop_name='flag_SB_envr'
       envr_props={0: 'B-S', 1: 'B-W', 2: 'B-V', 3: 'C-S', 4: 'C-W', 5: 'C-V'}
       data2print={'EAGLE': data}
       
       #12  subsamples 
       #evnr_prop_name='flag_sample' 
       #envr_props={-99: 'EAGLE',   0: 'BD-S', 1: 'BD-W', 2: 'BD-V',\
                                   # 3: 'BE-S', 4: 'BE-W', 5: 'BE-V',\
                                   # 6: 'CD-S', 7: 'CD-W', 8: 'CD-V',\
                                   # 9: 'CE-S', 10: 'CE-W', 11: 'CE-V'}
                                   
       data_LSB   = data[np.where(data['flag_SB_B']==1)[:][0]]
       data_LSB_S = data_LSB[np.where(data_LSB['envr']==3)[:][0]]
       data_LSB_W = data_LSB[np.where(data_LSB['envr']==2)[:][0]]
       data_LSB_V = data_LSB[np.where(data_LSB['envr']==1)[:][0]]
       
       data_HSB   = data[np.where(data['flag_SB_B']==0)[:][0]]
       data_HSB_S = data_HSB[np.where(data_HSB['envr']==3)[:][0]]
       data_HSB_W = data_HSB[np.where(data_HSB['envr']==2)[:][0]]
       data_HSB_V = data_HSB[np.where(data_HSB['envr']==1)[:][0]]                        
       
       data2print={'B':  data_LSB,
                   'B-S': data_LSB_S,
                   'B-W': data_LSB_W,
                   'B-V': data_LSB_V,
                   'C':  data_HSB,
                   'C-S': data_HSB_S,
                   'C-W': data_HSB_W,
                   'C-V': data_HSB_V}

       envr_props={-99: 'EAGLE'}
       
       print_props_table_format(data2print,
                               myprops=myprops,
                               envr_prop_name=evnr_prop_name,
                               envr_props=envr_props,
                               percentiles=[25.0,75.0],
                               order_myprops=True,
                               plot_median=True,
                               order_data_dict=True,
                               output_filename=mycomp+'anaconda/pro/myRun/plots/plots_table_stats/Hidden_Figures_I/')
       


       exit()

    


    def caseSwitcher(item):                        
             
        choose = {
            'crossMatch':           CrossMatch,
            'calcAssemblyHistory':  CalcAssemblyHistory,
            'filterData' :          FilterData,
            'assignSubsamples':     AssignSubsamples,
            'print2Table':          Print2Table,
            'print2TableBins':      Print2TableBins,            
            'generateBoxPlot':      GenerateBoxPlot
            }
            
        func = choose.get(item)
        return func(data)
    
    return  caseSwitcher(key)

def handle_SFH_Galacticus_data(data,
                              key=None):  

    def ReadDataProps():
        
        import pandas as pd
        dt = dt_MDPL_Galacticus_CMASS_v2()
                
        data = pd.read_csv(mycomp+'anaconda/pro/data/Galacticus_1Gpc/Galacticus_1Gpc_z_0.56_tarsel_new_mags_CMASS_down_sample3_more_props_sample_keys.txt',\
                           names=[dt[c][0] for c in xrange(len(dt))], skiprows=2, sep='  ')
            
        
        return df_to_sarray(data) 
    
    def FilterData():
        pass
    
    def CrossMatch(data):
        
        data = crossmatch_catalogs(data, key='CMASS_new_envr', redshift=0.56)
        
        return data
    
    
    def AssignSubsamples(data):
        
        data, data2print, odict = assign_SFH_samples(data,
                                                     'CMASS_sample_key')
        
        i=0
        while i<data.size:
            data['nr_sample_keys'][i] = data['CMASS_sample_key'][i].count(':')-1
            i+=1
                
        return data          


    def Print2Table(data):
        
        #data = ReadDataProps()
        #data = AssignSubsamples(data)
        #print data['orphan']
        myprops=[]
        for prop in data.dtype.names:
            prop_dtype = data[prop].dtype
            if prop_dtype!='int64' and prop_dtype!='int8' and prop.find('pos')==-1 and prop.startswith('zgasSF')==False:
                #print prop, prop_dtype
                myprops.append(prop)
                           
        evnr_prop_name  = 'env_1024' 
        envr_props      =  {3: 'knots', 2: 'filaments'}
        
        #envr_props      =  {-99: 'high-zcold-red'}
        
        #evnr_prop_name='flag_sample'
        #envr_props={-99: 'Gal-dens', 8: 'low-zcold-red', 9: 'high-zcold-red'}
            
        data=data[np.where(data['orphan']==0)[:][0]]
        print data.size
        data_lzr = data[np.where(data['flag_CMASS_lzr']==1)[:][0]]
        data_hzr = data[np.where(data['flag_CMASS_hzr']==1)[:][0]]
        
        data_lzr_knots = data[np.where((data['flag_CMASS_lzr']==1) & (data['env_1024']==3))[:][0]]
        data_lzr_filaments = data[np.where((data['flag_CMASS_lzr']==1) & (data['env_1024']==2))[:][0]]
        
        data_hzr_knots = data[np.where((data['flag_CMASS_hzr']==1) & (data['env_1024']==3))[:][0]]
        data_hzr_filaments = data[np.where((data['flag_CMASS_hzr']==1) & (data['env_1024']==2))[:][0]]
        
        #data2print={'knots': data_hzr_knots, 'filaments': data_hzr_filaments}
        data2print={'Gal-dens': data, 'low-zcold-red': data_lzr, 'high-zcold-red': data_hzr}
      
        myprops=['mhalo','mstar','SHMR','mbh','mcold', 'mhot', 'Mzgas','Mzstar', 'Mzhot_halo', 'cgf','zcold', 'zstar',\
                 'zcold_zstar', 'mAB_dA_total_cut_g_i','mAB_dA_total_cut_r_i', 'sfr', 'ssfr', 'Tcons','bheff', 'NFW_con','jbar','vbulge', 'angM_spheroid', 'jbulge']

        print_props_table_format(data2print,
                                myprops=myprops,
                                envr_prop_name=evnr_prop_name,
                                envr_props=envr_props,
                                percentiles=[31.8,68.2],
                                order_myprops=False,
                                order_data_dict=True,
                                order_envr_props_dict=False,
                                plot_median=False,
                                output_filename=mycomp+'anaconda/pro/myRun/plots/plots_table_stats/SFR_History_I/',
                                output_filename_key='_subsamples')

        exit()

    def GenerateBoxPlot(data):
        
        data=data[np.where(data['orphan']==0)[:][0]]
        
        #data_low = data[np.where(data['flag_CMASS_l']==1)[:][0]]
        #data_red = data[np.where(data['flag_CMASS_r']==1)[:][0]]
        #data_lzr = data[np.where(data['flag_CMASS_lzr']==1)[:][0]]
        #data_hzr = data[np.where(data['flag_CMASS_hzr']==1)[:][0]]
        #data_dict_start={'red': data_red, 'low': data_low, 'lowZ': data_lzr, 'highZ': data_hzr}
        
        #myprops = ['mhalo']
        myprops=['mstar', 'mbh', 'mcold', 'mhot', 'Mzgas','Mzstar', 'Mzhot_halo', 'cgf', 'zcold', 'zstar','zcold_zstar', 'mAB_dA_total_cut_g_i','mAB_dA_total_cut_r_i', 'sfr', 'ssfr', 'SHMR', 'Tcons','bheff', 'NFW_con']
        
        #myprops=myprops[::-1]
             
        #list_samples        = ['red', 'low', 'low-zcold', 'high-zcold']
        #list_sample_keys    = ['r', 'l', 'lzr', 'hzr']
        
        #list_samples        = ['all', 'low', 'high', 'passive', 'active', 'red', 'blue', 'low-zcold', 'high-zcold', 'low-zcold-red', 'high-zcold-red']
        #list_sample_keys    = ['' , 'l', 'h', 'p', 'a', 'r', 'b', 'lz','hz', 'lzr', 'hzr']
    
        list_samples        = ['all', 'low', 'passive', 'red', 'low-zcold-red', 'high-zcold-red']
        list_sample_keys    = ['' , 'l', 'p', 'r', 'lzr', 'hzr']
                   
        stats_test_samples(data,
                           myprops,
                           list_samples,
                           list_sample_keys,
                           print_stats_table=True,
                           prefix_flag_name='flag_CMASS_')
            


    def CalcStats(data):

        data=data[np.where(data['orphan']==0)[:][0]]    
        data_lzr = data[np.where(data['flag_CMASS_lzr']==1)[:][0]]
        data_hzr = data[np.where(data['flag_CMASS_hzr']==1)[:][0]]
        
        data_lzr_knots = data[np.where((data['flag_CMASS_lzr']==1) & (data['env_1024']==3))[:][0]]
        data_lzr_filaments = data[np.where((data['flag_CMASS_lzr']==1) & (data['env_1024']==2))[:][0]]
        
        data_hzr_knots = data[np.where((data['flag_CMASS_hzr']==1) & (data['env_1024']==3))[:][0]]
        data_hzr_filaments = data[np.where((data['flag_CMASS_hzr']==1) & (data['env_1024']==2))[:][0]]
        
        stats_test_samples(data)
            

    def caseSwitcher(item):                         
             
        choose = {
            'crossMatch':       CrossMatch,
            'assignSubsamples': AssignSubsamples,
            'print2Table':      Print2Table,
            'calcStats':        CalcStats,
            'generateBoxPlot':  GenerateBoxPlot
            }
 
        func = choose.get(item)
        return func(data)
   
    return caseSwitcher(key)               


def stats_test_samples(data,
                       myprops,
                       list_samples,
                       list_sample_keys,
                       print_stats_table=True,
                       prefix_flag_name='',
                       box_plot_output_key='',
                       is_flag_bool=True,
                       custom_bool_flag='',
                       my_hue=None,
                       my_flag_name=None,
                       mybin_flag=None):
    
    """
    input:
        data        is structured array
        myprops     list of properties to calc statistics from or boxplots generated
        print_stats_table   "True", it prints statistics, "other", generates boxplots
        
    
    """
    
    import pandas as pd
    myFuncs = mF.MyFunctions()   

    sample_name_dict = get_SFH_sample_name_dict()
    myprop_unit_map  = get_prop_unit_map_latex()
    #print_stats_table=False
    
    
    data_dict={}
    for prop in myprops:

        #print '\nprop:', prop, '\n==========\n'   
        for item, flag_name in zip(list_samples, list_sample_keys):
            print 'item:', item, 'flag_name:', flag_name, is_flag_bool, prefix_flag_name,
            #data_item = choose_random_sample(data[np.where(data['flag_CMASS_'+flag_name]==1)[:][0]], nr_tests)
            latex_name = sample_name_dict[item]['name_Latex']
            print 'latex_name:', latex_name
            try:
                if flag_name=='custom':
                    data_item = data[np.where(data[custom_bool_flag]==1)[:][0]]
                elif is_flag_bool==True:
                    data_item = data[np.where(data[prefix_flag_name+flag_name]==1)[:][0]]
                else:
                    data_item = data[np.where(data[prefix_flag_name]==flag_name)[:][0]]
            except:
                data_item = data
        
            if myprop_unit_map[prop][1]=='lin':
                data_to_test = data_item[prop]
            else:
                data_to_test = np.log10(data_item[prop])
                
            mask=np.where(np.isfinite(data_to_test))
            data_to_test=data_to_test[mask[:][0]]

            data_to_test=data_to_test[np.where(data_to_test!=-99.0)[:][0]]            
            
            data_dict.update({latex_name: data_to_test})
            
 
        if print_stats_table==True:

            print '\n\n\\newpage\n\n\\begin{table*}\n\t\\begin{center}\n\t\setlength{\\tabcolsep}{3pt}\n\t\t\\begin{tabular}{l||c|c|c|c|c|c|c}'
            print_table_stats_boxplot(data_dict,
                                      list_samples,
                                      prop)
            print' \t\t\t\hline\n\t\t\t (i)	\t & (ii)\t & (iii)\t & (iv)\t  & (v) \t & \t  (vi) \t  &  \t (vii) \t  & \t  (viii) \t break \n\t\t\end{tabular}\n\t\t\caption{\DS{Table just for checking, will not be published!}.}\label{tab:props_subsamples_stats_boxplots_'+prop+'}\n\t\end{center}\n\end{table*}\n\n'
    
            print '\\begin{figure*}\n\t\centering\n\t\includegraphics[width=0.75\\textwidth,angle=0]{plots/boxplots/SFH_box_plot_reduced_samples_'+prop+'.pdf}\n\t\caption{\DS{Table just for checking, will not be published!}.}\label{fig:box_plot_mhalo}\n\end{figure*}\n\n'

        else:
            
            if my_hue==None:            
                df = pd.DataFrame({sample_name_dict[list_samples[0]]['name_Latex']: data_dict[sample_name_dict[list_samples[0]]['name_Latex']]})     
                for item in list_samples[1::]:
                    df_add = pd.DataFrame({sample_name_dict[item]['name_Latex']: data_dict[sample_name_dict[item]['name_Latex']]})
                    df = pd.concat([df, df_add], axis=1)
            else:
                #my_hue: flag_SB_B LSB vs. HSB in different environments
                #data_to_test = data[[prop, 'envr', mybin_flag]]
                
                #my_hue: flag_SB_B LSB vs. HSB in different mass bins
                #data_to_test = data[['flag_SB_B', prop, 'flag_mass_bin']]
                
                #my_hue: envr in LSB and HSB bins
                data_to_test = data[['envr', prop, 'flag_SB_B']]
                
                if myprop_unit_map[prop][1]!='lin':
                    data_to_test[prop] = np.log10(data_to_test[prop])
                    mask=np.where(np.isfinite(data_to_test[prop]))
                    data_to_test=data_to_test[mask[:][0]]
                    
                data_to_test=data_to_test[np.where(data_to_test[prop]!=-99.0)[:][0]]



                x_axis_name='$M_{*}$'
                #x_axis='$M_{200c}$'
                #my_box_plot_order=['low', 'inter', 'high']
                #my_hue='flag_SB_B' 

                x_axis_name='Surface Brightness'
                x_prop='flag_SB_B'
                my_box_plot_order=['HSB','LSB']                
                my_hue='envr'
                
                # x_axis_name='Environment'
                # x_prop='envr'
                # my_box_plot_order=['S','W', 'V']                
                # my_hue='flag_SB_B'
 
                # x_axis_name='Surface Brightness & Mass Bin'
                # x_prop=mybin_flag
                # my_box_plot_order        = ['HSB-low', 'LSB-low', 'HSB-inter', 'LSB-inter', 'HSB-high', 'LSB-high']            
                # my_hue='envr'

                import numpy.lib.recfunctions as rcfuncs
                data_to_test = rcfuncs.append_fields([data_to_test], x_axis_name, data_to_test['envr'], usemask=False, dtypes='S10')
                df = pd.DataFrame({})
                for item, flag in zip(list_samples,list_sample_keys):
                    print 'item:', item, 'flag:', flag
                    data_to_test[x_axis_name][np.where(data_to_test[x_prop]==flag)[:][0]]=str(item)             
              

                #df=pd.DataFrame({'SB': data_to_test['flag_SB_B'], 'environment': data_to_test['envr'], prop: data_to_test[prop]})
                df=pd.DataFrame({my_hue: data_to_test[my_hue], x_axis_name: data_to_test[x_axis_name], prop: data_to_test[prop]})
                #df=pd.DataFrame({my_hue: data_to_test[my_hue], x_axis: data_to_test['flag_SB_B'], prop: data_to_test[prop]})                

                #df=pd.melt(df)
              
                #df_test = df_new.melt(id_vars=['SB'], value_vars=['envr'], var_name='environment', value_name=prop)
                print df.head()                 
                #print df_test

            generate_box_plot(df,
                              prop,
                              len(my_box_plot_order),
                              output_key=box_plot_output_key,
                              my_hue=my_hue,
                              x_axis=x_axis_name,
                              my_box_plot_order=my_box_plot_order)      

def perform_statistical_tests(data_dict,
                              nr_tests,
                              prop='mstar'):
        
    print 'perform F-Test -->\n'
    myFuncs.Ftest(data_dict, nr_tests=nr_tests)
    print 'DONE!\n'

    print 'perform ANOVA -->\n'     
    myFuncs.ANOVA(data_dict, prop)
    print 'DONE!\n'
    
def print_table_stats_boxplot(data_dict,
                              list_samples,
                              prop):
    
    myprop_unit_map  = get_prop_unit_map_latex()
    sample_name_dict = get_SFH_sample_name_dict()
    
    #print data_dict
    print list_samples
    print prop, myprop_unit_map[prop][0], myprop_unit_map[prop][1]
    
    stats_props_boxplot='\t\t\tsample \t & \t $N_{\mathrm{gal}}$ \t & \t n [$times$ $10^{-4}$ Mpc$^{-3}$] \t & \t'+myprop_unit_map[prop][0]+'$^{\mathrm{+Q3}}_{\mathrm{-Q1}}$ \t & \t IQR \t & \t'+myprop_unit_map[prop][0]+'$^{\mathrm{Q3+1.5}\\times\mathrm{IQR}}_{\mathrm{Q1-1.5}\\times\mathrm{IQR}}$ & \t $N_{\mathrm{outliers}}$ \t & \%\t & $\sigma_{\mathrm{\bar{x}}}$ \\\\ \n \t\t\t \hline\n\t\t\t\hline'
    for item in list_samples:
        #print item, '\n------------'
        list_outliers, count_outliers, perc_outliers = count_outliers_boxplot(data_dict[sample_name_dict[item]['name_Latex']])
        #print 'count ourtliers:', count_outliers, '/ ', "{0:.2f}".format(perc_outliers), '% \n'
        
        stats_props_boxplot+='\t\t\t'+sample_name_dict[item]['name_Latex']+'\t & \t'+str(data_dict[sample_name_dict[item]['name_Latex']].size)+'\t & \t'+str("{0:.3f}".format(data_dict[sample_name_dict[item]['name_Latex']].size/(1000.0/0.6778)**3*1e4))+'\t & \t'

        median  = np.nanmedian(data_dict[sample_name_dict[item]['name_Latex']])
        Q1   = np.nanmedian(data_dict[sample_name_dict[item]['name_Latex']])-np.nanpercentile(data_dict[sample_name_dict[item]['name_Latex']], 25)
        Q3   = np.nanpercentile(data_dict[sample_name_dict[item]['name_Latex']], 75)-np.nanmedian(data_dict[sample_name_dict[item]['name_Latex']])
        IQR     = abs(Q3 - Q1)
        lower = Q1 - 1.5 * IQR
        upper = Q3 + 1.5 * IQR    

        #print 'median:', "{0:.2f}".format(median), '-/+', "{0:.2f}".format(lower), '/', "{0:.2f}".format(upper)
        stats_props_boxplot+=str("{0:.2f}".format(median))+'$_{-'+str("{0:.2f}".format(Q1))+'}^{+'+str("{0:.2f}".format(Q3))+'}$\t & \t '+str("{0:.2f}".format(IQR))+'\t & \t '+str("{0:.2f}".format(median))+'$_{-'+str("{0:.2f}".format(lower))+'}^{+'+str("{0:.2f}".format(upper))+'}$\t & \t '+str(count_outliers)+'\t & \t'+str("{0:.2f}".format(perc_outliers))+'\t & \t'   
        
        #print bootstrap errors
        boot_res, conf_intervals_boot_res, dist_boot_res = mF.bootstrap_error(data_dict[sample_name_dict[item]['name_Latex']])

        stats_props_boxplot+=str("{0:.2f}".format(boot_res))+'$_{-'+str("{0:.2f}".format(conf_intervals_boot_res[0]))+'}^{+'+str("{0:.2f}".format(conf_intervals_boot_res[1]))+'}$'+' \t \\\\ \n \hline \n'
        
        
    print stats_props_boxplot
    
    return stats_props_boxplot
    
                   
    #exit()

def generate_box_plot(data,
                      prop,
                      nr_samples,
                      envr='',
                      output_key='',
                      my_hue=None,
                      x_axis=False,
                      my_box_plot_order=None):

    #Load plot related packages
    import seaborn as sns
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.ticker import MultipleLocator
    from matplotlib.ticker import NullFormatter
    nullfmt = NullFormatter()
    mpl.style.use('classic')            
    mpl.mathtext.use_cm = False
    mpl.rcParams['mathtext.fontset'] = 'custom'
    mpl.rcParams['mathtext.tt'] = 'Typewriter'
    mpl.mathtext.fallback_to_cm = True

    
    if envr=='':
        envr_space=''
    else: 
        envr_space='_'
    if output_key=='':
        output_key_space=''
    else: 
        output_key_space='_'

    
    myprop_unit_map  = get_prop_unit_map_latex()

    prop_name_latex = myprop_unit_map[prop][0]
    
    sample_name_dict = get_SFH_sample_name_dict()
      
    output_filename=mycomp+'anaconda/pro/myRun/plots/plots_table_stats/SFR_History_I/SFH_box_plot_'+output_key+output_key_space+'samples_'+prop+envr_space+envr+'.pdf'
    output_filename=mycomp+'anaconda/pro/myRun/plots/plots_table_stats/Hidden_Figures_I/SFH_box_plot_'+output_key+output_key_space+'samples_'+prop+envr_space+envr+'_reduced.pdf'

    
    if nr_samples<=1:
        colours = ['k', 'r']
        x_size_proxy=3.5
    elif nr_samples<=3:
        #colours = ['#ffab00', '#332288', '#CC6677'] #envr
        colours = ['#117733', '#e65100', '#CC6677'] #envr
        #colours = ['#ffab00', '#332288', '#696969'] #HSB vs LSB
        x_size_proxy=3.5
    elif nr_samples>3 and nr_samples<=6:        
        colours = [ '#117733', '#e65100', '#CC6677', '#ffab00', '#AA4499', '#332288']
        #colours = ['#ffab00', '#332288', '#696969'] #m200c bin flag
        #colours = ['#332288', '#88CCEE', '#44AA99'] #age bin flag
        #olours = ['#CC6677', '#AA4499', '#7439A6'] #SFE gas 1.5ropt
        x_size_proxy=2.5
    else:
        colours = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', '#ffab00', '#CC6677', '#AA4499', '#7439A6', '#661100', '#333333', '#d2691e']
        x_size_proxy=1
    frame_lw = 4.0     

    fig = plt.figure(figsize=(x_size_proxy*nr_samples,9), dpi=150)
    ax = fig.add_subplot(111)
   
    if my_hue==None:
        check_max=max(data.max(axis='columns'))
        check_min=min(data.min(axis='columns'))
    else:
        check_max=max(data[prop])
        check_min=min(data[prop])       
    
    print 'prop:', prop, 'envr:', envr, 'check_max:', check_max, 'check_min:', check_min, '--> max-min:', check_max-check_min,

    if check_max>=10.0 and (check_max-check_min)>10.0:       
        float_format= 0
        y_minor     = 10
        ax.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda y, p: "{0:.0f}".format(y)))

    elif check_max>=10.0 and (check_max-check_min)>5.0:
        float_format= 0
        y_minor     = 0.5
        ax.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda y, p: "{0:.0f}".format(y)))

    elif check_max>=10.0 and (check_max-check_min)>1.0:
        float_format= 1
        y_minor     = 0.1
        ax.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda y, p: "{0:.1f}".format(y)))
    elif check_max>=10.0 and (check_max-check_min)>0.5:
        float_format= 1
        y_minor     = 0.1
        ax.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda y, p: "{0:.1f}".format(y)))
    elif check_max<10.0 and (check_max-check_min)>0.5:
        float_format= 1
        y_minor     = 0.5
        ax.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda y, p: "{0:.1f}".format(y))) 
    else:
        float_format= 2
        y_minor     = 0.05
        ax.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda y, p: "{0:.2f}".format(y))) 

    if my_hue==None:
        Q1 = np.nanpercentile(data.stack(), 20)
        Q3 = np.nanpercentile(data.stack(), 80)
    else:
        Q1 = np.nanpercentile(data[prop], 20)
        Q3 = np.nanpercentile(data[prop], 80)
        
    IQR = Q3 - Q1

    ymin                = round(Q1 - 1.5 * IQR, float_format)
    ymax                = round(Q3 + 1.5 * IQR, float_format)
    # 
    if prop_name_latex.find('Gyr')!=-1:
        ymax=15
        ymin = 0
        minor_y = 5
        ax.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda y, p: "{0:.0f}".format(y)))

    if nr_samples<=2:
        left=0.35
    elif nr_samples<4:
        left=0.5
    elif nr_samples<6 and nr_samples>=4:
        left=0.01*nr_samples
           
    right=0.98
    if my_hue==None:
        bottom=0.10
    else:
        bottom=0.15
        left=0.16 #6 samples
        left=0.28 #3 samples
    
    myflierprops = dict(marker='x', markerfacecolor='r', markersize=7,
                  linestyle='none', markeredgecolor='k')
    
    fig.subplots_adjust(hspace=0.01, wspace=0.01, left=left, bottom=bottom, right=right, top=0.94) 

    if prop.startswith('DvT'):
        ymin=0.0
        ymax=1.1
        y_minor=0.05
    elif prop.startswith('delta'):
        ymin=-5
        ymax=5
        y_minor=1        
    elif prop.startswith('cgf'):
        ymin=-1.6
        ymax=0.7
        y_minor=0.1       
    elif prop.startswith('MAB_'):
        ymin=-24
        ymax=-16
        y_minor=0.5
        ax.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda y, p: "{0:.0f}".format(y)))
    elif prop.startswith('Mgas_'):
        ymin=7
        ymax=10.25
        y_minor=0.25
        ax.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda y, p: "{0:.1f}".format(y)))
    elif prop.startswith('mgas_half'):
        ymin=8.
        ymax=9.65
        y_minor=0.1
        ax.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda y, p: "{0:.1f}".format(y)))
    elif prop.startswith('mstar'):
        ymin=9.0
        ymax=11.0
        y_minor=0.1
        ax.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda y, p: "{0:.1f}".format(y)))         
    elif prop=='mhalo_200c':
         ymin=10.9
         ymax=12.5
         y_minor=0.1
    elif prop=='SHMR_1.5ropt':
         ymin=-2.3
         ymax=-1.35
         y_minor=0.1
    elif prop=='ssfr_1.5ropt':
         ymin=-11.0
         ymax=-9.3
         y_minor=0.1
    elif prop=='ratio_lastM':
         ymin=0
         ymax=0.95
         y_minor=0.1
    elif prop=='angM_norm_1.5ropt_stars':
         ymin=-0.4
         ymax=1.05
         y_minor=0.1
    elif prop=='zcold_gasSF':
         ymin=9.95
         ymax=10.65
         y_minor=0.1
    elif prop=='Vrot2Vdisp_30kpc_T19':
         ymin=-0.1
         ymax=2.5
         y_minor=0.1
    elif prop=='vmax':
         ymin=80
         ymax=190
         y_minor=5
    elif prop.startswith('reff'):
         ymin=1
         ymax=25
         y_minor=0.5
    elif prop=='rVmax':
         ymin=0
         ymax=55
         y_minor=5
    elif prop=='r200c':
         ymin=50
         ymax=300
         y_minor=10
    elif prop.startswith('age'):
         ymin=0
         ymax=10
         y_minor=0.25
    elif prop.startswith('last'):
         ymin=0
         ymax=10
         y_minor=0.25
    elif prop.startswith('NFW'):
         ymin=7
         ymax=10
         y_minor=0.25
         ax.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda y, p: "{0:.1f}".format(y)))         
    
    minor_y = MultipleLocator(y_minor)
    ax.yaxis.set_minor_locator(minor_y)

    print 'ylim:', ymin, '/', ymax, 'minor_y:', y_minor#, 'yticks:',   
    yticks=[]
    l=0
    while l<len(ax.get_yticks()[:]):
        #print 'i:', i, ax0.get_xticks()[k]
        yticks+=[round(ax.get_yticks()[:][l],float_format)]
        l+=1       
        
    #print yticks     
    ax.set_xticks(yticks)
    ax.set_xticklabels(yticks)
    
    plt.ylim(ymin, ymax)
        
    #plt.ylim(10, 14)
    plt.setp(ax.get_yticklabels(), visible=True)    

    if my_hue==None:
        ax = sns.boxplot(data=data, linewidth=frame_lw,
                         palette=colours, medianprops={"color": "r", "linewidth": 4}, width=0.4, flierprops=myflierprops, showfliers=False)
        ax.set_xlabel("sample names", fontsize=30, labelpad=10)

    else:
        ax = sns.boxplot(x=x_axis, y=prop, hue=my_hue, order=my_box_plot_order, data=data, linewidth=frame_lw,
                          palette=colours, medianprops={"color": "r", "linewidth": 4}, width=0.4, flierprops=myflierprops, showfliers=False)
    
        ax.set_xlabel(x_axis, fontsize=28, labelpad=10)
        
    ax.set_ylabel(prop_name_latex, fontsize=32, labelpad=20)
               
    ax.tick_params(axis='x', which='major', top='on', bottom='on', pad=10, labelsize=22, length=12, width=frame_lw, direction='in', zorder=20) 
    
    ax.tick_params(axis='y', which='minor', top='on', bottom='on', pad=10, labelsize=28, length=8, width=frame_lw, direction='in', zorder=20) 
    ax.tick_params(axis='y', which='major', left='on', right='on', pad=10, labelsize=28, length=12, width=frame_lw, direction='in', zorder=20) 

    handles, _ = ax.get_legend_handles_labels()          # Get the artists.
    #ax.legend(handles, ["HSB", "LSB"], frameon=False, fancybox=False, loc=[0.15,1], fontsize=20, ncol=2)
    ax.legend(handles, ["voids", "walls", "skeleton"], frameon=False, fancybox=False, loc=[0.5,1], fontsize=20, ncol=nr_samples) 

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(frame_lw)
    
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(output_filename)
    plt.savefig(pp, format='pdf', rasterized=True, dpi=50, pad_inches=0.05, transparent=True)
    pp.close()
    
def count_outliers_boxplot(data):
    
    Q1 = np.nanpercentile(data, 25)
    Q3 = np.nanpercentile(data, 75)    
    IQR = Q3 - Q1
    
    lower = Q1 - 1.5 * IQR
    upper = Q3 + 1.5 * IQR
    #print "{0:.2f}".format(Q1), "{0:.2f}".format(Q3), "{0:.2f}".format(IQR)    

    data_outliers_lower = data[np.where(data < lower)[:][0]]
    data_outliers_upper = data[np.where(data > upper)[:][0]]
    
    data_outliers = np.append(data_outliers_lower, data_outliers_upper)
    #print data_outliers    

    return data_outliers, data_outliers.size , 100.0/data.size*data_outliers.size
    

    

    
    
    