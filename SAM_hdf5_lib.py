import os
system_info = os.getcwd()
start = system_info.find('anaconda')
mycomp = system_info[:start]

import numpy as np
import arangeData as aD
import outputData as oD

myOutput = oD.OutputData(config=False)

SAM_hdf5_filestruct_map = {}

band_map = {'0': 'u', '1': 'g', '2': 'r', '3': 'i', '4': 'z', '5': 'NUV', '6': 'FUV'}
map_surveys = {'0': 'SDSS_', '1': 'Galex_'}

map_comoving = {'0': 'observed'}
map_gal_parts = {'0': 'total_', '1': 'disk_', '2': 'spheroid_', '3': ''}

map_mags = {'0': 'mag', '1': 'mAB', '2': 'mAs'}

map_dust_atlas = {'0': 'dA_', '1': '', 'dA_': ':dustAtlas', '': '', '':''}

lum_map = {}
mag_type_map = {}

i=0
count_surveys=0
while i<len(map_surveys):
    j=0
    count=0
    while j<len(map_gal_parts):
        lum_map.update({str(count_surveys+j): 'L_'+map_surveys[str(i)]+'_'+map_gal_parts[str(j)]+'_'})
       
        a=0      
        while a<len(map_mags):
            mag_type_map.update({str(count): map_mags[str(a)]+'_'+map_gal_parts[str(j)]})
            count+=1              
            a+=1         
        j+=1
     
    count_surveys+=1
    i+=1


def catascii(catname,
             snapid,
             path_to_directory=False):

    if catname.find('nifty')!=-1:
        SAM_hdf5_filestruct_map[catname+'_haloid']  = 0
        SAM_hdf5_filestruct_map[catname+'_mcold']   = 9
        SAM_hdf5_filestruct_map[catname+'_mstar']   = 11
        
    elif catname.find('SDSS_DR7')!=-1:
        SAM_hdf5_filestruct_map[catname+'_DEC']     = 1
        SAM_hdf5_filestruct_map[catname+'_RA']     = 0
        SAM_hdf5_filestruct_map[catname+'_Z']       = 2
        SAM_hdf5_filestruct_map[catname+'_weight']       = 3
        SAM_hdf5_filestruct_map[catname+'_mstar']   = 4

         
    return SAM_hdf5_filestruct_map

def ROCKSTAR_HDF5_halocat_filestruct(catname,
                                     snapid,
                                     subdirid=False,
                                     path_to_directory=False):

    myData = aD.ArangeData()
    data = myData.readAnyFormat(config=False, 
                                mydtype=np.str_, 
                                mypath=mycomp+'anaconda/pro/data/'+catname+'/'+catname+'_file_names.txt', 
                                data_shape='shaped', 
                                comment='#')
      
    i=0
    while i<data.size:
        if data.size==1:
            SAM_hdf5_filestruct_map[catname+'_filename'+str(i)]            = path_to_directory+str(data)
        else:
            SAM_hdf5_filestruct_map[catname+'_filename'+str(i)]            = path_to_directory+data[i]

        i+=1

    SAM_hdf5_filestruct_map['catname']                              = catname
    SAM_hdf5_filestruct_map[catname+'_nr_files']                    = 1   
    SAM_hdf5_filestruct_map[catname+'_path_to_data']                = '/snapshot_'+snapid

    SAM_hdf5_filestruct_map[catname+'_haloid']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/HaloID'
    SAM_hdf5_filestruct_map[catname+'_npros']            = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/NumProgenitors'
    SAM_hdf5_filestruct_map[catname+'_firstProgenitorID']= SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/FirstProgenitor'   
    SAM_hdf5_filestruct_map[catname+'_siblingIndex']     = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Sibling'
    #print SAM_hdf5_filestruct_map 
    return SAM_hdf5_filestruct_map   

def Galacticus_HDF5_halocat_filestruct(catname,
                                       subdirid=False,
                                       path_to_directory=False):
    
    
    SAM_hdf5_filestruct_map['catname']                      = catname
    SAM_hdf5_filestruct_map[catname+'_path_to_directory']   = path_to_directory
    
    myData = aD.ArangeData()
    data = myData.readAnyFormat(config=False, mydtype=np.str_, mypath=path_to_directory+'input_file_names.txt', data_shape='shaped', comment='#')
      
    i=0
    while i<data.size:
        if data.size==1:
            SAM_hdf5_filestruct_map[catname+'_filename'+str(i)]            = path_to_directory+str(data)
        else:
            SAM_hdf5_filestruct_map[catname+'_filename'+str(i)]            = path_to_directory+data[i]
        i+=1
        
    SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']        = ''
        
    SAM_hdf5_filestruct_map[catname+'_halo_haloid']   = SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']+'/forestHalos/nodeIndex' 
    SAM_hdf5_filestruct_map[catname+'_halo_hostid']   = SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']+'/forestHalos/hostIndex'
    SAM_hdf5_filestruct_map[catname+'_halo_angmom']   = SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']+'/forestHalos/angularMomentum'
    SAM_hdf5_filestruct_map[catname+'_halo_desIndex']    = SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']+'/forestHalos/descendentIndex'         
    
    SAM_hdf5_filestruct_map[catname+'_halo_x_pos']    = SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']+'/forestHalos/position'
    SAM_hdf5_filestruct_map[catname+'_halo_y_pos']   = SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']+'/forestHalos/position'
    SAM_hdf5_filestruct_map[catname+'_halo_z_pos']    = SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']+'/forestHalos/position'
    
    SAM_hdf5_filestruct_map[catname+'_halo_nodemass']    = SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']+'/forestHalos/nodeMass'
    SAM_hdf5_filestruct_map[catname+'_halo_redshift']    = SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']+'/forestHalos/redshift'

    SAM_hdf5_filestruct_map[catname+'_halo_x_vel_disp']    = SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']+'/forestHalos/velocityDispersion'
    SAM_hdf5_filestruct_map[catname+'_halo_y_vel_disp']    = SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']+'/forestHalos/velocityDispersion'
    SAM_hdf5_filestruct_map[catname+'_halo_z_vel_disp']    = SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']+'/forestHalos/velocityDispersion'
    
    SAM_hdf5_filestruct_map[catname+'_halo_x_vel']    = SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']+'/forestHalos/velocity'
    SAM_hdf5_filestruct_map[catname+'_halo_y_vel']    = SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']+'/forestHalos/velocity'
    SAM_hdf5_filestruct_map[catname+'_halo_z_vel']    = SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']+'/forestHalos/velocity'

    SAM_hdf5_filestruct_map[catname+'_halo_rsca']    = SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']+'/forestHalos/scaleRadius'  
    SAM_hdf5_filestruct_map[catname+'_halo_vmax']    = SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']+'/forestHalos/velocityMaximum'
    SAM_hdf5_filestruct_map[catname+'_halo_spin']    = SAM_hdf5_filestruct_map[catname+'_path_to_halocat_data']+'/forestHalos/spin'
    
        
    return SAM_hdf5_filestruct_map    
 
def Galcticus_HDF5_filestruct(
                            catname,
                            snapid,
                            subdirid=False,
                            path_to_directory=False):  
    
    SAM_hdf5_filestruct_map['catname']                      = catname
    SAM_hdf5_filestruct_map[catname+'_snapid']              = snapid
    SAM_hdf5_filestruct_map[catname+'_path_to_directory']   = path_to_directory
    
    myData = aD.ArangeData()
    try:
        data = myData.readAnyFormat(config=False, mydtype=np.str_, mypath=path_to_directory+'file_names.txt', data_shape='shaped', comment='#')
    except:
        data = myData.readAnyFormat(config=False, mydtype=np.str_, mypath=mycomp+'anaconda/pro/data/'+catname+'/'+catname+'_file_names.txt', data_shape='shaped', comment='#')
                                   
    i=0
    while i<data.size:
        if data.size==1:
            SAM_hdf5_filestruct_map[catname+'_filename'+str(i)]            = path_to_directory+str(data)
        else:
            SAM_hdf5_filestruct_map[catname+'_filename'+str(i)]            = path_to_directory+data[i]
        i+=1      

    SAM_hdf5_filestruct_map[catname+'_nr_files'] = i    
   
    SAM_hdf5_filestruct_map[catname+'_path_to_data']        = 'Outputs/Output'
    SAM_hdf5_filestruct_map[catname+'_path_to_snapid']      = 'Outputs/Output'+snapid
    SAM_hdf5_filestruct_map[catname+'_redshift_attribute']  = 'outputExpansionFactor'
    
    SAM_hdf5_filestruct_map[catname+'_redshift']            = ''
    SAM_hdf5_filestruct_map[catname+'_Z']                   = ''

    SAM_hdf5_filestruct_map[catname+'_ngalaxies']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid
    SAM_hdf5_filestruct_map[catname+'_hostid']              = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/nodeIndex'
    SAM_hdf5_filestruct_map[catname+'_satelliteNodeIndex']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/satelliteNodeIndex' 
    SAM_hdf5_filestruct_map[catname+'_haloid']              = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/nodeIndex'
    SAM_hdf5_filestruct_map[catname+'_nodeIndex']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/nodeIndex'
    SAM_hdf5_filestruct_map[catname+'_parentIndex']         = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/parentIndex'
    SAM_hdf5_filestruct_map[catname+'_satelliteIndex']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/satelliteIndex'    
    SAM_hdf5_filestruct_map[catname+'_siblingIndex']        = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/siblingIndex'    
        
    SAM_hdf5_filestruct_map[catname+'_orphan']              = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/satelliteStatus'
    SAM_hdf5_filestruct_map[catname+'_satelliteMergeTime']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/satelliteMergeTime'     
    SAM_hdf5_filestruct_map[catname+'_isolated']            = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/nodeIsIsolated'
    SAM_hdf5_filestruct_map[catname+'_timeLastIsolated']    = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/basicTimeLastIsolated'
    
    SAM_hdf5_filestruct_map[catname+'_mhalo']               = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/basicMass'
    SAM_hdf5_filestruct_map[catname+'_mbasic']              = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/basicMass' 
    SAM_hdf5_filestruct_map[catname+'_mhalo_sat']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/satelliteBoundMass' 

    SAM_hdf5_filestruct_map[catname+'_mhalo_cents']               = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/basicMass'
    SAM_hdf5_filestruct_map[catname+'_mbasic_cents']              = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/basicMass' 
    
    SAM_hdf5_filestruct_map[catname+'_spinParameter']       = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/spinSpin'

    SAM_hdf5_filestruct_map[catname+'_mean_age_stars']     = ''
    SAM_hdf5_filestruct_map[catname+'_mbh']                = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/blackHoleMass'
    SAM_hdf5_filestruct_map[catname+'_bh_acc_rate']        = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/blackHoleAccretionRate'    
    SAM_hdf5_filestruct_map[catname+'_bh_rad_eff']         = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/blackHoleRadiativeEfficiency'
     
    SAM_hdf5_filestruct_map[catname+'_bhspin']              = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/blackHoleSpin'
    SAM_hdf5_filestruct_map[catname+'_bhcount']              = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/blackHoleCount'
    
    SAM_hdf5_filestruct_map[catname+'_angM_disk']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/diskAngularMomentum'
    SAM_hdf5_filestruct_map[catname+'_angM_spheroid']       = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/spheroidAngularMomentum'

    SAM_hdf5_filestruct_map[catname+'_rdisk']               = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/diskRadius' 
    SAM_hdf5_filestruct_map[catname+'_rbulge']              = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/spheroidRadius' 
    SAM_hdf5_filestruct_map[catname+'_rhalf_mass']          = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/halfMassRadius' 
    
    SAM_hdf5_filestruct_map[catname+'_mstar']               = ''
    SAM_hdf5_filestruct_map[catname+'_mstar_disk']          = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/diskMassStellar'
    SAM_hdf5_filestruct_map[catname+'_mstar_spheroid']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/spheroidMassStellar'     
    
    SAM_hdf5_filestruct_map[catname+'_mcold']               = ''
    SAM_hdf5_filestruct_map[catname+'_mcold_disk']          = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/diskMassGas'
    SAM_hdf5_filestruct_map[catname+'_mcold_spheroid']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/spheroidMassGas'

    SAM_hdf5_filestruct_map[catname+'_mhot']                = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/hotHaloMass'
    SAM_hdf5_filestruct_map[catname+'_mhot_outflow']        = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/hotHaloOutflowedMass'

    SAM_hdf5_filestruct_map[catname+'_sfr']                 = ''
    SAM_hdf5_filestruct_map[catname+'_sfr_disk']            = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/diskStarFormationRate'
    SAM_hdf5_filestruct_map[catname+'_sfr_spheroid']        = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/spheroidStarFormationRate'

    #What has the coluumn name of 'abundance' is acutally mass of metals --> wrongly named! [Msunh-1]
    SAM_hdf5_filestruct_map[catname+'_Mzstar']              = ''
    SAM_hdf5_filestruct_map[catname+'_Mzstar_disk']          = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/diskAbundancesStellarMetals'
    SAM_hdf5_filestruct_map[catname+'_Mzstar_spheroid']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/spheroidAbundancesStellarMetals'    
    
    SAM_hdf5_filestruct_map[catname+'_Mzgas']               = ''
    SAM_hdf5_filestruct_map[catname+'_Mzgas_disk']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/diskAbundancesGasMetals'
    SAM_hdf5_filestruct_map[catname+'_Mzgas_spheroid']       = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/spheroidAbundancesGasMetals'    

    SAM_hdf5_filestruct_map[catname+'_Mzhot_halo']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/hotHaloAbundancesMetals'   
   

    SAM_hdf5_filestruct_map[catname+'_x_pos']               = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/positionPositionX'
    SAM_hdf5_filestruct_map[catname+'_y_pos']               = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/positionPositionY'
    SAM_hdf5_filestruct_map[catname+'_z_pos']               = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/positionPositionZ'

    SAM_hdf5_filestruct_map[catname+'_x_vel']               = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/positionVelocityX'
    SAM_hdf5_filestruct_map[catname+'_y_vel']               = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/positionVelocityY'
    SAM_hdf5_filestruct_map[catname+'_z_vel']               = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/positionVelocityZ'

    SAM_hdf5_filestruct_map[catname+'_L_RGO_dA_spheroid_B_part1'] = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/spheroidLuminositiesStellar:RGO_B:rest:z'
    SAM_hdf5_filestruct_map[catname+'_L_RGO_dA_disk_B_part1']     = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/diskLuminositiesStellar:RGO_B:rest:z'
    SAM_hdf5_filestruct_map[catname+'_L_RGO_dA_spheroid_B']       = ''    
    SAM_hdf5_filestruct_map[catname+'_L_RGO_dA_spheroid_B_part2'] = ''
    SAM_hdf5_filestruct_map[catname+'_L_RGO_dA_disk_B_part2']     = ''
    SAM_hdf5_filestruct_map[catname+'_L_RGO_dA_disk_B']           = '' 
    SAM_hdf5_filestruct_map[catname+'_L_RGO_dA_total_B_part1']          = ''
    SAM_hdf5_filestruct_map[catname+'_L_RGO_dA_total_B_part2']          = ''
    SAM_hdf5_filestruct_map[catname+'_L_RGO_dA_total_B']          = ''

    SAM_hdf5_filestruct_map[catname+'_L_OII_disk_3726_part1']    = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/diskLineLuminosity:oxygenII3726:z'  
    SAM_hdf5_filestruct_map[catname+'_L_OII_spheroid_3726_part1']= SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/spheroidLineLuminosity:oxygenII3726:z' 
    SAM_hdf5_filestruct_map[catname+'_L_OII_disk_3726_part2']    = ''
    SAM_hdf5_filestruct_map[catname+'_L_OII_spheroid_3726_part2']= ''
    SAM_hdf5_filestruct_map[catname+'_L_OII_disk_3726']          = ''
    SAM_hdf5_filestruct_map[catname+'_L_OII_spheroid_3726']      = ''
    
    SAM_hdf5_filestruct_map[catname+'_L_OII_total_3726_part1']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/totalLineLuminosity:oxygenII3726:z'
    SAM_hdf5_filestruct_map[catname+'_L_OII_total_3726_part2']   = ''
    SAM_hdf5_filestruct_map[catname+'_L_OII_total_3726']         = ''  
    
    SAM_hdf5_filestruct_map[catname+'_L_OII_disk_3729_part1']    = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/diskLineLuminosity:oxygenII3729:z'  
    SAM_hdf5_filestruct_map[catname+'_L_OII_spheroid_3729_part1']= SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/spheroidLineLuminosity:oxygenII3729:z' 
    SAM_hdf5_filestruct_map[catname+'_L_OII_disk_3729_part2']    = ''
    SAM_hdf5_filestruct_map[catname+'_L_OII_spheroid_3729_part2']= ''
    SAM_hdf5_filestruct_map[catname+'_L_OII_disk_3729']          = ''
    SAM_hdf5_filestruct_map[catname+'_L_OII_spheroid_3729']      = ''
    SAM_hdf5_filestruct_map[catname+'_L_OII_total_3729']         = ''

    SAM_hdf5_filestruct_map[catname+'_L_OII_disk_cont_part1']    = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/diskOxygenContinuumLuminosity:z' 
    SAM_hdf5_filestruct_map[catname+'_L_OII_spheroid_cont_part1']= SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/spheroidOxygenContinuumLuminosity:z'  
    SAM_hdf5_filestruct_map[catname+'_L_OII_disk_cont_part2']    = ''
    SAM_hdf5_filestruct_map[catname+'_L_OII_spheroid_cont_part2']= ''
    SAM_hdf5_filestruct_map[catname+'_L_OII_disk_cont']          = ''
    SAM_hdf5_filestruct_map[catname+'_L_OII_spheroid_cont']      = ''
    SAM_hdf5_filestruct_map[catname+'_L_OII_total_cont']         = '' 

        
    SAM_hdf5_filestruct_map[catname+'_Mag_dA_total_B']    = ''
    
    SAM_hdf5_filestruct_map[catname+'_vmax']               = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/nBodyVelocityMaximum'   
    SAM_hdf5_filestruct_map[catname+'_vdisp']              = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/nBodyVelocityDispersion'  
    SAM_hdf5_filestruct_map[catname+'_sfr_disk_int']       = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/ageStatisticsDiskIntegratedSFR' 
    SAM_hdf5_filestruct_map[catname+'_sfr_spheroid_int']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/ageStatisticsSpheroidIntegratedSFR'
 
    SAM_hdf5_filestruct_map[catname+'_age_sfr_disk_int']       = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/ageStatisticsDiskTimeWeightedIntegratedSFR' 
    SAM_hdf5_filestruct_map[catname+'_age_sfr_spheroid_int']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/ageStatisticsSpheroidTimeWeightedIntegratedSFR'

  
    
    c=0
    while c<len(map_surveys):   
        i=0
        while i<len(band_map):
            a=0
            while a<len(lum_map):
                b=0
                while b<len(map_gal_parts):
                    d=0
                    while d<len(map_dust_atlas)/2.:
                        e=0
                        while e<len(map_comoving):
                            #print 'L_'+map_surveys[str(c)]+map_gal_parts[str(b)]+band_map[str(i)]
                            SAM_hdf5_filestruct_map[catname+'_L_'+map_surveys[str(c)]+map_dust_atlas[str(d)]+map_gal_parts[str(b)]+band_map[str(i)]] = ''
                            
                            if map_gal_parts[str(b)]=='':
                                
                                SAM_hdf5_filestruct_map[catname+'_L_'+map_surveys[str(c)]+map_dust_atlas[str(d)]+map_gal_parts[str(b)]+band_map[str(i)]+'_part1'] = ''
                                SAM_hdf5_filestruct_map[catname+'_L_'+map_surveys[str(c)]+map_dust_atlas[str(d)]+map_gal_parts[str(b)]+band_map[str(i)]+'_part2'] = ''
                                
                            else:
                                #print SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/'+map_gal_parts[str(b)][0:len(map_gal_parts[str(b)])-1]+'LuminositiesStellar:'+map_surveys[str(c)]+band_map[str(i)]+':'+map_comoving[str(e)]+':z',map_dust_atlas[str(d)] 
                                SAM_hdf5_filestruct_map[catname+'_L_'+map_surveys[str(c)]+map_dust_atlas[str(d)]+map_gal_parts[str(b)]+band_map[str(i)]+'_part1'] = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/nodeData/'+map_gal_parts[str(b)][0:len(map_gal_parts[str(b)])-1]+'LuminositiesStellar:'+map_surveys[str(c)]+band_map[str(i)]+':'+map_comoving[str(e)]+':z'                            
                                SAM_hdf5_filestruct_map[catname+'_L_'+map_surveys[str(c)]+map_dust_atlas[str(d)]+map_gal_parts[str(b)]+band_map[str(i)]+'_part2'] = map_dust_atlas[map_dust_atlas[str(d)]]           

                            
                            e+=1
                        d+=1
                    b+=1
                a+=1        
            
            j=0
            while j<len(mag_type_map):                    
                d=0
                while d<len(map_dust_atlas)/2.:    
                    SAM_hdf5_filestruct_map[catname+'_'+mag_type_map[str(j)]+map_dust_atlas[str(d)]+'_'+band_map[str(i)]] = ''
                    d+=1
                j+=1
            
            map_magnitudes_and_colour_cuts(catname, SAM_hdf5_filestruct_map, mag_type_map, index=i)
                
            i+=1
        c+=1
 
    
    #print SAM_hdf5_filestruct_map
        
    return SAM_hdf5_filestruct_map

def map_magnitudes_and_colour_cuts(catname,
                                   filestruct_map,
                                   map1,
                                   index=0):
    
    colour_cut_map = {'0': 'r_i', '1': 'd_mesa', '2': 'g_r'}
       
    a=0
    while a<len(colour_cut_map):            
        filestruct_map[catname+'_'+map1[str(index)]+'cut_'+colour_cut_map[str(a)]]  = ''
        a+=1

     
def SAG_HDF5_filestruct(catname,
                        snapid,
                        subdirid=1,
                        path_to_directory=False): 


#    i=0
#    while i<int(subdirid):
#        SAM_hdf5_filestruct_map[catname+'_filename'+str(i)]            = path_to_directory
#        i+=1  
    
    if len(snapid)<=2: snapid=str(3-len(snapid)-1)+snapid
    SAM_hdf5_filestruct_map['catname']                      = catname
    SAM_hdf5_filestruct_map[catname+'_snapid']              = 'snapshot_'+snapid
    SAM_hdf5_filestruct_map[catname+'_path_to_directory']   = path_to_directory
    
    #print 'path to directory:', path_to_directory
    
    myData = aD.ArangeData()
    data = myData.readAnyFormat(config=False, 
                                mydtype=np.str_, 
                                mypath=mycomp+'anaconda/pro/data/'+catname+'/'+catname+'_file_names.txt', 
                                data_shape='shaped', 
                                comment='#')
    
    #print 'data:', data.size, snapid

    i=0
    a=0
    while i<data.size:
        #print 'i:', i, 'a:', a, data[i].find(snapid), catname+'_filename'+str(a), data[i]
        if data.size==1:
            SAM_hdf5_filestruct_map[catname+'_filename'+str(a)] = path_to_directory+SAM_hdf5_filestruct_map[catname+'_snapid'] +'/'+str(data)
            a+=1
        else:
            #Check if the snapid (e.g. 125 or 075) corresponds to the filename "gal_+snapid_xxxxx.hdf5" --> be careful with the numerating of the files (001 to 128),
            #If the number of the file is corresponding by coindidence to the snapid the files wrongly chosen: Therefore two criteria >-1 and <=8 since the SAG files
            #are named like this e.g. gal_125_xxxx_numerating_xxx.hdf5
            if data[i].find(snapid)==4 or data[i].find(snapid)==5 or data[i].find(snapid)==8:  
                if catname.find('sub')!=-1:
                    #print 'HERE: sub!'
                    SAM_hdf5_filestruct_map[catname+'_filename'+str(a)] = path_to_directory+'/'+data[i]                
                else:
                    #print 'HERE: normal!'
                    SAM_hdf5_filestruct_map[catname+'_filename'+str(a)] = path_to_directory+SAM_hdf5_filestruct_map[catname+'_snapid'] +'/'+data[i]
                
                #print 'chosen filename:', SAM_hdf5_filestruct_map[catname+'_filename'+str(a)]
                a+=1
        i+=1  

    SAM_hdf5_filestruct_map[catname+'_nr_files'] = a 
  
    SAM_hdf5_filestruct_map[catname+'_path_to_data']        = ''
    SAM_hdf5_filestruct_map[catname+'_path_to_snapid']      = ''
    SAM_hdf5_filestruct_map[catname+'_redshift_attribute']  = 'Redshift'

    if catname.find('sub')!=-1:
   
        SAM_hdf5_filestruct_map[catname+'_haloid']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/HaloID' 
        SAM_hdf5_filestruct_map[catname+'_hostid']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/GalaxyStaticID'


    elif catname.find('v2')!=-1:

        SAM_hdf5_filestruct_map[catname+'_haloid']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/HaloID' 
        SAM_hdf5_filestruct_map[catname+'_hostid']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/GalaxyID' 

        SAM_hdf5_filestruct_map[catname+'_OH_gas_disk']    = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/OH_gas_disk'  
        SAM_hdf5_filestruct_map[catname+'_OH_gas_bulge']    = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/OH_gas_bulge'

        SAM_hdf5_filestruct_map[catname+'_MAB_total_u']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/Magnitudes/Mag_id119_AB_tot_r'
        SAM_hdf5_filestruct_map[catname+'_MAB_total_g']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/Magnitudes/Mag_id120_AB_tot_r'
        SAM_hdf5_filestruct_map[catname+'_MAB_total_r']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/Magnitudes/Mag_id121_AB_tot_r'
        SAM_hdf5_filestruct_map[catname+'_MAB_total_i']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/Magnitudes/Mag_id122_AB_tot_r'
        SAM_hdf5_filestruct_map[catname+'_MAB_total_z']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/Magnitudes/Mag_id123_AB_tot_r'
             
        SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_u']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/Magnitudes/Mag_ext_id119_AB_tot_r'
        SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_g']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/Magnitudes/Mag_ext_id120_AB_tot_r'
        SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_r']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/Magnitudes/Mag_ext_id121_AB_tot_r'
        SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_i']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/Magnitudes/Mag_ext_id122_AB_tot_r'
        SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_z']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/Magnitudes/Mag_ext_id123_AB_tot_r'        

        SAM_hdf5_filestruct_map[catname+'_OII_3727_ext']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/EmissionLines/OII_3727_ext'
        SAM_hdf5_filestruct_map[catname+'_OII_3727']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/EmissionLines/OII_3727'
        SAM_hdf5_filestruct_map[catname+'_OII_3729_ext']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/EmissionLines/OII_3729_ext'
        SAM_hdf5_filestruct_map[catname+'_OII_3729']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/EmissionLines/OII_3729'    

        SAM_hdf5_filestruct_map[catname+'_OII_cont_3727_ext']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/EmissionLines/Cont_OII_3727_ext'
        SAM_hdf5_filestruct_map[catname+'_OII_cont_3727']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/EmissionLines/Cont_OII_3727'
        SAM_hdf5_filestruct_map[catname+'_OII_cont_3729_ext']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/EmissionLines/Cont_OII_3729_ext'
        SAM_hdf5_filestruct_map[catname+'_OII_cont_3729']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/EmissionLines/Cont_OII_3729'
        
    else:
    
        SAM_hdf5_filestruct_map[catname+'_ngalaxies']= SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/'    
        SAM_hdf5_filestruct_map[catname+'_haloid']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/MainHaloID' 
        SAM_hdf5_filestruct_map[catname+'_hostid']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/GalaxyHostID'
       
        SAM_hdf5_filestruct_map[catname+'_OH_gas_disk_bulge']    = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/OH_gas_disk_bulge' 

        SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_u']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Magnitudes/Mag_uS_dust'
        SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_g']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Magnitudes/Mag_gS_dust'
        SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_r']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Magnitudes/Mag_rS_dust'
        SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_i']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Magnitudes/Mag_iS_dust'
        SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_z']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Magnitudes/Mag_zS_dust'       
    
        SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_B']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Magnitudes/Mag_B_dust'      
    
        SAM_hdf5_filestruct_map[catname+'_MAB_total_u']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Magnitudes/Mag_uS'
        SAM_hdf5_filestruct_map[catname+'_MAB_total_g']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Magnitudes/Mag_gS'
        SAM_hdf5_filestruct_map[catname+'_MAB_total_r']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Magnitudes/Mag_rS'
        SAM_hdf5_filestruct_map[catname+'_MAB_total_i']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Magnitudes/Mag_iS'
        SAM_hdf5_filestruct_map[catname+'_MAB_total_z']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Magnitudes/Mag_zS'


    
    
    SAM_hdf5_filestruct_map[catname+'_ngalaxies']= SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/' 

    SAM_hdf5_filestruct_map[catname+'_orphan']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Galaxy_Type'
    SAM_hdf5_filestruct_map[catname+'_mhalo']    = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Halo/M200c'
    SAM_hdf5_filestruct_map[catname+'_mhalo_cents']    = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Halo/M200c'
    SAM_hdf5_filestruct_map[catname+'_r200c']    = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Halo/R200c'
    SAM_hdf5_filestruct_map[catname+'_vmax']     = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Halo/Vmax'
    SAM_hdf5_filestruct_map[catname+'_vpeak']    = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Halo/Vpeak'
    SAM_hdf5_filestruct_map[catname+'_NFW_con']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Halo/cNFW'

    SAM_hdf5_filestruct_map[catname+'_spinParameter']   = ''

    SAM_hdf5_filestruct_map[catname+'_mbh']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Mbh'       
    SAM_hdf5_filestruct_map[catname+'_bhspin']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Lambda'
    
    SAM_hdf5_filestruct_map[catname+'_mean_age_stars']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/T_stars' 

    SAM_hdf5_filestruct_map[catname+'_rhalf_bulge']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Rhalf_bulge'
    SAM_hdf5_filestruct_map[catname+'_rhalf_disk']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Rhalf_disk'        
 
    SAM_hdf5_filestruct_map[catname+'_mstar'] = ''
    SAM_hdf5_filestruct_map[catname+'_mstar_disk'] = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/M_star_disk'
    SAM_hdf5_filestruct_map[catname+'_mstar_spheroid'] = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/M_star_bulge'

    SAM_hdf5_filestruct_map[catname+'_mstar_IC']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/M_ICstars' 
   
    SAM_hdf5_filestruct_map[catname+'_mcold']    = ''
    SAM_hdf5_filestruct_map[catname+'_mcold_disk']    = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/M_gas_disk'
    SAM_hdf5_filestruct_map[catname+'_mcold_spheroid']    = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/M_gas_bulge'

    SAM_hdf5_filestruct_map[catname+'_mbh']     = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Mbh'  
    SAM_hdf5_filestruct_map[catname+'_mhot']    = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/M_hot'

    
    SAM_hdf5_filestruct_map[catname+'_sfr']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SFR'
    SAM_hdf5_filestruct_map[catname+'_sfr_disk']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SFR'
    SAM_hdf5_filestruct_map[catname+'_sfr_spheroid']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SFR_bulge'

    SAM_hdf5_filestruct_map[catname+'_sfr_spheroid_inst']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/InstantSfrBulge'
    SAM_hdf5_filestruct_map[catname+'_sfr_quies_inst']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/InstantSfrQuies'

    SAM_hdf5_filestruct_map[catname+'_zgas_disk']     = ''
    SAM_hdf5_filestruct_map[catname+'_zgas_spheroid']     = ''
    
    SAM_hdf5_filestruct_map[catname+'_zstar_disk']     = ''
    SAM_hdf5_filestruct_map[catname+'_zstar_spheroid']     = ''

    SAM_hdf5_filestruct_map[catname+'_zhot_halo']     = ''


    SAM_hdf5_filestruct_map[catname+'_Mzstar']     = ''
    SAM_hdf5_filestruct_map[catname+'_Mzstar_disk']     = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/MZ_stars_disk'
    SAM_hdf5_filestruct_map[catname+'_Mzstar_spheroid']     = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/MZ_stars_bulge'

    SAM_hdf5_filestruct_map[catname+'_Mzgas']     = ''
    SAM_hdf5_filestruct_map[catname+'_Mzgas_disk']     = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/MZ_gas_disk'
    SAM_hdf5_filestruct_map[catname+'_Mzgas_spheroid']     = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/MZ_gas_bulge'

    SAM_hdf5_filestruct_map[catname+'_Mzhot_halo']     = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/MZ_hot_halo'

 

    SAM_hdf5_filestruct_map[catname+'_x_pos']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/X'
    SAM_hdf5_filestruct_map[catname+'_y_pos']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Y'
    SAM_hdf5_filestruct_map[catname+'_z_pos']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Z'

    SAM_hdf5_filestruct_map[catname+'_x_vel']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Vx'
    SAM_hdf5_filestruct_map[catname+'_y_vel']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Vy'
    SAM_hdf5_filestruct_map[catname+'_z_vel']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Vz'
    
#    i=0
#    while i<len(band_map):
#        j=0
#        while j<len(mag_type_map):    
#            SAM_hdf5_filestruct_map[catname+'_'+mag_type_map[str(i)]+band_map[str(i)]] = SAM_hdf5_filestruct_map[catname+'_path_to_data']+snapid+'/Magnitudes/Mag_ext_id122_AB_tot_'+band_map[str(i)]
#            
#            a=0
#            while a<len(lum_map):           
#                SAM_hdf5_filestruct_map[catname+'_'+lum_map[str(a)]+band_map[str(i)]]   = ''            
#                a+=1
#            
#            j+=1
#
#        map_magnitudes_and_colour_cuts(catname, SAM_hdf5_filestruct_map, mag_type_map, index=i)
#         
#        i+=1 

    SAM_hdf5_filestruct_map[catname+'_mAB_total_u']   = ''
    SAM_hdf5_filestruct_map[catname+'_mAB_total_g']   = ''
    SAM_hdf5_filestruct_map[catname+'_mAB_total_r']   = ''
    SAM_hdf5_filestruct_map[catname+'_mAB_total_i']   = ''
    SAM_hdf5_filestruct_map[catname+'_mAB_total_z']   = ''

    SAM_hdf5_filestruct_map[catname+'_mAB_dA_total_u']   = ''
    SAM_hdf5_filestruct_map[catname+'_mAB_dA_total_g']   = ''
    SAM_hdf5_filestruct_map[catname+'_mAB_dA_total_r']   = ''
    SAM_hdf5_filestruct_map[catname+'_mAB_dA_total_i']   = ''
    SAM_hdf5_filestruct_map[catname+'_mAB_dA_total_z']   = ''    
        
    #print 'SAG filestructure:',  SAM_hdf5_filestruct_map     
    
    return SAM_hdf5_filestruct_map
    
def SAG_NIFTY_filestruct(catname,
                        snapid,
                        subdirid=1,
                        path_to_directory=False): 

#    i=0
#    while i<int(subdirid):
#        SAM_hdf5_filestruct_map[catname+'_filename'+str(i)]            = path_to_directory
#        i+=1  
    
    if len(snapid)<=2: snapid=str(3-len(snapid)-1)+snapid
    SAM_hdf5_filestruct_map['catname']                      = catname
    SAM_hdf5_filestruct_map[catname+'_snapid']              = 'snapshot_'+snapid
#    SAM_hdf5_filestruct_map[catname+'_path_to_directory']   = path_to_directory
#    
#    myData = aD.ArangeData()
#    data = myData.readAnyFormat(config=False, 
#                                mydtype=np.str_, 
#                                mypath=mycomp+'anaconda/pro/data/'+catname+'/'+catname+'_file_names.txt', 
#                                data_shape='shaped', 
#                                comment='#')
    
    #print 'data:', data.size, snapid
       
#    i=0
#    a=0
#    while i<data.size:
#        #print 'i:', i, 'a:', a, data[i].find(snapid), catname+'_filename'+str(a), data[i]
#        if data.size==1:
#            SAM_hdf5_filestruct_map[catname+'_filename'+str(a)] = path_to_directory+SAM_hdf5_filestruct_map[catname+'_snapid'] +'/'+str(data)
#            a+=1
#        else:
#            #Check if the snapid (e.g. 125 or 075) corresponds to the filename "gal_+snapid_xxxxx.hdf5" --> be careful with the numerating of the files (001 to 128),
#            #If the number of the file is corresponding by coindidence to the snapid the files wrongly chosen: Therefore two criteria >-1 and <=8 since the SAG files
#            #are named like this e.g. gal_125_xxxx_numerating_xxx.hdf5
#            if data[i].find(snapid)==4 or data[i].find(snapid)==5 :               
#                SAM_hdf5_filestruct_map[catname+'_filename'+str(a)] = path_to_directory+SAM_hdf5_filestruct_map[catname+'_snapid'] +'/'+data[i]
#                #print 'chosen filename:', SAM_hdf5_filestruct_map[catname+'_filename'+str(a)]
#                a+=1
#        i+=1  

    #SAM_hdf5_filestruct_map[catname+'_nr_files'] = a 
  
    SAM_hdf5_filestruct_map[catname+'_path_to_data']        = ''
    SAM_hdf5_filestruct_map[catname+'_path_to_snapid']      = ''
    SAM_hdf5_filestruct_map[catname+'_redshift_attribute']  = ''
    
    SAM_hdf5_filestruct_map[catname+'_haloid']   = 0 
    SAM_hdf5_filestruct_map[catname+'_hostid']   = 1
    SAM_hdf5_filestruct_map[catname+'_orphan']   = 2
    SAM_hdf5_filestruct_map[catname+'_mhalo']    = ''

    
    SAM_hdf5_filestruct_map[catname+'_mstar'] = 11   
    SAM_hdf5_filestruct_map[catname+'_mbh']      = 12  
    SAM_hdf5_filestruct_map[catname+'_mcold']    = 9

     
    SAM_hdf5_filestruct_map[catname+'_sfr']      = 16
    #SAG has no 'sfr_disk' coluumn! In order to make the programm smooth the 'sfr' coluumn will be used instead
    SAM_hdf5_filestruct_map[catname+'_sfr_spheroid']      = 17

    SAM_hdf5_filestruct_map[catname+'_zgas']     = ''
    SAM_hdf5_filestruct_map[catname+'_zgas_disk']     = 33
    SAM_hdf5_filestruct_map[catname+'_zgas_spheroid']     = 34
    
    SAM_hdf5_filestruct_map[catname+'_x_pos']   = 3
    SAM_hdf5_filestruct_map[catname+'_y_pos']   = 4
    SAM_hdf5_filestruct_map[catname+'_z_pos']   = 5

    SAM_hdf5_filestruct_map[catname+'_x_vel']   = 6
    SAM_hdf5_filestruct_map[catname+'_y_vel']   = 7
    SAM_hdf5_filestruct_map[catname+'_z_vel']   = 8
    

    SAM_hdf5_filestruct_map[catname+'_mAB_total_u']   = ''
    SAM_hdf5_filestruct_map[catname+'_mAB_total_g']   = ''
    SAM_hdf5_filestruct_map[catname+'_mAB_total_r']   = ''
    SAM_hdf5_filestruct_map[catname+'_mAB_total_i']   = ''
    SAM_hdf5_filestruct_map[catname+'_mAB_total_z']   = ''

    SAM_hdf5_filestruct_map[catname+'_mAB_dA_total_u']   = ''
    SAM_hdf5_filestruct_map[catname+'_mAB_dA_total_g']   = ''
    SAM_hdf5_filestruct_map[catname+'_mAB_dA_total_r']   = ''
    SAM_hdf5_filestruct_map[catname+'_mAB_dA_total_i']   = ''
    SAM_hdf5_filestruct_map[catname+'_mAB_dA_total_z']   = ''    
    
    SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_u']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/Magnitudes/Mag_ext_id119_AB_tot_r'
    SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_g']   = 46
    SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_r']   = 47
    SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_i']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/Magnitudes/Mag_ext_id122_AB_tot_r'
    SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_z']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/Magnitudes/Mag_ext_id123_AB_tot_r'

    SAM_hdf5_filestruct_map[catname+'_MAB_total_u']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/Magnitudes/Mag_id119_AB_tot_r'
    SAM_hdf5_filestruct_map[catname+'_MAB_total_g']   = 52
    SAM_hdf5_filestruct_map[catname+'_MAB_total_r']   = 53
    SAM_hdf5_filestruct_map[catname+'_MAB_total_i']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/Magnitudes/Mag_id122_AB_tot_r'
    SAM_hdf5_filestruct_map[catname+'_MAB_total_z']   = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SED/Magnitudes/Mag_id123_AB_tot_r'
  
    #print 'SAG filestructure:',  SAM_hdf5_filestruct_map     
    
    return SAM_hdf5_filestruct_map
    
def SAGE_HDF5_filestruct(catname,
                        snapid,
                        subdirid=1,
                        path_to_directory=False): 
  
    SAM_hdf5_filestruct_map['catname']  = catname
    
    myData = aD.ArangeData()
    data = myData.readAnyFormat(config=False, 
                                mydtype=np.str_, 
                                mypath=mycomp+'anaconda/pro/data/'+catname+'/'+catname+'_file_names.txt', 
                                data_shape='shaped', 
                                comment='#')
                        
    SAM_hdf5_filestruct_map[catname+'_filename0']=str(data)     
    
    SAM_hdf5_filestruct_map[catname+'_nr_files'] = 1
    SAM_hdf5_filestruct_map[catname+'_path_to_data']=''
    SAM_hdf5_filestruct_map[catname+'_snapid']              = ''
    SAM_hdf5_filestruct_map[catname+'_redshift_attribute']  = 'Redshift'

    SAM_hdf5_filestruct_map[catname+'_ngalaxies']       = SAM_hdf5_filestruct_map[catname+'_path_to_data']+''    
    SAM_hdf5_filestruct_map[catname+'_haloid']          = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Central_Galaxy_ID' 
    SAM_hdf5_filestruct_map[catname+'_hostid']          = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Galaxy_ID'
    SAM_hdf5_filestruct_map[catname+'_orphan']          = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Galaxy_Classification'
    SAM_hdf5_filestruct_map[catname+'_mhalo']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Mvir'
    SAM_hdf5_filestruct_map[catname+'_mhalo_cents']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Central_Galaxy_Mvir'
    SAM_hdf5_filestruct_map[catname+'_mhalo_sat']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Subhalo_Mvir_at_Infall'
    SAM_hdf5_filestruct_map[catname+'_rvir']            = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Rvir'
    SAM_hdf5_filestruct_map[catname+'_vvir']            = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Vvir'
    SAM_hdf5_filestruct_map[catname+'_vvir_sat']            = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Subhalo_Vvir_at_Infall'
    SAM_hdf5_filestruct_map[catname+'_vmax']            = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Maximum_circular_velocity_of_the_halo'
    SAM_hdf5_filestruct_map[catname+'_vmax_sat']            = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Subhalo_Vmax_at_Infall' 
    SAM_hdf5_filestruct_map[catname+'_vdisp']            = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Velocity_Dispersion'
 
    SAM_hdf5_filestruct_map[catname+'_rdisk']            = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Disk_Scale_Radius'    
    
    SAM_hdf5_filestruct_map[catname+'_mstar']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Total_Stellar_Mass'
    SAM_hdf5_filestruct_map[catname+'_mstar_disk']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Total_Stellar_Mass'
    SAM_hdf5_filestruct_map[catname+'_mstar_spheroid']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Bulge_Stellar_Mass'
    
    SAM_hdf5_filestruct_map[catname+'_mstar_IC']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Intracluster_Stars_Mass'

    SAM_hdf5_filestruct_map[catname+'_mbh']             = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Black_Hole_Mass'  
    SAM_hdf5_filestruct_map[catname+'_mcold_disk']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Cold_Gas_Mass'
    SAM_hdf5_filestruct_map[catname+'_mhot']            = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Hot_Gas_Mass'
    
    SAM_hdf5_filestruct_map[catname+'_sfr']             = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Total_Star_Formation_Rate'
    SAM_hdf5_filestruct_map[catname+'_mean_age_stars']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Mean_Age_of_Stars'

    SAM_hdf5_filestruct_map[catname+'_Mzstar']          = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Metals_Total_Stellar_Mass'
    SAM_hdf5_filestruct_map[catname+'_Mzstar_disk']          = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Metals_Total_Stellar_Mass'
    SAM_hdf5_filestruct_map[catname+'_Mzstar_spheroid'] = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Metals_Bulge_Mass'
    SAM_hdf5_filestruct_map[catname+'_Mzgas_disk']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Metals_Cold_Gas_Mass'
    SAM_hdf5_filestruct_map[catname+'_Mzhot_halo']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Metals_Hot_Gas_Mass'

    SAM_hdf5_filestruct_map[catname+'_Z']            = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Redshift_Observed'

    SAM_hdf5_filestruct_map[catname+'_DEC']            = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Declination'
    SAM_hdf5_filestruct_map[catname+'_RA']            = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Right_Ascension'
   
    SAM_hdf5_filestruct_map[catname+'_x_pos']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/X'
    SAM_hdf5_filestruct_map[catname+'_y_pos']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Y'
    SAM_hdf5_filestruct_map[catname+'_z_pos']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Z'

    SAM_hdf5_filestruct_map[catname+'_x_vel']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/X_Velocity'
    SAM_hdf5_filestruct_map[catname+'_y_vel']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Y_Velocity'
    SAM_hdf5_filestruct_map[catname+'_z_vel']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Z_Velocity'

    SAM_hdf5_filestruct_map[catname+'_mAB_dA_total_u']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SDSS_u_Apparent'
    SAM_hdf5_filestruct_map[catname+'_mAB_dA_total_g']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SDSS_g_Apparent'
    SAM_hdf5_filestruct_map[catname+'_mAB_dA_total_r']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SDSS_r_Apparent'
    SAM_hdf5_filestruct_map[catname+'_mAB_dA_total_i']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SDSS_i_Apparent'
    SAM_hdf5_filestruct_map[catname+'_mAB_dA_total_z']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SDSS_z_Apparent'    
    
    SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_u']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SDSS_u_Absolute'
    SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_g']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SDSS_g_Absolute'
    SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_r']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SDSS_r_Absolute'
    SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_i']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SDSS_i_Absolute'
    SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_z']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/SDSS_z_Absolute'
  
    #print 'SAGE filestructure:',  SAM_hdf5_filestruct_map     
    
    return SAM_hdf5_filestruct_map


def LGALAXIES_HDF5_Tree_filestruct(catname,
                            snapid,
                            subdirid=1,
                            path_to_directory=False): 

  
    SAM_hdf5_filestruct_map['catname']  = catname
    
    myData = aD.ArangeData()
    data = myData.readAnyFormat(config=False, 
                                mydtype=np.str_, 
                                mypath=mycomp+'anaconda/pro/data/'+catname+'/'+catname+'_file_names.txt', 
                                data_shape='shaped', 
                                comment='#')
                                    
    print data                                
                                
    SAM_hdf5_filestruct_map[catname+'_nr_files']            = 1
    SAM_hdf5_filestruct_map[catname+'_path_to_data']        = '/Arr'
    SAM_hdf5_filestruct_map[catname+'_index']               = '/Index'
    SAM_hdf5_filestruct_map[catname+'_path_to_scale_factor']= '/Snapshot/a'
    SAM_hdf5_filestruct_map[catname+'_path_to_look_back_time']= '/Snapshot/LBT'
    SAM_hdf5_filestruct_map[catname+'_path_to_hubble_const']= '/Snapshot/Hz'
    
    SAM_hdf5_filestruct_map[catname+'_filename0']           = path_to_directory+str(data)
    SAM_hdf5_filestruct_map[catname+'_snapid']              = ''
    SAM_hdf5_filestruct_map[catname+'_redshift_attribute']  = ''

    SAM_hdf5_filestruct_map[catname+'_halo_age']            = '/Central/Halo_Age'
    SAM_hdf5_filestruct_map[catname+'_haloid']              = ''   
    SAM_hdf5_filestruct_map[catname+'_hostid']              = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/GalaxyID'
    SAM_hdf5_filestruct_map[catname+'_orphan']              = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Type'
    SAM_hdf5_filestruct_map[catname+'_mhalo']               = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Mvir'
    SAM_hdf5_filestruct_map[catname+'_mstar']               = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/MStell'
    SAM_hdf5_filestruct_map[catname+'_sfr']                 = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/sfr'
    SAM_hdf5_filestruct_map[catname+'_Z']                   = ''     
    #rint 'LGALAXIES filestructure:',  SAM_hdf5_filestruct_map     
    
    return SAM_hdf5_filestruct_map 

def LGALAXIES_HDF5_filestruct(catname,
                            snapid,
                            subdirid=1,
                            path_to_directory=False): 

  
    SAM_hdf5_filestruct_map['catname']  = catname
    
    myData = aD.ArangeData()
    data = myData.readAnyFormat(config=False, 
                                mydtype=np.str_, 
                                mypath=mycomp+'anaconda/pro/data/'+catname+'/'+catname+'_file_names.txt', 
                                data_shape='shaped', 
                                comment='#')
                                    
    SAM_hdf5_filestruct_map[catname+'_nr_files']            = 1
    SAM_hdf5_filestruct_map[catname+'_path_to_data']        = '/Galaxy'
    SAM_hdf5_filestruct_map[catname+'_filename0']           = path_to_directory+str(data)
    SAM_hdf5_filestruct_map[catname+'_snapid']              = ''
    SAM_hdf5_filestruct_map[catname+'_redshift_attribute']  = ''

    #print 'path to data:', SAM_hdf5_filestruct_map[catname+'_path_to_data']



    SAM_hdf5_filestruct_map[catname+'_ngalaxies']       = SAM_hdf5_filestruct_map[catname+'_path_to_data']+''    
    SAM_hdf5_filestruct_map[catname+'_haloid']          = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/haloID' 
    SAM_hdf5_filestruct_map[catname+'_hostid']          = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/GalaxyID'
    SAM_hdf5_filestruct_map[catname+'_lastProID']       = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/lastProgenitorID'
    SAM_hdf5_filestruct_map[catname+'_fofID']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/fofID'
    SAM_hdf5_filestruct_map[catname+'_mainLeafID']      = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/mainLeafID'
    SAM_hdf5_filestruct_map[catname+'_orphan']          = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Type'
    SAM_hdf5_filestruct_map[catname+'_mhalo']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Mvir'
    SAM_hdf5_filestruct_map[catname+'_mhalo_cents']     = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Mvir'
    SAM_hdf5_filestruct_map[catname+'_rvir']            = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Rvir'
    SAM_hdf5_filestruct_map[catname+'_vvir']            = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Vvir'
    SAM_hdf5_filestruct_map[catname+'_vmax']            = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Vmax'
    SAM_hdf5_filestruct_map[catname+'_vmax_sat']        = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Vmax_infall' 

    SAM_hdf5_filestruct_map[catname+'_rdisk']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/stellarDiskRadius'  
    SAM_hdf5_filestruct_map[catname+'_rbulge']          = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/bulgeSize'    
    
    SAM_hdf5_filestruct_map[catname+'_mstar']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/MStell'
    SAM_hdf5_filestruct_map[catname+'_mcold']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/CGmass'
    
    SAM_hdf5_filestruct_map[catname+'_sfr']             = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/sfr'
    SAM_hdf5_filestruct_map[catname+'_ssfr']             = ''
    SAM_hdf5_filestruct_map[catname+'_zcold']             = ''
    SAM_hdf5_filestruct_map[catname+'_cgf']             = ''
   
    SAM_hdf5_filestruct_map[catname+'_x_pos']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Pos'
    SAM_hdf5_filestruct_map[catname+'_y_pos']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Pos'
    SAM_hdf5_filestruct_map[catname+'_z_pos']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Pos'

    SAM_hdf5_filestruct_map[catname+'_x_vel']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Vel'
    SAM_hdf5_filestruct_map[catname+'_y_vel']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Vel'
    SAM_hdf5_filestruct_map[catname+'_z_vel']           = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/Vel'   
    
    SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_u']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/uDust'
    SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_g']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/gDust'
    SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_r']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/rDust'
    SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_i']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/iDust'
    SAM_hdf5_filestruct_map[catname+'_MAB_dA_total_z']  = SAM_hdf5_filestruct_map[catname+'_path_to_data']+'/zDust'
      
    #print 'LGALAXIES filestructure:',  SAM_hdf5_filestruct_map     
    
    return SAM_hdf5_filestruct_map 