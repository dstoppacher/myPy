"""
@author: DS (doris.stoppacher@csic.es)                    last update 12/26/2017

Descriptions:
-----------------
This script reads user selected galaxy properties through a list OR optional all
galaxy properties column wise to a numpy structured array. Galaxy property names,
units and data type will automatically mapped.

The data type and units of the properties will be taken automatically from the
HDF5-file. You only need to enter your favorite galaxy properties in the list
'mycol' using the corresponding IDs from the dictionary 'gal_properties'.

When the reading was successfully - you will get a message printed on the screen -
then your data structure is ready to work with. Below a few examples how to deal
with it.

Input:
-----------------
List of integer numbers referring to the galaxy property names OR empty list (then all
properties will be read into a structure)

Output:
-----------------
numpy structured array with as many columns as galaxy properties chosen in the
input list, units and data types are automatically included.

name
-----------------

gal_properties: python dictionary {}, stores the names of all available galaxy
                properties given in the MultiDark-Galaxies HDF5 output file.
                Properties are alphabetically order. Please NOTE that not
                all semi-analytical models have all properties stored in their
                output files.
                   
mycol:          python list [], you enter the integer values of the names in
                the 'gal_properties' dictionary.
                example:  mycol=[1,2,5,18] for properties with 'name1', name2, 
                ...etc.
                In this example the following properties will be read into the structure
                    --> HaloMass, HostHaloID, LstarSDSSr and MstarDisk.
               
                if you let the 'mycol' list empty like:
                    mycol=[]
                then all available properties for a certain model we be read into
                the structure.
                   
                Information about the galaxy property and what the name stands for
                can be found in the
                      README_MultiDark-Galaxies.txt
                where the file is stored in the same folder as this script.
               
Examples:
-----------------

I want to address the 'HaloMass'
++++++++++++++++++++++++++++++++++++++++++++++++++++

    print data_struct['HaloMass']
   
    --> only write the name of the property in [] next to the data structure.

the unit of the 'HaloMass' can be found here and printed to screen

    print data_struct['HaloMass'].attrs['unit'] --> h-1Msun
   

I only want to print the first 10 'HaloMasses'
++++++++++++++++++++++++++++++++++++++++++++++++++++

    print data_struct['HaloMass'][0:9]
   
   
I want info about the column 'HaloMass'
++++++++++++++++++++++++++++++++++++++++++++++++++++

    print data_struct.dtype.fields['HaloMass']

All fields of the galaxy property 'HaloMass' is printed on the screen
        --> (dtype('float64'), 1, 'h-1Msun')

The output has always the same order: data type, ID number from the 'gal_properties'
dictionary, units of the properties

if you want all properties of your 'data_struct' to be printed use:
 
    for k in data_struct.dtype.names:
        print 'galaxy property:', k, 'unit:', data_struct.dtype.fields[k][2], \
              'data_type:', data_struct.dtype.fields[k][0], 'ID number:', \
               data_struct.dtype.fields[k][1]    

Get information about the MD-Galaxies output file of your choice
++++++++++++++++++++++++++++++++++++++++++++++++++++

    output_file_info = {}                          
    i=0
    while i<len(f.attrs.keys()):
        #in the keys of the file attributes information about name of the
        #semi-analytical model, the dark matter simulation, box size etc. can be read-in
        print f.attrs.keys()[i], f.attrs.values()[i]
        output_file_info[f.attrs.keys()[i]] = f.attrs.values()[i]
        i+=1
   
    #you can then set those attributes to proper name to use further in your
    #script e.g. the redshift of the snapshot
    redshift     = output_file_info['redshift']
    print 'redshift of this catalog is,' redshift
  
"""

import numpy as np
import h5py as hdf5

def readMDGal():        
    print '#####################################################################'
    print 'read MultiDark-Galaxies HDF5 file format\n'
    
    #Dictionary in alphabethical order of the galaxy properties available!
    #Please note that not every model has all properties available ...
    gal_properties={
                     'name0': 'GalaxyType',
                     'name1': 'HaloMass',
                     'name2': 'HostHaloID',
                     'name3': 'LstarSDSSg',
                     'name4': 'LstarSDSSi',
                     'name5': 'LstarSDSSr',
                     'name6': 'LstarSDSSu',
                     'name7': 'LstarSDSSz',                         
                     'name8': 'MagStarSDSSg',
                     'name9': 'MagStarSDSSi',
                     'name10': 'MagStarSDSSr',
                     'name11': 'MagStarSDSSu',
                     'name12': 'MagStarSDSSz',
                     'name13': 'MainHaloID',
                     'name14': 'Mbh',
                     'name15': 'McoldSpheroid',
                     'name16': 'MeanAgeStars',
                     'name17': 'Mhot',
                     'name18': 'MstarDisk',
                     'name19': 'MstarIC',
                     'name20': 'MstarSpheroid',
                     'name21': 'MZgasDisk',
                     'name22': 'MZgasSpheroid',
                     'name23': 'MZhotHalo',
                     'name24': 'MZstarDisk',
                     'name25': 'MZstarSpheroid',
                     'name26': 'NFWconcentration',
                     'name27': 'OH_gas_disk_bulge',
                     'name28': 'rbulge',
                     'name29': 'rdisk',
                     'name30': 'rhalf_bulge',
                     'name31': 'rhalf_disk',
                     'name32': 'rhalf_mass',
                     'name33': 'sfr_quies_inst',
                     'name34': 'sfr_spheroid_inst',                         
                     'name35': 'SFRdisk',
                     'name36': 'SFRspheroid',
                     'name37': 'SpinParameter',
                     'name38': 'Vmax',
                     'name39': 'Vpeak',
                     'name40': 'Vx',
                     'name41': 'Vy',
                     'name42': 'Vz',
                     'name43': 'X',
                     'name44': 'Y',
                     'name45': 'Z',
                     'name46': 'ZgasDisk', 
                     'name47': 'ZgasSpheroid'                 
            }

    #dfefault: all galaxy properties for a certain model will be read in!    
    mycol=[]
    #Enter the numbers of the galaxy properties you want to read in! -->
    #mycol=[0, 1, 2, 3, 4]
    
    if mycol==[]:
        mycol=[k for k in np.arange(len(gal_properties))]
    
    #enter filename you want to read in:
    common_path ='/store/erebos/doris/SkiesANDUniverses/'
    
    SAM_name    = 'Galacticus'
    redshift    = 0.00
    #example for path+filename to read-in:
    myfilename  = common_path+'MDPL2_'+SAM_name+'_z_'+str("{0:.2f}".format(redshift))+'.hdf5'
    
    f = hdf5.File(myfilename, "r")

    
    #construct data types for the your chosen galaxy properties 
    dt={}
    print 'the following galaxy properties have been selected for\n\tSAM:\t ', \
           SAM_name, '\n\tredshift: ', redshift, '\n\tfilename:', myfilename, '\n'
    class unit(object):
        """This class serve to trick numpy-structured arrays to accept two
        field 'titles' with the same name by assigning each an object with unit 
        identifier!"""
        def __init__(self,unit):
            self.unit = unit
            
    for col in mycol:
        try:
            print 'ID', col, '\tproperty name:', gal_properties['name'+str(col)].ljust(18),\
                  'unit:', str('['+f[gal_properties['name'+str(col)]].attrs['unit']+']').ljust(25), \
                  'data type:', f[gal_properties['name'+str(col)]].dtype
                  
            dt.update({gal_properties['name'+str(col)]: (\
                      np.dtype(f[gal_properties['name'+str(col)]].dtype),\
                      col,\
                      unit(f[gal_properties['name'+str(col)]].attrs['unit']))})
        except:
            print 'galaxy property is not available for this SAM!'
            pass
        
    #construct a structured array which will be filled with the chosen galaxy
    #properties
    data_struct=np.zeros((f[dt.items()[0][0]].size ,), dtype=dt)                  
 
    print '\nprocessing ...\n'         
    #read all chosen galaxy properties to the structured array
    for name in data_struct.dtype.names:                  
        print 'name', name.ljust(18), str('['+f[name].attrs['unit']+']').ljust(25), \
              'ngal:', f[name].size, 
        data_struct[name][0:f[name].size] = f[name]
        print 'successfully read!'
        
print 'EOF - END OF FUN'
print '-------------------------'