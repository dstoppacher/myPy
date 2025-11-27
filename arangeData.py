from __future__ import print_function
import numpy as np
from astropy.io import fits
import h5py as hdf5
from time import time
import SAM_hdf5_lib as hdf5lib
import myLib as mL
import outputData as oD
import arangeData as aD
import os, sys
from os import listdir
from os.path import isfile, join
#import read_data_Cholla as rdC

#myOutput = oD.OutputData(config=False)

system_info = os.getcwd()
start = system_info.find('anaconda')
mycomp = system_info[:start]

class	ArangeData:
                          
    def readAnyFormat(
                self,
                config=True,
                mypath=False,
                myfilename=False,
                mypath_softlink=False,
                data_format='ASCII',
                data_shape='shaped',
                filename=False,
                delim=None,
                nr_col=0,
                id_col=0,
                nr_rows=500000,
                comment='#',
                skiprow=0,
                mydat_off=0,
                mydtype=np.str_,
                mydtype_after_read=False,
                id_col_array=[],
                config_array=[],
                name_conv_map=[],                
                catname='default',
                snapid=False,
                show_content=False,
                value_is_redshift=True,
                file_count=0,
                halocat_code=False,
                halo_id_col_array=[],
                start_fileID=0,
                nr_files_snapshot=1,
                end_fileID=1,
                scale_factor_map=None):


        def create_structed_array(set_nr_rows=False):
            
            if set_nr_rows!=False:
                nr_rows=set_nr_rows                
            
            mytypes={}
            a=0
            while a<int(id_col_array['nr_entries']):
                mytypes.update({id_col_array['name'+str(a)]: id_col_array['data_type'+str(a)]})                    
                a+=1
 
            dt = np.dtype([(k, mytypes[k]) for k in mytypes.keys()])

            try:
                structured_array=np.zeros((nr_rows,), dtype=dt)
                nr_entries=int(id_col_array['nr_entries'])
            except:
                structured_array=np.zeros((nr_rows,), dtype=dt)
                nr_entries=nr_col
                
            return structured_array               

        def read_simple_hdf52struct_array(path,
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
            
            f=hdf5.File(path, 'r')
            print('\nSimple read-in hdf5 data to structured array! --> \npath:', path)
            
            redshift    =f.attrs['redshift']
            scale_factor=f.attrs['scaleFactor']
            
            if cols==[]:
                cols=f.keys()
        
            dt = np.dtype([(k, f[k].dtype) for k in cols])
        
            data = np.zeros((f[cols[0]].size,), dtype=dt)
           
            #print(np.info(data))
        
            for name in cols:
                #print('property:', name,)
                for dim in range(len(f[cols[0]].shape)):
                    #print('dim:', dim)
                    try:
                        data[name][dim,:]=f[name][dim]
                    except:
                       data[name]=f[name][:] 
                                  
            f.close()
                
            return data, redshift, scale_factor

        def readANDConvert(convert):            
            print('###################################################################################\nreading and converting data formats ...')

            def readFITS():            
                print('###################################################################################\nreading fits ...')

                #mypath='/home/doris/anaconda/pro/data/CMASS_SAMS/CMASS_CATALOGUES/LSS/cmass-dr12v4-allmasses-spall-complete.dat.fits'
                hdulist = fits.open(mypath, memmap=True)                
                #hdulist.info()
                fits_header_map={}
#                hdr = hdulist[0].header                   
#                print(repr(hdr)             
#                hdr = hdulist[1].header 
#                print(repr(hdr)             
#                hdr = hdulist[2].header 
#                print(repr(hdr)              
 
                #print(hdulist
                #print(hdulist[2].header
                #print(hdulist[2].data
                exit()
                i=0
                a=1
                for entry in hdulist[2].header:
#                    try:
#                        print('a:', a, hdulist[1].header[i], hdulist[1].data[hdulist[1].header[i]].shape                      
#                        check_data_shape=False
#                        fits_header_map.update({'name'+str(a-1): hdulist[1].header[i], hdulist[1].header[i]: hdulist[1].data[hdulist[1].header[i]]})
#                        a+=1
#                    except:
#                        try:
#                            check_data_shape=True
#                            fits_header_map.update({'name'+str(a-1): hdulist[1].header[i], hdulist[1].header[i]: hdulist[1].data[hdulist[1].header[i]]})
#                            a+=1
#                        except:
#                            print('failed!'
                    i+=1 
                fits_header_map.update({'nr_entries': a})
                
                #print(fits_header_map
                #print(my_cols2extract
             
                with hdf5.File(mycomp+'anaconda/pro/data/'+catname+'/'+catname+'_'+config_array[catname+'_output_filename_code']+'.hdf5', 'w') as hf:
                    i=0
                    while i<fits_header_map['nr_entries']-1:
                       #print('i:', i, fits_header_map['name'+str(i)]

                        if np.any(my_cols2extract==fits_header_map['name'+str(i)]):
                            print('i:', i, 'fits name:', fits_header_map['name'+str(i)], fits_header_map[fits_header_map['name'+str(i)]].shape, 'hdf5 name:', name_mapping[fits_header_map['name'+str(i)]], end='')                            
                            try:
                                print(fits_header_map[fits_header_map['name'+str(i)]][0,:])
                                #fits_header_map[fits_header_map['name'+str(i)]][0:2,0],
                                dset = hf.create_dataset(name_mapping[fits_header_map['name'+str(i)]], (len(fits_header_map[fits_header_map['name'+str(i)]]),), dtype=np.float32)                                
                                dset[0::] = fits_header_map[fits_header_map['name'+str(i)]][:,0]
                                dset.attrs['colName'] = name_mapping[fits_header_map['name'+str(i)]]                                
                                print('--> first element done!')                                
                                a=1
                                while a<len(fits_header_map[fits_header_map['name'+str(i)]][0,:]) and check_data_shape==True:
                                    print(' a:', a, 'name:', name_mapping[fits_header_map['name'+str(i)]+str(a)], fits_header_map[fits_header_map['name'+str(i)]][0:2,a])
    
                                    dset = hf.create_dataset(name_mapping[fits_header_map['name'+str(i)]+str(a)], (len(fits_header_map[fits_header_map['name'+str(i)]][:,a]),), dtype=np.float32)
                                    dset[0::] = fits_header_map[fits_header_map['name'+str(i)]][:,a]
                                    dset.attrs['colName'] = name_mapping[fits_header_map['name'+str(i)]+str(a)]                                         
                                    a+=1                               
                            except:
                                dset = hf.create_dataset(name_mapping[fits_header_map['name'+str(i)]], (len(fits_header_map[fits_header_map['name'+str(i)]]),), dtype=np.float32)
                                if name_mapping[fits_header_map['name'+str(i)]].find('mstar')!=-1 or name_mapping[fits_header_map['name'+str(i)]].find('ssfr')!=-1: 
                                    print('10**\n')
                                    if name_mapping[fits_header_map['name'+str(i)]]=='mstar_char': 
                                        dset[0::] = 10**(fits_header_map[fits_header_map['name'+str(i)]]-0.03925)
                                    else:
                                        dset[0::] = 10**fits_header_map[fits_header_map['name'+str(i)]]
                                    
                                    
                                else:                                    
                                    dset[0::] = fits_header_map[fits_header_map['name'+str(i)]]
                                    print('\n')
                                dset.attrs['colName'] = name_mapping[fits_header_map['name'+str(i)]]
                        i+=1
                
                hdulist.close()
                exit()
            def convert2HDF5():
                print('###################################################################################\nconverting to HDF5 ...')
    
                mytypes={}
                a=0
                while a<int(id_col_array['nr_entries']):
                    mytypes.update({id_col_array['name'+str(a)]: id_col_array['data_type'+str(a)]})                    
                    a+=1
     
                dt = np.dtype([(k, mytypes[k]) for k in mytypes.keys()])
    
                try:
                    self.data_array=np.zeros((nr_rows,), dtype=dt)
                    nr_entries=int(id_col_array['nr_entries'])
                except:
                    self.data_array=np.zeros((nr_rows,), dtype=dt)
                    nr_entries=nr_col
                
                print('myfilename:', myfilename)
                f = hdf5.File(myfilename, "r")


                while i<int(name_mapping['nr_entries']):
                
                        try:
                            if config_array[catname+'_UNIT_CODE']=='MD':
                                print('i:', i, name_conv_map[id_col_array['name'+str(i)]], '-->', end='')
                                name=name_conv_map[id_col_array['name'+str(i)]]
                            else:
                                name=id_col_array['name'+str(i)]
                        except:
                            name=id_col_array['name'+str(i)]
        
                        print(id_col_array['name'+str(i)], 'size:', f[name].size, 'i:', i, 'id_col:', id_col_array[id_col_array['name'+str(i)]+'_col_id'])
                           
                        self.data_array[id_col_array['name'+str(i)]][0:f[name].size] = f[name]
                        
            
            def read2Array(fits_header_map):
                i=0
                count=0
                while i<fits_header_map['nr_entries']-1:
                    if np.any(my_cols2extract==fits_header_map['name'+str(i)]):
                        print(fits_header_map['name'+str(i)], fits_header_map[fits_header_map['name'+str(i)]].shape, name_mapping[fits_header_map['name'+str(i)]])
                        if count==0:
                            self.data_array = fits_header_map[fits_header_map['name'+str(i)]]                        
                            count=1
                        else:
                            self.data_array = np.column_stack((self.data_array, fits_header_map[fits_header_map['name'+str(i)]]))  
    
                    i+=1
                
                print(self.data_array[0:10,:])

            def caseSwitcher(convert):

                choose = {
                    'FITS2HDF5': FITS2HDF5,
                    'HDF52HDF5': HDF52HDF5,
                    'READ2ARRAY': READ2ARRAY
                    }
                    
                func = choose.get(convert)
                return func()
                
            def FITS2HDF5():
                readFITS()
            
            def READ2ARRAY():
                read2Array()                
    
            def HDF52HDF5():
                convert2HDF5()
            
            my_cols2extract = mL.multicolTestAlgorithm(config_array[catname+'_fits_map_names'], delimiter=';')
            my_cols2extract_mapping = mL.multicolTestAlgorithm(config_array[catname+'_fits_map_names_mapping'], delimiter=';')
            name_mapping={}
            #print(my_cols2extract
            #print(my_cols2extract_mapping
            i=0
            while i<len(my_cols2extract):
                #print('i:', i, my_cols2extract[i], my_cols2extract_mapping[i]
                name_mapping.update({my_cols2extract[i]: my_cols2extract_mapping[i]})
                i+=1
            
            name_mapping.update({'nr_entries': i})
            #print(name_mapping
           
            caseSwitcher(convert)
                 
            exit()
           
        def readBINARY():
            print('###################################################################################\nreading BINARY FILES GALAXY CATALOUGE+MHALO ...')
 
            #prin 'mypath in the end:', mypath
            print('nr_entries:', int(id_col_array['nr_entries']))
            self.data_array=np.zeros((nr_rows, int(id_col_array['nr_entries'])), dtype=np.float64)
                       
            start = time()
            count = 0
            count_new = 0          
            start_byte=8
            add_bytes = 8
            block_size=224
            data_block_col_size = 53

            def caseSwitcherReadBinary(name):

                choose = {
                    'ngalaxies': Ngalaxies,
                    'haloid': Haloid,
                    'hostid': Hostid,
                    'orphan': Orphan,
                    'mcold': Mcold,
                    'mstar': Mstar,
                    'mbh': Mbh,
                    'zgas': Zgas,
                    'sfr': Sfr,
                    'sfr_disk': skipThat,
                    'sfr_spheroid': Sfrspheriod,
                    'mstar_spheroid': Mstarsph,
                    'mhalo': Mhalo,
                    'spin': skipThat,
                    'mstar_disk': skipThat,
                    'x_pos': X_pos,
                    'y_pos': Y_pos,
                    'z_pos': Z_pos,
                    'x_vel': X_vel,
                    'y_vel': Y_vel,
                    'z_vel': Z_vel,
                    'mAB_dA_total_r': skipThat,
                    'mAB_dA_total_g': skipThat,
                    'MAB_dA_total_r': MAB_dA_total_r,
                    'MAB_dA_total_g': MAB_dA_total_g,
                    'mAB_total_r': skipThat,
                    'mAB_total_g': skipThat,
                    'MAB_total_r': MAB_total_r,
                    'MAB_total_g': MAB_total_g
                    }
                    
                func = choose.get(name)
                #print('func:', func
                return func()
                
            def Haloid():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]     = haloid[0]
                
            def Ngalaxies():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]  = ngalaxies[0]         
        
            def Hostid():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]     = data_block["galaxyhostid"]
                
            def Orphan():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]     = data_block["galaxy_is_orphan"]
            
            def Mcold():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]      = data_block["cols3_56"][:,6]
            
            def Mstar():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]      = data_block["cols3_56"][:,8]
            
            def Mbh():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]        = data_block["cols3_56"][:,9]
            
            def Zgas():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]       = data_block["cols3_56"][:,10]
            
            def Sfr():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]        = data_block["cols3_56"][:,13]
            
            def Mstarsph():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]   = data_block["cols3_56"][:,22]
            
            def Mhalo():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]      = data_block["cols3_56"][:,52]
                
            def X_pos():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]      = data_block["cols3_56"][:,0]

            def Y_pos():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]      = data_block["cols3_56"][:,1]

            def Z_pos():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]      = data_block["cols3_56"][:,2]

            def X_vel():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]      = data_block["cols3_56"][:,3]

            def Y_vel():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]      = data_block["cols3_56"][:,4]

            def Z_vel():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]      = data_block["cols3_56"][:,5]

            def MAB_dA_total_g():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]      = data_block["cols3_56"][:,43]

            def MAB_dA_total_r():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]      = data_block["cols3_56"][:,44]

            def MAB_total_g():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]      = data_block["cols3_56"][:,49]

            def MAB_total_r():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]      = data_block["cols3_56"][:,50]
                            
            def Sfrspheriod():
                self.data_array[count:count_new, str(id_col_array['col_id'+str(a)])]      = data_block["cols3_56"][:,14] 

                        
            def skipThat():
                return
         
            print('mypath', mypath)
            openfile = open(mypath,'rb')
            nhalos = np.fromstring(openfile.read(8), dtype = np.uint64)
            #print('nhalos:', nhalos[0]
            
            dt = np.dtype([
                        ('galaxyhostid', np.uint64),
                        ('galaxy_is_orphan', np.int32),
                        ('cols3_56', np.float32, data_block_col_size)
                        ])           
            
            i=0
            while i<nhalos[0]:
                openfile.seek(start_byte)
                ngalaxies = np.fromstring(openfile.read(8), dtype = np.uint64)
                start_byte+=add_bytes 

                openfile.seek(start_byte)
                haloid = np.fromstring(openfile.read(8), dtype = np.uint64)
                start_byte+=add_bytes 
                
                openfile.seek(start_byte)
                data_block = np.fromfile(openfile, dtype=dt, count=ngalaxies[0])
               
                count_new+=ngalaxies[0]
                count_new = count_new.astype(np.uint32)
                #print('ngalaxies:', ngalaxies[0]
                #print('count_new:', count_new
                
                a=0
                while a<id_col_array['nr_entries']:
                    caseSwitcherReadBinary(str(id_col_array['name'+str(a)])) 
                    a+=1
                
                            
                count=count_new               
                start_byte+=np.uint32(ngalaxies[0]*block_size)
                
                i+=1

            openfile.close()

            self.data_array = self.data_array[0:count,:]
            
            self.redshift=None
            self.scale_factor=None

        def readBINARY_SAGE(mypath):

            print('###################################################################################\nreading BINARY FILES GALAXY CATALOUGE+MHALO ...')
 
            SAM_binary_filestruct_map = {}
            SAM_binary_filestruct_map['catname']                      = catname
            SAM_binary_filestruct_map[catname+'_snapid']              = snapid
            
            data = self.readAnyFormat(config=False, mydtype=np.str_, mypath=mycomp+'anaconda/pro/data/SAGE_1Gpc/SAGE_1Gpc_file_names.txt', data_shape='shaped', comment='#')

            mytypes={}
            d=0
            while d<int(id_col_array['nr_entries']):
                mytypes.update({id_col_array['name'+str(d)]: id_col_array['data_type'+str(d)]})                    
                d+=1
 
            dt = np.dtype([(k, mytypes[k]) for k in mytypes.keys()])

            self.data_array=np.zeros((nr_rows,), dtype=dt)

            i=0
            a=0
            while i<data.size:
                #print('i:', i, 'a:', a, data[i].find(snapid), catname+'_filename'+str(a), data[i]
                if data.size==1:
                    SAM_binary_filestruct_map[catname+'_filename'+str(a)] = mypath+'/'+str(data)
                    a+=1
                else:
                    if data[i].find(snapid)==7 or data[i].find(snapid)==8:               
                        SAM_binary_filestruct_map[catname+'_filename'+str(a)] = mypath+'/'+data[i]
                        #print('chosen filename:', SAM_binary_filestruct_map[catname+'_filename'+str(a)]
                        a+=1
                i+=1      
           
            SAM_binary_filestruct_map[catname+'_nr_files'] = a        

            dt =    [
                    ('SnapNum'                      , np.int32),                    
                    ('orphan'                       , np.int32),   #Type                
                    ('GalaxyIndex'                  , np.int64),   #GalaxyIndex               
                    ('CentralGalaxyIndex'           , np.int64),   #CentralGalaxyIndex                
                    ('hostid'                       , np.int64),   #CtreesHaloID
                    ('TreeIndex'                    , np.int32),   #TreeIndex
                    ('haloid'                       , np.int64),   #CtreesCentralID
                    ('mergeType'                    , np.int32),                       
                    ('mergeIntoID'                  , np.int32),                        
                    ('mergeIntoSnapNum'             , np.int32),     
                    ('dT'                           , np.float32),                    
                    ('Pos'                          , (np.float32, 3)),           
                    ('Vel'                          , (np.float32, 3)),             
                    ('spinParameter'                , (np.float32, 3)), #Spin             
                    ('Len'                          , np.int32),                    
                    ('mhalo_200c'                        , np.float32), #Mvir                 
                    ('mhalo_cents'                  , np.float32), #CentralMvir                 
                    ('rvir'                         , np.float32), #Rvir                 
                    ('vvir'                         , np.float32), #Vvir                 
                    ('vmax'                         , np.float32), #Vmax               
                    ('vdisp'                      , np.float32), #Vdisp                 
                    ('mcold_disk'                   , np.float32), #ColdGas                 
                    ('mstar'                        , np.float32),  #StellarMass                
                    ('mstar_spheroid'               , np.float32),  #BulgeMass               
                    ('mhot'                         , np.float32),  #HotGas                
                    ('EjectedMass'                  , np.float32),                  
                    ('mbh'                          , np.float32),  #BlackHoleMass                
                    ('mstar_IC'                      , np.float32),  #IntraClusterStars                
                    ('Mzgas_disk'                   , np.float32),   #MetalsColdGas               
                    ('Mzstar'                       , np.float32),   #MetalsStellarMass               
                    ('Mzstar_spheroid'              , np.float32),   #MetalsBulgeMass              
                    ('Mzhot_halo'                   , np.float32),   #MetalsHotGas               
                    ('MetalsEjectedMass'            , np.float32),                  
                    ('MetalsIntraClusterStars'      , np.float32),                  
                    ('sfr_disk'                     , np.float32), #SfrDisk                 
                    ('sfr_spheroid'                 , np.float32), #SfrBulge                
                    ('zgas_disk'                     , np.float32), #SfrDiskZ                
                    ('zgas_spheroid'                 , np.float32), #SfrBulgeZ                
                    ('rdisk'                       , np.float32), #DiskRadius                  
                    ('Cooling'                      , np.float32),                  
                    ('Heating'                      , np.float32),
                    ('mbh_acc_quasar'                , np.float32), #QuasarModeBHaccretionMass
                    ('lastMajorM'                   , np.float32), #TimeOfLastMajorMerger
                    ('lastMinorM'                   , np.float32), #TimeOfLastMinorMerger
                    ('OutflowRate'                  , np.float32),
                    ('mean_age_stars'               , np.float32), #MeanStarAge
                    ('mhalo_sat'                   , np.float32), #infallMvir
                    ('vvir_sat'                   , np.float32),#infallVvir
                    ('vmax_sat'                   , np.float32) #infallVmax
                    ]
            
            
            names = [dt[c][0] for c in xrange(len(dt))]
            formats = [dt[c][1] for c in xrange(len(dt))]
            dt = np.dtype({'names':names, 'formats':formats}, align=True)
#
#            for name in names:
#                print('name:', name
#                
#            exit()
        
            start = time()
            count = 0
            count_new = 0          
        
            i=0
            while i<SAM_binary_filestruct_map[catname+'_nr_files']:
                #print('i:', i, SAM_binary_filestruct_map[catname+'_filename'+str(i)]
                openfile = open(SAM_binary_filestruct_map[catname+'_filename'+str(i)],'rb')
                Ntrees = np.fromfile(openfile, np.dtype(np.int32), 1)  # Read number of trees in file
                ngalaxies = np.fromfile(openfile, np.dtype(np.int32), 1)[0]  # Read number of gals in file.
                GalsPerTree = np.fromfile(openfile, np.dtype((np.int32, Ntrees)),1) # Read the number of gals in each tree

                data_block = np.fromfile(openfile, dt, ngalaxies) # Read in the galaxy structures
#                b=0
#                while b<20:
#                    print('b:', b, data_block[0][b]
#                    b+=1


                
                count_new+=ngalaxies
                #print('i:', i, 'ngalaxies:', ngalaxies, 'data_block.shape:', data_block.shape, 'count:', count, 'count_new:', count_new
                a=0
                while a<id_col_array['nr_entries']:
                    #print('a:', a, 'name:', id_col_array['name'+str(a)], 'id_col:', id_col_array['col_id'+str(a)]
                    if id_col_array['name'+str(a)]=='Z'\
                        or id_col_array['name'+str(a)]=='ssfr'\
                        or id_col_array['name'+str(a)]=='Tcons'\
                        or id_col_array['name'+str(a)]=='cgf'\
                        or id_col_array['name'+str(a)]=='zcold'\
                        or id_col_array['name'+str(a)]=='fbar'\
                        or id_col_array['name'+str(a)]=='bheff'\
                        or id_col_array['name'+str(a)]=='SHMR'\
                        or id_col_array['name'+str(a)]=='regionName':
                        pass                        
                    elif id_col_array['name'+str(a)]=='sfr':
                        #print('here: sfr!'
                        self.data_array[id_col_array['name'+str(a)]][count:count_new] = self.data_array[id_col_array['name'+str(a)]+'_spheroid'][count:count_new]+self.data_array[id_col_array['name'+str(a)]+'_disk'][count:count_new]
                    elif id_col_array['name'+str(a)].find('x_pos')!=-1:
                        #print('x:', data_block['Pos'][:,0]
                        self.data_array[id_col_array['name'+str(a)]][count:count_new] = data_block['Pos'][:,0]
                    elif id_col_array['name'+str(a)].find('y_pos')!=-1:
                        #print('y:', data_block['Pos'][:,1]
                        self.data_array[id_col_array['name'+str(a)]][count:count_new] = data_block['Pos'][:,1]
                    elif id_col_array['name'+str(a)].find('z_pos')!=-1:
                        #print('z:', data_block['Pos'][:,2]
                        self.data_array[id_col_array['name'+str(a)]][count:count_new] = data_block['Pos'][:,2]
                    elif id_col_array['name'+str(a)].find('x_vel')!=-1:
                        #print('x:', data_block['Pos'][:,0]
                        self.data_array[id_col_array['name'+str(a)]][count:count_new] = data_block['Vel'][:,0]
                    elif id_col_array['name'+str(a)].find('y_vel')!=-1:
                        #print('y:', data_block['Pos'][:,1]
                        self.data_array[id_col_array['name'+str(a)]][count:count_new] = data_block['Vel'][:,1]
                    elif id_col_array['name'+str(a)].find('z_vel')!=-1:
                        #print('z:', data_block['Pos'][:,2]
                        self.data_array[id_col_array['name'+str(a)]][count:count_new] = data_block['Vel'][:,2]
                    elif id_col_array['name'+str(a)]=='spinParameter':
                        #print('z:', data_block['Pos'][:,2]
                        #print('spin[,[0,1,2]', data_block['spinParameter'][0:10,[0,1,2]]
                        self.data_array[id_col_array['name'+str(a)]][count:count_new] = (data_block['spinParameter'][:,0]**2+data_block['spinParameter'][:,1]**2+data_block['spinParameter'][:,2]**2)**0.5 / (2.0**0.5 * data_block['rvir'] * data_block['vvir'])
                        #print(self.data_array[id_col_array['name'+str(a)]][count:count_new]
#                    elif id_col_array['name'+str(a)]=='zgas_disk':
#                        self.data_array[id_col_array['name'+str(a)]][count:count_new] = data_block['Mzgas_disk']/data_block['mcold_disk']
                    elif id_col_array['name'+str(a)]=='mstar_disk':
                        self.data_array[id_col_array['name'+str(a)]][count:count_new] = data_block['mstar']
                    elif id_col_array['name'+str(a)]=='Mzstar_disk':
                        self.data_array[id_col_array['name'+str(a)]][count:count_new] = data_block['Mzstar']  
                    elif id_col_array['name'+str(a)]=='mstar+IC':
                        self.data_array[id_col_array['name'+str(a)]][count:count_new] = data_block['mstar']                        
                    else:
                        #print('default!'
                        self.data_array[id_col_array['name'+str(a)]][count:count_new] = data_block[str(id_col_array['name'+str(a)])]
                    a+=1

                count+=ngalaxies                          
                i+=1


            openfile.close()

            self.data_array = self.data_array[:count]

            try:
                #write redshift in a column if Z is selected
                self.data_array['Z']=format(float(snapid), '.2f')
            except:
                pass
            
            self.redshift=np.float32(snapid)
            self.scale_factor = mL.redshift_to_expfactor(self.redshift)
            
            print(self.data_array.shape, self.redshift, self.scale_factor)


           
        def readASCII():

#            print('###################################################################################'
#            print('read ASCII'
#            print(' '
        
                                     
            if data_shape=='shaped':  
              
                self.data_array = np.genfromtxt(mypath, dtype=mydtype, delimiter=delim, comments=comment, skip_header=skiprow, encoding=None)
                    
            else:

                self.data_array = np.zeros((nr_rows, nr_col), dtype=np.float32)                

                #print('DATA FORMAT:', data_format, 'nr col:', nr_col, 'id_col:', id_col, mypath, 'rows 2 read:', nr_rows)
                
                check_path = os.path.exists(mypath)
                
                if check_path==True:                
                
                    i = 0
                    count = 0
                    line_count = 0
                    a=0
                    with open(mypath, 'r') as f: 
                        
                        for line in f:
                                
                                #if mypath==self.config_array['mydir'][0]+'62.5Mpc_SAG_AHF.txt': print('line_count:', line_count
                                if line_count<nr_rows:
                                    data = np.fromstring(line, dtype=np.double, sep='\n')#, count=2)#,  missing_values="0", filling_values="0")#, usecols=(self.config_array['id_col1'], self.config_array['id_col2']))
                                    j = 0
                                    #print('i start:', i
                                    if data.size==mydat_off and i!=0: 
                                        self.data_array[i-1,:]=self.data_array[i,:]
                                        i-=1
                                        
                                    
                                    while (j<data.size-1):
                                        
                                        if data.size==mydat_off:
                                            self.data_array[i,j] = data[j]
                                            if j==data.size-1: count+=1
                                            
                                        else:
                                            self.data_array[i,0:mydat_off] = self.data_array[i-1,0:mydat_off]
                                            self.data_array[i-1,j+mydat_off] = data[j]
                                           
                                        j+=1
                               
                                    i+=1
                                    line_count+=1 
                                    a+=1

        def read_MDGalaxies():
            
            print('###################################################################################\nread MultiDark-Galaxies HDF5 file format\n')
            
            #Galaxy properties availalble
            gal_props =['GalaxyType','HaloMass','HostHaloID',\
                        'LstarSDSSg','LstarSDSSi','LstarSDSSr','LstarSDSSu','LstarSDSSz',\
                        'MZgasDisk','MZgasSpheroid','MZhotHalo','MZstarDisk','MZstarSpheroid',\
                        'MainHaloID','Mbh','McoldDisk','McoldSpheroid','Mhot','MhotOutflow','MstarDisk','MstarSpheroid',\
                        'SFRdisk','SFRspheroid','SpinParameter',\
                        'Vx','Vy','Vz','X','Y','Z',\
                        'rbulge','rdisk','rhalf_mass']      

            cols=[] 

            myfilename=config_array[catname+'_use_store_registerpath']+'MDPL2_Galacticus_z_'+str(format(scale_factor_map[catname+'_redshift'+str(file_count)], '0.2f'))+'.hdf5'
            print('myfilename:', myfilename)
                
            data,self.redshift, self.scale_factor = read_simple_hdf52struct_array(myfilename)  

            #print(np.info(data)
            self.data_array=create_structed_array(set_nr_rows=data.size)       
            
            
            for a, item in enumerate(self.data_array.dtype.names):
                try:
                    print('a:', a, 'item:', item, name_conv_map[item])
                    self.data_array[item] = data[name_conv_map[item]]
                except:
                    print('DOES NOT EXIST here set to default!', item)
                    if item.find('AB')!=-1:
                        print('--> magnitudes ... set to 0.0!')
                        self.data_array[item] = 0.0
                    else:
                        print('--> ... set to -99!')
                        self.data_array[item] = -99

                
            #print(name_conv_map
            print(np.info(self.data_array))                                                     

        def readHDF5(crossmatch=''):
            
            print('###################################################################################\nread HDF5')

            mytypes={}
            a=0
            while a<int(id_col_array['nr_entries']):  
                #print('a:', a, id_col_array['name'+str(a)], id_col_array['data_type'+str(a)])
                mytypes.update({id_col_array['name'+str(a)]: id_col_array['data_type'+str(a)]})                    
                a+=1

            #print(mytypes)
 
            dt = np.dtype([(k, mytypes[k]) for k in mytypes.keys()])

            try:
                self.data_array=np.zeros((nr_rows,), dtype=dt)
                nr_entries=int(id_col_array['nr_entries'])
            except:
                self.data_array=np.zeros((nr_rows,), dtype=dt)
                nr_entries=nr_col
            #print('HERE!'
            #myfilename='/store/erebos/doris/Galacticus_1Gpc_z_0.09_CUT3_Contreras+13_mcold.hdf5'
            print('myfilename:', myfilename)
            
            f = hdf5.File(myfilename+crossmatch, "r")
                

            self.cat_attributes = {}                           
            i=0
            while i<len(f.attrs.keys()):
                #print(list(f.attrs.keys())[i], list(f.attrs.values())[i])
                self.cat_attributes[list(f.attrs.keys())[i]] = list(f.attrs.values())[i]
                i+=1
            print('+++++++++++++++++++++++++++\n')
            try:
                self.redshift = self.cat_attributes['redshift']
                self.scale_factor = self.cat_attributes['scaleFactor']
            except:
                self.redshift=False
                self.scale_factor=False

            i=0

            while i<nr_entries:
                try:                   
                #     try:
                    if config_array[catname+'_UNIT_CODE']=='MD':
                        #print('i:', i, name_conv_map[id_col_array['name'+str(i)]], '-->', end=' ')
                        #print('('+name_conv_map[id_col_array['name'+str(i)]]+','+
                        name=name_conv_map[id_col_array['name'+str(i)]]
                    else:
                        name=id_col_array['name'+str(i)]
                    # except:
                    #     name=id_col_array['name'+str(i)]
    
                    print(id_col_array['name'+str(i)], 'size:', f[name].size, 'i:', i, 'id_col:', id_col_array[id_col_array['name'+str(i)]+'_col_id'], end='')

                    self.data_array[id_col_array['name'+str(i)]][0:f[name].size] = f[name]
                    print(' --> default!')
                    
                    mysave_colname=name
                
                except:
                    print('-->', id_col_array['name'+str(i)], '/', name, 'i:', i, 'not excisting ...')
                                                                      
                    if id_col_array['name'+str(i)]=='Z':
                        print('--> enter redshift in column: Z!')
                        self.data_array[id_col_array['name'+str(i)]] = self.redshift

                    elif id_col_array['name'+str(i)].find('AB')!=-1:
                        print('--> magnitudes ... set to 0.0!')
                        try:
                            self.data_array[id_col_array['name'+str(i)]][0::] = 0.0
                        except:
                            self.data_array[id_col_array['name'+str(i)]] = 0.0
                    elif id_col_array['name'+str(i)].find('sample_key')!=-1:
                        print('--> string ... DO NOTHING!')
                    elif id_col_array['name'+str(i)].find('flag')!=-1:
                        print('--> flag ... set to FALSE/0!')
                        self.data_array[id_col_array['name'+str(i)]]   = 0                     
                    else:
                        print('here set to -99!', id_col_array['name'+str(i)])
                        try:
                            self.data_array[id_col_array['name'+str(i)]][0::] = -99
                        except:
                            self.data_array[id_col_array['name'+str(i)]] = -99
                i+=1
            

            self.data_array = self.data_array[:f[mysave_colname].size]
            print('self.data_array after read-in:', self.data_array.shape, 'redshift:', self.redshift, 'scale factor:', self.scale_factor)
            #print(self.data_array[0:3]


        def readROCKSTAR_ASCII(halocat_code):
            print('###################################################################################\nread ROCkSTAR ASCII FORMAT for catname: ', catname, 'snapid:', snapid)

            path = mypath+'out_'+str(snapid)+'.list'
            print('path:', path)
    

            #Properties from Rockstar CHOLLA 50Mpc run
            #ID DescID Mvir Vmax Vrms Rvir Rs Np X Y Z VX VY VZ JX JY JZ Spin rs_klypin Mvir_all M200b M200c M500c M2500c Xoff Voff spin_bullock 
            #b_to_a c_to_a A[x] A[y] A[z] b_to_a(500c) c_to_a(500c) A[x](500c) A[y](500c) A[z](500c) T/|U| M_pe_Behroozi M_pe_Diemer Halfmass_Radius


        def readEAGLE_ASCII(key):
                    
            import pandas as pd
            
            if key=='full':
                print('###################################################################################read EAGLE ASCII FORMAT\n catname: ', catname, 'snapid:', snapid)
                mypath=mycomp+'anaconda/pro/data/EAGLE_100Mpc/original_catalogs/EAGLE_100Mpc_full_catalog.dat'
                print('mypath:', mypath)
    
                data=pd.read_csv(mypath, skiprows=1, \
                                  names=['fofID', 'haloid', 'hostid', 'galaxyID', 'nSubhalos',  'mhalo_fof', 'mhalo_200c', 'r200c',\
                                         'x_pos_subhalo', 'y_pos_subhalo', 'z_pos_subhalo', 'mPart_DM', 'rhalf_DM', 'rhalf_DM_2D', 'mPart_stars', 'mstar_30kpc',\
                                         'rhalf_stars',  'rhalf_stars_2D', 'mPart_gas', 'rhalf_gas', 'rhalf_gas_2D', 'mPart_bh', 'rhalf_bh', 'rhalf_bh_2D',\
                                         'vmax', 'rVmax', 'x_pos', 'y_pos', 'z_pos', 'x_vel', 'y_vel', 'z_vel', 'spinGas_x', 'spinGas_y', 'spinGas_z',\
                                         'mbh', 'vdisp', 'mstar_birth', 'sfr_total', 'Etot',  'Ekin', 'Etherm'] , sep=' ')
                    
                # names_props=['fofID', 'haloid', 'hostid', 'galaxyID', 'nSubhalos',  'mhalo_fof', 'mhalo_200c', 'r200c',\
                #            'x_pos_subhalo', 'y_pos_subhalo', 'z_pos_subhalo', 'mPart_DM', 'rhalf_DM', 'rhalf_DM_2D', 'mPart_stars', 'mstar_30kpc',\
                #            'rhalf_stars',  'rhalf_stars_2D', 'mPart_gas', 'rhalf_gas', 'rhalf_gas_2D', 'mPart_bh', 'rhalf_bh', 'rhalf_bh_2D',\
                #            'vmax', 'rVmax', 'x_vel', 'y_vel', 'z_vel', 'spinGas_x', 'spinGas_y', 'spinGas_z',\
                #            'vdisp', 'sfr_total']
                # data = data[names_props]
                

            elif key=='full_trees':
                
                print('###################################################################################')
                print('read EAGLE ASCII FULT MERGER TREE FORM DATA BASE\n catname: ', catname, 'snapid:', snapid)
                mypath=mycomp+'anaconda/pro/data/EAGLE_ev_100Mpc/EAGLE_ev_100Mpc_full_centrals_main_progenitor_trees.txt'
                print('mypath:', mypath)
    
                data= pd.read_csv(mypath, skiprows=1,\
                                  names=['fofID',  'progFofID',  'redshift',  'haloid',  'hostid',  'galaxyID',  'lastProgID',  'topLeafID',\
                                         'mhalo_200c',  'r200c',  'nSubhalos',  'x_pos',  'y_pos',  'z_pos',  'x_pos_subhalo',  'y_pos_subhalo',  'z_pos_subhalo',\
                                         'x_pos_cof',  'y_pos_cof',  'z_pos_cof',  'mPart_DM',  'mPart_bh',  'mPart_gas',  'mPart_stars',\
                                         'sfr',  'vmax',  'rVmax',  'vdisp',  'mstar_birth',  'mbh',  'bh_acc_rate',\
                                         'mstar_30kpc',  'sfr_30kpc',  'vdisp_30kpc',  'mgas_30kpc',  'mhalo_30kpc',\
                                         'rhalf_stars',  'rhalf_gas',  'mgas_SF',  'mgas_NSF',  'spinGasNSF_x',  'spinGasNSF_y',  'spinGasNSF_z',\
                                         'spinGasSF_x',  'spinGasSF_y',  'spinGasSF_z',  'spinStars_x',  'spinStars_y',  'spinStars_z',\
                                         'z_mean_birth_stars',  'age_mean_stars',  'x_vel',  'y_vel',  'z_vel'],\
                                         sep=',')

         
                 
            elif key=='Patricia1':
                """
                Catalog created by Patricia Tissera
                name:           REF100_generaltable_z0-header.dat
                sim:            EAGLE 100Mpc/h, fiducial model 1504^3 particles, Gadget-3, Planck cosmology, Chabrier IMF
                description:    catalog used for Yetli's project on void galaxies
                properties:     jsub, ndiskt, nopts, noptg, sfr, mstartotal, DT_stars, DT_gas, r80, M200, R200, ssfr, refftotal, MBH, 
                                MBHrate, idhalo ,dtgx,sfrgastota
                    
                    Property description
                    ================================
                    jsub:       internal index
                    ndiskt:     number of stellar particles in the discs within 1.5 r80 (or ropt)
                    nopts:      number of stars particles within 1.5 ropt
                    noptg:      number of gas particles withn 1.5 ropt
                    sfr:        star formation rate estimated by using stellar particles younger than 2e9 yrs (Msun/yr).
                    mstartotal: stellar mass within 1.5 ropt (Msun)
                    DT_stars:   disc-to-total stellar mass >0.4
                    DT_gas:     disc-to-total gas mass >0.4 most of the gas
                    ropt/r80:   radius that enclosed 83% of the stellar mass (Mpc)
                    M200:       virial mass in Msun
                    R200:       virial radius in Mpc
                    sSFR:       Specific star foramtion rate defined sfr/mstars (yr-1)
                    refftotal:  half-mass stellar radius (all particles within 1.5 ropt has been used) (Mpc)
                    MBH:        mass of the black hole (Msun).
                    MBHrate1:   rate of accretion (Msun/yr) 
                    idhalo:     ID of the FoF group
                    DTstars1:   the same as DT_stars but using a  circularity >0.5
                    sfrgastotal: estimating using the instantaneous star formation rate of gas particles  within 1.5 ropt (Msun/yr)
                     
    
            """           
                
                print('###################################################################################')
                print('read EAGLE ASCII FORMAT\n catname: ', catname, 'snapid:', snapid)
                mypath=mycomp+'anaconda/pro/data/EAGLE_100Mpc/original_catalogs/REF100_generaltable_z0-header.dat'
                print('mypath:', mypath)
                                       
                data= pd.read_csv(mypath, skiprows=19,\
                                   names=['jsub','np_disk_1.5ropt','np_stars_1.5ropt','np_gas_1.5ropt','sfr','mstar_1.5ropt','DvT_stars_c0.4','DvT_gas_c0.4','ropt','mhalo_200c','r200c',\
                                          'ssfr','rhalf_stars_1.5ropt', 'mbh', 'bh_acc_rate', 'fofID', 'DvT_stars_c0.5', 'sfr_int_gas'],
                                   sep=' ')           
                    
  
                    
            elif key=='Patricia2':
                 """
                 Catalog created by Patricia Tissera and sent March 2023, it differce a little form the first catalog she sent (DTgas circularity
                 is now c>0.5), 1/3 less stellar particles in the disk).
                 name:           REF25_generaltable_z0-header.dat
                 sim:            EAGLE 100Mpc/h, fiducial model 1504^3 particles, Gadget-3, Planck cosmology, Chabrier IMF
                 description:    catalog used for Yetli's project on void galaxies
                 properties:     jsub, ndiskt, nopts, noptg, sfr, mstartotal, DT_stars, DT_gas, r80, M200, R200, ssfr, refftotal, MBH, 
                                 MBHrate, idhalo ,dtgx,sfrgastota
                     
                     Property description
                     ================================
                     jsub:       internal index
                     ndiskt:     number of stellar particles in the discs within 1.5 r80 (or ropt) (1/3 less than in the catalog before!)
                     nopts:      number of stars particles within 1.5 ropt
                     noptg:      number of gas particles withn 1.5 ropt
                     sfr:        star formation rate estimated by using stellar particles younger than 2e9 yrs (Msun/yr).
                     mstartotal: stellar mass within 1.5 ropt (Msun)
                     DT_stars_c0.4:   disc-to-total stellar mass >0.4
                     DT_gas_c0.5:     disc-to-total gas mass >0.5 most of the gas (before was c>0.4)
                     ropt/r80:   radius that enclosed 83% of the stellar mass (Mpc)
                     M200:       virial mass in Msun
                     R200:       virial radius in Mpc
                     sSFR:       Specific star foramtion rate defined sfr/mstars (yr-1)
                     refftotal:  half-mass stellar radius (all particles within 1.5 ropt has been used) (Mpc)
                     MBH:        mass of the black hole (Msun).
                     MBHrate1:   rate of accretion (Msun/yr) 
                     idhalo:     ID of the FoF group
                     DTstars1:   the same as DT_stars but using a  circularity >0.5
                     sfrgastotal: estimating using the instantaneous star formation rate of gas particles  within 1.5 ropt (Msun/yr)
                     reffgas:    effective radius (Mpc)
                     jmodgas:    angular momentuio with 1.5 ropt normalized to Mgas
                     jmodstars:  with 1.5 ropt normalized to Mstars
            
             """           
                 
                 print('###################################################################################')
                 print('read EAGLE ASCII FORMAT\n catname: ', catname, 'snapid:', snapid)
                 mypath=mycomp+'anaconda/pro/data/EAGLE_100Mpc/REF25_generaltable_z0.dat'
                 print('mypath:', mypath)
            
                 data= pd.read_csv(mypath, skiprows=1,\
                                   names=['jsub','np_disk_1.5ropt','np_stars_1.5ropt','np_gas_1.5ropt','sfr','mstar_1.5ropt','DvT_stars_c0.4','DvT_gas_c0.5','ropt','mhalo_200c','r200c',\
                                          'ssfr','rhalf_stars_1.5ropt', 'mbh', 'bh_acc_rate', 'fofID', 'DvT_stars_c0.5', 'sfr_int_gas', 'ropt_gas', 'angM_norm_1.5ropt_stars', 'angM_norm_1.5ropt_gas'],
                                   sep=' ')
                
            elif key=='Yetli':
                 """Catalog with Mgas and environment send January 2023
             """           
                 
                 print('###################################################################################')
                 print('read EAGLE ASCII FORMAT\n catname: ', catname, 'snapid:', snapid)
                 mypath=mycomp+'anaconda/pro/data/EAGLE_100Mpc/Catalogue_Galaxies_Doris_January2023.dat'
                 print('mypath:', mypath)
            
                 data= pd.read_csv(mypath, skiprows=1,\
                                   names=['haloid', 'hostid', 'envr', 'lastMajorM', 'lastMinorM',\
                                          'Mgas_HI_30kpc_GK11', 'Mgas_HII_30kpc_GK11', 'Mgas_HI_30kpc_K13', 'Mgas_HII_30kpc_K13',\
                                          'Mgas_HI_60kpc_GK11', 'Mgas_HII_60kpc_GK11', 'Mgas_HI_60kpc_K13', 'Mgas_HII_60kpc_K13',\
                                           'MAB_total_u', 'MAB_total_i', 'MAB_total_r', 'MAB_total_z', 'MAB_total_g', 'mstar_30kpc'],
                                   sep=' ')
                     

                    
            elif key=='Silvio':
                """
                Catalog created by Patricia Tissera for Silvio (handed to me March 2023), reference: Van De Sande, J. 2019
                name:           EAGLE_L100N1504_ForSilvio
                sim:            EAGLE 100Mpc/h, fiducial model 1504^3 particles, Gadget-3, Planck cosmology, Chabrier IMF
                description:    catalog has kinematic information and Hydrogen Masses, 54 properties in total
                properties:     GroupNumber SubGroupNumber stellar_mass[Msun]
                                r50(2D)[pkpc] v_dispersion_2D[km/s, edge-on, r50_2D] V/S(r50_2D, edge-on) v_dispersion_2D[km/s, random, r50_2D] V/S(r50_2D, random)  ellipticity(r50_2D, edge-on) ellipticity(r50_2D, random
                                B/T_kinematics(Lagos17b) n_sersic(3D) lambdaR(r50_2D, edge-on) lambdaR(r50_2D, random)
                                MHI(60kpc,GK11) MHI(60kpc,K13) MH2(60kpc,GK11) MH2(60kpc,K13) 
                                v_dispersion_2D[km/s, edge-on, 2r50_2D] V/S(2r50_2D, edge-on)  v_dispersion_2D[km/s, random, 2r50_2D] V/S(2r50_2D, random)
                                ellipticity(2r50_2D, edge-on) ellipticity(2r50_2D, random) lambdaR(2r50_2D, edge-on) lambdaR(2r50_2D, random)
                                stellar_age(rband,r50_2D,[Gyr]) stellar_age(rband,2r50_2D,[Gyr]) 
                                SFR[Msun/yr] GroupMass[Msun] 
                                LambdaR_lasttimecentral(r50_2D, random) LambdaR_lasttimecentral(2r50_2D, random) SSFR_lastimecentral(yr^-1) 
                                SFR_lastimecentral(Msun yr^-1) stellar_mass_lastimecentral(Msun) V/S_lasttimecentral(r50_2D, random) 
                                V/S_lasttimecentral(2r50_2D, random) V/S_lasttimecentral(r50_2D, edge-on) V/S_lasttimecentral(2r50_2D, edge-on)
                                LookbackTime_LastCentral[Gyr]
                                Pos_x[Mpc] Pos_y[Mpc] Pos_z[Mpc] Vel_x[km/s] Vel_y[km/s] Vel_z[km/s]
                                DistToCentre/R200crit SurfaceDensN10[Mpc^-2] SurfaceDensN7[Mpc^-2] SurfaceDensN5[Mpc^-2]
                                LookbackTime_LastMerger[Gyr] MassRatio_LastMerger r50_SM(3D)[pkpc] r50_SFR(3)[pkpc]
                    
                    Property description
                    ================================
                    GroupNumber: FoFID
                    SubGroupNumber: Satellite ID
                    stellar_mass[Msun]
                    r50(2D)[pkpc]: half mass radius projected 2D
                    v_dispersion_2D[km/s,edge-on,r50_2D]
                    V/S(r50_2D,edge-on)
                    v_dispersion_2D[km/s,random,r50_2D]
                    V/S(r50_2D,random)
                    ellipticity(r50_2D,edge-on)
                    ellipticity(r50_2D,random)
                    B/T_kinematics(Lagos17b)
                    n_sersic(3D)
                    lambdaR(r50_2D,edge-on): angular momentum
                    lambdaR(r50_2D,random)
                    MHI(60kpc,GK11)
                    MHI(60kpc,K13)
                    MH2(60kpc,GK11)
                    MH2(60kpc,K13)
                    v_dispersion_2D[km/s,edge-on,2r50_2D]
                    V/S(2r50_2D,edge-on)
                    v_dispersion_2D[km/s,random,2r50_2D]
                    V/S(2r50_2D,random)
                    ellipticity(2r50_2D,edge-on)
                    ellipticity(2r50_2D,random)
                    lambdaR(2r50_2D,edge-on)
                    lambdaR(2r50_2D,random)
                    stellar_age(rband,r50_2D,[Gyr])
                    stellar_age(rband,2r50_2D,[Gyr])
                    SFR[Msun/yr]
                    GroupMass[Msun]: virial halo mass
                    LambdaR_lasttimecentral(r50_2D,random)
                    LambdaR_lasttimecentral(2r50_2D,random)
                    SSFR_lastimecentral(yr^-1)
                    SFR_lastimecentral(Msunyr^-1)
                    stellar_mass_lastimecentral(Msun)
                    V/S_lasttimecentral(r50_2D,random)
                    V/S_lasttimecentral(2r50_2D,random)
                    V/S_lasttimecentral(r50_2D,edge-on)
                    V/S_lasttimecentral(2r50_2D,edge-on)
                    LookbackTime_LastCentral[Gyr]
                    Pos_x[Mpc]
                    Pos_y[Mpc]
                    Pos_z[Mpc]
                    Vel_x[km/s]
                    Vel_y[km/s]
                    Vel_z[km/s]
                    DistToCentre/R200crit: distance to center vs. r200c
                    SurfaceDensN10[Mpc^-2]
                    SurfaceDensN7[Mpc^-2]
                    SurfaceDensN5[Mpc^-2]
                    LookbackTime_LastMerger[Gyr]
                    MassRatio_LastMerger
                    r50_SM(3D)[pkpc]
                    r50_SFR(3)[pkpc]                 
    
            """           
                
                print('###################################################################################')
                print('read EAGLE ASCII FORMAT KINEMATIC CATALOG\n catname: ', catname, 'snapid:', snapid)
                mypath=mycomp+'anaconda/pro/data/EAGLE_100Mpc/EAGLE_L100N1504_ForSilvio.dat'
                print('mypath:', mypath)
    
                data= pd.read_csv(mypath, skiprows=1,\
                                  names=['fofID','hostid','mstar','rhalf_stars_2D',\
                                         'vdisp_r502D_edgeOn', 'VvS_r502D_edgeOn', 'vdisp_r502D_random', 'VvS_r502D_random',\
                                         'ell_r502D_edgeOn','ell_r502D_random', 'BvT', 'Sersic_n',\
                                         'lambda_r502D_edgeOn','lambda_r502D_random', 'Mgas_HI_60kpc_GK11', 'Mgas_HI_60kpc_K13', 'Mgas_HII_60kpc_GK11', 'Mgas_HII_60kpc_K13',\
                                         'vdisp_2r502D_edgeOn', 'VvS_2r502D_edgeOn', 'vdisp_2r502D_random', 'VvS_2r502D_random', \
                                         'ell_2r502D_edgeOn', 'ell_2r502D_random',\
                                         'lambda_2r502D_edgeOn', 'lambda_2r502D_random', 'age_stars_rband_r502D', 'age_stars_rband_2r502D',\
                                         'sfr', 'mhalo',\
                                         'lambda_r502D_random_ltc', 'lambda_2r502D_random_ltc', 'ssfr_ltc', 'sfr_ltc', 'mstar_ltc',\
                                         'VvS_r502D_random_ltc', 'VvS_2r502D_random_ltc', 'VvS_r502D_edgeOn_ltc', 'VvS_2r502D_edgeOn_ltc', 'lastTimeCentral',\
                                         'x_pos', 'y_pos', 'z_pos', 'x_vel', 'y_vel', 'z_vel',\
                                         'rcenter2r200c', 'sDensity_N10', 'sDensity_N7', 'sDensity_N5', 'lastMerger', 'ratio_lastM', 'rhalf_stars_3D', 'rhalf_sfr_3D'],\
                                         sep=' ')
                    
            else:
                print('No key set! --> Go back to main!')
                exit()


            self.data_array = mL.df_to_sarray(data)

            #import numpy.lib.recfunctions as rcfuncs
            #self.data_array = rcfuncs.append_fields([self.data_array], ['orphan','nsats'] ,[np.zeros(self.data_array.size,),np.zeros(self.data_array.size,)], usemask=False)

            #self.data_array=self.data_array[np.where(self.data_array['redshift']<0.0001)[0][:]]
            self.data_array=self.data_array[np.where(self.data_array['hostid']==0)[0][:]]
            #print('ncentrals:', self.data_array.size
            #self.data_array['fofID']-=28000000000000
            #rint self.data_array
            print(np.info(self.data_array))
            #exit() 

        def readASCII_TTH():
            """300 Properties.txt
                """           
            import pandas as pd
            print('###################################################################################')
            print('read Galacticus 300 ASCII FORMAT\n catname: ', catname, 'snapid:', snapid)
            mypath=mycomp+'anaconda/pro/data/300/March2025/SFH_Properties/SFH_Properties_300_z_0.00__props.txt'
            print('mypath:', mypath)
       
            data= pd.read_csv(mypath, skiprows=2,\
                              names=['regionID', 'haloid', 'hostid', 'parentIndex', 'orphan', 'nsats', 'x_pos', 'y_pos', 'z_pos',\
                                     'mhalo', 'mstar', 'SHMR', 'mcold', 'mhot', 'mbh', 'mhot_outflow', 'Mzgas', 'Mzstar', 'Mzhot_halo', 
                                     'Mzhot_outflowHalo', 'zcold', 'zstar', 'sfr', 'ssfr', 'Tcons', 'bheff', 'rhalf_mass', 'rhotHalo',
                                     'cgf', 'jbar', 'jhotHalo', 'joutHotHalo', 'L_SDSS_dA_total_g', 'L_SDSS_dA_total_r', 'L_SDSS_dA_total_i',
                                     'mAB_dA_total_g', 'mAB_dA_total_r', 'mAB_dA_total_i',
                                     'MAB_dA_total_g', 'MAB_dA_total_r', 'MAB_dA_total_i', 
                                     'mAB_dA_total_cut_r_i', 'mAB_dA_total_cut_g_r', 'mAB_dA_total_cut_g_i'],
                              sep=' ')   
 
            self.data_array = mL.df_to_sarray(data)
            
            import numpy.lib.recfunctions as rcfuncs
            self.data_array = rcfuncs.append_fields([self.data_array], ['redshift','flag_Del', 'flag_Seb', 'haloidz0', 'DelSHMR', 'Delz', 'Delcol', 'z50Mhalo', 'z50Mstar', 'fbar'],\
                                                    [np.zeros(self.data_array.size,),self.data_array['orphan'],self.data_array['orphan'],self.data_array['haloid'],np.zeros(self.data_array.size,),np.zeros(self.data_array.size,),np.zeros(self.data_array.size,),np.zeros(self.data_array.size,),np.zeros(self.data_array.size,),np.zeros(self.data_array.size,)], usemask=False)


            self.data_array[['flag_Del','flag_Seb', 'haloidz0']] = -99
            #print(self.data_array.size, len(self.data_array[0]))              

        def readCHOLLAHDF5(halocat_code):
            
            def get_domain_block( proc_grid, box_size, grid_size ):
                np_x, np_y, np_z = proc_grid
                Lx, Ly, Lz = box_size
                nx_g, ny_g, nz_g = grid_size
                dx, dy, dz = Lx/np_x, Ly/np_y, Lz/np_z
                nx_l, ny_l, nz_l = nx_g//np_x, ny_g//np_z, nz_g//np_z
                nprocs = np_x * np_y * np_z
                domain = {}
                domain['global'] = {}
                domain['global']['dx'] = dx
                domain['global']['dz'] = dz
                for k in range(np_z):
                    for j in range(np_y):
                        for i in range(np_x):
                            pId = i + j*np_x + k*np_x*np_y
                            domain[pId] = { 'box':{}, 'grid':{} }
                            xMin, xMax = i*dx, (i+1)*dx
                            yMin, yMax = j*dy, (j+1)*dy
                            zMin, zMax = k*dz, (k+1)*dz
                            domain[pId]['box']['x'] = [xMin, xMax]
                            domain[pId]['box']['y'] = [yMin, yMax]
                            domain[pId]['box']['z'] = [zMin, zMax]
                            domain[pId]['box']['dx'] = dx
                            domain[pId]['box']['dy'] = dy
                            domain[pId]['box']['dz'] = dz
                            domain[pId]['box']['center_x'] = ( xMin + xMax )/2.
                            domain[pId]['box']['center_y'] = ( yMin + yMax )/2.
                            domain[pId]['box']['center_z'] = ( zMin + zMax )/2.
                            gxMin, gxMax = i*nx_l, (i+1)*nx_l
                            gyMin, gyMax = j*ny_l, (j+1)*ny_l
                            gzMin, gzMax = k*nz_l, (k+1)*nz_l
                            domain[pId]['grid']['x'] = [gxMin, gxMax]
                            domain[pId]['grid']['y'] = [gyMin, gyMax]
                            domain[pId]['grid']['z'] = [gzMin, gzMax]
                return domain

            def select_procid( proc_id, subgrid, domain, ids, ax ):
                domain_l, domain_r = domain
                subgrid_l, subgrid_r = subgrid
                
                if domain_l <= subgrid_l and domain_r > subgrid_l:
                    ids.append(proc_id)
                if domain_l >= subgrid_l and domain_r <= subgrid_r:
                    ids.append(proc_id)
                if domain_l < subgrid_r and domain_r >= subgrid_r:
                    ids.append(proc_id)

            def select_ids_to_load( subgrid, domain, proc_grid ):
                subgrid_x, subgrid_y, subgrid_z = subgrid
                nprocs = proc_grid[0] * proc_grid[1] * proc_grid[2]
                ids_x, ids_y, ids_z = [], [], []
                
                print('subrid x/y/z:', subgrid_x, subgrid_y, subgrid_z)
                
                for proc_id in range(nprocs):
                    domain_local = domain[proc_id]
                    
                    domain_x = domain_local['grid']['x']
                    domain_y = domain_local['grid']['y']
                    domain_z = domain_local['grid']['z']
                    
                    select_procid( proc_id, subgrid_x, domain_x, ids_x, 'x' )
                    select_procid( proc_id, subgrid_y, domain_y, ids_y, 'y' )
                    select_procid( proc_id, subgrid_z, domain_z, ids_z, 'z' )
                    
                set_x = set(ids_x)
                set_y = set(ids_y)
                set_z = set(ids_z)
                set_ids = (set_x.intersection(set_y)).intersection(set_z )
                return list(set_ids)

            def load_snapshot_data_distributed( nSnap, inDir, data_type, fields, subgrid,  precision, proc_grid,  box_size, grid_size, show_progess=True ):
                # Get the doamin domain_decomposition
                domain = get_domain_block( proc_grid, box_size, grid_size )
                
                # Find the ids to load 
                ids_to_load = select_ids_to_load( subgrid, domain, proc_grid )
                print('ids to load:', ids_to_load)
                print('Loading Snapshot: {0}'.format(nSnap))
                
                #Find the boundaries of the volume to load
                domains = { 'x':{'l':[], 'r':[]}, 'y':{'l':[], 'r':[]}, 'z':{'l':[], 'r':[]}, }
                
                for id in ids_to_load:
                    for ax in list(domains.keys()):
                        d_l, d_r = domain[id]['grid'][ax]
                        domains[ax]['l'].append(d_l)
                        domains[ax]['r'].append(d_r)
                
                boundaries = {}
                for ax in list(domains.keys()):
                    boundaries[ax] = [ min(domains[ax]['l']),  max(domains[ax]['r']) ]
                    
                # Get the size of the volume to load
                nx = int(boundaries['x'][1] - boundaries['x'][0])    
                ny = int(boundaries['y'][1] - boundaries['y'][0])    
                nz = int(boundaries['z'][1] - boundaries['z'][0]) 
                
                dims_all = [ nx, ny, nz ]
                data_out = {}
                data_out[data_type] = {}
                
                for field in fields:
                    
                    data_particels = False
                    if field in ['pos_x', 'pos_y', 'pos_z', 'vel_x', 'vel_y', 'vel_z']: data_particels = True 
                    if not data_particels: data_all = np.zeros( dims_all, dtype=precision )
                    else: data_all = []
                    
                    added_header = False
                    n_to_load = len(ids_to_load)
                    for i, nBox in enumerate(ids_to_load):
                        name_base = 'h5'
                        if data_type == 'particles': inFileName = '{0}_particles.{1}.{2}'.format(nSnap, name_base, nBox)
                        if data_type == 'hydro': inFileName = '{0}.{1}.{2}'.format(nSnap, name_base, nBox)
                        
                        inFile = hdf5.File( inDir + inFileName, 'r')
                        available_fields = inFile.keys()
                        head = inFile.attrs
                        
                        if added_header == False:
                            print(' Loading: '+inDir+inFileName) 
                            #print(' Available Fields:  ', available_fields
                            for h_key in list(head.keys()):
                                if h_key in ['dims', 'dims_local', 'offset', 'bounds', 'domain', 'dx', ]: continue
                                data_out[h_key] = head[h_key][0]
                                if h_key == 'current_z':
                                    print(' current_z: {0}'.format( data_out[h_key]))
                                    self.redshift = float(format(data_out[h_key], '0.6f'))
                                    #print('here: redshift', self.redshift
                                    
                                if h_key == 'current_a':
                                    #print(' current_a: {0}'.format( data_out[h_key])
                                    self.scale_factor = float(format(data_out[h_key], '0.6f'))
                                    #print('here: scale factor', self.scale_factor
                            added_header = True
                            
                        if show_progess:
                            terminalString  = '\r Loading File: {0}/{1}   {2}'.format(i, n_to_load, field)
                            sys.stdout. write(terminalString)
                            sys.stdout.flush()
                            
                        if not data_particels:
                            procStart_x, procStart_y, procStart_z = head['offset']
                            procEnd_x, procEnd_y, procEnd_z = head['offset'] + head['dims_local']
                            # Substract the offsets
                            procStart_x -= boundaries['x'][0]
                            procEnd_x   -= boundaries['x'][0]
                            procStart_y -= boundaries['y'][0]
                            procEnd_y   -= boundaries['y'][0]
                            procStart_z -= boundaries['z'][0]
                            procEnd_z   -= boundaries['z'][0]
                            procStart_x, procEnd_x = int(procStart_x), int(procEnd_x)
                            procStart_y, procEnd_y = int(procStart_y), int(procEnd_y)
                            procStart_z, procEnd_z = int(procStart_z), int(procEnd_z)
                            data_local = inFile[field][...]
                            data_all[ procStart_x:procEnd_x, procStart_y:procEnd_y, procStart_z:procEnd_z] = data_local
                            
                        else:
                            data_local = inFile[field][...]
                            data_all.append( data_local )
                        
                    if not data_particels:
                        # Trim off the excess data on the boundaries:
                        trim_x_l = subgrid[0][0] - boundaries['x'][0]
                        trim_x_r = boundaries['x'][1] - subgrid[0][1]  
                        trim_y_l = subgrid[1][0] - boundaries['y'][0]
                        trim_y_r = boundaries['y'][1] - subgrid[1][1]  
                        trim_z_l = subgrid[2][0] - boundaries['z'][0]
                        trim_z_r = boundaries['z'][1] - subgrid[2][1]  
                        trim_x_l, trim_x_r = int(trim_x_l), int(trim_x_r) 
                        trim_y_l, trim_y_r = int(trim_y_l), int(trim_y_r) 
                        trim_z_l, trim_z_r = int(trim_z_l), int(trim_z_r) 
                        data_output = data_all[trim_x_l:nx-trim_x_r, trim_y_l:ny-trim_y_r, trim_z_l:nz-trim_z_r,  ]
                        data_out[data_type][field] = data_output
                    else:
                        data_all = np.concatenate( data_all )
                        data_out[data_type][field] = data_all
                        if show_progess: print(' ')
                return data_out

            def read_and_write_attributes(file_struc_info, res):
                
                output_filename_grid=mycomp+'/anaconda/pro/data/Cholla'+res+'_50Mpc/Cholla'+res+'_grid.txt'
                myOutput.writeIntoFile(
                           output_filename_grid,
                           [''],
                           myheader='CHOLLA data resolution: '+res+' from Bruno \n(1) filename (2) a (3) z'\
                           +'(4) dt (5) dx (6) n_fields (6) n_step (7) n_fields (8) offset (9) t (10) dims (11) dims_local (12) bounds and/or n_procs',
                           append_mytext=False,
                           data_is_string=False,
                           data_format='%s')            
             
                i=0
                while i<int(file_struc_info[catname+'_nr_files']):           
                                    
                    if i==i_break and end_fileID!='False':
                       break
                    path = file_struc_info[catname+'_filename'+str(i)]
                    #path='/store/multidark/NewMD_3840_Planck1/Galacticus/latest/job0/the_trees_0_1000000_0_results.hdf5'
                    #path='/data3/users/abenson/the_trees_0_1000000_0_results.hdf5'               
                    f = hdf5.File(path, "r")
                    print(path)
                    if path.find('pchw18')!=-1:
                        stats=(path[path.find('pchw18')+7:len(path)]).ljust(9)+'\t'
                    else:
                        stats=(path[path.find('_files')+7:len(path)]).ljust(9)+'\t'
                    #stats=path+'\t'
    
            
    
    
                    for k in ['Current_a','Current_z', 'dt', 'dx', 'n_fields','n_step','offset','t','dims','dims_local','bounds','n_procs']:
                        print('k:', k, end='')
                        try:
                            print(f.attrs[k][0])
                            new_stats= str(f.attrs[k][0])+'\t'+str(f.attrs[k][1])+'\t'+str(f.attrs[k][2])+'\t'
                            new_stats= str(f.attrs[k])+'\t'
                        except:
                            try:
                                if k.startswith('n_') or k.startswith('dim'):
                                    new_stats= str(f.attrs[k][0])+'\t'
                                else:
                                    new_stats= str(format(f.attrs[k][0],'0.5f'))+'\t'
                            except:
                                new_stats='\t'
                               
                        stats+=new_stats
                    
                    stats=stats[:-1]
        
                    myOutput.writeIntoFile(
                               output_filename_grid,
                               stats+'\n',
                               append_mytext=True,
                               data_is_string=True,
                               data_format='%s')                
                    
                    i+=1
    
                exit()


            def scan_file_format():
                
                i=0
                while i<file_struc_info[catname+'_nr_files']:
                    path = file_struc_info[catname+'_filename'+str(i)]    

                    f = hdf5.File(path, "r")
                    
                    if i==0:
                        print('\n1) Attributes\n-------------------------')
                        for k in range(0,len(f.attrs.keys()),1):
                            print('\t', f.attrs.keys()[k], ':\t\t', f.attrs.values()[k][0])
                            

                    print('\n+++++++++++++++++++++\nSCANNING FILE FORMAT of ...', path,'\n')    
                    print('\n2) Data\n-------------------------')                
                    for name in f:
                        print('name:\t', name, end='')
                        try:
                            print('size:\t', f[name].size, end='') 
                            print(f[name], end='')
                            try:
                                print('min/max:\t', format(min(f[name][:]), '0.2f'), '/', format(max(f[name][:]), '0.2f'))
                            except:
                                for k in [0,1,2]:
                                    for g in [0,1,2]:
                                        #print('k:', k, 'g:', g,
                                        #print('\tdim: [:,'+str(k)+'][:,0]', format(min(f[name][:,k][:,0]), '0.5f'), '/', format(max(f[name][:,k][:,0]), '0.5f')
                                        #print('\tdim: [:,0][:,'+str(g)+']', format(min(f[name][:,0][:,g]), '0.5f'), '/', format(max(f[name][:,0][:,g]), '0.5f')
                                        print('\tsize:', f[name][:,k][:,g].size, '\tdim: [:,'+str(k)+'][:,'+str(g)+']', format(min(f[name][:,k][:,g]), '0.2f'), '/', format(max(f[name][:,k][:,g]), '0.2f'))
                                    
                        except:
                            print('--> failed!')
                    print('\n')
                    i+=1
                
                exit()


            print('###################################################################################')
            print('read CHOLLA HDF5')
            print('HDF5', 'catname:', catname)
            print('HDF5 begin: mypath:', mypath)


            file_struc_info= hdf5lib.CHOLLA_50Mpc_HDF5_filestruct(catname, snapid, path_to_directory=mypath)

            mytypes={}
            a=0
            while a<int(id_col_array['nr_entries']):
                mytypes.update({id_col_array['name'+str(a)]: id_col_array['data_type'+str(a)]})                    
                a+=1
 

            dt = np.dtype([(k, mytypes[k]) for k in mytypes.keys()])
            try:
                self.data_array=np.zeros((nr_rows,), dtype=dt)
                nr_entries=int(id_col_array['nr_entries'])
            except:
                self.data_array=np.zeros((nr_rows,), dtype=dt)
                nr_entries=nr_col
                           
            i=start_fileID*nr_files_snapshot

            if end_fileID!='False':
                i_break = start_fileID*nr_files_snapshot+(end_fileID-start_fileID+1)*nr_files_snapshot
            else:
                i_break='False'

            print('start_fileID:', start_fileID, 'end_fileID:', end_fileID, 'start i:', i,'i break:', i_break, 'nr files2read:', file_struc_info[catname+'_nr_files'])


            if catname.find('512')!=-1:
                res=512
            else:
                res=256



            ####################################################################################
            #   SCAN AND ANALYSE FILE STRUCTURE ONLY 

            #scan_file_format()

            #read_and_write_attributes(file_struc_info, str(res))

            ####################################################################################
            #   READ DATA CHOLLA FROM BRUNO ADJUSTED 
            
            
            #​dataDir = '/raid/bruno/data/'
            dataDir = '/data/groups/comp-astro/bruno/'
            dataDir = mycomp
            #inDir = dataDir + 'cosmo_sims/512_hydro_50Mpc/output_files_pchw18/'
            inDir = '/data/256_hydro_50Mpc/'
            #nDir = '/data/groups/comp-astro/bruno/cosmo_sims/512_hydro_50Mpc/output_files_pchw18/'
            
            n_snapshot = 72 #file_struc_info[catname+'_filename'+str(i)]
            
            data_type = 'hydro'
            data_type = 'particles'
  
            colname_map = { 'pos_x':'x_pos',
                            'pos_y':'y_pos',
                            'pos_z':'z_pos',
                            'vel_x':'x_vel',
                            'vel_y':'y_vel',
                            'vel_z':'z_vel'}

            fields=[]
            for col in colname_map:
                fields+= [col]
                
            #print(fields
            #fields=['density']
            precision = np.float32
            Lbox = 50000		#kpc/h
            #proc_grid = [ 4, 2, 2]
            proc_grid = [ 2, 2, 2]
            box_size = [ Lbox, Lbox, Lbox ]
            grid_size = [ res, res, res ] #Size of the simulation grid
            subgrid = [ [0, res], [0, res], [0, res] ] #Size of the volume to load
            data = load_snapshot_data_distributed(n_snapshot, inDir, data_type, fields, subgrid, precision, proc_grid,	box_size, grid_size, show_progess=True)

            print(data)

            for k in data[data_type].keys():
                print('k:', end='')
                try:
                    print(k, data[data_type][k].size)
                    self.data_array[colname_map[k]][0:data[data_type][k].size]=data[data_type][k]
                except:
                    print('--> wrong format for structured array!')
            
            self.data_array = self.data_array[:data[data_type][k].size]
            #print(np.info(self.data_array)
            #print(self.data_array[0:3]

            print('self.data_array after read-in:', self.data_array.shape, 'redshift:', self.redshift, 'scale factor:', self.scale_factor)

            #exit()

        def readHydroHDF5(key):      

            print('###################################################################################')
            print('read Hydro HDF5\n')
            print('HDF5', 'catname:', catname, 'mypath:', mypath)
            
            def caseSwitcher(item):

                choose = {
                    'EAGLE_': EAGLE,
                    'EAGLE_ev': EAGLE_ev,
                    'IllustrisTNG300': IllustrisTNG300
                    }
                    
                func = choose.get(item)
                return func()


            def EAGLE():
                file_struc_info= hdf5lib.EAGLE_HDF5_filestruct(catname, snapid, path_to_directory=mypath)
                return file_struc_info
 
            def EAGLE_ev():
                file_struc_info= hdf5lib.EAGLE_ev_HDF5_filestruct(catname, snapid, path_to_directory=mypath)
                return file_struc_info
            
            def IllustrisTNG300():
                file_struc_info= hdf5lib.IllustrisTNG300_HDF5_Subhalo_filestruct(catname, snapid, path_to_directory=mypath)
                return file_struc_info             

            from cosmolopy import cparam, cd, cc

            fidcosmo = cparam.PlanckEAGLE(flat=True, extras=False)

            file_struc_info, snapidzred = caseSwitcher(key)
            
            mytypes={}
            a=0
            while a<int(id_col_array['nr_entries']):
                mytypes.update({id_col_array['name'+str(a)]: id_col_array['data_type'+str(a)]})                    
                a+=1
 

            dt = np.dtype([(k, mytypes[k]) for k in mytypes.keys()])
            try:
                self.data_array=np.zeros((nr_rows,), dtype=dt)
                nr_entries=int(id_col_array['nr_entries'])
            except:
                self.data_array=np.zeros((nr_rows,), dtype=dt)
                nr_entries=nr_col   
                
            f = hdf5.File(file_struc_info[catname+'_filename0'], "r")

            print(f.values())
            
            for item in f.keys():
                print('env:', item, f[item].values()[0])


            snapidzred.sort(order=['snapid'], axis=0)
            
            rowNrs = [3,65,209,305,101,400]
            #env = 'WallGals'
            
            for env in f.keys():
                for rowNr in rowNrs:
                    print(env,'\t', rowNr, '\t', f[env+'/Group'][rowNr],'\t', format(f[env+'/M200'][rowNr], '0.6e'), format(f[env+'/Mstar30kpc'][rowNr], '0.6e'), format(f[env+'/TimeGalaxy50'][rowNr], '0.4f'), format(f[env+'/MBH'][rowNr], '0.6e'))
            
            #exit()
            
            # print('haloid:', f[env+'/Group'][rowNr], f[env+'/Subgroup'][rowNr], 'rVoid:', f[env+'/Rvoid'][rowNr], 'r2rVoid:', f[env+'/DtoVoid'][rowNr]
            # print('t50_stars:', f[env+'/TimeGalaxy50'][rowNr], 't70_stars:', f[env+'/TimeGalaxy70'][rowNr]
            # print('mstar_30kpc(z0):', f[env+'/Mstar30kpc'][rowNr], 'evol:', f[env+'/Mstar_ev'][rowNr]
            # print('m200c(z0):', f[env+'/M200'][rowNr], 'evol:', f[env+'/M200_ev'][rowNr]
            # print('mbh(z0):', f[env+'/MBH'][rowNr], 'evol:', f[env+'/BHmass_ev'][rowNr]
            # print('mgas_30kpc(z0):', f[env+'/Mgas30kpc'][rowNr], 'evol:',f[env+'/Mgas_ev'][rowNr] 
            
            # from scipy import interpolate
            # func = interpolate.interp1d(f[env+'/Mstar_ev'][rowNr]/f[env+'/Mstar30kpc'][rowNr], snapidzred['z'])
            
            # print('\t z 50%/70%:', format(func(0.5), '0.2f'), '/', format(func(0.7), '0.2f')
            # print('\t LBT t50/t70:', cd.lookback_time(func(0.5), z0=0.0, **fidcosmo)/3.1536e+16, '/', cd.lookback_time(func(0.7), z0=0.0, **fidcosmo)/3.1536e+16 #seconds to Gyrs     

            #print(snapidzred
            exit()
   
            #print(snapidzred['z']
            
            last_valid_item='haloid'
            size_before=0
            l=0
            for a, env in enumerate(file_struc_info[catname+'_myenv'].values()):
                
                for item in sorted(self.data_array.dtype.names, key=str.lower):

                    if item=='envr':
                        
                        print('SET ENVR FLAG:', a, 'env:', env, 'last_valid_item:', last_valid_item)
                        self.data_array[item][size_before:size_before+f[file_struc_info[catname+'_'+env+'_'+last_valid_item]][:].size] = a        
                 
                    else:
                        print('env:', env, 'a:', a, 'item:', item, file_struc_info[catname+'_'+env+'_'+item], 'size_before:', size_before, 'size:', f[file_struc_info[catname+'_'+env+'_'+item]][:].size) 
                        #print(f[file_struc_info[catname+'_'+env+'_'+item]][:]
                        self.data_array[item][size_before:size_before+f[file_struc_info[catname+'_'+env+'_'+item]][:].size] = f[file_struc_info[catname+'_'+env+'_'+item]][:]
                        last_valid_item=item                     

                size_before+=f[file_struc_info[catname+'_'+env+'_'+last_valid_item]][:].size
                
            
            self.data_array = self.data_array[:size_before]
            #print(np.info(self.data_array)
            #self.data_array.sort(order=['haloid'], axis=0)
            #print(self.data_array[['haloid', 'envr', 'mhalo_200c', 'mstar_30kpc', 't50_stars', 't70_stars', 'sfr_30kpc']][0:100]
                  
            
            self.redshift=False
            self.scale_factor=False
            
            #print(np.info(self.data_array)
            #exit()
            print('self.data_array after read-in:', self.data_array.shape, 'redshift:', self.redshift, 'scale factor:', self.scale_factor)            


        def readSAMHDF5(halocat_code):
            print('###################################################################################')
            print('read SAMHDF5')
            print('HDF5', 'catname:', catname,'HDF5 begin: mypath:', mypath)

            def caseSwitcher(catname):

                choose = {
                    'Galacticus_': Galacticus,
                    'SAG_': SAG,
                    'SAGE_': SAGE,
                    'LGALAXIES_': LGALAXIES,
                    'sussing_': ROCKSTAR
                    }
                    
                func = choose.get(catname)
                return func()

            def Galacticus():
                if halocat_code=='False' or halocat_code==False: 
                    file_struc_info= hdf5lib.Galcticus_HDF5_filestruct(catname, snapid, path_to_directory=mypath)
                else:
                    file_struc_info= hdf5lib.Galacticus_HDF5_halocat_filestruct(catname, path_to_directory=mypath)
                return file_struc_info
                
            def SAG():
                if halocat_code=='False' or halocat_code==False:
                    file_struc_info= hdf5lib.SAG_HDF5_filestruct(catname, snapid, path_to_directory=mypath)
                else:
                    file_struc_info= hdf5lib.SAG_HDF5_halocat_filestruct(catname, path_to_directory=mypath)
                return file_struc_info
            
            def SAGE():
                file_struc_info= hdf5lib.SAGE_HDF5_filestruct(catname, snapid, path_to_directory=mypath)
                return file_struc_info
            
            def LGALAXIES():
                if halocat_code=='False' or halocat_code==False: 
                    file_struc_info= hdf5lib.LGALAXIES_HDF5_filestruct(catname, snapid, path_to_directory=mypath)
                else:
                    file_struc_info= hdf5lib.LGALAXIES_HDF5_Tree_filestruct(catname, snapid, path_to_directory=mypath)
                return file_struc_info

            def ROCKSTAR():
                file_struc_info= hdf5lib.ROCKSTAR_HDF5_halocat_filestruct(catname, snapid, path_to_directory=mypath)
                return file_struc_info
            
                
            def other():
                file_struc_info= hdf5lib.default_HDF5_filestruct(catname, snapid)
                return file_struc_info

            def handle_multirow_dataset(nr_cols, a, size_before, new_size):
                print(' ... loading handle_multirow_dataset ... ', end='')
                try:
                    redshift = mL.expfactor_to_redshift(f[file_struc_info[file_struc_info['catname']+'_path_to_scale_factor']][:])
                except:
                    print('nothing to do ...')
                b=0
                while b<nr_cols:
                    print('b:', b, 'a:', a, 'name:', myid_col_array['name'+str(a)], 'nr_cols:', nr_cols)                               
                    if catname.startswith('LGALAXIES')!=-1:
                        new_size=size_before+f[file_struc_info[file_struc_info['catname']+'_index']][:,0].size
                        print('size_before', size_before, 'new_size:', new_size)
                        
                        if myid_col_array['name'+str(a)]=='haloid':
                            print('multirow haloid!')
                            self.data_array['haloid'][size_before:new_size]=f[file_struc_info[file_struc_info['catname']+'_index']][:,b]
                            print(self.data_array['haloid'][size_before:new_size].shape, 'min:', min(self.data_array['haloid'][size_before:new_size]), 'max:', max(self.data_array['haloid'][size_before:new_size]), 'nr of gal with id==-99,', len(np.where(self.data_array['haloid'][size_before:new_size]==-99)[:][0]))
                            print(self.data_array['haloid'][size_before:new_size][(np.where(self.data_array['haloid'][size_before:new_size]!=-99)[:][0])], 'length of list:', len(np.where(self.data_array['haloid'][size_before:new_size]!=-99)[:][0]))
                        elif myid_col_array['name'+str(a)]=='Z':  
                            print('multirow Z!')                              
                            self.data_array['Z'][size_before:new_size]=redshift[b]
                        else:
                            print('multirow default!')
                            self.data_array[myid_col_array['name'+str(a)]][size_before:new_size] = f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]][:,b]
                        size_before=new_size
                    else:
                        print('default!')
                        self.data_array[myid_col_array['name'+str(a)]][size_before:new_size] = f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]][:,b]
                    b+=1

            try:
                file_struc_info = caseSwitcher(catname[0:catname.find('_')+1])
            except:
                file_struc_info = caseSwitcher(catname)

            size = 0
            size_before=0
            new_size=0

            mytypes={}
            a=0
            while a<int(id_col_array['nr_entries']):
                mytypes.update({id_col_array['name'+str(a)]: id_col_array['data_type'+str(a)]})                    
                a+=1
 

            dt = np.dtype([(k, mytypes[k]) for k in mytypes.keys()])

            myid_col_array = id_col_array
            if halocat_code==False: 
                self.data_array=np.zeros((nr_rows,), dtype=dt)
            else:
                self.data_array=np.zeros((nr_rows*2,), dtype=dt)
                           
            i=start_fileID*nr_files_snapshot

            if end_fileID!='False':
                i_break = start_fileID*nr_files_snapshot+(end_fileID-start_fileID+1)*nr_files_snapshot
            else:
                i_break='False'

            print('start_fileID:', start_fileID, 'end_fileID:', end_fileID, 'start i:', i,'i break:', i_break, 'nr files2read:', file_struc_info[catname+'_nr_files'])

            while i<int(file_struc_info[catname+'_nr_files']):           
                                
                if i==i_break and end_fileID!='False':
                   break
                path = file_struc_info[catname+'_filename'+str(i)]
                #path='/store/multidark/NewMD_3840_Planck1/Galacticus/latest/job0/the_trees_0_1000000_0_results.hdf5'
                #path='/data3/users/abenson/the_trees_0_1000000_0_results.hdf5'
                f = hdf5.File(path, "r")
                #print('i:', i, 'filename:', path)

                #List snapshots for SMDPL-Galacticus 400Mpc
#                l=1
#                while l<=96:
#                    try:
#                        print('snapshot: ', l, '\ta:', format(f['/Outputs/Output'+str(l)].attrs.get('outputExpansionFactor'), '.4f'), '\tz:', format(mL.expfactor_to_redshift(f['/Outputs/Output'+str(l)].attrs.get('outputExpansionFactor')), '.2f')
#                    except:
#                        print('failed'
#                    l+=1
#                exit()
            
#                p=0
#                while p<(f['Arr/GalaxyID'][0].size):
#                    print(f['Arr/GalaxyID'][0][:]
#                    p+=1
#                exit()
#                print(f.keys()
#                for key in f.keys():
#                    print(key, '\n-----------------\n'
#
#
#
#                exit()

#                SAG                
#                 for name in f:
#                     count_att=0
#                     print('name:', name.ljust(20), 
#                     for att in f[name].attrs:
#                         print('description:', f[name].attrs.get('Description')
#                         count_att=1
#                     if count_att==0:
#                         count_att2=0
#                         print('\n',
#                         for key in f[name].keys():
#                             print('        ', key.ljust(17), 
#                             for att in f[name+'/'+key].attrs:
#                                 print'description:', f[name+'/'+key].attrs.get('Description')
#                                 count_att2=1
#                             if count_att2==0:
#                                 print('\n',
#                                 for key2 in f[name+'/'+key].keys():
#                                     print('           ', key2.ljust(14), 'description:', f[name+'/'+key+'/'+key2].attrs.get('Description')
# #                                    
#                 exit()                                    
                #Galacticus
                # path='Outputs/Output79/nodeData/'
                # print(f.keys()
                # for key in f.keys():
                #     print(key, '\n-----------------\n'
                #     for keykey in f[key].keys():
                #         print(keykey, ':\t', f[key].attrs.get(keykey)
                #         #path='Outputs/Output96/nodeData/'
                #     print('++++++++++++++++++++\n'
                    

# #
#                 for name in f[path]:
#                     print('name:', name.ljust(20)

#                 exit()                      
                            
              
#                 path='Parameters/'

                    
#                 for name in f[path].attrs.items():            
#                     print(name[0].ljust(60), name[1]

#                 for name in f[path].keys():            
#                     print(name
#                     for key in f[path+'/'+name].attrs.items():
#                         print('\t', key[0].ljust(20), key[1]                                             

#                 exit()
                #print(f['Outputs/Output79/nodeData/spheroidLuminositiesStellar:SDSS_i:observed:z0.0000:dustAtlas'])    
                    
                if halocat_code=='False' and file_struc_info['catname'].find('SAGE')==-1:
                    if file_struc_info['catname'].find('SAGE')!=-1:
                        self.scale_factor = None
                        self.redshift = snapid
                    elif file_struc_info['catname'].find('SAG_')==-1 and file_struc_info['catname'].find('LGALAXIES')==-1 and file_struc_info['catname'].find('Illustris')==-1 and file_struc_info['catname'].find('EAGLE')==-1:
                        self.scale_factor = f[file_struc_info[catname+'_path_to_snapid']].attrs.get(file_struc_info[catname+'_redshift_attribute'])
                        self.redshift = mL.expfactor_to_redshift(self.scale_factor)                     
                    else:
                        #print('else:'
                        if myid_col_array['name'+str(0)]!='Z':
                            try:
                                self.redshift = f.attrs.get(file_struc_info[catname+'_redshift_attribute'])[0]
                            except:
                                self.redshift=None
                        else:
                            try:
                                self.redshift = f.attrs.get(file_struc_info[catname+'_redshift_attribute'])[1]
                            except:
                                print('no redshift attribute found')
                        try:
                            self.scale_factor = mL.redshift_to_expfactor(self.redshift)
                        except:
                            self.scale_factor=None
                    
                    #print('self.redshift:', self.redshift, str(format(self.redshift, '.4f'))

                    #set "redshift" as part of the filename to read in Galacticus Luminosities
                    try:
                        redshift = str(format(self.redshift, '.4f'))
                    except:
                        redshift = False

                    #print('snapdid', snapid, 'self.i:', file_count, 'z:', self.redshift, 'scale_factor:', self.scale_factor)
                else:
                    if file_struc_info['catname'].find('LG')!=-1:
                        self.scale_factor=f[file_struc_info[file_struc_info['catname']+'_path_to_scale_factor']][62-(62-45)-file_count]
                        self.redshift=mL.expfactor_to_redshift(self.scale_factor)
                        print('file_count:', file_count, 'a:', self.scale_factor,'z:', self.redshift)
   
                    else:
                        self.redshift=False
                        self.scale_factor=False


                size_before=new_size

                if str(halocat_code)!='False':
                    size=f[file_struc_info[file_struc_info['catname']+'_'+myid_col_array['name'+str(1)]]][:,0].size
                    first_prop=myid_col_array['name'+str(0)]
                elif myid_col_array['name'+str(0)]!='Z':
                    try:
                        size=f[file_struc_info[file_struc_info['catname']+'_'+myid_col_array['name'+str(0)]]].size     
                        first_prop=myid_col_array['name'+str(0)]
                    except:
                        try:
                            size=f[file_struc_info[file_struc_info['catname']+'_'+myid_col_array['name'+str(1)]]].size
                            first_prop=myid_col_array['name'+str(1)]
                        except:
                            print('Des homma net!\nTake another property as the first to chose from the file, the first and second set are both not exiting!')
                            exit()
                else:
                    size=f[file_struc_info[file_struc_info['catname']+'_'+myid_col_array['name'+str(1)]]].size             
                
                new_size=size+size_before

                if size_before==0 and np.sum(f[file_struc_info[file_struc_info['catname']+'_'+str(first_prop)]].shape)>size and file_struc_info['catname'].startswith('LGALAXIES')==False:
                    print('expand array!')
                    self.data_array = np.expand_dims(self.data_array, axis=1)
                    
                
                #print('here:', str(myid_col_array['name'+str(0)]), np.sum(f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(0)])]].shape), size, np.sum(f[file_struc_info[file_struc_info['catname']+'_'+str(first_prop)]].shape)     
                check_size=size
                
                #print('size_before:', size_before, 'datasize:', size,  'size_new:', size_before+size, 'check_size:', new_size            

                a=0
                while a<int(myid_col_array['nr_entries']):
                    #print('a:', a, myid_col_array['name'+str(a)], 'col_id:', myid_col_array['col_id'+str(a)]#, np.prod(f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]].shape)
                    if myid_col_array['name'+str(a)].startswith('L_'):
                        #print('here 1337', myid_col_array['name'+str(a)], end=' ')
                        file_struc_info[file_struc_info['catname']+'_'+myid_col_array['name'+str(a)]] = file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])+'_part1']+redshift+file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])+'_part2']                   
                        #print('Lum:', file_struc_info[file_struc_info['catname']+'_'+myid_col_array['name'+str(a)]])
                    

                    if check_size<np.sum(f[file_struc_info[file_struc_info['catname']+'_'+str(first_prop)]].shape) and i==0 and catname.startswith('LG')==-1:
                        self.data_array = np.expand_dims(self.data_array, axis=1)
                        #print('new shape:', self.data_array.shape, 'check_size:', check_size, np.sum(f[file_struc_info[file_struc_info['catname']+'_mstarsph']].shape)
                        check_size+=1
                        
                    if myid_col_array['name'+str(a)] == 'ngalaxies':
                        #print('ngalaxies at i:', i, f[file_struc_info[file_struc_info['catname']+'_ngalaxies']]
                        self.data_array['name'+str(a)][size_before:new_size] = size

                    elif myid_col_array['name'+str(a)]=='x_pos' and (catname.startswith('LGALAXIES') or catname.startswith('Illustris')):
                        print('LGALAXIES/Illustris positions ...')
                        self.data_array[myid_col_array['name'+str(a)]][size_before:new_size]=f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]][:,0] 
                        self.data_array['y_pos'][size_before:new_size]=f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]][:,1]                       
                        self.data_array['z_pos'][size_before:new_size]=f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]][:,2]
                        
                        #print(self.data_array[myid_col_array['name'+str(a)]][0:10]
                        #print(self.data_array['y_pos'][0:10]
                        #print(self.data_array['z_pos'][0:10]
 
                    elif myid_col_array['name'+str(a)]=='x_vel' and catname.startswith('LGALAXIES'):
                        print('LGALAXIES velocities ...')
                        self.data_array[myid_col_array['name'+str(a)]][size_before:new_size]=f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]][:,0] 
                        self.data_array['y_vel'][size_before:new_size]=f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]][:,1]                       
                        self.data_array['z_vel'][size_before:new_size]=f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]][:,2]                       
 
                    elif catname.startswith('LGALAXIES'):
                        self.data_array[myid_col_array['name'+str(a)]][size_before:new_size] = f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]][size_before:new_size,file_count]

                    elif myid_col_array['name'+str(a)]=='Z' and catname.startswith('SAGE')==False:
                        #print('manage redshift!', myid_col_array['name'+str(a)], int(myid_col_array['col_id'+str(a)])
                        self.data_array[myid_col_array['name'+str(a)]][size_before:new_size]=format(float(self.redshift), '.2f')                                
                    else:
                                                                                
                        if str(myid_col_array['name'+str(a)])=='mstar' or str(myid_col_array['name'+str(a)])=='sfr' or str(myid_col_array['name'+str(a)])=='mcold' or str(myid_col_array['name'+str(a)])=='Mzgas' or str(myid_col_array['name'+str(a)])=='Mzstar':
                            
                            if myid_col_array['name'+str(a)].startswith('Mz') and catname.startswith('Galacticus'):
                                pass
                            elif catname=='SAG_1Gpc' or catname.find('run2')!=-1 or catname.find('v2')!=-1 or catname.find('Galacticus')!=-1: 
                                self.data_array[myid_col_array['name'+str(a)]][size_before:new_size] = f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])+'_disk']][:] + f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])+'_spheroid']][:]
                                #print(str(myid_col_array['name'+str(a)]), 'here disk+sph!'
                            else:
                                #print('else'
                                #print('size_before', size_before, 'new_size:', new_size, file_struc_info['catname'], file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]
                                self.data_array[myid_col_array['name'+str(a)]][size_before:new_size] = f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]]

                        elif str(myid_col_array['name'+str(a)]).startswith('MA') or str(myid_col_array['name'+str(a)]).startswith('mA') or str(myid_col_array['name'+str(a)]).startswith('mag') or str(myid_col_array['name'+str(a)]).startswith('Mag'):
                            #print('mag!',
                            if catname.startswith('Galacticus'):# or catname.startswith('LGALAXIES'):                             
                                pass
                            elif (catname.startswith('SAG_') or catname.startswith('LGALAXIES')) and str(myid_col_array['name'+str(a)]).startswith('MA'):
                                #print('MAB'
                                self.data_array[myid_col_array['name'+str(a)]][size_before:new_size] = f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]]
                                #print(min(self.data_array[myid_col_array['name'+str(a)]][size_before:new_size]), max(self.data_array[myid_col_array['name'+str(a)]][size_before:new_size])                        
                            
                            elif catname.startswith('SAGE'):
                                #print('here:!'
                                self.data_array[myid_col_array['name'+str(a)]][size_before:new_size] = f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]]                                
                        #elif str(myid_col_array['name'+str(a)]).startswith('L_') and str(myid_col_array['name'+str(a)]).startswith('isk', len(myid_col_array['name'+str(a)])-5, len(myid_col_array['name'+str(a)])-2)==False and str(myid_col_array['name'+str(a)]).startswith('oid', len(myid_col_array['name'+str(a)])-5, len(myid_col_array['name'+str(a)])-2)==False and str(myid_col_array['name'+str(a)]).startswith('tal', len(myid_col_array['name'+str(a)])-5, len(myid_col_array['name'+str(a)])-2)==True:
                        elif str(myid_col_array['name'+str(a)]).startswith('L_') and myid_col_array['name'+str(a)].find('total')!=-1:
                                filter_name = myid_col_array['name'+str(a)][len(myid_col_array['name'+str(a)])-myid_col_array['name'+str(a)][::-1].find('_')::]
                                try:
                                    self.data_array[myid_col_array['name'+str(a)]][size_before:new_size] = f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]][:]
                                except:
                                    #print('no total Luminosity found! try ...'
                                    try:
                                        #print('Lum sph+disk!', myid_col_array['name'+str(a)], 'filter_name:', filter_name, file_struc_info[file_struc_info['catname']+'_'+myid_col_array['name'+str(a)][0:myid_col_array['name'+str(a)].find('total')]+'disk_'+filter_name], '+', file_struc_info[file_struc_info['catname']+'_'+myid_col_array['name'+str(a)][0:myid_col_array['name'+str(a)].find('total')]+'spheroid_'+filter_name],
                                        #print('Lum sph+disk!', myid_col_array['name'+str(a)], 'filter_name:', filter_name#, file_struc_info['catname']+'_'+myid_col_array['name'+str(a)][0:myid_col_array['name'+str(a)].find('total')]+'disk_'+filter_name, '+', file_struc_info['catname']+'_'+myid_col_array['name'+str(a)][0:myid_col_array['name'+str(a)].find('total')]+'spheroid_'+filter_name,
                                        #print(file_struc_info['catname']+'_'+myid_col_array['name'+str(a)][0:myid_col_array['name'+str(a)].find('total')]+'disk_'+filter_name
                                        #print(f[file_struc_info[file_struc_info['catname']+'_'+myid_col_array['name'+str(a)][0:myid_col_array['name'+str(a)].find('total')]+'disk_'+filter_name]][0:10]
                                        #print(file_struc_info['catname']+'_'+myid_col_array['name'+str(a)][0:myid_col_array['name'+str(a)].find('total')]+'spheroid_'+filter_name
                                        #print(f[file_struc_info[file_struc_info['catname']+'_'+myid_col_array['name'+str(a)][0:myid_col_array['name'+str(a)].find('total')]+'disk_'+filter_name]][0:10], '+', f[file_struc_info[file_struc_info['catname']+'_'+myid_col_array['name'+str(a)][0:myid_col_array['name'+str(a)].find('total')]+'spheroid_'+filter_name]][0:10]

                                        self.data_array[myid_col_array['name'+str(a)]][size_before:new_size] = f[file_struc_info[file_struc_info['catname']+'_'+myid_col_array['name'+str(a)][0:myid_col_array['name'+str(a)].find('total')]+'disk_'+filter_name]][:] + f[file_struc_info[file_struc_info['catname']+'_'+myid_col_array['name'+str(a)][0:myid_col_array['name'+str(a)].find('total')]+'spheroid_'+filter_name]][:]
                                        #print('-->CHECK!'
                                    except:
                                        pass
                                        #print('i:', i, 'filename:', path, 'z=', self.redshift) 
                                        #print('no spheroid or disk galaxy parts found ...!\na:', a, myid_col_array['name'+str(a)], 'col_id:', myid_col_array['col_id'+str(a)])

                            
                        try:
                            if np.prod(f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]].shape)>np.sum(f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]].shape):
                                nr_cols = int(np.prod(f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]].shape)/f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]][:,0].size)
                                handle_multirow_dataset(nr_cols, a, size_before, new_size)
                        except:
                            if np.prod(f[file_struc_info[file_struc_info['catname']+'_'+str(first_prop)]].shape)>np.sum(f[file_struc_info[file_struc_info['catname']+'_'+str(first_prop)]].shape):
                                nr_cols = int(np.prod(f[file_struc_info[file_struc_info['catname']+'_'+str(first_prop)]].shape)/f[file_struc_info[file_struc_info['catname']+'_'+str(first_prop)]][:,0].size)
                                handle_multirow_dataset(nr_cols, a, size_before, new_size)
                            
                        else:
                            #print('name:', myid_col_array['name'+str(a)], 'size_before:', size_before, 'new_size:', new_size, 'datasize:', file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]
                            self.data_array[myid_col_array['name'+str(a)]][size_before:new_size] = f[file_struc_info[file_struc_info['catname']+'_'+str(myid_col_array['name'+str(a)])]]

                    a+=1
                i+=1

            #self.data_array = np.resize(self.data_array,(new_size, int(myid_col_array['nr_cols2read'])))
            self.data_array = self.data_array[:new_size]
            if catname.startswith('SAG_'): self.data_array = np.reshape(self.data_array, (len(self.data_array),))

            print('self.data_array after read-in:', self.data_array.shape)


        def caseSwitcher(data_format):

            choose = {
                'ASCII': ASCII,
                'ASCII_SPLITTED_CAT': ASCII,
                'CATASCII': ASCII,
                'BINARY': BINARY,
                'BINARY_SAGE': BINARY_SAGE,
                'HDF5': HDF5,
                'SAMHDF5': SAMHDF5,
                'HYDROHDF5': HYDROHDF5,
                'CHOLLAHDF5': CHOLLAHDF5,
                'MDHDF5': MDGalaxies, 
                'HDF52HDF5': HDF52HDF5,
                'READ2ARRAY': READ2ARRAY,
                'FITS2HDF5': FITS2HDF5,
                'ROCKSTAR_ASCII': ROCKSTAR_ASCII,
                'EAGLE_ASCII_Silvio': EAGLE_ASCII,
                'EAGLE_ASCII_Patricia1': EAGLE_ASCII,
                'EAGLE_ASCII_Patricia2': EAGLE_ASCII,
                'EAGLE_ASCII_Patricia_rgas': EAGLE_ASCII,
                'EAGLE_ASCII_Yetli': EAGLE_ASCII,
                'EAGLE_ASCII_full': EAGLE_ASCII,
                'EAGLE_ASCII_full_trees': EAGLE_ASCII,
                'EAGLE_ASCII': EAGLE_ASCII,
                'TTH_ASCII': TTH_ASCII
                }
                
            func = choose.get(data_format)
            return func()
            
        def ASCII():
            readASCII()          
    
        def BINARY():
            readBINARY()

        def BINARY_SAGE():
            readBINARY_SAGE(mypath) 
           
        def SAMHDF5():
            readSAMHDF5(halocat_code)
 
        def HYDROHDF5():
            key = catname[0:catname.find('ev')+2]
            print('key', key)  
            readHydroHDF5(key)
            
        def MDGalaxies():
            read_MDGalaxies()            

        def EAGLE_ASCII():
            key = data_format[data_format.find('ASCII')+6::]
            print('key', key)      
            readEAGLE_ASCII(key)
            
        def TTH_ASCII():    
            readASCII_TTH()

        def ROCKSTAR_ASCII():
            readROCKSTAR_ASCII(halocat_code) 
            
        def CHOLLAHDF5():
            readCHOLLAHDF5(halocat_code)            
        
        def HDF5():
            readHDF5()
 
        def HDF52HDF5():
            readANDConvert(data_format)
            
        def READ2ARRAY():
            readANDConvert(data_format)  

        def FITS2HDF5():
            readANDConvert(data_format)
            
        #print('mypath_softlink:', mypath_softlink
        #print('myfilename:', myfilename
        #print(data_format
        #print('mypath:', mypath
        check_path = os.path.exists(mypath)
        #print('check_path:', check_path
        #check_path_softlink = os.path.exists(mypath_softlink)
        
        if config!=False:

            if mypath!=mypath_softlink:
                mypath=mypath_softlink   
            elif data_format == 'SAMHDF5':
                mypath+=myfilename
            
            if data_format=='BINARY':
                mypath+=myfilename
                
            if data_format.find('FITS')!=-1:
                mypath+=myfilename

           
        #print('mypath:', mypath, 'myfilename:', myfilename        
                    
        if check_path==False:# and check_path_softlink==False:
            
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n', mypath, 'NOT EXISTS!\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')

        #print('config:', config)
     
        caseSwitcher(data_format)

        return self.data_array


    def selectData2Compute(
            self,
            data,
            single=False,
            selected_col=0,
            operator='>',
            condition=0.0):

        #print('data shape before selection:', data.shape
        if single==True:
            #print('single=yes'
            data = np.expand_dims(data, axis=1)
            print('new data shape: ', data.shape)
        #print('selected_col', selected_col    
        if selected_col=='*':
            #print('*', condition
            mask = np.where(data.any(axis=1) > condition)
            
        else:
            #print('condition:', 'operator:', operator, 'selected col:', selected_col, 'condition:', condition
            if operator =='>':
                #print(operator, '>'
                mask = np.where(data[selected_col] > condition)
            elif operator =='<':
                #print(operator, '<'
                mask = np.where(data[selected_col] < condition)
            elif operator =='>=':
                #print(operator, '>='
                mask = np.where(data[selected_col] >= condition)
            elif operator =='<=':
                #print(operator, '<='
                mask = np.where(data[selected_col] <= condition)
            elif operator =='==':
                #print(operator, '='
                mask = np.where(data[selected_col] == condition)
            elif operator =='!=':
                #print(operator, '='
                mask = np.where(data[selected_col] != condition)

        #self.selected_data = data[mask[:][0]]

        #print('data after selection:', data[mask[:][0]].shape

        return data[mask[:][0]]

              
 
      

