# Load packages
import config as conf
myConfig = conf.Configuration()
import scipy as scipy
import arangeData as aD
myData = aD.ArangeData()

import myLib as mL
import numpy as np
#from astropy.cosmology import Planck13 as cosmo
#import astropy.constants as const
#import cosmolopy as cosfunc

import os
system_info = os.getcwd()
start = system_info.find('anaconda')
mycomp = system_info[:start]

import outputData as oD

myOutput = oD.OutputData(config=False)

from time import time



class MyFunctions:
                        
    def normaliseHisto(
            self,
            histo,
            norm_x=1,
            norm_y=1,
            log_x=False,
            x_error=0,
            log_x_err=False,
            real_mean_x=False,
            log_y=False,
            y_error=0,
            y_error_type='Poisson',
            log_y_err=True):
                
        """
        Input:
        ------------------
        histo[:,0] = histo or binned data points x-axis
        histo[:,1] = histo or binned data points y-axis
        histo[:,2] = x-error (binsize)
        histo[:,3] = x-error (binsize)
        histo[:,4] = -dy
        histo[:,5] = +dy  
        histo[:,6] = number count of binned objects y-axis
        """
        print 'Normalisations y-Axis: NORM Y:', norm_y, 'error type:', y_error_type, 'NORM X:', norm_x
        #Choose x-axis, log or linear?

        #Normalisations x-Axis
        if log_x==True:
            histo[:,0] = np.log10(histo[:,0]/norm_x)
        else:
            histo[:,0]/=norm_x

        #Choose y-axis, log or linear?
        #Normalisations y-Axis    
        if log_y==True:
            histo[:,1] = np.log10(histo[:,1]/norm_y)  
        else:
            histo[:,1]/=norm_y

        if x_error!=0:
            #Choose x-error, log or linear?
            if log_x_err==True:
                histo[:,2] = np.log10(x_error/2.)
            else:
                histo[:,2] = x_error/2.

            histo[:,3] = histo[:,2]                
        
        #Choose error type for y-error
        if y_error_type=='Poisson':
            #y-error --> N/poisson(N**0.5)
            print 'log_y_err:', log_y_err
            if log_y_err==True:
                
                histo[:,4] = (histo[:,1]/(histo[:,6]**0.5))*0.43
            else:
                histo[:,4] = histo[:,1]/(histo[:,6]**0.5)

                
            histo[:,5] = histo[:,4]
        
        else:
            histo[:,[4,5]]/=norm_y
            
        
        return histo                                          

    def binUp(
        self,
        data,
        #names,
        nbins,
        log10bin=False,
        histo_min='min',
        histo_max='max',
        binup_2D=False,
        use_MAD=True,
        div_y2x=False,
        add_axis=False,
        norm_by_binsize=True,
        cumulative=False,
        weights=False,
        equal_bins=False):

        """
        Input:
        ------------------
        data[:,0] = data points x-axis (e.g. Mhalo)
        data[:,1] = data points y-axis input_col[1]/input_col[0] (e.g. Mstar/Mhalo)
        data[:,2] = data points x values for additional x-axis (e.g. Mstar in Mstar/Mhalo v. Mhalo plot)
        
        Output:
        ------------------
        binned_array[:,0] = results x-axis binned --> always, median (50%)
        binned_array[:,1] = results y-axis binned together (binup_2D=True) or N-count (binup_2D=False) --> always, median (50%)
        binned_array[:,2] = results x error classic standard deviation (std), or binsize (N-count), if binup_2D=True: 95%
        binned_array[:,3] = results x error classic standard deviation (std), or binsize (N-count), if binup_2D=True: 5%      
        binned_array[:,4] = results y error (dy-) std, MAD calculated here or Poission--> if 'normalise=True' calculated in normaliseHisto(), if binup_2D=True: 95%
        binned_array[:,5] = results y error (dy+) std, MAD calculated here or Poission--> if 'normalise=True' calculated in normaliseHisto(), if binup_2D=True: 5%
        binned_array[:,6] = number count of objects binned together in on bucket
        binned_array[:,7] = results x values for additional x-axis (e.g. Mstar in Mstar/Mhalo v. Mhalo plot)
        binned_array[:,8] = results x_add_axis error (dy+) std, MAD calculated here or Poission--> if 'normalise=True' calculated in normaliseHisto(), if binup_2D=True: 95%
        binned_array[:,9] = results x_add_axis error (dy-) std, MAD calculated here or Poission--> if 'normalise=True' calculated in normaliseHisto(), if binup_2D=True: 5%
        """
  
        #sort values to bins!
        def sortAndCalc():
            
            data2bin= {}
            i=0
            while (i<binned_array[:,0].size):
                #binned_array[i,0] = (data_min + binsize*(i+0.5))
                #print 'bucket example:', 'bin:', binned_array[i,0], 'Mhalo:', 9+i*0.3, (9+i*0.3-data_min)/binsize   
                key = 'x'+str(i)
                data2bin[key] = []
    
                key = 'y'+str(i)
                data2bin[key] = []
                
                key = 'x_add_axis'+str(i)
                data2bin[key] = []                           
                i+=1
               
            def calc_buckets(data1,data2,data3):
                

                data1=data1[np.where(data1>data_min)[:][0]]
                data1=data1[np.where(data1<=data_max)[:][0]]
                
                if weights!=False: 
                    data_weight=data2
                    
                else:
                    data_weight=np.ones((data1.size,), dtype=np.int8)
                                
                for i in range(len(data1)):                        
                    if log10bin==True:
                        bucket = np.floor((data1[i]-data_min)/binsize)
                    elif linbin==True:
                        bucket = np.floor((data1[i]+data_min)/binsize)
                    else:
                        bucket = np.floor(np.log10(data1[i]/data_min)/binsize)
                    
                    bucket = np.int_(bucket)
                    #print 'i:', i, 'bucket:', bucket, 'data[i]:', "{0:.10e}".format(10**data1[i]), "{0:.10e}".format(10**data_weight[i])
            
                    if data1[i]==data_max:
                        print 'data_max:', data_max,#, data[i,0]
                        bucket=nbins-1
                        print bucket
#                    elif data1[i]==data_min:
#                        print 'data_min:', data_min, data[i,0]
#                        bucket=0 
#                        print bucket
                        
            
                    if log10bin==True:
                        data1[i]=10**data1[i]
                        data2[i]=10**data2[i]
                        data3[i]=10**data3[i]
                             
                        data2bin['x'+str(bucket)].append(data1[i])
                        data2bin['y'+str(bucket)].append(data2[i])
                        data2bin['x_add_axis'+str(bucket)].append(data3[i])                              
                    
                    #bucket number count n
                    binned_array[bucket,6]+=1*data_weight[i]       

            
            #sort into the buckets    
            try:               
                calc_buckets(data[str(names[0])],data[str(names[1])],data[str(names[2])])    
            except:
                calc_buckets(data[:,0],data[:,1],data[:,2]) 

            #Calculate median/mean values of the buckets!              
            i=0            
            while i<binned_array[:,0].size:

                if names[1]=='ngal':
                        #print 'i:', i, 'n counts:', binned_array[i,6]
                        if len(data2bin['x'+str(i)])>0:
                            binned_array[i,0] = np.nanmedian(data2bin['x'+str(i)])
                            binned_array[i,2] = np.percentile(data2bin['x'+str(i)], 95)                
                        if len(data2bin['y'+str(i)])>0:
                            #print data2bin['y'+str(i)]
                            #print 'sum:', np.sum(data2bin['y'+str(i)]), '/nhalos:', np.sum(data2bin['y'+str(i)])/len(data2bin['y'+str(i)])#, '+/-', np.percentile(data2bin['y'+str(i)], 95)/len(data2bin['y'+str(i)]), '/', np.percentile(data2bin['y'+str(i)], 5)/len(data2bin['y'+str(i)])
                            binned_array[i,1] = np.sum(data2bin['y'+str(i)])/len(data2bin['y'+str(i)])
    
                            binned_array[i,3] = np.percentile(data2bin['y'+str(i)], 95)/binned_array[i,6]
                            binned_array[i,4] = np.percentile(data2bin['y'+str(i)], 5)/binned_array[i,6]                         
                        try:
                            if len(data2bin['x_add_axis'+str(i)])>0:
                                binned_array[i,7] = np.sum(data2bin['x_add_axis'+str(i)])
                                binned_array[i,8] = np.percentile(data2bin['x_add_axis'+str(i)], 95)/binned_array[i,6]
                                binned_array[i,9] = binned_array[i,8]                                
                        except:
                            pass
                else:
                    if use_MAD==True:
                        #print 'USE MEDIAN AND MAD!'
                        if len(data2bin['x'+str(i)])>0: 
                            binned_array[i,0] = np.nanmedian(data2bin['x'+str(i)])
                            #MAD is 75% percentile!
                            binned_array[i,2] = robust.mad(data2bin['x'+str(i)])
                            binned_array[i,3] = binned_array[i,2]                             
                        if len(data2bin['y'+str(i)])>0:
                            binned_array[i,1] = np.nanmedian(data2bin['y'+str(i)])
                            #MAD is 75% percentile!
                            binned_array[i,4] = robust.mad(data2bin['y'+str(i)])
                            binned_array[i,5] = binned_array[i,4]                            
                        try:
                            if len(data2bin['x_add_axis'+str(i)])>0:
                                binned_array[i,7] = np.median(data2bin['x_add_axis'+str(i)])
                                #MAD is 75% percentile!
                                binned_array[i,8] = robust.mad(data2bin['x_add_axis'+str(i)])
                                binned_array[i,9] = binned_array[i,8]                                 
                        except:
                            pass
    
                    else:
                        print 'USE MEDIAN AND PERCENTILES!'
                        for col_data,col_bin,perc in zip(['x','y','x','x','y','y','x_add_axis','x_add_axis','x_add_axis'],[0,1,2,3,4,5,7,8,9],[50,50,16,84,32,68,50,16,84]):
                            print 'col_bin:', col_bin, 'col_data:', col_data, 'perc:', perc
                            try:
                                binned_array[i,col_bin] =  np.percentile(data2bin[col_data+str(i)], perc)
                            except:
                                print 'bin:', col_bin, 'for:', col_data, 'not calculateable for ', perc                                                                                                                
                i+=1
      

        def binData():
            #binning the data if a number counted histogram is calculated (no cumulative, no bin2D)
            
            if names[1]!='ngal':         
                print 'binup_2D and cumulative: False --> normalise through!', binsize
                binned_array[:,1] = binned_array[:,6]/binsize

            if equal_bins==True:
                if log10bin==True:
                    bins=np.arange(data_min,data_max,binsize)
                    print bins
                    binned_array[:,0]=10**bins
                    
            else:
                if log10bin==True:
                    print 'log10bin: True --> log quantity , linear binning!'               
                    #print 'data_min:', data_min, 'data_max:', data_max, 'bin_size:', binsize
                    i=1
                    while i<nbins+1:
                        bucket_left = data_min + (i-1)*binsize
                        bucket_right = data_min + i*binsize
                                       
                        binned_array[i-1,0] = (10**bucket_left*10**bucket_right)**0.5
                        
                        #print 'i-1:', i-1, 'bucket_left:', bucket_left, 'bucket_right:', bucket_right, 'bucket_center:', binned_array[i-1,0], 'log10bin', (bucket_left*bucket_right)**0.5
                        i+=1                                      
                else:
                    print 'log10bin: False --> linear quantity , log binning!'                           
                    #print 'data_min:', data_min, 'data_max:', data_max, 'bin_size:', binsize                    
                    i=1
                    while i<nbins+1:
                        bucket_left = data_min*10**(binsize*(i-1))
                        bucket_right = data_min*10**(binsize*i)
                                       
                        binned_array[i-1,0] = (bucket_left+bucket_right)*0.5
                        #print 'i-1:', i-1, 'bucket_left:', bucket_left, 'bucket_right:', bucket_right, 'bucket_center:', binned_array[i-1,0], 'log10bin', (bucket_left*bucket_right)**0.5
                        i+=1 
          

        def cumHisto():
            print 'create cumulative histogramm!'
            cumhisto, limits, binsize, extrapoints=scipy.stats.cumfreq(data[:,0], defaultreallimits=(data_min, data_max), numbins=nbins)

            binned_array[:,1] = cumhisto[::-1]
            binned_array[:,2] = binsize
            binned_array[:,6] = cumhisto[::-1]
           
        
        #MAIN binup()
        #------------------------------------------------------------------
        from statsmodels import robust            
        start = time()
        #rearange input data array --> has to have this shape data.shape = (n-rows,3)

        try: 
            #try in case data is structured array!
            #print 'try!'

            names=[]        
            for f in data.dtype.names:
                names+=[f]
                
            print names
            if len(names)<3:       
                import numpy.lib.recfunctions as rcfuncs
                data = rcfuncs.append_fields([data], ['add_axis'] ,[data[str(names[0])]], usemask=False)
                names+=['add_axis']

            data_min = min(data[str(names[0])])
            data_max = max(data[str(names[0])])
      
            if div_y2x==True:
                print 'div_y2x:', div_y2x
                data[str(names[1])]/=data[str(names[0])]    
    
            if log10bin==True:
                for f in data.dtype.names:
                    data[f]=np.log10(data[f])
        
        except:
            names=[99,99]          
            #print 'except!'
            #data is not a structured array, treat as numpy ndarray!
            if div_y2x==False:
                myexpand_col = 0
            else:
                myexpand_col = 1
    
            #rearange input data array --> has to have this shape data.shape = (n-rows,3)       
            data = mL.dim_expander(data, 3, id_of_col_to_expand=myexpand_col) 
    
            if div_y2x==True:
                print 'div_y2x:', div_y2x
                data[:,1]/=data[:,0]

            if histo_min=='min':   
                data_min = min(data[:,0])
            else:
                data_min=float(histo_min)
            if histo_max=='max':                 
                data_max = max(data[:,0])
            else:
                data_max=float(histo_max)

            if log10bin==True:
                data = np.log10(data)
        

        print 'use MAD for calculating errorbars?', use_MAD, 'use log10 bining?', log10bin

        linbin=False
        if log10bin==True:
            print 'log quantity, linear binning!'
            
            if equal_bins==True:              
                                
                print 'equal bins!', histo_min
                data_min=float(histo_min)
                data_max=float(histo_max)
                                
            data_min = np.log10(data_min)
            data_max = np.log10(data_max)
            binsize = (data_max-data_min)/nbins
           
        else:
                   
            binsize = np.log10(data_max/data_min)/nbins

            if np.int_(binsize)<0:
                binsize = np.log10(data_max-data_min)/nbins
                print 'linear quantity, linear binning!'
                linbin=True
            else:                             
                print 'linear quantity, log binning!'
                
        if add_axis!=False:
            size_binned_array=10
        else:
            size_binned_array=8
            
        binned_array = np.zeros((nbins, size_binned_array), dtype=np.float64)
        
        print 'data_min:', data_min, 'data_max:', data_max, 'histo min/max:', histo_min, '/', histo_max 
                                                                                      
        print 'binsize:', binsize, 'shape of binned_array:', binned_array.shape 

        if cumulative==False:
            sortAndCalc()

        else:
            cumHisto()

        if use_MAD==True and binup_2D==False:
            binData()
            
        print 'Time 2 bin up with binUp:', (time()-start)/60.0, 'min/', (time()-start), 'sec'

#        print 'check_shape:', binned_array.shape
#        print binned_array

        return binned_array, binsize
