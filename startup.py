# Load packages
import os
system_info = os.getcwd()
start = system_info.find('anaconda')
if start==-1:
    mycomp='/home/dstoppacher/'
else:
    mycomp = system_info[:start]

#print 'mycomp:', mycomp, system_info, start
import compileall
compileall.compile_dir(mycomp+'anaconda/pro/myPy/', force=1)

import subprocess as subs

if mycomp.startswith('/homer'):
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat.sh', shell=True)
    subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_home.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_EAGLE_cross_history_assembly.sh', shell=True)
else:
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_RS.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_RS2.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_LG2.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_Gal2300.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_CMASS.sh', shell=True) 
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_Gal2.sh', shell=True)    
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_Gal400.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_TNG300.sh', shell=True)    
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_300_Gal.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_300_SAG.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_Gal_for_David.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_SAG_for_David.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_SAGE_for_David.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_300_SAGE.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_300.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_Cholla.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_EAGLE_full.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_EAGLE_trees.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_EAGLE_cross.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_EAGLE_cross_history.sh', shell=True)
    subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_EAGLE_cross_history_assembly.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_EAGLE_cross_history_fit_test.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_EAGLE_cross_reduced_props.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_EAGLE_cross_for_P.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_EAGLE_Patricia.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_EAGLE_Patricia2.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_EAGLE_Ben.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_EAGLE_ev_Yetli.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_EAGLE_Yetli.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_filter_EAGLE_Silvio.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_SFH_SDSS2CMASS.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_SFH.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_SFH_Gal400.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat_SFH_run2.sh', shell=True)    
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat2.sh', shell=True)
    #subs.call(mycomp+'anaconda/pro/myRun/gen_input_cat.sh', shell=True) 

#Check with this routine if a bad character is in a certain file!
#with open(mycomp+'anaconda/pro/myPy/myLib.py') as fp:
#    for i, line in enumerate(fp):
#        if "\xce" in line:
#            print 'Here is the nasty bastard:', i, repr(line)
#exit()

import config as conf
import doAnalysis as a

myConfig = conf.Configuration()
myConfig.SAMScaleFactorMapping()
myConfig.readMyRunConfig(myfile_format='', read_filenr_sequence=False)
myConfig.plotMappingArray()
myConfig.plotMarkerColorKewords()
myConfig.physicsSpecs()
myConfig.histoConfigs()
myConfig.filterDataConfig()
myConfig.nameConvMap()

#print myConfig.config_array
#print myConfig.snapid_array
#print myConfig.histo_config_array

myAnalysis = a.DoAnalysis(myConfig.config_array,
                         myConfig.snapid_array,
                         myConfig.plot_map_array,
                         myConfig.physics_specs,
                         myConfig.histo_config_array,
                         myConfig.mycond_config_array,
                         myConfig.SAM_scale_factor_map,
                         myConfig.name_conv_map)
                  
myAnalysis.MAIN()  
                
print 'EOF - END OF FUN'
print '-------------------------'                   

 
