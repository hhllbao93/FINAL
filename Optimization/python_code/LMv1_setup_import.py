# #########################################################
# ### Script to set up and import data for main scripts ###
# #########################################################

# ## OPTS
# ## # opts[0]: Initial (i.e., Niters==1) mu = opts[0] * max(np.diag(A))
# ## # opts[1]: Stop if ng <= opts[1]
# ## # opts[2]: Stop if nh <= opts[2]*nx
# ## # opts[3]: Max # iterations

## Import required modules

from netCDF4 import Dataset
import numpy as np
import glob
import re
import csv
import time
import shutil
from datetime import datetime as dt
import sys
import argparse
import timeit
    
import LMv1_ancillaries as ssr

# Name present working directory (add trailing slash if not present)
pwd = os.getcwd()
if pwd[len(pwd)-1] != '/':
    pwd = pwd + '/'

## Parse input arguments

def restrict_SSEpctile(x):
    x = float(x)
    if x < 0.0 or x > 100.0:
        raise argparse.ArgumentTypeError("Bad SSEpctile: %r not in range [0.0, 100.0]"%(x,))
    return x

parser = argparse.ArgumentParser(description='This parses input arguments.')
parser.add_argument('-p','--param_work_dir',type=str, help='Specify parameterization_work directory (relative path from bash script).')
parser.add_argument('-o','--output_dir',type=str, help='Specify directory with output data for all runs in this set (relative path from bash script).')
parser.add_argument('-v','--verbose', dest='verbose', type=int, choices=[0, 1, 2], help='Turns on verbose output.')
parser.add_argument('-n','--n_years',type=int,help='Use the last [1] years of output data.')
parser.add_argument('--constrain_AGBparams', dest='constrain_AGBparams', action='store_true', help='If present, constrain AGB parameters (Li2012-style only!).')
parser.add_argument('--constrain_Ia', dest='constrain_Ia', action='store_true', help='If present, constrain alpha and param1 for Ia to >=0.')
parser.add_argument('--constrain_popDsupp_NF', dest='constrain_popDsupp_NF', action='store_true', help='If present, constrain popDsupp_NF parameters.')
parser.add_argument('--constrain_RHparams', dest='constrain_RHparams', action='store_true', help='If present, constrain AGB parameters (Li2012-style only!).')
parser.add_argument('--constrain_duration', dest='constrain_duration', action='store_true', help='If present, constrain fire duration for each species to >=0.')
parser.add_argument('--constrain_ROSmax', dest='constrain_ROSmax', action='store_true', help='If present, constrain max rate of spread for each species to >=0.')
parser.add_argument('--constrain_thetaEsq', dest='constrain_thetaEsq', action='store_true', help='If present, constrain thetaEsquared to >0.')
parser.add_argument('--constrain_magicScalar', dest='constrain_magicScalar', action='store_true', help='If present, constrain magic_scalar to >0.')
parser.add_argument('--constrain_strict', dest='constrain_strict', action='store_true', help='If present, constrain Ia_param2 and RH_up to <=1.')
parser.add_argument('-d','--damping_style',type=str,choices=['additive','multiplicative','transtrum_minVal'], help='Specify damping style: ["additive"], "multiplicative", or "transtrum_minVal".' )
parser.add_argument('--f_rh_style',type=str,choices=['li2012','logistic','gompertz'],help='Specify f_rh_style: ["li2012"], "logistic", or "gompertz".')
parser.add_argument('--f_agb_style',type=str,choices=['li2012','logistic','gompertz'],help='Specify f_agb_style: ["li2012"] or "logistic", or "gompertz".')
parser.add_argument('--f_theta_style',type=str,choices=['li2012','logistic','gompertz'],help='Specify f_theta_style: ["li2012"] or "logistic", or "gompertz".')
parser.add_argument('-r','--rescale_params', dest='rescale_params', action='store_true', help='If present, adjust parameter values & derivatives to improve scaling.')
parser.add_argument('-O','--options', dest='opts', nargs=4, type=float)
parser.add_argument('--doBins',dest='doBins',action='store_true', help='If present, shift scnd_biomass_bins to account for adjusting biomass parameters. Ignored if biomass params not being optimized.')
parser.add_argument('-m', '--minVal_forTranstrum', dest='minVal_forTranstrum', type=float)
parser.add_argument('--obsBA_filename',type=str,help='Specify file for obsBA.')
parser.add_argument('-S','--SSEpctile',type=restrict_SSEpctile,help='Ignore cells whose SSE (summed across all months) is above 95th percentile.')
parser.add_argument('-F','--use_BF', dest='use_BF', action='store_true', help='If present, judge against BF instead of BA.')
parser.add_argument('--newest_unpacked', action='store_true', help='If present, use latest unpacking results.')
parser.add_argument('--flip_RH1', action='store_true', help='If present, reverse sign of RH param1.')
parser.add_argument('--flip_RH2', action='store_true', help='If present, reverse sign of RH param2.')
parser.add_argument('--flip_thetaE', action='store_true', help='If present, reverse sign of thetaE.')


parser.set_defaults(n_years=1,
                    verbose=0,
                    constrain_Ia=False,
                    constrain_AGBparams=False,
                    constrain_popDsupp_NF=False,
                    constrain_RHparams=False,
                    constrain_duration=False,
                    constrain_ROSmax=False,
                    constrain_thetaEsq=False,
                    constrain_magicScalar=False,
                    constrain_strict=False,
                    rescale_params=False,
                    damping_style='additive',
                    f_rh_style='li2012',
                    f_agb_style='li2012',
                    f_theta_style='li2012',
                    param_work_dir='/path/to/param_work_dir/',   # SUBSTITUTE
                    output_dir = '/path/to/output_dir/',   # SUBSTITUTE
                    opts = [1e-11, 1e-15, 1e-15, 100],
                    doBins=False,
                    minVal_forTranstrum = 5.0e3,
                    obsBA_filename = 'estCnstBA_GFED3s_global_2001-2009_LM3res_OTHR.nc',
                    SSEpctile = 100.0,
                    use_BF=False,
                    newest_unpacked=False,
                    flip_RH1=False,
                    flip_RH2=False,
                    flip_thetaE=False)
                    
args = parser.parse_args()
opts = args.opts
    
if args.newest_unpacked:
    obsBA_filename = 'estCnstBA_GFED3s_global_2001-2009_LM3res_OTHR.latest.nc'
else:
    obsBA_filename = 'estCnstBA_GFED3s_global_2001-2009_LM3res_OTHR.nc'
    
args.obsBA_filename = '../' + obsBA_filename
if args.verbose>0:
    print 'Using unpacked fire from ' + obsBA_filename
    
if args.verbose>0:
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('%%% Begin verbose output %%%')
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
if args.verbose>0:
    args_dict = vars(args)
    print('Options:')
    for i in args_dict:
        print('   ' + i + ': ' + str(args_dict[i]))
    print(' ')
    
param_work_dir = args.param_work_dir
if location=='gaea':
    # Change to specified param_work_dir    
    os.chdir(param_work_dir)
    param_work_dir = os.getcwd()
    if param_work_dir[len(param_work_dir)-1] != '/':
        param_work_dir = param_work_dir + '/'

txtFiles_active = param_work_dir + 'txtFiles_active/'
    
# Find working_dirs
working_dirs = args.output_dir
if args.verbose>0:
    print('Looking for output directory in ' + working_dirs)
    
# Add trailing slashes, if missing
txtFiles_orig = param_work_dir + 'txtFiles_orig/'
if working_dirs[len(working_dirs)-1] != '/':
    working_dirs = working_dirs + '/'

print(os.getcwd())


## Set up runs

# Find directory containing output of run
ssr.write_status('Finding runs...',txtFiles_active)

if thisRun == 'latest':
    runDirs_all = glob.glob(working_dirs + '*')
    runDirs_all.sort(key=lambda x: os.path.getmtime(x))
    runDir = runDirs_all[len(runDirs_all)-1]
    n = 1
    while runDir[-4:]=='.zip' or runDir[-4:]=='.tar' or 'paramwork' in runDir:
        n += 1
        runDir = runDirs_all[len(runDirs_all)-n]
    if runDir[-1:] != '/':   # Add trailing /, if not present
        runDir = runDir + '/'
    if args.verbose>0:
        print('Using output in ' + runDir)
    del n
else:
    runDir = working_dirs + thisRun + '/'
    
# What variables are being optimized?
ssr.write_status('Getting varList...',txtFiles_active)
varList = ssr.get_varList(txtFiles_orig,args)
Nvars = len(varList)
Nmonths = 12*args.n_years
if args.verbose>0:
    print('Using last ' + str(Nmonths) + ' months of output.')
    
if args.doBins:
    if 'AGBparam1' in varList or 'AGBparam2' in varList:
        bins_nonAGB = np.loadtxt(txtFiles_orig + 'bins_orig_notAGBparams.txt')
    else:
        bins = np.loadtxt(txtFiles_orig + 'bins_orig.txt')

    
## Import data

# Import run output
ssr.write_status('Importing run data...',txtFiles_active)
estBA_mc, estBA, derivs_mcv, derivs_xv, Ncells, cellDirs = \
   ssr.import_cellRuns_ntrlscndOnly(runDir,varList,args)

# Parse lat-lons
ssr.write_status('Parsing lat-lons...',txtFiles_active)
maps = ssr.import_latLons(runDir, cellDirs, Nmonths)

# Import observed BA
ssr.write_status('Importing observed BA...',txtFiles_active)
obsBA_mc, obsBA = ssr.import_obsBA(maps,Nmonths,args.obsBA_filename,args.use_BF)

if np.any(np.isnan(estBA)):
    raise NameError('At least one element of estBA is NaN.')
elif np.any(np.isnan(obsBA)):
    raise NameError('At least one element of obsBA is NaN.') 
    
# Restrict, if doing so
if args.SSEpctile < 100.0:
    ssr.write_status('Doing SSE restriction...',txtFiles_active)
    if os.path.isfile(txtFiles_active + 'okcells'):
        okcells = np.load(txtFiles_active + 'okcells.npy')
    else:   # I.e., this is the first run
        okcells = ssr.SSErestrict_getOK(obsBA_mc, estBA_mc, args, txtFiles_active)
    Ncells = len(okcells)
    obsBA_mc, obsBA, estBA_mc, estBA, derivs_mcv, derivs_xv = \
        ssr.SSErestrict_doIt(obsBA_mc, obsBA, estBA_mc, estBA, derivs_mcv, \
                             derivs_xv, okcells, Ncells, args)
