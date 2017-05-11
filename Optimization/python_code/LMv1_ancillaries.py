from netCDF4 import Dataset
import numpy as np
import os
import glob
import re
import csv
import time
import shutil
from datetime import datetime as dt
import matplotlib.pyplot as plt
import timeit


def import_cellRuns(runDir,varList,args):
    
    raise NameError('This is obsolete. Do not use!')
    
    cellDirs = glob.glob(runDir + '*[EW]*[NS]*')
    Ncells = len(cellDirs)
    Nmonths = 12*args.n_years
    
    if Ncells == 0:
        raise NameError('No cells found!')
        
    # Which variables are present?
    Nvars = len(varList)
    
    # Initialize
    BA_rate_mc = np.empty((Nmonths,Ncells,))
    BA_rate_mc[:] = np.NAN
    derivs_mcv = np.empty((Nmonths,Ncells,Nvars,))
    derivs_mcv[:] = np.NAN
    
    for c in range(Ncells):
        inFile = cellDirs[c]+'/land_month.nc'
        if not os.path.isfile(inFile):
           raise NameError('File not found: ' + inFile)
        modelOut_ds = Dataset(cellDirs[c]+'/land_month.nc',mode='r')
        if c == 1:
            # Get days_in_year
            if modelOut_ds.variables['average_T1'].units == 'days since 1501-01-01 00:00:00':
                thisYear = 1501 ;
                days_since = max(np.squeeze(modelOut_ds.variables['average_T1']))
                while days_since >= 0:
                    if thisYear%4 == 0:
                        days_since = days_since - 366
                    else:
                        days_since = days_since - 365
                    thisYear = thisYear + 1
                lastYear = thisYear - 1
                year_indices = range(args.n_years)
                days_in_month_list = list()   # Probably not proper since it gets turned into a numpy array.
                for y in year_indices:
                    thisYear_2 = lastYear - (1+year_indices[-(y+1)])
                    if thisYear_2%4==0:
                        days_in_month_list = np.append(days_in_month_list,[31,29,31,30,31,30,31,31,30,31,30,31])
                    else:
                        days_in_month_list = np.append(days_in_month_list,[31,28,31,30,31,30,31,31,30,31,30,31])
            else:
                raise NameError('Figure out how to determine year when average_T1.units = ' + modelOut_ds.variables['average_T1'].units + '.')
            days_in_month = np.array(days_in_month_list)
        BA_rate_mc[:,c] = np.squeeze(modelOut_ds.variables['BA_rate'])[-Nmonths:]
        for v in range(len(varList)):
            thisVar = 'BA_DERIVwrt_' + varList[v]
            derivs_mcv[:,c,v] = np.squeeze(modelOut_ds.variables[thisVar])[-Nmonths:]
        modelOut_ds.close
    
    if np.any(np.isnan(derivs_mcv)):
        raise NameError('At least one element of derivs_mcv is NaN before converting to per-month.')

    # Convert from per-day to per-month
    BA_rate_mc = BA_rate_mc * np.transpose(np.tile(days_in_month,(Ncells,1)))
    derivs_mcv = derivs_mcv * np.transpose(np.tile(days_in_month,(Nvars,Ncells,1)))
    if np.any(np.isnan(derivs_mcv)):
        raise NameError('At least one element of derivs_mcv is NaN after converting to per-month.')
    
    # Vectorize
    BA_rate = np.reshape(BA_rate_mc,Ncells*Nmonths) ; del BA_rate_mc
    derivs_xv = np.reshape(derivs_mcv,(Ncells*Nmonths,Nvars)) ; del derivs_mcv
    
    # Rescale parameters, if doing so
    if args.rescale_params:
        derivs_xv = rescale_params('in',varList,derivs_xv,args)
    
    return BA_rate, derivs_xv, Ncells, cellDirs


def import_cellRuns_ntrlscndOnly(runDir,varList,args):
    
    cellDirs = [f for f in glob.glob(runDir + '*[EW]*[NS]*') \
                if 'zzz' not in f]   # To exclude a cell, add 'zzz' to its dirName.
    Ncells = len(cellDirs)
    Nmonths = 12*args.n_years
    
    if Ncells == 0:
        raise NameError('No cells found!')
        
    # Which variables are present?
    Nvars = len(varList)
    
    # Initialize
    BA_rate_mc_ntrl = np.empty((Nmonths,Ncells,))
    BA_rate_mc_scnd = np.empty((Nmonths,Ncells,))
    derivs_mcv_ntrl = np.empty((Nmonths,Ncells,Nvars,))
    derivs_mcv_scnd = np.empty((Nmonths,Ncells,Nvars,))
    
    for c in range(Ncells):
        tic = timeit.default_timer()
        inFile = cellDirs[c]+'/land_month_ntrl.nc'
        if not os.path.isfile(inFile):
           raise NameError('File not found: ' + inFile)
        modelOut_ds = Dataset(cellDirs[c]+'/land_month_ntrl.nc',mode='r')
        if c == 1:
            # Get days_in_year
            if modelOut_ds.variables['average_T1'].units == 'days since 1501-01-01 00:00:00':
                thisYear = 1501 ;
                days_since = max(np.squeeze(modelOut_ds.variables['average_T1']))
                while days_since >= 0:
                    if thisYear%4 == 0:
                        days_since = days_since - 366
                    else:
                        days_since = days_since - 365
                    thisYear = thisYear + 1
                lastYear = thisYear - 1
                year_indices = range(args.n_years)
                days_in_month_list = list()   # Probably not proper since it gets turned into a numpy array.
                for y in year_indices:
                    thisYear_2 = lastYear - (1+year_indices[-(y+1)])
                    if thisYear_2%4==0:
                        days_in_month_list = np.append(days_in_month_list,[31,29,31,30,31,30,31,31,30,31,30,31])
                    else:
                        days_in_month_list = np.append(days_in_month_list,[31,28,31,30,31,30,31,31,30,31,30,31])
            else:
                raise NameError('Figure out how to determine year when average_T1.units = ' + modelOut_ds.variables['average_T1'].units + '.')
            days_in_month = np.array(days_in_month_list)
        BA_rate_Variable = modelOut_ds.variables['BA_rate_ntrl']
        missing_value = BA_rate_Variable.missing_value
        BA_rate_data = np.squeeze(BA_rate_Variable)[-Nmonths:]
        BA_rate_data[BA_rate_data==missing_value] = 0.0
        if args.use_BF:
            try:
                land_area = modelOut_ds.variables['land_area'][-Nmonths:]
            except KeyError:
                land_static = Dataset(cellDirs[c]+'/land_static.nc',mode='r')
                land_area = land_static.variables['land_area'][-Nmonths:]
            BA_rate_data = BA_rate_data / land_area
            BA_rate_data[land_area==0.0] = 0.0
        
        BA_rate_mc_ntrl[:,c] = BA_rate_data
        del BA_rate_Variable ; del missing_value ; del BA_rate_data
        for v in range(len(varList)):
            thisVar = 'BA_DERIVwrt_' + varList[v] + '_ntrl'
            thisVar_Variable = modelOut_ds.variables[thisVar]
            missing_value = thisVar_Variable.missing_value
            thisVar_data = np.squeeze(thisVar_Variable)[-Nmonths:]
            thisVar_data[thisVar_data==missing_value] = 0.0
            if args.use_BF:
                thisVar_data = thisVar_data / land_area
            derivs_mcv_ntrl[:,c,v] = thisVar_data
            del thisVar ; del thisVar_Variable ; del missing_value ; del thisVar_data
        modelOut_ds.close
        toc = timeit.default_timer()
        if args.verbose>2:
            print('NTRL: Cell ' + str(c) + ' took ' + str(round(toc-tic,1)) + ' secs')
        
        
    for c in range(Ncells):
        tic = timeit.default_timer()
        inFile = cellDirs[c]+'/land_month_scnd.nc'
        if not os.path.isfile(inFile):
           raise NameError('File not found: ' + inFile)
        modelOut_ds = Dataset(cellDirs[c]+'/land_month_scnd.nc',mode='r')
        if c == 1:
            # Get days_in_year
            if modelOut_ds.variables['average_T1'].units == 'days since 1501-01-01 00:00:00':
                thisYear = 1501 ;
                days_since = max(np.squeeze(modelOut_ds.variables['average_T1']))
                while days_since >= 0:
                    if thisYear%4 == 0:
                        days_since = days_since - 366
                    else:
                        days_since = days_since - 365
                    thisYear = thisYear + 1
                lastYear = thisYear - 1
                year_indices = range(args.n_years)
                days_in_month_list = list()   # Probably not proper since it gets turned into a numpy array.
                for y in year_indices:
                    thisYear_2 = lastYear - (1+year_indices[-(y+1)])
                    if thisYear_2%4==0:
                        days_in_month_list = np.append(days_in_month_list,[31,29,31,30,31,30,31,31,30,31,30,31])
                    else:
                        days_in_month_list = np.append(days_in_month_list,[31,28,31,30,31,30,31,31,30,31,30,31])
            else:
                raise NameError('Figure out how to determine year when average_T1.units = ' + modelOut_ds.variables['average_T1'].units + '.')
            days_in_month = np.array(days_in_month_list)
        BA_rate_Variable = modelOut_ds.variables['BA_rate_scnd']
        missing_value = BA_rate_Variable.missing_value
        BA_rate_data = np.squeeze(BA_rate_Variable)[-Nmonths:]
        BA_rate_data[BA_rate_data==missing_value] = 0.0
        if args.use_BF:
            BA_rate_data = BA_rate_data / land_area
            BA_rate_data[land_area==0.0] = 0.0
            
        BA_rate_mc_scnd[:,c] = BA_rate_data
        del BA_rate_Variable ; del missing_value ; del BA_rate_data
        for v in range(len(varList)):
            thisVar = 'BA_DERIVwrt_' + varList[v] + '_scnd'
            thisVar_Variable = modelOut_ds.variables[thisVar]
            missing_value = thisVar_Variable.missing_value
            thisVar_data = np.squeeze(thisVar_Variable)[-Nmonths:]
            thisVar_data[thisVar_data==missing_value] = 0.0
            if args.use_BF:
                thisVar_data = thisVar_data / land_area
            derivs_mcv_scnd[:,c,v] = thisVar_data
            del thisVar ; del thisVar_Variable ; del missing_value ; del thisVar_data
        modelOut_ds.close
        toc = timeit.default_timer()
        if args.verbose>2:
            print('SCND: Cell ' + str(c) + ' took ' + str(round(toc-tic,1)) + ' secs')
        
    # Account for missing values
    if np.any(np.isnan(derivs_mcv_ntrl)):
        raise NameError('At least one element of derivs_mcv_ntrl is NaN before converting to per-month.')
    elif np.any(np.isnan(derivs_mcv_scnd)):
        raise NameError('At least one element of derivs_mcv_scnd is NaN before converting to per-month.')

    # Combine NTRL + SCND
    BA_rate_mc = BA_rate_mc_ntrl + BA_rate_mc_scnd
    derivs_mcv = derivs_mcv_ntrl + derivs_mcv_scnd
    
    # Flip sign, if doing so
    if args.flip_RH1:
        if 'RHparam1' in varList:
            print('Flipping sign of derivatives w/r/t RHparam1.')
            derivs_mcv[:,:,varList.index('RHparam1')] = -1.0 * derivs_mcv[:,:,varList.index('RHparam1')]
        else:
            print('flip_RH1==TRUE but RHparam1 not in varList. Ignoring.')
    if args.flip_RH2:
        if 'RHparam2' in varList:
            print('Flipping sign of derivatives w/r/t RHparam2.')
            derivs_mcv[:,:,varList.index('RHparam2')] = -1.0 * derivs_mcv[:,:,varList.index('RHparam2')]
        else:
            print('flip_RH2==TRUE but RHparam2 not in varList. Ignoring.')
    if args.flip_thetaE:
        if 'thetaE' in varList:
            print('Flipping sign of derivatives w/r/t thetaE.')
            derivs_mcv[:,:,varList.index('thetaE')] = -1.0 * derivs_mcv[:,:,varList.index('thetaE')]
        else:
            print('flip_thetaE==TRUE but thetaE not in varList. Ignoring.')
    
    # Convert from per-day to per-month
    BA_rate_mc = BA_rate_mc * np.transpose(np.tile(days_in_month,(Ncells,1)))
    derivs_mcv = derivs_mcv * np.transpose(np.tile(days_in_month,(Nvars,Ncells,1)))
    
    if np.any(np.isnan(derivs_mcv)):
        raise NameError('At least one element of derivs_mcv is NaN after converting to per-month.')
    
    # Vectorize
    BA_rate = np.reshape(BA_rate_mc,Ncells*Nmonths)
    derivs_xv = np.reshape(derivs_mcv,(Ncells*Nmonths,Nvars))
    
    # Rescale parameters, if doing so
    if args.rescale_params:
        derivs_xv = rescale_params('in',varList,derivs_xv,args)
            
    return BA_rate_mc, BA_rate, \
           derivs_mcv, derivs_xv, \
           Ncells, cellDirs
                    
        
def import_latLons(runDir, cellDirs, Nmonths):
    
    Ncells = len(cellDirs)
    
    lons = np.empty((Ncells,1,)) ; lons[:] = np.NAN
    lats = np.empty((Ncells,1,)) ; lats[:] = np.NAN
    for c in range(Ncells):
        LatLon_str = cellDirs[c].replace(runDir,'')
        lon_str = re.search('(.+?)x',LatLon_str).group(1)   # Extract everything before 'x'
        lat_str = re.search('x(.+)',LatLon_str).group(1)    # Extract everything after 'x'
        if lon_str[-1:] == '/':      # Remove trailing /, if necessary
            lon_str = lon_str[:-1]
        if lat_str[-1:] == '/':      # Remove trailing /, if necessary
            lat_str = lat_str[:-1]
        lon_str = lon_str.replace('z','',2)
        lat_str = lat_str.replace('z','',2)
        lon_tmp = float(lon_str[:len(lon_str)-1])
        lat_tmp = float(lat_str[:len(lat_str)-1])
        if lon_str[len(lon_str)-1] == 'W':
            lon_tmp = lon_tmp * -1.
        if lat_str[len(lat_str)-1] == 'S':
            lat_tmp = lat_tmp * -1.
        lons[c] = lon_tmp
        lats[c] = lat_tmp
        
    vec_lon = np.arange(-178.75,180.,2.5)
    vec_lat = np.arange(-89.,90.,2.)
    map_lat = np.flipud(np.transpose(np.tile(vec_lat,(144,1))))
    map_lon = np.tile(vec_lon,(90,1))
    map_lon_Marray = np.transpose(np.tile(map_lon,(Nmonths,1,1)),(1,2,0))
    map_lat_Marray = np.transpose(np.tile(map_lat,(Nmonths,1,1)),(1,2,0))
    
    return dict(lons=lons, lats=lats, 
                lon_Marray=map_lon_Marray,
                lat_Marray=map_lat_Marray)
    

def import_obsBA(maps,Nmonths,obsBA_filename,use_BF):
    
    Ncells = len(maps['lats'])
    
    if not os.path.isfile(obsBA_filename):
        raise NameError(obsBA_filename + 'does not exist!')
        
    if obsBA_filename[-4:]=='.mat':
        mat = scipy.io.loadmat(obsBA_filename)
        if use_BF:
            tmp = mat['estBF_othr_LM3_YXm']
        else:
            tmp = mat['estBA_othr_LM3_YXm']
    elif obsBA_filename[-3:]=='.nc':
        obsBA_file = Dataset(obsBA_filename,mode='r')
        tmp  = obsBA_file.variables['BA_othr'][:]
        obsBA_file.close()
    else:
        raise NameError('What kind of file is ' +  obsBA_filename + '?')
    
    
    if tmp.shape==(144,90,108):
        obsBA_tmp = np.flipud(np.transpose(tmp,(1, 0, 2)))
    elif tmp.shape==(108,90,144):
        obsBA_tmp = np.flipud(np.transpose(tmp,(1, 2, 0)))
    elif tmp.shape==(90,144,108):
        obsBA_tmp = tmp
    else:
        raise NameError('What shape is this?')
    
    
    obsBA_mc = np.empty((Nmonths,Ncells,)) ; obsBA_mc[:] = np.NAN
    for c in range(Ncells):
        obsBA_mc[:,c] = np.squeeze(obsBA_tmp[np.where((maps['lon_Marray']==maps['lons'][c]) & (maps['lat_Marray']==maps['lats'][c]))])
    
    obsBA = np.reshape(obsBA_mc,Ncells*Nmonths)
    return obsBA_mc, obsBA
    

def checkStop_ng(g,opts,Niters,txtFiles_active):
    ng = np.linalg.norm(g , ord=np.inf)
    stop = False
    if ng <= opts[1]:
        stop = True
        save_bashDo([str(-1), str(ng), str(opts[1]), dt.now().strftime("%Y%m%d_%H%M%S_%f")],Niters,txtFiles_active)
    return stop
    
    
def calcA(J):
    A = np.dot( np.transpose(J) , J )

    if np.any(np.isnan(A)):
        raise NameError('At least one element of A is NaN.')

    return A
    
    
def calcAdamped(A,mu,Nvars,varList,args):
    if args.damping_style=='additive':
        # # Additive damping can have issues with badly-scaled problems
        # # but protects against a rank-deficient Jacobian (Lampton, 1997).
        A_damped = A + mu*np.eye(Nvars)
    elif args.damping_style=='multiplicative':
        A_damped = np.copy(A)
        np.fill_diagonal(A_damped,np.diag(A)*(1+mu))
    elif args.damping_style=='transtrum_minVal':
        # Just like "multiplicative," but with minimum values specified for
        # diagonal. Using 1 as minimum value for now, but Transtrum (email)
        # notes that some scaling can be done based on what I know about scaling
        # of parameters.
        A_damped_tmp = np.copy(A)
        np.fill_diagonal(A_damped_tmp,np.diag(A)*(1+mu))
        A_damped = np.copy(A_damped_tmp)
        if any(np.diag(A_damped)<0):
            raise NameError('At least one element of diag(A_damped) is negative. HOW?? Fix minVal to use absolute value.')
        for v in range(Nvars):
            A_damped[v,v] = max(A_damped[v,v],args.minVal_forTranstrum)
    elif args.damping_style=='transtrum_minVal_v2':
        # Just like "multiplicative," but with minimum values specified for
        # diagonal. With this second version, the minimum value varies by
        # parameter, after Transtrum (email).
        A_damped_tmp = np.copy(A)
        np.fill_diagonal(A_damped_tmp,np.diag(A)*(1+mu))
        A_damped = np.copy(A_damped_tmp)
        minVal_list = args.minVal_forTranstrum * np.power(rescale_params('out',varList,np.ones_like(np.diag(A_damped)),args),2)
        if any(np.diag(A_damped)<0):
            raise NameError('At least one element of diag(A_damped) is negative. HOW?? Fix minVal to use absolute value.')
        for v in range(Nvars):
            thisMinVal = minVal_list[v]
            A_damped[v,v] = max(A_damped[v,v],thisMinVal)
    else:
        raise NameError('Invalid value for args.damping_style.')
        
    if args.verbose==2:
        print('A = \n' + str(A))
        if args.damping_style=='transtrum_minVal' or args.damping_style=='transtrum_minVal_v2':
            print('A_damped_tmp = \n' + str(A_damped_tmp))
        if args.damping_style=='transtrum_minVal_v2':
            print('minVal_list = \n' + str(minVal_list))
        print('diag(A) = ')
        print(np.diag(A))
        print('diag(A_damped, using ' + args.damping_style + ') = ')
        print(np.diag(A_damped))
        
    if args.verbose>=1:
        if np.array_equal(A,A_damped):
            print('No damping actually occurred.')
        else:
            print('A did get damped.')
    
    return A_damped
    
    
def calcBins(bins_nonAGB,AGBparam1,AGBparam2):
    
    if AGBparam1 == AGBparam2:
        if AGBparam1 == 0.0:
            newBins = np.append(np.append(bins_nonAGB,9999.2),9999.3)
        else:
            newBins = np.append(np.append(bins_nonAGB,AGBparam1),9999.2)
    elif AGBparam1 == 0.0:
        newBins = np.append(np.append(bins_nonAGB,AGBparam2),9999.2)
    else:
        newBins = np.append(np.append(bins_nonAGB,AGBparam1),AGBparam2)
        
    newBins.sort()
        
    return newBins
    
        
def calcH(J,g,mu,nu,x,opts,Niters,txtFiles_active,args,varList):
    Nvars = len(g)
    
    # Calculate A(J)
    write_status('Calculating A...',txtFiles_active)
    A = calcA(J)
        
    # Calculate delta to next guess
    write_status('Calculating A_damped...',txtFiles_active)
    A_damped = calcAdamped(A,mu,Nvars,varList,args)
    write_status('Calculating h...',txtFiles_active)
    h = np.dot( np.linalg.inv(A_damped) , -1*g)
        
    if args.rescale_params:
        h = rescale_params('out',varList,h,args)
    if np.any(np.isnan(h)):
        raise NameError('At least one element of h is NaN.')
        
    if args.verbose>0:
        if args.verbose==2:
            print('inv(A_damped):' + str(np.linalg.inv(A_damped)))
        print('h = '+ str(h))
        print('Frac. change = '+ str(h/x))
        
    # Calculate xnew
    write_status('Calculating xnew...',txtFiles_active)
    xnew = x + h
    
    # Constrain! (If doing so.)
    write_status('Constraining xnew...',txtFiles_active)
    xnew = constrain_params(xnew,args,varList)
    h = xnew - x
    
    # If necessary, update h (and xnew) until dL>0
    dL = np.dot( np.transpose(h) , mu*h-g ) / 2
    if dL <= 0:
        write_status('Updating h, xnew until dL>0...',txtFiles_active)
    if args.verbose>0:
        print('dL = ' + str(dL))
        if dL <= 0:
            print('% % % % % % % % % % % % %')
            print('% dL <=0, so repeating % ')
            print('% % % % % % % % % % % % %')
    while dL <= 0:
        mu = mu * nu
        nu = nu * 2
        A_damped = calcAdamped(A,mu,Nvars,varList,args)
        h = np.dot( np.linalg.inv(A_damped) , -1*g)
                
        if args.rescale_params:
            h = rescale_params('out',varList,h)
        xnew = x + h
        xnew_preconstrained = np.copy(xnew)
        
        # Constrain! (If doing so.)
        xnew = constrain_params(xnew,args,varList)
        h = xnew - x
        
        dL = np.dot( np.transpose(h) , mu*h-g ) / 2
        
        if args.verbose>0:
            print('% mu = ' + str(mu))
            print('% nu = ' + str(nu))
            if args.verbose==2:
                print('% A_damped = \n %' + str(A_damped))
            print('% h = ' + str(h))
            print('% xnew = ' + str(xnew_preconstrained))
            print('% xnew_constrained = ' + str(xnew))
            print('% dL = ' + str(dL))
            
        if np.any(np.isnan(h)):
            raise NameError('At least one element of h is NaN.')
        elif np.isnan(dL):
            raise NameError('dL is NaN.')
        elif np.isinf(dL):
            raise NameError('dL is +/-Inf.')
            
        if dL>0 and args.verbose>0:
            print('% % % % % % % % % % %')
            print('% Ending this loop % ')
            print('% % % % % % % % % % %')
    
    # Check stopping criteria
    write_status('Checking stop criteria in calcH...',txtFiles_active)
    nh = np.linalg.norm(h,2)
    nx = np.linalg.norm(x,2) + opts[2]
    stop = False
    if nh <= opts[2] * nx:
        stop = True
        save_bashDo([str(-2), str(nh), str(opts[2]*nx), dt.now().strftime("%Y%m%d_%H%M%S_%f")],Niters,txtFiles_active)
        print('Stopping: nh <= opts[2] * nx')
    elif nh >= nx / np.spacing(1):
        stop = True
        save_bashDo([str(-4), str(nh), str(nx/np.spacing(1)), dt.now().strftime("%Y%m%d_%H%M%S_%f")],Niters,txtFiles_active)
        print('Stopping: nh >= nx / ' + str(np.spacing(1)))
    if args.verbose>0:
        print('nh = '+ str(nh))
        print('nx = '+ str(nx))
        
        
    if not stop:        
        # Save dL
        write_status('Saving dL...',txtFiles_active)
        if Niters==1: # If first iteration, erase existing file.
            dL_record = open(txtFiles_active + '/dL_record.txt','wb')
        else:
            dL_record = open(txtFiles_active + '/dL_record.txt','ab')
        writer = csv.writer(dL_record, delimiter=' ', quoting=csv.QUOTE_NONE)
        writer.writerow([dL, dt.now().strftime("%Y%m%d_%H%M%S_%f")])
        dL_record.close() ; del writer
        np.save(txtFiles_active + 'lastRun_dL',dL)
        
    
    # Return
    return A, A_damped, h, xnew, dL, mu, nu, stop
    
    
def calcSSE(estBA,obsBA):
    
    # Calculate errors. Needs to be estBA-obsBA... Here's why:
    # ### Let's say BA increases with Param, so that dBA/dParam > 0. For a cell
    # ### with too much BA, SSE should increase with Param. So we want
    # ### dSSE/dParam = 2*f*dBA/dParam to be > 0 as well. Since we said that
    # ### dBA/dParam > 0, then f must also be > 0.
    f = estBA - obsBA
    
    # Calculate sum of squared errors. Divide by 2 for some reason?
    SSEhalf = np.dot(np.transpose(f),f) / 2
    
    return f, SSEhalf
    
    
def readLatest_bashDo(txtFiles_active):
    try:
        bash_do = open(txtFiles_active + 'bash_do.txt','r')
        bash_do_lineList = bash_do.readlines()
        bash_do.close() ; del bash_do
        Niters = len(bash_do_lineList)
        last_bashDo = int(bash_do_lineList[Niters-1].split()[0])
    except IOError:
        # bash_do.txt does not yet exist, which means that this is happening
        # after the very first model run
        last_bashDo = 0
    return last_bashDo
    
    
def readLatest_fJ_acc(txtFiles_active,args):
    
    f = np.load(txtFiles_active + 'lastAcc_f.npy')
    J = np.load(txtFiles_active + 'lastAcc_J.npy')
    
    return f, J


def readLatest_muNu(txtFiles_active):
    mu = np.load(txtFiles_active + 'last_mu.npy')
    nu = np.load(txtFiles_active + 'last_nu.npy')
    
    return mu, nu

    
def readLatest_params_run(txtFiles_active,Nvars):
    
    x = np.load(txtFiles_active + 'lastRun_params.npy')
    
    return x
    
    
def readLatest_params_acc(txtFiles_active,Nvars,Niters):
    
    if Niters > 1:
        x = np.load(txtFiles_active + 'lastAcc_params.npy')
    else:
        params_run = open(txtFiles_active + 'params_acc.txt','r')
        params_run_lineList = params_run.readlines()
        params_run.close() ; del params_run
        Niters = len(params_run_lineList)
        x = np.empty(Nvars)
        for v in range(Nvars):
            x[v] = float(params_run_lineList[Niters-1].split()[v])
        del params_run_lineList
    
    return x
    
    
def get_Niters(txtFiles_active):
    params_run = open(txtFiles_active + 'params_run.txt','r')
    params_run_lineList = params_run.readlines()
    params_run.close() ; del params_run
    Niters = len(params_run_lineList)
    return Niters
    

def save_bashDo(thisRow, Niters, txtFiles_active):
    if Niters==1: # If first iteration, erase existing file. 
        bash_do = open(txtFiles_active + 'bash_do.txt','wb')
    else:
        bash_do = open(txtFiles_active + 'bash_do.txt','ab')
    writer = csv.writer(bash_do, delimiter=' ', quoting=csv.QUOTE_NONE)
    writer.writerow(thisRow)
    bash_do.close() ; del writer
    
    
def save_fJ(f,J,txtFiles_active,args):
    
    np.save(txtFiles_active + 'lastAcc_f',f)
    np.save(txtFiles_active + 'lastAcc_J',J)
    
    
def save_muNu(mu, nu, Niters,txtFiles_active):
    if Niters==1: # If first iteration, erase existing file.
        mu_nu = open(txtFiles_active + 'mu_nu.txt','wb')
    else:
        mu_nu = open(txtFiles_active + 'mu_nu.txt','ab')
    writer = csv.writer(mu_nu, delimiter=' ', quoting=csv.QUOTE_NONE)
    writer.writerow([mu, nu, dt.now().strftime("%Y%m%d_%H%M%S_%f")])
    mu_nu.close() ; del writer
    
    np.save(txtFiles_active + 'last_mu',mu)
    np.save(txtFiles_active + 'last_nu',nu)


def save_nextGuess(xnew,txtFiles_active,varList,newBins,args):
    
    # Save next-guess parameter set to lastRun_params.npy
    np.save(txtFiles_active + 'lastRun_params',xnew)
    
    # Save next-guess parameter set to params_run.txt
    params_run = open(txtFiles_active + 'params_run.txt','ab')
    writer = csv.writer(params_run, delimiter=' ', quoting=csv.QUOTE_NONE)
    row2write = map(str,xnew)
    row2write.append(dt.now().strftime("%Y%m%d_%H%M%S_%f"))
    writer.writerow(row2write)
    params_run.close() ; del writer
    
    # Save next-guess bins to vegn_data.nml
    if args.doBins:
        shutil.copyfile('txtFiles_orig/vegn_data_incomplete.nml',txtFiles_active + '/vegn_data.nml')
        vegn_data_nml = open(txtFiles_active + 'vegn_data.nml','ab')
        newBins_str = ''
        for b in range(len(newBins)):
            newBins_str = newBins_str + str(newBins[b]) + ' '
        vegn_data_nml.write('         scnd_biomass_bins = ' + newBins_str + '\n')
        vegn_data_nml.write('/')
        vegn_data_nml.close()
    
    # Save next-guess parameter set to fire.nml
    shutil.copyfile('txtFiles_orig/fire_incomplete.nml',txtFiles_active + '/fire.nml')
    fire_nml = open(txtFiles_active + 'fire.nml','ab')
#    fire_nml.write('\n')
    for v in range(len(varList)):
        thisVar = varList[v]
        if thisVar=='AGBparam1':
            if args.f_agb_style=='li2012':
                fire_nml.write('      agb_lo = ' + str(xnew[v]) + '\n')
            elif args.f_agb_style=='logistic':
                fire_nml.write('      agb_psi2 = ' + str(xnew[v]) + '\n')
            elif args.f_agb_style=='gompertz':
                fire_nml.write('      agb_gom2 = ' + str(xnew[v]) + '\n')
        elif thisVar=='AGBparam2':
            if args.f_agb_style=='li2012':
                fire_nml.write('      agb_up = ' + str(xnew[v]) + '\n')
            elif args.f_agb_style=='logistic':
                fire_nml.write('      agb_psi3 = ' + str(xnew[v]) + '\n')
            elif args.f_agb_style=='gompertz':
                fire_nml.write('      agb_gom3 = ' + str(xnew[v]) + '\n')
        elif thisVar=='RHparam1':
            if args.f_rh_style=='li2012':
                fire_nml.write('      rh_lo = ' + str(xnew[v]) + '\n')
            elif args.f_rh_style=='logistic':
                fire_nml.write('      rh_psi2 = ' + str(xnew[v]) + '\n')
            elif args.f_rh_style=='gompertz':
                fire_nml.write('      rh_gom2 = ' + str(xnew[v]) + '\n')
        elif thisVar=='RHparam2':
            if args.f_rh_style=='li2012':
                fire_nml.write('      rh_up = ' + str(xnew[v]) + '\n')
            elif args.f_rh_style=='logistic':
                fire_nml.write('      rh_psi3 = ' + str(xnew[v]) + '\n')
            elif args.f_rh_style=='gompertz':
                fire_nml.write('      rh_gom3 = ' + str(xnew[v]) + '\n')
        elif thisVar=='alphaM':
            fire_nml.write('      Ia_alpha_monthly = ' + str(xnew[v]) + '\n')
        elif thisVar=='IaParam1':
            fire_nml.write('      Ia_param1 = ' + str(xnew[v]) + '\n')
        elif thisVar=='IaParam2':
            fire_nml.write('      Ia_param2 = ' + str(xnew[v]) + '\n')
        elif thisVar=='popDsupp_NF_eps1':
            fire_nml.write('      popD_supp_eps1 = ' + str(xnew[v]) + '\n')
        elif thisVar=='popDsupp_NF_eps2':
            fire_nml.write('      popD_supp_eps2 = ' + str(xnew[v]) + '\n')
        elif thisVar=='popDsupp_NF_eps3':
            fire_nml.write('      popD_supp_eps3 = ' + str(xnew[v]) + '\n')
        elif thisVar=='thetaE':
            if args.f_theta_style=='logistic' or args.f_theta_style=='gompertz':
                raise NameError('thisVar==thetaE but f_theta_style==' + args.f_theta_style)
            fire_nml.write('      theta_extinction = ' + str(xnew[v]) + '\n')
        elif thisVar=='THETAparam1':
            if args.f_theta_style=='li2012':
                raise NameError('thisVar==THETAparam1 but f_theta_style==li2012')
            elif args.f_theta_style=='logistic':
                fire_nml.write('      theta_psi2 = ' + str(xnew[v]) + '\n')
            elif args.f_theta_style=='gompertz':
                fire_nml.write('      theta_gom2 = ' + str(xnew[v]) + '\n')
        elif thisVar=='THETAparam2':
            if args.f_theta_style=='li2012':
                raise NameError('thisVar==THETAparam2 but f_theta_style==li2012')
            elif args.f_theta_style=='logistic':
                fire_nml.write('      theta_psi3 = ' + str(xnew[v]) + '\n')
            elif args.f_theta_style=='gompertz':
                fire_nml.write('      theta_gom3 = ' + str(xnew[v]) + '\n')
        elif thisVar=='fireDur_c4':
            fire_nml.write('      fire_duration_ave_c4 = ' + str(xnew[v]) + '\n')
        elif thisVar=='fireDur_c3':
            fire_nml.write('      fire_duration_ave_c3 = ' + str(xnew[v]) + '\n')
        elif thisVar=='fireDur_dt':
            fire_nml.write('      fire_duration_ave_dt = ' + str(xnew[v]) + '\n')
        elif thisVar=='fireDur_tt':
            fire_nml.write('      fire_duration_ave_tt = ' + str(xnew[v]) + '\n')
        elif thisVar=='fireDur_et':
            fire_nml.write('      fire_duration_ave_et = ' + str(xnew[v]) + '\n')
        elif thisVar=='ROSmax_c4':
            fire_nml.write('      ROS_max_C4GRASS = ' + str(xnew[v]) + '\n')
        elif thisVar=='ROSmax_c3':
            fire_nml.write('      ROS_max_C3GRASS = ' + str(xnew[v]) + '\n')
        elif thisVar=='ROSmax_dt':
            fire_nml.write('      ROS_max_TEMPDEC = ' + str(xnew[v]) + '\n')
        elif thisVar=='ROSmax_tt':
            fire_nml.write('      ROS_max_TROPICAL = ' + str(xnew[v]) + '\n')
        elif thisVar=='ROSmax_ts' or thisVar=='ROSmax_tsav':
            fire_nml.write('      ROS_max_TROPSAV = ' + str(xnew[v]) + '\n')
        elif thisVar=='ROSmax_tshr':
            fire_nml.write('      ROS_max_TROPSHR = ' + str(xnew[v]) + '\n')
        elif thisVar=='ROSmax_et':
            fire_nml.write('      ROS_max_EVERGR = ' + str(xnew[v]) + '\n')
        elif thisVar=='magicScalar':
            fire_nml.write('      magic_scalar = ' + str(xnew[v]) + '\n')
        else:
            fire_nml.write('      ' + varList[v] + ' = ' + str(xnew[v]) + '\n')
    fire_nml.write('/')
    fire_nml.close()
    
    # Concatenate input_incomplete.nml, vegn_data.nml, and fire.nml
    if args.doBins:
        filenames = ['txtFiles_orig/input_incomplete.nml', txtFiles_active + '/vegn_data.nml', txtFiles_active + '/fire.nml']
    else:
        filenames = ['txtFiles_orig/input_incomplete.nml', txtFiles_active + '/fire.nml']
    with open(txtFiles_active + '/input.nml', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                outfile.write(infile.read())
    

def save_txtFiles_active(runDir,working_dirs,txtFiles_active):
    
    histDir = txtFiles_active + '../txtFiles_active_history/'
    histDir_this = histDir + 'after_' + runDir.replace(working_dirs,'')
        
    # Copy files. Based on 
    # http://pythoncentral.io/how-to-recursively-copy-a-directory-folder-in-python/
    try:
        shutil.copytree(txtFiles_active, histDir_this)
    # Directories are the same
    except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
    # Any error saying that the directory doesn't exist
    except OSError as e:
        print('Directory not copied. Error: %s' % e)
        

def write_lineBreak(file):
    file_handle = open(file,'ab')
    file_handle.write('\n')
    file_handle.close()
    
def get_varList(txtFiles_orig,args):
    
    varList = open(txtFiles_orig + 'parameters_present.txt','r').read().splitlines()
    
    # Remove any spaces
    varList = [v.replace(' ','') for v in varList]
    
    return varList
    
def constrain_params(xnew_in,args,varList):
    
    if args.verbose>0:
        print('xnew = ' + str(xnew_in))
    did_constrain = False ;

    if any(np.isnan(xnew_in)):
        raise NameError('At least one xnew value is NaN before constraining.')

    xnew_out = np.copy(xnew_in)
    
    if args.constrain_AGBparams:
        did_constrain = True
        if args.f_agb_style=='li2012':
            # AGB parameter 1 (agb_lo) can't be < 0
            if 'AGBparam1' in varList and xnew_out[varList.index('AGBparam1')] < 0.0:
                xnew_out[varList.index('AGBparam1')] = 0.0
            # AGB parameter 2 (agb_up) can't be < agb_lo
            if 'AGBparam2' in varList and xnew_out[varList.index('AGBparam2')] < xnew_out[varList.index('AGBparam1')]:
                xnew_out[varList.index('AGBparam2')] = xnew_out[varList.index('AGBparam1')]
        elif args.f_agb_style=='logistic':
            raise NameError('How do I constrain when f_agb_style==logistic?')
        elif args.f_agb_style=='gompertz':
            # Don't let f_agb(0)>0.001
            if 'AGBparam1' in varList and np.exp(-xnew_out[varList.index('AGBparam1')]) > 0.001:
                xnew_out[varList.index('AGBparam1')] = -np.log(0.001)
            
    if args.constrain_Ia:
        did_constrain = True
        # alpha_m for Ia can't be < 0
        if 'alphaM' in varList and xnew_out[varList.index('alphaM')] < 0.0:
            xnew_out[varList.index('alphaM')] = 0.0
        # IaParam1 can't be < 0
        if 'IaParam1' in varList and xnew_out[varList.index('IaParam1')] < 0.0:
            xnew_out[varList.index('IaParam1')] = 0.0
        # IaParam2 can't be > 1, or else Ia DECREASES with popD---that's the job
        # of popDsupp_NF! But only if STRICTLY constraining.
        if 'IaParam2' in varList and xnew_out[varList.index('IaParam2')] > 1.0 and args.constrain_strict:
            xnew_out[varList.index('IaParam2')] = 1.0
        # IaParam2 can't be < 0, or else you can get Ia near Infinity for
        # realistic values of popD. ("Runaway Ia.")
        if 'IaParam2' in varList and xnew_out[varList.index('IaParam2')] < 0.0:
            xnew_out[varList.index('IaParam2')] = 0.0
            
    if args.constrain_popDsupp_NF:
        did_constrain = True
        # Fraction of fires suppressed in densely-populated regions... 
        if 'popDsupp_NF_eps1' in varList:
            # ...can't be > 1
            if xnew_out[varList.index('popDsupp_NF_eps1')] > 1.0:
                xnew_out[varList.index('popDsupp_NF_eps1')] = 1.0
            # ...can't be < 0
            elif xnew_out[varList.index('popDsupp_NF_eps1')] < 0.0:
                xnew_out[varList.index('popDsupp_NF_eps1')] = 0.0
        # Fraction of fires suppressed in uninhabited regions...
        if 'popDsupp_NF_eps2' in varList:
            # ...can't be > popDsupp_NF_eps1
            if xnew_out[varList.index('popDsupp_NF_eps2')] > xnew_out[varList.index('popDsupp_NF_eps1')]:
                xnew_out[varList.index('popDsupp_NF_eps2')] = xnew_out[varList.index('popDsupp_NF_eps1')]
            # ...can't be < 0
            elif xnew_out[varList.index('popDsupp_NF_eps2')] < 0.0:
                xnew_out[varList.index('popDsupp_NF_eps2')] = 0.0
        # eps3 can't be negative, because that will result in negative fire 
        if 'popDsupp_NF_eps3' in varList:
            if xnew_out[varList.index('popDsupp_NF_eps3')] < 0.0:
                xnew_out[varList.index('popDsupp_NF_eps3')] = 0.0
                
    if args.constrain_RHparams:
        did_constrain = True
        if args.f_rh_style=='li2012':
            # RH parameter 1 (rh_lo)...
            if 'RHparam1' in varList:
                # ...can't be < 0
                if xnew_out[varList.index('RHparam1')] < 0.0:
                    xnew_out[varList.index('RHparam1')] = 0.0
                # ...can't be > 1
                elif xnew_out[varList.index('RHparam1')] > 1.0:
                    xnew_out[varList.index('RHparam1')] = 1.0
            # RH parameter 2 (rh_up)...
            if 'RHparam2' in varList:
                if 'RHparam1' in varList:
                    # ...can't be <= rh_lo
                    if xnew_out[varList.index('RHparam2')] < xnew_out[varList.index('RHparam1')]:
                        # Bump up to RH_lo + max(0.01,0.01*RH_lo)
                        xnew_out[varList.index('RHparam2')] = xnew_out[varList.index('RHparam1')] + max(0.01,0.01*xnew_out[varList.index('RHparam1')])
                # ...can't be > 1 if STRICTLY constraining
                if xnew_out[varList.index('RHparam2')] > 1.0 and args.constrain_strict:
                    xnew_out[varList.index('RHparam2')] = 1.0
        elif args.f_rh_style=='logistic':
            if 'RHparam1' in varList:
                # ...can't be < 0
                if xnew_out[varList.index('RHparam1')] < 0.0:
                    xnew_out[varList.index('RHparam1')] = 0.0
        elif args.f_rh_style=='gompertz':
            print('f_rh_style==gompertz, so not constraining. Maybe you should?')
                
    if args.constrain_duration:
        did_constrain = True
        if 'fireDur_c4' in varList and xnew_out[varList.index('fireDur_c4')] < 0.0:
            xnew_out[varList.index('fireDur_c4')] = 0.0
        if 'fireDur_c3' in varList and xnew_out[varList.index('fireDur_c3')] < 0.0:
            xnew_out[varList.index('fireDur_c3')] = 0.0
        if 'fireDur_dt' in varList and xnew_out[varList.index('fireDur_dt')] < 0.0:
            xnew_out[varList.index('fireDur_dt')] = 0.0
        if 'fireDur_tt' in varList and xnew_out[varList.index('fireDur_tt')] < 0.0:
            xnew_out[varList.index('fireDur_tt')] = 0.0
        if 'fireDur_et' in varList and xnew_out[varList.index('fireDur_et')] < 0.0:
            xnew_out[varList.index('fireDur_et')] = 0.0
            
    if args.constrain_ROSmax:
        did_constrain = True
        if 'ROSmax_c4' in varList and xnew_out[varList.index('ROSmax_c4')] < 0.0:
            xnew_out[varList.index('ROSmax_c4')] = 0.0
        if 'ROSmax_c3' in varList and xnew_out[varList.index('ROSmax_c3')] < 0.0:
            xnew_out[varList.index('ROSmax_c3')] = 0.0
        if 'ROSmax_dt' in varList and xnew_out[varList.index('ROSmax_dt')] < 0.0:
            xnew_out[varList.index('ROSmax_dt')] = 0.0
        if 'ROSmax_tt' in varList and xnew_out[varList.index('ROSmax_tt')] < 0.0:
            xnew_out[varList.index('ROSmax_tt')] = 0.0
        if 'ROSmax_ts' in varList and xnew_out[varList.index('ROSmax_ts')] < 0.0:
            xnew_out[varList.index('ROSmax_ts')] = 0.0
        if 'ROSmax_tshr' in varList and xnew_out[varList.index('ROSmax_tshr')] < 0.0:
            xnew_out[varList.index('ROSmax_tshr')] = 0.0
        if 'ROSmax_tsav' in varList and xnew_out[varList.index('ROSmax_tsav')] < 0.0:
            xnew_out[varList.index('ROSmax_tsav')] = 0.0
        if 'ROSmax_et' in varList and xnew_out[varList.index('ROSmax_et')] < 0.0:
            xnew_out[varList.index('ROSmax_et')] = 0.0
    
    if args.f_theta_style=='li2012':
        if args.constrain_thetaEsq:
            did_constrain = True
            # thetaE (theta_extinction) must be positive
            if 'thetaE' in varList and xnew_out[varList.index('thetaE')] <= 0.001:
                xnew_out[varList.index('thetaE')] = 0.001
    elif args.f_theta_style=='gompertz':
        if 'THETAparam1' in varList:
            # ...can't be < 0
            if xnew_out[varList.index('THETAparam1')] < 0.0:
                xnew_out[varList.index('THETAparam1')] = 0.0
    else:
        print('f_theta_style==' + args.f_theta_style + ', so skipping constrain.')
                
    if args.constrain_magicScalar:
        if 'magicScalar' in varList and xnew_out[varList.index('magicScalar')] <= 0.001:
            did_constrain = True
            xnew_out[varList.index('magicScalar')] = 0.001
            
    if args.verbose>0 and did_constrain:
        print('xnew_constrained = ' + str(xnew_out))

    if any(np.isnan(xnew_out)):
        raise NameError('At least one xnew value is NaN after constraining.')

    return xnew_out
  

def rescale_params(in_or_out,varList,data,args):
    # in_or_out: "in" if converting derivatives, "out" if converting xnew.
    
    Nvars = len(varList)
    
    # Define dictionary for conversions (multiply xnew by this to get actual new
    # parameter guesses). Brings all (initial) parameter values to between 
    # 10^0 and 10^1.
    convFactors = {'alphaM':           1000.0,
                   'IaParam1':            1.0,
                   'IaParam2':           10.0,
                   'popDsupp_NF_eps1':   10.0,
                   'popDsupp_NF_eps2':   10.0,
                   'popDsupp_NF_eps3':  100.0,
                   'fireDur_c4':      10000.0,
                   'fireDur_c3':      10000.0,
                   'fireDur_dt':      10000.0,
                   'fireDur_tt':      10000.0,
                   'fireDur_et':      10000.0,
                   'ROSmax_c4':          10.0,
                   'ROSmax_c3':          10.0,
                   'ROSmax_dt':          10.0,
                   'ROSmax_tt':          10.0,
                   'ROSmax_ts':          10.0,
                   'ROSmax_tshr':        10.0,
                   'ROSmax_tsav':        10.0,
                   'ROSmax_et':          10.0,}
                   
    if args.f_agb_style == 'li2012':
        convFactors['AGBparam1'] = 10.0
        convFactors['AGBparam2'] = 1.0
    elif args.f_agb_style == 'gompertz':
        convFactors['AGBparam1'] = 1.0
        convFactors['AGBparam2'] = 1.0
    else:
        raise NameError('convFactor not defined for AGBparams when f_agb_style==' + args.f_agb_style + '!')
        
    if args.f_rh_style == 'li2012':
        convFactors['RHparam1'] = 10.0
        convFactors['RHparam2'] = 10.0
    elif args.f_rh_style == 'gompertz':
        convFactors['RHparam1'] = 1000.0
        convFactors['RHparam2'] = 1.0
    else:
        raise NameError('convFactor not defined for RHparams when f_rh_style==' + args.f_rh_style + '!')
    
    if args.f_theta_style == 'li2012':
        convFactors['thetaE'] = 10.0
    else:
        raise NameError('convFactor not defined for theta params when f_theta_style==' + args.f_theta_style + '!')
        
    if in_or_out=='in':
        for var in varList:
            try:
                data[:,varList.index(var)] = data[:,varList.index(var)] / convFactors[var]
            except KeyError:
                raise KeyError('convFactor not defined for ' + var + '!')
    elif in_or_out=='out':
        for var in varList:
            try:
                data[varList.index(var)] = data[varList.index(var)] * convFactors[var]
            except KeyError:
                raise KeyError('convFactor not defined for ' + var + '!')
    else:
        raise NameError('Invalid value for in_or_out. Use either "in" or "out".')
        
    return data
    
    
def trace_results(txtFiles_active,Nvars,varList):
    record = np.loadtxt(txtFiles_active + '/SSEhalf_record.txt',usecols=range(0,Nvars+1))
    Nruns = record.shape[0]
    
    # Get record of ACCEPTED runs
    SSEhalf = record[:,0]
    record_acc = record[0,:]
    SSEhalf_lastAcc = SSEhalf[0]
    for r in [x+1 for x in range(0,Nruns-1)]:
        if SSEhalf[r] < SSEhalf_lastAcc:
            record_acc = np.vstack([record_acc,record[r,:]])
            SSEhalf_lastAcc = SSEhalf[r]
        else:
            record_acc = np.vstack([record_acc,np.zeros(Nvars+1)+np.nan])
    record = record[:,1:]
    SSEhalf_acc = record_acc[:,0]
    record_acc = record_acc[:,1:]
    accMask = np.isfinite(SSEhalf_acc)
    
    # Set up figure
    thisFig = plt.figure(1)
    Nsubplots_x = np.floor(np.sqrt(Nvars+1))
    Nsubplots_y = np.ceil(np.sqrt(Nvars+1))
    while Nsubplots_x*Nsubplots_y < Nvars+1:
        if Nsubplots_x <= Nsubplots_y:
            Nsubplots_x = Nsubplots_x + 1
        else:
            Nsubplots_y = Nsubplots_y + 1
    
    # Make figure
    ax = thisFig.add_subplot(Nsubplots_y,Nsubplots_x,1,axisbg='yellow')
    ax.margins(0.05)   # Makes it so points on edges of X-axis aren't cut off by plot box
    ax.tick_params(axis='both', which='major', labelsize=10)
    plt.plot(np.arange(1,Nruns+1),SSEhalf,'rs',
            np.arange(1,Nruns+1)[accMask],SSEhalf_acc[accMask],'b',
            np.arange(1,Nruns+1)[accMask],SSEhalf_acc[accMask],'bs')
    plt.ylabel('SSE (halved)')
    for v in range(Nvars):
        ax = thisFig.add_subplot(Nsubplots_y,Nsubplots_x,v+2)
        ax.margins(0.05)   # Makes it so points on edges of X-axis aren't cut off by plot box
        ax.tick_params(axis='both', which='major', labelsize=10)
        plt.plot(np.arange(1,Nruns+1),record[:,v],'rs',
                np.arange(1,Nruns+1)[accMask],record[:,v][accMask],'b',
                np.arange(1,Nruns+1)[accMask],record[:,v][accMask],'bs')
        plt.ylabel(varList[v])  
    plt.tight_layout()
    
    # Save figure
    plt.savefig(txtFiles_active + '/record_trace.eps')
        
    
def SSErestrict_getOK(obsBA_mc, estBA_mc, args, txtFiles_active):
        
    SSE_mc = np.power(estBA_mc - obsBA_mc,2)
    SSE_sumM = np.squeeze(np.sum(SSE_mc,axis=0))
    okcells = np.squeeze(np.nonzero(SSE_sumM <= np.percentile(SSE_sumM,args.SSEpctile)))   # np.nonzero() finds indices where TRUE
    
    # Save for reference in next run
    np.save(txtFiles_active + 'okcells',okcells)
    
    return okcells
    
    
def SSErestrict_doIt(obsBA_mc, obsBA, estBA_mc, estBA, derivs_mcv, derivs_xv, okcells, Ncells, args):
    if args.verbose>0:
        print('Restricting to cells with total SSE <= ' + str(args.SSEpctile) + 'th percentile.')
        
    # Restrict
    Nmonths = np.shape(obsBA_mc)[0]
    obsBA_mc = obsBA_mc[:,okcells]
    obsBA = np.reshape(obsBA_mc,Ncells*Nmonths)
    estBA_mc = estBA_mc[:,okcells]
    estBA = np.reshape(estBA_mc,Ncells*Nmonths)
    if len(np.shape(derivs_mcv))==3:
        derivs_mcv = derivs_mcv[:,okcells,:]
        derivs_xv = np.reshape(derivs_mcv,(Ncells*Nmonths,np.shape(derivs_mcv)[2]))
    elif len(np.shape(derivs_mcv))==2:
        derivs_mcv = derivs_mcv[:,okcells]
        derivs_xv = np.reshape(derivs_mcv,Ncells*Nmonths)
    
    return obsBA_mc, obsBA, estBA_mc, estBA, derivs_mcv, derivs_xv
            
            
def write_status(msg,txtFiles_active):
    with open(txtFiles_active + '/status.txt','ab') as f:
        f.write('        ' + time.strftime("%H:%M:%S") + '  ' + msg + '\n')           
            
