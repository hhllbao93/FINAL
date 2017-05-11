# ################################
# ### LM for LM3, v1: Script 2 ###
# ################################

# This script comes after all model runs other than the first. It first
# calculates the difference in the sum of squared errors between the
# latest model run and the last accepted model run. If the former is
# less than the latter, this script saves the latest parameter set guess
# as "accepted," reduces mu (and resets nu), and generates the next guess.
# Otherwise, this script saves the latest parameter guess as "rejected,"
# increases mu (and nu), and generates the next guess (based on the last
# ACCEPTED parameter set).

## Set up and import data

import os
import glob
# Change directory to one with scripts
if os.getcwd()[:7] == '/lustre' or os.getcwd()[:7] == '/lustre/':
    # I.e., on Gaea
    location = 'gaea'
    thisDir = os.getcwd()
    if thisDir[len(thisDir)-1] != '/':
        thisDir = thisDir + '/'
    
    # Change (temporarily) to python_scripts
    os.chdir('python_scripts_v12/')
else:
    raise NameError('What platform are you on??')
execfile('./LMv1_setup_import.py')


## Calculate f, dF, etc.

# Read previous ACCEPTED and GUESSED parameter sets from file, and get Niters
ssr.write_status('Reading previous param sets...',txtFiles_active)
Niters = ssr.get_Niters(txtFiles_active)
x_lastRun = ssr.readLatest_params_run(txtFiles_active,Nvars)
x_lastAcc = ssr.readLatest_params_acc(txtFiles_active,Nvars,Niters)

if args.verbose>0:
    print('Niters = ' + str(Niters))
    print('x_lastRun = ' + str(x_lastRun))
    print('x_lastAcc = ' + str(x_lastAcc))

# Calculate f, SSEhalf_new
ssr.write_status('Calculating new run SSE...',txtFiles_active)
f, SSEhalf_new = ssr.calcSSE(estBA,obsBA)

# Save SSEhalf and parameters from this run to SSEhalf_record.txt
ssr.write_status('Saving SSEhalf and params from this run...',txtFiles_active)
SSEhalf_record = open(txtFiles_active + '/SSEhalf_record.txt','a')
writer = csv.writer(SSEhalf_record, delimiter=' ', quoting=csv.QUOTE_NONE)
row2write = [str(SSEhalf_new)] ;
[row2write.append(i) for i in map(str,x_lastRun)]
row2write.append(dt.now().strftime("%Y%m%d_%H%M%S_%f"))
writer.writerow(row2write)
SSEhalf_record.close() ; del writer ; del row2write

# Update trace of SSE and parameters
ssr.write_status('Making trace plot...',txtFiles_active)
ssr.trace_results(txtFiles_active,Nvars,varList)

# Read SSEhalf of most recently accepted parameter set from file.
ssr.write_status('Reading previous accepted SSEhalf...',txtFiles_active)
SSEhalf_old = np.load(txtFiles_active + 'lastAcc_SSEhalf.npy')

# Calculate dF
dF = SSEhalf_old - SSEhalf_new
if args.verbose>0:
    print('SSEhalf_old = ' + str(SSEhalf_old))
    print('SSEhalf_new = ' + str(SSEhalf_new))
    print('dF = ' + str(dF))


## Revise mu & nu in preparation for calculating next guess

# Read latest mu and nu from file
ssr.write_status('Reading latest mu & nu...',txtFiles_active)
mu, nu = ssr.readLatest_muNu(txtFiles_active)

if dF > 0:   # I.e., SSEhalf did improve with new parameters.

    # Read dL (from BEFORE the model run) from file.
    ssr.write_status('Loading pre-run dL...',txtFiles_active)
    # Needed to reduce mu.
    dL_preModelRun = np.load(txtFiles_active + 'lastRun_dL.npy')
    
    # Reset mu and nu for calculation of next parameter set
    ssr.write_status('Resetting mu & nu...',txtFiles_active)
    mu = mu * max(1./3, 1. - (2.*dF/dL_preModelRun - 1.)**3)
    nu = 2.0
    if args.verbose>0:
        print('SSE improved with new parameters. New mu & nu:')
        print('mu = ' + str(mu) + ', nu = ' + str(nu))
        
    # Save SSEhalf and parameters from this run to SSEhalf_acc.txt
    ssr.write_status('Saving things from this run...',txtFiles_active)
    SSEhalf_acc = open(txtFiles_active + '/SSEhalf_acc.txt','a')
    writer = csv.writer(SSEhalf_acc, delimiter=' ', quoting=csv.QUOTE_NONE)
    row2write = [str(SSEhalf_new)] ;
    [row2write.append(i) for i in map(str,x_lastRun)]
    row2write.append(dt.now().strftime("%Y%m%d_%H%M%S_%f"))
    writer.writerow(row2write)
    SSEhalf_acc.close() ; del writer ; del row2write
    np.save(txtFiles_active + 'lastAcc_SSEhalf',SSEhalf_new)
    
    # Save these parameter values as "accepted," along with dF
    params_acc = open(txtFiles_active + '/params_acc.txt','a')
    writer = csv.writer(params_acc, delimiter=' ', quoting=csv.QUOTE_NONE)
    row2write = map(str,x_lastRun)
    row2write.append(dF)
    row2write.append(dt.now().strftime("%Y%m%d_%H%M%S_%f"))
    writer.writerow(row2write)
    params_acc.close() ; del writer ; del row2write
    np.save(txtFiles_active + 'lastAcc_params',x_lastRun)
    
    # Calculate J from this run
    J = derivs_xv
    
    # Save these f and J arrays
    ssr.write_status('Saving f & J...',txtFiles_active)
    ssr.save_fJ(f,J,txtFiles_active,args)
    
    # Replace x_lastAcc with x_lastRun; delete x_lastRun to avoid any possible confusion
    x_lastAcc = x_lastRun
    del x_lastRun
    
else:   # I.e., SSEhalf didn't improve with new parameters
    
    # Load f and J from run with last accepted parameter set
    ssr.write_status('Reading last accepted f and J...',txtFiles_active)
    del f
    f, J = ssr.readLatest_fJ_acc(txtFiles_active,args)
    
    # Increase mu and nu
    mu = mu * nu
    nu = nu * 2.0
    if args.verbose>0:
        print('SSE did not improve with new parameters. New mu & nu:')
        print('   mu = ' + str(mu))
        print('   nu = ' + str(nu))
    
    
## Calculate next guess
    
# Calculate g(f,J)
ssr.write_status('Calculating g...',txtFiles_active)
g = np.dot( np.transpose(J) , f )
if args.verbose==2:
    print('g = \n' + str(g))
if np.any(np.isnan(g)):
    raise NameError('At least one element of g is NaN.')
    
# Check for stopping based on ng. (Triggered if close to a local minimum.)
ssr.write_status('Checking for stop based on ng...',txtFiles_active)
stop_ng = ssr.checkStop_ng(g,opts,Niters,txtFiles_active)

if not stop_ng:

    # Calculate h, xnew; save new dL
    A, A_damped, h, xnew, dL, mu, nu, stop_h = ssr.calcH(J,g,mu,nu,x_lastAcc,opts,Niters,txtFiles_active,args,varList)
        
    if not stop_h:
                                                        
        # Calculate bins, if necessary
        if args.doBins:
            ssr.write_status('Calculating new bins...',txtFiles_active)
            if 'AGBparam1' in varList or 'AGBparam2' in varList:
                newBins = ssr.calcBins(bins_nonAGB,xnew[varList.index('AGBparam1')],xnew[varList.index('AGBparam2')])
            else:
                newBins = bins
        else:
            newBins = []
                            
        # Save next-guess parameter set for scripts and model
        ssr.write_status('Saving next-guess parameter set...',txtFiles_active)
        ssr.save_nextGuess(xnew,txtFiles_active,varList,newBins,args)
        
        # Save mu and nu for next call of Script 2
        ssr.write_status('Saving mu & nu...',txtFiles_active)
        ssr.save_muNu(mu,nu,Niters,txtFiles_active)
        
        # Tell the model we're not at a stop condition.
        ssr.write_status('Saving bashDo...',txtFiles_active)
        if dF > 0:
            ssr.save_bashDo([str(1), np.nan, np.nan, dt.now().strftime("%Y%m%d_%H%M%S_%f")],Niters,txtFiles_active)
        else:
            ssr.save_bashDo([str(2), np.nan, np.nan, dt.now().strftime("%Y%m%d_%H%M%S_%f")],Niters,txtFiles_active)
        
        # Save record of txtFiles_active at this point
        if location == 'gaea':
            ssr.save_txtFiles_active(runDir,working_dirs,txtFiles_active)


if args.verbose>0:
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('%%% End Script2 verbose output %%%')
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

ssr.write_status('All done?...',txtFiles_active)