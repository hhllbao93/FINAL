# ################################
# ### LM for LM3, v1: Script 1 ###
# ################################

# This script only comes after the first model run. Its only purpose
# is to calculate the SSE of the first model run and then generate
# a new guess.

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


## Do copying/deleting as necessary, since it's the first run

ssr.write_status('Copying/deleting files...',txtFiles_active)
if args.verbose>0:
    print('Copying ' + txtFiles_orig + 'params_run_orig.txt to files ' + txtFiles_active + 'params_[run,acc].txt')
shutil.copyfile(txtFiles_orig + 'params_run_orig.txt',txtFiles_active + 'params_run.txt')
shutil.copyfile(txtFiles_orig + 'params_run_orig.txt',txtFiles_active + 'params_acc.txt')

for f in glob.glob(txtFiles_active + '*.npy'):
    os.remove(f)


## Read previous parameters and calculate SSE

# Calculate f, SSE
ssr.write_status('Calculating first run SSE...',txtFiles_active)
f, SSEhalf = ssr.calcSSE(estBA,obsBA)
if args.verbose>0:
    print('SSEhalf = ' + str(SSEhalf))
    
# Read original parameter set from file
ssr.write_status('Reading original param set...',txtFiles_active)
x_init = ssr.readLatest_params_acc(txtFiles_active,Nvars,1)
if args.verbose>0:
    print('x_init = ' + str(x_init))

# Save to various files
ssr.write_status('Saving to various files...',txtFiles_active)
SSEhalf_acc = open(txtFiles_active + '/SSEhalf_acc.txt','wb')
writer = csv.writer(SSEhalf_acc, delimiter=' ', quoting=csv.QUOTE_NONE)
row2write = [str(SSEhalf)] ;
[row2write.append(i) for i in map(str,x_init)]
row2write.append(dt.now().strftime("%Y%m%d_%H%M%S_%f"))
writer.writerow(row2write)
SSEhalf_acc.close() ; del writer
SSEhalf_record = open(txtFiles_active + '/SSEhalf_record.txt','wb')
writer = csv.writer(SSEhalf_record, delimiter=' ', quoting=csv.QUOTE_NONE)
writer.writerow(row2write)
SSEhalf_record.close() ; del writer
params_run = open(txtFiles_active + '/params_run.txt','wb')
writer = csv.writer(params_run, delimiter=' ', quoting=csv.QUOTE_NONE)
row2write = map(str,x_init)
row2write.append(dt.now().strftime("%Y%m%d_%H%M%S_%f"))
writer.writerow(row2write)
params_run.close() ; del writer
params_acc = open(txtFiles_active + '/params_acc.txt','wb')
writer = csv.writer(params_acc, delimiter=' ', quoting=csv.QUOTE_NONE)
writer.writerow(row2write)
params_acc.close() ; del writer
np.save(txtFiles_active + 'lastAcc_SSEhalf',SSEhalf)
np.save(txtFiles_active + 'lastAcc_params',x_init)


## Generate next parameter set guess

# Calculate J, g(f,J)
ssr.write_status('Calculating J and g...',txtFiles_active)
J = derivs_xv
g = np.dot( np.transpose(J) , f )
if args.verbose==2:
    print('g = \n' + str(g))
if np.any(np.isnan(g)):
    raise NameError('At least one element of g is NaN.')


# Check for stopping based on ng
stop_ng = ssr.checkStop_ng(g,opts,1,txtFiles_active)
if not stop_ng:
    
    # Initialize mu and nu
    ssr.write_status('Initializing mu and nu...',txtFiles_active)
    A = ssr.calcA(J)
    mu = opts[0] * max(np.diag(A))
    nu = 2.0
    if args.verbose>0:
        print('Starting mu = ' + str(mu))
        print('         nu = ' + str(nu))
        
    # Calculate h, xnew; save dL.
    A, A_damped, h, xnew, dL, mu, nu, stop_h = ssr.calcH(J,g,mu,nu,x_init,opts,1,txtFiles_active,args,varList)
    
    if not stop_h:
        
        # Calculate bins, if necessary
        if args.doBins:
            if 'AGBparam1' in varList or 'AGBparam2' in varList:
                newBins = ssr.calcBins(bins_nonAGB,xnew[varList.index('AGBparam1')],xnew[varList.index('AGBparam2')])
            else:
                newBins = bins
        else:
            newBins = []
                            
        # Save next-guess parameter set for scripts and model
        ssr.write_status('Saving next guess...',txtFiles_active)
        ssr.save_nextGuess(xnew,txtFiles_active,varList,newBins,args)
        
        # Save mu and nu for next call of Script 2
        ssr.write_status('Saving mu & nu...',txtFiles_active)
        ssr.save_muNu(mu,nu,1,txtFiles_active)
        
        # Save f and J for Script 2
        ssr.write_status('Saving f & J...',txtFiles_active)
        ssr.save_fJ(f,J,txtFiles_active,args)
        
        # Tell the model we're not at a stop condition.
        ssr.write_status('Saving bashDo...',txtFiles_active)
        ssr.save_bashDo([str(1), np.nan, np.nan, dt.now().strftime("%Y%m%d_%H%M%S_%f")],1,txtFiles_active)
        
        # Save record of txtFiles_active at this point
        if location == 'gaea':
            histDir = txtFiles_active + '../txtFiles_active_history/'
            os.makedirs(histDir)
            ssr.save_txtFiles_active(runDir,working_dirs,txtFiles_active)
        
        
        
if args.verbose>0:
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('%%% End Script1 verbose output %%%')
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

ssr.write_status('All done?',txtFiles_active)


