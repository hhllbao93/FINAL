#!/bin/bash
#PBS -d /lustre/f1/unswept/Sam.Rabin/param_est/testing_12month_parameterization_20150511
#PBS -o /lustre/f1/unswept/Sam.Rabin/param_est/testing_12month_parameterization_20150511_results
#PBS -N main_job_script
#PBS -l walltime=57600
#PBS -W umask=026
#PBS -S /bin/bash
#PBS -r y
#PBS -q batch
#PBS -l partition=c1:c2
#PBS -m abe
#PBS -j oe
#PBS -E
#PBS -l size=250
#PBS -l qos=norm
#PBS -A gfdl_b

set -e # quit on errors
set -u # quit on undefined variables

### What directory are we in?
baseDir=/lustre/f1/unswept/Sam.Rabin/param_est/testing_12month_parameterization_20150511
cd $baseDir

### What directories for restart data?
restartDir_1=/lustre/f1/unswept/Sam.Rabin/param_est/initConds__tikal_lad2_20141008_slm_ulmFix/working_dirs/spinup_CORPSE_250_hill_v2/19811990
restartDir_2=/lustre/f1/unswept/Sam.Rabin/param_est/spinup_run/1948-1990

### Which parameters are we running, and what are their starting guesses?
agb_gom2=7.3157
agb_gom3=4.1100
alphaM=0.0035
popD_supp_eps3=0.025
rh_gom2=0.0062
rh_gom3=-9.1912
theta_gom2=0.0750
theta_gom3=-6.3741
ROSmax_tt=0.3
ROSmax_c4=0.4

# define number of processors: if PBS_NP is undefined, use 1 processor
npes=${PBS_NP:-1}

# initialize environment modules
source /opt/modules/3.2.6.6/init/bash
module use -a /ncrc/home2/fms/local/modulefiles
module unload PrgEnv-pgi PrgEnv-pathscale PrgEnv-intel PrgEnv-gnu PrgEnv-cray
module unload netcdf fre fre-commands
module load PrgEnv-intel
#module swap intel intel/12.0.5.220
module load fre/bronx-8
module load python
module load python_netcdf4
module load python_matplotlib/1.4.3
module list

# Get the name of this script
#mainScript=`basename "$0"`
mainScript=${PBS_JOBNAME##*/}

# Get the name of this set of runs
thisRun=$mainScript.`date +%Y%m%d_%H%M%S`

# I don't know what these vars do, this is pasted from the standard FRE runscript
export F_UFMTENDIAN=big
export KMP_STACKSIZE=2g
export MPICH_ENV_DISPLAY 
export MPICH_CPUMASK_DISPLAY
export MPICH_GNI_LOCAL_CQ_SIZE=131072

# input
namelists="$baseDir/input_1991-2009_CORPSE_hill.nml"
diagTable="$baseDir/diagTable_fireMin_byLC_fireDur"
siteTable="$baseDir/points.txt"
siteDir="$baseDir/points"

# path to executables
executable="$baseDir/code_this_branch/exec/fms_code.x"
  makeGrid="/autofs/na1_home1/Sergey.Malyshev/bin/make-1x1-gridspec"
makeRivers="/autofs/na1_home1/Sergey.Malyshev/bin/make-dummy-rivers"
taskfarmer="/lustre/f1/unswept/Sam.Rabin/slm_for_ssr/taskfarmer/src/taskfarmer"

# We assemble all data into the dataDir, and then hard-link it for each of the points.
# We do this instead of soft-linking to the files directly for each point because that
# stopped working for some reason. SSR modified this so that the data only need to 
# get copied once for all runs in this directory.
dataDir="$baseDir/data"

# Prepare directory structure in dataDir, parameterization_work, and outDir_root
mkdir -p $dataDir/{INPUT,RESTART} || exit 1
param_workDir=$baseDir/paramwork.$thisRun
mkdir -p $param_workDir/txtFiles_orig || exit 1
txtFiles_orig_dir=$param_workDir/txtFiles_orig
outDir_root=$baseDir/working_dirs/output.$thisRun
mkdir $outDir_root

# Copy Python scripts to parameterization working directory
#cp $baseDir/python_scripts/*.py $param_workDir
cp $baseDir/python_scripts_v12/*.py $param_workDir

# Copy input data
. ./do_copyData.sh

# prepare forcing data
ls -1 /lustre/f1/unswept/Sergey.Malyshev/DATA/Sheffield/2013-03-19/[0-9][0-9][0-9][0-9].nc > $dataDir/file_table

# copy field table
cat > $dataDir/field_table <<EOF
# co2
 "TRACER", "atmos_mod",    "co2"
           "longname",     "carbon dioxide"
           "units",        "kg/kg"
           "convection",   "all"
           "profile_type","fixed","surface_value = 434.5626E-06"/
 "TRACER", "land_mod",     "co2"
           "longname",     "carbon dioxide"
           "units",        "kg/kg" /
# specific humidity for moist runs
 "TRACER", "atmos_mod", "sphum"
           "longname",     "specific humidity"
           "units",        "kg/kg"
       "profile_type", "fixed", "surface_value=3.e-6" /
 "TRACER", "land_mod", "sphum"
           "longname",     "specific humidity"
           "units",        "kg/kg" /
EOF

# create data table
### For historical runs (1860-), will need to change as I did for XMLs.
cat > $dataDir/data_table <<EOF
"ATM", "co2_bot", "co2", "INPUT/co2_gblannualdata.nc", .false., 1.0e-6
EOF

# diag table
cat $diagTable > $dataDir/diag_table 

# Which parameters are we running? These should be named after the suffixes to
# the BA_derivWRT_* variables corresponding to the parameters being optimized.
# (Python scripts will use these for importing.)
cat > $txtFiles_orig_dir/parameters_present.txt <<EOF
AGBparam1
AGBparam2
alphaM
popDsupp_NF_eps3
RHparam1
RHparam2
THETAparam1
THETAparam2
ROSmax_tt
ROSmax_c4
EOF

# Make starting version of params_run.txt
touch $txtFiles_orig_dir/params_run_orig.txt
cat >> $txtFiles_orig_dir/params_run_orig.txt <<EOF
$agb_gom2 $agb_gom3 $alphaM $popD_supp_eps3 $rh_gom2 $rh_gom3 $theta_gom2 $theta_gom3 $ROSmax_tt $ROSmax_c4
EOF

# Start incomplete namelist
#####cat $namelists > $dataDir/input.nml
#####cp $dataDir/input.nml $txtFiles_orig_dir/input_incomplete.nml
cat $namelists > $txtFiles_orig_dir/input_incomplete.nml
cp $txtFiles_orig_dir/input_incomplete.nml $txtFiles_orig_dir/input.nml

# Add basic fire stuff to incomplete namelist
cat > $txtFiles_orig_dir/fire_incomplete.nml <<EOF

   &fire_nml
      fire_to_use = 'unpacked'
      do_calc_derivs = .TRUE.,
      f_agb_style = 'gompertz'
      f_rh_style = 'gompertz'
      f_theta_style = 'gompertz'
      min_BA_to_split = 1.0
      do_fire_fragmentation = .TRUE.
      frag_incl_nonveg = .TRUE.
      frag_incl_PAST = .FALSE.
      zero_TMF = .TRUE.
      use_FpopD_ba = .FALSE.
      use_Fgdp_nf = .FALSE.
      use_Fgdp_ba = .FALSE.
      theta_ROSeffect_asFnTheta = .TRUE.
      C_beta_params_likeRH = .FALSE.
      lock_ROSmaxC3_to_ROSmaxC4 = .TRUE.
EOF
#####cat $txtFiles_orig_dir/fire_incomplete.nml >> $dataDir/input.nml
cat $txtFiles_orig_dir/fire_incomplete.nml >> $txtFiles_orig_dir/input.nml
#####cat >> $dataDir/input.nml <<EOF
cat >> $txtFiles_orig_dir/input.nml <<EOF
      agb_gom2 = $agb_gom2
      agb_gom3 = $agb_gom3
      Ia_alpha_monthly = $alphaM
      popD_supp_eps3 = $popD_supp_eps3
      rh_gom2 = $rh_gom2
      rh_gom3 = $rh_gom3
      theta_gom2 = $theta_gom2
      theta_gom3 = $theta_gom3
      ROS_max_TROPICAL = $ROSmax_tt
      ROS_max_C4GRASS  = $ROSmax_c4
/
EOF

####################################
### First model run and Script 1 ###
####################################

# Make new txtFiles_active directory
txtFiles_active_dir=$param_workDir/txtFiles_active
mkdir $txtFiles_active_dir
# Put starting values of things in there
cp $txtFiles_orig_dir/params_run_orig.txt $txtFiles_active_dir/params_run.txt
cp $txtFiles_orig_dir/parameters_present.txt $txtFiles_active_dir
touch $txtFiles_active_dir/status.txt

Nruns=1
echo " "
echo "########## MODEL RUN ##########"
date
echo "###############################"
echo " "
date >> $txtFiles_active_dir/status.txt
echo "|___Model run $Nruns..." >> $txtFiles_active_dir/status.txt
. ./doModel.sh
rm $executable_tmp

echo " "
echo "########## SCRIPT 1 ##########"
date
echo "##############################"
echo " "
date >> $txtFiles_active_dir/status.txt
echo "|___Script 1..." >> $txtFiles_active_dir/status.txt
echo $param_workDir/LMv1_script1.py
python $param_workDir/LMv1_script1.py --SSEpctile 100.0 -n 9 --constrain_thetaEsq --constrain_ROSmax --constrain_Ia --constrain_duration --constrain_AGBparams --constrain_popDsupp_NF --constrain_RHparams --damping_style transtrum_minVal --f_rh_style gompertz --f_agb_style gompertz --f_theta_style gompertz -v 1 -p $param_workDir -o $outDir_root -O 1e-13 1e-15 1e-15 100 --newest_unpacked

# Check bash_do
string=$(tail -n 1 $txtFiles_active_dir/bash_do.txt)
IFS=' ' read -ra bash_do_array <<< "$string"
bash_do=${bash_do_array[0]}

##############################
### Subsequent model runs  ###
##############################

while [ "$bash_do" -gt 0 ]
do

   # Run model
   let Nruns=$Nruns+1
   echo " "
   echo "########## MODEL RUN $Nruns ##########"
   date
   echo "######################################"
   echo " "
   date >> $txtFiles_active_dir/status.txt
   echo "|___Model run $Nruns..." >> $txtFiles_active_dir/status.txt
   . ./doModel.sh   # The "." at the beginning tells bash to run the script in this shell.
   rm $executable_tmp
   
   # Run Script 2
   echo " "
   echo "########## SCRIPT 2 ##########"
   date
   echo "##############################"
   date >> $txtFiles_active_dir/status.txt
   echo "|___Script 2..." >> $txtFiles_active_dir/status.txt
   python $param_workDir/LMv1_script2.py --SSEpctile 100.0 --constrain_thetaEsq -n 9 --constrain_ROSmax --constrain_Ia --constrain_duration --constrain_AGBparams --constrain_popDsupp_NF --constrain_RHparams --damping_style transtrum_minVal --f_rh_style gompertz --f_agb_style gompertz --f_theta_style gompertz -v 1 -p $param_workDir -o $outDir_root -O 1e-13 1e-15 1e-15 100 --newest_unpacked
   
   # Check bash_do
   string=$(tail -n 1 $txtFiles_active_dir/bash_do.txt)
   IFS=' ' read -ra bash_do_array <<< "$string"
   bash_do=${bash_do_array[0]}
   if [ "$bash_do" -eq 2 ] ; then
       echo "    SSE increased :-(" >> $txtFiles_active_dir/status.txt
   elif [ "$bash_do" -eq 1 ] ; then
       echo "    SSE decreased! :-D" >> $txtFiles_active_dir/status.txt
   fi 

done

##############
### Finish ###
##############

echo " "
echo " " >> $txtFiles_active_dir/status.txt
if   [ "$bash_do" -eq -1 ]; then
   echo "########## END: NG <= OPTS[1] ##########"
   echo "########## END: NG <= OPTS[1] ##########" >> $txtFiles_active_dir/status.txt
elif [ "$bash_do" -eq -2 ]; then
   echo "########## END: PARAMETER STEP TOO SMALL ##########" 
   echo "########## END: PARAMETER STEP TOO SMALL ##########" >> $txtFiles_active_dir/status.txt 
elif [ "$bash_do" -eq -3 ]; then
   echo "########## END: TOO MANY ITERATIONS ##########" 
   echo "########## END: TOO MANY ITERATIONS ##########" >> $txtFiles_active_dir/status.txt 
elif [ "$bash_do" -eq -4 ]; then
   echo "########## END: (ALMOST) SINGULAR ##########" 
   echo "########## END: (ALMOST) SINGULAR ##########" >> $txtFiles_active_dir/status.txt 
fi



exit 0
