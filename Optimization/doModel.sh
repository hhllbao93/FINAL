#!/bin/bash

set -e # quit on errors
set -u # quit on undefined variables

# define working directory name
outDir_this=$outDir_root/`date +%Y%m%d_%H%M%S`
mkdir $outDir_this

# Make a copy of the executable to avoid the problem (?) of multiple tasks
# hitting fms_code.x at the same time
executable_tmp=$executable.$thisRun.`date +%Y%m%d%H%M%S`
cp $executable $executable_tmp

# loop over points, clone dataDir for each point, and add a model run command line 
# to the task list file
### Uses input from $sitetable (see "done" line). read returns FALSE when no more lines left.
echo Setting up each point...
while read lon lat REST
do
   # convert geographical notation to coordinates
   ### Geog. notation used for nice dir names
   case $lon in
      *W) x=-${lon%?} ;;
      *E) x=${lon%?}  ;;  
      *)  x=${lon}    ;;  
   esac
   case $lat in
      *S) y=-${lat%?} ;;
      *N) y=${lat%?}  ;;  
      *)  y=${lat}    ;;  
   esac

   #create working directory for the point   
   name=${lon}x${lat}
   echo $name
   workDir=$outDir_this/$name
   mkdir $workDir
   # clone pre-created data directory for the points
   mkdir -p $workDir/{INPUT,RESTART}
   ln -s $dataDir/INPUT/* $workDir/INPUT/
   cp $dataDir/*table $workDir

   # Copy original or Python-modified namelist file
   if [ $Nruns -eq 1 ]; then
      cp $txtFiles_orig_dir/input.nml $workDir
   else
      cp $txtFiles_active_dir/input.nml $workDir
   fi

   # Copy restart data
   if [ ! -f $restartDir_1/$name/RESTART/land.res.nc ]; then
       if [ ! -f $restartDir_2/$name/RESTART/land.res.nc ]; then
           echo "No data in $restartDir_2/$name/RESTART/ !!!"
           exit 42
       else
            cp $restartDir_2/$name/RESTART/* $workDir/INPUT/
       fi

   else
       cp $restartDir_1/$name/RESTART/* $workDir/INPUT/
   fi

   # create grid spec and dummy river file    
   if [ ! -f $siteDir/$name/grid_spec.nc ]; then
       echo "^Making grid spec and dummy river file...^"
      mkdir -p $siteDir/$name
      $makeGrid -- $x $y $siteDir/$name/grid_spec.nc 1> /dev/null
      $makeRivers $siteDir/$name/grid_spec.nc $siteDir/$name/river_data.nc 1> /dev/null
   fi
   cp $siteDir/$name/* $workDir/INPUT/

   # build the list of tasks for taskfarmer
   echo "cd $workDir ; $executable_tmp 2>&1 >fms.out" >> $outDir_this/tasks.txt
done < $siteTable

echo $outDir_root
echo $outDir_this 

# run the multiple copies of the model in taskfarmer
# taskfarmer gradually consumes the task file, as the jobs get executed,
# so that at any time the taskfile contains the jobs that has not been started yet.
echo Done. Running model...
aprun -n $npes $taskfarmer -v -f $outDir_this/tasks.txt
