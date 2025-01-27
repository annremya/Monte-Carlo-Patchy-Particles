#!/bin/bash
#@ output     = test.out
#@ error      = test.err
#@ job_type   = serial
#@ class      = Long
#@ environment = COPY_ALL
#@ queue
Jobid=`echo $LOADL_STEP_ID | cut -f 6 -d .`
tmpdir=$HOME/scratch/50_p0.01_2_$Jobid
mkdir -p $tmpdir; cd $tmpdir
cp -R $LOADL_STEP_INITDIR/* $tmpdir
time ./patchy
mv ../50_p0.01_2_$Jobid $LOADL_STEP_INITDIR 
