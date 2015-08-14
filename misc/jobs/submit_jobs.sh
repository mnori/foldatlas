#!/bin/bash

# Create all the job scripts, same folder as this one
for i in `seq 0 257`;
do
        bsub -q NBI-Test128 -R"rusage[mem=4096]" "./job_$i.sh"
done
echo "Jobs submitted!"
