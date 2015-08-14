#!/bin/bash

# Create all the job scripts, same folder as this one
for i in `seq 0 257`;
do
        echo "CHUNK=$i" > tmp1.sh
        cat tmp1.sh _job_tpl.sh > job_$i.sh
        rm tmp1.sh
done

chmod 777 job_*
echo "Jobs created"
