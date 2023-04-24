#!/bin/bash

# --- Check if the Hough-Transform training file already exists:
if [ -f "HTAdata.txt" ]
then
    echo "Removing HTAdata.txt ..."
    /bin/rm HTAdata.txt
fi

# --- Recompile the executables
echo "Recompiling the executables ..."
make

echo "Starting the Hough-Transform array training."

for phase in {1..3};
do

    echo -n "Training PHASE $phase ... "
    ./HTArrayTraining >& out_HTtraining_$phase.log
    echo "DONE"

done


