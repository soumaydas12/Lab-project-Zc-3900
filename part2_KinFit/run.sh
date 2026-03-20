#!/bin/bash

#=============================================================================
# This function called "run" does:
# 1. open root with "root -l"
# 2. compile ZcAnalysis.C with ".L"
# 3. create a ZcAnalysis object "ana"
#    - the first argument "$1" is the path to the data
#    - the second argument "$2" is the save path for the plots
#    - the "Loop()"" function automatically runs
# 4. quit root with ".q"
#=============================================================================

run () {
root -l << INPUT
    .L ZcAnalysis.C
    ZcAnalysis* ana = new ZcAnalysis("$1", "$2")
    .q
INPUT
}

#=============================================================================
# The next lines execute the "run" function from above with two arguments.
# You can run the analyses with different paths in parallel with "&".
# Comment out the ones you don't need!
#=============================================================================

run "../NTuples/mc_pippimJpsi.root" "plots_pippimJpsi/" &
run "../NTuples/data_4260.root" "plots_data/" &
run "../NTuples/mc_pipmZc3900.root" "plots_pipmZc3900/" &

#=============================================================================
# Wait for all "run" functions to finish
#=============================================================================

wait
