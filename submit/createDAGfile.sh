#!/bin/sh

s='nominal'
LLH='effective'

LEmasses=(40 63 100 158 251 398 631 1000 1585)
HEmasses=(1000 1585 2512 3981 6310 10000 15850 25120 39810)

DM_processes=('annihilation' 'decay')
DM_profiles=('NFW' 'Burkert')
FINAL_STATES=('nue' 'b' 'W' 'tau' 'mu')

for c in "${FINAL_STATES[@]}"; do

    for p in "${DM_profiles[@]}"; do
        
        for t in "${DM_processes[@]}"; do
        
            for m in "${LEmasses[@]}"; do

                if [ $c == 'nue' ]; then
                    oversampling='100'
                else
                    oversampling='-1'
                fi

                JOBIDSENS=$c-$m-$p-$t-LE-Sens
                echo JOB $JOBIDSENS Sensitivity.submit
                echo VARS $JOBIDSENS JOBNAME=\"$JOBIDSENS\" TYPE=\"$t\" CHANNEL=\"$c\" PROFILE=\"$p\" SYST=\"$s\" LECUT=\"0.15\" HECUT=\"0.2\" MASS=\"$m\" OVERSAMPLING=\"$oversampling\" REBINE=\"2\" REBINPSI=\"5\"  CL=\"90\" LLH=\"$LLH\"
            done

            for m in "${HEmasses[@]}"; do

                if [ $c == 'nue' ]; then
                    oversampling='200'
                else
                    oversampling='-1'
                fi
                
                JOBIDSENS=$c-$m-$p-$t-HE-Sens
                echo JOB $JOBIDSENS Sensitivity.submit
                echo VARS $JOBIDSENS JOBNAME=\"$JOBIDSENS\" TYPE=\"$t\" CHANNEL=\"$c\" PROFILE=\"$p\" SYST=\"$s\" LECUT=\"-1.0\" HECUT=\"0.3\" MASS=\"$m\" OVERSAMPLING=\"$oversampling\" REBINE=\"2\" REBINPSI=\"5\"  CL=\"90\" LLH=\"$l\"

            done
        done
    done
done