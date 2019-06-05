# This the section that is run in martinate.sh if DAFT is set


# This section is pretty much obsolete, as the DAFT workflow 
# invokes martinate with a complete CG structure and topology

echo Splitting molecules in energy groups for DAFT stuff

# Make a new energy group index file
# It will be added to the master index later on
echo > EnergyGroups.ndx

# To make every molecule (so far) its own energygrp we have to

#  a. Identify molecules (reusing code from above, STEP 1A, section IV.6.f-i)
#  b. Make an index group for it
#  c. List it
MOLECULES=($(awk "${awk_top_get_moleculelist}" $base-cg.top))
echo MOLECULES: ${MOLECULES[@]}

ITPFILES=($($SED -n -e "${sed_top_get_includes}" $base-cg.top))

# We only take the FIRST moleculetype from the itp file
MOLTYPES=($(for i in  ${ITPFILES[@]}; do awk "${awk_itp_get_moltype}" $i; done))

MOLITP=()
COUNT=1
K=1
EnergyGroups=()
echo "Molecules and corresponding topology files:"
for ((i=0; i<${#MOLECULES[@]}; i+=2))
do
    for ((j=0; j<${#MOLTYPES[@]}; j++))
    do
        if [[ ${MOLECULES[$i]} == ${MOLTYPES[$j]} ]]
        then
            MOLITP[$((i/2))]=${ITPFILES[$j]}
            # Get the number of atoms for this moleculetype
            N=$(awk '/^ *\[ *atoms *\]/,/^ *$/{if ($1 ~ /^ *[0-9]/) N++}END{print N}' ${ITPFILES[$j]})
            # For each molecule of this kind, write an index group
            for ((k=0; k<${MOLECULES[$((i+1))]}; k++))
            do
                EnergyGroups[${#EnergyGroups[@]}]=${MOLECULES[$i]}_$K
                echo -e "[ ${MOLECULES[$i]}_$K ]\n$(printf "$fmt\n" `SEQ $COUNT $((COUNT+N-1))` | $SED 's/ 0//g')" >> EnergyGroups.ndx
                : $((COUNT+=N)) $((K++))
            done
        fi
    done
    printf "    %20s : %20s\n" ${MOLECULES[$i]}: ${MOLITP[$((i/2))]}
done
__mdp_cg__energygrps=$($SED 's/ /,/g' <<< "${EnergyGroups[@]} Membrane Solvent")
echo @@ $__mdp_cg__energygrps
