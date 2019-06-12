# Residue groups used for classifying atoms in the structure file.
# Ions are typically considered positioned after solvent.
# Membrane is the complementary group to the rest.
# The structure file is assumed to have the following composition:
#  - Protein
#  - Nucleic acids
#  - Membrane
#  - Solvent
#  - Ions
# All groups are optional (as long as there are some)
amino_acids=(ALA CYS ASP GLU PHE GLY HIS HIH ILE LYS LEU MET ASN PRO HYP GLN ARG SER THR VAL TRP TYR)
nucleic_acids=(DG DA DC DT G A C U)
solvent_names=(W WF PW BMW SOL HOH)


## Format for printing numbers to index files
fmt=" %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d"


# Amino acids
function sed_protein()
{
  for a in ${amino_acids[@]}
  do 
    echo '-e /^.....\s*'$a/$1
  done
}


# Nucleic acids
function sed_nucleic()
{
  for a in ${nucleic_acids[@]}
  do 
    echo '-e /^.....\s*'$a/$1
  done
}


# Solvent
function sed_solvent()
{
  for a in ${solvent_names[@]}
  do 
    echo '-e /^.....\s*'$a/$1
  done
}


function INDEX()
{
  [[ -n $2 ]] && fn=$2 || fn=basic.ndx
 
  exec 6>&1 && exec >$fn

  # Whatever we do not have explicitly must be membrane

  # Total number of atoms
  N=$(awk '{getline; print; exit}' $1)
  echo "[ System ]"
  printf "$fmt\n" `SEQ 1 $N` | $SED 's/ 0//g'

  # Solvent atoms (including ions, etc, listed after 'SOL/W/...')
  SOL=$(( $($SED -n $(sed_solvent '{=;q;}') $1) - 2 ))
  echo "[ Solvent ]"
  printf "$fmt\n" `SEQ $SOL $N` | $SED 's/ 0//g'
    
  # Base system: solute and membrane, if present
  echo "[ Base ]"
  printf "$fmt\n" `SEQ 1 $((SOL - 1))` | $SED 's/ 0//g'

  # Protein - Get the first and last atom of protein stuff
  PROTSTART=$(( $($SED -n $(sed_protein '{=;q;}') $1) - 2 ))
  PROTEND=$(( $($SED -n -e 1,2d $(sed_protein d) -e '{=;q;}' $1) - 3 )) 
  echo "[ Protein ]"
  [[ $PROTSTART -gt 0 ]] && printf "$fmt\n" `SEQ $PROTSTART $PROTEND` | $SED 's/ 0//g'

  # Protein - Get the first and last atom of protein stuff
  NUCSTART=$(( $($SED -n $(sed_nucleic '{=;q;}') $1) - 2 ))
  NUCEND=$(( $($SED -n -e 1,2d $(sed_nucleic d) -e '{=;q;}' $1) - 3 )) 
  echo "[ Nucleic ]"
  [[ $NUCSTART -gt 0 ]] && printf "$fmt\n" `SEQ $NUCSTART $NUCEND` | $SED 's/ 0//g'

  # Membrane, if any
  MEMBRANE=$(( $($SED -n -e '1,2d' $(sed_protein d) $(sed_nucleic d) -e '{=;q;}' $1) - 2 ))
  echo '[ Membrane ]'
  if [[ -n $MEMBRANE && $MEMBRANE -gt 0 && $MEMBRANE != $SOL ]]
  then
    printf "$fmt\n" `SEQ $MEMBRANE $((SOL - 1))` | $SED 's/ 0//g'
  fi

  echo '[ Solute ]'
  if [[ $PROTSTART -gt 0 && $NUCSTART -gt 0 ]]
  then
    printf "$fmt\n" `SEQ $PROTSTART $NUCEND` | $SED 's/ 0//g'
  elif [[ $PROTSTART -gt 0 ]]
  then
    printf "$fmt\n" `SEQ $PROTSTART $PROTEND` | $SED 's/ 0//g'
  elif [[ $NUCSTART -gt 0 ]]
  then
    printf "$fmt\n" `SEQ $NUCSTART $NUCEND` | $SED 's/ 0//g'
  fi

  exec 1>&6 6>&-

  return 0
}

