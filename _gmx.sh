echo
echo "# GROMACS "
echo "#========="

# Check and set the gromacs related stuff  
GMXRC=$(find_program_function gmxrc)
echo "# Gromacs RC file: $GMXRC"

# Source the gromacs RC file if one was found
# Otherwise the script will rely on the active gromacs commands 
[[ $GMXRC ]] && source $GMXRC 

# Find out which Gromacs version this is
GMXVERSION=$(mdrun -h 2>&1 | $SED -n '/^.*VERSION \([^ ]*\).*$/{s//\1/p;q;}')
# From version 5.0.x on, the commands are gathered in one 'gmx' program
# The original commands are aliased, but there is no guarantee they will always remain
[[ -z $GMXVERSION ]] && GMXVERSION=$(gmx -h 2>&1 | $SED -n '/^.*VERSION \([^ ]*\).*$/{s//\1/p;q;}')
# Version 2016 uses lower case "version", which is potentially ambiguous, so match carefully
[[ -z $GMXVERSION ]] && GMXVERSION=$(gmx -h 2>&1 | $SED -n '/^GROMACS:.*gmx, version \([^ ]*\).*$/{s//\1/p;q;}')
ifs=$IFS; IFS="."; GMXVERSION=($GMXVERSION); IFS=$ifs

# Set the directory for binaries
[[ $GMXVERSION -gt 4 ]] && GMXBIN=$(which gmx) || GMXBIN=$(which mdrun)
# Extract the directory
GMXBIN=${GMXBIN%/*}
# Set the directory to SCRIPTDIR if GMXBIN is empty 
GMXBIN=${GMXBIN:-$SCRIPTDIR}

# Set the command prefix and set the GMXLIB variable to point to 
# the force field data and such.
# In some cases, 'gromacs' is part of $GMXDATA
if [[ $GMXVERSION -gt 4 ]]
then
  GMX="$GMXBIN/gmx " 
  GMXLIB=
else
  GMX=$GMXBIN/
  export GMXLIB=${GMXDATA}/gromacs/top
  [[ -d $GMXLIB ]] || export GMXLIB=${GMXDATA%/gromacs*}/gromacs/top
fi

# Now finally, test a command and see if it works
# otherwise raise a fatal error.
${GMX}grompp -h >/dev/null 2>&1 || executable_not_found_error "GROMACS (GMXRC)"

# Check if Gromacs can handle RTC (if so needed)
# This only works for GMX <5
if [[ $RotationalConstraints == "rtc" ]]
then
  echo Testing rotational constraints
  echo "comm_mode = RTC" > gromacs_rtc_test.mdp
  if ${GMX}grompp -f gromacs_rtc_test.mdp 2>&1 | grep -q "Invalid enum 'RTC'"
  then
    rm gromacs_rtc_test.mdp
    FATAL "Roto-translational constraints requested (comm_mode = RTC), but not supported by ${GMX}grompp"
  fi
  rm gromacs_rtc_test.mdp mdout.mdp
fi

echo "# - binaries:  $GMXBIN"
echo "# - data:      $GMXDATA"
echo "# - libraries: $GMXLIB"
echo 
