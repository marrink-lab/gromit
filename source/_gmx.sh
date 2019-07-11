function load_gromacs() {    
    ## Contents:
    ##
    ##  I. Determining the gromacs version
    ## II. Expressions for sed and awk to handle topologies and structures
    ##
    
    
    #------------------------------#
    # I. Determine gromacs version #
    #------------------------------#
    
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
    GMXBIN=${GMXBIN:-$SDIR}
    
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
}
    
    
#---------------------------------------------------------------------#
# II. Expressions for sed and awk to handle topologies and structures #
#---------------------------------------------------------------------#

# Extract the moleculetype names from an itp file (or only the first)
#    - at the line matching 'moleculetype' 
#      read in the next line
#      continue reading next lines until one is not beginning with ;
#      print the first field
awk_itp_get_moltype='/moleculetype/{getline; while ($0 ~ /^ *;/) getline; print $1; exit}'
awk_itp_get_moltypes='/moleculetype/{getline; while ($0 ~ /^ *;/) getline; print $1}'
# Usage: MTP=($(for i in ${ITP[@]}; do awk "${awk_itp_get_moltypes}" $i; done))
    
    
# Get the current list of molecules from a top file 
#    and find the corresponding itp files
#    Some of these may have been changed to remove duplicate listings
#    This list should match the list of chains
#    awk expression
#    - if M is nonzero and the line does not start with ';', print
#    - at the line with the tag [ molecules ] set M to 1
#   Following each molecule name is the number of instances
awk_top_get_moleculelist='(M && ! /^ *;/); /\[ *molecules *\]/{M=1}'
# Usage: MOLECULES=($(awk "${awk_top_get_moleculelist}" $base-aa.top))
    
    
# Get the list of (itp) files #included in a top file
#    - at lines starting with #include and not containing a slash in the file name
#      substitute the line with the name of the included file and print
sed_top_get_includes='/^#include .*"\([^/]*\)"/s//\1/p'
# Usage: ITPFILES=($($SED -n -e "${sed_top_get_includes}" $base-aa.top))
