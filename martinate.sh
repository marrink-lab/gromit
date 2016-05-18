#!/bin/bash

PROGRAM=martinate.sh
VERSION=0.1
VERSTAG=devel-160502-0800-TAW
AUTHOR="Tsjerk A. Wassenaar, PhD"
YEAR="2016"
AFFILIATION="
University of Groningen
The Netherlands"

CMD="$0 $@"


DESCRIPTION=$(cat << __DESCRIPTION__
This is a convenience script to set up, equilibrate and run coarse 
grained (and multiscaled) systems, using the Martini force field. 
It is built as a wrapper around martinize.py and insane.py, allowing 
automated processing of membrane proteins, with full control of membrane 
composition. If no input file is provided, only a membrane is built.
Otherwise, if a protein is given and -L is not set, the protein is 
solvated and run, but if -L is set, a membrane is built in addition.
Options given that do not match an option in this script are passed to
martinize.py.
__DESCRIPTION__
)


# These will be looked for before running, and can be set from the cmdline, e.g.:
#    -gmxrc /usr/local/gromacs-5.1/bin/GMXRC
# If not set, the default name will be searched for in
#    1. the environment (if PROGEVAR is given)
#    2. the directory where this calling script (martinate) is located
#    3. the PATH 
DEPENDENCIES=( dssp  gmxrc  martinize     insane   )
PROGEXEC=(     dssp  GMXRC  martinize.py  insane.py)
PROGEVAR=(     DSSP  GMXRC)


# The Gromacs RC file for use on the WeNMR GRID
# The existence of this file is checked later to
# to see if this is a GRID run.
# Note that the GMXRC file can also be specified with -gmxrc
# This is just a convenience hack.
GRIDRC="${VO_ENMR_EU_SW_DIR}/BCBR/gromacs/4.5.3-rtc/bin/GMXRC.bash"
GRID=false


# Starting step 
STEPS=(AA CG SOLVENT EM NVT-PR NPT PREPRODUCTION TPR PRODUCTION ANALYSIS END)
get_step_fun() { for ((i=0; i<${#STEPS[@]}; i++)) do [[ ${STEPS[$i]} =~ ^$1 ]] && return $i; done; }
STEP=AA
STOP=PRODUCTION
# Macros to echo stuff only if the step is to be executed -- more about steps further down
RED='\x1b[1;31m'
YEL='\x1b[1;33m'
OFF='\x1b[0m'
LINES() { sed -e 's/^/ /;s/$/ /;h;s/./-/g;p;x;p;x;' <<< "$@"; }
SHOUT() { [[ $STEP == $NOW ]] && echo && LINES "$@" | sed 's/^/#/'; }
WARN()  { [[ $STEP == $NOW ]] && echo -e "$RED" && LINES "$@" | sed 's/^/# WARNING: /' && echo -e "$OFF"; }
FATAL() { [[ $STEP == $NOW ]] && echo -e "$RED" && LINES "$@" | sed 's/^/# FATAL: /' && echo -e "$OFF"; exit 1; }
NOTE()  { [[ $STEP == $NOW ]] && echo -e "$YEL" && LINES "$@" | sed 's/^/# NOTE: /' && echo -e "$OFF"; }
LINE()  { [[ $STEP == $NOW ]] && echo && sed -e 's/^/ /;s/$/ /p;s/./-/g;' <<< "$@" | sed 's/^/#/'; }
ECHO()  { [[ $STEP == $NOW ]] && echo "$@"; }


# MARTINI Force field parameters
MARTINI=martini22
DRY=


# Solvents
# NOTE: with standard Martini water, exclude AA/CG interactions
SOLVENTS=( dry    polarizable      PW          pw         bmw       W     simple standard martini )
SOLNAMES=(  NO         PW          PW          PW         BMW       W       W       W        W   )
SOLTYPE=(   NO    polarizable polarizable polarizable polarizable plain   plain   plain    plain )
EPSR_CG=(   15        2.50        2.50        2.50        1.30     15      15      15       15   ) # CG-CG dielectric constant
EPSR_AA=(    0        1.45        1.45        1.45        1.00      0       0       0        0   ) # CG-AA dielectric constant
SOLVFF=(   dry          p           p           p          bmw                                   ) # MARTINI subversion and solvent tag


# Options:

# Downstream programs:

# - membrane and periodic boundary conditions (insane options):
INSANE=()

# - martinizing
MARTINIZE=()


# This program:

# - protein:
PDB=
TOP=
NDX=
NOHETATM=true
VSITE=false

# - multiscaling
ForceField=gromos45a3
MULTI=()
ALL=false
M=false
NCH=0
SOL=
HybridIons=false

# - nonbonded interactions
EPSR=
EPSRF=78
LJDP=6
LJRP=12
LJSW=-1
RC=1.4

# - other
DIR=.
NP=1
FORCE=false
MAXH=-1       # Maximum duration of run
NOEXEC=
JUNK=()
PROGOPTS=()
MDPOPTS=()
DAFT=
SCRATCH=      # Scratch directory
KEEP=false

# - simulation parameters
DELT=0.020       # Time step in picoseconds
TIME=0           # Nanoseconds
AT=0.5           # Write output every 0.5 ns
EMSTEPS=500      # Number of steps for EM
Temperature=310  # Degree Kelvin
Pressure=1       # Bar
Salinity=0.1536  # Concentration NaCl
SEED=$$


# Directories and executables
# This gives the directory where the script is located
SDIR=$( [[ $0 != ${0%/*} ]] && cd ${0%/*}; pwd )


OPTIONS=$(cat << __OPTIONS__
# File options: 
    -f  	 Input PDB file                                                *FILE:  None 
    -top         Input topology file                                           *FILE:  None
    -cg          Coarse grained force field                                    *STR:   $MARTINI
    -ndx         Index file                                                    *STR:   None
# Multiscaling: 
    -m           Chains for multiscaling. Other chains will be coarse grained. *STR:   None
    -M           Multiscale all chains, except those specified with -m         *BOOL:  False
    -ff          Atomistic force field for multiscaling                        *STR:   $ForceField
# Conditions:
    -sol         Solvent to use                                                *STR:   $SOL
    -epsr        Dielectric constant of the medium                             *FLOAT: $EPSR
    -T           Temperature                                                   *FLOAT: $TEMP (K) 
    -P           Pressure                                                      *FLOAT: $PRES (bar)
    -salt        Salt (NaCl) Concentration                                     *FLOAT: $Salinity (M)
# Simulation control: 
    -time        Simulation length in ns                                       *FLOAT: $TIME (ns)
    -at          Output resolution                                             *FLOAT: $AT (ns)
    -np          Number of processors to run on                                *INT:   $NP
    -vsite       Use virtual sites                                             *BOOL:  False
# Script control options: 
    -step        Step at which to start or resume the setup                    *STR:   $STEP
    -stop        Step at which to stop execution                               *STR:   $STOP
    -force       Force overwriting existing files                              *BOOL:  False
    -noexec      Do not actually simulate                                      *BOOL:  False
    -keep        Keep all rubbish generated; do not clean up                   *BOOL:  False
# Other options:
    -daft        Index file with DAFT energy groups. Invokes DAFT run.         *FILE:  None
    -dry         Use dry martini                                               *STR:   None

Arguments not listed are passed through to martinize.py
Check '$MARTINIZE -h' for a listing

__OPTIONS__
)


USAGE ()
{
    cat << __USAGE__

$PROGRAM version $VERSION:

$DESCRIPTION

$OPTIONS

(c)$YEAR $AUTHOR
$AFFILIATION

__USAGE__
}


BAD_OPTION ()
{
  echo
  echo "Unknown option "$1" found on command-line"
  echo "It may be a good idea to read the usage:"

  USAGE

  exit 1
}


# Functions for handling argument lists for downstream programs
# 1. Expanding options from '-opt=val1\,val2' to '-opt val1 val2', not changing long options (--long-opt=val)
function expandOptList()   { for i in $@; do [[ $i =~ --+ ]] && echo $i || (j=${i/=/ }; echo ${j//\\,/ }); done; }
# 2. Condensing options from '-opt1=val1,val2 -opt2=val3,val4' to '-opt1=val1\,val2,-opt2=val3\,val4' 
function condenseOptList() { echo $(sed 's/,/\,/g;s/  */,/g;' <<< $@); }
# 3. Reading condensed options from the command line
function readOptList() { sed "s/\\\,/##/g;s/,/ /g;s/##/,/g;s/--[^{]\+{\(.*\)}/\1/;" <<< $1; }


# Collect errors, warnings and notes to (re)present to user at the end
# Spaces are replaced by the unlikely combination QQQ to keep the 
# messages together.
errors_array=()
store_error_fun() { a="$@"; errors_array+=(${x// /QQQ}); FATAL "$@"; }
warnings_array=()
store_warning_fun() { a=$@; warnings_array+=(${x// /QQQ}); WARN "$@"; }
notes_array=()
store_note_fun() { a=$@; notes_array+=(${x// /QQQ}); NOTE "$@"; }


while [ -n "$1" ]; do

  # Check for program option
  depset=false
  NDEP=${#DEPENDENCIES[@]}
  for ((i=0; i<$NDEP; i++))
  do
      if [[ $1 == "-${DEPENDENCIES[$i]}" ]]
      then
          PROGEXEC[$i]=$2
          shift 2
          depset=true
      fi
  done
  # If we set a dependency, skip to the next cycle
  $depset && continue

  # Check for other options
  case $1 in
         -h)       USAGE                             ; exit 0 ;;
     # File options
         -f)                   PDB=$2                ; shift 2; continue ;;
        -cg)               MARTINI=$2                ; shift 2; continue ;;
       -top)                   TOP=$2                ; shift 2; continue ;;
       -ndx)                   NDX=$2                ; shift 2; continue ;;
    -hetatm)              NOHETATM=false             ; shift 1; continue ;;
       -sol)                   SOL=$2                ; shift 2; continue ;;
      -epsr)                  EPSR=$2                ; shift 2; continue ;;
     -epsrf)                 EPSRF=$2                ; shift 2; continue ;;
      -ljdp)                  LJDP=$2                ; shift 2; continue ;;
      -ljrp)                  LJRP=$2                ; shift 2; continue ;;
      -ljsw)                  LJSW=$2                ; shift 2; continue ;;
        -rc)                    RC=$2                ; shift 2; continue ;;
         -m)     MULTI[$((NCH++))]=$2; M=true        ; shift 2; continue ;;
         -M)                   ALL=true; M=true      ; shift 1; continue ;;
        -ff)            ForceField=$2                ; shift 2; continue ;;
         -T)           Temperature=$2                ; shift 2; continue ;;
         -P)              Pressure=$2                ; shift 2; continue ;;
      -salt)              Salinity=$2                ; shift 2; continue ;;
        -dt)                  DELT=$2                ; shift 2; continue ;;
      -time)                  TIME=$2                ; shift 2; continue ;;
        -at)                    AT=$2                ; shift 2; continue ;;
        -em)               EMSTEPS=$2                ; shift 2; continue ;;
       -dir)                   DIR=$2                ; shift 2; continue ;;
#     -gmxrc)                 GMXRC=$2                ; shift 2; continue ;;
#      -dssp)                  DSSP=$2                ; shift 2; continue ;;
      -step)                  STEP=$2                ; shift 2; continue ;;
      -stop)                  STOP=$2                ; shift 2; continue ;;
        -np)                    NP=$2                ; shift 2; continue ;;
      -maxh)                  MAXH=$2                ; shift 2; continue ;;
     -vsite)                 VSITE=true              ; shift 1; continue ;;
     -force)                 FORCE=true              ; shift 1; continue ;;
      -daft)                  DAFT=$2                ; shift 2; continue ;;
       -dry)                   DRY=$2                ; shift 2; continue ;;
      -keep)                  KEEP=true              ; shift 1; continue ;;
    -noexec)                NOEXEC=echo              ; shift 1; continue ;;
    --mdp-*)     MDPOPTS[${#MDPOPTS[@]}]=${1#--mdp-} ; shift 1; continue ;;

    # Options for downstream programs
    # If the options are given on the command line, they are expanded and each
    # option will be formatted as --program-opt=val
    --martinize-*) MARTINIZE+=(${1#--martinize})     ; shift  ; continue;;
    --insane-*)       INSANE+=(${1#--insane})        ; shift  ; continue;;

    # If the options are passed by another program, they will be formatted like
    #   --program{-opt1=val1,-opt2=val2\,val3}
    # In this case the option needs to be parsed explicitly:
    --martinize*)  MARTINIZE+=($(readOptList $1))    ; shift  ; continue;;
    --insane*)        INSANE+=($(readOptList $1))    ; shift  ; continue;;

    # Other program-specific options
    --*)   PROGOPTS[${#PROGOPTS[@]}]=$1              ; shift 1; continue ;;

    # All options should be covered above. Anything else raises an error here.
      *)         BAD_OPTION "$1";;

  esac
done


cat << __RUNINFO__

$PROGRAM version $VERSION:

(c)$YEAR $AUTHOR
$AFFILIATION

Now executing...

$CMD

__RUNINFO__

echo $CMD > cmd.log

NOW=$STEP

#--------------------------------------------------------------------
#---GROMACS AND RELATED STUFF
#--------------------------------------------------------------------


## 0. Finding programs

dependency_not_found_error()
{
    FATAL The required dependency $@ was not found.
}

NDEP=${#DEPENDENCIES[@]}
find_program_fun()
{
    for ((i=0; i<$NDEP; i++)); do
        if [[ ${DEPENDENCIES[$i]} == "$1" ]] 
        then
            progr=${PROGEXEC[$i]}
            envvar=${PROGEVAR[$i]}
        fi
    done

    # Check if the program is in the environment
    [[ -n $envvar ]] && [[ -f ${!envvar} ]] && echo ${!envvar} && return 0

    # Check if the program is in the directory of this script
    [[ -f $SDIR/$progr ]] && echo $SDIR/$progr && return 0

    # Check if the program is in the PATH
    which $progr 2>/dev/null && return 0 || return 1
}


##  1. GROMACS  ##

# Check and set the gromacs related stuff  
GMXRC=$(find_program_fun gmxrc)
echo Gromacs RC file: $GMXRC

# Source the gromacs RC file if one was found
# Otherwise the script will rely on the active gromacs commands 
[[ $GMXRC ]] && source $GMXRC 

# Find out which Gromacs version this is
GMXVERSION=$(mdrun -h 2>&1 | sed -n '/^.*VERSION \([^ ]*\).*$/{s//\1/p;q;}')
# From version 5.0.x on, the commands are gathered in one 'gmx' program
# The original commands are aliased, but there is no guarantee they will always remain
[[ -z $GMXVERSION ]] && GMXVERSION=$(gmx -h 2>&1 | sed -n '/^.*VERSION \([^ ]*\).*$/{s//\1/p;q;}')
ifs=$IFS; IFS="."; GMXVERSION=($GMXVERSION); IFS=$ifs

# Set the directory for binaries
[[ $GMXVERSION -gt 4 ]] && GMXBIN=$(which gmx) || GMXBIN=$(which mdrun)
# Extract the directory
GMXBIN=${GMXBIN%/*}
# Set the directory to SCRIPTDIR if GMXBIN is empty 
GMXBIN=${GMXBIN:-$SCRIPTDIR}

# Make binaries executable if they are not
# (This may be required for Grid processing)
[[ -f $GMXBIN/grompp && ! -x $GMXBIN/grompp ]] && chmod +x $GMXBIN/grompp
[[ -f $GMXBIN/mdrun  && ! -x $GMXBIN/mdrun  ]] && chmod +x $GMXBIN/mdrun
[[ -f $GMXBIN/gmx    && ! -x $GMXBIN/gmx    ]] && chmod +x $GMXBIN/gmx

# Set the command prefix
[[ $GMXVERSION -gt 4 ]] && GMX="$GMXBIN/gmx " || GMX=$GMXBIN/

# Set the GMXLIB variable to point to the force field data and such
# In some cases, 'gromacs' is part of $GMXDATA
export GMXLIB=${GMXDATA}/gromacs/top
[[ -d $GMXLIB ]] || export GMXLIB=${GMXDATA%/gromacs*}/gromacs/top
echo Gromacs data directory: $GMXLIB

# Now finally, test a command and see if it works
# otherwise raise a fatal error.
${GMX}grompp -h >/dev/null 2>&1 || executable_not_found_error "GROMACS (GMXRC)"


## 2. DSSP ##

# Search the DSSP binary, from environment, from path, or guess
# Only required if we have an input file
echo -n 'Checking DSSP binary (for martinizing proteins)... '
DSSP=$(find_program_fun dssp)
if [[ $? == 1 ]]
then
    warn="DSSP binary not found - Will martinize without secondary structure :S"
    store_warning_fun "$warn"
else
    echo "$DSSP"
    MARTINIZE+=(-dssp=$DSSP)
fi


## 3. PATH ##

# Add the script location to the PATH
export PATH="${PATH}:${SDIR}"


## 4. Echo mdp options specified on the command line

# These options are formatted like --mdp-param=value
# This will add 'param = value' to the MDP file for 
# all runs following energy minimization. 
# *NOTE*: Options specified on the command line take
# precedence over internal parameters and over those
# read from an mdp file, provided as value to option
# -mdp 

if [[ -n $MDPOPTS ]]
then
    echo 'Simulation parameters specified on command line (note how flexible!):'
    for ((i=0; i<${#MDPOPTS[@]}; i++)); do echo ${MDPOPTS[$i]}; done
    echo ===
fi


## 5. Locate insane if STEP lies before SOLVENT and STOP lies after.

SOLSTEP=$(get_step_fun SOLVENT)
if [[ $STEP -le $SOLSTEP && $STOP -ge $SOLSTEP ]]
then
    INSA=$(find_program_fun insane)
    if [[ $? != 0 ]]
    then
	FATAL "Dependency (insane) required for building solvent/membrane, but not found."
    fi
fi

## 5. Set the correct sed version for multi-platform use
# Also try to avoid GNU specific sed statements for the
# poor bastards that are stuck with one of those Mac things
SED=$(which gsed || which sed)


#--------------------------------------------------------------------
#---TIMING (GRID STUFF)
#--------------------------------------------------------------------

# If UNTIL is not set, then we run until all work is done, or until we crash
UNTIL=


# If we have a proxy, then we adhere to it, assuming to be running on the GRID
# The allowed runtime is set slightly lower to allow the script to finish
which voms-proxy-info && LEFT=$(( $(voms-proxy-info -timeleft &> /dev/null) - 300 )) || LEFT=


# If LEFT is set, MAXH may not be set negative
[[ $LEFT ]] && (( MAXH < 0 )) && MAXH=48


# Maximum time in seconds
if [[ $MAXH =~ ":" ]]
then
    # Format HH:MM:SS
    ifs=$IFS; IFS=":"; MAXS=($MAXH); IFS=$ifs
    MAXS=$((3600*MAXS[0] + 60*MAXS[1] + MAXS[2]))
else
    # Format x.y HH
    # BASH floating point arithmetics
    int=${MAXH%.*}
    frac=${MAXH#$int}
    frac=${frac/./}
    MAXS=$(( ${MAXH/./} * 3600 / (10**${#frac})  ))
fi


[[ -n $LEFT ]] && (( MAXS > LEFT )) && MAXS=$LEFT


if (( MAXS > 0 ))
then
    UNTIL=$(( $(date +%s) + MAXS ))
    echo "# $PROGRAM will run until $(date --date=@$UNTIL), or until run has finished"
else
    echo "# No maximum runtime set. Will run until finished or until crash."
fi

# This variable will be reset to the time needed for the last run run
# A following run will usually take longer.
LASTRUN=0


#--------------------------------------------------------------------
#---WARMING UP VARIABLE GYMNASTICS
#--------------------------------------------------------------------


## 1. Format for printing numbers to index files
fmt=" %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d"


## 2. WORKING DIRECTORY AND SOURCE DIRECTORY ##
SRCDIR=$(pwd)
[[ ! -d $DIR ]] && mkdir -p $DIR; pushd $DIR >/dev/null


## 3. START/STOP FLOW CONTROL ##
for ((i=0; i<${#STEPS[@]}; i++)); do [[ ${STEPS[$i]} == ${STEP}* ]] && STEP=$i && break; done
for ((i=0; i<${#STEPS[@]}; i++)); do [[ ${STEPS[$i]} == ${STOP}* ]] && STOP=$i && break; done
NOW=$STEP
echo Will run from step ${STEPS[$STEP]} until ${STEPS[$STOP]}

# Macro to do stuff only if the step is to be executed
DO() { [[ $STEP == $NOW ]] && echo "$@" && $@; }

# Sed wrapper - echo the command line before running it
SED() { echo sed "$@" 1>&2; sed "$@"; }

# Sequence generation; seq might not be available
SEQ() { for ((i=$1; i<=$2; i++)); do echo $i; done; }


## 4. SET THE SOLVENT  ##

# Select the solvent to use
# Default solvent is martini water
if [[ -z $SOL ]]
then
    if [[ $MARTINI == *p ]]
    then
	SOL=polarizable
    fi
fi
SID=; for ((i=0; i<${#SOLVENTS[@]}; i++)); do [[ ${SOLVENTS[$i]} == ${SOL}* ]] && SID=$i; done

# Override if option -dry is given... Serving Dry Martini
[[ -n $DRY ]] && SID=0

# Check whether we found a matching solvent model
[[ -z $SID ]] && echo Unknown solvent model \"$SOL\" specified. && exit

# Check whether the solvent type is polarizable
[[ ${SOLTYPE[$SID]} == polarizable ]] && POLARIZABLE=true || POLARIZABLE=false


## 5. FORCE FIELD ##

#     a. COARSE GRAINED: MARTINI/ELNEDYN
#        Set a pointer to the correct forcefield generating script
#        (of the form: martini_2.1.py or martini_2.1_P.py)
FFMARTINIPY=$SDIR/${MARTINI}${SOLVFF[$SID]}.py
if [[ ! -f $FFMARTINIPY ]]
then
    if [[ -f $SDIR/${MARTINI}.py ]]
    then
        # If martini22p was specified in stead of martini22 with PW,
	# then we end up here, setting the script to martini22p.py
	FFMARTINIPY=$SDIR/${MARTINI}.py
    else
	# Unclear dependency, but shipped with the scripts...
	FATAL Forcefield script $FFMARTINIPY does not exist, nor does $SDIR/${MARTINI}.py
    fi
fi


#     b. ATOMISTIC
#        Set pointers to ffnonbonded.itp and ffbonded.itp if we do multiscaling
$M && ffnb=$GMXLIB/$ForceField.ff/ffnonbonded.itp || ffnb=
$M && ffbn=$GMXLIB/$ForceField.ff/ffbonded.itp    || ffbn=

#     c. GENERATE martini.itp FOR COARSE GRAINED, MULTISCALE or DRY
if [[ -n $DRY ]]
then
    $FFMARTINIPY "$DRY" > martini.itp
else
    $FFMARTINIPY $ffnb $ffbn > martini.itp
fi

#     d. UPDATE martini.itp FOR DUMMIES
#        Replace the #include statement for ff_dum.itp for atomistic force fields by the contents of it
$M && sed -i -e "/#include \"ff_dum.itp\"/r$GMXLIB/$ForceField.ff/ff_dum.itp" -e "/#include \"ff_dum.itp\"/d" martini.itp


## 6. ELECTROSTATICS AND TABLES ##

EPSR_CG=${EPSR_CG[$SID]}
EPSR_AA=${EPSR:-${EPSR_AA[$SID]}}
$M && TABLES=-tables || TABLES=


#--------------------------------------------------------------------
#---INPUT CHECKING, SPLITTING, TRIMMING, GROOMING
#--------------------------------------------------------------------

if [[ -n $PDB ]]
then

    # If the input file is not found, check whether it was given without extension.
    # If that is not the case, then fetch the file from the PDB repository.
    if [[ ! -f $PDB ]]
    then
	if [[ ! -f $PDB.pdb ]]
	then
	    echo Input file $PDB not found... Trying server

            # Try fetching it from the PDB    
            pdb=$(tr [A-Z] [a-z] <<< ${PDB%.pdb})
            RCSB="http://files.rcsb.org/download/${pdb%%.*}.pdb.gz"
            echo "# Input file not found, but will try fetching it from the PDB ($RCSB)"
            # Use wget or curl; one of them should work
            wget $RCSB 2>/dev/null || curl -O "$RCSB" 
	    gunzip $pdb.pdb.gz
            [[ -n $SCRATCH ]] && cp $pdb.pdb $DIR
	    PDB=$pdb
	fi
	PDB=$PDB.pdb
    fi


    # If the input file is missing now, raise an errorr
    if [[ ! -f $PDB ]]
    then
	echo Input file $PDB not found and fetching from PDB server failed.
	exit 1
    fi


    # Check whether the input file is here or in another directory.
    # In the latter case, copy it here
    [[ $PDB == ${PDB##*/} || $PDB == ./${PDB##*/} ]] || cp $PDB .   


    pdb=${PDB##*/}          # Filename
    base=${pdb%.*}          # Basename
    ext=${pdb##*.}          # Extension
    dirn=${PDB%$pdb}        # Directory
    [[ $dirn ]] || dirn="." 
    dirn=`cd $dirn && pwd`  # Full path to input file directory


    if [ $dirn != `pwd` ]
    then
	NOTE "The run is performed here (`pwd`), while the input file is elsewhere ($dirn)."
    fi

    if [[ -z $TOP && $ext == "pdb" ]]
    then
	if $NOHETATM -a $(grep -q HETATM $PDB)
	then
	    NOTE Removing HETATM entries from PDB file
	    sed -i '' -e /^HETATM/d $PDB
	fi

        # Extract a list of chains from PDB file
	CHAINS=( $(grep '^\(ATOM\|HETATM\)' $PDB | cut -b 22 | uniq) )

        # Unpack lists of chains to multiscale separated by commas
        MULTI=( $(for i in ${MULTI[@]}; do echo ${i//,/ }; done) )

        # Residues defined in martinize.py
        AA=(ALA CYS ASP GLU PHE GLU HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR)
        # Sed query for residues (separated by \|):
	SED_AA=$(sed 's/ \+/\\\|/g' <<< ${AA[@]})
        # ATOM selection (martinizable residues)
	ATOM='/^\(ATOM  \|HETATM\)/{/.\{17\} *'$SED_AA' */p;}'
        # HETATM selection (non-martinizable residues)
	HETATM='/^\(ATOM  \|HETATM\)/{/.\{17\} *'$SED_AA' */!p;}'

        # Split the pdb file in stuff that can be processed with pdb2gmx 
        # and stuff that cannot be processed with it
        #  Extract the names of building blocks defined in the rtp file 
        # of the force field used. These blocks are defined as '[ ALA ]',
        # so we match a name in square brackets at the start of a line.
        # The name is appended to the hold space.
        RTPENTRIES='/^\[ *\(...\).*\].*/{s//\1/;H;}'
        # At the end of the list, the building block names are reformatted
        # to make a regular expression string matching each word. First,
        # the hold space is swapped with the pattern space, the first bit is 
        # removed and then all embedded newlines are replaced by '\|'
        FORMAT='${x;s/\n...\n//;s/\n/\\\|/g;p;}'
        # Finally sed is called processing all rtp files of the force field
        DEF=$(sed -n -e "$RTPENTRIES" -e "$FORMAT" $GMXLIB/$ForceField.ff/*.rtp)
 
        # Now we can split the input PDB file into a processable and a non-processable part
	echo $DEF
        #sed -e '/^\(TER\|MODEL\|ENDMDL\)/p' -e '/^\(ATOM  \|HETATM\)/{/.\{17\} *\('$DEF'\) */p;}' $dirn/$base.pdb  > $base-def.pdb
        #sed -e '/^\(TER\|MODEL\|ENDMDL\)/p' -e '/^\(ATOM  \|HETATM\)/{/.\{17\} *\('$DEF'\) */!p;}' $dirn/$base.pdb > $base-ndef.pdb  
    fi
else
    base=
fi


if [[ $PBC == retain && -n $PDB ]]
then
    PBC="-pbc keep"
else
    PBC="-pbc $PBC"
fi


echo Done checking



SED() { echo sed "$@" 1>&2; sed "$@"; }


## SED stuff

# Extract the moleculetype name
SED_MOLECULETYPE='/moleculetype/{
:a
n
/^;/ba
s/\s\+.*$//p
q
}
'

echo Done gymnastics


#--------------------------------------------------------------------
#---SIMULATION PARAMETERS--          
#--------------------------------------------------------------------

## OT N ## For every parameter not defined the default is used
## NOTE ## This is probably fine for equilibration, but check the defaults to be sure
## E OT ## The list as is was set up for gromacs 4.5

# This function lists the mdp options requested based on a preceding tag
mdp_options ()
{
    echo MDP tags: $@ 1>&2
    for tag in $@
    do
        # Find variables declared with specified tag
        for opt in `set | grep ^__mdp_${tag}__`
        do
            var=${opt%%=*}
            val=${opt#*=}
            # Replace the tag and redeclare in local space
            # If the variable was already declared it will
            # be overridden.
            local ${var/mdp_$tag/mdp}=$val
        done
    done
    # Find all variables starting with __mdp__ and echo them
    set | sed -n '/^__mdp__/{s/__mdp__//;s/,/ /g;p;}'
}

#--------------------------------------------------------------------
# Global parameters

__mdp_cg__nstlist=10
if [[ -n $DELT ]]
then
    __mdp_cg__dt=$DELT
elif $ALL -o [[ -n $MULTI ]]
then
    if $VSITE
    then
	__mdp_cg__dt=0.004
	__mdp_cg__nstlist=4
    else
	__mdp_cg__dt=0.002
    fi	
else
    __mdp_cg__dt=0.020
fi

# Output parameters
TIME=$(python -c "print int(1000*$TIME/$__mdp_cg__dt + 0.5 )") 
AT=$(python -c "print int(1000*$AT/$__mdp_cg__dt + 0.5)") 
__mdp_cg__nsteps=$TIME
__mdp_cg__nstxout=$AT
__mdp_cg__nstvout=0 
__mdp_cg__nstfout=0
__mdp_cg__nstlog=$AT
__mdp_cg__nstenergy=$AT
[[ $GMXVERSION -gt 4 ]] && __mdp_cg__nstxout_compressed=$AT || __mdp_cg__nstxtcout=$AT

# Nonbonded interactions
if [[ $GMXVERSION -gt 4 ]]
then
    # GMX 5 parameters from De Jong, Baoukina and Marrink
    __mdp_cg__cutoff_scheme=Verlet
    __mdp_cg__coulombtype=Cut-off
    __mdp_cg__coulomb_modifier=Potential-shift
    __mdp_cg__rcoulomb=1.1
    __mdp_cg__epsilon_r=$EPSR_CG
    __mdp_cg__vdw_type=Cut-off
    __mdp_cg__vdw_modifier=Potential-shift
    __mdp_cg__rvdw=1.1
    __mdp_cg__dispcorr=No
    __mdp_cg__nstlist=20
else
    __mdp_cg__coulombtype=Shift
    __mdp_cg__rcoulomb=1.2
    __mdp_cg__rcoulomb_switch=0.0
    __mdp_cg__epsilon_r=$EPSR_CG
    __mdp_cg__rlist=1.2
    __mdp_cg__vdw_type=Shift
    __mdp_cg__rvdw=1.2
    __mdp_cg__rvdw_switch=0.9
    __mdp_cg__dispcorr=No
fi

# Coupling - depending on presence of protein/membrane
# Pressure coupling semiisotropic for membranes,
# set to isotropic if no membrane is present
__mdp_cg__tcoupl=v-rescale
__mdp_cg__nsttcouple=$__mdp_cg__nstlist

__mdp_cg__pcoupl=no # Overridden at step NpT
__mdp_cg__nstpcouple=$__mdp_cg__nstlist

__mdp_cg__pcoupltype=Isotropic
__mdp_cg__compressibility=3e-4
__mdp_cg__tau_p=3.0 
__mdp_cg__ref_p=$Pressure

__mdp_mem__pcoupltype=Semiisotropic
__mdp_mem__compressibility=3e-4,3e-4
__mdp_mem__tau_p=3.0,3.0 # Overridden at step NpT with GMX5
__mdp_mem__ref_p=$Pressure,$Pressure

# Gromacs can handle empty groups, so this is fine
# These are set based on groups in the index file 
#__mdp_cg__tc_grps=Solute,Membrane,Solvent
#__mdp_cg__tau_t=1.0,1.0,1.0
#__mdp_cg__ref_t=$Temperature,$Temperature,$Temperature

# Other
__mdp_cg__constraints=none
__mdp_cg__energygrps=Solute,Membrane,Solvent
__mdp_cg__comm_mode=Linear
__mdp_cg__comm_grps=System
__mdp_cg__nstcomm=$__mdp_cg__nstlist


#--------------------------------------------------------------------
# Multiscale
__mdp_ms__rlist=$RC
__mdp_ms__coulombtype=user
__mdp_ms__rcoulomb=$RC
__mdp_ms__vdw_type=user 
__mdp_ms__rvdw=$RC
__mdp_ms__energygrps=AA,CG,VZ
__mdp_ms__epsilon_r=1
if $POLARIZABLE
then
    __mdp_ms__energygrp_table=AA,AA,AA,CG
    __mdp_ms__energygrp_excl=AA,VZ,VZ,VZ
else
    __mdp_ms__energygrp_table=AA,AA
    __mdp_ms__energygrp_excl=AA,VZ,AA,CG,VZ,VZ
fi

#--------------------------------------------------------------------
# Energy minimization

__mdp_em__integrator=steep
__mdp_em__nsteps=$EMSTEPS
__mdp_em__emstep=0.001
__mdp_em__lincs_order=6
__mdp_em__lincs_iter=8
__mdp_em__lincs_warnangle=90
__mdp_em__define=-DFLEXIBLE

#--------------------------------------------------------------------
# Equilibration runs: position restraints, NVT, NPT

# Position restraints are relieved at step 8
__mdp_equil__define=-DPOSRES
__mdp_equil__dt=0.002
__mdp_equil__nsteps=5000
__mdp_equil__nstlist=1
__mdp_equil__nstlog=8
__mdp_equil__nstenergy=8
[[ $GMXVERSION -gt 4 ]] && __mdp_cg__nstxout_compressed= || __mdp_cg__nstxtcout=0
__mdp_equil__tcoupl=v-rescale
__mdp_equil__lincs_order=6
__mdp_equil__lincs_iter=8
__mdp_equil__lincs_warnangle=90

__mdp_equil__tau_p=10

#--------------------------------------------------------------------
# User specified parameters


if [[ -n $MDP ]]
then
    # Check if the MDP file specified exists and exit if not
    [[ ! -f $MDP ]] && echo "MDP file $MDP specified, but not found" && exit
    
    # Gather options
    # 1. Delete empty lines
    # 2. Delete comment lines
    # 3. Remove whitespace on either side of the equality sign
    # 4. Remove trailing whitespace
    # 5. Replace spaces with commas
    USER_MDP=( $( sed '/^ *$/d;/^ *;/d;s/\s*=\s*/=/;s/\s*$//;s/ /,/g;' $MDP ) )
    
    for i in ${USER_MDP[@]}
    do 
        opt=${i%%=*}
        val=${i#$opt}
        eval "__mdp_usr__${opt//-/_}$val"
    done    
fi


if [[ -n $MDPOPTS ]]
then
    for i in ${MDPOPTS[@]}
    do
        # --mdp-energygrps=bla,bla,bla
        opt=${i%%=*}
        val=${i#$opt}
        eval "__mdp_usr__${opt//-/_}$val"
    done
fi


#--------------------------------------------------------------------


#--------------------------------------------------------------------
#---SUBROUTINES--
#--------------------------------------------------------------------


ERROR=0

archive ()
{
    if [[ -n $ARCHIVE ]]
    then
        tar cfz $ARCHIVE.tmp.tgz `ls --ignore=$ARCHIVE`
        mv $ARCHIVE.tmp.tgz $ARCHIVE
    fi
}
trap "archive" 2 9 15


exit_clean()
{
    [[ -f RUNNING ]] && rm -f RUNNING

    touch DONE

    if ! $KEEP
    then
        echo Deleting redundant files:
        printf "%-25s %-25s %-25s %-25s %-25s\n" ${JUNK[@]} \#*\#
        for i in ${JUNK[@]} \#*\#; do [[ -f $i ]] && rm $i; done
    fi

    exit 0
}


exit_error()
{

  echo "**"
  echo "** Something went wrong in running script martinize.sh from"
  echo "** `pwd`/:"
  echo "**"
  echo "** exit code: $1"
  echo "** $MSG"
  echo "**"

  [[ -f RUNNING ]] && rm -f RUNNING
  touch ERROR

  archive

  exit 1
}


# Trashing files
trash()
{
    for item in $@; do JUNK[${#JUNK[@]}]=$item; done
}


# Check for the existence of all arguments as files
all_exist()
{
  local have=true 
  for f in ${OUTPUT[@]} 
  do 
      [[ -e $f ]] || have=false
  done
  echo $have
}


# Parsing program specific options
program_options()
{    
    local OPTS=
    for opt in ${PROGOPTS[@]}
    do
        if [[ $opt =~ --$1 ]]
        then
            OPTS="$OPTS $(sed 's/--[^-]*//;s/=/ /' <<< $opt)"
        fi
    done
    echo $OPTS
}

# Amino acids
amino_acids=(ALA CYS ASP GLU PHE GLY HIS HIH ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR)
function sed_protein()
{
    for a in ${amino_acids[@]}
    do 
	echo '-e /'$a/$1
    done
}

# Nucleic acids
nucleic_acids=(DG DA DC DT)
function sed_nucleic()
{
    for a in ${nucleic_acids[@]}
    do 
	echo '-e /'$a/$1
    done
}

# Solvent
solvent_names=(W\ \+W PW BMW SOL)
function sed_solvent()
{
    for a in ${solvent_names[@]}
    do 
	echo '-e /'$a/$1
    done
}

function INDEX()
{
    [[ -n $2 ]] && fn=$2 || fn=basic.ndx
 
    exec 6>&1 && exec >$fn

    fmt="%5d %5d %5d %5d %5d %5d %5d %5d %5d %5d"

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


MDRUNNER ()
{
    local NP=1
    local fnOUT=
    local fnNDX=
    local FRC=false
    local SPLIT=
    local TABLES=
    while test -n "$1"; do
        case $1 in
            -f)      local  fnMDP=$2        ; shift 2; continue;;
            -c)      local   fnIN=$2        ; shift 2; continue;;
            -n)      local  fnNDX=$2        ; shift 2; continue;;
            -o)      local  fnOUT=$2        ; shift 2; continue;;
            -p)      local  fnTOP=$2        ; shift 2; continue;;
            -l)      local  fnLOG=$2        ; shift 2; continue;;
            -force)  local    FRC=$2        ; shift 2; continue;;
            -np)     local     NP=$2        ; shift 2; continue;;
            -split)  local  SPLIT=-noappend ; shift  ; continue;;
	    -tables) local TABLES="-table table.xvg -tablep tablep.xvg"; shift; continue;;
            *)  echo "PANIC!: Internal Argument Error ($1) in routine MDRUNNER"; exit;;
        esac
    done
    local TPR=

    # Check input
    [[ -f $fnIN  ]] || exit_error "Input structure not found in routine MDRUNNER ($fnIN)"
    [[ -f $fnTOP ]] || exit_error "Input topology not found in routine MDRUNNER ($fnTOP)"
    [[ -n $fnNDX && -f $fnNDX ]] || INDEX $fnIN $fnNDX

    # Infer basename from output file name
    baseOUT=${fnOUT%.*}

    #if [[ -n $SCRATCH ]]
    #then
    #    # Make sure to copy the TPR file if we already have one
    #    [[ -f $DIR/$baseOUT.tpr ]] && cp $DIR/$baseOUT.tpr .    
    #
    #    # Make sure we deal well with checkpointing
    #    [[ -f $DIR/$baseOUT.cpt ]] && cp $DIR/$baseOUT.cpt .
    #fi

    # Check if there are parts of runs (always in the RUN directory)
    if [[ -n $SPLIT ]]
    then
        # A neat way to check for the last of a set of numbered files
        # The work is done in the while test, incrementing the number 
        # before checking if the file exists. The statement is empty 
        # (using the true command ':'). The while loop ends with z being
        # equal to the first number not used.
        local z=0; while [[ -f $DIR/$baseOUT.part$(printf "%04d" $((++z))).log ]]; do :; done
        # The last existing file has number one less than z
        last=$baseOUT.part$(printf "%04d" $((z-1))).log
        # The logfile to write -- THIS IS PROBABLY REDUNDANT! (it should already be written like that)
        log=$baseOUT.part$(printf "%04d" $z).log
    else        
        log=$baseOUT.log 
        last=$DIR/$log
    fi

    # Check whether the log file actually exists
    step=0
    if [[ -e $last ]]
    then
        # Get the last Step/Time/Lambda listed in the last log file
        # Nifty: At each line containing Step/Time/Lambda, 
        # read a next line (n) and put it in the hold space (h)
        # At the end of the file, switch the hold space and the 
        # pattern space (x) and print (p)
        step=($($SED  -n -e '/^ *Step *Time *Lambda/{n;h;}' -e '${x;p;}' $last))
    fi
        
    local nsteps=$(sed -n '/nsteps/s/^.*=\([^;]*\).*$/\1/p' $fnMDP)


    # Check whether we need to do anything
    if $FRC
    then
        removed=()
        for z in ${baseOUT}.*; do [[ -f $z ]] && rm $z && removed[${#removed[@]}]=$z; done
        echo "# Forced execution."
        [[ -n $removed ]] && echo "# Removed files:" && echo ${removed[@]}
    elif [[ -f $fnOUT || -f $DIR/$fnOUT ]]
    then
        echo "# Output found ($fnOUT). Skipping step."
        return 0
    elif [[ $step -gt $nsteps ]]
    then
        echo "# A log file exists which reports having run $step of $STEPS steps ($last)"
        return 0
    fi


    echo "# $(date): STARTING MDRUNNER"    

    # Set the options for the parameter and index files
    fnMDP="-f $fnMDP -po ${fnMDP%.mdp}-out.mdp"
    [[ -n $fnNDX ]] && fnNDX="-n $fnNDX"


    # Skip generation of run input file if it exists
    if [[ ! -e $baseOUT.tpr ]]
    then
        # Build command
        # The warnings are pretty much controlled. We usually get a warning because of using a plain cut-off for EM.
        # With custom ligand topologies we may get warnings about overriding parameters. This should be fine? :S
        GROMPP="${GMX}grompp $fnMDP -c $fnIN -p $fnTOP $fnNDX -o $baseOUT.tpr $(program_options grompp) -maxwarn -1"

        echo "$GROMPP" | tee $LOG
        echo ": << __GROMPP__" >>$LOG
        $GROMPP >>$LOG 2>&1 || exit_error "Execution of grompp failed in routine MDRUNNER. More information in $LOG."
        echo "__GROMPP__" >>$LOG
    fi


    # If the output file extension is 'tpr' we should be done here
    [[ "${fnOUT#$baseOUT}" == ".tpr" ]] && return 0


    # If we extend a partial run, mention it
    [[ ! -e $fnOUT && -e $baseOUT.cpt ]] && echo $(date): FOUND PARTIAL ${STEPS[$STEP]} RESULTS... CONTINUING


    # Check the time
    # MAXH needs to be updated after every cycle
    if [[ -n $UNTIL ]]
    then
        MAXS=$(( UNTIL - $(date +%s) ))
        if (( $MAXS < 0 ))
        then
            exit_error "Somehow we violated the allowed run time... Exiting NOW!"
        fi
        MAXH=$(bc <<< "scale=3;$MAXS/3600") 
        
        # May want to check if MAXS makes sense for the run...
        if ((MAXS < LASTRUN))
        then
            echo "# The remaining time is $MAXS seconds, which does not seem sufficient for this part."
            echo "# The previous part took $LASTRUN seconds."
            exit_error "INSUFFICIENT TIME LEFT. RUN INCOMPLETE. RESUBMIT."
        fi
    fi

    
    # Run the run
    MDRUN="${GMX}mdrun -nice 0 -deffnm $baseOUT -c $fnOUT -cpi $baseOUT.cpt -nt $NP $SPLIT $(program_options mdrun) -maxh $MAXH"
    echo "$MDRUN" | tee -a $fnLOG
    echo ": << __MDRUN__" >>$LOG
    $MDRUN >>$fnLOG 2>&1 || exit_error "Execution of mdrun failed in routine MDRUNNER. More information in $LOG."
    echo "__MDRUN__" >>$LOG


    # If we split then we have to do some more work to see if we have finished
    if [[ -n $SPLIT ]]
    then
        step=($SED -n -e '/Step *Time *Lambda/{n;h;}' -e '${x;p;}' $fnLOG)
        [[ $step == $FIN ]] && cp ${fnLOG%.log}.gro $fnOUT
    fi

    
    # If $fnOUT exists then we finished this part
    if [[ -e $fnOUT ]]
    then
        echo "# $(date): FINISHED MDRUNNER (STEP ${STEPS[$STEP]})" | tee -a $fnLOG
        return 0
    else
        echo "# $(date): MDRUN EXITED (STEP ${STEPS[$STEP]}), BUT RUN NOT COMPLETE" | tee -a $fnLOG
        return 1
    fi
}


TABLE ()
{
    ARGS=$(sed 's/ /,/g' <<< $@)

    # Create a table for the Coulomb/RF interaction
    python - << __PYTHON__ 

epsR, epsRF, RC, LD, LR, LC, LS = $ARGS

Krf = (epsRF-epsR)/(2*epsRF+epsR)/(RC**3)
Crf = 3*epsRF/(2*epsRF+epsR)/RC

# Shifted power
def vf_fun(a,rc,rs):
    if rs < 0: 
        return lambda r: r**-a, lambda r: a*r**-(a+1)
    A, B = -((a+4)*rc-(a+1)*rs)/((rc-rs)**2)/(rc**(a+2)), ((a+3)*rc-(a+1)*rs)/((rc-rs)**3)/(rc**(a+2))
    C    = (rc**-a)-(a*A*(rc-rs)**3/3)-(a*B*(rc-rs)**4/4)
    V    = lambda r: r**-a - C - (r>rs and a*(A*(r-rs)**3/3+B*(r-rs)**4/4))
    F    = lambda r: a*r**-(a+1) + (r>rs and a*(A*(r-rs)**2+B*(r-rs)**3))
    return V, F

VD,FD = vf_fun(LD,LC,LS)
VR,FR = vf_fun(LR,LC,LS)

r  = [ i/1e3 or 0 for i in range(2,3001,2)]
print """#
# Coulomb cut-off/reaction-field: epsRF = %f, epsR = %f, RC = %f
# Lennard-Jones dispersion: power=%f, cutoff=%f, %sshifted %s
# Lennard-Jones repulsion:  power=%f, cutoff=%f, %sshifted %s
""" % (epsRF, epsR, RC, 
       LD, LC, LS<0 and "not " or "", LS>=0 and "(rshift=%f)"% LS or "",
       LR, LC, LS<0 and "not " or "", LS>=0 and "(rshift=%f)"% LS or ""),

f0 = [ (1/i+Krf*i*i-Crf)/epsR for i in r ]
f1 = [ (1/(i*i)-2*Krf*i)/epsR for i in r ]

g0 = [ i<LC and -VD(i) for i in r ]
g1 = [ i<LC and -FD(i) for i in r ]

h0 = [ i<LC and  VR(i) for i in r ]
h1 = [ i<LC and  FR(i) for i in r ]

table = [(0,0,0,0,0,0,0)]+zip(r,f0,f1,g0,g1,h0,h1)

print "\n".join([(7*"%.10e ")%i for i in table])
__PYTHON__
}



# Always ECHO the first line
NOW=$STEP

#--------------------------------------------------------------------
SHOUT "---= THIS IS WHERE WE START =---"
#--------------------------------------------------------------------

NOW=0

                
if [[ -n $DAFT ]] && ($ALL || [[ -n $MULTI ]])
then
    echo "Currently, DAFTly splitting molecules in energy groups and running multiscaled is not possible."
    echo "Will run using default multiscale splitting."
    # Simply setting DAFT to false should disable any attempt to do splitting
    DAFT=
fi	



# Set the basename to 'membrane' if no input structure is given
[[ -n $PDB ]] || base=membrane


#---------------------------------------------------------------------
SHOUT "---STEP 1A: GENERATE ATOMISTIC STRUCTURE AND TOPOLOGY"
#---------------------------------------------------------------------


# Skip the atomistic topology if we do not multiscale
if [[ $STEP == $NOW ]]
then
    echo "Skipping step... (not multiscaling)"
    $M  || : $((STEP++))
fi


# Skip this step if we run DAFT
if [[ $STEP == $NOW ]]
then
    echo DAFT run. Skipping step.
    [[ -n $DAFT  ]] && : $((STEP++))
fi


if [[ $STEP == $NOW ]]
then
    # Output for this section:
    OUT=$base-aa.pdb
    TOP=$base-aa.top
    LOG=01-TOPOLOGY-AA.log

    OUTPUT=($OUT $TOP)

    # Delete existing output if we force this step
    $FORCE && rm ${OUTPUT[@]}

    if $(all_exist ${OUTPUT[@]})
    then
	echo Found $TOP and $OUT. Skipping structure/topology building with pdb2gmx...
	: $((STEP++))
    fi
fi


MDPMS=
ForceFieldAA=


if [[ $STEP == $NOW ]]
then

    # I. Set the mdp tag
    MDPMS=ms


    # II. Atomistic force field for multiscaling

    #    1. List force fields available
    AAFF=($(ls -d $GMXLIB/*.ff | sed 's#.*/\(.*\)\.ff#\1#'))

    #    2. Try complete matching
    for i in ${AAFF[@]}; do [[ "$i" == "$ForceField" ]] && ForceFieldAA=$i; done

    #    3. Try partial matching if no complete match was found
    match=0
    [[ -n $ForceFieldAA ]] || for i in ${AAFF[@]}; do [[ $i =~ $ForceField ]] && ForceFieldAA=$i && : $(( match++ )); done
    if [[ $match -eq 0 && -z $ForceFieldAA ]]
    then
	echo No matching atomistic forcefield found for $ForceField... Bailing out.
	exit 1
    elif [[ $match -gt 1 ]]
    then
	echo Ambiguous selection for atomistic force field... Bailing out.
	exit 1
    fi
  
    
    # III. Interaction tables   
    #      epsilon_r epsilon_rf cutoff LJ_dispersion LJ_repulsion LJ_cutoff LJ_switch
    TABLE  $EPSR_CG   $EPSRF     $RC      $LJDP         $LJRP         1.2       0.9   > table.xvg
    TABLE  1          $EPSRF     $RC      6             12            $RC      -1     > table_AA_AA.xvg 
    TABLE  1          $EPSRF     $RC      6             12            $RC      -1     > tablep.xvg 
    if $POLARIZABLE
    then
	TABLE $EPSR_AA $EPSRF    $RC      6             12            $RC      -1     > table_AA_CG.xvg 
    fi

    
    # IV. pdb2gmx

    #     1. Basic stuff
    PDB2GMX="${GMX}pdb2gmx -v -f $dirn/$base.pdb -o $OUT -p $TOP"
    PDB2GMX="$PDB2GMX -i $base-posre.itp -posrefc 200 -ignh -ff $ForceFieldAA -water none"

    #     2. Virtual sites
    $VSITE && PDB2GMX="$PDB2GMX -vsite hydrogens" 

    #     3. Add program options specified on command line
    PDB2GMX="$PDB2GMX $(program_options pdb2gmx)"

    #     4. Specification of protonation states
    if [[ -e $dirn/$base.tit || -e $dirn/titratables.dat ]]
    then
	[[ -e $dirn/$base.tit ]] && TITR=$base.tit || TITR=titratables.dat
	echo SETTING PROTONATION STATES FROM $dirn/$TITR
        # Acidic residues ASP and GLU: deprotonated=0 protonated=1
	ACID='/\(ASP\|GLU\)/{s/.*[Hh0]\s*$/1/;s/.*\-\s*/0/}'
        # Basic residue LYS: deprotonated=0 protonated=1
	LYS='/LYS/{s/.*0\s*$/0/;s/.*[Hh+]\s*$/1/}'
        # Histidine: 
	HIS='/HIS/{s/.*[DdAa]\s*$/0/;s/.*[EeBb]\s*$/1/;s/.*[Hh\+]\s*$/2/}'
        # N-terminal
	NTER='/NTER/{s/.*\+\s*$/0/;s/.*0\s*$/1/;s/.*N\s*/2/}'
        # C-terminal
	CTER='/CTER/{s/.*\-\s*$/0/;s/.*0\s*$/1/;s/.*N\s*/2/}'
    
	sed -e "$NTER" -e "$CTER" -e "$ACID" -e "$LYS" -e "$HIS"  $dirn/$TITR > pdb2gmx.query
	trash pdb2gmx.query
    fi

    #     5. Check for gmxquery and, if found, set interactive mode 
    if [ -e "$dirn/pdb2gmx.query" ]; then
	GMXQUERY=$(cat $dirn/pdb2gmx.query)
	PDB2GMX="$PDB2GMX -ter -lys -asp -glu -his"
    else
	GMXQUERY=
    fi

    #     6. Process
    if [[ -n $NOEXEC ]] || $(all_exist ${OUTPUT[@]})
    then
        # Skipping verbosely; showing what would have been run
	echo -e "Skipping: \n$PDB2GMX"   
	echo "Also skipping atomistic topology acrobatics"
    else

	LOG=01-TOPOLOGY-AA.log
	
	echo "echo $GMXQUERY | $PDB2GMX" | tee -a $LOG
	echo $GMXQUERY | $PDB2GMX >>$LOG 2>&1 || exit_error 1


        #-- Lots of bookkeeping:
        #-- Simplifying topology if identical chains are found

        # a. First get the list of moleculetype itp files
	ITP=($(ls ${base}-aa_*.itp 2>/dev/null) ${base}-aa.top)

        # b. Extract the moleculetype names from the itp files
        #    awk expression:
        #    - at the line matching 'moleculetype' 
        #      read in the next line
	#      continue reading next lines until one is not beginning with ;
        #      print the first field
	awkexpr='/moleculetype/{getline; while ($0 ~ /^ *;/) getline; print $1}'
	MTP=($(for i in ${ITP[@]}; do awk "$awkexpr" $i; done))

        # c. Compare each pair of itp files
	for ((i=1; i<${#ITP[@]}; i++))
	do
            for ((j=0; j<$i; j++))
            do
                # Check how many lines are different, excluding includes and comments
                # For identical moleculetypes only the moleculetype name should differ
                # That will give four lines of output from diff
		if [[ $(diff -I "^\(;\|#include\)" ${ITP[$i]} ${ITP[$j]} | wc -l) == 4 ]]
		then
                    echo Removing duplicate moleculetype definition in ${ITP[$i]}
                    # 1. remove the #include statement for that itp file
                    # 2. rename the moleculetype under [ system ]
                    SED -i -e "/${ITP[$i]}/d" -e "/[\s*system\s*]/,\$s/${MTP[$i]}/${MTP[$j]}/" $base-aa.top
                    # List the file for removal
                    trash ${ITP[$i]} $(sed 's/_/-posre_/' <<< ${ITP[$i]})
                    break
		fi
            done
	done        


        # d. Extract the chain identifiers from the pdb2gmx log, listed as:
        # 
        #   chain  #res #atoms
        #  1 'B'     7     53  
        #  2 'A'     7     53  
        #  3 'C'     8     67  
        #
        # In case there is one chain without identifier, set the identifier
	CHAINS=(`sed -n '/chain  #res #atoms/,/^\s*$/s/^\s\+[0-9]\+..\(.\).*$/\1/p' 01-TOPOLOGY-AA.log`)
	[[ -n $CHAINS ]] || CHAINS=(A)
	MS=()
	N=0;
	for i in ${CHAINS[@]} 
	do 
	    MS[$N]=$ALL 
	    for j in ${MULTI[@]} 
	    do 
		[[ $i == $j ]] && MS[$N]=true
	    done
	    : $((N++)) 
	done 

        # e. If there is only one chain, gromacs writes a single topology file,
        #    which will make the protocol choke. That is fixed here.
        #    The moleculetype definition is taken from the topology file and
        #    put in a separate include topology.
        #    In addition, the numbers of residues and atoms are not listed 
        #    in the output from pdb2gmx for single chains as they are for 
        #    multiple chains. So for single chains we just take the atom
        #    count from the coordinate file. 
	if [[ ${#CHAINS[@]} -eq 1 ]]
	then
	    ITP=${base}-aa_molecule_chain_${CHAINS[0]}.itp
	    if [[ ! -f $ITP ]] # || [[ $ITP -ot $base-aa.top]]
	    then
                # Ahh, sed magic. Removing lines from one file, 
	        # while writing them to another.
		#@@@
		SED -i.bck \
		    -e "/^ *\[ *moleculetype *\] */{h;s/.*/#include \"$ITP\"/p;x;}" \
		    -e '/moleculetype/,/^#endif/w'$ITP \
		    -e '/moleculetype/,/^#endif/d' \
		    $base-aa.top
	    fi
	fi

        # f. Get the current list of molecules from the top file 
        #    and find the corresponding itp files
        #    Some of these may have been changed to remove duplicate listings
        #    This list should match the list of chains
	#    awk expression
	#    - if M is nonzero and the line does not start with ';', print
	#    - at the line with the tag [ molecules ] set M to 1
        #   Following each molecule name is the number of instances
	awkexpr='(M && ! /^ *;/); /\[ *molecules *\]/{M=1}'
	MOLECULES=($(awk "$awkexpr" $base-aa.top))
	echo MOLECULES: ${MOLECULES[@]}	


        # g. Get the list of itp files #included in the top file
        #    One sed expression:
        #    - at lines starting with #include and not containing a slash in the file name
        #      substitute the line with the name of the included file and print
	sedexpr='/^#include .*"\([^/]*\)"/s//\1/p'
	ITPFILES=($(sed -n -e "$sedexpr" $base-aa.top))


        # h. Get the list of moleculetypes from the itp files
        #    One awk expression:
        #    - at the line matching 'moleculetype' 
        #      read in the next line
	#      continue reading next lines until one is not beginning with ;
        #      print the first field
	awkexpr='/moleculetype/{getline; while ($0 ~ /^ *;/) getline; print $1}'
	MOLTYPES=($(for i in  ${ITPFILES[@]}; do awk "$awkexpr" $i; done))


        # i. Match the molecules with the itp files
        #    For each molecule in the list, check which moleculetype it is and 
        #    and fetch the corresponding itp file
        #    Mind that the molecule list also has a count per molecule
        #    following each name
	MOLITP=()
	echo "Molecules and corresponding topology files:"
	for ((i=0; i<${#MOLECULES[@]}; i+=2)) 
	do 
	    for ((j=0; j<${#MOLTYPES[@]}; j++))
	    do 
		[[ ${MOLECULES[$i]} == ${MOLTYPES[$j]} ]] && MOLITP[$((i/2))]=${ITPFILES[$j]}
	    done
	    printf "    %20s : %20s\n" ${MOLECULES[$i]}: ${MOLITP[$((i/2))]}
	done	
    fi 

    #     7. Set multiscale arguments for Martinize
    $ALL && M_MULTI="-multi all" || M_MULTI=$(sed 's/\(^\| \)/ -multi /g' <<< ${MULTI[@]})	

    PDB=$OUT
fi


# END OF ATOMISTIC TOPOLOGY

[[ $STOP ==   $NOW     ]] && exit_clean
[[ $STEP == $((NOW++)) ]] && : $((STEP++))


#---------------------------------------------------------------------
SHOUT "---STEP 1B: GENERATE COARSE GRAINED STRUCTURE AND TOPOLOGY"
#---------------------------------------------------------------------

# Skip this step if we run DAFT. In that case GRO equals PDB.
GRO=$PDB
if [[ $STEP == $NOW ]]
then
    [[ -n $DAFT  ]] && : $((STEP++))
fi

if [[ $STEP == $NOW ]]
then
    # Output for this section:
    GRO=$base-cg.gro
    TOP=$base-cg.top
    NDX=$base-cg.ndx
    LOG=01-TOPOLOGY-CG.log

    OUTPUT=($GRO $TOP $NDX)

    # Delete existing output if we force this step
    $FORCE && rm ${OUTPUT[@]}

    if $(all_exist ${OUTPUT[@]})
    then
	echo Found $TOP and $GRO. Skipping step with martinize/insane...
	: $((STEP++))
    fi
fi

if [[ $STEP == $NOW ]]
then
    NPROT=0

    # If a protein is given, build the martinize command line and echo it
    # regardless of whether we actually run.
    if [[ -n $pdb ]]
    then
        # Convert structure to coarse grained and build topology using martinize.py
	MART=$(find_program_fun martinize)
	if [[ $? != 0 ]]
	then
	    FATAL "Coarse graining PDB file ($pdb), but martinize was not found."
	fi
	MARTINIZE="$MART $(expandOptList ${MARTINIZE[@]})"
	MARTINIZE="$MARTINIZE -f $pdb -o $TOP -x $base-mart.pdb -n $base-mart.ndx -ff $MARTINI $M_MULTI"
	echo $MARTINIZE
    fi

    # Only if we have a pdb file and we actually run (NOEXEC is not set) 
    # then this block is executed
    if [[ -n $pdb && -z $NOEXEC ]]
    then
        # Executing the command. 
        # We need to know the moleculetype names (and itp files) after martinizing
        # 1. The command is executed
        # 2. stdout and stderr are swapped
        # 3. stdout (originally stderr) is also written to stderr using 'tee'
        # 4. stdout is parsed by sed, extracting moleculetype names for each chain
        # 5. The result is stored as array MARMOLS
        #          |-1-| |------2-----|  |-------3------|   |--------------------------4---------------------|
	MARMOLS=($($MARTINIZE 3>&1 1>&2 2>&3 | tee /dev/stderr | sed -n '/^INFO *[0-9]\+-> \+\([^ ]\+\).*$/s//\1/p'))

        # If we use polarizable or BMW water, we need to change AC1/AC2 to C1/C2
	if $POLARIZABLE
	then
	    sed -i -e 's/AC[12]/ C[12]/' $base-mart.pdb
	fi

        # Check for chains to multiscale
	if $ALL -o [[ -n $MULTI ]]
	then
            # Have to have constraints.itp for multiscale topologies to work
	    if [[ ! -e $GMXLIB/$ForceField.ff/constraints.itp && ! -e ./constraints.itp ]]
	    then
                # Generate from ffbonded.itp
		awk '/#define *gb_/{sub("gb_","gc_"); print $1, $2, $3}' $GMXLIB/$ForceField.ff/ffbonded.itp > constraints.itp 

                # Add specific constrainttypes:
		S_S='S      S       1    0.2040'
		NR_FE='NR     FE      1    0.1980'
		S_CR1='S      CR1     1    0.1830'
		SED -i -e "/angletypes/s/^/[ constrainttypes ]\n${S_S}\n${NR_FE}\n${S_CR1}\n\n/" martini.itp
	    fi

            # #include the constraintsfile in the master topology if it is not already
	    grep -q '#include *"constraints.itp"' $TOP || SED -i '/#include *"martini.itp"/s/$/\n\n#include "constraints.itp"\n\n/' $TOP

	    for ((i=0,j=0; i<${#MOLECULES[@]}; i+=2,j++))
	    do
		if ${MS[$j]}
		then
   		    # Adapt the atomistic topology for multiscaling
		    awk '/^\[ *bonds *\]/{print "#ifdef FLEXIBLE"; bonds=1} 
                         sub("^ *"P,M) || 1
                         bonds {
                             B[N++]=$0 
                             if ($0 ~ /^ *$/) {
                                 print "#else\n[ constraints ]"
                                 bonds=0
                                 for (i=1; i<N-1; i++) {
                                     if (split(B[i],Q)==3) 
                                         printf "%5d %5d\n", Q[1], Q[2]
                                     else {
                                         sub(/2 *gb_/,"1 gc_",B[i])
                                         print B[i]
                                     }
                                 } 
                                 print "#endif\n"
                             }
                         }' P=${MOLECULES[$i]} M=${MARMOLS[$j]} ${MOLITP[$j]} > ${MARMOLS[$j]}_MS.itp

      		    # Add the virtual sites generated by martinize.py
		    cat ${MARMOLS[$j]}.itp >> ${MARMOLS[$j]}_MS.itp

    		    # Update the moleculetype #include file in the master topology
		    SED -i -e '/^\(#include \+"\)'${MARMOLS[$j]}.itp'/s//\1'${MARMOLS[$j]}_MS.itp'/' $TOP
		fi
	    done
	elif [[ -n $DAFT ]]
	then
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
            awkexpr='(M && ! /^ *;/); /\[ *molecules *\]/{M=1}'
            MOLECULES=($(awk "$awkexpr" $base-cg.top))
            echo MOLECULES: ${MOLECULES[@]}

            sedexpr='/^#include .*"\([^/]*\)"/s//\1/p'
            ITPFILES=($(sed -n -e "$sedexpr" $base-cg.top))

            # We only take the FIRST moleculetype from the itp file
            awkexpr='/moleculetype/{getline; while ($0 ~ /^ *;/) getline; print $1; exit}'
            MOLTYPES=($(for i in  ${ITPFILES[@]}; do awk "$awkexpr" $i; done))

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
                            echo -e "[ ${MOLECULES[$i]}_$K ]\n$(printf "$fmt\n" `SEQ $COUNT $((COUNT+N-1))` | sed 's/ 0//g')" >> EnergyGroups.ndx
                            : $((COUNT+=N)) $((K++))
			done
                    fi
		done
		printf "    %20s : %20s\n" ${MOLECULES[$i]}: ${MOLITP[$((i/2))]}
            done  	                   
            __mdp_cg__energygrps=$(sed 's/ /,/g' <<< "${EnergyGroups[@]} Membrane Solvent")
            echo @@ $__mdp_cg__energygrps
	fi

	touch residuetypes.dat elements.dat
	trash residuetypes.dat elements.dat
	${GMX}editconf -f $base-mart.pdb -o $base-mart.gro -d 100 -noc >/dev/null 2>&1

	NPROT=$(awk '{getline;  print; exit}' $base-mart.gro)

        # Have to energy minimize output structure for stability
	#__mdp_mart__pbc=no
	MDP=em-mart.mdp
	OPT=(em $MDPMS mart)
	OUT=$base-mart-EM.gro
	
	mdp_options ${OPT[@]} > $MDP
	MDRUNNER -f $MDP -c $base-mart.gro -p $TOP -o $OUT -n $base-mart.ndx -np 1 -l $LOG -force $FORCE $TABLES
	echo 0 | ${GMX}trjconv -s $base-mart.gro -f $OUT -pbc nojump -o $base-mart-nj.gro
	mv $base-mart-nj.gro $OUT

        # trash $base-mart.gro
	GRO=$OUT
    fi
fi


# If $GRO is not set, we set it equal to $pdb
[[ -z $GRO ]] && GRO=$pdb


# END OF COARSE GRAINING

[[ $STOP ==   $NOW     ]] && exit_clean
[[ $STEP == $((NOW++)) ]] && : $((STEP++))


#---------------------------------------------------------------------
SHOUT "---STEP 1C: SOLVATE"
#---------------------------------------------------------------------

# At this point, we need to know how many atoms there are 
# ...


# Skip this step if we run DAFT
if [[ $STEP == $NOW ]]
then
    [[ -n $DAFT  ]] && : $((STEP++))
fi


if [[ $STEP == $NOW ]]
then
    # Output for this section:
    OUT=$base-cg.gro
    TOP=$base-cg.top
    NDX=$base-cg.ndx
    LOG=01-SOLVENT.log

    OUTPUT=($OUT $TOP $NDX)

    # Delete existing output if we force this step
    $FORCE && rm ${OUTPUT[@]}

    if $(all_exist ${OUTPUT[@]})
    then
	echo Found $TOP and $OUT. Skipping step with martinize/insane...
	: $((STEP++))
    fi
fi


if [[ $STEP == $NOW ]]
then 

    if [[ -z $INSA ]]
    then
	FATAL "Dependency not found: insane (building membrane/solvent)."
    fi

    INSANE="$INSA $(expandOptList ${INSANE[@]})"
    if [[ "$INSANE" =~ " -l " ]]
    then
	echo "# Calling with insane -l, indicating a membrane is included."
    fi

    if [[ -n $SOL && ! "$INSANE" =~ " -sol " ]]
    then
	note_msg="No solvent included in insane statement... Adding \"-sol $SOL\""
	store_note_fun "$note_msg"
	INSANE="$INSANE -sol $SOL"
    fi

    [[ -n $pdb ]] && INSANE="$INSANE -o $OUT -f $base-mart-EM.gro" || INSANE="$INSANE -o $OUT -p $TOP"
    
    echo "$INSANE"

    if [[ -n $PDB ]]
    then
	$NOEXEC $INSANE 2>&1 | tee -a $TOP
    else
	$NOEXEC $INSANE
    fi

    N='\'$'\n'

    # Sugar hack
    if true
    then
        # Uncomment the include topology for sugars
	cp $SDIR/martini_v2.0_sugars.itp ./
	grep -q sugars.itp $TOP && CARBFIX='/sugars.itp/s/^; //' || CARBFIX='/\[ *system *\]/{s/^/#include "martini_v2.0_sugars.itp"'"$N$N"'/;}'
	SED -i"" -e "$CARBFIX" $TOP
    fi

    if [[ $INSANE =~ "-l " ]]
    then
        # Uncomment the include topology for lipids
	cp $SDIR/martini_v2.0_lipids.itp ./
	grep -q lipids.itp $TOP && LIPFIX='/lipids.itp/s/^; //' || LIPFIX='/\[ *system *\]/{s/^/#include "martini_v2.0_lipids.itp"'"$N$N"'/;}'
	SED -i"" -e "$LIPFIX" $TOP
    else
	echo "$INSANE"
    fi

    # When multiscaling, the ions can interact with the protein on the atomistic level
    if $ALL -o [[ -n $MULTI ]]
    then
	IONSITP=$ForceField.ff/ions.itp
    else
	IONSITP=martini_v2.0_ions.itp
	cp $SDIR/$IONSITP ./
    fi

    grep -q ions.itp $TOP && IONFIX='/ions.itp/s/^; //' || IONFIX='/\[ *system *\]/{s,^,#include "'$IONSITP'"'"$N$N"',;}'
    SED -i -e "$IONFIX" $TOP 

    GRO=$OUT
fi


# Bookkeeping

NTOT=$(awk '{getline; print; exit}' $GRO)                           
NSOL=$(grep "^..... *${SOLNAMES[$SID]} " $GRO | wc -l)
NION=$(grep '^..... *\(NA[+]* \+\|CL[-]* \+\)' $GRO | wc -l)
NMEM=$((NTOT-NSOL-NION-NPROT)) 
[[ $NMEM -gt 0 ]] && MEMBRANE=true || MEMBRANE=false

echo $GRO: NTOT=$NTOT NPROT=$NPROT NSOL=$NSOL NMEM=$NMEM NION=$NION


if [[ $STEP = $NOW ]]
then
    if [[ -n $PDB ]]
    then	
	# Add membrane and solvent to CG list in $base-mart.ndx
	# The structure has Protein/DNA/RNA Membrane Solvent Ions
	cp $base-mart.ndx $base-cg.ndx
	echo -e "$(printf "$fmt\n" `SEQ $((NPROT + 1)) $((NPROT + NMEM + NSOL))` | sed 's/ 0//g')" >> $base-cg.ndx

	# Add the ions to the CG group of we do not want hybrid ions
	if ! $HybridIons && [[ $NION -gt 0 ]]
	then
	    echo -e "$(printf "$fmt\n" `SEQ $((NPROT+NMEM+NSOL+1)) $((NPROT+NMEM+NSOL+NION))` | sed 's/ 0//g')" >> $base-cg.ndx
	fi

	echo -e "[ Solute ]\n$(printf "$fmt\n" `SEQ 1 $NPROT` | sed 's/ 0//g')" >> $base-cg.ndx
    else
	# Everything is in the CG group
	echo -e "[ CG ]\n$(printf "$fmt\n" `SEQ 1 $NTOT`      | sed 's/ 0//g')"  > $base-cg.ndx
        # And we have an empty Solute group
        echo '[ Solute ]' >> $base-cg.ndx
    fi

    # The membrane tag is always added to the index file; it may be an empty group
    echo '[ Membrane ]' >> $base-cg.ndx
    if [[ $NMEM -gt 0 ]] 
    then
	echo -e "$(printf "$fmt\n" `SEQ $((NPROT+1)) $((NPROT+NMEM))` | sed 's/ 0//g')" >> $base-cg.ndx
    fi

    # Solvent 
    echo -e "[ Solvent ]\n$(printf "$fmt\n" `SEQ $((NPROT+NMEM+1)) $((NPROT+NMEM+NSOL+NION))` | sed 's/ 0//g')" >> $base-cg.ndx

    if $HybridIons && [[ $NION -gt 0 ]]
    then
	# Make ions group (hybrid character) for multiscaling
	echo -e "[ Ions ]\n$(printf "$fmt\n" `SEQ $((NPROT+NMEM+NSOL+1)) $((NPROT+NMEM+NSOL+NION))` | sed 's/ 0//g')" >> $base-cg.ndx
    fi

    # The whole system... easy.
    echo -e "[ System ] \n$(printf "$fmt\n" `SEQ 1 $NTOT` | sed 's/ 0//g')" >> $base-cg.ndx

    # Energy groups, if we have them
    [[ -n $DAFT ]] && [[ -f $DAFT ]] && cat $DAFT >> $base-cg.ndx

    NDX=$base-cg.ndx
else
    echo "Skipping coarsegrained topology acrobatics"
fi


if ( $ALL || [[ -n $MULTI ]] ) && grep -q '\[ Ions \]' $base-cg.ndx
then
    # We need to change the moleculetype names for 
    # the ions. CG names are NA+, CL-, etc, whereas AA names
    # are NA, CL, etc
    sed -i '/\[ *molecules *\]/,$s/\(NA\|K\|CA\|ZN\|CU\|CL\)2*[+-]/\1 /' $base-cg.top
fi


if ( $ALL || [[ -n $MULTI ]] ) && grep -q '\[ Ions \]' $base-cg.ndx
then
    # Also (re)set tables and interactions for multiscaling
    cp table_AA_AA.xvg table_AA_Ions.xvg
    __mdp_ms__energygrp_table=$__mdp_ms__energygrp_table,Ions,AA
    __mdp_ms__energygrps=$__mdp_ms__energygrps,Ions
    __mdp_ms__energygrp_excl=$__mdp_ms__energygrp_excl,Ions,VZ
fi


trash *.ssd 


# END OF COARSE GRAINED/MULTISCALE TOPOLOGY

[[ $STOP ==   $NOW     ]] && exit_clean
[[ $STEP == $((NOW++)) ]] && : $((STEP++))


# END OF STEP TOPOLOGY

#[[ $STOP ==   $NOW     ]] && exit_clean
#[[ $STEP == $((NOW++)) ]] && : $((STEP++))


#--------------------------------------------------------------------
SHOUT "---STEP 2: ENERGY MINIMIZATION"
#--------------------------------------------------------------------

# DAFT runs will usually start here.
# Just to make sure we have the correct energy groups set when running with -daft:
# The index file obtained from daft.py sets everything, but also includes
# System and Solute
if [[ -f $DAFT ]]
then
    NDX=$DAFT
    cp $SDIR/martini_v2.0_{ions,lipids}.itp .
    __mdp_cg__energygrps=$(sed 's/ /,/g' <<< `sed '/^ *[^[]/d;s/\[ *\(.*\) *\]/\1/;/Solute/d;/System/d;' $DAFT`)
fi

trash $base-EM.{tpr,edr,trr} em-out.mdp

MDP=em-sol.mdp
OPT=(em $MDPMS)
OUT=$base-EM.gro
LOG=02-EM.log

mdp_options ${OPT[@]} > $MDP
MD="MDRUNNER -f $MDP -c $GRO -p $TOP -o $OUT -n $NDX -np 1 -l $LOG -force $FORCE $TABLES"
echo $MD

[[ $STEP ==   $NOW     ]] && $NOEXEC $MD && : $((STEP++)) && GRO=$OUT && archive
[[ $STOP == $((NOW++)) ]] && exit_clean

# Here we check whether energy minimization was successful. If not, the potential energy or
# force will probably have gone infinite. In that case the run is terminated to avoid 
# hanging of the following run.
grep "^\(Pot\|Max\|Norm\).*inf" $LOG && echo Energy minimization failed && exit_error ${STEP}5

# To avoid clutter we delete the 'step' files, rather than trashing them
rm step*.pdb


#--------------------------------------------------------------------
SHOUT "---STEP 3: POSITION RESTRAINED NVT"
#--------------------------------------------------------------------

# From here on temperature coupling is needed
GRPS=
TTAU=
TTMP=
for i in Solute Membrane Solvent
do
    grep -q "\[ *$i *\]" $NDX && GRPS="$GRPS,$i" && TTAU="$TTAU,1.0" && TTMP="$TTMP,$Temperature"
done 
__mdp_cg__tc_grps=$GRPS
__mdp_cg__tau_t=$TTAU
__mdp_cg__ref_t=$TTMP

# Base MDP options for all cycles
OPT=(cg $MDPMS equil)

# Running with PW on more than 4 cores during equilibration is problematic
[[ $NP -gt 65535 ]] && PRNP=4 || PRNP=$NP

# This step is skipped for DAFT runs with normal water.
echo This step is skipped if DAFT is set and normal water is used. DAFT is set to "${DAFT:-nothing, so this step is executed}."
if [[ -z $DAFT ]] || $POLARIZABLE
then
    LOG=03-PR-NVT.log
    OUT=$base-PR-NVT.gro
    MDP=pr-nvt.mdp
    mdp_options ${OPT[@]} > $MDP
    MD="MDRUNNER -f $MDP -c $GRO -p $TOP -o $OUT -n $NDX -np $PRNP -l $LOG -force $FORCE $TABLES"
    $NOEXEC $MD 
    GRO=$OUT
    trash $base-PR-NVT.{cpt,tpr} 
fi

[[ $STEP ==   $NOW     ]] && : $((STEP++)) && archive 
[[ $STOP == $((NOW++)) ]] && exit_clean


#--------------------------------------------------------------------
SHOUT "---STEP 4: POSITION RESTRAINED NpT"
#--------------------------------------------------------------------

# From here on pressure coupling is needed
__mdp_cg__pcoupl=Berendsen


# Checking what to do
if [[ -n $MDPMS ]]
then
    # Multiscaled - Using virtual sites or not
    $VSITE && DT=(0.0005 0.001 0.002 0.003 0.004) || DT=(0.0005 0.001 0.002)
else
    # Coarsegrained 
    # "With DAFT runs, we just skip this step (PR-NpT)"
    [[ -n $DAFT ]] && DT=() || DT=(0.005 0.010 0.020)
fi


# Base MDP options for all cycles
$MEMBRANE && OPT=(cg mem $MDPMS equil npt) || OPT=(cg $MDPMS equil npt)
# Isotropic pressure coupling for membranes:
$MEMBRANE && __mdp_equil__tau_p=$__mdp_equil__tau_p,$__mdp_equil__tau_p
for __mdp_equil__dt in ${DT[@]}
do
    LOG=04-PR-NPT-$__mdp_equil__dt.log
    OUT=$base-PR-NPT-$__mdp_equil__dt.gro
    MDP=pr-npt-$__mdp_equil__dt.mdp
    mdp_options ${OPT[@]} > $MDP
    MD="MDRUNNER -f $MDP -c $GRO -p $TOP -o $OUT -n $NDX -np $NP -l $LOG -force $FORCE $TABLES"
    $NOEXEC $MD 
    GRO=$OUT
    trash $base-PR-NPT-$__mdp_equil__dt.{cpt,tpr} 
done

[[ $STEP ==   $NOW     ]] && : $((STEP++)) && archive
[[ $STOP == $((NOW++)) ]] && exit_clean


#--------------------------------------------------------------------
SHOUT "---STEP 5: UNRESTRAINED NpT"
#--------------------------------------------------------------------

$MEMBRANE && OPT=(cg mem $MDPMS equil usr) || OPT=(cg $MDPMS equil usr)
GRO=$OUT
OUT=$GRO
LOG=05-NPT.log
MDP=npt.mdp

# Clear position restraints
__mdp_equil__define=

# Checking what to do
if [[ -n $MDPMS ]]
then
    # Multiscaled - Using virtual sites or not
    $VSITE && DT=(0.004 $DELT) || DT=(0.002 $DELT)
else
    # Coarsegrained
    DT=(0.020 $DELT)
fi

for __mdp_equil__dt in ${DT[@]}
do
    LOG=05-NPT-$__mdp_equil__dt.log

    OUT=$base-NPT-$__mdp_equil__dt.gro
    MDP=md-init-$__mdp_equil__dt.mdp
    mdp_options ${OPT[@]} > $MDP
    nsteps=$(sed -n '/nsteps/s/.*=//p' $MDP)
    before=$(date +%s)
    MD="MDRUNNER -f $MDP -c $GRO -p $TOP -o $OUT -n $NDX -np $NP -l $LOG -force $FORCE $TABLES"
    $NOEXEC $MD 
    timing=$(( $(date +%s) - before ))
    echo "Ran $nsteps steps in $timing seconds"
    GRO=$OUT
    trash $base-NPT-$__mdp_equil__dt.{cpt,tpr} 
done


# : $((STEP++))

[[ $STEP ==   $NOW     ]] && : $((STEP++)) && archive
[[ $STOP == $((NOW++)) ]] && exit_clean


#--------------------------------------------------------------------
SHOUT "---STEP 6: PRODUCTION RUN"
#--------------------------------------------------------------------

$MEMBRANE && OPT=(cg mem $MDPMS usr) || OPT=(cg $MDPMS usr)
GRO=$OUT
OUT=$base-MD.gro
LOG=06-MD.log
MDP=md.mdp

# If the run length is set to zero, write an mdp file 
# with nsteps=-1, generate the .tpr, and make a clean exit
if [[ $TIME == 0 ]]
then
    __mdp_cg__nsteps=-1
    OUT=$base-MD.tpr
fi

# Set pressure coupling to PR if running with GMX5
if [[ $GMXVERSION -gt 4 ]]
then
    __mdp_cg__pcoupl=Parrinello-Rahman
    __mdp_cg__tau_p=12.0
    __mdp_mem__tau_p=12.0,12.0
fi

echo ${STEPS[$STEP]} ${STEPS[$STOP]} ${STEPS[$NOW]}

# If the Stop-Step is set to TPR only write a tpr and do not run
[[ $STEP ==   $NOW     ]] && : $((STEP++))
[[ $STOP == $((NOW++)) ]] && OUT=$base-MD.tpr

echo @@@ $OUT

mdp_options ${OPT[@]} > $MDP

mdsteps=$(sed -n '/nsteps/s/.*=//p' $MDP)
estimate=$(( (timing*(1000*mdsteps)/nsteps)/1000 ))
estfin=$(( $(date +%s) + estimate ))
echo "Expected runtime for $mdsteps step: $(( estimate/3600 ))H:$(( (estimate%3600)/60 ))M:$(( estimate%60 ))S (until $(date -r $estfin))"

MD="MDRUNNER -f $MDP -c $GRO -p $TOP -o $OUT -n $NDX -np $NP -l $LOG -split -force $FORCE $TABLES"
$NOEXEC $MD 

# : $((STEP++))

[[ $STEP ==   $NOW     ]] && : $((STEP++)) && archive
[[ $STOP == $((NOW++)) ]] && exit_clean

exit_clean
