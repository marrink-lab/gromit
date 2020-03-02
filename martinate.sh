#!/bin/bash

PROGRAM=martinate.sh
VERSION=0.99
VERSTAG=devel-180428-1200-TAW
AUTHORS="Tsjerk A. Wassenaar, PhD"
YEAR="2018"
AFFILIATION="
University of Groningen
Nijenborgh 7
9747AG Groningen
The Netherlands"


: << __NOTES_FOR_USE_AND_EDITING__

IF YOU CHANGE THE PARAMETERS AND/OR WORKFLOW, PLEASE RENAME THE PROGRAM AND
STATE THE NATURE AND PURPOSE OF THE CHANGES MADE.

This has grown to be a rather complicated bash script. It is intended to
work through the MD process as a person would, issuing shell commands and reading
and editing files. Bash feels more natural for this than a Python/C wrapper.
It is advised to (get to) know about bash loops and variable substitution, as
these are used plenty. In addition, since there are many occassions where files
need to be read and edited, there are a lot of calls to sed, with quite a
few less ordinary commands.

To keep the code manageable, it is structured in sections and every section is
ordered, preferrably by numbered chunks. In addition, there is extensive
documentation. Every statement should be clear, either by itself or by a
preceding explanation. In case advanced bash/sed/... features are used, they
ought to be explained. That will keep the program manageable and make it a nice
place for learning tricks :)

Oh, and please note that usual copyright laws apply...

TAW - 20120718
(copied from gromit.sh)

__NOTES_FOR_USE_AND_EDITING__


DESCRIPTION=$(cat << __DESCRIPTION__

$PROGRAM $VERSION is a versatile wrapper built around
GROMACS, insane, and martinize, for setting up and running
MARTINI COARSE GRAIN molecular dynamics simulations of solvents,
membranes, proteins and/or nucleic acids in any combination.

It is built to allow automated processing of membrane proteins,
with full control of membrane composition. If no input file is provided,
only a membrane and/or solvent is built.
Options given that do not match an option in this script are passed to
martinize.py.

The script contains a complete and flexible workflow, consisting of the
following steps:

    1.   Generate topology from input structure
         A. Generate atomistic topology             (AA)
         B. Generate MARTINI CG/multiscale topology (CG)
    2.   Solvation and adding ions                  (SOLVENT)
    5.   Energy minimization                        (EM)
    6.   Position restrained NVT equilibration      (NVT-PR)
    7.   Unrestrained NpT equilibration             (NPT)
    8.   Equilibration under run conditions         (PREPRODUCTION)
    9.   Production simulation
         A. Run input file                          (TPR)
         B. Simulation (possibly in parts)          (PRODUCTION)

The program allows running only part of the workflow by specifying the
start and end step (-step/-stop), using an argument uniquely matching
one of the tags given between parentheses.

This program requires a working installation of Gromacs. To link
the program to the correct version of Gromacs, it should be placed in the
Gromacs binaries directory or the Gromacs GMXRC file should be passed as
argument to the option -gmxrc

The workflow contained within this program corresponds to a standard protocol
that should suffice for routine CG molecular dynamics simulations of proteins
and/or nucleic acids in aqueous solution with or without a membrane.
It follows the steps that are commonly taken in MD tutorials
(e.g. those at http://cgmartini.nl).

This program is designed to enable high-throughput processing of CG molecular
dynamics simulations in which specific settings are varied systematically. These
settings include protein/nucleic acid, ligand, temperature, and pressure, as well
as many others.


## -- IMPORTANT -- ##

Molecular dynamics simulations are complex, with many contributing factors.
The workflow in this program has been tested extensively and used many times.
Nonetheless, it should not be considered failsafe. No MD protocol ever is.
Despite careful set up, simulations may crash, and the possibility that a crash
is encountered is larger when many simulations are run. If the run crashes,
the intermediate results will be kept and can be investigated to identify the
source of the problem.

If the run finishes to completion, this does not automatically imply that the
results are good. The results from the simulations should always be subjected
to integrity and quality assurance checks to assert that they are correct within
the objectives of the study.

__DESCRIPTION__
)

#--------------------------------------------------------------------
#---Parsing COMMAND LINE ARGUMENTS AND DEPENDENCIES--
#--------------------------------------------------------------------

CMD="$0 $@"
echo "$CMD" | tee CMD


# Directory where this script is
SDIR=$( [[ $0 != ${0%/*} ]] && cd ${0%/*}; pwd )
SRCDIR="$SDIR/source"
FFDIR="$SDIR/forcefield"

# Sourcing modules
source "$SRCDIR"/_logging.sh
source "${SRCDIR}"/_optionhandling.sh
source "${SRCDIR}"/_functions.sh
source "${SRCDIR}"/_mdp_martinate.sh
source "${SRCDIR}"/_mdp.sh
source "${SRCDIR}"/_mdrunner.sh
source "${SRCDIR}"/_pdb.sh
source "${SRCDIR}"/_martinate_index.sh
source "${SRCDIR}"/_martinate_multiscale.sh
source "${SRCDIR}"/_martinate_daft.sh
source "${SRCDIR}"/_gmx.sh

trap "archive" 2 9 15 # function archive in _functions.sh


# These will be looked for before running, and can be set from the cmdline, e.g.:
#    -gmxrc /usr/local/gromacs-5.1/bin/GMXRC
# If not set, the default name will be searched for in
#    1. the environment (if PROGEVAR is given)
#    2. the directory where this calling script (martinate) is located
#    3. the PATH
DEPENDENCIES=( dssp  gmxrc  martinize     insane     liptop           squeeze)
PROGEXEC=(     dssp  GMXRC  martinize.py  insane     $FFDIR/liptop.py squeeze)
PROGEVAR=(     DSSP  GMXRC)


# Run control
MONALL=       # Monitor all steps
CONTROL=
CHECKTIME=300 # Run control every five minutes


# Stepping stuff
STEPS=(AA CG SOLVENT EM NVT-PR NPT PREPRODUCTION TPR PRODUCTION ANALYSIS END)
get_step_fun() { for ((i=0; i<${#STEPS[@]}; i++)) do [[ ${STEPS[$i]} =~ ^$1 ]] && echo $i; done; }
STEP=AA
STOP=PRODUCTION

# MARTINI Force field parameters
MARTINI=martini22
DRY=
FFITP=
FFTAG=v2.0
USRITP=() # Additional topologies to include

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
VirtualSites=false

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

# - run control and files
DIR="."       # Directory to run and write
TPR=          # Run input file... skip to production or to STEP
NAME=         # Run name
FETCH=        # Try to fetch PDB file from web
MSGFILE=/dev/stdout      # Master log file (stdout)
ERRFILE=/dev/stderr      # Error log file (stderr)
EXEC=         # Execute run
NP=1          # Number of processors
MDP=          # User-provided MDP file
MDARGS=       # User-provided arguments to mdrun
MAXH=-1       # Maximum duration of run
JUNK=()       # Garbage collection
SCRATCH=      # Scratch directory
ARCHIVE=      # Archive file name
FORCE=false   # Overwrite existing run data
KEEP=false    # Keep intermediate rubbish (except junk)
DAFT=

# - system setup
Salinity=0.1536  # Salt concentration of solvent


# - group definitions
NATOMS=0                 # Total number of atoms
Biomol=()                # Biomolecules (protein, nucleic acid, lipid)
Solute=()                # Solute molecule (protein or so)
Membrane=()              # Lipids, excluding protein
Solvent=()               # Solvent, including ions
Ligand=()                # Ligands
Ligenv=()                # Ligand environment to consider for LIE contributions
CoupleGroups=()          # Groups for temperature coupling
EnergyGroups=()          # Groups for energy calculations
LIE=false                # LIE run


# - simulation parameters
TIME=0           # Production run time (ns)
AT=0.5           # Output frequency for positions, energy and log (ns)
DELT=0.020       # Time step in picoseconds
EMSTEPS=500      # Number of steps for EM
Temperature=310  # Degree Kelvin
Tau_T=0.1                # ps
Pressure=1               # Bar
Tau_P=1.0                # ps
SEED=$$                  # Random seed for velocity generation
RotationalConstraints=   # Use rotational constraints, which is mandatory with NDLP


# User defined gromacs program options and simulation parameters (way flexible!)
PROGOPTS=()              # User-defined program options (--program-option=value)
MDPOPTS=()               # User-defined mdp parametesrs (--mdp-option=value)

hlevel=0
olevel=0

# Collect errors, warnings and notes to (re)present to user at the end
# Spaces are replaced by the unlikely combination QQQ to keep the
# messages together.
errors_array=()
store_error_fun() { a="$@"; errors_array+=(${x// /QQQ}); FATAL "$@"; }
warnings_array=()
store_warning_fun() { a=$@; warnings_array+=(${x// /QQQ}); WARN "$@"; }
notes_array=()
store_note_fun() { a=$@; notes_array+=(${x// /QQQ}); NOTE "$@"; }


##>> OPTIONS

if [[ -z "$1" ]]; then
  echo "No command line arguments give. Please read the program usage:"
  USAGE 1
  exit
fi


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
    #=0
    #=0 OPTIONS
    #=0 =======
    #=0
    -h       ) USAGE 0                              ; exit 0 ;; #==0 Display help
    --help   ) USAGE 0                              ; exit 0 ;; #==1 Display help
    -hlevel  ) hlevel=$2                            ; shift 2; continue ;; #==1 Set level of help (use before -h/--help)
    -olevel  ) olevel=$2                            ; shift 2; continue ;; #==1 Set level of options to display

    #=1
    #=1 File options
    #=1 ------------
    #=1
    -f       ) PDB=$2                               ; shift 2; continue ;; #==0 Input PDB file
    -g       ) MSGFILE=$2                           ; shift 2; continue ;; #==9 Standard output log file (default: /dev/stdout)
    -e       ) ERRFILE=$2                           ; shift 2; continue ;; #==9 Standard error log file (default: /dev/stderr)
    -name    ) NAME=$2                              ; shift 2; continue ;; #==9 Name of project
    -top     ) TOP=$2                               ; shift 2; continue ;; #==1 Input topology file
    -ndx     ) NDX=$2                               ; shift 2; continue ;; #==1 Input index file
    -mdp     ) MDP=$2                               ; shift 2; continue ;; #==9 MDP (simulation parameter) file
    -scratch ) SCRATCH=$2                           ; shift 2; continue ;; #==9 Scratch directory to perform simulation
    -fetch   ) FETCH=$2                             ; shift 2; continue ;; #==9 Database to fetch input structure from
    -hetatm  ) NOHETATM=false                       ; shift 1; continue ;; #==1 Keep HETATM records

    #=1
    #=1 Overall control options
    #=1 -----------------------
    #=1
    -step    ) STEP=$2                              ; shift 2; continue ;; #==1 Step to start protocol
    -stop    ) STOP=$2                              ; shift 2; continue ;; #==1 Step to stop protocol
    -keep    ) KEEP=true                            ; shift 1; continue ;; #==2 Whether or not to keep intermediate data
    -dir     ) DIR=$2                               ; shift 2; continue ;; #==2 Directory for running simulation
    -np      ) NP=$2                                ; shift 2; continue ;; #==1 Number of cores (processes) to use
    -maxh    ) MAXH=$2                              ; shift 2; continue ;; #==2 Maximum run time
    -force   ) FORCE=true                           ; shift 1; continue ;; #==2 Whether or not to force redoing parts already run
    -noexec  ) EXEC=echo                          ; shift 1; continue ;; #==2 Whether or not to actually execute the commands

    #=1
    #=1 Forcefield control options
    #=1 --------------------------
    #=1
    -cg      ) MARTINI=$2                           ; shift 2; continue ;; #==1 Coarse grain force field
    -sol     ) SOL=$2                               ; shift 2; continue ;; #==1 Solvent type to use
    -ffitp   ) FFITP=$2                             ; shift 2; continue ;; #==2 Coarse-grain force field definition
    -ffdir   ) FFDIR=$2                             ; shift 2; continue ;; #==2 Directory for force field files
    -fftag   ) FFTAG=$2                             ; shift 2; continue ;; #==2 Tag for force field files (v3.0 -> martini_v3.0_lipids.itp)
    -itp     ) USRITP+=($2)                         ; shift 2; continue ;; #==2 User-provided ITP file
    -dry     ) DRY=$2                               ; shift 2; continue ;; #==2 Use dry martini from file definition

    #=1
    #=1 Simulation control options
    #=1 --------------------------
    #=1
    -T       ) Temperature=$2                       ; shift 2; continue ;; #==1 Temperature
    -P       ) Pressure=$2                          ; shift 2; continue ;; #==1 Pressure
    -salt    ) Salinity=$2                          ; shift 2; continue ;; #==1 Salt concentration
    -dt      ) DELT=$2                              ; shift 2; continue ;; #==2 Integration time step
    -time    ) TIME=$2                              ; shift 2; continue ;; #==1 Production simulation time
    -at      ) AT=$2                                ; shift 2; continue ;; #==1 Output sampling frequency
    -em      ) EMSTEPS=$2                           ; shift 2; continue ;; #==2 Number of steps for EM
#   -gmxrc   ) GMXRC=$2                             ; shift 2; continue ;;
#   -dssp    ) DSSP=$2                              ; shift 2; continue ;;
    -rtc     ) RotationalConstraints=rtc            ; shift  ; continue ;; #==2 Whether or not to use rotational constraints

    #=2
    #=2 Multiscale options
    #=2 ------------------
    #=2
    -m       ) MULTI[$((NCH++))]=$2; M=true         ; shift 2; continue ;; #==2 Chains for multiscaling
    -M       ) ALL=true; M=true                     ; shift 1; continue ;; #==2 Multiscale all chains
    -ff      ) ForceField=$2                        ; shift 2; continue ;; #==2 Atomistic force field for multiscaling
    -vsite   ) VirtualSites=true                    ; shift 1; continue ;; #==2 Use virtual sites in multiscaling
    -epsr    ) EPSR=$2                              ; shift 2; continue ;; #==2 Dielectric constant of vacuum
    -epsrf   ) EPSRF=$2                             ; shift 2; continue ;; #==2 Dielectric constant of Reaction-Field
    -ljdp    ) LJDP=$2                              ; shift 2; continue ;; #==2 Lennard-Jones dispersion
    -ljrp    ) LJRP=$2                              ; shift 2; continue ;; #==2 Lennard-Jones repulsion
    -ljsw    ) LJSW=$2                              ; shift 2; continue ;; #==2 Lennard-Jones switch radius
    -rc      ) RC=$2                                ; shift 2; continue ;; #==2 Cut-off for non-bonded terms

    #=1
    #=1 Monitor options
    #=1 ---------------
    #=1
    -monall  ) MONALL=-monitor                      ; shift 1; continue ;; #==2 Monitor all steps using control script
    -control )                                                             #==2 Simulation monitor script
      while [[ -n $2 && $2 != ';' ]]
      do
        CONTROL="$CONTROL $2"
        shift
      done
      shift 2
      echo MONITOR COMMAND: $CONTROL
      continue;;
    -ctime   ) CHECKTIME=$2                         ; shift 2; continue ;; #==2 Time for running monitor

    #=2
    #=2 A control process is either a program, script or command
    #=2 that monitors the production run and terminates it
    #=2 upon a certain condition, indicated by its exit code.

    #=2
    #=2 Protocol options
    #=2 ----------------
    #=2
    -daft    ) DAFT=$2; NDX=$2                      ; shift 2; continue ;; #==2 Run martinate in DAFT pipeline

    #=1
    #=1 Advanced control options
    #=1 ------------------------
    #=1
    #=2 This program allows specifying options for advanced control of
    #=2 program invocation and simulation parameters. These options are
    #=2 described below.
    #=2

    # The first one is the template/dummy for the help system
    --mdp-option=value) olevel=2; hlevel=2; USAGE 1; continue;; #==2 Command-line specified simulation parameters
    --mdp-*  ) MDPOPTS+=(${1#--mdp-})               ; shift  ; continue ;;
    #=2
    #=2 This will add 'option = value' to the MDP file for all simulations
    #=2 following energy minimization. MDP options specified on the command line
    #=2 take precedence over those specified in an input file (-mdp), which take
    #=2 precedence over parameters defined in this script.
    #=2 If the option takes multiple arguments, then 'value' should be a
    #=2 comma separated list.
    #=2 The STEP/STOP controls can be used to set parameters for (pre)production
    #=2 simulations selectively.
    #=2

    # Options for downstream programs
    # If the options are given on the command line, they are expanded and each
    # option will be formatted as --program-opt=val
    --martinize-option=value) olevel=2; hlevel=2; USAGE 1; continue;; #==2 Parameters for martinize
    --martinize-*) MARTINIZE+=(${1#--martinize})    ; shift  ; continue;;

    --insane-option=value) olevel=2; hlevel=2; USAGE 1; continue;; #==2 Parameters for insane
    --insane-*)       INSANE+=(${1#--insane})       ; shift  ; continue;;

    # If the options are passed by another program, they will be formatted like
    #   --program{-opt1=val1,-opt2=val2\,val3}
    # In this case the option needs to be parsed explicitly:
    --martinize*)  MARTINIZE+=($(readOptList $1))    ; shift  ; continue;;
    --insane*)        INSANE+=($(readOptList $1))    ; shift  ; continue;;

    # Other program-specific options
    --*)   PROGOPTS[${#PROGOPTS[@]}]=$1              ; shift 1; continue ;;

    #=0
    #=0

    # All options should be covered above. Anything else raises an error here.
    *)         BAD_OPTION "$1";;

  esac
done
##<< OPTIONS


#--------------------------------------------------------------------
#---GLOBAL PARAMETERS AND STUFF--
#--------------------------------------------------------------------

exec 3>&1 4>&2
[[ -n $MSGFILE ]] && exec 1>$MSGFILE
[[ -n $ERRFILE ]] && exec 2>$ERRFILE

cat << __RUNINFO__

$PROGRAM version $VERSION:

(c)$YEAR $AUTHOR
$AFFILIATION

Now executing...

$CMD

__RUNINFO__

echo $CMD > cmd.log

# Time. To keep track of the remaining run time
START=$(date +%s)


# START/STOP FLOW CONTROL
for ((i=0; i<${#STEPS[@]}; i++)); do [[ ${STEPS[$i]} == ${STEP}* ]] && STEP=$i && break; done
for ((i=0; i<${#STEPS[@]}; i++)); do [[ ${STEPS[$i]} == ${STOP}* ]] && STOP=$i && break; done


# Set the scratch directory, if any:
#    scratch directory, user name, random number
if [[ -n $SCRATCH ]]
then
    # The scratch directory can be specified as
    # (escaped) variable, like \$TMPDIR. This
    # variable will be expanded at runtime.
    # That may be handy on clusters, where the
    # $TMPDIR is set for every node.
    if [[ ${SCRATCH:0:1} == '$' ]]
    then
	tmp=${SCRATCH:1}
	SCRATCH=${!tmp}
    fi

    # To ensure that there is no further rubbish
    # the scratch directory is extended with the
    # username, the data and the process ID. There
    # the run will be performed.
    SCRATCH=$SCRATCH/$(date +%F).$USER.$$
    if ! mkdir -p $SCRATCH
    then
	echo Scratch directory $SCRATCH is not available... exiting
	exit
    fi

    echo $SCRATCH > SCRATCH
fi


#--------------------------------------------------------------------
#---Sed and awk
#--------------------------------------------------------------------

# Awk expression for extracting moleculetype
#    - at the line matching 'moleculetype'
#      read in the next line
#      continue reading next lines until one is not beginning with ;
#      print the first field
AWK_MOLTYPE='/moleculetype/{getline; while ($0 ~ /^ *;/) getline; print $1}'


#--------------------------------------------------------------------
#---GROMACS AND RELATED STUFF
#--------------------------------------------------------------------

## 0. Finding programs

dependency_not_found_error()
{
  FATAL The required dependency $@ was not found.
}

NDEP=${#DEPENDENCIES[@]}
find_program_function()
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
  # Python scripts may be available as 'binaries' (martinize/insane)
  which $progr 2>/dev/null && return 0
  which ${progr%.py} 2>/dev/null && return 0 || return 1
}


##  1. GROMACS  ##

load_gromacs

## 2. DSSP ##

# Search the DSSP binary, from environment, from path, or guess
# Only required if we have an input file
SOLSTEP=$(get_step_fun CG)
if [[ $STEP -le $SOLSTEP && $STOP -ge $SOLSTEP ]]
then
  echo -n '# Checking DSSP binary (for martinizing proteins)... '
  DSSP=$(find_program_function dssp)
  if [[ $? == 1 ]]
  then
    warn="DSSP binary not found - Will martinize without secondary structure :S"
    store_warning_fun "$warn"
  else
    echo "$DSSP"
    MARTINIZE+=(-dssp=$DSSP)
  fi
fi



## 4. Echo mdp/program options specified on the command line

MSG="Program options specified on command line:"
echo_additional_options ${PROGOPTS[@]}

MSG="MDP options specified on command line (note how flexible!):"
echo_additional_options ${MDPOPTS[@]}

## 5. Locate insane if STEP lies before SOLVENT and STOP lies after.

SOLSTEP=$(get_step_fun SOLVENT)
if [[ $STEP -le $SOLSTEP && $STOP -ge $SOLSTEP ]]
then
    INSA=$(find_program_function insane)
    if [[ $? != 0 ]]
    then
	STEP=$NOW
	FATAL "Dependency (insane) required for building solvent/membrane, but not found."
    fi
fi


#--------------------------------------------------------------------
#---TIMING
#--------------------------------------------------------------------

# Maximum time in seconds
if [[ $MAXH =~ ":" ]]
then
  # Format HH:MM:SS
  ifs=$IFS; IFS=":"; MAXS=($MAXH); IFS=$ifs
  MAXS=$((3600*MAXS[0] + 60*MAXS[1] + MAXS[2]))
else
  # Format x.y HH
  MAXS=$(awk '{printf "%d\n", $1*3600}' <<< $MAXH )
fi

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

## 2. WORKING DIRECTORY AND SOURCE DIRECTORY ##
# SRCDIR=$(pwd)
[[ ! -d $DIR ]] && mkdir -p $DIR; pushd $DIR >/dev/null


## 3. START/STOP FLOW CONTROL ##
NOW=$STEP
echo "# Will run from step ${STEPS[$STEP]} until ${STEPS[$STOP]}"


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

#     a. ATOMISTIC
#        Set pointers to ffnonbonded.itp and ffbonded.itp if we do multiscaling
$M && ffnb=$GMXLIB/$ForceField.ff/ffnonbonded.itp || ffnb=
$M && ffbn=$GMXLIB/$ForceField.ff/ffbonded.itp    || ffbn=

#     b. COARSE GRAINED (MULTISCALE) ITP
#        i. If a forcefield ITP is given, use that
#        ii. If a forcefield ITP generating script is available use that
#        iii. If $FFDIR contains a suitable forcefield ITP use that (and warn)
#        iv. Raise an error
FFMARTINIPY=
if [[ -n $FFITP ]]
then
  # i.
  cp $FFITP ./martini.itp
else
  FFMARTINIPY=$FFDIR/${MARTINI}${SOLVFF[$SID]}.py
  if [[ ! -f $FFMARTINIPY ]]
  then
    if [[ -f $FFDIR/${MARTINI}.py ]]
    then
      # If martini22p was specified in stead of martini22 with PW,
      # then we end up here, setting the script to martini22p.py
      FFMARTINIPY=$FFDIR/${MARTINI}.py
    fi
  fi

  if [[ -f $FFMARTINIPY ]]
  then
    # ii.
    if [[ -n $DRY ]]
    then
      $FFMARTINIPY "$DRY" > martini.itp
    else
      $FFMARTINIPY $ffnb $ffbn > martini.itp
      # UPDATE martini.itp FOR DUMMIES
      # Replace the #include statement for ff_dum.itp for atomistic force fields by the contents of it
      $M && $SED -i -e "/#include \"ff_dum.itp\"/r$GMXLIB/$ForceField.ff/ff_dum.itp" -e "/#include \"ff_dum.itp\"/d" martini.itp
    fi
  else
    # iii.
    # Check if an FF include exists in $FFDIR
    FFITP=$(
      for itp in $FFDIR/*itp
      do
        tags=$(grep '\[' $itp)
        [[ "$tags" =~ defaults && "$tags" =~ atomtypes && "$tags" =~ nonbond_params ]] && echo $itp
      done
    )
    if [[ -n $FFITP ]]
    then
      NOTE Found force field include file $FFITP
      cp $FFITP martini.itp
    fi
  fi
fi

if [[ ! -f martini.itp ]]
then
  FATAL Could not find forcefield itp file or generating script.
fi

# TODO - get something smarter for the FFTAG
if [[ "$MARTINI" =~ "martini3" ]]
then
    FFTAG=v3.0
fi

## 6. ELECTROSTATICS AND TABLES ##

EPSR_CG=${EPSR_CG[$SID]}
EPSR_AA=${EPSR:-${EPSR_AA[$SID]}}
$M && TABLES=-tables || TABLES=



#--------------------------------------------------------------------
#---SIMULATION PARAMETERS--
#--------------------------------------------------------------------

## OT N ## For every parameter not defined the default is used
## NOTE ## This is probably fine for equilibration, but check the defaults to be sure
## E OT ## The list as is was set up for gromacs 4.5 and 5.1

init_mdp_parameters
read_mdp_file
read_mdp_options

#--------------------------------------------------------------------


#--------------------------------------------------------------------
#---SUBROUTINES--
#--------------------------------------------------------------------


ERROR=0

# Always ECHO the first line
NOW=$STEP

#--------------------------------------------------------------------
SHOUT "---= THIS IS WHERE WE START =---"
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#---INPUT CHECKING, SPLITTING, TRIMMING, GROOMING
#--------------------------------------------------------------------

if [[ -n $PDB ]]
then
  # If the input file is not found, check whether it was given without extension.
  # If that is not the case, then fetch the file from the PDB repository.
  if [[ ! -f $PDB ]]
  then
    PDB=$PDB.pdb
    if [[ ! -f $PDB ]]
    then
      # Try fetching it from the PDB
      pdb=$(tr [A-Z] [a-z] <<< ${PDB%.pdb})
      fetch_structure $pdb $FETCH
      [[ -n $SCRATCH ]] && cp $pdb.pdb $DIR
    fi
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
      NOTE Removing HETATM entries from PDB file $PDB
      $SED -i'' -e "/^HETATM/d" "$PDB"
    fi

    # Extract a list of chains from PDB file
    CHAINS=( $(grep '^\(ATOM\|HETATM\)' $PDB | cut -b 22 | uniq) )

    # Unpack lists of chains to multiscale separated by commas
    MULTI=( $(for i in ${MULTI[@]}; do echo ${i//,/ }; done) )

    # Residues defined in martinize.py
    AA=(ALA CYS ASP GLU PHE GLU HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR)
    # Sed query for residues (separated by \|):
    SED_AA=$($SED 's/ \+/\\\|/g' <<< ${AA[@]})
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
    DEF=$($SED -n -e "$RTPENTRIES" -e "$FORMAT" $GMXLIB/$ForceField.ff/*.rtp)

    # Now we can split the input PDB file into a processable and a non-processable part
    echo "# Defined building blocks in $ForceField RTP files:"
    echo "# $DEF"
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


echo "# Done checking"
echo "# Done gymnastics"


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
if ! $M && [[ $STEP == $NOW ]]
then
    echo "Skipping step... (not multiscaling)"
    : $((STEP++))
fi


# Skip this step if we run DAFT
if [[ -n $DAFT ]] && [[ $STEP == $NOW ]]
then
    echo DAFT run. Skipping step.
    : $((STEP++))
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
    AAFF=($(ls -d $GMXLIB/*.ff | $SED 's#.*/\(.*\)\.ff#\1#'))

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
    TABLE="${SRCDIR}"/_table.py

    #      epsilon_r epsilon_rf cutoff LJ_dispersion LJ_repulsion LJ_cutoff LJ_switch
    $TABLE  $EPSR_CG   $EPSRF     $RC      $LJDP         $LJRP         1.2       0.9   > table.xvg
    $TABLE  1          $EPSRF     $RC      6             12            $RC      -1     > table_AA_AA.xvg
    $TABLE  1          $EPSRF     $RC      6             12            $RC      -1     > tablep.xvg
    if $POLARIZABLE
    then
	$TABLE $EPSR_AA $EPSRF    $RC      6             12            $RC      -1     > table_AA_CG.xvg
    fi
fi




## I. pdb2gmx

# 1. Basic stuff
PDB2GMX="${GMX}pdb2gmx -v -f $dirn/$base.pdb -o $OUT -p $TOP"

# 2. Position restraints
#    * The position restraint fc (-posrefc) is bogus and
#      intended to allow easy replacement with sed.
#      These will be placed under control of a #define
PDB2GMX="$PDB2GMX -i $base-posre.itp -posrefc 200 -ignh -ff $ForceFieldAA -water none"

# 3. Virtual sites
$VirtualSites && PDB2GMX="$PDB2GMX -vsite hydrogens"

# 4. Add program options specified on command line (--pdb2gmx-option=value)
PDB2GMX="$PDB2GMX $(program_options pdb2gmx)"

# 5. Specification of protonation states
if [[ $STEP == $NOW ]]
then
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
	
	$SED -e "$NTER" -e "$CTER" -e "$ACID" -e "$LYS" -e "$HIS" $dirn/$TITR > pdb2gmx.query
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
    if [[ -n $EXEC ]] || $(all_exist ${OUTPUT[@]})
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
	MTP=($(for i in ${ITP[@]}; do awk "${awk_itp_get_moleculetypes}" $i; done))

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
                    LSED -i -e "/${ITP[$i]}/d" -e "/[\s*system\s*]/,\$s/${MTP[$i]}/${MTP[$j]}/" $base-aa.top
                    # List the file for removal
                    trash ${ITP[$i]} $($SED 's/_/-posre_/' <<< ${ITP[$i]})
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
	CHAINS=($($SED -n '/chain  #res #atoms/,/^\s*$/s/^\s\+[0-9]\+..\(.\).*$/\1/p' 01-TOPOLOGY-AA.log))
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
		LSED -i.bck \
		    -e "/^ *\[ *moleculetype *\] */{h;s/.*/#include \"$ITP\"/p;x;}" \
		    -e '/moleculetype/,/^#endif/w'$ITP \
		    -e '/moleculetype/,/^#endif/d' \
		    $base-aa.top
	    fi
	fi

	MOLECULES=($(awk "${awk_top_get_moleculelist}" $base-aa.top))
	ITPFILES=($($SED -n -e "${sed_top_get_includes}" $base-aa.top))
	MOLTYPES=($(for i in  ${ITPFILES[@]}; do awk "${awk_itp_get_moltypes}" $i; done))
	echo MOLECULES: ${MOLECULES[@]}	

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
    $ALL && M_MULTI="-multi all" || M_MULTI=$($SED 's/\(^\| \)/ -multi /g' <<< ${MULTI[@]})	

    PDB=$OUT
fi

# Here we may have $OUT (if multiscaling)

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

NPROT=
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
	MART=$(find_program_function martinize)
	if [[ $? != 0 ]]
	then
	    FATAL "Coarse graining PDB file ($pdb), but martinize was not found."
	fi
	# If multiscaling, let martinize generate the index
	# martinize2 cannot generate it, but cannot do multiscaling anyway
	[[ -n $M_MULTI ]] && MNDX="-n $base-mart.ndx"
	MARTINIZE="$MART $(expandOptList ${MARTINIZE[@]})"
	MARTINIZE="$MARTINIZE -f $pdb -o $TOP -x $base-mart.pdb $MNDX -ff $MARTINI $M_MULTI"
	echo $MARTINIZE
    fi

    # Only if we have a pdb file and we actually run (EXEC is not set)
    # then this block is executed
    if [[ -n $pdb && -z $EXEC ]]
    then
        # Executing the command.
        # We need to know the moleculetype names (and itp files) after martinizing
        # 1. The command is executed
        # 2. stdout and stderr are swapped
        # 3. stdout (originally stderr) is also written to stderr using 'tee'
        # 4. stdout is parsed by sed, extracting moleculetype names for each chain
        # 5. The result is stored as array MARMOLS
        #          |-1-| |------2-----|  |-------3------|   |--------------------------4---------------------|
	MARMOLS=($($MARTINIZE 3>&1 1>&2 2>&3 | tee /dev/stderr | $SED -n '/^INFO *[0-9]\+-> \+\([^ ]\+\).*$/s//\1/p'))

	# If martinize did not generate the index file, we need to do it now
	[[ -z $M_MULTI ]] && { echo "[ CG ]"; seq $(grep -c "^ATOM" $base-mart.pdb); } > $base-mart.ndx

        # If we use polarizable or BMW water, we need to change AC1/AC2 to C1/C2
	if $POLARIZABLE
	then
	    $SED -i -e 's/AC[12]/ C[12]/' $base-mart.pdb
	fi

        # Check for chains to multiscale
	if $ALL -o [[ -n $MULTI ]]
	then
            multiscale_topology
	elif [[ -n $DAFT ]]
	then
            daft_groups
	fi

	touch residuetypes.dat elements.dat
	trash residuetypes.dat elements.dat
	${GMX}editconf -f $base-mart.pdb -o $base-mart.gro -box 100 100 100 -noc >/dev/null 2>&1

        # Have to energy minimize output structure for stability
	#__mdp_mart__pbc=no
	MDP=em-mart.mdp
	OPT=(em $MDPMS mart)
	OUT=$base-mart-EM.gro
	
	mdp_options ${OPT[@]} > $MDP
	MDRUNNER -f $MDP -c $base-mart.gro -p $TOP -o $OUT -n $base-mart.ndx -np 1 -l $LOG -force $FORCE $TABLES
	echo 0 | ${GMX}trjconv -s $base-mart.gro -f $OUT -pbc nojump -o $base-mart-nj.gro >/dev/null 2>&1
	mv $base-mart-nj.gro $OUT

        # trash $base-mart.gro
	GRO=$OUT
    fi
fi

# If $GRO is not set, we set it equal to $pdb
[[ -z $GRO ]] && GRO=$pdb

# NPROT may be set, but if we skipped this step it's not
# If it is not set, then the input could be solute (prot/nucl) and/or membrane
NPROT=$(LSED -n '2{p;q;}' "$GRO")

# END OF COARSE GRAINING
[[ $STOP ==   $NOW     ]] && exit_clean
[[ $STEP == $((NOW++)) ]] && : $((STEP++))


#---------------------------------------------------------------------
SHOUT "---STEP 1C: SOLVATE"
#---------------------------------------------------------------------

# At this point, we need to know how many atoms there are
# ...


# If this step is skipped, the structure provided has
# membrane and/or solvent. Then an index file must be
# provided, which we parse for the groups.


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

    OUTPUT=($OUT $NDX)

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
	# Add a comment at the end of the TOP file to avoid adding stuff to existing lines
	echo ';' >> $TOP
	$EXEC $INSANE 2>&1 | tee -a $TOP
    else
	$EXEC $INSANE 2>insane.stderr
	cat insane.stderr | tee -a $TOP
    fi


    # For each lipid built with insane (-al* options), make a topology:
    LIPTOP=$(find_program_function liptop)
    [[ $? != 0 ]] && NOLIPTOP=true || NOLIPTOP=false
    INS=($INSANE)
    USRDEF=()
    USRLIP=()
    for ((i=1; i<${#INS[@]}; i++))
    do
	case ${INS[$i]} in
	    -alhead) USRDEF+=(-he ${INS[$((++i))]}); continue;;
	    -allink) USRDEF+=(-li ${INS[$((++i))]}); continue;;
	    -altail) USRDEF+=(-ta ${INS[$((++i))]}); continue;;
	    -alname) USRDEF+=(-name ${INS[$((++i))]} -o ${INS[$i]}.itp); USRLIP+=(${INS[$i]}); continue;;
	esac
    done
    echo ${USRDEF[@]}
    for ((i=0; i<${#USRDEF[@]}; i+=10))
    do
	LIPDEF=${USRDEF[@]:$i:10}
	echo $LIPTOP ${LIPDEF[@]}
	$LIPTOP ${LIPDEF[@]}
    done
    echo "# There are ${#USRLIP[@]} user defined lipids: ${USRLIP[@]}"

    N='\'$'\n'

    # Add custom lipid topologies:
    USRFIX=
    for lip in ${USRLIP[@]}
    do
	USRFIX="${USRFIX}"'#include "'$lip'.itp"'"$N"
    done
    if [[ -n $USRFIX ]]
    then
	LSED -i"" -e '/\[ *system *\]/{s/^/'"$USRFIX"'/;}' "$TOP"
    fi

    # Sugar hack
    if true
    then
        # Uncomment the include topology for sugars
	cp $FFDIR/martini_${FFTAG}_sugars.itp ./ || touch martini_${FFTAG}_sugars.itp
	grep -q sugars.itp $TOP && CARBFIX='/sugars.itp/s/^; //' || CARBFIX='/\[ *system *\]/{s/^/#include "martini_'$FFTAG'_sugars.itp"'"$N$N"'/;}'
	LSED -i"" -e "$CARBFIX" "$TOP"
    fi

    if [[ $INSANE =~ "-l " ]]
    then
        # Uncomment the include topology for lipids
	cp $FFDIR/martini_${FFTAG}_lipids.itp ./ || touch martini_${FFTAG}_lipids.itp
	# Check which lipids are defined and which (custom) lipid topologies to add
	lipids=($($SED -n '/^ *\[ *moleculetype/{n;/^ *;/n;p;}' martini_${FFTAG}_lipids.itp))
	for lip in ${alname[@]}
	do
	    isdefined=false
	    for lip2 in ${lipids[@]}; do [[ $lip == $lip2 ]] && isdefined=true && break; done
	    $isdefined && break
	    echo >> martini_${FFTAG}_lipids.itp
	    cat $lip.itp >> martini_${FFTAG}_lipids.itp
	done
	grep -q lipids.itp $TOP && LIPFIX='/lipids.itp/s/^; //' || LIPFIX='/\[ *system *\]/{s/^/#include "martini_'$FFTAG'_lipids.itp"'"$N$N"'/;}'
	LSED -i"" -e "$LIPFIX" "$TOP"
    else
	echo "$INSANE"
    fi

    # When multiscaling, the ions can interact with the protein on the atomistic level
    if $ALL -o [[ -n $MULTI ]]
    then
	IONSITP=$ForceField.ff/ions.itp
    else
	IONSITP=martini_${FFTAG}_ions.itp
	cp $FFDIR/$IONSITP ./ || touch martini_${FFTAG}_ions.itp
    fi

    grep -q ions.itp $TOP && IONFIX='/ions.itp/s/^; //' || IONFIX='/\[ *system *\]/{s,^,#include "'$IONSITP'"'"$N$N"',;}'
    LSED -i -e "$IONFIX" "$TOP"

    # Include user defined topologies
    USRFIX=
    for itp in ${USRITP[@]}
    do
	USRFIX="${USRFIX}"'#include "'$itp'"'"$N"
    done
    USRFIX='/\[ *system *\]/{s,^,'${USRFIX}"$N"',;}'
    LSED -i -e "$USRFIX" "$TOP"

    # We made a topology... extract groups
    NSOL=($(grep '; NDX Solvent' $TOP))
    NSOL=$(( NSOL[4]-NSOL[3]+1 ))
    [[ $NSOL -lt 0 ]] && NSOL=0
    NMEM=($(grep '; NDX Membrane' $TOP))
    NMEM=$(( NMEM[4]-NMEM[3]+1 ))
    [[ $NMEM -lt 0 ]] && NMEM=0
    NPROT=($(grep '; NDX Solute' $TOP))
    NPROT=$(( NPROT[4]-NPROT[3]+1 ))
    [[ $NPROT -lt 0 ]] && NPROT=0

    GRO=$OUT
else
    # We did not make a topology, but there should be an index file.
    NPROT=$(LSED -n '/\[ *Solute/,/\[/{/\[/d;p;}' "$NDX" | wc -w)
    NMEM=$(LSED -n '/\[ *Membrane/,/\[/{/\[/d;p;}' "$NDX" | wc -w)
    NSOL=$(LSED -n '/\[ *Solvent/,/\[/{/\[/d;p;}' "$NDX" | wc -w)
fi

[[ $NMEM -gt 0 ]] && MEMBRANE=true || MEMBRANE=false

NTOT=$(LSED -n '2{p;q;}' "$GRO")
# Assuming that ions come after solvent
NION=$(( NTOT - $(LSED -n $(sed_solvent '{p;d;}') "$GRO" | wc -l) - NMEM - NPROT ))
NSOL=$(( NSOL - NION ))
echo $GRO: NTOT=$NTOT NPROT=$NPROT NSOL=$NSOL NMEM=$NMEM NION=$NION


if [[ $STEP = $NOW ]]
then
    if [[ -n $PDB ]]
    then	
	# Add membrane and solvent to CG list in $base-mart.ndx
	# The structure has Protein/DNA/RNA Membrane Solvent Ions
	cp $base-mart.ndx $base-cg.ndx
	echo -e "$(printf "$fmt\n" `SEQ $((NPROT + 1)) $((NPROT + NMEM + NSOL))` | $SED 's/ 0//g')" >> $base-cg.ndx

	# Add the ions to the CG group of we do not want hybrid ions
	if ! $HybridIons && [[ $NION -gt 0 ]]
	then
	    echo -e "$(printf "$fmt\n" `SEQ $((NPROT+NMEM+NSOL+1)) $((NPROT+NMEM+NSOL+NION))` | $SED 's/ 0//g')" >> $base-cg.ndx
	fi

	echo -e "[ Solute ]\n$(printf "$fmt\n" `SEQ 1 $NPROT` | $SED 's/ 0//g')" >> $base-cg.ndx
    else
	# Everything is in the CG group
	echo -e "[ CG ]\n$(printf "$fmt\n" `SEQ 1 $NTOT`      | $SED 's/ 0//g')"  > $base-cg.ndx
        # And we have an empty Solute group
        echo '[ Solute ]' >> $base-cg.ndx
    fi

    # The membrane tag is always added to the index file; it may be an empty group
    echo '[ Membrane ]' >> $base-cg.ndx
    if [[ $NMEM -gt 0 ]]
    then
	echo -e "$(printf "$fmt\n" `SEQ $((NPROT+1)) $((NPROT+NMEM))` | $SED 's/ 0//g')" >> $base-cg.ndx
    fi

    # Solvent
    echo -e "[ Solvent ]\n$(printf "$fmt\n" `SEQ $((NPROT+NMEM+1)) $((NPROT+NMEM+NSOL+NION))` | $SED 's/ 0//g')" >> $base-cg.ndx

    if $HybridIons && [[ $NION -gt 0 ]]
    then
	# Make ions group (hybrid character) for multiscaling
	echo -e "[ Ions ]\n$(printf "$fmt\n" `SEQ $((NPROT+NMEM+NSOL+1)) $((NPROT+NMEM+NSOL+NION))` | $SED 's/ 0//g')" >> $base-cg.ndx
    fi

    # The whole system... easy.
    echo -e "[ System ] \n$(printf "$fmt\n" `SEQ 1 $NTOT` | $SED 's/ 0//g')" >> $base-cg.ndx

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
    $SED -i '/\[ *molecules *\]/,$s/\(NA\|K\|CA\|ZN\|CU\|CL\)2*[+-]/\1 /' $base-cg.top
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
    cp $FFDIR/martini_${FFTAG}_{ions,lipids}.itp .
    grps=($($SED '/^ *[^[]/d;s/\[ *\(.*\) *\]/\1/;/Solute/d;/System/d;' $DAFT))
    __mdp_cg__energygrps=$($SED 's/ /,/g' <<< ${grps[@]})
fi

trash $base-EM.{tpr,edr,trr} em-out.mdp

MDP=em-sol.mdp
OPT=(em $MDPMS)
OUT=$base-EM.gro
LOG=02-EM.log

mdp_options ${OPT[@]} > $MDP
MD="MDRUNNER -f $MDP -c $GRO -p $TOP -o $OUT -n $NDX -np 1 -l $LOG -force $FORCE $TABLES"
echo $MD

[[ $STEP ==   $NOW     ]] && $EXEC $MD && : $((STEP++)) && GRO=$OUT && archive
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
    MD="MDRUNNER -f $MDP -c $GRO -p $TOP -o $OUT -n $NDX -np $PRNP -l $LOG -force $FORCE $TABLES $MONALL"
    $EXEC $MD
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
    $VirtualSites && DT=(0.0005 0.001 0.002 0.003 0.004) || DT=(0.0005 0.001 0.002)
else
    # Coarsegrained
    # "With DAFT runs, we just skip this step (PR-NpT)"
    [[ -n $DAFT ]] && DT=() || DT=(0.005 0.010 0.020)
fi


# Base MDP options for all cycles
$MEMBRANE && OPT=(cg mem $MDPMS equil npt $RotationalConstraints) || OPT=(cg $MDPMS equil npt $RotationalConstraints)
# Set pressure coupling to PR if running with GMX5
if [[ $GMXVERSION -le 4 ]]
then
    $MEMBRANE && __mdp_equil__tau_p=$__mdp_equil__tau_p,$__mdp_equil__tau_p
fi
run=1
for __mdp_equil__dt in ${DT[@]}
do
    LOG=04-PR-NPT-$__mdp_equil__dt-$run.log
    OUT=$base-PR-NPT-$__mdp_equil__dt-$run.gro
    MDP=pr-npt-$__mdp_equil__dt-$run.mdp
    mdp_options ${OPT[@]} > $MDP
    MD="MDRUNNER -f $MDP -c $GRO -p $TOP -o $OUT -n $NDX -np $NP -l $LOG -force $FORCE $TABLES $MONALL"
    $EXEC $MD
    GRO=$OUT
    : $((run++))
    trash $base-PR-NPT-$__mdp_equil__dt.{cpt,tpr}
done

[[ $STEP ==   $NOW     ]] && : $((STEP++)) && archive
[[ $STOP == $((NOW++)) ]] && exit_clean


#--------------------------------------------------------------------
SHOUT "---STEP 5: UNRESTRAINED NpT"
#--------------------------------------------------------------------

$MEMBRANE && OPT=(cg mem $MDPMS equil $RotationalConstraints usr) || OPT=(cg $MDPMS equil $RotationalConstraints usr)
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
    $VirtualSites && DT=(0.004 $DELT) || DT=(0.002 $DELT)
else
    # Coarsegrained
    DT=(0.020 $DELT)
fi

run=1
for __mdp_equil__dt in ${DT[@]}
do
    LOG=05-NPT-$__mdp_equil__dt-$run.log

    OUT=$base-NPT-$__mdp_equil__dt-$run.gro
    MDP=md-init-$__mdp_equil__dt-$run.mdp
    mdp_options ${OPT[@]} > $MDP
    nsteps=$($SED -n '/nsteps/s/.*=//p' "$MDP")
    before=$(date +%s)
    echo ==--- $GRO $OUT
    MD="MDRUNNER -f $MDP -c $GRO -p $TOP -o $OUT -n $NDX -np $NP -l $LOG -force $FORCE $TABLES $MONALL"
    $EXEC $MD
    timing=$(( $(date +%s) - before ))
    echo "Ran $nsteps steps in $timing seconds"
    GRO=$OUT
    : $((run++))
    trash $base-NPT-$__mdp_equil__dt-$run.{cpt,tpr}
done


# : $((STEP++))

[[ $STEP ==   $NOW     ]] && : $((STEP++)) && archive
[[ $STOP == $((NOW++)) ]] && exit_clean


#--------------------------------------------------------------------
SHOUT "---STEP 6: PRODUCTION RUN"
#--------------------------------------------------------------------

$MEMBRANE && OPT=(cg mem $MDPMS $RotationalConstraints usr) || OPT=(cg $MDPMS $RotationalConstraints usr)
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
    #__mdp_mem__tau_p=12.0,12.0
    __mdp_mem__tau_p=12.0
fi

echo ${STEPS[$STEP]} ${STEPS[$STOP]} ${STEPS[$NOW]}

# If the Stop-Step is set to TPR only write a tpr and do not run
[[ $STEP ==   $NOW     ]] && : $((STEP++))
[[ $STOP == $((NOW++)) ]] && OUT=$base-MD.tpr

echo @@@ $OUT

mdp_options ${OPT[@]} > $MDP

mdsteps=$($SED -n '/nsteps/s/.*=//p' "$MDP")
estimate=$(( (timing*(1000*mdsteps)/nsteps)/1000 ))
estfin=$(( $(date +%s) + estimate ))
echo "Expected runtime for $mdsteps step: $(( estimate/3600 ))H:$(( (estimate%3600)/60 ))M:$(( estimate%60 ))S (until $(date -r $estfin))"

MD="MDRUNNER -f $MDP -c $GRO -p $TOP -o $OUT -n $NDX -np $NP -l $LOG -split -force $FORCE $TABLES -monitor"
$EXEC $MD

# : $((STEP++))

[[ $STEP ==   $NOW     ]] && : $((STEP++)) && archive

#--------------------------------------------------------------------
SHOUT "---DONE SIMULATING"
#--------------------------------------------------------------------

# We stop here if this was the last bit to do
[[ $STOP == $((NOW++)) ]] && exit_clean

# Otherwise we just continue with some more fun stuff
# If we have stuff to do, actually.

#-----------------------------------------------------------------------------------
SHOUT "Huh!? You made it all the way to the end of the script. This can't be right."
#-----------------------------------------------------------------------------------

archive
exit_clean
