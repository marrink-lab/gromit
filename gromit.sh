#!/bin/bash

PROGRAM=${0##*/}
VERSION=2.5      # 20170706.1200
GMXVERSION=5.x.y # But 4.x.y and 2016 are supported
HISTORY="\
"

AUTHORS="Tsjerk A. Wassenaar, PhD"
YEAR="2017"
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

__NOTES_FOR_USE_AND_EDITING__



DESCRIPTION=$(cat << __DESCRIPTION__

$PROGRAM $VERSION is a versatile wrapper for setting up and running
molecular dynamics simulations of proteins and/or nucleic acids in solvent.
The script contains a complete and flexible workflow, consisting of the 
following steps:

    1.   Generate topology from input structure 
         A. Process structure against force field (TOPOLOGY)
         B. Add ligands                           (LIGANDS)
    2.   Set up periodic boundary conditions      (BOX)
    3.   Energy minimize system in vacuum         (EMVACUUM)
    4.   Solvation and adding ions                (SOLVATION)
    5.   Energy minimization                      (EMSOL)
    6.   Position restrained NVT equilibration    (NVT-PR)
    7.   Unrestrained NpT equilibration           (NPT)
    8.   Equilibration under run conditions       (PREPRODUCTION)
    9.   Production simulation                    
         A. Run input file                        (TPR)
         B. Simulation (possibly in parts)        (PRODUCTION)

The program allows running only part of the workflow by specifying the
start and end step (-step/-stop), using an argument uniquely matching 
one of the tags given between parentheses.

This program requires a working installation of Gromacs. To link 
the program to the correct version of Gromacs, it should be placed in the 
Gromacs binaries directory or the Gromacs GMXRC file should be passed as 
argument to the option -gmxrc

The workflow contained within this program corresponds to a standard protocol
that should suffice for routine molecular dynamics simulations of proteins 
and/or nucleic acids in aqueous solution. It follows the steps that are 
commonly taken in MD tutorials (e.g. http://md.chem.rug.nl/~mdcourse/molmod2012/).

This program is designed to enable high-throughput processing of molecular
dynamics simulations in which specific settings are varied systematically. These
settings include protein/nucleic acid, ligand, temperature, and pressure, as well
as many others.

The program supports a number of specific protocols, including Linear Interaction 
Energy (LIE), Transition Path Sampling (TPS), and GRID-distributed MD.


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
#---EXTERNAL/USER STUFF
#--------------------------------------------------------------------

# Some variables can be set here to provide defaults 
# for example for GRID-specific settings.
LFC_ROOT=/grid/enmr.eu
export LFC_HOME=$LFC_ROOT/gromacs
export LFC_HOST=lfcserver.cnaf.infn.it
GRID_SE=se-enmr.cerm.unifi.it

# Structure databases
RCSB=(http://files.rcsb.org/download/ pdb.gz)
BIO=(http://www.rcsb.org/pdb/files/ pdb1.gz)
OPM=(http://opm.phar.umich.edu/pdb/ pdb)

# Automated topology builder
ATB=http://atb.uq.edu.au
ATBSEARCHPDB="$ATB/search.py?"
ATBSEARCHMOL="$ATB/search.py?"
ATBGET="$ATB/download.py?outputType=top&ffVersion=54a7"
ATBITP="$ATBGET&file=rtp_uniatom&molid"
ATBPDB="$ATBGET&file=pdb_uniatom_optimised&molid"
ATBPDBU="$ATBGET&file=pdb_uniatom_unoptimised&molid"
ATBGROMOS=https://atb.uq.edu.au/forcefield_files/atb_gromacs/5/gromos54a7_atb.ff.tar.gz


#--------------------------------------------------------------------
#---Parsing COMMAND LINE ARGUMENTS AND DEPENDENCIES--
#--------------------------------------------------------------------

# Directory where this script is
SDIR=$( [[ $0 != ${0%/*} ]] && cd ${0%/*}; pwd )

# These will be looked for before running, and can be set from the cmdline, e.g.:
#    -gmxrc /usr/local/gromacs-5.1/bin/GMXRC
# If not set, the default name will be searched for in
#    1. the environment (if PROGEVAR is given)
#    2. the directory where this calling script (gromit) is located
#    3. the PATH 
DEPENDENCIES=( gmxrc squeeze)
PROGEXEC=(     GMXRC squeeze)
PROGEVAR=(     GMXRC)


# Run control and files
DIR="."       # Directory to run and write           
fnIN=         # Input file name
TPR=          # Run input file... skip to production or to STEP
NAME=         # Run name
TOP=          # Topology file
FETCH=        # Try to fetch PDB file from web
MSGFILE=/dev/stdout      # Master log file (stdout)
ERRFILE=/dev/stderr      # Error log file (stderr)
ARCHIVE=      # Archive file name
FORCE=false   # Overwrite existing run data
EXEC=         # Execute run
NP=1          # Number of processors
MDP=          # User-provided MDP file
MDARGS=       # User-provided arguments to mdrun
GMXRC=        # File for loading GMX version
SQUEEZE=      # Squeeze executable
ANALYSIS=()   # Analysis tags 
MAXH=-1       # Maximum duration of run
JUNK=()       # Garbage collection
SCRATCH=      # Scratch directory
KEEP=false    # Keep intermediate rubbish (except junk)
HETATM=true   # Keep HETATM records in input PDB file
GRID=false    # Whether this is a grid-enabled run
CMD="$0 $@"   # The command as it was issued
echo $CMD > CMD


# Run control
MONALL=       # Monitor all steps
CONTROL=
CHECKTIME=300 # Run control every five minutes


# Stepping stuff
STEPS=(TOPOLOGY LIGANDS BOX EMVACUUM SOLVATION EMSOLVENT NVT-PR NPT PREPRODUCTION TPR PRODUCTION ANALYSIS END)
get_step_fun() { for ((i=0; i<${#STEPS[@]}; i++)) do [[ ${STEPS[$i]} =~ ^$1 ]] && echo $i; done; }
STEP=
STOP=PRODUCTION

# Macros to echo stuff only if the step is to be executed -- more about steps further down
RED='\x1b[1;31m'
YEL='\x1b[1;33m'
GRN='\x1b[1;32m'
CYA='\x1b[1;36m'
OFF='\x1b[0m'
LINES() { sed -e 's/^/ /;s/$/ /;h;s/./-/g;p;x;p;x;' <<< "$@"; }
SHOUT() { [[ $STEP == $NOW ]] && echo && LINES "$@" | sed 's/^/#/'; }
WARN()  { [[ $STEP == $NOW ]] && echo -e "$RED" && LINES "$@" | sed 's/^/# WARNING: /' && echo -e "$OFF"; }
FATAL() { [[ $STEP == $NOW ]] && echo -e "$RED" && LINES "$@" | sed 's/^/# FATAL: /' && echo -e "$OFF"; exit 1; }
NOTE()  { [[ $STEP == $NOW ]] && echo -e "$YEL" && LINES "$@" | sed 's/^/# NOTE: /' && echo -e "$OFF"; }
LINE()  { [[ $STEP == $NOW ]] && echo && sed -e 's/^/ /;s/$/ /p;s/./-/g;' <<< "$@" | sed 's/^/#/'; }
ECHO()  { [[ $STEP == $NOW ]] && echo "$@"; }


# Force field
ForceFieldFamilies=(gromos  charmm  amber   opls )
ForceFieldSolvents=(spc     tip3p   tip3p   tip4p)
SolventFiles=(      spc216  spc216  spc216  tip4p)
ForceField=gromos45a3
WaterModel=
SolModel=
SolName=SOL
SolFile=
LIGANDS=()
VirtualSites=false
AtomTypes=()
MoleculeTypes=()


# System setup
PBCDIST=2.25             # Minimal distance between periodic images
BOXTYPE=dodecahedron     # Box type
Salt=NA,CL               # Salt species
SaltCharge=1,-1          # Charges of salt species
Salinity=0.1539976       # Salt concentration of solvent
CHARGE=                  # Net charge to set on system
NDLP=false               # Use optimal simulation cell (near-densest-lattice-packing)


# Simulation parameters
TIME=0                   # Production run time (ns)
EquilTime=0.1            # Equilibration run time (ns)          
PreTime=0.5              # Preproduction run time (ns)
EMSteps=500              # Number of steps for EM
AT=0.05                  # Output frequency for positions, energy and log (ns)
Temperature=200,300      # Degree Kelvin
Tau_T=0.1                # ps
Pressure=1               # Bar
Tau_P=1.0                # ps
PosreFC=200,200          # Position restraint force constant(s)
Electrostatics=          # Electrostatics scheme to use if not opting for default
RotationalConstraints=   # Use rotational constraints, which is mandatory with NDLP
SEED=$$                  # Random seed for velocity generation


# GROUP DEFINITIONS
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


# User defined gromacs program options and simulation parameters (way flexible!)
PROGOPTS=()              # User-defined program options (--program-option=value)
MDPOPTS=()               # User-defined mdp parametesrs (--mdp-option=value)


hlevel=1
olevel=1
# Function for parsing help from option list
help_fun()
{
  local desc='/^.*#=[0-'${hlevel}']/s//    /p'
  local item='/#==[0-'${olevel}']/{s/).*#==./  /p;d;}'
  sed -n '/##>> OPTIONS/,/##<< OPTIONS/{'"${item};${desc};}" $SDIR/$PROGRAM
}


# Function for displaying USAGE information
function USAGE ()
{
    cat << __USAGE__

$PROGRAM version $VERSION:

$DESCRIPTION

(c)$YEAR $AUTHORS
$AFFILIATION

__USAGE__

help_fun
}


if [[ -z "$1" ]]; then
  echo "No command line arguments give. Please read the program usage:"
  USAGE
  exit
fi


BAD_OPTION ()
{
    echo
    echo "Unknown option "$1" found on command-line"
    echo "It may be a good idea to read the usage:"

    USAGE

    echo " /\                                               /\ " 
    echo "/||\  Unknown option "$1" found on command-line  /||\ "
    echo " ||   It may be a good idea to read the usage     || "
    echo " ||                                               || "

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
error_function() { a="$@"; errors_array+=(${x// /QQQ}); FATAL "$@"; }
warnings_array=()
warning_function() { a=$@; warnings_array+=(${x// /QQQ}); WARN "$@"; }
notes_array=()
note_function() { a=$@; notes_array+=(${x// /QQQ}); NOTE "$@"; }


##>> OPTIONS
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

    case $1 in
	#=0
	#=0 OPTIONS
        #=0 =======
        #=0
	-h       ) USAGE                                ; exit 0 ;; #==0 Display help
	--help   ) USAGE                                ; exit 0 ;; #==1 Display help
	-hlevel  ) hlevel=$2                            ; shift 2; continue ;; #==1 Set level of help (use before -h/--help)
	-olevel  ) olevel=$2                            ; shift 2; continue ;; #==1 Set level of options to display

	#=1
        #=1 File options
	#=1 ------------
        #=1
	-f       ) fnIN=$2                              ; shift 2; continue ;; #==0 Input coordinate file (PDB)
        -g       ) MSGFILE=$2                           ; shift 2; continue ;; #==1 Standard output log file (default: /dev/stdout)
        -e       ) ERRFILE=$2                           ; shift 2; continue ;; #==1 Standard error log file (default: /dev/stderr)
        -tpr     ) TPR=$2                               ; shift 2; continue ;; #==1 Run input file 
	-name    ) NAME=$2                              ; shift 2; continue ;; #==1 Name of project
	-top     ) TOP=$2                               ; shift 2; continue ;; #==1 Input topology file
	-atp     ) AtomTypes+=($2)                      ; shift 2; continue ;; #==2 Additional atom type definitions (force field file)
	-itp     ) MoleculeTypes+=($2)                  ; shift 2; continue ;; #==2 Additional molecule type definitions 
	-l       ) LIGANDS+=($2)                        ; shift 2; continue ;; #==2 Ligands to include (topology or structure,topology)
	-mdp     ) MDP=$2                               ; shift 2; continue ;; #==2 MDP (simulation parameter) file
	-scratch ) SCRATCH=$2                           ; shift 2; continue ;; #==2 Scratch directory to perform simulation
	-fetch   ) FETCH=$2                             ; shift 2; continue ;; #==1 Database to fetch input structure from 
        -rmhet   ) HETATM=false                         ; shift  ; continue ;; #==2 Whether or not to remove HETATM records

	#=1
	#=1 Overall control options
	#=1 -----------------------
	#=1
	-step    ) STEP=$2                              ; shift 2; continue ;; #==1 Step to start protocol
	-stop    ) STOP=$2                              ; shift 2; continue ;; #==1 Step to end protocol
        -grid    ) GRID=true                            ; shift 1; continue ;; #==2 GRID-enabled run
        -keep    ) KEEP=true                            ; shift  ; continue ;; #==2 Whether or not to keep intermediate data
	-dir     ) DIR=$2                               ; shift 2; continue ;; #==2 Directory where to perform simulation (make if required)
	-np      ) NP=$2                                ; shift 2; continue ;; #==1 Number of processors/threads to use
        -maxh    ) MAXH=$2                              ; shift 2; continue ;; #==2 Maximum time to run (in hours)
	-archive ) ARCHIVE=${2%.tgz}.tgz                ; shift 2; continue ;; #==2 Archive file name to save data in
	-force   ) FORCE=true                           ; shift  ; continue ;; #==2 Whether or not to force redoing parts already run
        -noexec  ) EXEC=echo                            ; shift  ; continue ;; #==2 Whether or not to actually execute the commands

	#=1
	#=1 Simulation control options
	#=1 --------------------------
	#=1
	-rtc     ) RotationalConstraints=rtc            ; shift  ; continue ;; #==2 Whether or not to use rotational constraints
	-ndlp    ) NDLP=true; RotationalConstraints=rtc ; shift  ; continue ;; #==2 Whether or not to use NDLP (molecular shaped) PBC
        -bt      ) BOXTYPE=$2                           ; shift 2; continue ;; #==2 Box type to use
	-salt    ) Salt=$2                              ; shift 2; continue ;; #==2 Salt to use (NA,CL)
	-conc    ) Salinity=$2                          ; shift 2; continue ;; #==2 Salt concentration
	-sq      ) SaltCharge=$2                        ; shift 2; continue ;; #==2 Charge of ions from salt (1,-1)
	-charge  ) CHARGE=$2                            ; shift 2; continue ;; #==2 
	-t       ) Temperature=$2                       ; shift 2; continue ;; #==1 Temperature
	-ttau    ) Tau_T=$2                             ; shift 2; continue ;; #==2 Temperature coupling constant
        -p       ) Pressure=$2                          ; shift 2; continue ;; #==1 Pressure
	-ptau    ) Tau_P=$2                             ; shift 2; continue ;; #==2 Pressure coupling constant
	-d       ) PBCDIST=$2                           ; shift 2; continue ;; #==1 Distance between images over PBC
	-prfc    ) PosreFC=$2                           ; shift 2; continue ;; #==2 Position restraint force constant
	-time    ) TIME=$2                              ; shift 2; continue ;; #==1 Time to run production simulation (ns)
	-at      ) AT=$2                                ; shift 2; continue ;; #==1 Output resolution
        -em      ) EMSteps=$2                           ; shift 2; continue ;; #==2 Number of steps in EM
        -equil   ) EquilTime=$2                         ; shift 2; continue ;; #==2 Equilibration run time
        -pre     ) PreTime=$2                           ; shift 2; continue ;; #==2 Time for preproduction run
	-elec    ) Electrostatics=$2                    ; shift 2; continue ;; #==2 Electrostatics treatment (cutoff/RF/PME)
	-ff      ) ForceField=$2                        ; shift 2; continue ;; #==1 Force field to use
        -vsite   ) VirtualSites=true                    ; shift  ; continue ;; #==2 Whether or not to use virtual sites
	-seed    ) SEED=$2                              ; shift 2; continue ;; #==2 Seed for random number generator
	-solvent ) SolModel=$2                          ; shift 2; continue ;; #==2 Solvent model name(s); can be itp file(s)
	-solfile ) SolFile=$2                           ; shift 2; continue ;; #==2 Solvent (configuration) file

	#=1
	#=1 Analysis options
	#=1 ----------------
	#=1
	-lie     ) LIE=true                             ; shift  ; continue ;; #==2 Whether or not to use LIE setup and analysis 
        -analysis) ANALYSIS+=($2)                       ; shift 2; continue ;; #==2 Analysis protocols to run

        #=2
        #=2 Perform analysis specified (provided it is implemented) *STR:   None
        #=2 Analysis routines can be added to the end of the script
        #=2 They should be tagged to allow being called through the 
        #=2 command line. 
        #=2 
        #=2 Currently available routines:
        #=2
        #=2     * LIE
        #=2       ---
        #=2       Extract ligand - environment interaction energies
        #=2       for post-hoc LIE calculations. This analysis 
        #=2       requires running with the option -lie and is then
        #=2       selected automatically.

	#=1
	#=1 Monitor options
	#=1 ---------------
	#=1
	-monall  ) MONALL=-monitor                      ; shift 1; continue ;; #==2 Monitor all steps using control script
	-control ) CONTROL=$2                           ; shift 2; continue ;; #==2 Simulation monitor script
	-ctime   ) CHECKTIME=$2                         ; shift 2; continue ;; #==2 Time for running monitor

	#=2
	#=2 A control process is either a program, script or command 
        #=2 that monitors the production run and terminates it
        #=2 upon a certain condition, indicated by a zero exit code.
        #=2 The control process gets the following command-line
        #=2 options added from mdrun:
	#=2
        #=2         -s ${NAME}-MD.tpr 
        #=2         -f ${NAME}-MD.part#.trr 
        #=2         -x ${NAME}-MD.part#.xtc 
        #=2         -e ${NAME}-MD.part#.edr
        #=2         -g ${NAME}-MD.part#.log 

	#=1
	#=1 Advanced control options
	#=1 ------------------------
	#=1
	#=2 This program allows specifying options for advanced control of 
	#=2 program invocation and simulation parameters. These options are 
	#=2 described below.
	#=2

	#=2    --mdp-option=value         Command-line specified simulation parameters
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

	#=2    --program-option=value     Command-line specified program parameters
	--*      ) PROGOPTS+=($1)                       ; shift  ; continue ;;
	#=2
	#=2 This will add "-option value" to the command line of the call to 'program'.
	#=2 Note that this does not allow overriding options specified internally. 
	#=2 Trying to do so will result in an error due to double specification of the 
	#=2 option. If the option takes multiple arguments, then 'value' should be a 
	#=2 comma separated list.


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


# Time. To keep track of the remaining run time
START=$(date +%s)


# Set the script directory
SCRIPTDIR=$(cd ${0%${0##*/}}; pwd)


# Directory where command was issued
BDIR=$(pwd)


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

# Set the correct sed version for multi-platform use
# Also try to avoid GNU specific sed statements for the
# poor bastards that are stuck with one of those Mac things
SED=$(which gsed || which sed)


# Awk expression for extracting moleculetype
#    - at the line matching 'moleculetype' 
#      read in the next line
#      continue reading next lines until one is not beginning with ;
#      print the first field
AWK_MOLTYPE='/moleculetype/{getline; while ($0 ~ /^ *;/) getline; print $1}'


#--------------------------------------------------------------------
#---GROMACS STUFF
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
    [[ -f $progr ]] && echo $progr && return 0
    
    # Check if the program is in the PATH
    which $progr 2>/dev/null && return 0 || return 1
}


##  1. GROMACS  ##

# Check and set the gromacs related stuff  
GMXRC=$(find_program_function gmxrc)
echo Gromacs RC file: $GMXRC

# Source the gromacs RC file if one was found
# Otherwise the script will rely on the active gromacs commands 
[[ $GMXRC ]] && source $GMXRC 

# Find out which Gromacs version this is
GMXVERSION=$(mdrun -h 2>&1 | sed -n '/^.*VERSION \([^ ]*\).*$/{s//\1/p;q;}')
# From version 5.0.x on, the commands are gathered in one 'gmx' program
# The original commands are aliased, but there is no guarantee they will always remain
[[ -z $GMXVERSION ]] && GMXVERSION=$(gmx -h 2>&1 | sed -n '/^.*VERSION \([^ ]*\).*$/{s//\1/p;q;}')
# Version 2016 uses lower case "version", which is potentially ambiguous, so match carefully
[[ -z $GMXVERSION ]] && GMXVERSION=$(gmx -h 2>&1 | sed -n '/^GROMACS:.*gmx, version \([^ ]*\).*$/{s//\1/p;q;}')
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

[[ $GMXVERSION -gt 4 ]] && 

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
    echo Gromacs data directory: $GMXLIB
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

# 2. SQUEEZE/NDLP for minimal-volume simulation.
# - Requires Gromacs with RTC support
# - Requires SQUEEZE executable:
if $NDLP
then
    SQUEEZE=$(find_program_function squeeze)
    if [[ $? != 0 ]]
    then
        FATAL "NDLP SETUP REQUESTED, BUT SQUEEZE EXECUTABLE NOT FOUND"	
    fi
    if ! $SQUEEZE -h >/dev/null 2>&1
    then
        echo 
	echo "Squeeze was probably compiled for a different version of Gromacs."
        FATAL "NDLP SETUP REQUESTED, BUT SQUEEZE EXECUTABLE FAILED" 
    fi
fi


METHOD=$(cat << __METHOD__
Simulations were performed using GROMACS version $GMXVERSION [1] through the automated workflow $PROGRAM [2].


[1] $GMXREF
[2] $THISREF

__METHOD__
)

#--------------------------------------------------------------------
#---TIMING (GRID STUFF)
#--------------------------------------------------------------------

# If UNTIL is not set, then we run until all work is done, or until we crash
UNTIL=


# If we have a proxy, then we adhere to it, assuming to be running on the GRID
# The allowed runtime is set slightly lower to allow the script to finish
$GRID && which voms-proxy-info && LEFT=$(( $(voms-proxy-info -timeleft &> /dev/null) - 300 )) || LEFT=


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
#---INITIAL CHECKS AND LOGGING
#--------------------------------------------------------------------


## 1. Expand options that can take multiple, comma-separated values

# This concerns options for equilibration, such as temperature,
# pressure and position restraint force constants. These will
# be used to set up cycles of equilibration. Position restraints Fc
# and temperature are cycled together in STEP 6 (NVT), followed by 
# pressure equilibration in STEP 7 (NPT).

# Store the Internal Field Separator (IFS)
ifs=$IFS
IFS=","
Temperature=($Temperature)
Tau_T=($Tau_T)
Pressure=($Pressure)
Tau_P=($Tau_P)
PosreFC=($PosreFC) 
Salt=($Salt)
SaltCharge=($SaltCharge)
# Restore the field separator
IFS=$ifs


## 2. Echo options for direct modulation of program calls

# These options are formatted like --program-option=value
# This will add '-option value' on the command line of the 
# specific program. If an option takes multiple arguments
# then these should be given comma separated. The commas 
# will be replaced by spaces.
# *NOTE*: Some options are defined internally and can not be 
# overriden. Attempting to do so will result in an error, as
# options will be doubly specified.

if [[ -n $PROGOPTS ]]
then
    echo "# Program options specified on command line:"
    for ((i=0; i<${#PROGOPTS[@]}; i++)); do echo "#     ${PROGOPTS[$i]}"; done
    echo "#==="
fi


## 3. Echo mdp options specified on the command line

# These options are formatted like --mdp-param=value
# This will add 'param = value' to the MDP file for 
# all runs following energy minimization. 
# *NOTE*: Options specified on the command line take
# precedence over internal parameters and over those
# read from an mdp file, provided as value to option
# -mdp 

if [[ -n $MDPOPTS ]]
then
    echo "# Simulation parameters specified on command line (note how flexible!):"
    for ((i=0; i<${#MDPOPTS[@]}; i++)); do echo "#     ${MDPOPTS[$i]}"; done
    echo "#==="
fi


#--------------------------------------------------------------------
#---WARMING UP VARIABLE GYMNASTICS
#--------------------------------------------------------------------


# Parse input file names - expand to full path
pdb=${fnIN##*/}                                  # Filename
dirn=$(cd ${fnIN%${fnIN##*/}}./; pwd)            # Directory
ext=${pdb##*.}                                   # Extension
[[ "$ext" == "tpr" && -z $TPR ]] && TPR=$pdb     # Got a run input file as input file
[[ -n $fnIN ]] && fnIN=$dirn/$pdb                # Input file with full path
topdir=$(cd ${TOP%${TOP##*/}}./; pwd)            # Full path to directory of topology file
[[ -n $TOP ]]  && TOP=$topdir/${TOP##*/}         # Topology file with full path


# Set the name. If one is given, that one is used. Otherwise, 
# if an input structure is given, the name is derived from it. 
# If there is no input structure, but there is a ligand, then
# "ligand" is used. If nothing is set, then just use
# "wallace" 
base=${pdb%.*}                                   # Set basename
[[ -n $NAME ]] && base=$NAME                     # Override base name if name is given
[[ -z $base && -n $LIGANDS ]] && base=ligand     # Still unset.., set to "ligand" if we have a ligand
[[ -z $base ]] && base=wallace                   # Okay, then just go with this


# Change working directory, creating one if necessary
[[ ! -d $DIR ]] && mkdir -p $DIR
[[ -n $TPR && ! -f $DIR/${TPR##*/} ]] && cp $fnIN $DIR
cd $DIR
DIR=$(pwd)


echo
echo "# $( date +"%a %b %d %H:%M:%S %Y" ) "
echo "# Command:"
echo $CMD
echo

[[ -n $fnIN    ]] && echo "# Input structure:     $fnIN"
[[ -n $TOP     ]] && echo "# Input topology:      $TOP"
                     echo "# Base name:           $base"
                     echo "# Source directory:    $dirn"
                     echo "# Command issued from: $BDIR"
                     echo "# Run directory:       $DIR"
[[ -n $SCRATCH ]] && echo "# Scratch directory:   $SCRATCH"

echo
echo "# GROMACS "
echo "#========="
echo "# - binaries:  $GMXBIN"
echo "# - data:      $GMXDATA"
echo "# - libraries: $GMXLIB"
echo 


# Copy topology stuff if we specify a topology
[[ -n $TOP ]] && cp $topdir/*itp ./


# Set a trap for signals
function archive ()
{
    if [[ -n $ARCHIVE ]]
    then
	tar cfz $DIR/$ARCHIVE.tmp.tgz $(ls | grep -v $ARCHIVE)
	mv $DIR/$ARCHIVE.tmp.tgz $DIR/$ARCHIVE
    fi
    exit
}
trap "archive" 2 9 15


# Set the forcefield tag
case $ForceField in
    gromos*) ForceFieldFamily=gromos;;
    amber*)  ForceFieldFamily=amber;;
    charmm*) ForceFieldFamily=charmm;;
    opls*)   ForceFieldFamily=opls;;
esac


#--------------------------------------------------------------------
#---SET THE SOLVENT MODEL
#--------------------------------------------------------------------

# If it is not specified, the solvent is the water model 
# most appropriate for the force field. This can be overriden
# by explicitly stating an alternative water model or by
# providing a solvent model topology. 

# The solvent model needs to be accompanied by a solvent
# structure file, which is to be used for the solvation.
# For the water models these files are available in the
# distribution.

if [[ -z $SolModel ]]
then
    # Get the right solvent model for the force field selected
    for ((i=0; i<${#ForceFieldFamilies[@]}; i++))
    do
	if [[ $ForceFieldFamily == ${ForceFieldFamilies[$i]} ]]
	then
	    WaterModel=${ForceFieldSolvents[$i]}
	    SolModel=$WaterModel
	    SolventTopology=$ForceField.ff/$WaterModel.itp
	    [[ -z $SolFile ]] && SolFile=${SolventFiles[$i]}
	fi
    done
else
    # If the name is known (in the ForceFieldSolvents list)
    # use the corresponding file
    found=false
    # Get the right solvent model for the force field selected
    for ((i=0; i<${#ForceFieldFamilies[@]}; i++))
    do
	if [[ $SolModel == ${ForceFieldSolvents[$i]} ]]
	then
	    WaterModel=$SolModel
	    SolventTopology=$WaterModel.itp
	    [[ -z $SolFile ]] && SolFile=${SolventFiles[$i]}
	    found=true
	fi
    done

    if ! $found
    then
	# Solvent model was specified, but not found
	if [[ -f $SolModel ]]
	then
	    SolName=$(awk "$AWK_MOLTYPE" $SolModel)
	    SolventTopology=$SolModel
	    NOTE Using solvent model $SolName from $SolModel with file $SolFile
	else
	    # Don't know what to do...
	    FATAL The solvent model should be a registered name or provided as topology
	fi
    fi
fi


#--------------------------------------------------------------------
#---STEPPING STUFF
#--------------------------------------------------------------------

# Set the starting/stopping step
# For LIE calculations also execute the analysis part
$LIE && ANALYSIS+=(LIE)

# If we have a TPR as input file run only production
[[ -n $TPR ]] && STEP=PRODUCTION

# If there is analysis to be done, set the stopping step accordingly
[[ -n $ANALYSIS ]] && [[ PRODUCTION == ${STOP}* || END == ${STOP}* ]] && STOP=ANALYSIS

# Step up to the step-in step
for ((i=0; i<${#STEPS[@]}; i++)); do [[ ${STEPS[$i]} == ${STEP}* ]] && STEP=$i && break; done
# Step up to the stop-step: stop if the step stepped up to is the step to stop at
for ((i=0; i<${#STEPS[@]}; i++)); do [[ ${STEPS[$i]} == ${STOP}* ]] && STOP=$i && break; done

echo "# Starting at step $STEP:    ${STEPS[$STEP]}"
echo "# Stopping at step $STOP:    ${STEPS[$STOP]}"

# Remove flags from previous runs
[[ -e DONE  ]] && rm DONE
[[ -e ERROR ]] && echo "# Found ERROR flag, probably from previous run. Trying again." && rm ERROR


#--------------------------------------------------------------------
#---INFORMATIVE OUTPUT--          
#--------------------------------------------------------------------


echo "# Starting MD protocol for $fnIN"

if [[ -z $TPR ]]
then
    echo "# Using $ForceFieldFamily force field $ForceField with $WaterModel water model"

    [[ -n $Electrostatics ]] \
	&& echo "# Using $Electrostatics for treatment of long range coulomb interactions" \
	|| echo "# Inferring electrostatics treatment from force field family (check mdp files)"

    $VirtualSites \
	&& echo "# Using virtual sites" \
	|| echo "# Not using virtual sites"
    
    $NDLP \
	&& echo "# Simulations will be performed using a near-densest lattice packing unit cell" \
	|| echo "# Simulations will be performed in a rhombic dodecahedron unit cell" 
fi


#--------------------------------------------------------------------
#---SIMULATION PARAMETERS--          
#--------------------------------------------------------------------

## OT N ## For every parameter not defined the default is used
## NOTE ## This is probably fine for equilibration, but check the defaults to be sure
## E OT ## The list as is was set up for gromacs 4.5 and 5.1

# This function lists the mdp options requested based on a preceding tag
function mdp_options ()
{    
    for tag in $@
    do
	# Find variables declared with specified tag
        for opt in `set | grep ^__mdp_${tag}__`
        do
	    # Strip everything from the first = to the end 
	    # to get the variable name
            var=${opt%%=*}
	    # Strip everything up to the first = to get
	    # the value
            val=${opt#*=}
	    # Replace the tag and redeclare in local space
	    # If the variable was already declared it will
	    # be overridden.
            local ${var/mdp_$tag/mdp}=$val
        done
    done
    # Find all variables starting with __mdp__ and echo them
    set | $SED -n '/^__mdp__/{s/__mdp__//;s/,/ /g;p;}'
}


#--------------------------------------------------------------------
# Global parameters

__mdp_md__dt=0.002
[[ $GMXVERSION -gt 4 ]] && __mdp_md__nstlist=10 || __mdp_md__nstlist=5
[[ $GMXVERSION -gt 4 ]] && __mdp_md__cutoff_scheme=Verlet

# Virtual site specific
if $VirtualSites; then
  __mdp_md__dt=0.004
  __mdp_md__lincsorder=5
  __mdp_md__nstlist=2
fi

# Output parameters
TIME=$(python -c "print int(1000*$TIME/$__mdp_md__dt + 0.5 )") 
AT=$(python -c "print int(1000*$AT/$__mdp_md__dt + 0.5)") 
__mdp_md__nsteps=$TIME
__mdp_md__nstxout=$AT
__mdp_md__nstvout=0 
__mdp_md__nstfout=0
__mdp_md__nstlog=$AT
__mdp_md__nstenergy=$AT
[[ $GMXVERSION -gt 4 ]] && __mdp_md__nstxout_compressed=$AT || __mdp_md__nstxtcout=$AT

# Coupling
# Listed here are the overall controls. Specific controls
# (temperature,pressure,time constants) are controlled through
# the command line interface and taken care of in steps 6/7
# The temperature is set later on
__mdp_md__tcoupl=v-rescale
__mdp_md__nsttcouple=$__mdp_md__nstlist
__mdp_md__nstpcouple=$__mdp_md__nstlist
__mdp_md__compressibility=4.5e-5

# Nonbonded interactions
__mdp_md__coulombtype=PME
__mdp_md__fourierspacing=0.125
__mdp_md__rcoulomb=0.9
__mdp_md__rlist=0.9
__mdp_md__rvdw=1.4
__mdp_md__pme_order=4

# Other
__mdp_md__constraints=all-bonds
__mdp_md__comm_mode=Linear
__mdp_md__comm_grps=System
__mdp_md__nstcomm=$__mdp__nstlist


#--------------------------------------------------------------------
# Rotational constraints

__mdp_rtc__comm_mode=RTC
__mdp_rtc__comm_grps=Solute

#--------------------------------------------------------------------
# Force field specific parameters

# AMBER, all versions 
# (parameters taken from Acpype protein-ligand tutorial by Alan Wilter ...)
__mdp_amber__coulombtype=PME
__mdp_amber__coulomb_modifier=Potential-shift-Verlet
__mdp_amber__rcoulomb=1
__mdp_amber__rlist=1
__mdp_amber__vdwtype=cut-off
__mdp_amber__rvdw=1.0
__mdp_amber__DispCorr=EnerPres
__mdp_amber__lincs_iter=2
#__mdp_amber__optimize_fft=yes
# Additional options specified in named tutorial:
# __mdp_amber__pcoupltype=Parinello-Rahman
# __mdp_amber__tau_p=0.5

# GROMOS96, all versions
__mdp_gromos__coulombtype=Reaction-Field
__mdp_gromos__rcoulomb=1.4
__mdp_gromos__epsilon_rf=54

# CHARMM
__mdp_charmm__coulombtype=Switch
__mdp_charmm__rcoulomb_switch=1.0
__mdp_charmm__rcoulomb=1.2
__mdp_charmm__vdwtype=Switch
__mdp_charmm__rvdw_switch=1.0
__mdp_charmm__rvdw=1.2

# OPLS/AA
__mdp_opls__coulombtype=PME
__mdp_opls__rcoulomb=1.0
__mdp_opls__rlist=1.0
__mdp_opls__rvdw=1.0
__mdp_opls__fourierspacing=0.135
__mdp_opls__dt=0.001

#--------------------------------------------------------------------
# Energy minimization

__mdp_em__define=-DFLEXIBLE
__mdp_em__integrator=steep
__mdp_em__nsteps=$EMSteps
__mdp_em__pbc=xyz
__mdp_em__constraints=none

#--------------------------------------------------------------------
# Equilibration runs: position restraints, NVT, NPT

# Position restraints are relieved at step 7
TIME=$(python -c "print int(1000*$EquilTime/$__mdp_md__dt + 0.5 )")
__mdp_equil__define=-DPOSRES
__mdp_equil__nsteps=$TIME
__mdp_equil__nstlog=10
__mdp_equil__nstenergy=10
__mdp_equil__nstxtcout=0

# Velocities are only generated once
# After the first NVT/PR cycle 'genvel' is set to no
# and the other two options are ignored
__mdp_equil__genvel=yes
__mdp_equil__gen_seed=$SEED
__mdp_equil__gen_temp=${Temperature[0]}

#--------------------------------------------------------------------
# User specified parameters
# These will be used in PR-NVT, PR-NPT, MD-INIT, MD-PRE and MD
# To specify parameters for a single step, make cunning use of
# -step and -stop.

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
    USER_MDP=( $( $SED '/^ *$/d;/^ *;/d;s/\s*=\s*/=/;s/\s*$//;s/ /,/g;' $MDP ) )
    
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


function writeToLog() 
{  
    # Write message to log as formatted string in the same way as the server cron scripts do
  
    local MSG_STATUS=$2
    local LOGMSG="# $( date +"%a %b %d %H:%M:%S %Y" ) MDS $$ $MSG_STATUS: $1"
    echo "$LOGMSG"
}


function print_messages() 
{
    [[ -z $3 ]] && return 0
    N=$#
    echo -e $1 There were $((N-2)) $2:
    for ((i=2; i<=$N; i++))
    do 
        echo [$i] ${!i//QQQ/ }
    done
}


function uploadToSE() 
{
    # Upload results to a storage element
   
    # These variables should be available on the command line or
    # at least higher up in the script in a user-modifiable section
    LFC_ROOT='/grid/enmr.eu/'
    export LFC_HOME=$LFC_ROOT/gromacs
    export LFC_HOST='lfcserver.cnaf.infn.it'
  
    if [[ -e $ARCHIVE ]]
    then
        lcg-cr -v --vo enmr.eu -d se-enmr.cerm.unifi.it -l $ARCHIVE $ARCHIVE
    else
        echo "ERROR: no archive $ARCHIVE found. Cannot upload to SE"
    fi
}


function archive ()
{
  if [[ -n $ARCHIVE ]]; then
    
    FILES=$( ls )
    EXCLUDE="$ARCHIVE gmx45setup.sh"
    for EX in $EXCLUDE; do
      FILES=$( echo $FILES | $SED "s/ $EX / /g;" ) 
    done
    
    tar cfz gmx-server-tmp-archive.tgz $FILES
    [[ -e gmx-server-tmp-archive.tgz ]] && \mv gmx-server-tmp-archive.tgz $ARCHIVE
    
  fi
}


function exit_clean()
{
    # Finished... Removing flag for running
    [[ -f RUNNING ]] && rm -f RUNNING

    # Set flag indicating job done
    touch DONE

    # Add message to log file
    if [[ -n $1 ]]
    then
	writeToLog "$@"
    fi
    writeToLog "Wrapping up and exiting."

    # Clean up
    echo "# Deleting redundant files:"
    printf "%-25s %-25s %-25s %-25s %-25s\n" ${JUNK[@]} \#*\#
    for i in ${JUNK[@]} \#*\# *.bck; do [[ -f $i ]] && rm $i; done

    [[ -n $SCRATCH ]] && cd $DIR # && rm -rf $SCRATCH

    # Archive
    archive
    #uoploadToSE

    writeToLog "$PROGRAM finished successfully"

    exit 0
}


function exit_error()
{
  touch ERROR
 
  # Give error message
  writeToLog "$PROGRAM STEP $STEP ${STEPS[$STEP]} : $@" ERROR

  # Set flags
  [[ -f RUNNING ]] && rm -f RUNNING
 
  [[ -n $SCRATCH ]] && cp * $DIR && cd $DIR # && rm -rf $SCRATCH

  # Archive everything
  archive

  # Upload stuff to SE
  # uploadToSE

  [[ -n $MSGFILE ]] && exec 1>&3
  [[ -n $ERRFILE ]] && exec 2>&4

  print_messages $CYA NOTES    ${notes_array[@]}
  print_messages $YEL WARNINGS ${warnings_array[@]}
  print_messages $RED ERRORS   ${errors_array[@]}

  FATAL "$@"
}


function terminate()
{
    echo
    writeToLog "SIGNAL CAUGHT... "
    while [[ -n $1 ]]
    do
	if kill -0 $1 >/dev/null 2>&1 
	then
	    writeToLog "Terminating \"$(ps -o command= $1)\" (PID $1)"
	    kill -15 $1
	else
	    writeToLog "Attempt to terminate already terminated process \"$(ps -o command= $1)\" (PID $1)"
	fi
	shift
    done
    
    writeToLog "Exiting"
    
    touch ERROR
    touch TERMINATED

    exit 1
}


function pdb2gmx_error()
{
  sed -n -e '/^Fatal error/,/^----/p' $1
  exit_error "Converting structure with pdb2gmx failed."
}


# Trashing files
function trash()
{
    for item in $@; do JUNK[${#JUNK[@]}]=$item; done
}


# Routine for gathering program specific options
# In the version edited by MvD (UU) this routine 
# uses numbers: for ((opt=0; opt<${#...
# Do not understand why.
function program_options()
{    
    local OPTS=
    for opt in ${PROGOPTS[@]}
    do
	if [[ $opt =~ --$1- ]]
	then
            OPTS="$OPTS $($SED 's/--[^-]*//;s/=/ /;' <<< $opt)"
	fi
    done
    echo $OPTS
}


# Check for the existence of all arguments as files
function all_exist()
{
    for f in ${OUTPUT[@]} 
    do 
        [[ -e $f || -e $dirn/$f ]] || return 1
    done
}


# This is a function to generate a sequence of numbers
# On Linux platforms there usually is 'seq', but we 
# can not depend on it.
function SEQ(){ for ((i=$1; i<=$2; i++)); do echo $i; done };


# Routine for extracting the charge from a TPR file
# a. Extract and set the active molecule name 
AWK_TPR_MOLNAME='/^ *name=/{sub(".*=","",$0); T=$0}'
# b. Extract the moleculetype name
AWK_TPR_MOLTYPE='/^ *moltype *=/{M=$4}'
# c. Extract the number of instances of the current moleculetype
#    and initialize charge for moleculetype  
AWK_TPR_MOLNUM='/#molecules *=/{N[M]=$3; C[M]=0}'
# d. Extract and accumulate charge
AWK_TPR_CHARGE='/^ *atom.*q=/{sub(".*q=","",$0); sub(",.*","",$0); C[T]+=$0}'
# e. Finalize
AWK_TPR_END='END{S=0; for (i in C) {S+=N[i]*C[i]}; if (S<0) S-=0.5; else S+=0.5; printf "%d\n", S}'
function getCharge ()
{
    gmxdump -s $1 2>&1 | awk "$AWK_TPR_MOLNAME $AWK_TPR_MOLTYPE $AWK_TPR_MOLNUM $AWK_TPR_CHARGE $AWK_TPR_END"
}


# Routine for generating a simple index file
function INDEX()
{
    [[ -n $2 ]] && fn=$2 || fn=basic.ndx
 
    exec 6>&1 && exec >$fn

    fmt="%5d %5d %5d %5d %5d %5d %5d %5d %5d %5d"

    # Total number of atoms
    N=$(awk '{getline; print; exit}' $1)
    echo "[ System ]"
    printf "$fmt\n" `SEQ 1 $N` | $SED 's/ 0//g'
  
    # Solvent atoms (including ions, etc, listed after 'SOL')
    SOL=$(( $(sed -n '/'$SolName'/{=;q;}' $1) - 2 ))
    echo "[ Solvent ]"
    printf "$fmt\n" `SEQ $SOL $N` | $SED 's/ 0//g'

    # Base system: solute and membrane, if present
    echo "[ Base ]"
    printf "$fmt\n" `SEQ 1 $((SOL - 1))` | $SED 's/ 0//g'

    # Membrane, if any
    MEMBRANE=$($SED -n '/\(POP\|DPP\|DMP\|DOP\|PPC\)/{=;q;}' $1)
    if [[ -n $MEMBRANE ]]
    then
        echo '[ Membrane ]'
        printf "$fmt\n" `SEQ $MEMBRANE $((SOL - 1))` | $SED 's/ 0//g'
    else
        MEMBRANE=SOL
    fi

    echo '[ Solute ]'
    printf "$fmt\n" `SEQ 1 $((MEMBRANE - 1))` | $SED 's/ 0//g'

    exec 1>&6 6>&-

    return 0
}


function MDRUNNER ()
{
    writeToLog "MDRUNNER"    
    
    local NP=1
    local fnOUT=
    local fnNDX=
    local FRC=false
    local SPLIT=
    local TABLES=
    local MONITOR=false
    while test -n "$1"; do
	case $1 in
	    -f)       local fnMDP=$2        ; shift 2; continue;;
            -c)       local  fnIN=$2        ; shift 2; continue;;
	    -n)       local fnNDX=$2        ; shift 2; continue;;
            -o)       local fnOUT=$2        ; shift 2; continue;;
	    -p)       local fnTOP=$2        ; shift 2; continue;;
	    -l)       local fnLOG=$2        ; shift 2; continue;;
	    -force)   local   FRC=$2        ; shift 2; continue;;
	    -np)      local    NP=$2        ; shift 2; continue;;
	    -split)   local SPLIT=-noappend ; shift  ; continue;;
	    -monitor) local MONITOR=true    ; shift  ; continue;;
	    -tables) local TABLES="-table table.xvg -tablep tablep.xvg"; shift; continue;;
            *)  echo "PANIC!: Internal Argument Error ($1) in routine MDRUNNER"; exit;;
        esac
    done
    
    
    if [[ "${fnIN##*.}" == "tpr" ]]
    then
	local TPR=$fnIN
	local fnOUT=${TPR%.tpr}.gro
    else
        # Check input
	[[ -f $fnIN  ]] || exit_error "Input structure not found in routine MDRUNNER ($fnIN)"
	[[ -f $fnTOP ]] || exit_error "Input topology not found in routine MDRUNNER ($fnTOP)"
	[[ -f $fnMDP ]] || exit_error "Input parameter file not found in routine MDRUNNER ($fnMDP)"
	[[ -n $fnNDX && -f $fnNDX ]] || INDEX $fnIN $fnNDX
	local TPR=
    fi
    
    
    # Infer basename from output file name
    baseOUT=${fnOUT%.*}
    
    
    if [[ -n $SCRATCH ]]
    then
        # Make sure to copy the TPR file if we already have one
	[[ -f $DIR/$baseOUT.tpr ]] && cp $DIR/$baseOUT.tpr .    
	
        # Make sure we deal well with checkpointing
	[[ -f $DIR/$baseOUT.cpt ]] && cp $DIR/$baseOUT.cpt .
    fi
    
    
    if $FRC
    then
	removed=()
	for z in ${baseOUT}.*; do [[ -f "$z" && "$z" != "$TPR" ]] && rm $z && removed[${#removed[@]}]=$z; done
	echo "# Forced execution."
        [[ -n $removed ]] && writeToLog "Removed files:" && echo ${removed[@]}
    fi
    
    
    if [[ -f $fnOUT || -f $DIR/$fnOUT ]]
    then
	writeToLog "Output found ($fnOUT). Skipping step ${STEPS[$STEP]}"
	return 0
    fi
    
    
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
	local -i step=($($SED  -n -e '/^ *Step *Time *Lambda/{n;h;}' -e '${x;p;}' $last))
	
	# Check the number of steps for this cycle
	local -i RUNSTEPS=$(sed -n '/nsteps/s/^.*=\([^;]*\).*$/\1/p' $fnMDP)
	
	if [[ $step -gt 0 ]]
	then
	    writeToLog "A log file exists which reports having run $step of $RUNSTEPS steps ($last)"
	fi	
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
    [[ ! -e $fnOUT && -e $baseOUT.cpt ]] && writeToLog "Found partial ${STEPS[$STEP]} results... Continuing"
    
    
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
    
    
    # Set up the mdrun command and start it in the background 
    MDRUN="${GMX}mdrun -nice 0 -deffnm $baseOUT -c $fnOUT -cpi $baseOUT.cpt -nt $NP $SPLIT $(program_options mdrun) -maxh $MAXH"
    echo
    echo "$MDRUN" | tee -a $fnLOG
    echo
    $MDRUN >>$fnLOG 2>&1 &
    local -i MDPID=$!
    writeToLog "MDRUN PID: $MDPID"
    
    
    # Start the monitor (if this is a production run)
    if $MONITOR && [[ -n $CONTROL ]]
    then
	MON=$CONTROL
	MON=${MON/@PID/$MDPID}
	MON=${MON/@TPR/$baseOUT.tpr}
	MON=${MON/@TRR/$baseOUT.trr}
	MON=${MON/@XTC/$baseOUT.xtc}
	MON=${MON/@EDR/$baseOUT.edr}
	MON=${MON/@LOG/$baseOUT.log}
	MON=${MON/@TOP/$fnTOP}
	MON=${MON/@NDX/$fnNDX}
	MON=${MON/@MDP/$fnMDP}
	writeToLog "Monitoring run"
	echo
	echo "$MON"
	echo
	# Mark the output
	# Make sure that the exit code from the monitor is exitcode of the process. 
	($MON | sed 's/^/MONITOR: /'; exit ${PIPESTATUS[0]}) &

	local -i MONID=$!
	writeToLog "Monitor PID: $MONID"
	
	# In most cases, the monitor finishes before the run.
	# The exit code then tells what to do (see below).
	trap "terminate $MDPID $MONID" SIGHUP SIGINT SIGTERM SIGCHLD

    	wait $MONID
	MONEXIT=$?
	echo MONITOR EXIT: $MONEXIT
    else
	trap "terminate $MDPID" SIGHUP SIGINT SIGTERM SIGCHLD
    fi


    # Wait for the MD run to finish
    wait $MDPID
    MDEXIT=$?


    echo "__MDRUN__" >>$LOG


    case $MONEXIT in
	0) writeToLog "Monitor exited after run finished without taking action (PASS)";;
	1) exit_error "Monitor exited after run finished: condition not met, while it should be (FAIL)";;
	2) writeToLog "Monitor terminated process, because condition was met (PASS)";;
	3) exit_error "Monitor terminated process, because condition was met, but shouldn't be (FAIL)";;
	4) writeToLog "Monitor killed process, because condition was met (PASS)";;
	5) exit_error "Monitor killed process, because condition was met, but shouldn't be (FAIL)";;
    esac


    if [[ $MDEXIT != 0 ]]
    then
	exit_error "MDRUN exited with non-zero exit code ($MDEXIT) ... Exiting $PROGRAM."
    fi


    # If we split then we have to do some more work to see if we have finished
    if [[ -n $SPLIT ]]
    then
	step=($SED -n -e '/Step *Time *Lambda/{n;h;}' -e '${x;p;}' $fnLOG)
	[[ $step == $RUNSTEPS ]] && cp ${fnLOG%.log}.gro $fnOUT
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


# Macros to echo stuff only if the step is to be executed
function SHOUT() { [[ $STEP == $NOW ]] && echo && L=$($SED -e 's/./-/g;s/^/#/;' <<< " $@ ") && echo $L && sed 's/^/#/' <<< "$@ " && echo $L ; }
function LINE()  { [[ $STEP == $NOW ]] && echo && $SED -e 's/^/ /;s/$/ /p;s/./-/g;' <<< "$@" | sed 's/^/#/'; }
function ECHO()  { [[ $STEP == $NOW ]] && echo "$@"; }


# Macro to do stuff only if the step is to be executed
function DO() { [[ $STEP == $NOW ]] && echo "$@" && $@; }


# Sed wrapper (loud sed) - echo the command line before running it
function LSED() { echo $SED "$@" 1>&2; $SED "$@"; }


# Always ECHO the first line
NOW=$STEP


#--------------------------------------------------------------------
SHOUT "---= THIS IS WHERE WE really START =---"
#--------------------------------------------------------------------

NOW=0

# We may be working in two directories:
#     1. The RUN directory $DIR
#     2. The SCRATCH directory $SCRATCH
# If we are working in the scratch directory, we
# copy back files after every step to keep the 
# run directory in sync.
# What we copy back depends on the value of $KEEP.
# If that is 'true' everything is copied.
# We always go to the scratch directory if we have one.
[[ -n $SCRATCH ]] && cd $SCRATCH && echo $DIR > SOURCE


# If we do not have an input file, definitely skip the first step
# This may happen if we set up a ligand-in-solvent simulation, or
# just a box of solvent, or maybe a membrane...
[[ $STEP == $NOW && -z $fnIN ]] && : $((STEP++)) && echo "# Skipping step 1A: No structure to generate topology for."


# If we are at this step, but have a top file already, increase the STEP
if [[ $STEP == $NOW && -n $TOP && -n $fnIN ]]
then
    echo "# Topology file ${TOP##*/} provided for structure ${fnIN##*/}. Skipping step 1A."
    # If the accompanying GRO file does not exist, convert the PDB file
    [[ -e $base.gro ]] || ${GMX}editconf -f $fnIN -o $base.gro &>/dev/null
    : $(( STEP++ ))
fi


# If we are still at this step, do some checking of input and 
# maybe fetch a file from the PDB.
if [[ $STEP == $NOW ]]
then

    # Check the input file
    if [[ ! -f $dirn/$pdb ]]
    then
	[[ -f $dirn/$pdb.pdb    ]] && pdb=$pdb.pdb
	[[ -f $dirn/$pdb.pdb.gz ]] && pdb=$pdb.pdb.gz
	[[ -f $dirn/$pdb.gz     ]] && pdb=$pdb.gz
    fi
    PDB=$dirn/$pdb

    if [[ -f $PDB ]]
    then
        [[ -n $SCRATCH ]] && cp $PDB .
    elif [[ -n $FETCH ]]
    then
        # Try fetching it from the repository
        # Biological units use MODEL/ENDMDL statements to separate chains
        # These need to be removed.
        case $FETCH in
            pdb  | PDB)  
                GET=${RCSB[0]}${pdb%%.*}.${RCSB[1]}
                RMMOD=false 
                ;;
            bio* | BIO*)
                # Check if a number is given for selecting the BIO assembly
                [[ $FETCH =~ : ]] && BIO[1]=pdb${FETCH##*:} 
                GET=${BIO[0]}${pdb%%.*}.${BIO[1]}
                RMMOD=true  
                ;;
            opm  | OPM)  
                GET=${OPM[0]}${pdb%%.*}.${OPM[1]}
                RMMOD=false
                ;;
        esac

	echo "# Input file not found, but will try fetching it ($GET)"
        # Use wget or curl; one of them should work
	wget $GET 2>/dev/null || curl -O $GET
        # Remove the URL part of what we just got
        pdb=${GET##*/}
        [[ $GET =~ \.gz$ ]] && gunzip $pdb && pdb=${pdb%.gz} 
        # With biological units, Gromacs chokes on the pdb1 extension
        [[ $pdb =~ \.pdb1$ ]] && mv $pdb ${pdb%1} && pdb=${pdb%1}
        PDB=$pdb

        if $RMMOD
        then
            # Removing the MODEL records and replacing ENDMDL by TER
            sed -i '' -e /^MODEL/d -e s/ENDMDL/TER/ $PDB
        fi

	[[ -n $SCRATCH ]] && cp $pdb $DIR
    else
	echo "# Input file $dirn/$pdb specified, but not found."
	echo "# For automated download of a PDB file, add -fetch to the command line."
	exit 1
    fi


    # Allow feeding zipped files - temporarily unzipping them
    if [[ "${PDB: -2}" == "gz" ]]
    then
	gz=$PDB
	pdb=$base.pdb
	[[ -z $NAME ]] && base=${base%.*}
	gunzip -c $gz > $pdb && trash $pdb
	ls
    fi

    # Remove HETATM records
    if ! $HETATM
    then
        sed -i'.het' /^HETATM/d $pdb
    fi
fi


#--------------------------------------------------------------------
#---CHECK AND FIX THE INPUT STRUCTURE
#--------------------------------------------------------------------

## I. Protein

# a. Sequence
# b. Side chains
# c. PTMs

## II. Ligands and cofactors 

# a. Ions
# b. Ligands/Factors

## III. Nucleic acids

# a. DNA
# b. RNA

## IV. Lipids

## V. Solvent


if [[ $STEP == $NOW ]]
then
  :
fi


#---------------------------------------------------------------------
SHOUT "---STEP 1A: GENERATE STRUCTURE AND TOPOLOGY FOR INPUT PDB FILE"
#---------------------------------------------------------------------

# Output for this section:
OUTPUT=($base.top $base.gro)


# Delete existing output if we force this step
# If we are in $SCRATCH, we just refrain from copying
[[ $STEP == $NOW ]] && $FORCE && rm ${OUTPUT[@]}


# If we are here now, we should generate a top file
[[ $STEP == $NOW ]] && TOP=$base.top


## I. Build command 


# 1. Basic stuff
PDB2GMX="${GMX}pdb2gmx -v -f $dirn/$pdb -o $base.gro -p $base.top -ignh -ff $ForceField -water $WaterModel"


# 2. Position restraints
#    * The position restraint fc (-posrefc) is bogus and 
#      intended to allow easy replacement with sed.
#      These will be placed under control of a #define
PDB2GMX="$PDB2GMX -i $base-posre.itp -posrefc 999"


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
    
    if [ -e "$dirn/pdb2gmx.query" ]; then
	GMXQUERY=$(cat $dirn/pdb2gmx.query)
	PDB2GMX="$PDB2GMX -ter -lys -asp -glu -his"
    else
	GMXQUERY=
    fi
fi


## II. Process (or not)


# Skipping because of existing output (and not forcing)
if [[ $STEP == $NOW ]] && ! $FORCE && $(all_exist ${OUTPUT[@]})
then
    echo "# Output found, skipping topology generation" 
    : $((STEP++))
fi


# Skipping verbosely; showing what would have been run
if [[ $STEP == $NOW && -n $EXEC ]]
then
    echo -e "# Skipping:\n@ $PDB2GMX" 
    : $((STEP++))
fi


# Execute
if [[ $STEP == $NOW ]]
then
    LOG=01-PDB2GMX.log

    # The following lines echo the command to the log 
    # and capture the results as a here-document
    echo "echo $GMXQUERY | $PDB2GMX" | tee -a $LOG
    echo ": << __PDB2GMX__" >>$LOG
    ERROR="Generation of topology with pdb2gmx failed. Probably the input structure has issues. Check $LOG for more information."
    echo $GMXQUERY | $PDB2GMX >>$LOG 2>&1 # || exit_error $ERROR
    EXITCODE=$?
    echo "__PDB2GMX__" >>$LOG
    [[ $EXITCODE == 0 ]] || pdb2gmx_error $LOG

    NATOMS=$(awk '{getline; print; exit}' $base.gro)

    for itp in ${base}-posre*.itp
    do
	# Use sed to insert an #ifndef .. #endif block
	# and to subsitute 999 for POSREFC
	# We introduce newlines with sed by processing the raw string
	# with \n using the bash $ operator, which substitutes 
	# backslash escape characters with control characters.
	# This is necessary for introducing newlines with BSD sed (OS X).
        $SED -i.bck -e $'1s/^/#ifndef POSREFC\\\n  #define POSREFC 200\\\n#endif\\\n/' \
                    -e 's/999  999  999/POSREFC POSREFC POSREFC/' $itp
    done

    #-- Simplifying topology if identical chains are found

    # a. First get the list of moleculetype itp files
    ITP=($(ls ${base}_*.itp 2>/dev/null) ${base}.top)

    # b. Extract the moleculetype names from the itp files
    MTP=($(for i in ${ITP[@]}; do awk "$AWK_MOLTYPE" $i; done))

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
		echo "# Removing duplicate moleculetype definition in ${ITP[$i]}"
		# 1. remove the #include statement for that itp file
		# 2. rename the moleculetype under [ system ]
		LSED -i.bck -e "/${ITP[$i]}/d" -e "/[\s*system\s*]/,\$s/${MTP[$i]}/${MTP[$j]}/" $base.top
		# List the file for removal
		trash ${ITP[$i]} $($SED 's/_/-posre_/' <<< ${ITP[$i]})
		
		break
	    fi
	done
    done	

    # Save what is needed to the RUN directory
    [[ -n $SCRATCH ]] && cp $LOG $base.gro $base.top $base*.itp $DIR

    echo "# $(date): Finished topology"
fi


## BOOKKEEPING ##

# Bookkeeping is always done, also if the step is not executed #

# Index groups, Coupling Groups, Energy Groups
if [[ -e $DIR/$base.gro ]]
then
    NATOMS=$(awk '{getline; print; exit}' $DIR/$base.gro)

    # Set the solute start and end atom:
    Biomol=(1 $NATOMS)

    # Set the solute start and end atom:
    Solute=(1 $NATOMS)

    # If we add ligands, the solute is also in the ligand environment
    Ligenv=(1 $NATOMS)

    # List the solute as coupling group and as energy group
    CoupleGroups=(Solute)
    EnergyGroups=(Solute)
fi


# Set the topology. Best use the local local one :p
# Both will be good if we are running in $DIR
# The bottom one will be best if we run in $SCRATCH
[[ -e $DIR/$base.top ]] && TOP=$DIR/$base.top
[[ -e $base.top      ]] && TOP=$base.top


# If we have a gro file here, it is named $base.gro.
# We won't have one if we did not have an input file,
# like when setting up a box of solvent
# Just maybe the gro file is not exactly here.
if [[ -n $fnIN ]]
then
    GRO=$base.gro
    [[ -e $GRO ]] || GRO=$DIR/$GRO
else
    GRO=
fi


# End of step
[[ $STOP ==   $NOW     ]] && exit_clean
[[ $STEP == $((NOW++)) ]] && : $((STEP++))


# If there are no ligands to add, increase the step one more
[[ -z $LIGANDS && $STEP == $NOW ]] && : $((STEP++))


#--------------------------------------------------------------------
SHOUT "---STEP 1B: ADD LIGANDS SPECIFIED ON COMMAND LINE"
#--------------------------------------------------------------------

# ADD LIGANDS

LOG=01-LIGANDS.log


# Skip this step if there are no ligands to add
[[ $STEP == $NOW && ${#LIGANDS[@]} == 0 ]] && : $((STEP++))


# Output for this section:
OUTPUT=($base-lig.top $base-lig.gro)


# Delete existing output if we force this step
[[ $STEP == $NOW ]] && $FORCE && rm ${OUTPUT[@]}


# Check output
if [[ $STEP == $NOW ]] && ! $FORCE && $(all_exist ${OUTPUT[@]})
then
    echo Found output... Skipping addition of ligand\(s\).
    : $((STEP++))
fi


# Execute this step
if [[ $STEP == $NOW && ${#LIGANDS[@]} -gt 0 ]]
then
    echo Adding Ligands


    ## 1. Input GRO file

    # This step requires a GRO file, but we allow stepping in with a
    # PDB file, as long as the TOP file is present. The PDB file is
    # then converted to GRO format here.
    [[ -z $GRO && -f $fnIN ]] && ${GMX}editconf -f $fnIN -o $base.gro >/dev/null 2>&1

    # If $GRO is not set here, it must be $base.gro
    # This can also be written GRO=${GRO:-$base.gro}
    [[ -z $GRO ]] && GRO=$base.gro

    # If $GRO is not available, it may be in $DIR
    [[ ! -e $GRO && -e $DIR/$GRO ]] && cp $DIR/$GRO ./


    ## 2. Setting up coordinate file containing ligands

    #     Overwrite the file if it already exists.
    #     a. Retain title
    #        Here we use a trick to write the title of the new .gro file.
    #        If we already have $GRO the first line is copied,
    #        otherwise a new line is written. To make sure that whichever
    #        goes to standard out and into the new file, the commands are 
    #        grouped together, using curly braces. Mind the semicolon.
    #        The redirect captures the stdout of the compound command.
    #        An empty line is added to be replaced by the number of atoms.
    { [[ -f $GRO ]] && head -1 $GRO || echo Ligand; } > $base-lig.gro
    echo >> $base-lig.gro

    #     b. Note the current number of atoms
    #        Again, base the results on the existence of $base.gro
    [[ -f $GRO ]] && natoms=$(awk '{getline; print; exit}' $GRO) || natoms=0

    #     c. Store the line number of the box definition
    box=$((natoms+3)) 

    #     d. Write the coordinates to the new file, if there are any
    [[ $natoms -gt 0 ]] && $SED -ne "3,$((natoms+2))p" $GRO >> $base-lig.gro


    ## 3. Listing of ligands

    #     a. Moleculetype names listing
    mols=()
    
    #     b. Moleculetype definitions
    itps=()

    #     c. Add ligand coordinates to the structure file
    #
    #        Ligands are specified like:
    #
    #            -l structure[,topology[,name]]
    #
    #        The structure file should contain the coordinates
    #        for a single ligand. The topology file should 
    #        contain a [ moleculetype ] definition of the ligand,
    #        but may be omitted if the definition is already available,
    #        e.g. from the default libraries or from a definition in the
    #        master topology file. 
    #        If no file is specified containing the moleculetype 
    #        definition, then the moleculetype name can be specified.
    #        Otherwise, the name is set to the residuename of the first 
    #        atom. 
    #
    for lig in ${LIGANDS[@]}
    do
	# The ligand specification must contain a coordinate file,
        # but the topology is optional:
        #   -l file.gro[,file.itp]
	# Split the structure file name from the rest.
	struc=${lig%%,*}       # structure file
	rest=${lig#$struc}     # topology file or empty

	# If the structure has an absolute path, leave it.
	# Otherwise set a full path to allow running from $SCRATCH
	[[ ${struc:0:1} == "/" ]] || struc=$BDIR/$struc

	# Check the file format and convert to GRO if necessary.
	# If stripping of ".gro" gives the same variable, it will 
	# be a PDB file.
	if [[ ${struc%.gro} == $struc ]]
	then	    
	    s=${struc##*/}.gro
	    [[ ! -e $s ]] && ${GMX}editconf -f $struc -o $s
	    struc=$s
	    trash $s
	fi

	# Get the atom count for the ligand
	latoms=$(awk '{getline; print; exit}' $struc)

	# Add the atoms to the GRO file
	$SED -n 3,$((latoms+2))p $struc >> $base-lig.gro

	# Add the number to the total
	: $((natoms += latoms))

	mtd=
	# Check whether there is more
	if [[ -z $rest ]]
	then
	    # If nothing else is specified, take the name from the
	    # residue name of the first atom in the GRO file
	    # (the characters five to ten on the line).
	    lname=$($SED -n '3{s/^.....\(.....\).*/\1/p;q;}' $struc)
	else
	    # Take off the preceding comma
	    rest=${rest#,}
	    # Check what is given: try to take off whatever matches
	    # up to and including the first comma (may be nothing)
	    lname=${rest#*,}
	    rest=${rest%$lname}
	    # If there is a rest still, then it should be a file
	    # with a moleculetype definition (mtd)
	    # Otherwise lname contains a filename or the molecule name	    
	    if [[ -n $rest ]]
	    then
		mtd=${rest%,}
                [[ ${mtd:0:1} == "/" ]] || mtd=$BDIR/$mtd
	    elif [[ -e $lname || -e $BDIR/$lname ]]
	    then
		mtd=$lname
                [[ ${mtd:0:1} == "/" ]] || mtd=$BDIR/$mtd
		# Extract the first moleculetype name
		lname=($(awk "$AWK_MOLTYPE" $mtd))
	    fi
	fi

	# Add the ligand name to the list
	mols[${#mols[@]}]=$lname

        if [[ -n $mtd ]]
	then
            # If it is not an absolute path, it is relative to $BDIR
	    [[ ${mtd:0:1} == "/" ]] || mtd=$BDIR/$mtd

            # Add the file containing the definition to the list
	    [[ -n $mtd ]] && itps[${#itps[@]}]=$mtd

	    # If the file is not here, it needs to be copied
	    [[ -e ${mtd##*/} ]] || cp $mtd ./
	    
	    # If the file is not there, it needs to be copied
	    [[ -e $DIR/${mtd##*/} ]] || cp $mtd $DIR
	fi
    done

    
    ## 4. Finalize GRO file

    #     a. Add the original box definition
    { [[ -f $GRO ]] && $SED -ne "${box}p" $GRO || echo 0 0 0; } >> $base-lig.gro
    
    #     b. Update the atom count on the second line
    $SED -i.bck "2s/^.*/$(printf %5d $natoms)/" $base-lig.gro


    ## 5. Topology

    #     a. Set up the topology file containing ligands
    #        adding #include statements before [ system ] directive
    if [[ -n $TOP ]]
    then
	# We already have a topology.
        # First copy the topology up to [ system ]:
        #  * Suppress printing (not to print the directive)
        #  * At the directive, quit
        #  * Otherwise, print the line
	$SED -n "/^\[ system \]/q;p;" $TOP > $base-lig.top
    else
	# We don't have a topology yet.
	# Build one!
	echo -e "#include \"$ForceField.ff/forcefield.itp\"\n#include \"$SolventTopology\"" > $base-lig.top
    fi

    #     b. Add include statements for ligands. Only add each file once
    #        Use sed to replace spaces by newlines, and feed the results 
    #        to 'sort -u' to uniqify. Loop over the resulting entries, 
    #        #including each in the topology 
    for mtd in $($SED $'s/ /\\\n/g' <<< ${itps[@]} | sort -u); do echo "#include \"$mtd\""; done >> $base-lig.top

    #     c. Add the original system definition (after adding a newline)
    echo >> $base-lig.top
    if [[ -n $TOP ]]
    then
        # i.  Print from the [ system ] directive up to and including the [ molecules ] directive
	EXPR1='/^\[ *system *\]/,/^\[ *molecules *\]/;' 
        # ii. At the [ molecules ] directive, start a loop, printing non-blank lines
	#     NOTE: MOLECULES ARE ONLY ADDED IF WE ALREADY HAD AN INPUT STRUCTURE
	[[ -f $GRO ]] && EXPR2='/^\[ *molecules *\]/{while (getline) {if ($0 !~ /^ *$/) print}}' || EXPR2=""
	awk "$EXPR1 $EXPR2" $TOP >> $base-lig.top 
    else
	echo -e '[ system ]\nLigand in solvent\n\n[ molecules ]' >> $base-lig.top
    fi

    #     d. Then add the molecule list for the ligands
    #        Add specifications of identical molecules together
    #mols=( $($SED $'s/ /\\\n/g' <<< ${mols[@]} | uniq -c) )
    mols=( $(for i in ${mols[@]}; do echo $i; done | uniq -c) )
    for ((i=1; i<${#mols[@]}; i+=2)); do printf "%-5s %10d \n" ${mols[$i]} ${mols[$((i-1))]} ; done >> $base-lig.top


    [[ -n $SCRATCH ]] && cp $base-lig.gro $base-lig.top $DIR
    echo "# $(date): Added ${#LIGANDS[@]} ligands to coordinate and topology file"
fi


## BOOKKEEPING ##

# Bookkeeping is always done, also if the step is not executed 
# In this case only if we actually had ligands

if [[ ${#LIGANDS[@]} -gt 0  ]] 
then    

    [[ -n $GRO && ! -f $GRO ]] && GRO=$DIR/$GRO

    # The number of atoms of the base system (without ligand)
    # NATOMS corresponds to the first non-base-system atom.
    [[ -f $GRO ]] && NATOMS=$(( $(awk '{getline; print; exit}' $base.gro) + 1 )) || NATOMS=1

    # Since ligands were aded, there must be $base-lig.gro somewhere
    [[ -f $base-lig.gro ]] && LGRO=$base-lig.gro || LGRO=$DIR/$base-lig.gro

    # The number of atoms including ligands
    LATOMS=$(awk '{getline; print; exit}' $LGRO)

    # The groups 
    Solute=(1 $LATOMS)
    Ligand=($NATOMS $LATOMS)

    # Set the correct structure/topology file as master structure/topology
    GRO=$DIR/$base-lig.gro
    TOP=$DIR/$base-lig.top
fi


if [[ -n $TOP ]]
then
    # Add atomtypes and moleculetypes if given
    sed /^#include.*forcefield.itp/q $TOP > $base-usr.top
    for atp in ${AtomTypes[@]}
    do
	echo
	sed -n -e '/\[ *atomtypes *\]/p' -e '/\[ *atomtypes *\]/,/^ *\[/{/^ *\[/d;p;}' $atp
    done >> $base-usr.top
    for mtp in ${MoleculeTypes[@]}
    do
	echo
	sed -n -e '/\[ *moleculetype *\]/,/\[ *system *\]/{/\[ *system *\]/d;p;}' $mtp
    done >> $base-usr.top
    sed -n -e 'H' -e '/^#include.*forcefield.itp/{n;x;}' -e '${x;p;}' $TOP >> $base-usr.top
    TOP=$DIR/$base-usr.top
fi


# End of step
[[ $STOP ==   $NOW     ]] && exit_clean
[[ $STEP == $((NOW++)) ]] && : $((STEP++))


#--------------------------------------------------------------------
#SHOUT "---STEP 1C: (GUESS... :))"
#--------------------------------------------------------------------

# Well... This should be setting up a membrane.
# However, that is far from trivial, and the suggested route
# is using a combined workflow 
#
#    1. gromit    (AA structure+topology of protein) 
#    2. martinate (CG equilibrated structure)
#    3. backward  (Reverse transformation)
#    4. gromit    (AA simulation)
#
# All these tools are available from http://www.cgmartini.nl/
#

#--------------------------------------------------------------------
SHOUT "---STEP 2: SET PERIODIC BOUNDARY CONDITIONS"
#--------------------------------------------------------------------


# Output for this section:
OUTPUT=($base-pbc.gro)


# Log file
LOG=02-PBC.log


# Delete existing output if we force this step
[[ $STEP == $NOW ]] && $FORCE && rm ${OUTPUT[@]}


# If there is no structure file at this point, build a 
# rhombic dodecahedron with diameter $PBCDIST
if [[ $STEP == $NOW && ! -f $GRO ]]
then
    echo "# Building rhombic dodecahedron solvent box with diameter $PBCDIST"
    echo -e "Solvent box\n    0" > $base-pbc.gro    
    python -c "print \"\".join([\"%10.5f\"%($PBCDIST*i) for i in [1,1,.7071068,0,0,0,0,.5,.5]])" >> $base-pbc.gro

    [[ -n $SCRATCH ]] && cp $base-pbc.gro $DIR

    # Skip the rest of this step
    # Also skip the following step (little sense in EM in vacuum if there is nothing in the vacuum)
    : $((STEP+=2))

    # Stop here if stop was set to PBC or EMVAC
    [[ $STEP -ge $STOP ]] && exit_clean
fi


## I. Build command

## 1. Select the program
if $NDLP
then
    PBC="$SQUEEZE -f $GRO -s $GRO -o $base-pbc.gro -d $PBCDIST -c"
    tag=squeeze
else
    PBCDIST=$(python -c "print 0.5*$PBCDIST")
    PBC="${GMX}editconf -f $GRO -o $base-pbc.gro -bt $BOXTYPE -d $PBCDIST -c"
    tag=editconf
fi

## 2. Add program options specified on command line
PBC="$PBC $(program_options $tag)"

## 3. Process
if [[ $STEP == $NOW && -n $EXEC ]]
then
    # Not running, but showing what would have been run
    echo -e "# Skipping: \n$PBC"   
elif [[ $STEP == $NOW ]] && ! $FORCE && $(all_exist ${OUTPUT[@]})
then
    echo "# Output found. Skipping setting up PBC."
elif [[ $STEP == $NOW ]]
then
    # Echo the command and capture the output nicely
    echo "echo 0 0 0 | $PBC" | tee $LOG
    echo ": << __PBC__" >>$LOG
    ERROR="Setting up periodic boundary conditions failed. Check $LOG for more information."
    echo 0 0 0 | $PBC >>$LOG 2>&1 || exit_error $ERROR
    echo "__PBC__" >>$LOG


    # Check the box. The minimum vector length should be 3
    # We only check the first value, which is the length of
    # the first vector. For a rhombic dodecahedron it is the 
    # length of all vectors. For an NDLP unit cell it gives
    # the length of the longest vector. If that is less than
    # 3, the cell is obviously too small. In both cases, the 
    # cell is replaced by a rhombic dodecahedron with length 3.
    box=($(tail -1 $base-pbc.gro))
    if [[ ${box%%.*} -lt 3 ]]
    then
	cat << __WARNING__

WARNING:

The box set up according to the method ($tag) 
and the distance criterion ($PBCDIST) is too small:

u: ${box[0]} ${box[3]} ${box[4]}
v: ${box[1]} ${box[5]} ${box[6]}
w: ${box[2]} ${box[7]} ${box[8]}

Replacing box with a rhombic dodecahedron of radius 3.

__WARNING__

	box='   3.00000   3.00000   2.12132   0.00000   0.00000   0.00000   0.00000   1.50000   1.50000'
        awk -vbox="$box" 'NR==2{X=$1+3} NR==X{print box; exit} 1' $base-pbc.gro > tmp.gro
	mv tmp.gro $base-pbc.gro
	
    fi

    # rm aminoacids.dat

    [[ -n $SCRATCH ]] && cp $LOG $base-pbc.gro $DIR
fi

## BOOKKEEPING ##


# bookkeeping is always done, also if the step is not executed 

GRO=$DIR/$base-pbc.gro


# End of step
[[ $STOP ==   $NOW     ]] && exit_clean
[[ $STEP == $((NOW++)) ]] && : $((STEP++))


#--------------------------------------------------------------------
SHOUT "---STEP 3: RUN EM IN VACUUM"
#--------------------------------------------------------------------

$KEEP || trash $base-EMv.{tpr,edr,log,trr} em-vac-out.mdp

MDP=em-vac.mdp
OPT=(em)
OUT=$base-EMv.gro
LOG=03-EMv.log
STRT=$(date +%s)

# Command for running MD
MD="MDRUNNER -f $MDP -c $GRO -p $TOP -o $OUT -np 1 -l $LOG -force $FORCE"

if [[ $STEP == $NOW ]]
then
    # Generate mdp file
    mdp_options ${OPT[@]} > $MDP

    # Execute command
    $EXEC $MD && : $((STEP++))

    [[ -n $SCRATCH && -e $OUT ]] && cp $LOG $OUT $DIR >/dev/null 2>&1

    # Set name of GRO file to output of EM
    GRO=$DIR/$OUT
fi

# Check for exit
[[ $STOP == $((NOW++)) ]] && exit_clean

LASTRUN=$(( $(date +%s) - STRT ))


#--------------------------------------------------------------------
SHOUT "---STEP 4: SOLVATION AND ADDING IONS"
#--------------------------------------------------------------------


OUTPUT=($base-sol.gro $base-sol.top $base-sol.ndx)


# If we execute this step and we force regeneration of files
# then delete any output files that may be present
[[ $STEP == $NOW ]] && $FORCE && rm ${OUTPUT[@]}


# Check output
if [[ $STEP == $NOW ]] && ! $FORCE && $(all_exist ${OUTPUT[@]})
then
    echo "# Found output... Skipping addition of solvent(s)."
    : $((STEP++))
fi


# If there is no topology file yet, we have to set one up
if [[ -z $TOP ]]
then
    TOP=$base.top
    cat << __TOP__ > $base.top
#include "$ForceField.ff/forcefield.itp

[ system ]
Box of solvent, maybe with ions

[ molecules ]
__TOP__

    # Add atomtypes and moleculetypes if given
    sed /^#include.*forcefield.itp/q $TOP > $base-usr.top
    for atp in ${AtomTypes[@]}
    do
	echo
	sed -n -e '/\[ *atomtypes *\]/p' -e '/\[ *atomtypes *\]/,/^ *\[/{/^ *\[/d;p;}' $atp
    done >> $base-usr.top
    for mtp in ${MoleculeTypes[@]}
    do
	echo
	sed -n -e '/\[ *moleculetype *\]/,/\[ *system *\]/{/\[ *system *\]/d;p;}' $mtp
    done >> $base-usr.top
    sed -n -e 'H' -e '/^#include.*forcefield.itp/{n;x;}' -e '${x;p;}' $TOP >> $base-usr.top
    TOP=$DIR/$base-usr.top
fi


if [[ $STEP == $NOW ]]
then

    LOG=04-SOLVATION.log

    # Mark redundant files generated in this step
    trash $base-sol-b4ions.{gro,top,tpr} sol-b4ions.ndx genion_node0.log empty.mdp defaults.mdp


    ## 1. Solvation

    # a. Basic stuff
    [[ $GMXVERSION -gt 4 ]] && SOLVATE="$GMXBIN/gmx solvate" || SOLVATE=$GMXBIN/genbox
    SOLVATE="$SOLVATE -cp $GRO -cs -o $base-sol-b4ions.gro"

    # b. Add program specific options from command line
    SOLVATE="$SOLVATE $(program_options genbox)"
    SOLVATE="$SOLVATE $(program_options solvate)"

    # c. Make noise
    echo $SOLVATE | tee $LOG
    
    # d. Execute
    $EXEC $SOLVATE >>$LOG 2>&1

    # e. Extract number of atoms before adding solvent and amount of solvent added from log
    SED_SOL="/^Generated solvent/s/^.*in \(.*\) residues/\1/p"
    SED_ATOMS="/^Containing/s/^Containing \(.*\) atoms.*/\1/p"
    #NATOMS=($(LSED -n -e "$SED_SOL" -e "$SED_ATOMS" $LOG))
    #NSOL=${NATOMS[1]}
    NSOL=$(LSED -n -e "$SED_SOL" $LOG)

    # f. Update topology: add the number of solvent molecules added
    cp $TOP $base-sol-b4ions.top
    printf "$SolName %17d ; B4IONS\n" $NSOL >> $base-sol-b4ions.top

    # g. Add solvent model include file if it is not present yet
    #    First check if there is a moleculetype named 'SOL'
    #    This is a bit awkward and may give wrong results if we use 
    #    solvent other than water, but it is necessary to prevent 
    #    redefining the SOL moleculetype in certain (eTox) cases.
    MOLTYPES="$(awk "$AWK_MOLTYPE" $TOP)"
    if [[ ! $MOLTYPES =~ $SolName ]] 
    then
	# check if a file is #included with the solvent model
	if ! grep -q '#include.*'$SolventTopology $TOP 
	then
	    # Check if the topology for the solvent is here or there
	    #[[ -f $SolventTopology ]] || SolventTopology=$ForceField.ff/$SolventTopology
	    $SED -i.bck '/^\[ *system *\]/s,^,#include "'$SolventTopology$'"\\\n\\\n,' $base-sol-b4ions.top
	fi
    fi

    # h. Make some more noise
    echo "# Solvent added: $NSOL molecules"
	

    ## 2. Adding ions

    #
    # NOTES
    # 
    # One could argue that the following is taking bash over the top. And that is true.
    # It would be easier and neater just using python. In addition, it could even be left
    # to 'genion', which has options for neutralization and specification of the 
    # concentration. 
    # Yet there are several reasons for doing this here. First of all, the reason for 
    # doing it ourselves is that genion uses the box volume to calculate the number of
    # ions to add, given the volume. This assumes that the solvent in the box has 
    # equilibrium density already, and that the space occupied by solute counts as space
    # to consider for ions. The latter neglects the difference between macroscopic and 
    # microscopic solutions of macromolecules and the fact that such macromolecules are 
    # solvated in an isotonic solution, rather than brought to isotonicity afterwards,
    # based on the volume of the macromolecule solution. 
    # Another reason for doing it here is ... because we can :) Yes, it is showing off
    # to some extent, but it also demonstrates some features of bash, in particular 
    # related to integer math, which may come in handy.
    # Admittedly, this is not the most efficient way to handle this. But compared to the
    # simulation bits, the overhead is limited.
    #
    # On to adding ions...
    #


    # a. First to extract net charge
    #    Combine the topology and structure, yielding a full listing of atoms and charges
    echo nsteps=1 > empty.mdp
    tag="sol-b4ions"
    GROMPP="${GMX}grompp -v -f empty.mdp -c $base-$tag.gro -p $base-$tag.top -o $base-$tag.tpr -po defaults.mdp -maxwarn -1"


    # b. Tell what is happening
    echo $GROMPP | tee -a $LOG


    # c. Execute
    $NOEXEC $GROMPP >>$LOG 2>&1


    # d. Get the charge of the system
    #    For LIE calculations only calculate the charge excluding the ligand
    #    At this point, that corresponds to the 'ligand environment' group
    #    since we have not added the solvent to it.
    #    So first we create an index group, if needed.
    NDX=
    if $LIE
    then
	echo "# LIEing: Excluding ligand from charge calculation"
	# Of course, this only makes sense if we have something other than ligand
	if [[ -n $Ligenv ]]
	then
	    printf "%5d %5d %5d %5d %5d\n" $(SEQ ${Ligenv[@]}) | sed $'s/ 0//g;1s/^/[ check ]\\\n/' > charge.ndx
	    NDX="-n charge.ndx"
	    [[ $GMXVERSION -gt 4 ]] && TPRCONV="$GMXBIN/gmx tpr-convert" || TPRCONV="$GMXBIN/tpbconv"
	    $TPRCONV -s $base-sol-b4ions.tpr -o $base-sol-b4ions-noligand.tpr -n charge.ndx >/dev/null 2>&1
            NCHARGE=$(getCharge $base-sol-b4ions-noligand.tpr)
	    trash charge.ndx $base-sol-b4ions-noligand.tpr
	else
	    NCHARGE=0
	fi
    else
	NCHARGE=$(getCharge $base-sol-b4ions.tpr)
    fi
    echo "# Net charge of system: $NCHARGE"


    # e. Check if we should neutralize the system ... and whether we do
    [[ -n $CHARGE ]] && echo "# Setting charge to $CHARGE, while system charge is $NCHARGE" && NCHARGE=$CHARGE
    [[ $NCHARGE != 0 && ${Salinity:0:1} == - ]] && NCHARGE=0 && echo "# Not adding counterions, despite charge."


    # f. Calculate NPOS and NNEG given the charge and salinity, correcting for water removed

    #    i.    Salt: aX^+m, bY^n
    #          Infer stoichiometry from names - Check if salt name starts with a number
    [[ ${Salt[0]} =~ ^[0-9] ]] && a=${Salt[0]%%[^0-9]*} || a=1
    [[ ${Salt[1]} =~ ^[0-9] ]] && b=${Salt[1]%%[^0-9]*} || b=1

    #    ii.   Set names for ions - Strip the number in front (if any)
    PNAM=${Salt[0]#$a}
    NNAM=${Salt[1]#$b}    
	
    #    iii.  Now if Q is the number of ions to be added, and S is the number of solvent molecules, then
    #              Q = PS / 55.4(1+P)
    #          if the solvent is water and has a molarity of 55.4, and P = C(a+b) is the number of
    #          ions per liter, with C being the concentration of the salt.
    #          We cheat by calling python for floating point arithmetics.
    #          Yes, we could have used bc...
    Q=$(python -c "print int(0.5+$Salinity*($a+$b)*$NSOL/(55.4+$Salinity*($a+$b)))")

    #    iv.   Let U and V denote the number of positive and negative ions, respectively, then for
    #          a system with a net charge equal to Z we have
    #             Z = Um + Vn           (u and v integer)
    #             Q = U  + V
    #             U = (Qn - Z)/(n - m)
    #             V = (Z - Qn)/(m - n)
    m=${SaltCharge[0]}
    n=${SaltCharge[1]}
    U=$(python -c "print int(0.5+( $Q*$n+$NCHARGE)/($n-$m))")
    V=$(python -c "print int(0.5+(-$Q*$m-$NCHARGE)/($n-$m))")

    #    v.    Now U and V may still be slightly off, so we correct iteratively
    #          If the net charge is odd and the difference between charges is
    #          even, then it will never converge...
    prev=(999 999)
    i=0
    while [[ $((-U*m-V*n)) != $NCHARGE ]] 
    do
	# If the charge is too low, increase the number of positives
	[[ $((-U*m-V*n)) -lt $NCHARGE ]] && : $((U--))
	# If the number is still too low, decrease the number of negatives
	[[ $((-U*m-V*n)) -lt $NCHARGE ]] && : $((V++))
	# If the number is too high, increase the number of negatives
	[[ $((-U*m-V*n)) -gt $NCHARGE ]] && : $((V--))
	# If the number is still too high, decrease the number of positives
	[[ $((-U*m-V*n)) -gt $NCHARGE ]] && : $((U++))
	# Store the net charge with this configuration
	prev[${#prev[@]}]=$((-U*m-V*n))
	# Check if we are in an endless loop
	[[ $((-U*m-V*n)) == ${prev[$((i++))]} ]] && echo "# Breaking out of endless loop" && break
	ERROR="Problem calculating ion numbers to compensate net charge. More than 100 iterations run. Check the salt specification."
	[[ $i -gt 100 ]] && exit_error $ERROR
    done

    #    vi.   Check if one of the two is negative, and correct if so
    [[ $U -lt 0 ]] && : $((V-=U)) && U=0
    [[ $V -lt 0 ]] && : $((U-=V)) && V=0


    #NCL=$(python -c "print max(min(0,$NCHARGE),int(0.5+0.5*($Salinity*$NSOL/(27.7+$Salinity)+$NCHARGE)))")
    #NNA=$((NCL-NCHARGE))
    echo "# Replacing $(( U + V )) solvent molecules with $U $PNAM(+$m) and $V $NNAM($n) ions."

    # - Make an index file for the solvent added for genion
    echo "[ $SolName ]" > sol-b4ions.ndx

    # Make a listing of the solvent added with respect to the output from EM
    N1=$(( $(awk '{getline; print; exit}' $GRO) + 1 ))
    N2=$(awk '{getline; print; exit}' $base-sol-b4ions.gro)
    echo "[ $SolName ]" > sol.ndx
    printf "%5d %5d %5d %5d %5d\n" $(SEQ $N1 $N2) | $SED 's/ 0//g' >> sol.ndx
    
    # Only call genionif ions should be added
    if [[ $((U + V)) -gt 0 ]]
    then
        # - Then call genion
	#   Earlier versions of GMX use a -random flag. Check if this one does.
	RND=$(genion -h 2>&1 | sed -n 's/^\(-random\) .*/\1/p')
        GENION="${GMX}genion -s $base-sol-b4ions.tpr -o $base-sol.gro -n sol.ndx"
        GENION="$GENION -pname $PNAM -nname $NNAM -np $U -nn $V -pq $m -nq $n -rmin 0.5 $RND"

        echo $GENION | tee -a $LOG

        $GENION >>$LOG 2>&1
    
        # - And finally, update the topology
	#   Do check if 'ions.itp' is included or whether the ions are listed in the moleculetypes
	#   Either have all ions defined, or have none, and use ions.itp
	
	IONSITP=""
	if [[ $MOLTYPES =~ $PNAM && $MOLTYPES =~ $NNAM ]]
	then
	    echo "Moleculetype definitions for $PNAM and $NNAM found in topology file"
	elif [[ $MOLTYPES =~ $PNAM ]]
	then
	    exit_error "Moleculetype definition found in $TOP for ion $PNAM, but none found for $NNAM" 
	elif [[ $MOLTYPES =~ $NNAM ]]
	then
	    exit_error "Moleculetype definition found in $TOP for ion $NNAM, but none found for $PNAM" 
	else
	    N='\'$'\n'
	    grep -q '#include.*ions.itp' $base-sol-b4ions.top || IONSITP=$'/^\[ *system *\]/s,^,#include "'$ForceField.ff/'ions.itp"'$N$N','
	fi
        LSED -e "/B4IONS/d;$IONSITP;" $base-sol-b4ions.top > $base-sol.top
        printf "$SolName %17d\n$PNAM %17d\n$NNAM %17d" $(( NSOL - U - V )) $U $V >> $base-sol.top	
    else
        LSED "s/B4IONS//" $base-sol-b4ions.top > $base-sol.top
        cp $base-sol-b4ions.gro $base-sol.gro
    fi
	
    [[ -n $SCRATCH ]] && cp $LOG $base-sol.top $base-sol.gro $DIR

elif [[ $STEP == $STOP ]]
then
    $(all_exist ${OUTPUT[@]}) && echo -e "Output files present...\nSkipping: $GENBOX"
    $EXEC && echo -e "# Not executing:\n$GENBOX"
fi


## BOOKKEEPING ##

# bookkeeping is always done, also if the step is not executed 

# If we need the ligand interaction energy, set the energy groups to
# ligand and environment, otherwise set to solute and solvent
$LIE && EnergyGroups=(Ligand Ligenv) || EnergyGroups=(Solute Solvent)

# The temperature coupling groups are solute (including ligands)
# and solvent (including ions)
CoupleGroups=(Solute Solvent)

if [[ -e $DIR/$base-sol.gro ]]
then
    # Now see how much solvent was added in total and list in the right group
    # Also add to the ligand environment
    N2=$(awk '{getline; print; exit}' $DIR/$base-sol.gro)
    Solvent[${#Solvent[@]}]=${N1:-1}
    Solvent[${#Solvent[@]}]=$N2
    Ligenv[${#Ligenv[@]}]=${N1:-1}
    Ligenv[${#Ligenv[@]}]=$N2

    if [[ ! -e $DIR/$base-sol.ndx ]]
    then
        ## WRITE MASTER INDEX FILE ##
    
        # Here the whole system has been built. Time to make an index file with all the definitions needed

        # First a basic one
	echo q | ${GMX}make_ndx -f $DIR/$base-sol.gro -o $base-sol.ndx >/dev/null 2>&1

        # Add the Solute and Solvent (and Membrane?) groups
	fmt="%5d %5d %5d %5d %5d %5d %5d %5d %5d %5d"    
	echo "[ Solute ]" >> $base-sol.ndx
	for ((i=0; i<${#Solute[@]}; ))
	do
	    A=${Solute[$((i++))]}
	    B=${Solute[$((i++))]}
	    printf "$fmt\n" `SEQ $A $B` | $SED 's/ 0//g' >> $base-sol.ndx
	done
	echo "[ Membrane ]" >> $base-sol.ndx
	for ((i=0; i<${#Membrane[@]}; ))
	do
	    A=${Membrane[$((i++))]}
	    B=${Membrane[$((i++))]}
	    printf "$fmt\n" `SEQ $A $B` | $SED 's/ 0//g' >> $base-sol.ndx
	done
	echo "[ Solvent ]" >> $base-sol.ndx
	for ((i=0; i<${#Solvent[@]}; ))
	do
	    A=${Solvent[$((i++))]}
	    B=${Solvent[$((i++))]}
	    printf "$fmt\n" `SEQ $A $B` | $SED 's/ 0//g' >> $base-sol.ndx
	done
        # Finally add the ligand and ligenv groups if any ligands were added
	if [[ ${#LIGANDS[@]} -gt 0 ]]
	then
	    echo "[ Ligand ]" >> $base-sol.ndx
	    for ((i=0; i<${#Ligand[@]}; ))
	    do
		A=${Ligand[$((i++))]}
		B=${Ligand[$((i++))]}
		printf "$fmt\n" `SEQ $A $B` | $SED 's/ 0//g' >> $base-sol.ndx
	    done	
	    echo "[ Ligenv ]" >> $base-sol.ndx
	    for ((i=0; i<${#Ligenv[@]}; ))
	    do
		A=${Ligenv[$((i++))]}
		B=${Ligenv[$((i++))]}
		printf "$fmt\n" `SEQ $A $B` | $SED 's/ 0//g' >> $base-sol.ndx
	    done	
	fi

        [[ -n $SCRATCH ]] && cp $base-sol.ndx $DIR
    fi
fi

[[ $STOP ==   $NOW     ]] && exit_clean
[[ $STEP == $((NOW++)) ]] && : $((STEP++))


#--------------------------------------------------------------------
SHOUT "---STEP 5: ENERGY MINIMIZATION IN SOLVENT"
#--------------------------------------------------------------------

trash $base-EMs.{tpr,edr,log,trr} em-sol-out.mdp

# Turn on PBC for EM
__mdp_em__pbc=xyz

MDP=em-sol.mdp
OPT=(em)
TOP=$DIR/$base-sol.top
GRO=$DIR/$base-sol.gro
OUT=$base-EMs.gro
LOG=05-EMs.log
STRT=$(date +%s)

MD="MDRUNNER -f $MDP -c $GRO -p $TOP -o $OUT -np $NP -l $LOG -force $FORCE"

[[ $STEP ==   $NOW     ]] && mdp_options ${OPT[@]} > $MDP
[[ $STEP ==   $NOW     ]] && $EXEC $MD && : $((STEP++)) && archive

[[ -n $SCRATCH && -e $OUT ]] && cp $LOG $OUT ${OUT%.gro}.edr $DIR

[[ $STOP == $((NOW++)) ]] && exit_clean

LASTRUN=$(( $(date +%s) - STRT ))


#----------------------------------------------------------------------------------
SHOUT "---STEP 6: POSITION RESTRAINT MD, NVT -- CYCLE THROUGH PRFC AND TEMP/TAU_T"
#----------------------------------------------------------------------------------

if [[ $Electrostatics == PME ]]
then
    # Do not use a force field tag; just use mdp defaults
    ForceFieldMDP=
else
    ForceFieldMDP=$ForceFieldFamily
fi

__mdp_md__tc_grps=$(    $SED 's/ /,/g' <<< ${CoupleGroups[@]})
__mdp_md__energygrps=$( $SED 's/ /,/g' <<< ${EnergyGroups[@]})

# Only print stuff if we actually execute this step

if [[ $STEP == $NOW ]]
then
    echo "#"
    echo "#            Equilibration (NVT/PR):"
    echo "#                Coupling:                ${CoupleGroups[@]}"
    echo "#                Temperatures:            ${Temperature[@]}"
    echo "#                Coupling times:          ${Tau_T[@]}"
    echo "#                Position restraint Fcs:  ${PosreFC[@]}"
    echo "#"
fi

NDX=$DIR/$base-sol.ndx
TOP=$DIR/$base-sol.top
GRO=$DIR/$base-EMs.gro

# To avoid too much overhead, identify files with a position restraint tag
posre=($(grep -l POSREFC $DIR/*top $DIR/*itp))

# Counters
i=0
j=0
k=0
m=1
# Numbers
I=$((${#Temperature[@]}-1))
J=$((${#Tau_T[@]}-1))
K=$((${#PosreFC[@]}-1))
while :
do

    # Set the temperature and the position restraint force for this cycle
    T=${Temperature[$i]}
    tau=${Tau_T[$j]}
    F=${PosreFC[$k]}
    STRT=$(date +%s)

    # Specify temperature for each group (Solute/Solvent/Membrane)
    # Note that this (re)sets the master temperature control
    __mdp_md__ref_t=$($SED 's/ /,/g' <<< $(for i in ${CoupleGroups[@]}; do echo $T;   done))
    __mdp_md__tau_t=$($SED 's/ /,/g' <<< $(for i in ${CoupleGroups[@]}; do echo $tau; done))

    if [[ $STEP == $NOW ]]
    then
        # Make some noise
	LINE "NVT Equilibration at $T Kelvin (tau_t=$tau) with Position restraint force $F"

        # Modify position restraint definitions
	LSED -i -e "s/^\( *#define \+POSREFC\).*\$/\1 $F/" ${posre[@]}
    fi

    # Set the file names and such for this cycle
    MDP=pr-$F-nvt-$T-$tau.mdp
    OPT=(md $ForceFieldMDP equil usr)
    OUT=$base-PR-$F-NVT-$T-$tau.gro
    LOG=06-PR-NVT-$((m++)).log

    # Build the command
    MD="MDRUNNER -f $MDP -c $GRO -p $TOP -o $OUT -n $NDX -np $NP -l $LOG -force $FORCE $MONALL"

    # Execute
    [[ $STEP ==   $NOW     ]] && mdp_options ${OPT[@]} > $MDP
    [[ $STEP ==   $NOW     ]] && $EXEC $MD 

    # Mark the time
    LASTRUN=$(( $(date +%s) - STRT ))

    [[ -n $SCRATCH && -e $OUT ]] && cp $LOG $OUT ${OUT%.gro}.edr $DIR

    # Set current structure to last output structure
    GRO=$DIR/$OUT

    # Disable generation of velocities after first cycle
    __mdp_equil__genvel=no

    # Break the loop if all values are at maximum
    [[ $i == $I && $j == $J && $k == $K ]] && break

    # Increment the counters if there are more entries for it
    # A note for whoever is reading this: ':' is a bash command
    # which is just 'true', and can be used conveniently for 
    # doing in-place math operations using $(( )).
    [[ $i -lt $I ]] && : $(( i++ ))
    [[ $j -lt $J ]] && : $(( j++ ))
    [[ $k -lt $K ]] && : $(( k++ ))
  
done


# Store intermediate results if we just finished a step
[[ $STEP ==   $NOW     ]] && : $((STEP++)) && archive

# Exit if this step is the stop step
[[ $STOP == $((NOW++)) ]] && exit_clean


#----------------------------------------------------------------------------
SHOUT "---STEP 7: UNRESTRAINED MD 20 ps NPT -- CYCLE THROUGH PRESSURE/TAU_P"
#----------------------------------------------------------------------------

# Turning on the pressure
__mdp_md__pcoupl=Berendsen

# Relieve position restraints
__mdp_equil__define=


if [[ $STEP == $NOW ]]
then
    echo "#"
    echo "#            Equilibration (NpT):"
    echo "#                Pressures:       ${Pressure[@]}"
    echo "#                Coupling times:  ${Tau_P[@]}"
    echo "#"
fi

TOP=$DIR/$base-sol.top
GRO=$DIR/$OUT 

# Counters
i=0
j=0
m=1
# Numbers
I=$((${#Pressure[@]}-1))
J=$((${#Tau_P[@]}-1))
while :
do
    # Set the pressure and coupling time for this cycle
    P=${Pressure[$i]}
    tau=${Tau_P[$j]}
    STRT=$(date +%s)


    # Specify pressure
    # Note that this (re)sets the master pressure control
    __mdp_md__ref_p=$P
    __mdp_md__tau_p=$tau


    # Make some noise if we are actually executing this step
    [[ $STEP == $NOW ]] && LINE "NpT Equilibration at $P bar (tau_p=$tau)"


    # Set the file names and such for this cycle
    MDP=npt-$P-$tau.mdp
    OPT=(md $ForceFieldMDP equil $RotationalConstraints usr)
    OUT=$base-NPT-$P-$tau.gro
    LOG=07-NPT-$((m++)).log


    # Build the command
    MD="MDRUNNER -f $MDP -c $GRO -p $TOP -o $OUT -n $NDX -np $NP -l $LOG -force $FORCE $MONALL"


    # Execute
    [[ $STEP ==   $NOW     ]] && mdp_options ${OPT[@]} > $MDP
    [[ $STEP ==   $NOW     ]] && $EXEC $MD 

    # Mark the time
    LASTRUN=$(( $(date +%s) - STRT ))

    [[ -n $SCRATCH && -e $OUT ]] && cp $LOG $OUT ${OUT%.gro}.edr $DIR

    # Set current structure to last output structure
    GRO=$DIR/$OUT

    # Break the loop if all values are at maximum
    if [[ $i == $I && $j == $J ]]
    then
	break
    fi

    # Increment the counters if there are more entries for it
    # A note for whoever is reading this: ':' is a bash command
    # which is just 'true', and can be used conveniently for 
    # doing in-place math operations using $(( )).
    [[ $i -lt $I ]] && : $(( i++ ))
    [[ $j -lt $J ]] && : $(( j++ ))
  
done


# Store intermediate results if we just finished a step
[[ $STEP ==   $NOW     ]] && : $((STEP++)) && archive

# Exit if this step is the stop step
[[ $STOP == $((NOW++)) ]] && exit_clean


#--------------------------------------------------------------------
SHOUT "---STEP 8: SHORT RUN UNDER PRODUCTION CONDITIONS"
#--------------------------------------------------------------------

MDP=md-pre.mdp
OPT=(md $ForceFieldMDP equil $RotationalConstraints usr)
TOP=$DIR/$base-sol.top
GRO=$DIR/$OUT
OUT=$base-MD-PRE.gro
LOG=08-MD-PRE.log
STRT=$(date +%s)

MD="MDRUNNER -f $MDP -c $GRO -p $TOP -o $OUT -n $NDX -np $NP -l $LOG -force $FORCE"

[[ $STEP ==   $NOW     ]] && mdp_options ${OPT[@]} > $MDP
[[ $STEP ==   $NOW     ]] && $EXEC $MD && : $((STEP++)) && archive

[[ -n $SCRATCH && -e $OUT ]] && cp $LOG $OUT ${OUT%.gro}.edr $DIR

[[ $STOP == $((NOW++)) ]] && exit_clean

LASTRUN=$(( $(date +%s) - STRT ))

#--------------------------------------------------------------------
SHOUT "---STEP 9: PRODUCTION RUN"
#--------------------------------------------------------------------

MDP=md-prod.mdp
OPT=(md $ForceFieldMDP $RotationalConstraints usr)
TOP=$DIR/$base-sol.top
GRO=$DIR/$OUT
OUT=$base-MD.gro
LOG=09-MD.log

# If the run length is set to zero, write an mdp file 
# with nsteps=-1, generate the .tpr, and make a clean exit
if [[ $TIME == 0 ]]
then
    __mdp_cg__nsteps=-1
    OUT=$base-MD.tpr
fi

# If the Stop-Step is set to TPR only write a tpr and do not run
# Do set the STOP to the next STEP, as we still need to run MDRUNNER
[[ $STEP ==   $NOW     ]] && : $((STEP++))
[[ $STOP == $((NOW++)) ]] && OUT=$base-MD.tpr && : $((STOP++))

# If we have a TPR file for this step 
[[ $STEP == $NOW && -n $TPR ]] && GRO=$TPR && OUT=${TPR%.tpr}.gro && TPR=

MD="MDRUNNER -f $MDP -c $GRO -p $TOP -o $OUT -n $NDX -np $NP -l $LOG -force $FORCE -split -monitor"

[[ $STEP == $NOW ]] && mdp_options ${OPT[@]} > $MDP

# The following is like the usual mantra for running the simulations, but
# with a minor addition to see if the run finished to completion or not.
# If the first part is true (STEP==NOW) the expression in curly braces is
# executed, running the simulation and setting DONE to the exit code of
# MDRUNNER (0 if complete, 1 if not). The assignment is always true, so the
# step is increased and the stuff is archived, regardless of the result of 
# the run.
[[ $STEP == $NOW ]] && { $EXEC $MD; DONE=$?; } && archive

# Here we need to copy the whole run TPR/CPT/*.part*.*/GRO
# Since it is impossible to have multiple runs in one scratch
# directory, we can just do
[[ -n $SCRATCH ]] && cp $LOG ${OUT%.*}* $DIR

# We stop here if the run has not finished to completion
[[ $DONE == 1 ]] && touch INCOMPLETE && exit_clean

#--------------------------------------------------------------------
SHOUT "---DONE SIMULATING"
#--------------------------------------------------------------------

# We stop here if this was the last bit to do
[[ $STOP == $((NOW++)) ]] && exit_clean

# Otherwise we just continue with some more fun stuff
# If we have stuff to do, actually.

# If there is no analysis to run we skip the last bit
[[ -n $ANALYSIS ]] || : $((STEP++))


#--------------------------------------------------------------------
SHOUT "---STEP 10: ANALYSIS"
#--------------------------------------------------------------------

: << __NOTES__

This section contains some analyses to run on the data obtained. 
In general this should be restricted to relatively lightweight 
analysis and preferrably be controlled by the command line. 
Initially 

__NOTES__


# Convert the array of analyses into a string for matching
ANALYSIS=$(sed 's/ /./g' <<< .${ANALYSIS[@]}.)


if [[ $ANALYSIS =~ .LIE. ]]
then
    # Process the energy file(s) to extract the temperature and 
    # interaction energy between the ligand and the environment.

    # These are the energy terms to extract from the energy file
    terms=(
	Temperature
	Coul-SR:Ligand-Ligenv
	LJ-SR:Ligand-Ligenv
	Coul-LR:Ligand-Ligenv
	LJ-LR:Ligand-Ligenv
	0
    )
    # Process with sed to add newlines
    terms=$(sed $'s/ /\\\n/g' <<< ${terms[@]})

    # The production MD is run in parts; Process each part
    ls $base-MD.*edr >/dev/null 2>&1 || exit_error "LIE ANALYSIS ERROR: No energy \(.edr\) files found."

    for edr in *-MD.part*.edr 
    do
        [[ $GMXVERSION -gt 4 ]] && ENE="$GMXBIN/gmx energy" || ENE=$GMXBIN/g_energy
	echo $terms | $ENE -f $edr -o ${edr%.edr}.xvg 2>/dev/null | sed '/^Energy/,/^ *$/{/^ *$/q}' > ${edr%.edr}.lie
    done
fi


#-----------------------------------------------------------------------------------
SHOUT "Huh!? You made it all the way to the end of the script. This can't be right."
#-----------------------------------------------------------------------------------

archive
exit_clean

