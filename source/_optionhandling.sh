hlevel=1
olevel=1
# Function for parsing help from option list
help_fun()
{
  local desc='/^.*#=[0-'${hlevel}']/s//    /p'
  local item='/#==[0-'${olevel}']/{s/^ */        /;s/).*#==./  /p;d;}'
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

  [[ $1 -lt 0 ]] || exit $1
}


BAD_OPTION ()
{
  echo
  echo "Unknown option "$1" found on command-line"
  echo "It may be a good idea to read the usage:"

  USAGE -1

  echo " /\                                               /\ " 
  echo "/||\  Unknown option "$1" found on command-line  /||\ "
  echo " ||   It may be a good idea to read the usage     || "
  echo " ||                                               || "

  exit 1
}


# Functions for handling argument lists for downstream programs

# Expanding options from '-opt=val1\,val2' to '-opt val1 val2', not changing long options (--long-opt=val)
function expandOptList()
{
    for i in $@
    do
	[[ $i =~ --+ ]] && echo $i || (j=${i/=/ }; echo ${j//\\,/ })
    done
}

# Condensing options from '-opt1=val1,val2 -opt2=val3,val4' to '-opt1=val1\,val2,-opt2=val3\,val4' 
function condenseOptList()
{
    echo $(sed 's/,/\,/g;s/  */,/g;' <<< $@)
}

# Reading condensed options from the command line
function readOptList()
{
    sed "s/\\\,/##/g;s/,/ /g;s/##/,/g;s/--[^{]\+{\(.*\)}/\1/;" <<< $1
}

# Echo options from an array (for programs or MDP file)
function echo_additional_options() {
    if [[ -n $1 ]]
    then
	echo
	echo "# $MSG"
	while [[ -n $1 ]]
	do
	    echo "#     $1"
	    shift
	done
	echo "# ==="
    fi
}

