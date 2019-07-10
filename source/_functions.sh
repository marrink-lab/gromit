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


function archive()
{
    if [[ -n $ARCHIVE ]]
    then
        tar cfz $ARCHIVE.tmp.tgz `ls --ignore=$ARCHIVE`
        mv $ARCHIVE.tmp.tgz $ARCHIVE
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
    if ! $KEEP
    then
        echo "# Deleting redundant files:"
        printf "%-25s %-25s %-25s %-25s %-25s\n" ${JUNK[@]} \#*\#
        for i in ${JUNK[@]} \#*\# *.bck; do [[ -f $i ]] && rm $i; done
    fi
    [[ -n $SCRATCH ]] && cd $DIR # && rm -rf $SCRATCH

    # Archive
    archive

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


# Trashing files
function trash()
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


# Routine for gathering program specific options
# In the version edited by MvD (UU) this routine 
# uses numbers: for ((opt=0; opt<${#...
# Do not understand why.
function program_options()
{
    local OPTS=
    for opt in ${PROGOPTS[@]}
    do
        if [[ $opt =~ --$1 ]]
        then
            OPTS="$OPTS $($SED 's/--[^-]*//;s/=/ /' <<< $opt)"
        fi
    done
    echo $OPTS
}


# This is a function to generate a sequence of numbers
# On Linux platforms there usually is 'seq', but we 
# can not depend on it.
function SEQ(){ for ((i=$1; i<=$2; i++)); do echo $i; done };



