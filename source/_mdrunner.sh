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
    local -i RUNSTEPS=$($SED -n '/nsteps/s/^.*=\([^;]*\).*$/\1/p' $fnMDP)

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
    GROMPP="${GMX}grompp $fnMDP -c $fnIN -r $fnIN -p $fnTOP $fnNDX -o $baseOUT.tpr $(program_options grompp) -maxwarn -1"

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
    MAXH=$(awk '{print $1/3600}' <<< "$MAXS") 
    
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
    ($MON | $SED 's/^/MONITOR: /'; exit ${PIPESTATUS[0]}) &

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
