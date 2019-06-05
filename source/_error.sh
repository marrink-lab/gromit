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

