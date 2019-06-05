# This function lists the mdp options requested based on a preceding tag
function mdp_options ()
{
  echo MDP tags: $@ 1>&2
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

