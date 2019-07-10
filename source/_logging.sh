## Macros to echo stuff only if the step is to be executed -- more about steps further down

# Terminal coloring
RED='\x1b[1;31m'
YEL='\x1b[1;33m'
GRN='\x1b[1;32m'
CYA='\x1b[1;36m'
OFF='\x1b[0m'

## Set the correct sed version for multi-platform use
# Also try to avoid GNU specific sed statements for the
# poor bastards that are stuck with one of those Mac things
SED=$(which gsed || which sed)

## Sed wrapper (loud sed) - echo the command line before running it
function LSED() { echo $SED "$@" 1>&2; $SED "$@"; }

## Functions for messaging

# Just echo the line (if the step asks for it)
function ECHO()  { [[ $STEP == $NOW ]] && echo "$@"; }

# Print with a line below
function LINE()  { [[ $STEP == $NOW ]] && echo && $SED -e 's/^/ /;s/$/ /p;s/./-/g;' <<< "$@" | $SED 's/^/#/'; }

# Add whitespace before and after and print with a line above and below
function LINES() { $SED -e 's/^/ /;s/$/ /;h;s/./-/g;p;x;p;x;' <<< "$@"; }

# Print with lines and preceding '#'
function SHOUT() { [[ $STEP == $NOW ]] && echo && LINES "$@" | $SED 's/^/#/'; }

# Print in red with lines and preceding '# WARNING: '
function WARN()  { [[ $STEP == $NOW ]] && echo -e "$RED" && LINES "$@" | $SED 's/^/# WARNING: /' && echo -e "$OFF"; }

# Print in red with lines and preceding '# FATAL: ' and then exit (hey, that's not nice...)
function FATAL() { [[ $STEP == $NOW ]] && echo -e "$RED" && LINES "$@" | $SED 's/^/# FATAL: /' && echo -e "$OFF"; exit 1; }

# Print in yellow with lines and preceding '# NOTE: '
function NOTE()  { [[ $STEP == $NOW ]] && echo -e "$YEL" && LINES "$@" | $SED 's/^/# NOTE: /' && echo -e "$OFF"; }

## Execution

# Macro to do stuff only if the step is to be executed
function DO() { [[ $STEP == $NOW ]] && echo "$@" && $@; }

