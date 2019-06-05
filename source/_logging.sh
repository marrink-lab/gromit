# Macros to echo stuff only if the step is to be executed -- more about steps further down
RED='\x1b[1;31m'
YEL='\x1b[1;33m'
GRN='\x1b[1;32m'
CYA='\x1b[1;36m'
OFF='\x1b[0m'

# Set the correct sed version for multi-platform use
# Also try to avoid GNU specific sed statements for the
# poor bastards that are stuck with one of those Mac things
SED=$(which gsed || which sed)

LINES() { $SED -e 's/^/ /;s/$/ /;h;s/./-/g;p;x;p;x;' <<< "$@"; }
SHOUT() { [[ $STEP == $NOW ]] && echo && LINES "$@" | $SED 's/^/#/'; }
WARN()  { [[ $STEP == $NOW ]] && echo -e "$RED" && LINES "$@" | $SED 's/^/# WARNING: /' && echo -e "$OFF"; }
FATAL() { [[ $STEP == $NOW ]] && echo -e "$RED" && LINES "$@" | $SED 's/^/# FATAL: /' && echo -e "$OFF"; exit 1; }
NOTE()  { [[ $STEP == $NOW ]] && echo -e "$YEL" && LINES "$@" | $SED 's/^/# NOTE: /' && echo -e "$OFF"; }
LINE()  { [[ $STEP == $NOW ]] && echo && $SED -e 's/^/ /;s/$/ /p;s/./-/g;' <<< "$@" | $SED 's/^/#/'; }
LINE()  { [[ $STEP == $NOW ]] && echo && $SED -e 's/^/ /;s/$/ /p;s/./-/g;' <<< "$@" | $SED 's/^/#/'; }
ECHO()  { [[ $STEP == $NOW ]] && echo "$@"; }

# Macro to do stuff only if the step is to be executed
DO() { [[ $STEP == $NOW ]] && echo "$@" && $@; }

# Sed wrapper (loud sed) - echo the command line before running it
LSED() { echo $SED "$@" 1>&2; $SED "$@"; }

