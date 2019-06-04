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
