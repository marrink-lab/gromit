# Structure databases
RCSB=(http://files.rcsb.org/download/ pdb.gz)
BIO=(http://www.rcsb.org/pdb/files/ pdb1.gz)
OPM=(http://opm-assets.storage.googleapis.com/pdb/ pdb)

function fetch_structure() # structure source
{
    pdb=$1
    FETCH=$2
    
    RMMOD=false
    RMDUM=false

    # Try fetching it from the repository
    # Biological units use MODEL/ENDMDL statements to separate chains
    # These need to be removed.
    case $FETCH in
        pdb  | PDB)  
            GET=${RCSB[0]}${pdb%%.*}.${RCSB[1]}
            ;;
        bio* | BIO*)
            # Check if a number is given for selecting the BIO assembly
            [[ $FETCH =~ : ]] && BIO[1]=pdb${FETCH##*:} 
            GET=${BIO[0]}${pdb%%.*}.${BIO[1]}
            RMMOD=true		
            ;;
        opm  | OPM)  
            GET=${OPM[0]}${pdb%%.*}.${OPM[1]}
            RMDUM=true
            ;;
    esac

    echo "# Input file not found, but will try fetching it ($GET)"
    # Use wget or curl; one of them should work
    ERRMSG="Could not retrieve input structure ($GET) from $FETCH database."
    wget $GET 2>/dev/null || curl -O $GET || exit_error "$ERRMSG"
    # Remove the URL part of what we just got
    pdb=${GET##*/}
    [[ $GET =~ \.gz$ ]] && gunzip $pdb && pdb=${pdb%.gz} 
    # With biological units, Gromacs chokes on the pdb1 extension
    [[ $pdb =~ \.pdb1$ ]] && mv $pdb ${pdb%1} && pdb=${pdb%1}
    PDB=$pdb

    if $RMMOD
    then
        # Removing the MODEL records and replacing ENDMDL by TER
        $SED -i '' -e /^MODEL/d -e s/ENDMDL/TER/ $PDB
    fi

    if $RMDUM
    then
	# Removing DUM particles, which are added by OPM to mark the membrane
	$SED -i '' -e '/^HETATM.\{11\}DUM/d' $PDB
    fi
}
