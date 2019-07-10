function init_mdp_parameters() {
    # This file contains the basis parameters for running simulations with Martini
    # (with or without multiscaling)
    
    # Common setting for time step and neighborlist update
    __mdp_cg__nstlist=10
    __mdp_cg__dt=0.020
    
    # Override for multiscale simulations
    if $ALL -o [[ -n $MULTI ]]
    then
        if $VSITE
        then
            __mdp_cg__dt=0.004
            __mdp_cg__nstlist=4
        else
            __mdp_cg__dt=0.002
        fi
    fi
    
    # Override if timestep is set explicitly
    if [[ -n $DELT ]]
    then
      __mdp_cg__dt=$DELT
    fi
    
    # Output parameters
    TIME=$(python -c "print int(1000*$TIME/$__mdp_cg__dt + 0.5 )") 
    AT=$(python -c "print int(1000*$AT/$__mdp_cg__dt + 0.5)") 
    __mdp_cg__nsteps=$TIME
    __mdp_cg__nstxout=$AT
    __mdp_cg__nstvout=0 
    __mdp_cg__nstfout=0
    __mdp_cg__nstlog=$AT
    __mdp_cg__nstenergy=$AT
    [[ $GMXVERSION -gt 4 ]] && __mdp_cg__nstxout_compressed=$AT || __mdp_cg__nstxtcout=$AT
    
    # Nonbonded interactions
    if [[ $GMXVERSION -gt 4 ]]
    then
      # GMX 5 parameters from De Jong, Baoukina and Marrink
      __mdp_cg__cutoff_scheme=Verlet
      __mdp_cg__coulombtype=Cut-off
      __mdp_cg__coulomb_modifier=Potential-shift
      __mdp_cg__rcoulomb=1.1
      __mdp_cg__epsilon_r=$EPSR_CG
      __mdp_cg__vdw_type=Cut-off
      __mdp_cg__vdw_modifier=Potential-shift
      __mdp_cg__rvdw=1.1
      __mdp_cg__dispcorr=No
      __mdp_cg__nstlist=20
    else
      __mdp_cg__coulombtype=Shift
      __mdp_cg__rcoulomb=1.2
      __mdp_cg__rcoulomb_switch=0.0
      __mdp_cg__epsilon_r=$EPSR_CG
      __mdp_cg__rlist=1.2
      __mdp_cg__vdw_type=Shift
      __mdp_cg__rvdw=1.2
      __mdp_cg__rvdw_switch=0.9
      __mdp_cg__dispcorr=No
    fi
    
    # Coupling - depending on presence of protein/membrane
    # Pressure coupling semiisotropic for membranes,
    # set to isotropic if no membrane is present
    __mdp_cg__tcoupl=v-rescale
    __mdp_cg__nsttcouple=$__mdp_cg__nstlist
    
    __mdp_cg__pcoupl=no # Overridden at step NpT
    __mdp_cg__nstpcouple=$__mdp_cg__nstlist
    
    __mdp_cg__pcoupltype=Isotropic
    __mdp_cg__compressibility=3e-4
    __mdp_cg__tau_p=3.0 
    __mdp_cg__ref_p=$Pressure
    
    __mdp_rtc__comm_mode=RTC
    __mdp_rtc__comm_grps=Solute
    
    __mdp_mem__pcoupltype=Semiisotropic
    __mdp_mem__compressibility=3e-4,3e-4
    __mdp_mem__tau_p=3.0,3.0 # Overridden at step NpT with GMX5
    __mdp_mem__ref_p=$Pressure,$Pressure
    
    # Gromacs can handle empty groups, so this is fine
    # These are set based on groups in the index file 
    #__mdp_cg__tc_grps=Solute,Membrane,Solvent
    #__mdp_cg__tau_t=1.0,1.0,1.0
    #__mdp_cg__ref_t=$Temperature,$Temperature,$Temperature
    
    # Other
    __mdp_cg__constraints=none
    __mdp_cg__energygrps=Solute,Membrane,Solvent
    __mdp_cg__comm_mode=Linear
    __mdp_cg__comm_grps=System
    __mdp_cg__nstcomm=$__mdp_cg__nstlist
    
    
    #--------------------------------------------------------------------
    # Multiscale
    __mdp_ms__rlist=$RC
    __mdp_ms__coulombtype=user
    __mdp_ms__rcoulomb=$RC
    __mdp_ms__vdw_type=user 
    __mdp_ms__rvdw=$RC
    __mdp_ms__energygrps=AA,CG,VZ
    __mdp_ms__epsilon_r=1
    if $POLARIZABLE
    then
      __mdp_ms__energygrp_table=AA,AA,AA,CG
      __mdp_ms__energygrp_excl=AA,VZ,VZ,VZ
    else
      __mdp_ms__energygrp_table=AA,AA
      __mdp_ms__energygrp_excl=AA,VZ,AA,CG,VZ,VZ
    fi
    
    #--------------------------------------------------------------------
    # Energy minimization
    
    __mdp_em__integrator=steep
    __mdp_em__nsteps=$EMSTEPS
    __mdp_em__emstep=0.001
    __mdp_em__lincs_order=6
    __mdp_em__lincs_iter=8
    __mdp_em__lincs_warnangle=90
    __mdp_em__define=-DFLEXIBLE
    
    #--------------------------------------------------------------------
    # Equilibration runs: position restraints, NVT, NPT
    
    # Position restraints are relieved at step 8
    __mdp_equil__define=-DPOSRES
    __mdp_equil__dt=0.002
    __mdp_equil__nsteps=5000
    __mdp_equil__nstlist=1
    __mdp_equil__nstlog=8
    __mdp_equil__nstenergy=8
    [[ $GMXVERSION -gt 4 ]] && __mdp_cg__nstxout_compressed= || __mdp_cg__nstxtcout=0
    __mdp_equil__tcoupl=v-rescale
    __mdp_equil__lincs_order=6
    __mdp_equil__lincs_iter=8
    __mdp_equil__lincs_warnangle=90
    
    __mdp_equil__tau_p=10
}
