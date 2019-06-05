#--------------------------------------------------------------------
# Global parameters

__mdp_md__dt=0.002
[[ $GMXVERSION -gt 4 ]] && __mdp_md__nstlist=10 || __mdp_md__nstlist=5
[[ $GMXVERSION -gt 4 ]] && __mdp_md__cutoff_scheme=Verlet

# Virtual site specific
if $VirtualSites; then
  __mdp_md__dt=0.004
  __mdp_md__lincsorder=5
  __mdp_md__nstlist=2
fi

# Output parameters
TIME=$(python -c "print int(1000*$TIME/$__mdp_md__dt + 0.5 )") 
AT=$(python -c "print int(1000*$AT/$__mdp_md__dt + 0.5)") 
__mdp_md__nsteps=$TIME
__mdp_md__nstxout=$AT
__mdp_md__nstvout=0 
__mdp_md__nstfout=0
__mdp_md__nstlog=$AT
__mdp_md__nstenergy=$AT
[[ $GMXVERSION -gt 4 ]] && __mdp_md__nstxout_compressed=$AT || __mdp_md__nstxtcout=$AT

# Coupling
# Listed here are the overall controls. Specific controls
# (temperature,pressure,time constants) are controlled through
# the command line interface and taken care of in steps 6/7
# The temperature is set later on
__mdp_md__tcoupl=v-rescale
__mdp_md__nsttcouple=$__mdp_md__nstlist
__mdp_md__nstpcouple=$__mdp_md__nstlist
__mdp_md__compressibility=4.5e-5

# Nonbonded interactions
__mdp_md__coulombtype=PME
__mdp_md__fourierspacing=0.125
__mdp_md__rcoulomb=0.9
__mdp_md__rlist=0.9
__mdp_md__rvdw=1.4
__mdp_md__pme_order=4

# Other
__mdp_md__constraints=all-bonds
__mdp_md__comm_mode=Linear
__mdp_md__comm_grps=System
__mdp_md__nstcomm=$__mdp__nstlist


#--------------------------------------------------------------------
# Rotational constraints

__mdp_rtc__comm_mode=RTC
__mdp_rtc__comm_grps=Solute

#--------------------------------------------------------------------
# Force field specific parameters

# AMBER, all versions 
# (parameters taken from Acpype protein-ligand tutorial by Alan Wilter ...)
__mdp_amber__coulombtype=PME
__mdp_amber__coulomb_modifier=Potential-shift-Verlet
__mdp_amber__rcoulomb=1
__mdp_amber__rlist=1
__mdp_amber__vdwtype=cut-off
__mdp_amber__rvdw=1.0
__mdp_amber__DispCorr=EnerPres
__mdp_amber__lincs_iter=2
#__mdp_amber__optimize_fft=yes
# Additional options specified in named tutorial:
# __mdp_amber__pcoupltype=Parinello-Rahman
# __mdp_amber__tau_p=0.5

# GROMOS96, all versions
__mdp_gromos__coulombtype=Reaction-Field
__mdp_gromos__rcoulomb=1.4
__mdp_gromos__epsilon_rf=54

# CHARMM
if [[ $GMXVERSION -gt 4 ]]
then
  __mdp_charmm__constraints=h-bonds
  __mdp_charmm__cutoff_scheme=Verlet
  __mdp_charmm__vdwtype=cutoff
  __mdp_charmm__vdw_modifier=force-switch
  __mdp_charmm__rlist=1.2
  __mdp_charmm__rvdw=1.2
  __mdp_charmm__rvdw_switch=1.0
  __mdp_charmm__coulombtype=PME
  __mdp_charmm__rcoulomb=1.2
  __mdp_charmm__DispCorr=no
else
  __mdp_charmm__coulombtype=Switch
  __mdp_charmm__rcoulomb_switch=1.0
  __mdp_charmm__rcoulomb=1.2
  __mdp_charmm__vdwtype=Switch
  __mdp_charmm__rvdw_switch=1.0
  __mdp_charmm__rvdw=1.2
fi

# OPLS/AA
__mdp_opls__coulombtype=PME
__mdp_opls__rcoulomb=1.0
__mdp_opls__rlist=1.0
__mdp_opls__rvdw=1.0
__mdp_opls__fourierspacing=0.135
__mdp_opls__dt=0.001

#--------------------------------------------------------------------
# Energy minimization

__mdp_em__define=-DFLEXIBLE
__mdp_em__integrator=steep
__mdp_em__nsteps=$EMSteps
__mdp_em__pbc=xyz
__mdp_em__constraints=none

#--------------------------------------------------------------------
# Equilibration runs: position restraints, NVT, NPT

# Position restraints are relieved at step 7
TIME=$(python -c "print int(1000*$EquilTime/$__mdp_md__dt + 0.5 )")
__mdp_equil__define=-DPOSRES
__mdp_equil__nsteps=$TIME
__mdp_equil__nstlog=10
__mdp_equil__nstenergy=10
[[ $GMXVERSION -gt 4 ]] && __mdp_equil__nstxout_compressed=0 || __mdp_equil__nstxtcout=0

# Velocities are only generated once
# After the first NVT/PR cycle 'genvel' is set to no
# and the other two options are ignored
__mdp_equil__genvel=yes
__mdp_equil__gen_seed=$SEED
__mdp_equil__gen_temp=${Temperature[0]}

