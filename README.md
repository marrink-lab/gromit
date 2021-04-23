# gromit and martinate

Auxiliary tools for automated atomistic (gromit) and coarse-grained (martinate) molecular dynamics simulations using GROMACS

## Synopsis

Molecular dynamics simulations have complex workflows, including the generation of a model, setting up the environment, relaxation of the system and finally the production simulation. Despite the intrinsic complexity, the steps of 
the process are well-defined. For simulations of protein and/or DNA in solution, with or without ligand and with or without ions standard protocols are available. Gromit and martinate are versatile wrappers providing such 
protocols for atomistic (gromit) and coarse-grain (martinate) simulations, using the molecular simulation package Gromacs and, for martinate, the coarse grain Martini force field.

## Example

## Motivation

## Installation

Martinate requires GROMACS, insane.py script provided on martini website and the python package vermouth. Additionally, Martini3 forcefield needs to be downloaded to create systems in this.

You can install gromit, martinize2, insane.py and martinize2 (vermouth) by running:
```bash
git clone https://github.com/marrink-lab/gromit.git

# Put insane in gromit path
wget http://www.cgmartini.nl/images/tools/insane/insane.py \
    -O gromit/insane
chmod +x gromit/insane

# Install martinize2
pip install vermouth

# Download forcefield and mappings as given there:
# http://cgmartini.nl/index.php/force-field-parameters/particle-definitions
wget http://www.cgmartini.nl/images/martini_v300.zip
unzip -d martini_v300 martini_v300.zip
rm martini_v300.zip
```
## Test

## Contributors

## License
