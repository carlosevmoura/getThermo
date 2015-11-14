#!/bin/bash

### install.sh - A script to install getThermo program
### Partition Functions and Gibbs Free Energy from selected Normal Modes from GAMESS
### Grupo de Espectroscopia Teórica e Modelagem Molecular - GET (DFq/IQ/UFRJ)

### gFORTRAN Compiler
FCOMP="gfortran"

$FCOMP -ffree-line-length-none getThermo.f90 -o getThermo.x
mv getThermo.x ../bin/
