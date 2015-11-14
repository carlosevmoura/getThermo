#!/bin/bash

### runexamples.sh - A script to execute examples of getThermo Fortran program
### A program to calculate Partition Functions and Gibbs Free Energy from selected Normal Modes from GAMESS
### Grupo de Espectroscopia Teórica e Modelagem Molecular - GETMM (IQ/UFRJ)

echo -e "\033[01;32mgetThermo examples in execution...\033[00m"  	    

# Example 1
cd example1/
../../bin/getThermo example1 > /dev/null

#Example 2
cd ../example2
../../bin/getThermo example2 > /dev/null

echo -e "\033[01;32mgetThermo examples were sucesfully executed. Compare with results in output folders.\033[00m"