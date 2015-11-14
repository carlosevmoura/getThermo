# getThermo
A tool to obtain thermodynamic properties of chemical systems selecting normal modes which describes vibrational movement.

**getThermo** works with [GAMESS-US](http://www.msg.ameslab.gov/gamess/index.html) output files from Hessian calculation type and [Gaussian09](http://www.gaussian.com/g_prod/g09.htm) output files from Frequency calculation.

## Installation
To install **getThermo**, you should compile getThermo.f90 code, with your favourite FORTRAN compiler, and move the binary file to internal *bin* directory. We recommend to use [GFortran](https://gcc.gnu.org/fortran/).

The script *install.sh* can do this procedure automatically.

## Usage
To use **getThermo**, you just need to run the script *getThermo* in the same folder that contains the aim output files, setting the input file to be used.

Example: ``` bash $your-getThermo-dir/bin/getThermo $InputFile ```

Read the documentation PDF file *User's Guide* to obtain more informations about it.

## About authors
This application is developed by Grupo de Espectroscopia Teórica e Modelagem Molecular (GET) from Instituto de Química of Universidade Federal do Rio de Janeiro (IQ/UFRJ), coordinated by [professor Alexandre Rocha](mailto:rocha@iq.ufrj.br). 

The FORTRAN code and shell-scripts are developed by [Carlos Eduardo de Moura](mailto:carlosevmoura@iq.ufrj.br) and [Ricardo Oliveira](mailto:rrjunior@iq.ufrj.br).
Fell encouraged to contact us with any question about getThermo!

Distribution of **getThermo** is given on our [GitHub repository](https://github.com/carlosevmoura/getThermo), under [MIT License](https://opensource.org/licenses/MIT).
