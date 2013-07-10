# TDPF #
**T**emperature- **D**ependent **P**ower **F**low implemented in MATLAB

Copyright &copy; 2013  Stephen Frank, Jason Sexauer, and Salman Mohagheghi

## What Is TDPF? ##
TDPF is a proof-of-concept Temperature-Dependent Power Flow algorithm implemented as a collection of MATLAB scripts. TDPF augments conventional power flow by integrating an estimate of branch temperatures with the conventional power flow equations. For the technical details, please read:

1.	S. Frank, J. Sexauer, and S. Mohagheghi, "Temperature-dependent power flow," *IEEE Transactions on Power Systems*, 2013, to be published. Available: [http://dx.doi.org/10.1109/TPWRS.2013.2266409](http://dx.doi.org/10.1109/TPWRS.2013.2266409  "Temperature-Dependent Power Flow by S. Frank, J. Sexauer, and S. Mohagheghi").

For those without access to an *IEEE Transactions on Power Systems subscription*, a pre-publication manuscript of the paper is available [here](http://files.stevefrank.info/pub/TDPF.pdf "Temperature-Dependent Power Flow by S. Frank, J. Sexauer, and S. Mohagheghi (Authors' manuscript)").

## How To Use TDPF ##
TDPF is implemented as a collection of MATLAB functions:

*	`importCaseData()` - Import power system data from CSV, [IEEE Common Data Format](http://dx.doi.org/10.1109/TPAS.1973.293571 "Description of IEEE Common Data Format"), or [MATPOWER](http://www.pserc.cornell.edu/matpower/ "MATPOWER Website")
*	`makeYBus()` - Create the system admittance matrix for TDPF
*	`evalJacobian()` - Evaluate the Jacobian matrix for TDPF
*	`evalMismatch()` - Evaluate power and temperature mismatches for TDPF
*	`PF()` - Execute a conventional Newton-Raphson or fast decoupled power flow
*	`FC_TDPF()` - Execute fully coupled temperature-dependent power flow
*	`PD_TDPF()` - Execute partially decoupled temperature-dependent power flow
*	`FD_TDPF()` - Execute fast decoupled temperature-dependent power flow
*	`SD_TDPF()` - Execute sequentially decoupled temperature-dependent power flow

The `.m` file for each function provides documentation of the function inputs, outputs, and options. The script `tdpf_example.m`, which provides a step-by-step example for executing a temperature dependent power flow for the 39 bus New England test system, is the best place to start.

## Citing TDPF ##
If you use TDPF in your research, we would appreciate it if you would cite the article referenced above.

## License ##
This program is free software: you may redistribute it and/or modify it under the terms of the GNU General Public License (GPL) as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

Text and HTML copies of the GNU General Public License should be distributed with these scripts in the files `gpl-3.0.txt` and `gpl-3.0.html`. If not, please visit [the GPL website](http://www.gnu.org/licenses/ "GNU General Public License").