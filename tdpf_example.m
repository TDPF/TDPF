%% TDPF Example
% This script provides a simple example for how to execute a temperature
% dependent power flow using the TDPF scripts for MATLAB. To execute this
% example, the TDPF function .m files must be either in the current MATLAB
% working directory or elsewhere in the MATLAB path.
%
% For the technical details of TDPF, please read:
%   S. Frank, J. Sexauer, and S. Mohagheghi, "Temperature-dependent power
%   flow," 2013, submitted for publication.
%   Available: http://files.stevefrank.info/pub/TDPF.pdf
%
% For more information regarding the various TDPF functions, please read
% the documentation within the .m files that define each function.


%% Data Import
% This example uses data from the MATPOWER 39-bus New England test system.
% For convenience, this data has already been converted to TDPF format and
% saved in the file '39bus_data.mat'. Therefore, MATPOWER does not need to
% be installed to run the example.

% Load the case data.
% NOTE: The file '39bus_data.mat' must be in the MATLAB working directory.
load('39bus_data.mat');

% If MATPOWER is installed and present in the MATLAB path, at your option
% you may instead import the case data directly from MATPOWER by executing
% the commands:
% >> casename = 'case39';
% >> [bus, branch, SBase, TBase] = importCaseData(casename,'MATPOWER');


%% Execute Conventional Power Flow
% This section executes conventional PF and fast-decoupled PF and displays
% the resulting mismatches.
%
% Outputs:
%   V = Voltage magnitudes
%   d = Voltage angles

% Execute conventional PF
[V1,d1] = PF(bus,branch);

% Execute fast-decoupled PF
[V2,d2] = PF(bus,branch,'FD',true);

% Compute mismatches in PF results [P, Q]
% (This requires first setting up P, Q, and H sets for evalMismatch() and
% computing the complex admittance matrix Y.)
sets.P = bus.id((bus.type == 0) | (bus.type == 1) | (bus.type == 2));
sets.Q = bus.id((bus.type == 0) | (bus.type == 1));
sets.H = [];

Y = makeYBus(bus,branch);

mm0 = evalMismatch(sets, bus.V_mag, bus.V_angle * (pi/180), ...
        [], Y, bus, branch); % Original case data
mm1 = evalMismatch(sets, V1, d1  * (pi/180), ...
        [], Y, bus, branch); % Computed results - Conventional PF
mm2 = evalMismatch(sets, V2, d2  * (pi/180), ...
        [], Y, bus, branch); % Computed results - Fast-decoupled PF
    
% Display computed mismatches
semilogy(abs([mm0, mm1, mm2]));
xlabel('State Vector Index');
ylabel('Mismatch Magnitude');
title('Mismatches - Conventional PF');
legend('Original Data','Conventional PF','Fast Decoupled PF');


%% Execute TDPF
% This section executes TDPF using each of the four possible variants:
%   1. FC_TDPF() - Fully Coupled TDPF
%   2. PD_TDPF() - Partially Decoupled TDPF
%   3. FD_TDPF() - Fast Decoupled TDPF
%   4. SD_TDPF() - Sequentially Decoupled TDPF
% and displays the resulting mismatches.
%
% Outputs:
%   V = Voltage magnitudes
%   d = Voltage angles
%   T = Branch temperatures

% Execute FC-TDPF
[V1,d1,T1,bus1,branch1] = FC_TDPF(bus, branch);

% Execute PD-TDPF
[V2,d2,T2,bus2,branch2] = PD_TDPF(bus, branch);

% Execute FD-TDPF
[V3,d3,T3,bus3,branch3] = FD_TDPF(bus, branch);

% Execute SD-TDPF
[V4,d4,T4,bus4,branch4] = SD_TDPF(bus, branch);

% Set up sets
sets.P = bus.id((bus.type == 0) | (bus.type == 1) | (bus.type == 2));
sets.Q = bus.id((bus.type == 0) | (bus.type == 1));
sets.H = branch.id((branch.type == true));

% Compute final admittance matrices
Y1 = makeYBus(bus1,branch1);
Y2 = makeYBus(bus2,branch2);
Y3 = makeYBus(bus3,branch3);
Y4 = makeYBus(bus4,branch4);

% Compute final mismatches [P, Q, H]
mm1 = evalMismatch(sets, V1, d1  * (pi/180), T1, Y1, bus1, branch1);
mm2 = evalMismatch(sets, V2, d2  * (pi/180), T2, Y2, bus2, branch2);
mm3 = evalMismatch(sets, V3, d3  * (pi/180), T3, Y3, bus3, branch3);
mm4 = evalMismatch(sets, V4, d4  * (pi/180), T4, Y4, bus4, branch4);

% Display computed mismatches
semilogy(abs([mm1, mm2, mm3, mm4]));
xlabel('State Vector Index');
ylabel('Mismatch Magnitude');
title('Mismatches - Temperature-Dependent PF');
legend('FC-TDPF','PD-TDPF','FD-TDPF','SD-TDPF');


%% History and Timings
% This section illustrates how to obtain history and timing data) from the
% TDPF functions, using FD_TDPF as an example.

% Execute FD-TDPF, enabling display of timings
[~] = FD_TDPF(bus, branch, 'timing', true);

% Execute FD-TDPF, logging the history
[~,~,~,~,~,hist] = FD_TDPF(bus, branch, 'history', true);

% The 'hist' structure contains state and mismatch data from each iteration
% of the power flow. For example, the following plots the maximum mismatch
% at each iteration using the history:
semilogy(hist.iter, hist.maxErr, 'ks-');
xlabel('Iteration');
ylabel('Maximum Mismatch');
title('FD-TDPF Convergence');


%% Branch Loss Comparison
% This section provides an example comparison of branch losses as computed
% by conventional PF and TDPF.

% Execute conventional PF
[V1,d1] = PF(bus,branch);
bus1 = bus; branch1 = branch;

% Execute FC-TDPF
[V2,d2,T2,bus2,branch2] = FC_TDPF(bus, branch);

% Compute branch losses
branch1 = evalBranchLoss(bus1,branch1);
branch2 = evalBranchLoss(bus2,branch2);

% Extract branch losses
branchloss1 = real(branch1.S_loss);
branchloss2 = real(branch2.S_loss);

% Round to 6 decimal places
% (This avoids numerical errors in pct. differences)
branchloss1 = round( branchloss1 .* 1e6 ) ./ 1e6;
branchloss2 = round( branchloss2 .* 1e6 ) ./ 1e6;

% Plot branch losses
plot(branch.id, [branchloss1; branchloss2]);
xlabel('Branch ID');
ylabel('Loss (pu)');
title('Branch loss comparison');
legend('Conventional PF','TDPF');

% Plot absolute difference
plot(branch.id, branchloss2 - branchloss1);
xlabel('Branch ID');
ylabel('Difference in Loss (pu)');
title('Absolute Difference in Branch Loss (TDPF vs. PF)');

% Plot percent difference
% (Note: 0/0 yields NaN, which is then replaced by 0% difference)
pctdiff = (branchloss2 - branchloss1) ./ branchloss1 * 100;
pctdiff( isnan(pctdiff) ) = 0;
plot(branch.id, pctdiff);
xlabel('Branch ID');
ylabel('Difference in Loss (%)');
title('Percent Difference in Branch Loss (TDPF vs. PF)');

