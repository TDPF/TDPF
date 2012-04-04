%% Scratchpad: Tests version 2 of the TDPF algorithms
%% Setup
clear;
clc;

%% Generate 3-Bus Test System
% System bases
SBase = 10;     % Power base (MW)
TBase = 1;      % Temperature base (deg C)

% System info:
Slack = [1];
PV = [2];
PQ = [3];

% Set bus data:
bus.id = [1 2 3];
bus.type = [3 2 0];
bus.G = [0 0 0];
bus.B = [0 0 0];
bus.V_mag = [1.0 0.99 0.9651];
bus.V_angle = [0 -0.03886 -0.04225];
bus.P_load = [0 5 12] ./ SBase;
bus.Q_load = [0 2 6] ./ SBase;
bus.P_gen = [13.415 4 0] ./ SBase;
bus.Q_gen = [0.384 8.317 0] ./ SBase;
bus.P_net = bus.P_gen - bus.P_load;
bus.Q_net = bus.Q_gen - bus.Q_load;

% Set branch data
branch.id = [1 2 3];
branch.type = logical([1 1 1]);
branch.from = [1 1 2];
branch.to = [2 3 3];
branch.R = [0.0300963 0.0314789 0.0243657];
branch.X = [0.0521283 0.0675068 0.0264888];
branch.B = [0 0 0];
branch.tap = [1 1 1];
branch.T_amb = [25 25 25] ./ TBase;
branch.T = branch.T_amb;
branch.R_ref = branch.R;
branch.T_ref = branch.T;
branch.T_f = [228.1 228.1 228.1];
branch.T_rrise = [25 25 45] ./ TBase;
branch.P_rloss = [.0702884 .0925866 .0196559];
branch.R_therm = branch.T_rrise ./ branch.P_rloss;

% Dimensions
N = 3;
M = 1;
R = 3;

% Bus and branch sets for the variables and mismatch equations
sets.P = [2 3];
sets.Q = [3];
sets.H = [1 2 3];
sets.V = [3];
sets.delta = [2 3];
sets.T = [1 2 3];

%% Test 3-bus Test System Admittance Matrix
% Generate admittance matrix
[Y,G,B,~,~] = makeYBus(bus,branch);

%% Test 3-bus Test System Jacobian -- Iteration 1 (from Mathematica)
% Flat-start voltage and angle vectors
% (accounting for voltage-controlled buses)
V = [bus.V_mag(Slack) bus.V_mag(PV) ones(1,M)];
delta = [bus.V_angle(Slack) zeros(1,N-1)];

% Starting temperature info.
T = branch.T;

% Evaluate the full Jacobian matrix for FC-TDPF (flat start)
J = evalJacobian(1,sets,V,delta,T,G,B,branch);

%% Test 3-bus Test System Mismatches -- Iteration 1 (from Mathematica)
% Flat-start states were defined above
mm = evalMismatch(sets,V,delta,T,Y,bus,branch);

%% Test 3-bus Test System Jacobian -- Iteration 2 (from Mathematica)
% Uses Iteration 2 data from Mathematica work
% New states...
V = [V(Slack) V(PV) 0.966193];
delta = [delta(Slack) -0.0369742 -0.0401394];
T = [25.2956 25.0 0.15426];

% New branch resistances...
branch.R = branch.R_ref .* ...
	( (T + branch.T_f) ./ (branch.T_ref + branch.T_f) );

% New YBus...
[Y,G,B,~,~] = makeYBus(bus,branch);

% And new Jacobian!
J = evalJacobian(1,sets,V,delta,T,G,B,branch);

%% Test 3-bus Test System Mismatches -- Iteration 2 (from Mathematica)
% Updated states were defined above
mm = evalMismatch(sets,V,delta,T,Y,bus,branch);

%% Algorithm Tests
% Test, in verbose mode, the 4 algorithms

% Full-Coupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2] = FC_TDPF(bus,branch,'verbose',true);

% Partially-Decoupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2] = PD_TDPF(bus,branch,'verbose',true);

% Fast-Decoupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2] = FD_TDPF(bus,branch,'verbose',true);

% Sequentially-Decoupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2] = SD_TDPF(bus,branch,'verbose',true);


%% MATPOWER Test -- Data Import
% Run an FC-TDPF where temperature dependence of all branches has been
% disabled using the MATPOWER test cases and compare to MATPOWER's
% results.

% NOTES:
%   Strange cases:
%       case300     Seems to have extra generators connected to buses
%                   labeled outside the system?? Doesn't work.
%       case2383wp  FC_TDPF results in a voltage profile about 0.1 higher
%                   overall than the MATPOWER solution!
%       case2736sp  Many branches are status 0, that is, offline. Our
%                   algorithm isn't set up to handle online/offline
%                   branches, so it fails to create YBus correctly.
%                   (This is also true of some of the other large test
%                   systems.)

casename = 'case2383wp';        % Set test case here

% Import MATPOWER test case (needs MATPOWER in the path) and check for
% open lines.
casedata = loadcase(casename);      
if any( casedata.branch(:,11) == 0 )
    warning(['Open lines detected in MATPOWER case data. Our ' ...
             'algorithm can''t handle that yet.']);
end

% Convert casedata to our format
[bus,branch,SBase,TBase] = importCaseData(casename,'MATPOWER');

% For visualization: Scatterplot rated losses vs. resistance values
% loglog(branch.rating, branch.P_rloss,'.')

% Reference YBus
YBusRef = makeYbus(SBase, casedata.bus, casedata.branch);

% Calculated YBus
[Y,G,B,~,~] = makeYBus(bus,branch);

% Difference (should always equal 0 within working precision)
YDiff = YBusRef - Y;
max(abs(YDiff(:)))

%% MATPOWER Test -- Conventional PF
% MATPOWER settings
mpopt = mpoption();
mpopt(2) = 1e-8;        % Set tolerance
mpopt(6) = 0;           % Don't enforce Q-limits on generation
mpopt(31) = 1;          % Turn on verbose mode

% Perform power flows w/ MATPOWER and FC_TDPF
results = runpf(casedata,mpopt);    % MATPOWER
[V,delta] = PF(bus,branch);         % Conventional PF

% Check deviations from MATPOWER values
norm(results.bus(:,8) - V,inf)
norm(results.bus(:,9) - delta,inf)

% Plots of deviations
plot(bus.id, results.bus(:,8) - V)
plot(bus.id, results.bus(:,9) - delta)

% Comparison plots
plot(bus.id, [results.bus(:,8), V])
plot(bus.id, [results.bus(:,9), delta])

% Oddly, these deviations can be very large in some systems even though the
% mismatches come out within tolerance!

% Check mismatches of MATPOWER results vs. our algorithm results, using the
% evalMismatch(). (Should ferret out any major differences in mismatch
% computations.)
sets.P = bus.id((bus.type == 0) | (bus.type == 1) | (bus.type == 2));
sets.Q = bus.id((bus.type == 0) | (bus.type == 1));
sets.H = [];
mm1 = evalMismatch(sets, results.bus(:,8), results.bus(:,9) * (pi/180), ...
        [], Y, bus, branch); % MATPOWER results
mm2 = evalMismatch(sets, V, delta  * (pi/180), ...
        [], Y, bus, branch); % Computed results
norm(mm1,inf)
norm(mm2,inf)
plot([mm1, mm2]);

% For most systems tested, both results are within tolerance

%% MATPOWER Test -- Histories
% Full-Coupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2,hist] = FC_TDPF(bus,branch,'history',true);
maxErrFC = hist.maxErr;

% Partially-Decoupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2,hist] = PD_TDPF(bus,branch,'history',true);
maxErrPD = hist.maxErr;

% Fast-Decoupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2,hist] = FD_TDPF(bus,branch,'history',true);
maxErrFD = hist.maxErr;

% Sequentially-Decoupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2,hist] = SD_TDPF(bus,branch,'history',true);
maxErrSD = hist.maxErr;

% Plot error histories
semilogy( 1:length(maxErrFC), maxErrFC, ...
          1:length(maxErrPD), maxErrPD, ...
          1:length(maxErrFD), maxErrFD, ...
          1:length(maxErrSD), maxErrSD )

%% MATPOWER Test -- Timings
% Full-Coupled Temperature Dependant Power Flow
tic
for tmp = 1:10
    [~] = FC_TDPF(bus,branch);
end
x = toc;
x/10

% Partially-Decoupled Temperature Dependant Power Flow
tic
for tmp = 1:10
    [~] = PD_TDPF(bus,branch);
end
x = toc;
x/10

% Fast-Decoupled Temperature Dependant Power Flow
tic
for tmp = 1:10
    [~] = FD_TDPF(bus,branch);
end
x = toc;
x/10

% Sequentially-Decoupled Temperature Dependant Power Flow
tic
for tmp = 1:10
    [~] = SD_TDPF(bus,branch);
end
x = toc;
x/10

% Put a space in the console
% (Good for comparing times when repeatedly evaluating this block)
disp(' ');
      
      
%% IEEE 30 Bus System Data Import
% Import from file
[bus,branch,SBase,TBase] = importCaseData( ...
    {'IEEE_30_BUS.csv','IEEE_30_BRANCH.csv'}, 'CSV', 'TBase', 100);

% Generate N, L, and PV/PQ/Slack bus sets
N = length(bus.id);
L = length(branch.id);
PQ = sort(bus.id(bus.type == 0 | bus.type == 1));
PV = sort(bus.id( bus.type == 2));
Slack = sort(bus.id( bus.type == 3));

%% IEEE 30 Bus System Test -- Conventional PF
% Same kind of test but with published IEEE 30 bus test case
branch1 = branch; branch1.type = false;
[V,delta,~,bu2,br2,hist] = FC_TDPF(bus,branch1,'history',true);

% Check deviations from published values
norm(bus.V_mag' - V,inf)
norm(bus.V_angle' - delta,inf)

% Plots of deviations
plot(bus.id, bus.V_mag' - V)
plot(bus.id, bus.V_angle' - delta)

% TO DO: Investigate these errors. what happened???
% Things I've checked:
% 1. My sparsity changes to the Jacobian matrix have no effect on this
%    error.

% 2. No weird errors where V and delta aren't matching the values placed
%    back into the 'bus' and 'branch' structures.
norm(bu2.V_mag' - V,inf)
norm(bu2.V_angle' - delta,inf)

% 3. No errors in non-PQ bus voltage magnitudes, so no PV buses have been
%    converted to PQ (gen limit reached).
plot(bus.id(bus.type ~= 0), bus.V_mag(bus.type ~= 0)' - V(bus.type ~= 0))

% 4. Using the evalMismatch() function, we see that the mismatches *appear*
%    better using our computed results than with the published values.
%    Perhaps there is an error in the way mismatches are calculated? Or
%    perhaps the published values are just wrong??
Y = makeYBus(bus,branch);
sets.P = bus.id((bus.type == 0) | (bus.type == 1) | (bus.type == 2));
sets.Q = bus.id((bus.type == 0) | (bus.type == 1));
sets.H = [];
mm1 = evalMismatch(sets, bus.V_mag', bus.V_angle' * (pi/180), ...
        [], Y, bu2, br2); % Published results
mm2 = evalMismatch(sets, V, delta  * (pi/180), ...
        [], Y, bu2, br2); % Computed results
norm(mm1,inf)
norm(mm2,inf)
plot([mm1, mm2]);

% 5. MATPOWER's IEEE 30-bus test case seems to have errors in the published
%    values as well, compared to MATPOWER's computed values. Perhaps it is
%    simply an error in the master data set?
casedata = case_ieee30;
results = runpf(casedata);
plot(casedata.bus(:,1), casedata.bus(:,8) - results.bus(:,8))   % V
plot(casedata.bus(:,1), casedata.bus(:,9) - results.bus(:,9))   % delta

%% IEEE 30 Bus System Test -- Histories
% Full-Coupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2,hist] = FC_TDPF(bus,branch,'history',true);
maxErrFC = hist.maxErr;

% Partially-Decoupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2,hist] = PD_TDPF(bus,branch,'history',true);
maxErrPD = hist.maxErr;

% Fast-Decoupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2,hist] = FD_TDPF(bus,branch,'history',true);
maxErrFD = hist.maxErr;

% Sequentially-Decoupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2,hist] = SD_TDPF(bus,branch,'history',true);
maxErrSD = hist.maxErr;

% Plot error histories
semilogy( 1:length(maxErrFC), maxErrFC, ...
          1:length(maxErrPD), maxErrPD, ...
          1:length(maxErrFD), maxErrFD, ...
          1:length(maxErrSD), maxErrSD )
      
%% IEEE 30 Bus System Test -- Timings
% NOTE: These no longer give much useful info. since the import process
% removes the temperature-dependence of lines in the IEEE 30 bus test
% system.

% Full-Coupled Temperature Dependant Power Flow
tic
for tmp = 1:10
    [~] = FC_TDPF(bus,branch);
end
x = toc;
x/10

% Partially-Decoupled Temperature Dependant Power Flow
tic
for tmp = 1:10
    [~] = PD_TDPF(bus,branch);
end
x = toc;
x/10

% Fast-Decoupled Temperature Dependant Power Flow
tic
for tmp = 1:10
    [~] = FD_TDPF(bus,branch);
end
x = toc;
x/10

% Sequentially-Decoupled Temperature Dependant Power Flow
tic
for tmp = 1:10
    [~] = SD_TDPF(bus,branch);
end
x = toc;
x/10

% Put a space in the console
% (Good for comparing times when repeatedly evaluating this block)
disp('- - -');
