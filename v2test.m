%% Scratchpad: Tests version 2 of the TDPF algorithms
%% Setup
clear all;
close all;
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
%   Normal cases:
%       case39          New England 39-bus test system
%       case30          Derived from IEEE 30-bus test system
%       case2383wp      Polish 1999-2000 winter peak
%       case2736sp      Polish 2004 summer peak
%
%   Strange cases:
%       case3012wp      MATPOWER refuses to load this data (throws error)
%
%       case3375wp      Generates dimension mismatch errors when trying to
%                       compare MATLAB's solution with ours. (??)

casename = 'case30';        % Set test case here

% Case NYISO - outside data source
% load('test_cases/nyiso_converted.mat');
% casename = nyiso;

% Import MATPOWER test case (needs MATPOWER in the path).
% 
% To add MATPOWER to your path, do
%   addpath('C:\PATH\TO\MATPOWER');
casedata = loadcase(casename);

% Convert casedata to our format
[bus,branch,SBase,~] = importCaseData(casename,'MATPOWER');

% For visualization: Scatterplot rated losses vs. resistance values
% loglog(branch.rating, branch.P_rloss,'.')

% Reference YBus
casedata2 = ext2int(casedata);
YBusRef = makeYbus(SBase, casedata.bus, casedata.branch);
clear('casedata2');

% Calculated YBus
[Y,~,~,~,~] = makeYBus(bus,branch);

% Difference (should always equal 0 within working precision)
YDiff = YBusRef - Y;
max(abs(YDiff(:)))

%% MATPOWER Test -- Conventional PF
% MATPOWER settings
mpopt = mpoption();
mpopt(2) = 1e-8;        % Set tolerance
mpopt(6) = 0;           % Don't enforce Q-limits on generation
mpopt(31) = 1;          % Enable verbose mode

% Perform power flows w/ MATPOWER and FC_TDPF
results = runpf(casedata,mpopt);                    % MATPOWER
[V,d] = PF(bus,branch);                             % Conventional PF
[V2,d2] = PF(bus,branch,'FD',true,'maxIter',1000);	% Fast-decoupled PF

% Check deviations from MATPOWER values
norm(results.bus(:,8) - V,inf)
norm(results.bus(:,9) - d,inf)

% Plots of deviations
plot(bus.id, results.bus(:,8) - V)
plot(bus.id, results.bus(:,9) - d)

% Comparison plots
plot(bus.id, [results.bus(:,8), V, V2]);
legend('MATPOWER','Conventional PF','FDPF');

plot(bus.id, [results.bus(:,9), d, d2])
legend('MATPOWER','Conventional PF','FDPF');

% Oddly, these deviations can be very large in some systems even though the
% mismatches come out within tolerance!

% Check mismatches of MATPOWER results vs. our algorithm results, using the
% evalMismatch(). (Should ferret out any major differences in mismatch
% computations.)
sets.P = bus.id((bus.type == 0) | (bus.type == 1) | (bus.type == 2));
sets.Q = bus.id((bus.type == 0) | (bus.type == 1));
sets.H = [];
mm0 = evalMismatch(sets, bus.V_mag, bus.V_angle * (pi/180), ...
        [], Y, bus, branch); % Original case data
mm1 = evalMismatch(sets, results.bus(:,8), results.bus(:,9) * (pi/180), ...
        [], Y, bus, branch); % MATPOWER results
mm2 = evalMismatch(sets, V, d  * (pi/180), ...
        [], Y, bus, branch); % Computed results - Conventional PF
mm3 = evalMismatch(sets, V2, d2  * (pi/180), ...
        [], Y, bus, branch); % Computed results - Fast-decoupled PF
norm(mm0,inf)
norm(mm1,inf)
norm(mm2,inf)
norm(mm3,inf)

% Plot magnitude of mismatches
semilogy(abs([mm0, mm1, mm2, mm3]));
legend('Original Data','MATPOWER','Conventional PF','FDPF');

% For most systems tested, MATPOWER and conventional PF are within
% tolerance. (Fast-decoupled PF often times out.)

%% MATPOWER Test -- Histories
% Initialize to solved power flow
V_init = V;
delta_init = d;

% Full-Coupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2,hist] = FC_TDPF(bus, branch, ...
    'history', true, 'V', V_init, 'delta', delta_init );
maxErrFC = hist.maxErr;

% Partially-Decoupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2,hist] = PD_TDPF(bus, branch, ...
    'history', true, 'V', V_init, 'delta', delta_init );
maxErrPD = hist.maxErr;

% Fast-Decoupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2,hist] = FD_TDPF(bus, branch, ...
    'history', true, 'V', V_init, 'delta', delta_init );
maxErrFD = hist.maxErr;

% Sequentially-Decoupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2,hist] = SD_TDPF(bus, branch, ...
    'history', true, 'V', V_init, 'delta', delta_init );
maxErrSD = hist.maxErr;

% Plot error histories
semilogy( 1:length(maxErrFC), maxErrFC, ...
          1:length(maxErrPD), maxErrPD, ...
          1:length(maxErrFD), maxErrFD, ...
          1:length(maxErrSD), maxErrSD )
legend('FC','PD','FD','SD');

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

%% Test of Timing Features
% Run manually - test timings for each type of power flow
if false
    clc;
    [~] = PF(bus,branch,'timing',true);
    [~] = PF(bus,branch,'FD',true,'maxIter',1000,'timing',true);
    [~] = FC_TDPF(bus,branch,'timing',true);
    [~] = PD_TDPF(bus,branch,'timing',true);
    [~] = FD_TDPF(bus,branch,'timing',true);
    [~] = SD_TDPF(bus,branch,'timing',true);
end
      
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

%% Jacobian Matrix Test -- Off-nominal + Phase Shifting Branches
% Tests the Jacobian matrix elements for dP/dT and dQ/dT for small changes
% in temperature for off-nominal and phase-shifting branches. (The math for
% these branches is very difficult and I'm trying to verify it
% numerically.)

% Case data import
% casename = 'case39';
casename = 'case2383wp';

[bus,branch,SBase,TBase] = importCaseData(casename,'MATPOWER');

% Sets
    % Power mismatch -> PQ, PV buses
    sets.P = bus.id((bus.type == 0) | (bus.type == 1) | (bus.type == 2));
    % Reactive power mismatch -> PQ buses
    sets.Q = bus.id((bus.type == 0) | (bus.type == 1));
    % Temperature mismatch -> Temp. dependant lines
    sets.H = branch.id((branch.type == true));
    % Voltage angle variables -> PQ, PV buses
    sets.delta = sets.P;
    % Voltage magnitude variables -> PQ buses
    sets.V = sets.Q;
    % Temperature variables ->  Temp. dependant lines
    sets.T = sets.H;
    
% Sizes
    N = length(bus.id);
    M = length(sets.Q);
    L = length(sets.H);

% % Find a branch with off-nominal tap
% find((abs(branch.tap) ~= 1) & (branch.type ~= 0))	% Let's pick 33
% kn = find(branch.id == 33);
    
% Find a branch with off-nominal phase shift
find((imag(branch.tap) ~= 0) & (branch.type ~= 0))	% Let's pick 374
kn = find(branch.id == 374);

abs( branch.tap(kn) )                    % Display tap magnitude
angle( branch.tap(kn) ) * 180 / pi       % Display tap voltage angle [deg]

% Indices corresponding to this branch
k = branch.from(kn);
n = branch.to(kn);

% Rows and columns of Jacobian matrix to check
rowPf = find( sets.P == k );
rowPt = find( sets.P == n );
col = find( sets.H == kn) + (N+M-1);


% Solve this test case with sequentially decoupled TDPF
[V,delta,T,bus,branch,hist] = ...
    PD_TDPF(bus,branch,'history',true,'maxIter',100);
maxErrFD = hist.maxErr;

% Find Admittance Matrix for solved power flow
[Y,G,B,~,~] = makeYBus(bus,branch);

% Find sensitivities (1st derivatives)
J = evalJacobian(1, sets, V, delta * pi/180, T, G, B, branch);


% Check these entries in the Jacobian
disp( J(rowPf, col) );
disp( J(rowPt, col) );

% Evaluate the as-solved solved mismatches
mm = evalMismatch(sets, V, delta * pi/180, T, Y, bus,branch);

% Make a small change to T on that branch
deltaT = 0.01;
branch.T(kn) = branch.T(kn) + deltaT;
branch.R(kn) = branch.R_ref(kn) .* ...
    ( (branch.T(kn) + branch.T_f(kn)) ./ ...
    (branch.T_ref(kn) + branch.T_f(kn)) );
branch.g(kn) = branch.R(kn) ./ ...
    (branch.R(kn).^2 + branch.X(kn).^2);
branch.b(kn) = -branch.X(kn) ./ ...
    (branch.R(kn).^2 + branch.X(kn).^2);
[Y2,~,~,~,~] = makeYBus(bus,branch);

% Evaluate the new mismatches
mm2 = evalMismatch(sets,V,delta * pi/180,T,Y2,bus,branch);

% Check the difference at the branches of interest
disp('Actual P sensitivity, From/To')
disp( mm(rowPf) - mm2(rowPf) );
disp( mm(rowPt) - mm2(rowPt) );

disp('Computed P sensitivity, From/To')
disp( J(rowPf, col) * deltaT);
disp( J(rowPt, col) * deltaT);

%%
% Actual P sensitivity, From/To
%  -5.4567e-005
% 
%   7.0971e-005
% 
% Computed P sensitivity, From/To
%    (1,1)     1.6978e-004
% 
%    (1,1)    -1.5479e-004

