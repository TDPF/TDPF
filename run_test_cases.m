%% Run Test Cases
% This script executes the test cases presented in the paper.
clear all; close all; clc

% NOTE: MATPOWER must be installed and in MATLAB's path.
% To add MATPOWER to your path, execute the MATLAB command
% >> addpath('C:\PATH\TO\MATPOWER');

%% Select Test Case
% Choose one of the MATPOWER test cases below.

% Case 30 - IEEE 30 bus test system
%casename = 'case30';

% Case 39 - A small representation of the New England transmission system
casename = 'case39';

% Case 2383wp - Power flow data for Polish system: winter 1999-2000 peak
%casename = 'case2383wp';

% Case 2736sp - Power flow data for Polish system: summer 2004 peak
%casename = 'case2736sp';

%% Script Options
% Only do temperature calcs on the top 5% of lines (by loading)
tempCalcsMostLoaded = false;

% Initialize TDPF to solved PF - use for larger systems
initTDPFtoPF = true;

%% Import Case Information
% Import MATPOWER test case and check for open lines.
casedata = loadcase(casename);      

% Convert casedata to TDPF format
[bus, branch, SBase, TBase] = importCaseData(casename,'MATPOWER');

%% More Setup
% Initialize to solved PF if applicable
if initTDPFtoPF
    [~, ~, bus, branch] = PF(bus, branch);
end

% If applicable, subset temperature calculations to most loaded lines
if tempCalcsMostLoaded
    % Gets rid of temp calcs on lower 95% of lines that have finite
    % thermal resistance...
    
    % Run a standard LF to determine line loadings
    [~, ~, ~, branch2, hist] = ...
        PF(bus, branch, 'history', true, 'lineLoad', true);
    [~, iSortedLoading] = sort(hist.lineLoad, 'ascend');
    
    % Filter to lines with finite thermal resistance
    iSortedLoading = ...
        iSortedLoading( isfinite(branch2.R_therm(iSortedLoading)) );
    topIndex = floor(length(iSortedLoading) * 0.95);
    
    % Make branches not temp dependent
    branch.type(iSortedLoading(1:topIndex)) = 0;
end

% Flat start for smaller systems
if any(strcmp(casename, {'case30','case39'}))
    bus.V_mag(bus.type ~= 2 & bus.type ~= 3) = 1;
    bus.V_angle(bus.type ~= 3) = 1;
end

%% Histories
% Comparison table
compTable = {'PF'; 'FD-PF'; 'FC-TDPF'; 'PD-TDPF'; 'FD-TDPF'; 'SD-TDPF'};

% History cell array
history = cell(6, 1);

% Conventional Power Flow
[~, ~, ~, ~, history{1}] = PF(bus, branch, 'history', true);

% Fast Decoupled Power Flow
[~, ~, ~, ~, history{2}] = PF(bus, branch, 'history', true, ...
    'FD', true, 'maxIter', 100);

% Full-Coupled Temperature Dependant Power Flow
[~, ~, ~, ~, ~, history{3}] = FC_TDPF(bus, branch, 'history', true);

% Partially-Decoupled Temperature Dependant Power Flow
[~, ~, ~, ~, ~, history{4}] = PD_TDPF(bus, branch, 'history', true);

% Fast-Decoupled Temperature Dependant Power Flow
[~, ~, ~, ~, ~, history{5}] = FD_TDPF(bus, branch, 'history', true);

% Sequentially-Decoupled Temperature Dependant Power Flow
[~, ~, ~, ~, ~, history{6}] = SD_TDPF(bus, branch, 'history', true);

% Parse # iterations
for i = 1:5
    compTable{i,2} = max(history{i}.iter);
end
compTable{6,2} = [num2str(length(history{6}.subIter)-1) ' inner, ' ...
                  num2str(max(history{6}.iter)) ' outer'];

% Plot error histories
semilogy( 0:length(history{1}.maxErr)-1, history{1}.maxErr, 'k+-', ...
          0:length(history{2}.maxErr)-1, history{2}.maxErr, 'k.-', ...
          0:length(history{3}.maxErr)-1, history{3}.maxErr, 'k*-', ...
          0:length(history{4}.maxErr)-1, history{4}.maxErr, 'kd-', ...
          0:length(history{5}.maxErr)-1, history{5}.maxErr, 'ko-', ...
          0:length(history{6}.maxErr)-1, history{6}.maxErr, 'kx-' );
legend('Conventional PF', 'Fast Dec. PF', ...
    'FC-TDPF', 'PD-TDPF', 'FD-TDPF', 'SD-TDPF', ...
    'location', 'SouthEast');
title('Convergence History');
xlabel('Iteration');
ylabel('Largest Mismatch (log)');
xlim([0 20]);   % Truncate FD_TDPF for readability


%% Save Results
% Maximum length
maxlen = 1;
for i = 1:6
    maxlen = max(maxlen, length(history{i}.maxErr) );
end

% Create matrix
allerr = NaN(maxlen, 6);

% Fill matrix
for i = 1:6
    % Get history
    x = history{i}.maxErr(:);
    
    % Store it
    allerr(1:length(x), i) = x;
end

% Save to file
% Order is:
%	Conventional PF
%   Fast Decoupled PF
%	FC-TDPF
%	PD-TDPF
%	FD-TDPF
%	SD-TDPF
if tempCalcsMostLoaded
    csvwrite([casename ' - 95th Pct Lines.csv'], allerr)
else
    csvwrite([casename ' - All Lines.csv'], allerr)
end

%% Display Maximum Errors of Various Sorts
% Compute results for PF vs. TDPF (using FC_TDPF)
[V1, delta1, bus1, branch1, hist] = ...
    PF(bus, branch, 'history', true, 'lineLoad', true);
lineLoad = hist.lineLoad;

[V2, delta2, T2, bus2, branch2] = ...
    FC_TDPF(bus, branch);

% Max. difference in voltage magnitude
Vdiff = V1-V2;
norm(Vdiff,Inf)                         % Absolute
norm(Vdiff(abs(V1) > eps) ./ ...        % Pct. Relative to conventional PF
    V1(abs(V1) > eps),Inf)*100          

% Max. difference in voltage angle
ddiff = delta1-delta2;
norm(ddiff,Inf)                       	% Absolute
norm(ddiff(abs(delta1) > eps) ./ ...    % Pct. Relative to conventional PF
    delta1(abs(delta1) > eps),Inf)*100  

% Max. difference in branch resistance
Rdiff = branch1.R-branch2.R;
norm(Rdiff,Inf)                         % Absolute
norm(Rdiff(abs(branch1.R) > eps) ./ ...	% Pct. Relative to conventional PF
    branch1.R(abs(branch1.R) > eps),Inf)*100 

% Compute branch losses
branch1 = evalBranchLoss(bus1,branch1);
branch2 = evalBranchLoss(bus2,branch2);

% Extract branch losses
branchloss1 = real(branch1.S_loss);
branchloss2 = real(branch2.S_loss);
    
% Near zero -> zero
branchloss1( branchloss1 < 10*eps ) = 0;
branchloss2( branchloss2 < 10*eps ) = 0;

% Max. difference in line loss: only for lines loaded above 5%
ll1 = branchloss1( lineLoad > 0.05 );
ll2 = branchloss2( lineLoad > 0.05 );
lldiff = ll1-ll2;
norm(lldiff,Inf)                                % Absolute
norm(lldiff(ll1 > 0) ./ ll1(ll1 > 0),Inf)*100   % Pct. Relative to conventional PF

% Max. difference in total system loss
sysloss1 = sum(bus1.P_net);
sysloss2 = sum(bus2.P_net);
abs(sysloss1 - sysloss2)              	% Absolute
abs((sysloss1 - sysloss2)/sysloss1)*100 % Pct. Relative to conventional PF

% Plot relative difference in line loss vs. line loading
[~, ii] = sort(lineLoad(lineLoad > 0.05));
lldiffrel = (ll2 - ll1)./(ll1);
plot( lldiffrel(ii) );

%% Timings
% Select number of runs
if any(strcmp(casename, {'case30','case39'}))
    numRuns = 25;
else
    numRuns = 5;
end

% Conventional Power Flow
tic
for tmp = 1:numRuns
    [~] = PF(bus,branch);
end
x = toc;
compTable{1,3} = x/numRuns;

% Fast Decoupled Power Flow
tic
for tmp = 1:numRuns
    [~] = PF(bus,branch,'FD',true,'maxIter',100);
end
x = toc;
compTable{2,3} = x/numRuns;

% Full-Coupled Temperature Dependant Power Flow
tic
for tmp = 1:numRuns
    [~] = FC_TDPF(bus,branch);
end
x = toc;
compTable{3,3} = x/numRuns;

% Partially-Decoupled Temperature Dependant Power Flow
tic
for tmp = 1:numRuns
    [~] = PD_TDPF(bus,branch);
end
x = toc;
compTable{4,3} = x/numRuns;

% Fast-Decoupled Temperature Dependant Power Flow
tic
for tmp = 1:numRuns
    [~] = FD_TDPF(bus,branch);
end
x = toc;
compTable{5,3} = x/numRuns;

% Sequentially-Decoupled Temperature Dependant Power Flow
tic
for tmp = 1:numRuns
    [~] = SD_TDPF(bus,branch);
end
x = toc;
compTable{6,3} = x/numRuns;

% Compute timings relative to conventional PF
for i = 1:6
    compTable{i,4} = compTable{i,3} / compTable{1,3};
end


%% Display Timings
% Add in a header row
compTable(2:end+1, :) = compTable;
compTable(1,:) = {'Algorithm' 'Iterations' 'Exe Time (sec)' ...
                  'Relative to PF'};

% Print
disp(['Case "' casename '"'])
if tempCalcsMostLoaded,
    disp('(Only top 5% of most loaded lines calculated)');
end
disp(compTable)






