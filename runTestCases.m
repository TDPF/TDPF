%% Run Test Cases
% This test includes a representation of the New England transmission
% system.

clear all
close all
clc
%% Script Options
% Choose one of the test cases below

% Case 39
% This test includes a representation of the New England transmission
% system.
casename = 'case39';

% Case 2383wp
% Power flow data for Polish system - winter 1999-2000 peak.
%casename = 'case2383wp';

%
% Options:
%
% Only do temperature calcs on the top 5% of lines (by loading)
tempCalcsMostLoaded = false;


%% Import Case Information
% Import MATPOWER test case (needs MATPOWER in the path) and check for
% open lines.
% 
% To add MATPOWER to your path, do
%   addpath('C:\PATH\TO\MATPOWER');
%   addpath('C:\PATH\TO\MATPOWER\t');
casedata = loadcase(casename);      
if any( casedata.branch(:,11) == 0 )
    error(['Open lines detected in MATPOWER case data. Our ' ...
             'algorithm can''t handle that yet.']);
end

% Convert casedata to our format
[bus,branch,SBase,TBase] = importCaseData(casename,'MATPOWER');

if tempCalcsMostLoaded,
    % Get rid of temp calcs on lower 95% of lines
    
    % Run a standard LF to determine line loadings
    [V,delta,bus2,hist] = PF(bus,branch,'history',true,'lineLoad',true);
    [sortedLoading, iSortedLoading] = sort(hist.lineLoad, 'ascend');
    topIndex = floor(length(sortedLoading)*0.95);
    
    % Make branches not temp dependent
    branch.type(iSortedLoading(1:topIndex)) = 0;    
end

%% Histories
compTable = {'NR'; 'FC-TDPF'; 'PD-TDPF'; 'FD-TDPF'; 'SD-TDPF'};

% Conventional Power Flow
[V,delta,bus2,hist] = PF(bus,branch,'history',true);
maxErrPF = hist.maxErr;
compTable{1,2} = max(hist.iter);

% Full-Coupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2,hist] = FC_TDPF(bus,branch,'history',true);
maxErrFC = hist.maxErr;
compTable{2,2} = max(hist.iter);

% Partially-Decoupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2,hist] = PD_TDPF(bus,branch,'history',true);
maxErrPD = hist.maxErr;
compTable{3,2} = max(hist.iter);

% Fast-Decoupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2,hist] = FD_TDPF(bus,branch,'history',true);
maxErrFD = hist.maxErr;
compTable{4,2} = max(hist.iter);

% Sequentially-Decoupled Temperature Dependant Power Flow
[V,delta,T,bus2,branch2,hist] = SD_TDPF(bus,branch,'history',true);
maxErrSD = hist.maxErr;
compTable{5,2} = [num2str(max(hist.subIter)) ' inner, ' ...
                  num2str(max(hist.iter)) ' outer'];

% Plot error histories
semilogy( 1:length(maxErrPF), maxErrPF, 'k+-', ...
          1:length(maxErrFC), maxErrFC, 'k*-', ...
          1:length(maxErrPD), maxErrPD, 'kd-', ...
          1:length(maxErrFD), maxErrFD, 'ko-', ...
          1:length(maxErrSD), maxErrSD, 'kx-' );
legend('NR', 'FC-TDPF', 'PD-TDPF', 'FD-TDPF', 'SD-TDPF', ...
    'location', 'SouthEast');
title('Convergence History');
xlabel('Iteration');
ylabel('Largest Mismatch (log)');
xlim([1 20]);   % Truncate FD_TDPF for readability

%% Timings
% Conventional Power Flow
tic
for tmp = 1:10
    [trash] = PF(bus,branch);
end
x = toc;
compTable{1,3} = x/10;

% Full-Coupled Temperature Dependant Power Flow
tic
for tmp = 1:10
    [trash] = FC_TDPF(bus,branch);
end
x = toc;
compTable{2,3} = x/10;

% Partially-Decoupled Temperature Dependant Power Flow
tic
for tmp = 1:10
    [trash] = PD_TDPF(bus,branch);
end
x = toc;
compTable{3,3} = x/10;

% Fast-Decoupled Temperature Dependant Power Flow
tic
for tmp = 1:10
    [trash] = FD_TDPF(bus,branch);
end
x = toc;
compTable{4,3} = x/10;

% Sequentially-Decoupled Temperature Dependant Power Flow
tic
for tmp = 1:10
    [trash] = SD_TDPF(bus,branch);
end
x = toc;
compTable{5,3} = x/10;

% Add in a header row
compTable(2:end+1, :) = compTable;
compTable(1,:) = {'Algorithm' 'Iterations' 'Exe Time (sec)'};

disp(['Case "' casename '"'])
if tempCalcsMostLoaded,
    disp('(Only top 5% of most loaded lines calculated)');
end
disp(compTable)






