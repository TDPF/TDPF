%% Run Test Cases
% This test includes a representation of the New England transmission
% system.

clear all
close all
clc
%% Script Options
% Choose one of the test cases below

% Case 30
% IEEE 30 bus test system
casename = 'case30';

% Case 39
% This test includes a representation of the New England transmission
% system.
%casename = 'case39';

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
    [V,delta,bus2,branch2,hist] = ...
        PF(bus,branch,'history',true,'lineLoad',true);
    [sortedLoading, iSortedLoading] = sort(hist.lineLoad, 'ascend');
    topIndex = floor(length(sortedLoading)*0.95);
    
    % Make branches not temp dependent
    branch.type(iSortedLoading(1:topIndex)) = 0;    
end

%% Histories
compTable = {'NR'; 'FC-TDPF'; 'PD-TDPF'; 'FD-TDPF'; 'SD-TDPF'};

% All of the following now properly count iteration 0 as the initial
% conditions and increment thereafter on each update of the state
% variables, whether it is partial or complete. (I changed the algorithm
% code to fix the discrepancies. -SF)

% Conventional Power Flow
[V,delta,bus2,branch2,hist] = PF(bus,branch,'history',true);
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
semilogy( 0:length(maxErrPF)-1, maxErrPF, 'k+-', ...
          0:length(maxErrFC)-1, maxErrFC, 'k*-', ...
          0:length(maxErrPD)-1, maxErrPD, 'kd-', ...
          0:length(maxErrFD)-1, maxErrFD, 'ko-', ...
          0:length(maxErrSD)-1, maxErrSD, 'kx-' );
legend('NR', 'FC-TDPF', 'PD-TDPF', 'FD-TDPF', 'SD-TDPF', ...
    'location', 'SouthEast');
title('Convergence History');
xlabel('Iteration');
ylabel('Largest Mismatch (log)');
xlim([0 20]);   % Truncate FD_TDPF for readability

%% Maximum Errors of Various Sorts
% Compute results for PF vs. TDPF (using FC_TDPF)
[V1,delta1,bus1,branch1] = PF(bus,branch);
[V2,delta2,T2,bus2,branch2] = FC_TDPF(bus,branch);

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

% Compute line loss
V1p = V1 .* exp(1j.*delta1.*pi./180);
V2p = V2 .* exp(1j.*delta2.*pi./180);
lineloss1 = zeros(size(V1));
lineloss2 = zeros(size(V2));
for ik = 1:length(branch.id)
    % Branch 1 loss
    % Get appropriate indices
    i = branch1.from(ik);
    k = branch1.to(ik);
    
    y1 = 1 / (branch1.R(ik) + 1j*branch1.X(ik));
    lineloss1(ik) = real( ...
        V1p(i) * conj(y1 * (V1p(i) - V1p(k))) + ...
        V1p(k) * conj(y1 * (V1p(k) - V1p(i))) ...
        );
    
    % Branch 2 loss
    % Get appropriate indices
    i = branch2.from(ik);
    k = branch2.to(ik);
    
    y2 = 1 / (branch2.R(ik) + 1j*branch2.X(ik));
    lineloss2(ik) = real( ...
        V2p(i) * conj(y2 * (V2p(i) - V2p(k))) + ...
        V2p(k) * conj(y2 * (V2p(k) - V2p(i))) ...
        );
end
lineloss1( lineloss1 < 10*eps ) = 0;
lineloss2( lineloss2 < 10*eps ) = 0;

% Max. difference in line loss
lldiff = lineloss1-lineloss2;
norm(lldiff,Inf)                      	% Absolute
norm(lldiff(abs(lineloss1) > eps) ./ ...% Pct. Relative to conventional PF
    lineloss1(abs(lineloss1) > eps),Inf)*100

% Max. difference in total system loss
sysloss1 = sum(bus1.P_net);
sysloss2 = sum(bus2.P_net);
abs(sysloss1 - sysloss2)              	% Absolute
abs((sysloss1 - sysloss2)/sysloss1)*100 % Pct. Relative to conventional PF

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






