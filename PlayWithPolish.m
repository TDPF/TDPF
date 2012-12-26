%% Play around with Polish case
clear all
close all
clc
% Case 2383wp
% Power flow data for Polish system - winter 1999-2000 peak.
casename = 'case2383wp';
%casename = 'case39';

casedata = loadcase(casename); 

% Get rid of some buses
% 18 is slack bus (31 in case 39)
buses = [1:2383];

% Change slack bus
casedata.bus(31,2) = 1;
casedata.bus(35,2) = 3;

n = size(casedata.bus,1);
rembuses = 1:n;    % The buses to get rid of
rembuses(buses) = [];

% Remove buses
casedata.bus(rembuses,:) = [];

% Remove affected lines
frmBus = casedata.branch(:,1);
toBus = casedata.branch(:,2);
linesToDel = 0;
for i = 1:length(rembuses),
    linesToDel = [linesToDel find(frmBus==rembuses(i))' ...
        find(toBus==rembuses(i))'];
end
linesToDel(1) = [];
casedata.branch(linesToDel,:) = [];

% Remove affected gens
bus = casedata.gen(:,1);
genToDel = 0;
for i = 1:length(rembuses),
    genToDel = [genToDel find(bus==rembuses(i))'];
end
genToDel(1) = [];
casedata.gen(genToDel,:) = [];
casedata.gencost(genToDel,:) = [];
casedata.areas = [];

% Renumber buses as needed
newbranch = casedata.branch;
newgen = casedata.gen;
for i = 1:length(buses),
   % Bus id does not match sequence
   newbranch(casedata.branch(:,1)==buses(i),1) = i;
   newbranch(casedata.branch(:,2)==buses(i),2) = i;
   newgen(casedata.gen(:,1)==buses(i),1) = i;
   casedata.bus(i,1) = i;
end
casedata.branch = newbranch;
casedata.gen = newgen;

matpowerResults = runpf(casedata);
[bus,branch,SBase,TBase] = importCaseData(casedata,'MATPOWER');
[V,delta,bus2,hist] = PF(bus,branch,'history',true,'lineLoad',true);

subplot(2,1,1); plot(matpowerResults.bus(:,8));
subplot(2,1,2); plot(V);

