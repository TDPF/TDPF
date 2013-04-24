%% Load NY ISO test system
% Clear
clear; clc;

% From file
% (Download the file from http://www.pserc.cornell.edu/tcc/tcc.md?SYSTEM=3)
load('nyiso.mat');

%% Setup MATPOWER structure
% Number of buses, branches
N = length(Bus.NUMBER);
L = length(Branch.I);

% General fields
d.version = '2';
d.baseMVA = 100;        % Assumed!

%% Buses
% Populate field by field
% Things marked with a (?) are either uncertain or are made up data

% Bus number (?)
%d.bus(:,1) = Bus.NUMBER(:);        
d.bus(:,1) = 1:length(Bus.NUMBER);  % Changed to new numbers straight off

% Bus type
d.bus(:,2) = Bus.type;

% Load real power
d.bus(:,3) = real( Bus.Spuload ) * d.baseMVA;

% Load reactive power
d.bus(:,4) = imag( Bus.Spuload ) * d.baseMVA;

% Shunt conductance (?)
d.bus(:,5) = 0;

% Shunt susceptance (?)
d.bus(:,6) = 0;

% Area
d.bus(:,7) = Bus.AREA;

% Voltage magnitude
d.bus(:,8) = abs( Bus.V );

% Voltage angle
d.bus(:,9) = angle( Bus.V ) * 180 / pi;

% Min. voltage
d.bus(:,12) = Bus.Vlow;

% Max. voltage (?)
d.bus(:,13) = 1.10;

%% Generators
% What buses have generators?
g = find( (Bus.Pgmax ~= 0) | (Bus.Pgmin ~= 0) );

% Populate field by field
% Things marked with a (?) are either uncertain or are made up data

% Generator bus
d.gen(:,1) = d.bus(g,1);

% Real power
d.gen(:,2) = real( Bus.Spugen(g) ) * d.baseMVA;

% Reactive power
d.gen(:,3) = imag( Bus.Spugen(g) ) * d.baseMVA;

% Max. reactive power
d.gen(:,4) = Bus.Qgmax(g) * d.baseMVA;

% Min. reactive power
d.gen(:,5) = Bus.Qgmin(g) * d.baseMVA;

% Voltage setpoint (?)
d.gen(:,6) = abs( Bus.V(g) );

% MVA base
d.gen(:,7) = d.baseMVA;

% Generator status (?)
d.gen(:,8) = 1;

% Max. real power
d.gen(:,9) = Bus.Pgmax(g) * d.baseMVA;

% Min. real power
d.gen(:,10) = Bus.Pgmin(g) * d.baseMVA;

%% Branches
% Populate field by field
% Things marked with a (?) are either uncertain or are made up data

% From bus (?)
d.branch(:,1) = Branch.I;

% To bus (?)
d.branch(:,2) = Branch.J;

% Resistance
d.branch(:,3) = real( Branch.Z );

% Reactance
d.branch(:,4) = imag( Branch.Z );

% Line charging
d.branch(:,5) = real( Branch.B );

% Line ratings
d.branch(:,6:8) = Branch.RATES(:, 1:3);

% Tap magnitude
d.branch(:,9) = abs( Branch.TAP );

% Phase shift
d.branch(:,10) = angle( Branch.TAP ) * 180 / pi;

% Branch status (?)
d.branch(:,11) = 1;

% %% Renumber
% % MATPOWER has problems with the current way buses are numbered in this
% % data set, primarily because it skips numbers in the sequence. This rather
% % tedious set of loops performs a full renumbering of all bus information.
% 
% % Old -> New
% oldBus = d.bus(:,1);
% newBus = (1:size(d.bus, 1))';
% 
% % Loop and replace
% for i = 1:length(oldBus)
%     d.bus( d.bus(:,1) == oldBus(i), 1 ) = -newBus(i);      	% Buses
%     d.gen( d.gen(:,1) == oldBus(i), 1 ) = -newBus(i);      	% Generators
% end
% 
% % Undo negation
% d.bus(:,1) = -d.bus(:,1);
% d.gen(:,1) = -d.gen(:,1);
% 
% % NOTE: For some bizarre reason, branches already seem to be indexed to bus
% % positions in the vector, rather than numbers?? Therefore, they do not
% % need correction.

%% Save 
% Change name
nyiso = d;

% Save
save('nyiso_converted', 'nyiso');

%% Test
% Attempt to run a power flow on this data using MATPOWER
results = runpf(nyiso);

% The power flow runs, but terminates unsucessfully






