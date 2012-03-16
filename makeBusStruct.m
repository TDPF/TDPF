%% function makeBusStruct: Create bus data structure
% Creates the bus data structure for a power system using a .CSV
% input file with the same column order as the IEEE Common Data Format.
% 
% SYNTAX:
%   [bu,N,PV,PQ,Slack] = makeBusStruct(filename,varargin)
%
% INPUTS:
%   filename =	Filename of CSV file with IEEE test case formatted bus data
%   varargin =  (Optional) Additional arguments passed as name-value pairs
%               (see OPTIONAL INPUTS below)
%
% OPTIONAL INPUTS:
%   The following Optional inputs may be passed as name-value pairs
%   following 'filename':
%
%   'SBase', [val]  Specify the power base of the system. Default = 100 MW.
%
%   Note that the names of these optional arguments are not case sensitive.
%
% OUTPUTS:
%   bu =	Structure containing bus data, as follows:
%               .id         Numerical ID for the bus; will match the .CSV
%                           if the bus data is given in order. (Is always a
%                           vector with entries 1:N.)
%               .type       Bus type: 0 = PQ, 1 = PQ (Gen. limit reached),
%                           2 = PV, 3 = Slack
%               .V_mag      (Final) bus voltage magnitude [pu]
%               .V_angle    (Final) bus voltage angle [deg]
%               .P_load     Load power [pu]
%               .Q_load     Load reactive power [pu]
%               .P_gen      (Final) Generator power [pu] 
%               .Q_gen      (Final) Generator reactive power [pu]
%               .B          Bus shunt conductance [pu]
%               .G          Bus shunt susceptance [pu]
%               .P_net      Net injected real power [pu]
%               .Q_net      Net injected reactive power [pu]
%   N =     Number of buses in system
%   PQ =    Vector indicating the indices of the PQ buses
%   PV =    Vector indicating the indices of the PV buses    
%   Slack = Vector indicating the indices of the slack buses
%
% TO DO:
%   1. Do something special with type 1 buses?
%   2. Import bus ID numbers directly from the file? (I was having trouble
%      with this.)
%

function [bu,N,PV,PQ,Slack] = makeBusStruct(filename,varargin)
    % Import data in IEEE common data format, starting at row 2, column 4
    raw_data = dlmread(filename, ',', 1, 3);
    N = length(raw_data(:,1));  % Number of buses in system
    
    % Assign unique bus ID numbers in order of import
    % NOTE: ID's are in order and therefore can be used as an index
    bu.id = 1:N;
    
    % Set power base [MW]
    SBase = 100;        % Default = 100 MW
    
    % Override default temperature data with any 'varargin' arguments
	while ~isempty(varargin)
		name = lower( varargin{1} );
		switch name
			case {'sbase'}
                SBase = varargin{2};    % Power base [MW]
            otherwise
                warning([ ...
                    'Optional argument ''' varargin{1} ''' is not ' ...
                    'recognized and has therefore been ignored.']);
		end
		
		% Clear these two entries from 'varargin'
		varargin(1:2) = [];
    end
    
    % Assign imported data to bus structure
    % Note: equivalent CSV column = 'raw_data' column + 3
    bu.type = raw_data(:,3);    % Bus type:
                                % 3 = Slack, 2 = PV, 0 = PQ, 
                                %   1 = VAR lim of gen reached
    bu.V_mag = raw_data(:,4);   % (Final) bus voltage magnitude
    bu.V_angle = raw_data(:,5); % (Final) bus voltage angle
    bu.P_load = raw_data(:,6);  % Load power
    bu.Q_load = raw_data(:,7);  % Load reactive power
    bu.P_gen = raw_data(:,8);   % (Final) Generator power
    bu.Q_gen = raw_data(:,9);   % (Final) Generator reactive power
    bu.G = raw_data(:,14);      % Bus shunt conductance
    bu.B = raw_data(:,15);      % Bus shunt susceptance
                                
    % NOTE:
    %   For PV buses:
    %       V_mag is input data
    %       V_angle is (reference) output data
    %       P_load, Q_load are input data
    %       P_gen is input data
    %       Q_gen is (reference) output data
    %       P_net is a fixed quantity
    %   For PQ buses:
    %       V_mag is (reference) output data
    %       V_angle is (reference) output data
    %       P_load, Q_load are input data
    %       P_gen, Q_gen are input data
    %       P_net, Q_net are fixed quantities
    %   For Slack buses:
    %       V_mag is input data
    %       V_angle is input data
    %       P_load, Q_load are input data
    %       P_gen, Q_gen are (reference) output data
    
    % Power is imported in MW; convert to per-unit.
    bu.P_load = bu.P_load ./ SBase;
    bu.Q_load = bu.Q_load ./ SBase;
    bu.P_gen  = bu.P_gen  ./ SBase;
    bu.Q_gen  = bu.Q_gen  ./ SBase;
    
    % Compute net real and reactive power for all buses
    bu.P_net = bu.P_gen - bu.P_load;
    bu.Q_net = bu.Q_gen - bu.Q_load;
    
    % Generate vectors indicating which buses are PV, PQ, Slack
    % NOTE: PQ buses include PV buses w/ Q at the limit (i.e. type 1)
    PQ = sort(bu.id(bu.type == 0 | bu.type == 1));
    PV = sort(bu.id( bu.type == 2));
    Slack = sort(bu.id( bu.type == 3));
    
    % Flip any column vectors to rows...
    name = fieldnames(bu);
    for i = 1:length(name)
        [m,n] = size( bu.(name{i}) );
        if m > n
            bu.(name{i}) = bu.(name{i}).';
        end
    end
end
