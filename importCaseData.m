%% WORK IN PROGRESS

%% function importCaseData: Imports bus, branch, and system base data
% Imports system name, system base, bus data, and branch data from one of
% several possible data sources: IEEE Common Data Format text file, CSV
% (2 files), or MATPOWER case data. The source type must be explicitly
% given and match a known source type or the import will fail with an error
% message.
%
% IEEE Common Data Format:
%   The IEEE common data format is a standard way to exchange power flow
%   and optimal power flow data in text format. Selecting this option will
%   import from an IEEE CDF text file as per the column specification. This
%   can sometimes cause problems if the source file isn't perfectly
%   formatted; see the 'COMMENTS' section.
%
% CSV:
%   Comma Separated Values (.CSV) data is assumed to be in two files, one
%   each for the bus and branch data, with the columns given in the order 
%   specified in the IEEE common data format. Column headers will be
%   ignored and the data imported by column number. NOTE: The system base
%   should be specified when using CSV data if it is other than the default
%   value of 100 MVA.
%
% MATPOWER Case Data:
%   Selecting this option will import case data by running the specified
%   filename (casename) and extracting the relevant fields from the
%   results. Thus, the filename should be an executable script which
%   returns data in a MATPOWER case structure; MATPOWER's built-in test
%   cases (installed automatically with MATPOWER) satisfy this requirement.
%   See the MATPOWER documentation for more information.
% 
% SYNTAX:
%   [bus,branch,SBase,TBase] = importCaseData(filename,srcflag,varargin)
%
% INPUTS:
%   filename =	Name of input data file, or, for .CSV data, a cell array
%               containing 2 filenames: the first for bus data and the
%               second for branch data.
%   srcflag =   Text string indicated data source. Must be one of:
%                   'IEEECDF'   for IEEE common data format text data
%                   'CSV'       for .CSV formatted data
%                   'MATPOWER'  for MATPOWER case structure data
%               This input is not case sensitive.
%   varargin =  (Optional) Additional arguments passed as name-value pairs
%               (see OPTIONAL INPUTS below)
%
% OPTIONAL INPUTS:
%   The following optional inputs may be passed as name-value pairs
%   following 'srcflag':
%
%   'SBase', [val]  Specify the power base of the system. Default = 100 MW
%                   or as imported from case data.
%   'brtype', [val] Specify explicitly the type for each branch:
%                   TRUE = temp. dependant, FALSE = not temp. dependant.
%                   Must be either a scalar or a vector of length L = # of
%                   branches. If not specified, the correct type is
%                   inferred from the line parameters.
%   'TBase', [val]  Specify the temperature base of the system.
%                   Default = 100 deg C.
%   'TAmb', [val]   Specifiy the branch ambient temperature(s) (deg C).
%                   Must be either a scalar or a vector of length L = # of
%                   branches. Default = 25 deg C.
%   'TInit', [val]  Specifiy the branch initial temperature(s) (deg C).
%                   Must be either a scalar or a vector of length L = # of
%                   branches. This value is also used to correct R from the
%                   input data to RRef. Default = TAmb.
%   'TRef', [val]   Specifiy the branch reference temperature(s) (deg C).
%                   Must be either a scalar or a vector of length L = # of
%                   branches. Default = TAmb.
%   'Tf', [val]     Specifiy the branch temperature coefficient(s) (deg C).
%                   Must be either a scalar or a vector of length L = # of
%                   branches. Default = 228.1 deg C (aluminum).
%   'TRRise', [val] Specifiy the branch rated temperature rise(s) (deg C).
%                   Must be either a scalar or a vector of length L = # of
%                   branches. Default = 25 deg C over ambient.
%
%   The names of these optional arguments are not case sensitive.
%
%   Any input arguments given in vector form will be matched to the
%   branches in the order of import from the original data file. (IE, there
%   is no internal sorting prior to mapping the data to the branches.)
%
% OUTPUTS:
%   bus =       Structure containing bus data, as follows:
%               .id         Numerical ID for the bus; will match the .CSV
%                           if the bus data is given in order. (Is always a
%                           vector with entries 1:N, where N is the total
%                           number of buses. This ID may not match the
%                           original bus numbers if the input data is not
%                           in order.)
%               .type       Bus type: 0 = PQ, 1 = PQ (Gen. limit reached),
%                           2 = PV, 3 = Slack
%               .V_mag      (Final) bus voltage magnitude [pu]
%               .V_angle    (Final) bus voltage angle [deg]
%               .P_load     Load power [pu]
%               .Q_load     Load reactive power [pu]
%               .P_gen      (Final) Generator power [pu] 
%               .Q_gen      (Final) Generator reactive power [pu]
%               .P_net      Net injected real power [pu]
%               .Q_net      Net injected reactive power [pu]
%               .B          Bus shunt conductance [pu]
%               .G          Bus shunt susceptance [pu]
%   branch =    Structure containing branch element data, as follows:
%               .id         Numerical ID for the branch. (Is always a
%                           vector with entries 1:L, where L is the total
%                           number of branches.)
%               .type       Branch type: TRUE = temp. dependant,
%                           FALSE = not temp. dependant
%               .from       From bus, or tap bus (i in R_ij)
%               .to         To bus, or Z bus (j in R_ij)
%               .R_ref      Reference resistance at reference temperature
%               .R          Initial branch resistance
%               .X          Branch reactance
%               .B          Branch charging suseptance
%               .tap        Complex tap ratio (1 + j0 for nominal)
%               .rating     Power rating of the branch (pu)
%                           (0 indicates no rating present)
%               .T          Initial branch temperature
%               .T_amb      Ambient temperature for the branch
%               .T_ref      Branch reference temperature
%               .T_f        Thermal coefficent based on conductor material 
%               .T_rrise    Rated temperature rise at P_rloss
%               .P_rloss    Loss in conductor at rated temperature rise
%                           (Computed from conductor resistance at rated
%                           temperature rise times MVA rating squared)
%               .R_therm    Thermal resistance for each branch
%	SBase =		Power base of the system (MVA)
%	TBase =		Temperature base of the system (deg C)
%
% COMMENTS:
% 1. For the IEEE CDF column specifications, see:
%       Working Group on a Common Format for the Exchange of Solved Load
%       Flow Data, "Common Data Format for the Exchange of Solved Load Flow
%       Data", IEEE Transactions on Power Apparatus and Systems, Vol.
%       PAS-92, No. 6, November/December 1973, pp. 1916-1925.
% 
%    A partial specification (possibly with a few minor errors) is also
%    available at:
%       http://www.ee.washington.edu/research/pstca/formats/cdf.txt
%
%    CAUTION: Some available IEEE CDF data in text format has slight errors
%    in the column alignment which can cause the import to fail. (In
%    particular, the last two columns of the branch data seem to sometimes
%    divide at 118/119 instead of 119/120.)If errors occur during the
%    import process, try checking the alignment of the columns in the
%    source file. There's a ruler at the bottom of this script for that
%    purpose.
%
% 2. The bus ID's for data imported via this function WILL NOT MATCH the
%    original data if the original buses are not numbered consecutively
%    starting with 1. This is because the TDPF algorithms as written
%    require consecutive bus numbers 1:N for an N bus system (a limitation
%    of the original code). Keep this in mind when making comparisons.
%
% 3. When importing from MATPOWER, the generator voltage setpoint WILL
%    OVERRIDE any other specified bus voltage at the point where the
%    generator connects. This is for consistency with how MATPOWER
%    interprets the case data.
%
% 4. When importing from MATPOWER, if a bus is specified as PV but all
%    connected generators are disabled, the bus is converted to a PQ bus.
%    This is consistent with how MATPOWER treats bus types.
function [bus,branch,SBase,TBase] = importCaseData(filename,srcflag,varargin)
    %% Validation
    % Check for valid source type
    if ~any(strcmpi(srcflag,{'ieeecdf','csv','matpower'}))
        error(['Invalid source type specified. Please specify one of ' ...
               '''IEEECDF'', ''CSV'', or ''MATPOWER''.'])
    end
    
    % For CSV input, check that 'filename' is a cell array of size 2
    if strcmpi(srcflag,'csv') && ...
        (~iscell(filename) || (length(filename) ~= 2))
        error(['For CSV data, ''filename'' must be a cell array ' ...
               'of exactly two filenames as strings.'])
    end
    
    %% Import Data
    % Import the data from the specified location
    srcflag = lower(srcflag);
    switch srcflag
        
case 'ieeecdf'
%% START -- IEEE Common Data Format %%
% Open file
fid = fopen(filename,'r');

% Find boundaries of bus and branch data, designated by '-999' in the
% source file. 
D = textscan(fid,'%4s %*s','whitespace','','delimiter','\n');
D = D{1}; % 1st cell
endlines = find(strcmp(D,'-999'));

% Location of bus data
beginBus = 3;
N = endlines(1) - beginBus;     % Number of buses

% Location of branch data
beginBranch = endlines(1) + 2;
L = endlines(2) - beginBranch;  % Number of branches

% Rewind the file to the beginning
frewind(fid)

% Import general case data
tline = fgetl(fid);

% Extract case name
% CaseName = tline(46:end);	% Columns 46-end of 1st line (unused)

% Extract system base
SBase = str2double(tline(32:37));	% Columns 32-37 of 1st line

% Skip a line
fgetl(fid);

% Scan in BUS data in text format
% (Note: We're now starting on line 3 = line 'beginBus')
format = [ ...  %COLUMNS    FIELD                               FIELD NO.
    '%4s ' ...  %1- 4       Bus number                          1
    '%*1s ' ... %5          SKIP 
    '%12s ' ... %6-17       Name                                2
    '%*1s ' ... %18         SKIP
    '%2s ' ...  %19-20      Load flow area number               3
    '%3s ' ...  %21-23      Loss zone number                    4
    '%*1s ' ... %24         SKIP
    '%2s ' ...  %25-26      Type                                5
    '%*1s ' ... %27         SKIP
    '%6s ' ...  %28-33      Final voltage, p.u.                 6
    '%7s ' ...  %34-40      Final angle, degrees                7
    '%9s ' ...  %41-49      Load MW                             8
    '%10s ' ... %50-59      Load MVAR                           9
    '%8s ' ...  %60-67      Generation MW                       10
    '%8s ' ...  %68-75      Generation MVAR                     11
    '%*1s ' ... %76         SKIP
    '%7s ' ...  %77-83      Base KV                             12
    '%*1s ' ... %84         SKIP
    '%6s ' ...  %85-90      Desired volts (pu)                  13
    '%8s ' ...  %91-98      Maximum MVAR or voltage limit       14
    '%8s ' ...  %99-106     Minimum MVAR or voltage limit       15
    '%8s ' ...  %107-114    Shunt conductance G (pu)            16
    '%8s ' ...  %115-122    Shunt susceptance B (pu)            17
    '%*1s ' ... %123        SKIP
    '%4s ' ];   %124-127	Remote controlled bus number        18
        
busD = textscan(fid,format,N,'whitespace','','delimiter','\n');

% Skip two more lines
fgetl(fid);
fgetl(fid);

% Scan in BRANCH data in text format
% (Note: We're now starting on line 'beginBranch')
format = [ ...  %COLUMNS    FIELD                           FIELD NO.
    '%4s ' ...  %1- 4       Tap bus number                      1
    '%*1s ' ... %5          SKIP
    '%4s ' ...  %6- 9       Z bus number                        2
    '%*1s ' ... %10         SKIP
    '%2s ' ...  %11-12      Load flow area                      3
    '%2s ' ...  %13-14      Loss zone                           4
    '%*2s ' ... %15-16      SKIP (?)
    '%1s ' ...  %17         Circuit                             5
    '%*1s ' ... %18         SKIP
    '%1s ' ...  %19         Type                                6
    '%10s ' ... %20-29      Branch resistance R, per unit       7
    '%11s ' ... %30-40      Branch reactance X, per unit        8
    '%10s ' ... %41-50      Line charging B, per unit           9
    '%5s ' ...  %51-55      Line MVA rating No 1                10
    '%*1s ' ... %56         SKIP
    '%5s ' ...  %57-61      Line MVA rating No 2                11
    '%*1s ' ... %62         SKIP
    '%5s ' ...  %63-67      Line MVA rating No 3                12
    '%*1s ' ... %68         SKIP
    '%4s ' ...  %69-72      Control bus number                  13
    '%*1s ' ... %73         SKIP
    '%1s ' ...  %74         Side                                14
    '%*2s ' ... %75-76      SKIP
    '%6s ' ...  %77-82      Transformer final turns ratio       15
    '%*1s ' ... %83         SKIP
    '%7s ' ...  %84-90      Phase shifter final angle           16
    '%7s ' ...  %91-97      Minimum tap or phase shift          17
    '%7s ' ...  %98-104     Maximum tap or phase shift          18
    '%*1s ' ... %105        SKIP
    '%6s ' ...  %106-111    Step size                           19
    '%*1s ' ... %112        SKIP
    '%7s ' ...  %113-119 	Minimum voltage, MVAR or MW limit   20
    '%7s ' ];   %120-126    Maximum voltage, MVAR or MW limit   21

branchD = textscan(fid,format,L,'whitespace','','delimiter','\n');

% Close the file
fclose(fid);

% Extract bus data from text
bus.id = 1:N;
bus.type = cellfun(@str2double, busD{5});   % Bus type:
                                            % 3 = Slack, 2 = PV, 0 = PQ, 
                                            %   1 = VAR lim of gen reached
bus.V_mag = cellfun(@str2double, busD{6});  % (Final) bus voltage magnitude
bus.V_angle = cellfun(@str2double, busD{7});% (Final) bus voltage angle
bus.P_load = cellfun(@str2double, busD{8}); % Load real power
bus.Q_load = cellfun(@str2double, busD{9}); % Load reactive power
bus.P_gen = cellfun(@str2double, busD{10}); % (Final) Gen. real power
bus.Q_gen = cellfun(@str2double, busD{11}); % (Final) Gen. reactive power
bus.G = cellfun(@str2double, busD{16});     % Bus shunt conductance
bus.B = cellfun(@str2double, busD{17});     % Bus shunt susceptance

% Extract branch data from text
branch.id = 1:L;
branch.from = ...
    cellfun(@str2double, branchD{1});	% Branch from bus (or tap bus)
branch.to = ...
    cellfun(@str2double, branchD{2});	% Branch to bus (or Z bus)
branch.R = ...
    cellfun(@str2double, branchD{7});	% Branch resistance
branch.X = ...
    cellfun(@str2double, branchD{8}); 	% Branch reactance
branch.B = ...
    cellfun(@str2double, branchD{9});	% Branch shunt charging susceptance
branch.tap = ...
    cellfun(@str2double, branchD{15});	% Off-nominal tap ratio (real part)
pshift = ...
    cellfun(@str2double, branchD{16});  % Phase shift [deg]
branch.rating = ...
    cellfun(@str2double, branchD{10});  % Branch normal operation MVA rating
										% (Used to calculate rated loss)
%%% END %%

case 'csv'
%% START -- CSV Files %%
% Retrieve filenames for bus and branch data
[busFN, branchFN] = filename{1:2};

% Import raw bus data (skipping 1st row and 1st 2 columns)
% NOTE: 1st columns are skipped to avoid non-numeric data, which otherwise
% breaks dlmread()
busD = dlmread(busFN, ',', 1, 2);
N = length(busD(:,1));  				% Number of buses in system

% Extract bus data
% NOTE: CSV column is adjusted by 2 for skipped columns
bus.id = 1:N;
bus.type = busD(:,3);					% Bus type:
										% 3 = Slack, 2 = PV, 0 = PQ, 
										%   1 = VAR lim of gen reached
bus.V_mag = busD(:,4);					% (Final) bus voltage magnitude
bus.V_angle = busD(:,5);				% (Final) bus voltage angle
bus.P_load = busD(:,6);					% Load power
bus.Q_load = busD(:,7);					% Load reactive power
bus.P_gen = busD(:,8);					% (Final) Generator power
bus.Q_gen = busD(:,9);					% (Final) Generator reactive power
bus.G = busD(:,14);						% Bus shunt conductance
bus.B = busD(:,15);						% Bus shunt susceptance

% Import raw bus data (skipping 1st row)
branchD = dlmread(branchFN, ',', 1, 0);
L = length(branchD(:,1));  				% Number of branches in system

% Extract branch data
branch.id = 1:L;
branch.from = branchD(:,1);				% Branch from bus (or tap bus)
branch.to = branchD(:,2);				% Branch to bus (or Z bus)
branch.R = branchD(:,7);				% Branch resistance
branch.X = branchD(:,8);				% Branch reactance
branch.B = branchD(:,9);				% Branch shunt charging susceptance
branch.tap = branchD(:,15);				% Off-nominal tap ratio (real part)
pshift = branchD(:,16);					% Phase shift [deg]
branch.rating = branchD(:,10);          % Branch normal operation MVA rating
										% (Used to calculate rated loss)
%%% END %%

case 'matpower'
%% START -- MATPOWER Case Structure %%
% Import case data in MATPOWER format from file, either using MATPOWER's
% loadcase() function (if MATPOWER is installed) or else by running the file
% directly (assuming it will return the proper casedata structure)

% Detect to see if we're getting passed the casedata already (don't have to
% read it in from file)
if isstruct(filename),
    casedata = filename;
else
    % Check for MATPOWER loadcase() function
    if exist('loadcase','file') == 2	
        casedata = loadcase(filename);
    else
        casedata = eval(filename);
    end
end

% Compute bus and branch sizes
N = size(casedata.bus,1);
L = sum( casedata.branch(:,11) > 0);

% Import MVA base
SBase = casedata.baseMVA;

% Import bus data
bus.id = (1:N)';
bus.type = casedata.bus(:,2);
	bus.type( bus.type == 1 ) = 0;  % Map MATPOWER's PQ number to ours
bus.P_load = casedata.bus(:,3);
bus.Q_load = casedata.bus(:,4);
bus.G = casedata.bus(:,5);
bus.B = casedata.bus(:,6);
bus.V_mag = casedata.bus(:,8);
bus.V_angle = casedata.bus(:,9);

% Initialize generator data
bus.P_gen = zeros(N,1);
bus.Q_gen = zeros(N,1);

% Import generator data (also goes into bus data)
% -- Drops offline generators from the data set
% -- Aggregates all generators at each bus
% -- Ignores reactive power limits
ongen = find( casedata.gen(:,8) > 0 );
for g = ongen'                          % <-- Note the transpose
    % New bus ID for this generator
    i = find( casedata.gen(g,1) == casedata.bus(:,1), 1);

    % Real power
    bus.P_gen( i ) = bus.P_gen( i ) + casedata.gen(g,2);

    % Reactive power
    bus.Q_gen( i ) = bus.Q_gen( i ) + casedata.gen(g,3);
    
    % Voltage setpoint
    bus.V_mag( i ) = casedata.gen(g,6);
end

% For off generators, change associated bus type from PV -> PQ
offgen = find( casedata.gen(:,8) == 0 );
for g = offgen'                         % <-- Note the transpose
    % New bus ID for this generator
    i = find( casedata.gen(g,1) == casedata.bus(:,1), 1);
    
    % Break if an enabled generators is connected to this bus
    if ismember(i, casedata.gen(ongen,1))
        continue;
    end

    % Otherwise, swap PV -> PQ
    if bus.type(i) == 2
        bus.type(i) = 0;
    end
end

% Import Branch data
% -- Drops out of service branches
onb = find( casedata.branch(:,11) );    % Only in-service branches

branch.id = (1:L)';
branch.from = casedata.branch(onb,1);	% Branch from bus (or tap bus)
branch.to = casedata.branch(onb,2);		% Branch to bus (or Z bus)
branch.R = casedata.branch(onb,3);		% Branch resistance
branch.X = casedata.branch(onb,4);		% Branch reactance
branch.B = casedata.branch(onb,5);		% Branch shunt charging susceptance
branch.tap = casedata.branch(onb,9);	% Off-nominal tap ratio (real part)
pshift = casedata.branch(onb,10);		% Phase shift [deg]
branch.rating = casedata.branch(onb,6); % Branch normal operation MVA rating
										% (Used to calculate rated loss)
                                        
% Translate old bus ID's to new bus ID's
for l = 1:L
    % Old from bus -> New from bus
    branch.from(l) = find( branch.from(l) == casedata.bus(:,1), 1);
    
    % Old to bus -> New to bus
    branch.to(l) = find( branch.to(l) == casedata.bus(:,1), 1);
end
                                        
%%% END %%
    end
    
    %% Process 'varargin' arguments
    % Now that data has been imported, override any defaults using the
    % user-specified values.
	while ~isempty(varargin)
		name = lower( varargin{1} );
		switch name
			case {'sbase'}
                SBase = varargin{2};    % Power base [MW]
			case {'tbase'}
                TBase = varargin{2};
			case {'brtype'}
                x = varargin{2};
				if length(varargin{2}) == 1
                    branch.type = zeros(1,L);
                    branch.type(1:L) = x;
                elseif length(varargin{2}) == L
                    branch.type = x;
                else
                    warning([ ...
                        'Optional argument ''brtype'' is not a scalar ' ...
                        'or a vector of length L = # of lines ' ...
                        'and has therefore been ignored.']);
                end
			case {'tamb'}
                x = varargin{2};
				if length(varargin{2}) == 1
                    branch.T_amb(1:L) = x;
                elseif length(varargin{2}) == L
                    branch.T_amb = x;
                else
                    warning([ ...
                        'Optional argument ''TAmb'' is not a scalar ' ...
                        'or a vector of length L = # of lines ' ...
                        'and has therefore been ignored.']);
                end
			case {'tinit'}
                x = varargin{2};
				if length(varargin{2}) == 1
                    branch.T(1:L) = x;
                elseif length(varargin{2}) == L
                    branch.T = x;
                else
                    warning([ ...
                        'Optional argument ''TInit'' is not a scalar ' ...
                        'or a vector of length L = # of lines ' ...
                        'and has therefore been ignored.']);
                end
			case {'tref'}
                x = varargin{2};
				if length(varargin{2}) == 1
                    branch.T_ref(1:L) = x;
                elseif length(varargin{2}) == L
                    branch.T_ref = x;
                else
                    warning([ ...
                        'Optional argument ''TRef'' is not a scalar ' ...
                        'or a vector of length L = # of lines ' ...
                        'and has therefore been ignored.']);
                end
			case {'tf'}
                x = varargin{2};
				if length(varargin{2}) == 1
                    branch.T_f(1:L) = x;
                elseif length(varargin{2}) == L
                    branch.T_f = x;
                else
                    warning([ ...
                        'Optional argument ''Tf'' is not a scalar ' ...
                        'or a vector of length L = # of lines ' ...
                        'and has therefore been ignored.']);
                end
			case {'trrise'}
                x = varargin{2};
				if length(varargin{2}) == 1
                    branch.T_rrise(1:L) = x;
                elseif length(varargin{2}) == L
                    branch.T_rrise = x;
                else
                    warning([ ...
                        'Optional argument ''TRRise'' is not a scalar ' ...
                        'or a vector of length L = # of lines ' ...
                        'and has therefore been ignored.']);
                end
            otherwise
                warning([ ...
                    'Optional argument ''' varargin{1} ''' is not ' ...
                    'recognized and has therefore been ignored.']);
		end
		
		% Clear these two entries from 'varargin'
		varargin(1:2) = [];
    end
	
    %% Set Defaults
    % Default values for anything that wasn't imported or user-specified
    
    % Power Base
    if ~exist('SBase','var')
        SBase = 100;
    end
    
    % Temperature Base
    if ~exist('TBase','var')
        TBase = 100;
    end
	
    % Initialize temperature data to default values
    branch.T = 25 * ones(1,L);          % 25 deg C
    branch.T_amb = branch.T;            % TAmb
    branch.T_ref = branch.T;            % TAmb
    branch.T_f = 228.1 * ones(1,L);     % Assume aluminum (Al) conductors
    branch.T_rrise = 25 * ones(1,L);    % 25 deg C over ambient
	
	% Infer the correct type for each branch if 'type' has not been
    % specified explicitly.
    if isfield(branch, 'type')
		branch.type = logical(branch.type);		% Coerce to logical
    else
        % Compute branch type
        branch.type = (branch.R(:) ~= 0) & ...
			(branch.T_rrise(:) ~= 0) & (branch.rating(:) ~= 0);
    end
    
	%% Per-unit conversions
	% Convert bus power data to per-unit using 'SBase'
    bus.P_load = bus.P_load ./ SBase;
    bus.Q_load = bus.Q_load ./ SBase;
    bus.P_gen  = bus.P_gen  ./ SBase;
    bus.Q_gen  = bus.Q_gen  ./ SBase;
    bus.G      = bus.G      ./ SBase;   % MW   -> p.u. admittance
    bus.B      = bus.B      ./ SBase;   % MVar -> p.u. admittance
	
	% Covert temperature data to per-unit using 'TBase'
    branch.T = branch.T ./ TBase ;
    branch.T_amb = branch.T_amb ./ TBase ;
    branch.T_ref = branch.T_ref ./ TBase ;
    branch.T_f = branch.T_f ./ TBase ;
    branch.T_rrise = branch.T_rrise ./ TBase ;
    
    % Convert line ratings using 'SBase'
    branch.rating = branch.rating ./ SBase;
    
    %% Compute Required Quantities
	% Compute net real and reactive power for all buses
    bus.P_net = bus.P_gen - bus.P_load;
    bus.Q_net = bus.Q_gen - bus.Q_load;
	
	% Nominal tap ratios contain zeroes in the 'tap' field; find these
	% and convert to 1.0 (nominal)
	branch.tap( branch.tap == 0 ) = 1;
	
    % Compute complex branch tap ratios
	% NOTE: Nominal phase shift is 0, so no replacements necessary
    branch.tap = branch.tap .* ( cosd(pshift) + 1j * sind(pshift) );
    
    % Compute corrected reference branch resistances based on specified
    % initial temperatures, initial resistances, and reference
    % temperatures.
    branch.R_ref = branch.R(:) .* (branch.T_ref(:) + branch.T_f(:)) ./ ...
        (branch.T(:) + branch.T_f(:));
    
    % Compute rated loss as rated current square times the branch
    % resistance at rated temperature rise.
    % NOTE: Rated current = MVA rating, assuming V=1.0 pu
    Rhot = branch.R_ref(:) .* ...
        (branch.T_ref(:) + branch.T_rrise(:) + branch.T_f(:)) ./ ...
        (branch.T_ref(:) + branch.T_f(:));
    branch.P_rloss = Rhot .* (branch.rating(:)).^2;
    
    % Compute branch thermal resistance
    branch.R_therm = branch.T_rrise(:) ./ branch.P_rloss(:);

	
	%% Data Sorting
	% Transpose any column vectors to rows
    name = fieldnames(bus);
    for i = 1:length(name)
        [m,n] = size( bus.(name{i}) );
        if m > n
            bus.(name{i}) = bus.(name{i}).';
        end
    end
    name = fieldnames(branch);
    for i = 1:length(name)
        [m,n] = size( branch.(name{i}) );
        if m > n
            branch.(name{i}) = branch.(name{i}).';
        end
    end
	
	% Sort fields for consistent ordering
    bus = orderfields(bus, ...
        {'id','type','V_mag','V_angle','P_load','Q_load', ...
         'P_gen','Q_gen','P_net','Q_net','G','B'});
    branch = orderfields(branch, ...
        {'id','type','from','to','R','R_ref','X','B','tap','rating', ...
         'T','T_amb','T_ref','T_f','T_rrise','P_rloss','R_therm'});
end

% Useful line rulers for counting columns
% Note that the 1st '%' counts as a character as part of the ruler
%        10        20        30        40        50        60        70        80        90        100       110       120       130       140       150
%   .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |    .    |
