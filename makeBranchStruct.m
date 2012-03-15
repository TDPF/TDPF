%% function makeBranchStruct: Create branch element data structure
% Creates the branch element data structure for a power system using a .CSV
% input file with the same column order as the IEEE Common Data Format.
% 
% SYNTAX:
%   br = makeBranchStruct(filename)
%
% INPUTS:
%   filename =	Filename of CSV file with IEEE test case formatted branch
%               data
%   varargin =  (Optional) Additional arguments passed as name-value pairs
%               (see OPTIONAL INPUTS below)
%
% OPTIONAL INPUTS:
%   The following Optional inputs may be passed as name-value pairs
%   following 'filename':
%
%   'type', [val]   Specify explicitly the type for each branch:
%                   TRUE = temp. dependant, FALSE = not temp. dependant.
%                   Must be either a scalar or a vector of length L = # of
%                   branches. If not specified, the correct type is
%                   inferred from the line parameters.
%   'TAmb', [val]   Specifiy the ambient temperature for the branches. 
%                   Must be either a scalar or a vector of length L = # of
%                   branches. Default = 25 deg C.
%   'TInit', [val]  Specifiy the initial temperature for the branches.
%                   Must be either a scalar or a vector of length L = # of
%                   branches. This value is also used to correct R from the
%                   IEEE test case data to RRef. Default = TAmb.
%   'TRef', [val]   Specifiy the reference temperature for the branches. 
%                   Must be either a scalar or a vector of length L = # of
%                   branches. Default = TAmb.
%   'Tf', [val]     Specifiy the temperature coefficient for the branches.
%                   Must be either a scalar or a vector of length L = # of
%                   branches. Default = 228.1 deg C (aluminum).
%   'TRRise', [val] Specify the rated temperature rise for each branch.
%                   Must be either a scalar or a vector of length L = # of
%                   branches. Default = 25 deg C over ambient.
%   'TBase', [val]  Specify the temperature base of the system.
%                   Default = 100 deg C.
%
%   Note that the names of these optional arguments are not case sensitive.
%
% OUTPUTS:
%   br =	Structure containing branch element data, as follows:
%               .id         Numerical ID for the branch
%               .type       Branch type: TRUE = temp. dependant,
%                           FALSE = not temp. dependant
%               .from       From bus, or tap bus (i in R_ij)
%               .to         To bus, or Z bus (j in R_ij)
%               .R          Initial branch resistance
%               .X          Branch reactance
%               .B          Branch charging suseptance
%               .tap        Complex tap ratio (1 + j0 for nominal)
%               .T          Initial branch temperature
%               .T_amb      Ambient temperature for the branch
%               .T_ref      Branch reference temperature
%               .T_f        Thermal coefficent based on conductor material 
%               .T_rrise    Rated temperature rise at P_rloss
%               .R_ref      Reference resistance at reference temperature
%               .P_rloss    Loss in condcutor at rated temperature rise
%               .R_therm    Thermal resistance for each branch
%
% TO DO:
%   1. Overrides and/or better computation for rated loss
%   2. Start the branch ID's counting from a user-specified number?
%
function br = makeBranchStruct(filename,varargin)
    % Import data in IEEE common data format
    raw_data = dlmread(filename, ',', 1, 0);
    L = length(raw_data(:,1));  % Number of branches in system
    
    % Assign unique branch ID numbers in order of import
    % NOTE: ID's are in order and therefore can be used as an index
    br.id = 1:L;
    
    % Initialize branch types to an empty array (for now)
    br.type = [];
    
    % Assign imported data to branch element structure
    br.from = raw_data(:, 1)';      % Branch from bus (or tap bus)
    br.to = raw_data(:, 2)';        % Branch to bus (or Z bus)
    br.R = raw_data(:,7)';          % Branch resistance
    br.X = raw_data(:,8)';          % Branch reactance
    br.B = raw_data(:,9)';          % Branch shunt charging susceptance
    br.tap = raw_data(:,15)';       % Off-nominal tap ratio (real part)
                                    % (Imaginary part processed later)
    
    % Initialize temperature data to default values
    br.T = 25 * ones(1,L);          % 25 deg C
    br.T_amb = br.T;                % TAmb
    br.T_ref = br.T;                % TAmb
    br.T_f = 228.1 * ones(1,L);     % Assume aluminum (Al) conductors
    br.T_rrise = 25 * ones(1,L);    % 25 deg C over ambient
    
    % Initialize default temperature base
    TBase = 100;                    % 100 deg C
    
    % Override default temperature data with any 'varargin' arguments
	while ~isempty(varargin)
		name = lower( varargin{1} );
		switch name
			case {'type'}
                x = varargin{2};
				if length(varargin{2}) == 1
                    br.type = zeros(1,L);
                    br.type(1:L) = x;
                elseif length(varargin{2}) == L
                    br.type = x;
                else
                    warning([ ...
                        'Optional argument ''type'' is not a scalar ' ...
                        'or a vector of length L = # of lines ' ...
                        'and has therefore been ignored.']);
                end
			case {'tbase'}
                TBase = varargin{2};
			case {'tamb'}
                x = varargin{2};
				if length(varargin{2}) == 1
                    br.T_amb(1:L) = x;
                elseif length(varargin{2}) == L
                    br.T_amb = x;
                else
                    warning([ ...
                        'Optional argument ''TAmb'' is not a scalar ' ...
                        'or a vector of length L = # of lines ' ...
                        'and has therefore been ignored.']);
                end
			case {'tinit'}
                x = varargin{2};
				if length(varargin{2}) == 1
                    br.T(1:L) = x;
                elseif length(varargin{2}) == L
                    br.T = x;
                else
                    warning([ ...
                        'Optional argument ''TInit'' is not a scalar ' ...
                        'or a vector of length L = # of lines ' ...
                        'and has therefore been ignored.']);
                end
			case {'tref'}
                x = varargin{2};
				if length(varargin{2}) == 1
                    br.T_ref(1:L) = x;
                elseif length(varargin{2}) == L
                    br.T_ref = x;
                else
                    warning([ ...
                        'Optional argument ''TRef'' is not a scalar ' ...
                        'or a vector of length L = # of lines ' ...
                        'and has therefore been ignored.']);
                end
			case {'tf'}
                x = varargin{2};
				if length(varargin{2}) == 1
                    br.T_f(1:L) = x;
                elseif length(varargin{2}) == L
                    br.T_f = x;
                else
                    warning([ ...
                        'Optional argument ''Tf'' is not a scalar ' ...
                        'or a vector of length L = # of lines ' ...
                        'and has therefore been ignored.']);
                end
			case {'trrise'}
                x = varargin{2};
				if length(varargin{2}) == 1
                    br.T_rrise(1:L) = x;
                elseif length(varargin{2}) == L
                    br.T_rrise = x;
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
    
    % Infer the correct type for each branch if 'type' has not been
    % specified explicitly.
    if isempty(br.type)
        br.type = (br.R ~= 0) & (br.T_rrise ~= 0);
    else
        br.type = logical(br.type);
    end
    
    % Assume 0's in the import mean nominal turns ratios (not zero turns
    % ratios)
    br.tap( br.tap == 0 ) = 1;
    
    % Compute complex tap ratios
    pshift = raw_data(:,16)';       % Phase shift [deg]
    br.tap = br.tap .* ( cosd(pshift) + 1j * sind(pshift) );
    
    % Correct all temperature-related data to per-unit based on TBase
    br.T = br.T ./ TBase ;
    br.T_amb = br.T_amb ./ TBase ;
    br.T_ref = br.T_ref ./ TBase ;
    br.T_f = br.T_f ./ TBase ;
    br.T_rrise = br.T_rrise ./ TBase ;
    
    % Compute corrected reference branch resistances based on specified
    % initial temperatures, initial resistances, and reference
    % temperatures.
    br.R_ref = br.R .* (br.T_ref + br.T_f) ./ (br.T + br.T_f);
    
    % P_loss should be proportional to 1/R so this is my best guess
    %  for a function right now
    br.P_rloss = 0.005 ./ br.R;
    
    % Compute thermal resistance
    br.R_therm = br.T_rrise ./ br.P_rloss;
end

    