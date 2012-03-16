%% function PD_TDPF:
% Partially Decoupled Temperature-Dependent Newton-Rhapson Power Flow
% Algorithm. Runs a Newton-Rhapson type power flow augmented with
% temperature data. This is a partially decoupled algorithm with seperate
% updates for temperature and voltage magnitudes + angles (2 Jacobian
% matrices).
% 
% SYNTAX:
%   [V,delta,T,bus,branch,hist] = PD_TDPF(bus,branch,varargin)
%
% INPUTS:
%	bus =       Bus structure (as returned from makeBusStruct())
%	branch =    Branch structure (as returned from makeBranchStruct())
%   varargin =  (Optional) Additional arguments passed as name-value pairs
%               (see OPTIONAL INPUTS below)
%
% OPTIONAL INPUTS:
%   The following Optional inputs may be passed as name-value pairs
%   following 'branch':
%
%   'V', [val]          Specify starting voltage magnitudes. Must be a
%                       scalar or a vector of length N (corresponding to
%                       each bus). If specified as a scalar, all voltage
%                       magnitude variables will be initialized to this
%                       value. If specified as a vector:
%                           (a) Voltages magnitudes at fixed buses will 
%                               override the values from 'V_mag' in 'bus',
%                           (b) Voltages magnitudes at variable buses will
%                               be initialized to the values in this
%                               vector.
%   'delta', [val]      Specify starting voltage angles. Must be a 
%                       scalar or a vector of length N (corresponding to
%                       each bus). If specified as a scalar, all voltage
%                       angle variables will be initialized to this value.
%                       If specified as a vector:
%                           (a) Voltage angles at fixed buses will override 
%                               the values from 'V_angle' in 'bus',
%                           (b) Voltage angles at variable buses will be
%                               initialized to the values in this vector.
%                       Units = degrees
%   'T', [val]          Specify starting temperature vector. Must be a 
%                       scalar or a vector of length L (corresponding to
%                       each branch). If specified as a scalar, all
%                       temperature variables will be initialized to this
%                       value. If specified as a vector:
%                           (a) Resistances corresponding to fixed
%                               branches will be updated using the
%                               temperature values in this vector,
%                           (b) Temperatures at variable branches will be
%                               initialized using the temperature values 
%                               in this vector.
%   'tol', [val]        Specify mismatch tolerance, in per-unit.
%                       Default = 1e-08
%   'maxIter', [val]    Specify maximum number of iterations. Default = 30
%   'verbose', [val]    Enable verbose mode (lots of reporting)? TRUE/FALSE
%                       Default = FALSE
%   'history', [val]    Enable history mode (saves and returns state
%                       variables at each iteration)? TRUE/FALSE
%                       Default = FALSE
%
%   Note that the names of these optional arguments are not case sensitive.
%
% OUTPUTS:
%   V =         Final voltage magnitudes at all buses
%   delta =     Final voltage angles at all buses
%   T =         Final temperatures at all branches
%	bus =       Modified bus structure. Voltage magnitude and angle
%               variables will be set to their final values as determined
%               by the power flow algorithm.
%	branch =    Modified branch structure. Temperature variables and branch
%               resistances will be set to their final values as determined
%               by the power flow algorithm.
%   hist =      History structure of states and mismatches at each
%               iteration. Non-empty only for optional argument
%               'history' == TRUE. State and mismatch matrix columns
%               correspond to the individual iterations in the 'iter'
%               field. NOTE: For consistency with the internal function
%               math, angles in the history are in RADIANS, not degrees.
%
% COMMENTS:
%   1. All the relevant data for determining variables must be contained in
%      the 'type' fields of 'bus' and 'branch'. The function will create
%      the appropriate variable sets automatically. Variable values are
%      initialized to their corresponding values in 'bus' and 'branch', so
%      set these accordingly for a flat start (or whatever).
%   2. If the maximum number of iterations is exceeded, a warning will be
%      issued and the latest values returned.
%   3. Verbose mode turns on all kinds of reporting and is useful for
%      troubleshooting.
%
% TO DO:
%   1. Add ability to create a history matrix
function [V,delta,T,bus,branch,hist] = PD_TDPF(bus,branch,varargin)
    %%Setup	
    % Default control parameters
	tol = 1e-08;		% PQH mismatch tolerance
	maxIter = 30;		% Maximum number of iterations
	verbose = false;	% Set to TRUE to spit out a bunch of diagnostics
    history = false;    % Set to TRUE to save and return iteration history
	
    % Number of buses and branches
    N = length(bus.id);
    L = length(branch.id);
    
	% Parse optional arguments
	while ~isempty(varargin)
		name = lower( varargin{1} );
		switch name
            % Control parameters
			case {'tol'}
				tol = varargin{2};
			case {'maxiter'}
				maxIter = varargin{2};
			case {'verbose'}
				verbose = varargin{2};
			case {'history'}
				history = varargin{2};
            % Initial states, etc.
            case {'v'}
                x = varargin{2};
                if (length(x) == 1) || (length(x) == N)
                    VInput = x;
                else
                    warning([ ...
                        'Optional argument ''V'' is not a scalar ' ...
                        'or a vector of length N = # of buses ' ...
                        'and has therefore been ignored.']);
                end
            case {'delta'}
                x = varargin{2};
                if (length(x) == 1) || (length(x) == N)
                    deltaInput = x;
                else
                    warning([ ...
                        'Optional argument ''delta'' is not a scalar ' ...
                        'or a vector of length N = # of buses ' ...
                        'and has therefore been ignored.']);
                end
            case {'t'}
                x = varargin{2};
				if (length(x) == 1) || (length(x) == L)
                    TInput = x;
                else
                    warning([ ...
                        'Optional argument ''T'' is not a scalar ' ...
                        'or a vector of length L = # of branches ' ...
                        'and has therefore been ignored.']);
                end
		end
		
		% Clear these two entries
		varargin(1:2) = [];
	end
    
    % Storage structure for history (unused if 'history' == FALSE)
    hist = struct;
    
    %% Determine appropriate variable and mismatch sets
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
    
	% Determine more relevant dimensions
    M = length(sets.Q);
    R = length(sets.H);
	
    % Ensure only 1 slack bus
    if (length(sets.P) ~= (N-1))
        error('There can be only 1 slack bus (type == 3).');
    end
    
    %% Initialize States
    % Unknowns:
	%	delta	-> N-1
	%	V       -> M
    %   T       -> R
    
    % Set up Voltage Magnitude, Voltage Angle, Temperature vectors
    % Voltage magnitude
    if exist('VInput','var')
        if length(VInput) == 1
            V = bus.V_mag;
            V(sets.V) = VInput;
        else
            V = VInput;
        end  
    else
        V = bus.V_mag;
        V(sets.V) = 1.0;
    end
    
    % Voltage angle
    if exist('deltaInput','var')
        if length(deltaInput) == 1
            delta = bus.V_angle;
            delta(sets.delta) = deltaInput;
        else
            delta = deltaInput;
        end  
    else
        delta = bus.V_angle;
        delta(sets.delta) = 0.0;
    end
    delta = delta * pi / 180;   % Convert to radians
    
    % Temperature
    if exist('TInput','var')
        if length(TInput) == 1
            T = branch.T;
            T(sets.T) = TInput;
        else
            T = TInput;
        end  
    else
        T = branch.T;
    end
    
    % Ensure column vectors
    if (size(V,2) > size(V,1))              % If # columns > # rows
        V = V';
    end
    if (size(delta,2) > size(delta,1))      % If # columns > # rows
        delta = delta';
    end
    if (size(T,2) > size(T,1))              % If # columns > # rows
        T = T';
    end
    
    % Calculate branch resistances and conductances according to starting
    % temperatures
    branch.T = T';
    branch.R = branch.R_ref .* ...
        ( (branch.T + branch.T_f) ./ (branch.T_ref + branch.T_f) );
    branch.g = branch.R ./ (branch.R.^2 + branch.X.^2);
	
    %% Conventional Power Flow - 1 Iteration
    % Perform a single iteration of conventional power flow in order to
    % improve the initial estimate of line temperatures.
    
    % Initialize YBus
    [Y,G,B,~,~] = makeYBus(bus,branch);
    
    % Evaluate Jacobian
    J = evalJacobian(4,sets,V,delta,[],G,B,branch);

    % Record initial states in history
    if history
        % Iteration -1 = starting conditions
        hist.iter = -1;
        
        % States
        hist.states.delta = delta;
        hist.states.V = V;
        hist.states.T = T;
        
        % Mismatches
        mm = evalMismatch(sets,V,delta,T,Y,bus,branch);
        hist.mismatches.P = mm(1:(N-1));
        hist.mismatches.Q = mm((N-1)+(1:M));
        hist.mismatches.H = mm((N-1+M)+(1:R));
    end
    
    % Evaluate mismatches
    s2 = sets;
    s2.H = [];
    mm = evalMismatch(s2,V,delta,[],Y,bus,branch);
    
    % Set up state vector
    x = [delta(sets.delta); V(sets.V)];
    
    % Perform update
    xnew = x + J \ mm;
    delta(sets.delta) = xnew(1:(N-1));
    V(sets.V)         = xnew((N-1)+(1:M));
    
    %% Initial Temperature Estimate
    % Using the results of the power flow iteration above, compute better
    % guesses for initial starting temperatures.
    for ij = sets.T
        % Get appropriate indices
		i = branch.from(ij);
		j = branch.to(ij);
        
        % Compute power loss in branch
        PLoss = branch.g(ij) * ( V(i)^2 + V(j)^2 - ...
				  2 * V(i) * V(j) * cos(delta(i) - delta(j)) );
              
        % Calculate corresponding temperature rise and update the
        % temperature estimate
        T(ij) = branch.T_amb(ij) + PLoss * branch.R_therm(ij);
    end
        
    %% Power Flow Algorithm
	% Perform iteration until convergence (or max. iterations)
	iter = 0;			% Iteration counter
	while ( iter < maxIter )
		% Display iteration info
        if (verbose)
			disp( ['Iteration: ' int2str(iter)] );
			disp( 'Voltages Magnitudes:' ); disp(V);
			disp( 'Voltages Angles:' ); disp(delta);
            disp( 'Branch Temperatures:' ); disp(T);
        end
        
        % Evaluate YBus
        [Y,G,B,~,~] = makeYBus(bus,branch);
        
        % Evaluate Jacobian matrices
        J = evalJacobian(2,sets,V,delta,T,G,B,branch);
        JPQ = J{1};  % J1 - J4
        JH = J{2};  % J9

        % Evaluate mismatches
        mm = evalMismatch(sets,V,delta,T,Y,bus,branch);
        
        % Display Mismatches
		if (verbose)
			disp('Real Power Mismatch:')
				disp( mm(1:(N-1)) )
			disp('Reactive Power Mismatch:')
				disp( mm((N-1)+(1:M)) )
			disp('Temperature Mismatch:')
				disp( mm((N-1+M)+(1:R)) )
		end
 
        % Record this iteration in history
        if history
            % Iteration
            hist.iter = [hist.iter, iter];

            % States
            hist.states.delta = [hist.states.delta, delta];
            hist.states.V = [hist.states.V, V];
            hist.states.T = [hist.states.T, T];

            % Mismatches
            hist.mismatches.P = [hist.mismatches.P, mm(1:(N-1))];
            hist.mismatches.Q = [hist.mismatches.Q, mm((N-1)+(1:M))];
            hist.mismatches.H = [hist.mismatches.H, mm((N-1+M)+(1:R))];
        end
                
		% Update state vector
		xVd = [delta(sets.delta); V(sets.V)];
        xT = T(sets.T);
        
		% Calculate maximum error / perform update
		err = norm(mm, inf);
        if ( err <= tol )
            if (verbose)
            	disp('PD-TDPF Power Flow Converged.')
            end
            break;
        else
            % Perform update
            xVdnew = xVd + JPQ \ mm(1:(N+M-1));   	% V and delta
            xTnew = xT + JH \ mm((N+M-1)+(1:R));	% T 
            delta(sets.delta) = xVdnew(1:(N-1));
            V(sets.V)         = xVdnew((N-1)+(1:M));
            T(sets.T)         = xTnew;
        end

        % Update resistances using new temperatures
        % (Only needs to occur for temperature-dependant branches)
        branch.T(sets.T) = T(sets.T);
        branch.R(sets.T) = branch.R_ref(sets.T) .* ...
            ( (branch.T(sets.T) + branch.T_f(sets.T)) ./ ...
            (branch.T_ref(sets.T) + branch.T_f(sets.T)) );
        branch.g(sets.T) = branch.R(sets.T) ./ ...
            (branch.R(sets.T).^2 + branch.X(sets.T).^2);
        
		% Increment iteration
		iter = iter + 1;
		if (verbose), disp(''), end
      
        % Create history matrix
        % TO DO: Get back to this
	end	% End While
    
    % Warn if maximum iterations exceeded
    if iter >= maxIter
        warning(['Maximum number of power flow iterations exceeded ' ...
                 'prior to convergence within specified tolerance.']);
    end
    
    %% Parse Results of Power Flow
	% Convert delta back to degrees
	delta = delta * 180 / pi;
    
    % Store final results back to 'bus' and 'branch'
    bus.V_mag = V;
    bus.V_angle = delta;
    branch.T = T';
    
    % Compute maximum error at each iteration and store to history
    if history
        hist.maxErr = max( [ ...
                        max( abs( hist.mismatches.P ) ); ...
                        max( abs( hist.mismatches.Q ) ); ...
                        max( abs( hist.mismatches.H ) ) ...
                        ] );
    end
    
end	% End Function