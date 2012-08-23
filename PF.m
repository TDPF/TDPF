%% function PF:
% Runs a conventional Newton-Rhapson type power flow.
% 
% SYNTAX:
%   [V,delta,bus,hist] = PF(bus,branch,varargin)
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
%   'tol', [val]        Specify mismatch tolerance, in per-unit.
%                       Default = 1e-08
%   'maxIter', [val]    Specify maximum number of iterations. Default = 10
%   'verbose', [val]    Enable verbose mode (lots of reporting)? TRUE/FALSE
%                       Default = FALSE
%   'history', [val]    Enable history mode (saves and returns state
%                       variables at each iteration)? TRUE/FALSE
%                       Default = FALSE
%   'lineLoad', [val]   When history mode is enabled, saves the
%                       line loadings (as a percent of line rating) into
%                       the history structure.
%                       Default = FALSE
%
%   Note that the names of these optional arguments are not case sensitive.
%
% OUTPUTS:
%   V =         Final voltage magnitudes at all buses
%   delta =     Final voltage angles at all buses
%	bus =       Modified bus structure. Voltage magnitude and angle
%               variables will be set to their final values as determined
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
function [V,delta,bus,hist] = PF(bus,branch,varargin)
    %%Setup	
    % Default control parameters
	tol = 1e-08;		% PQH mismatch tolerance
	maxIter = 30;		% Maximum number of iterations
	verbose = false;	% Set to TRUE to spit out a bunch of diagnostics
    history = false;    % Set to TRUE to save and return iteration history
    calcLineLoadings = false;  % Set to TRUE to return line loadings
	
    % Number of buses and branches
    N = length(bus.id);
    
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
            case {'lineload'}
                calcLineLoadings = varargin{2};
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
		end
		
		% Clear these two entries
		varargin(1:2) = [];
	end
    
    % Storage structure for history (unused if 'history' == FALSE)
    hist = struct;
    if history
        % Iteration
        hist.iter = [];
        
        % States
        hist.states.delta = [];
        hist.states.V = [];
        
        % Mismatches
        hist.mismatches.P = [];
        hist.mismatches.Q = [];
    end
    
    %% Determine appropriate variable and mismatch sets
    % Power mismatch -> PQ, PV buses
    sets.P = bus.id((bus.type == 0) | (bus.type == 1) | (bus.type == 2));
    % Reactive power mismatch -> PQ buses
    sets.Q = bus.id((bus.type == 0) | (bus.type == 1));
    % No temperature mismatch set
    sets.H = [];
    % Voltage angle variables -> PQ, PV buses
    sets.delta = sets.P;
    % Voltage magnitude variables -> PQ buses
    sets.V = sets.Q;
    % No temperature variables
    sets.T = [];
    
	% Determine more relevant dimensions
    M = length(sets.Q);
	
    % Ensure only 1 slack bus
    if (length(sets.P) ~= (N-1))
        error('There can be only 1 slack bus (type == 3).');
    end
    
    %% Initialize States
    % Unknowns:
	%	delta	-> N-1
	%	V       -> M
    
    % Set up Voltage Magnitude, Voltage Angle
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
    
    % Ensure column vectors
    if (size(V,2) > size(V,1))              % If # columns > # rows
        V = V';
    end
    if (size(delta,2) > size(delta,1))      % If # columns > # rows
        delta = delta';
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
        end
        
        % Evaluate YBus
        [Y,G,B,trash,trash] = makeYBus(bus,branch);
        
        % Evaluate Jacobian (Conventional Power Flow)
        J = evalJacobian(4,sets,V,delta,[],G,B,branch);

        % Evaluate mismatches
        mm = evalMismatch(sets,V,delta,[],Y,bus,[]);
        
        % Display Mismatches
		if (verbose)
			disp('Real Power Mismatch:')
				disp( mm(1:(N-1)) )
			disp('Reactive Power Mismatch:')
				disp( mm((N-1)+(1:M)) )
		end
 
        % Record this iteration in history
        if history
            % Iteration
            hist.iter = [hist.iter, iter];

            % States
            hist.states.delta = [hist.states.delta, delta];
            hist.states.V = [hist.states.V, V];

            % Mismatches
            hist.mismatches.P = [hist.mismatches.P, mm(1:(N-1))];
            hist.mismatches.Q = [hist.mismatches.Q, mm((N-1)+(1:M))];
        end
                
        % Update state vector
        x = [delta(sets.delta); V(sets.V)];

        % Calculate maximum mismatch / perform update
        err = norm(mm, inf);
        if ( err <= tol )
            if (verbose)
                disp('Conventional Power Flow Converged.')
                disp('')
            end
            break;
        else
            % Perform update
            xnew = x + J \ mm;
            delta(sets.delta) = xnew(1:(N-1));
            V(sets.V)         = xnew((N-1)+(1:M));
        end
        
		% Increment iteration
		iter = iter + 1;
		if (verbose), disp(''), end
      
	end	% End While
    
    % Warn if maximum iterations exceeded
    if iter >= maxIter
        warning(['Maximum number of power flow iterations exceeded ' ...
                 'prior to convergence within specified tolerance.']);
    end
    
    %% Parse Results of Power Flow
	% Compute maximum error at each iteration and store to history
    if history
        hist.maxErr = max( [ ...
                        max( abs( hist.mismatches.P ) ); ...
                        max( abs( hist.mismatches.Q ) ) ...
                        ] );
    end
    
    % Compute line loadings as a percent of line rating if requested
    if calcLineLoadings && history,
        branch.y = 1 ./ (branch.R + 1j*branch.X);
        V_phasor = V .* exp(1j*delta);
        
        % Compute power loss in branch
        for ik = 1:length(branch.y),
            % Get appropriate indices
            i = branch.from(ik);
            k = branch.to(ik);

            % Compute power in branch
            S(ik) = V_phasor(i) * conj(branch.y(ik)*...
                    (V_phasor(i)-V_phasor(k)));
            % Compute as a percent of rating
            LineLoadings(ik) = abs(S(ik))/branch.rating(ik);
        end
              
        hist.lineLoad = LineLoadings;
    end
    
    % Convert delta back to degrees
	delta = delta * 180 / pi;
    
    % Store final results back to 'bus'
    bus.V_mag = V';
    bus.V_angle = delta';
    
end	% End Function