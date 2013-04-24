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
%   'timing', [val]     Display timing information? TRUE/FALSE
%                       Default = FALSE
%   'history', [val]    Enable history mode (saves and returns state
%                       variables at each iteration)? TRUE/FALSE
%                       Default = FALSE
%   'lineLoad', [val]   When history mode is enabled, saves the
%                       line loadings (as a percent of line rating) into
%                       the history structure.
%                       Default = FALSE
%   'FD', [val]         Use fast-decoupled algorithm? TRUE/FALSE
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
%      troubleshooting. Enabling timing mode will calculate and display
%      timing information for various algorithm tasks.
function [V,delta,bus,branch,hist] = PF(bus,branch,varargin)
    %% Setup	
    % Default control parameters
	tol = 1e-08;		% PQH mismatch tolerance
	maxIter = 30;		% Maximum number of iterations
	verbose = false;	% Set to TRUE to spit out a bunch of diagnostics
    timing = false;     % Set to TRUE to display timing info
    history = false;    % Set to TRUE to save and return iteration history
    calcLineLoadings = false;   % Set to TRUE to return line loadings
    fastDecoupled = false;      % Set to TRUE to use fast-decoupled alg.
    
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
            case {'timing'}
				timing = varargin{2};
			case {'history'}
				history = varargin{2};
            case {'lineload'}
                calcLineLoadings = varargin{2};
            case {'fd'}
                fastDecoupled = varargin{2};
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
    
    % Timing
    if timing
        % Structure for timing
        runtime = struct();
        runtime.Total = 0;
        runtime.Setup = 0;
        runtime.YBus = 0;
        runtime.Jacobian = 0;
        runtime.MM = 0;
        runtime.Update = 0;
        
        % Start timing
        TotalTIC = tic;
        SetupTIC = TotalTIC;
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
    V = bus.V_mag;
    if exist('VInput','var')
        if length(VInput) == 1
            V(sets.V) = VInput;
        else
            V(sets.V) = VInput(sets.V);
        end  
    end
    
    % Voltage angle
    delta = bus.V_angle;
    if exist('deltaInput','var')
        if length(deltaInput) == 1
            delta(sets.delta) = deltaInput;
        else
            delta(sets.delta) = deltaInput(sets.delta);
        end  
    end
    delta = delta * pi / 180;   % Convert to radians
    
    % Ensure column vectors
    if (size(V,2) > size(V,1))              % If # columns > # rows
        V = V';
    end
    if (size(delta,2) > size(delta,1))      % If # columns > # rows
        delta = delta';
    end
    
    % Timing
    if timing, runtime.Setup = toc(SetupTIC); end
    
	%% Power Flow Algorithm - Conventional Newton Raphson
	% Perform iteration until convergence (or max. iterations)
	if ~fastDecoupled
        % Evaluate YBus
        if timing, YBusTIC = tic; end
        [Y,G,B,~,~] = makeYBus(bus,branch);
        if timing, runtime.YBus = runtime.YBus + toc(YBusTIC); end
        
        iter = 0;			% Iteration counter
        while ( iter < maxIter )
            % Display iteration info
            if (verbose)
                disp( ['Iteration: ' int2str(iter)] );
                disp( 'Voltages Magnitudes:' ); disp(V);
                disp( 'Voltages Angles:' ); disp(delta);
            end

            % Evaluate Jacobian (Conventional Power Flow)
            if timing, JacobTIC = tic; end
            J = evalJacobian(4,sets,V,delta,[],G,B,branch);
            if timing, runtime.Jacobian = runtime.Jacobian + toc(JacobTIC); end

            % Evaluate mismatches
            if timing, MismatchTIC = tic; end
            mm = evalMismatch(sets,V,delta,[],Y,bus,[]);
            if timing, runtime.MM = runtime.MM + toc(MismatchTIC); end

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
            if timing, UpdateTIC = tic; end
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
            if timing, runtime.Update = runtime.Update + toc(UpdateTIC); end

            % Increment iteration
            iter = iter + 1;
            if (verbose), disp(''), end

        end	% End While

        % Warn if maximum iterations exceeded
        if iter >= maxIter
            warning(['Maximum number of power flow iterations exceeded ' ...
                     'prior to convergence within specified tolerance.']);
        end
	end
    
	%% Power Flow Algorithm - Fast Decoupled
	% Perform iteration until convergence (or max. iterations)
	if fastDecoupled
        % Evaluate YBus
        if timing, YBusTIC = tic; end
        [Y,G,B,~,~] = makeYBus(bus,branch);
        if timing, runtime.YBus = runtime.YBus + toc(YBusTIC); end
        
        % Jacobian matrices are computed and inverted only once.
        if timing, JacobTIC = tic; end
        J = evalJacobian(5,sets,V,delta,[],G,B,branch);
        JP = J{1};              % J1
        JQ = J{2};              % J4
        if timing, runtime.Jacobian = runtime.Jacobian + toc(JacobTIC); end
        
        iter = 0;			% Iteration counter
        while ( iter < maxIter )
            % Display iteration info
            if (verbose)
                disp( ['Iteration: ' int2str(iter)] );
                disp( 'Voltages Magnitudes:' ); disp(V);
                disp( 'Voltages Angles:' ); disp(delta);
            end

            % Evaluate mismatches
            if timing, MismatchTIC = tic; end
            mm = evalMismatch(sets,V,delta,[],Y,bus,[]);
            if timing, runtime.MM = runtime.MM + toc(MismatchTIC); end

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

            % Calculate maximum mismatch / perform update
            if timing, UpdateTIC = tic; end
            err = norm(mm, inf);
            if ( err <= tol )
                if (verbose)
                    disp('Fast Decoupled Power Flow Converged.')
                    disp('')
                end
                break;
            else
                % Perform update
                delta(sets.delta) = delta(sets.delta) + JP \ mm(1:(N-1));
                V(sets.V)         = V(sets.V) + JQ \ mm((N-1)+(1:M));
            end
            if timing, runtime.Update = runtime.Update + toc(UpdateTIC); end

            % Increment iteration
            iter = iter + 1;
            if (verbose), disp(''), end

        end	% End While

        % Warn if maximum iterations exceeded
        if iter >= maxIter
            warning(['Maximum number of power flow iterations exceeded ' ...
                     'prior to convergence within specified tolerance.']);
        end
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
    
    % Compute final power injections
    % (Uses last computed values of V, delta, and Y)
    Vp = V .* ( cos(delta) + 1j*sin(delta) );	% Phasor voltages
    S = Vp .* conj(Y * Vp);                     % Complex power injections
    P = real(S);								% Real power injections
    Q = imag(S);								% Reactive power injections
    
    % Buses which have P and Q as variables
    sets.varP = setdiff(bus.id,sets.P);
    sets.varQ = setdiff(bus.id,sets.Q);
    
    % Convert delta back to degrees
	delta = delta * 180 / pi;
    
    % Store final results back to 'bus',
    % including P and Q at relevant buses
    bus.V_mag = V';
    bus.V_angle = delta';
    bus.P_net(sets.varP) = P(sets.varP);
    bus.Q_net(sets.varQ) = Q(sets.varQ);
    bus.P_gen(sets.varP) = bus.P_net(sets.varP) + bus.P_load(sets.varP);
    bus.Q_gen(sets.varQ) = bus.Q_net(sets.varQ) + bus.Q_load(sets.varQ);
    
    % Display timings
    if timing
        % Final
        runtime.Total = toc(TotalTIC);
        runtime.Overhead = runtime.Total - runtime.Setup - runtime.YBus ...
            - runtime.Jacobian - runtime.Update - runtime.MM;
        
        % Display
        disp( '---Run Time---' );
        fprintf(1, '%0d PF Iterations\n', iter);
        fprintf(1, 'Setup:                 \t%9.1f ms\t(%0.2f%%)\n', ...
            runtime.Setup * 1000, ...
            100 * runtime.Setup / runtime.Total );
        fprintf(1, 'Calculating Y Bus:     \t%9.1f ms\t(%0.2f%%)\n', ...
            runtime.YBus * 1000, ...
            100 * runtime.YBus / runtime.Total );
        fprintf(1, 'Calculating Jacobian:  \t%9.1f ms\t(%0.2f%%)\n', ...
            runtime.Jacobian * 1000, ...
            100 * runtime.Jacobian / runtime.Total );
        fprintf(1, 'Calculating Mismatches:\t%9.1f ms\t(%0.2f%%)\n', ...
            runtime.MM * 1000, ...
            100 * runtime.MM / runtime.Total );
        fprintf(1, 'Calculating Updates:   \t%9.1f ms\t(%0.2f%%)\n', ...
            runtime.Update * 1000, ...
            100 * runtime.Update / runtime.Total );
        fprintf(1, 'Overhead:              \t%9.1f ms\t(%0.2f%%)\n', ...
            runtime.Overhead * 1000, ...
            100 * runtime.Overhead / runtime.Total );
        fprintf(1, 'TOTAL:                 \t%9.1f ms\t(%0.2f%%)\n', ...
            runtime.Total * 1000, ...
            100 * runtime.Total / runtime.Total );
    end
    
end	% End Function