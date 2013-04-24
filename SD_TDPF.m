%% function SD_TDPF:
% Sequnetially Decoupled Temperature-Dependent Newton-Rhapson Power Flow
% Algorithm. Runs a Newton-Rhapson type power flow algorithm and a set of
% temperature corrections in sequence until convergence, modifying the
% admittance matrix after each temperature correction.

% SYNTAX:
%   [V,delta,T,bus,branch,hist] = SD_TDPF(bus,branch,varargin)
%
% INPUTS:
%	bus =       Bus structure (as returned from importCaseData())
%	branch =    Branch structure (as returned from importCaseData())
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
%   'maxIter', [val]    Specify maximum number of outer iterations.
%                       Default = 10
%   'maxSubIter', [val] Specify maximum number of sub-iterations, or inner
%                       iterations, per outer iteration.
%                       Default = 10
%   'verbose', [val]    Enable verbose mode (lots of reporting)? TRUE/FALSE
%                       Default = FALSE
%   'history', [val]    Enable history mode (saves and returns state
%                       variables at each iteration)? TRUE/FALSE
%                       Default = FALSE
%   'timing', [val]     Display timing information? TRUE/FALSE
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
%      troubleshooting. Enabling timing mode will calculate and display
%      timing information for various algorithm tasks.
function [V,delta,T,bus,branch,hist] = SD_TDPF(bus,branch,varargin)
    %% Setup	
    % Default control parameters
	tol = 1e-08;		% PQH mismatch tolerance
	maxIter = 10;		% Maximum number of outer iterations
    maxSubIter = 10;    % Maximum number of sub-iterations
	verbose = false;	% Set to TRUE to spit out a bunch of diagnostics
    history = false;    % Set to TRUE to save and return iteration history
    timing = false;     % Set to TRUE to display timing info
	
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
			case {'maxsubiter'}
				maxSubIter = varargin{2};
			case {'verbose'}
				verbose = varargin{2};
			case {'history'}
				history = varargin{2};
            case {'timing'}
				timing = varargin{2};
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
    if history
        % Iteration/sub-iteration
        hist.iter = [];
        hist.subIter = [];
        
        % States
        hist.states.delta = [];
        hist.states.V = [];
        hist.states.T = [];
        
        % Mismatches
        hist.mismatches.P = [];
        hist.mismatches.Q = [];
        hist.mismatches.H = [];
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
        runtime.Temp = 0;
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
    branch.b = -branch.X ./ (branch.R.^2 + branch.X.^2);
    
    % Timing
    if timing, runtime.Setup = toc(SetupTIC); end
        
    %% Power Flow Algorithm
    % Sets without temperature
    setsPQ = sets;
    setsPQ.H = [];
    
	% Perform iteration until convergence (or max. iterations)
	iter = 0;			% Outer iteration counter
	while ( iter < maxIter )
		% Display iteration info
        if (verbose)
			disp( ['Outer Iteration: ' int2str(iter)] );
        end
        
        % Evaluate YBus
        if timing, YBusTIC = tic; end
        [Y,G,B,~,~] = makeYBus(bus,branch);
        if timing, runtime.YBus = runtime.YBus + toc(YBusTIC); end
 
        subIter = 0;        % Inner iteration counter
        while ( subIter < maxSubIter )
            % Display iteration info
            if (verbose)
                disp( ['Sub-iteration: ' int2str(subIter)] );
                disp( 'Voltages Magnitudes:' ); disp(V);
                disp( 'Voltages Angles:' ); disp(delta);
            end
            
            % Evaluate Jacobian (Conventional Power Flow)
            if timing, JacobTIC = tic; end
            J = evalJacobian(4,setsPQ,V,delta,[],G,B,branch);
            if timing, runtime.Jacobian = runtime.Jacobian + toc(JacobTIC); end

            % Evaluate mismatches
            if timing, MismatchTIC = tic; end
            mm = evalMismatch(setsPQ,V,delta,[],Y,bus,branch);
            if timing, runtime.MM = runtime.MM + toc(MismatchTIC); end

            % Display Power Mismatches
            if (verbose)
                disp('Real Power Mismatch:')
                    disp( mm(1:(N-1)) )
                disp('Reactive Power Mismatch:')
                    disp( mm((N-1)+(1:M)) )
            end
 
            % Record this sub-iteration in history
            if history
                % Iteration
                hist.iter = [hist.iter, iter];
                hist.subIter = [hist.subIter, subIter];

                % States
                hist.states.delta = [hist.states.delta, delta];
                hist.states.V = [hist.states.V, V];
                hist.states.T = [hist.states.T, T];

                % Mismatches
                mm2 = evalMismatch(sets,V,delta,T,Y,bus,branch);
                hist.mismatches.P = [hist.mismatches.P, mm2(1:(N-1))];
                hist.mismatches.Q = [hist.mismatches.Q, mm2((N-1)+(1:M))];
                hist.mismatches.H = [hist.mismatches.H, mm2((N-1+M)+(1:R))];
            end
                
            % Update state vector
            x = [delta(sets.delta); V(sets.V)];
        
            % Calculate maximum P/Q error / perform update
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
            
            % Increment sub-iteration
            subIter = subIter + 1;
            if (verbose), disp(''), end
        end
        
        % Warn if maximum iterations exceeded
        if subIter >= maxSubIter
            warning(['Maximum number of sub-iterations exceeded ' ...
                     'prior to convergence within specified tolerance.']);
        end
        
        % Perform an update of system temperatures and resistances
        if timing, TempTIC = tic; end
        for ij = sets.T
            % Get appropriate indices
            i = branch.from(ij);
            j = branch.to(ij);

            % Voltage magnitudes and angles
            % 'from' bus must be modified by tap ratio.
            Vi = V(i) / abs(branch.tap(ij));
            deltai = delta(i) - angle(branch.tap(ij));

            % 'to' bus is used directly
            Vj = V(j);
            deltaj = delta(j);
            
            % Compute power loss in branch
            PLoss = branch.g(ij) * ( Vi^2 + Vj^2 - ...
                      2 * Vi * Vj * cos(deltai - deltaj) );

            % Calculate corresponding temperature rise and update the
            % temperature estimate
            T(ij) = branch.T_amb(ij) + PLoss * branch.R_therm(ij);
        end
        if timing, runtime.Temp = runtime.Temp + toc(TempTIC); end
        
        % Display new temperatures
        if (verbose)
            disp( 'New Branch Temperatures:' ); disp(T);
        end
        
        % Update resistances using new temperatures
        % (Only needs to occur for temperature-dependant branches)
        if timing, TempTIC = tic; end
        branch.T(sets.T) = T(sets.T);
        branch.R(sets.T) = branch.R_ref(sets.T) .* ...
            ( (branch.T(sets.T) + branch.T_f(sets.T)) ./ ...
            (branch.T_ref(sets.T) + branch.T_f(sets.T)) );
        branch.g(sets.T) = branch.R(sets.T) ./ ...
            (branch.R(sets.T).^2 + branch.X(sets.T).^2);
        branch.b(sets.T) = -branch.X(sets.T) ./ ...
            (branch.R(sets.T).^2 + branch.X(sets.T).^2);
        if timing, runtime.Temp = runtime.Temp + toc(TempTIC); end
        
        % Evaluate mismatches
        if timing, MismatchTIC = tic; end
        mm = evalMismatch(sets,V,delta,T,Y,bus,branch);
        if timing, runtime.MM = runtime.MM + toc(MismatchTIC); end
        
        % Display Temperature Mismatches
        if (verbose)
            disp('Temperature Mismatch:')
				disp( mm((N-1+M)+(1:R)) )
        end
            
        % Calculate maximum error, including temperatures, to ensure
        % convergance of entire algorithm
		err = norm(mm, inf);
        if ( err <= tol )
            if (verbose)
            	disp('SD-TDPF Power Flow Converged.')
            end
            break;
        end
        
		% Increment iteration
		iter = iter + 1;
		if (verbose), disp(''), end
       
	end	% End While
    
    % Warn if maximum iterations exceeded
    if iter >= maxIter
        warning(['Maximum number of outer iterations exceeded ' ...
                 'prior to convergence within specified tolerance.']);
    end
    
    %% Parse Results of Power Flow
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
    
    % Store final results back to 'bus' and 'branch',
    % including P and Q at relevant buses
    bus.V_mag = V';
    bus.V_angle = delta';
    bus.P_net(sets.varP) = P(sets.varP);
    bus.Q_net(sets.varQ) = Q(sets.varQ);
    bus.P_gen(sets.varP) = bus.P_net(sets.varP) + bus.P_load(sets.varP);
    bus.Q_gen(sets.varQ) = bus.Q_net(sets.varQ) + bus.Q_load(sets.varQ);
    branch.T = T';
    
    % Compute maximum error at each iteration and store to history
    if history
        hist.maxErr = max( [ ...
                        max( abs( hist.mismatches.P ) ); ...
                        max( abs( hist.mismatches.Q ) ); ...
                        max( abs( hist.mismatches.H ) ) ...
                        ] );
    end
    
    % Display timings
    if timing
        % Final
        runtime.Total = toc(TotalTIC);
        runtime.Overhead = runtime.Total - runtime.Setup - runtime.YBus ...
            - runtime.Jacobian - runtime.Update - runtime.MM ...
            - runtime.Temp;
        
        % Display
        disp( '---Run Time---' );
        fprintf(1, '%0d (Outer) SD-TDPF Iterations\n', iter);
        fprintf(1, 'Setup:                  \t%9.1f ms\t(%0.2f%%)\n', ...
            runtime.Setup * 1000, ...
            100 * runtime.Setup / runtime.Total );
        fprintf(1, 'Calculating Y Bus:      \t%9.1f ms\t(%0.2f%%)\n', ...
            runtime.YBus * 1000, ...
            100 * runtime.YBus / runtime.Total );
        fprintf(1, 'Calculating Jacobian:   \t%9.1f ms\t(%0.2f%%)\n', ...
            runtime.Jacobian * 1000, ...
            100 * runtime.Jacobian / runtime.Total );
        fprintf(1, 'Calculating Mismatches: \t%9.1f ms\t(%0.2f%%)\n', ...
            runtime.MM * 1000, ...
            100 * runtime.MM / runtime.Total );
        fprintf(1, 'Calculating Updates:    \t%9.1f ms\t(%0.2f%%)\n', ...
            runtime.Update * 1000, ...
            100 * runtime.Update / runtime.Total );
        fprintf(1, 'Updating Temperatures:  \t%9.1f ms\t(%0.2f%%)\n', ...
            runtime.Temp * 1000, ...
            100 * runtime.Temp / runtime.Total );
        fprintf(1, 'Overhead:               \t%9.1f ms\t(%0.2f%%)\n', ...
            runtime.Overhead * 1000, ...
            100 * runtime.Overhead / runtime.Total );
        fprintf(1, 'TOTAL:                  \t%9.1f ms\t(%0.2f%%)\n', ...
            runtime.Total * 1000, ...
            100 * runtime.Total / runtime.Total );
    end
    
end	% End Function