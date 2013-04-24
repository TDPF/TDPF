%% function FD_TDPF:
% Fast Decoupled Temperature-Dependent Newton-Rhapson Power Flow
% Algorithm. Runs a Newton-Rhapson type power flow augmented with
% temperature data. This is a fast decoupled algorithm with seperate
% updates for temperature, voltage magnitudes, and voltage angles (3
% Jacobian matrices). Additionally, Jacobians are computed and inverted
% only once.
% 
% SYNTAX:
%   [V,delta,T,bus,branch,hist] = FD_TDPF(bus,branch,varargin)
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
%   'maxIter', [val]    Specify maximum number of iterations. Default = 100
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
%
% LICENSE:
%   This file is part of the Temperature Dependent Power Flow (TDPF) script
%   collection for MATLAB: see http://github.com/TDPF/TDPF.
%   
%   TDPF is free software: you may redistribute it and/or modify it under
%   the terms of the GNU General Public License (GPL) as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version.
%   
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%   
%   Text and HTML copies of the GNU General Public License should be
%   distributed with these scripts in the files 'gpl-3.0.txt' and
%   'gpl-3.0.html'. If not, plase visit <http://www.gnu.org/licenses/>
%
function [V,delta,T,bus,branch,hist] = FD_TDPF(bus,branch,varargin)
    %% Setup	
    % Default control parameters
	tol = 1e-08;		% PQH mismatch tolerance
	maxIter = 100;		% Maximum number of iterations
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
    
    %% Compute Jacobian Matrices
    % Initialize YBus
    if timing, YBusTIC = tic; end
    [Y,G,B,~,~] = makeYBus(bus,branch);
    if timing, runtime.YBus = runtime.YBus + toc(YBusTIC); end
    
    % Jacobian matrices are computed and inverted only once.
    if timing, JacobTIC = tic; end
    J = evalJacobian(3,sets,V,delta,T,G,B,branch);
    invJP = inv(J{1});      % J11 inverse
    invJQ = inv(J{2});      % J22 inverse
    invJH = inv(J{3});      % J33 inverse
%    JP = J{1};              % J11
%    JQ = J{2};              % J22
%    JH = J{3};              % J33
    if timing, runtime.Jacobian = runtime.Jacobian + toc(JacobTIC); end

    % Note: There are two possible ways to handle the matrix inversions in
    % the fast decoupled method:
    %   (1) Invert the matrices once and then use matrix-vector
    %       multiplication thereafter in each iteration. (Uncommented)
    %   (2) Do not invert the matrices at all and use Gaussian elimination
    %       in each iteration. (Commented)
    % Method (1) is O(n^3) for the matrix inversion and then O(n^2) at each 
    % iteration to perform the multiplication. Method (2) is O(n^3) at each
    % iteration to perform the Gaussian elimination. Thus, method (1) is
    % (theoretically) faster, although method (2) is more stable
    % numerically with respect to badly conditioned matrices. (See MATLAB
    % documentation for the inv() function.)
    %
    % For the IEEE 30-bus test system, a quick test showed a 7% speedup for
    % method (1) compared to method (2).
    
    %% Fast-Decoupled Power Flow - 1 Iteration
    % Perform a single iteration of fast-decoupled power flow in order to
    % improve the initial estimate of line temperatures.

    % Record initial states in history
    if history
        % Iteration 0 = starting conditions
        hist.iter = 0;
        
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
    
    if timing, MismatchTIC = tic; end
    mm = evalMismatch(s2,V,delta,[],Y,bus,branch);
    if timing, runtime.MM = runtime.MM + toc(MismatchTIC); end
    
    % Perform update
    if timing, UpdateTIC = tic; end
    delta(sets.delta) = delta(sets.delta) + invJP * mm(1:(N-1));
    V(sets.V)         = V(sets.V)         + invJQ * mm((N-1)+(1:M));
 %   delta(sets.delta) = delta(sets.delta) + JP \ mm(1:(N-1));
 %   V(sets.V)         = V(sets.V)         + JQ \ mm((N-1)+(1:M));
	if timing, runtime.Update = runtime.Update + toc(UpdateTIC); end
 
    %% Initial Temperature Estimate
    % Using the results of the power flow iteration above, compute better
    % guesses for initial starting temperatures.
    if timing, TempTIC = tic; end
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
    if timing, runtime.Temp = runtime.Temp + toc(TempTIC); end
        
    %% Power Flow Algorithm
	% Perform iteration until convergence (or max. iterations)
	iter = 1;			% Iteration counter
                        % (Starts at 1 to account for initial iteration
                        %  of conventional power flow)
	while ( iter < maxIter )
		% Display iteration info
        if (verbose)
			disp( ['Iteration: ' int2str(iter)] );
			disp( 'Voltages Magnitudes:' ); disp(V);
			disp( 'Voltages Angles:' ); disp(delta);
            disp( 'Branch Temperatures:' ); disp(T);
        end
        
        % Evaluate YBus
        if timing, YBusTIC = tic; end
        [Y,~,~,~,~] = makeYBus(bus,branch);
        if timing, runtime.YBus = runtime.YBus + toc(YBusTIC); end

        % Evaluate mismatches
        if timing, MismatchTIC = tic; end
        mm = evalMismatch(sets,V,delta,T,Y,bus,branch);
        if timing, runtime.MM = runtime.MM + toc(MismatchTIC); end
        
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
                
		% Calculate maximum error / perform update
        if timing, UpdateTIC = tic; end
		err = norm(mm, inf);
        if ( err <= tol )
            if (verbose)
            	disp('FD-TDPF Power Flow Converged.')
            end
            break;
        else
            % Perform update
            delta(sets.delta) = delta(sets.delta) + invJP * mm(1:(N-1));
            V(sets.V)         = V(sets.V) + invJQ * mm((N-1)+(1:M));
            T(sets.T)         = T(sets.T) + invJH * mm((N-1+M)+(1:R));
%            delta(sets.delta) = delta(sets.delta) + JP \ mm(1:(N-1));
%            V(sets.V)         = V(sets.V) + JQ \ mm((N-1)+(1:M));
%            T(sets.T)         = T(sets.T) + JH \ mm((N-1+M)+(1:R));
        end
        if timing, runtime.Update = runtime.Update + toc(UpdateTIC); end

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
        fprintf(1, '%0d FD-TDPF Iterations\n', iter);
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