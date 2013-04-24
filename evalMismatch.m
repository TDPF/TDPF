%% function evalMismatch: Evaluate TDPF Mismatch Equations
% Evaluates the power (P), reactive power (Q), and temperature (H) mismatch
% equations for temperature-dependent power flow (TDPF)
% 
% SYNTAX:
%   mm = evalMismatch(sets,V,delta,T,Y,bus,branch)
%
% INPUTS:
%   sets =      Structure giving the bus/branch sets required to
%               compute the dimensions and entries of the Jacobian
%               submatrices. Must have the following elemenets:
%                   .P = Set of buses for real power mismatch equations
%                   .Q = Set of buses for reactive power mismatch equations
%                   .H = Set of branches for temperature mismatch equations
%               Setting any of these to an empty vector [] will effectively
%               remove that type of mismatch from the calculation.
%   V =         Voltage magnitude state vector (length = # of buses)
%   delta =     Voltage angle state vector (length = # of buses)
%   T =         Temperature state vector (length = # of branches)
%               (May be an empty vector if H mismatches are not required.)
%   Y =         Complex bus admittance matrix (G + jB)
%   bus =       Structure containing bus load data, in ID order, as
%               follows:
%                   .P_net      Net injected real power [pu]
%                   .Q_net      Net injected reactive power [pu]
%               The net powers need only be specified at the appropriate
%               buses from sets P and Q. However, elements for each bus
%               must be present to keep ordering correct. A structure
%               returned from makeBusStruct() will satisfy these
%               requirements.
%   branch =    Structure giving data about each temperature-dependent
%               branch element in ID order. Must have the following
%               elements:
%                   .from = Index of from bus (bus 'i' in pair [i,j])
%                   .to = Index of to bus (bus 'j' in pair [i,j])
%                   .R = Series resistance of each branch
%                   .X = Series reactance of each branch
%                   .T_amb = Ambient temperature for each branch element
%                   .R_therm = Thermal resistance of each branch element
%               A structure returned from makeBranchStruct() will satisfy
%               these requirements. For an empty H set, this argument is
%               ignored (and may safely be replaced by an empty matrix).
%
% OUTPUTS:
%   mm =        Mismatch (column) vector for the system: [DP, DQ, DH]
%
% COMMENTS:
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
function mm = evalMismatch(sets,V,delta,T,Y,bus,branch)
    %% Compute required vector/matrix dimensions
    N = length(V);          % Total number of buses
    M = length(sets.Q);     % Number of PQ buses
    R = length(sets.H);     % Number of temperature-dependent branches
    
    %% Compute some necessary quantities for branch elements
    % (Required for non-empty temperature equation set)
    if R > 0
       % Compute series conductance of branch elements
       branch.g = branch.R ./ (branch.R.^2 + branch.X.^2);
       branch.b = -branch.X ./ (branch.R.^2 + branch.X.^2);
    end
    
    %% Check and correct dimensions if needed
	% Column vectors are required for multiplication by 'Y'
    if (size(V,2) > size(V,1))
        V = V';
    end;
    if (size(delta,2) > size(delta,1))
        delta = delta';
    end;
    
	%% Compute power injections
    Vp = V .* exp(1j .* delta);         % Phasor voltages
    S = Vp .* conj(Y * Vp);             % Complex power injections
    P = real(S);                        % Real power injections
    Q = imag(S);						% Reactive power injections
	
    %% Compute Mismatches
    % (Each is evaluated only for the specified set elements.)
    
    % The mismatch column vector is as follows:
    %
    % /  DP   \
    % | (N-1) |
    % |       |
    % |  DQ   |
    % |  (M)  |
    % |       |
    % |  DH   |
    % \  (R)  /
    
    % DP
    DP = zeros(N-1,1);
    for ii = 1:(N-1)
		i = sets.P(ii);
        DP(ii) = bus.P_net(i) - P(i);
    end
    
    % DQ
    DQ = zeros(M,1);
    for ii = 1:M
		i = sets.Q(ii);
        DQ(ii) = bus.Q_net(i) - Q(i);
    end
    
    % DH
	% TO DO: Implement coefficients other than 1.0 for the temperature rise
    DH = zeros(R,1);
    for ii = 1:R
        % Indices
        ij = sets.H(ii);
		i = branch.from(ij);
		j = branch.to(ij);
        
        % Voltage magnitudes and angles
        % 'from' bus must be modified by tap ratio.
        Vi = V(i) / abs(branch.tap(ij));
        deltai = delta(i) - angle(branch.tap(ij));
        
        % 'to' bus is used directly
        Vj = V(j);
        deltaj = delta(j);
        
        % Compute mismatch
        DH(ii) = ( branch.T_amb(ij) + branch.R_therm(ij) * ...
				  branch.g(ij) * ( Vi^2 + Vj^2 - ...
				  2 * Vi * Vj * cos(deltai - deltaj) ) ) ...
				- T(ij);
    end

    %% Return mismatch vector
    mm = [DP; DQ; DH];
end




