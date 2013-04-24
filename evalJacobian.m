%% function evalJacobian: Evaluate Jacobian Matrix
% Evaluates the Jacobian matrix for temperature-dependent power flow. The
% structure of the Jacobian(s) returned depends on the type of TDPF being
% performed, as specified by the user when the function is called.
% 
% SYNTAX:
%   J = evalJacobian(type,sets,V,delta,T,G,B,branches)
%
% INPUTS:
%   type =      Indicates type of Jacobian(s) to return, based on type of
%               TDPF being performed. Options are:
%               1 -- Jacobian matrix for fully-coupled TDPF (FC-TDPF)
%               2 -- Jacobian submatrices J11-J22 and J33 for
%                    partially-decoupled TDPF (PD-TDPF)
%               3 -- Jacobian submatrices J11, J22, and J33 for
%                    fast-decoupled TDPF (FD-TDPF)
%               4 -- Jacobian matrix for conventional power flow (J11-J22);
%                    needed for conventional power flow or for
%                    sequentially decoupled TDPF (SD-TDPF)
%               5 -- Jacobian submatrices J11 and J22 for conventional
%                    fast-decoupled power flow
%   sets =      Structure giving the bus/branch sets required to
%               compute the dimensions and entries of the Jacobian
%               submatrices. Must have the following elements:
%                   .P = Set of buses for real power mismatch equations
%                   .Q = Set of buses for reactive power mismatch equations
%                   .H = Set of branches for temperature mismatch equations
%                   .V = Set of buses with unknown voltage magnitudes
%                   .delta = Set of buses with unknown voltage angles
%                   .T = Set of branches with unknown temperatures
%   V =         Voltage magnitude state vector (length = # of buses)
%   delta =     Voltage angle state vector (length = # of buses)
%   T =         Temperature state vector (length = # of branches)
%   G =         Bus conductance matrix (real part of YBus)
%   B =         Bus susceptance matrix (imaginary part of YBus)
%   branch =    Structure giving data about each temperature-dependent
%               branch element in ID order. Must have the following
%               elements:
%                   .from = Index of from bus (bus 'i' in pair [i,j])
%                   .to = Index of to bus (bus 'j' in pair [i,j])
%                   .R = Series resistance of each branch
%                   .X = Series reactance of each branch
%                   .R_ref = Reference series resistance of each branch
%                   .T_ref = Reference temperature for each branch
%                   .T_f = Thermal coefficentfor each branch 
%                   .R_therm = Thermal resistance of each branch element
%               A structure returned from makeBranchStruct() will satisfy
%               these requirements.
%
% OUTPUTS:
%   J =         Jacobian matrix (or matrices) with structure dependent on
%               the input variable 'type' as follows:
%               1 -- J is a single matrix representing the entire Jacobian
%                    J11-J33 for use with FC-TDPF
%               2 -- J is a cell structure containing two Jacobian
%                    submatrices: J11-J22 and J33 (in that order).
%                    For use with PD-TDPF.
%               3 -- J is a cell structure containing three Jacobian
%                    submatrices: J11, J22 and J33 (in that order).
%                    For use with FD-TDPF.
%               4 -- J is single matrix representing the Jacobian J11-J22
%                    for conventional power flow.
%                    For use with conventional power flow or SD-TDPF.
%               5 -- J is a cell structure containing two Jacobian
%                    submatrices: J11 and J22 (in that order).
%               	 For use with fast-decoupled power flow.
%
% COMMENTS:
%   1. The Jacobian matrices are numbered here as in the paper:
%      / J11 = dP/dd     J12 = dP/dV     J13 = dP/dT \
%      | (N-1)x(N-1)     (N-1)x(M)       (N-1)x(R)   |
%      |                                             |
%      | J21 = dQ/dd     J22 = dQ/dV     J23 = dQ/dT |
%      | (M)x(N-1)       (M)x(M)         (M)x(R)     |
%      |                                             |
%      | J31 = dH/dd     J32 = dH/dV     J33 = dH/dT |
%      \ (R)x(N-1)       (R)x(M)         (R)x(R)     /
%
%      N = Number of buses
%      M = Number of PQ buses
%      R = Number temperature-dependent branches
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
function J = evalJacobian(type,sets,V,delta,T,G,B,branch)
    %% Compute required vector/matrix dimensions
    N = length(V);          % Total number of buses
    M = length(sets.V);     % Number of PQ buses
    L = length(branch.id);  % Total number of branches
    R = length(sets.T);     % Number of temperature-dependent branches

    %% Compute some necessary quantities for branch elements
    % (Required for all cases except type==4 and type==5)
    if (type ~= 4) && (type ~= 5)
       % Compute series admittance of branch elements
       branch.g = branch.R ./ (branch.R.^2 + branch.X.^2);
       branch.b = -branch.X ./ (branch.R.^2 + branch.X.^2);
    end
    
    %% Check and correct dimensions if needed
	% Row vectors are required for element-wise multiplication with rows of
    % 'G' and 'B' matrices.
    if (size(V,1) > size(V,2))
        V = V';
    end;
    if (size(delta,1) > size(delta,2))
        delta = delta';
    end;
    if (size(T,1) > size(T,2))
        T = T';
    end;
    
    %% Compute Sparsity Information
    % It's helpful to precompute some sparsity information about YBus to
    % make it easier to compute the Jacobian matrices efficiently.
    
    % Generates a sparse matrix indicating non-zero elements of Y
    % (Note that diagonals are considered non-zero even if they happen to
    % be zero for some odd reason)
    nzY = (G ~= 0) | (B ~= 0) | speye(N);
    
    % Generates a sparse connection matrix mapping from...
    %   Branch elements by row, that are connected to...
    %   Bus elements by column
    % Each row has exactly 2 corresponding nonzero columns.
    nzL = sparse( [1:L,1:L], [branch.from,branch.to], ones(1,2*L), L, N);
    
    
    %% Compute Jacobian submatrices
    % (Each is evaluated only if needed)
    
    % J11 = dP/dd    (N-1) x (N-1)
    % Needed for all types
        % Find indices of all nonzero YBus elements + diagonals where:
        %   The row    exists in sets.P
        %   The column exists in sets.delta
        [jRow, jCol]  = find( nzY( sets.P, sets.delta ) );
        
        % NOTE: The indexing is now relative to J11 rather than YBus. This
        % is important when translating back to V and delta values!
        
        % Storage for J11 elements
        jElem = zeros( size(jRow) );
        
        % Generate each nonzero element of J11
        for count = 1:length(jRow)
            % Find indices with respect to J11
            ii = jRow(count);
            kk = jCol(count);
            
            % Find indices with respect to YBus
            i = sets.P(ii);
            k = sets.delta(kk);            

            % Same bus
            if (k == i)
                jElem(count) = ...
                    V(i) * sum( V .* ( ...
                    -G(i,:) .* sin(delta(i) - delta) ...
                    +B(i,:) .* cos(delta(i) - delta) ...
                    )) - V(i)^2 * B(i,i);
            % Different buses
            else
                jElem(count) = ...
                    V(i)*V(k)*( ...
                    G(i,k)*sin(delta(i)-delta(k)) - ...
                    B(i,k)*cos(delta(i)-delta(k)) );
            end
        end
        
        % Generate J11
        J11 = sparse(jRow, jCol, jElem, N-1, N-1);
    
    
    % J12 = dP/dV    (N-1) x (M)
    if any( type == [1 2 4] )
        % Find indices of all nonzero YBus elements + diagonals where:
        %   The row    exists in sets.P
        %   The column exists in sets.V
        [jRow, jCol]  = find( nzY( sets.P, sets.V ) );
        
        % NOTE: The indexing is now relative to J12 rather than YBus. This
        % is important when translating back to V and delta values!
        
        % Storage for J12 elements
        jElem = zeros( size(jRow) );
        
        % Generate each nonzero element of J12
        for count = 1:length(jRow)
            % Find indices with respect to J12
            ii = jRow(count);
            kk = jCol(count);
            
            % Find indices with respect to YBus
            i = sets.P(ii);
            k = sets.V(kk);            

            % Same bus
            if (k == i)
                jElem(count) = V(i) * G(i,i) + ...
                    sum( V .* ( ...
                    G(i,:) .* cos(delta(i) - delta) + ...
                    B(i,:) .* sin(delta(i) - delta) ...
                    ));
            % Different buses
            else
                jElem(count) = V(i) * ( ...
                    G(i,k) * cos(delta(i) - delta(k)) + ...
                    B(i,k) * sin(delta(i) - delta(k)) );
            end
        end
        
        % Generate J12
        J12 = sparse(jRow, jCol, jElem, N-1, M);
    end
    
    
    % J21 = dQ/dd    (M) x (N-1)
    if any( type == [1 2 4] )
        % Find indices of all nonzero YBus elements + diagonals where:
        %   The row    exists in sets.Q
        %   The column exists in sets.delta
        [jRow, jCol]  = find( nzY( sets.Q, sets.delta ) );
        
        % NOTE: The indexing is now relative to J21 rather than YBus. This
        % is important when translating back to V and delta values!
        
        % Storage for J21 elements
        jElem = zeros( size(jRow) );
        
        % Generate each nonzero element of J21
        for count = 1:length(jRow)
            % Find indices with respect to J21
            ii = jRow(count);
            kk = jCol(count);
            
            % Find indices with respect to YBus
            i = sets.Q(ii);
            k = sets.delta(kk);            

            % Same bus
            if (k == i)
                jElem(count) = ...
                    V(i) * sum( V .* ( ...
                    G(i,:) .* cos(delta(i) - delta) + ...
                    B(i,:) .* sin(delta(i) - delta) ...
                    )) - V(i)^2 * G(i,i);
            % Different buses
            else
                jElem(count) = ...
                    V(i)*V(k)*( ...
                    -G(i,k)*cos(delta(i)-delta(k)) - ...
                    B(i,k)*sin(delta(i)-delta(k)) );
            end
        end
        
        % Generate J21
        J21 = sparse(jRow, jCol, jElem, M, N-1);
    end
    
    
    % J22 = dQ/dV	(M) x (M)
    % Needed for all types
        % Find indices of all nonzero YBus elements + diagonals where:
        %   The row    exists in sets.Q
        %   The column exists in sets.V
        [jRow, jCol]  = find( nzY( sets.Q, sets.V ) );
        
        % NOTE: The indexing is now relative to J22 rather than YBus. This
        % is important when translating back to V and delta values!
        
        % Storage for J22 elements
        jElem = zeros( size(jRow) );
        
        % Generate each nonzero element of J22
        for count = 1:length(jRow)
            % Find indices with respect to J22
            ii = jRow(count);
            kk = jCol(count);
            
            % Find indices with respect to YBus
            i = sets.Q(ii);
            k = sets.V(kk);            

            % Same bus
            if (k == i)
                jElem(count) = -V(i) * B(i,i) + ...
                    sum( V .* ( ...
                    G(i,:) .* sin(delta(i) - delta) - ...
                    B(i,:) .* cos(delta(i) - delta) ...
                    ));
            % Different buses
            else
                jElem(count) = V(i) * ( ...
                    G(i,k) * sin(delta(i) - delta(k)) - ...
                    B(i,k) * cos(delta(i) - delta(k)) );
            end
        end
        
        % Generate J22
        J22 = sparse(jRow, jCol, jElem, M, M);
    
    
    % J31 = dH/dd    (R) x (N-1)
    if type == 1
        % Find indices of all nonzero branch-to-bus connections where
        %   The row    exists in sets.H
        %   The column exists in sets.delta
        [jRow, jCol]  = find( nzL( sets.H, sets.delta ) );
        
        % NOTE: The indexing is now relative to J31.
        
        % Storage for J31 elements
        jElem = zeros( size(jRow) );
        
        % Generate each nonzero element of J31
        for count = 1:length(jRow)
            % Find indices with respect to J31
            ii = jRow(count);       % Branch-related
            kk = jCol(count);       % Bus-related
            
            % Find indices with respect to YBus and branch
            ij = sets.H(ii);        % Branch id
            i = branch.from(ij);    % From bus
            j = branch.to(ij);      % To bus
            k = sets.delta(kk);     % Bus corresponding to 'delta'
            
            % Voltage magnitudes and angles
            % 'from' bus must be modified by tap ratio.
            Vi = V(i) / abs(branch.tap(ij));
            deltai = delta(i) - angle(branch.tap(ij));

            % 'to' bus is used directly
            Vj = V(j);
            deltaj = delta(j);

            % Calculate result
            if (k == i)
                jElem(count) = -2 * branch.R_therm(ij) * ...
                               branch.g(ij) * ...
                               Vi * Vj * sin(deltai - deltaj);
            elseif (k == j)
                jElem(count) = 2 * branch.R_therm(ij) * ...
                               branch.g(ij) * ...
                               Vi * Vj * sin(deltai - deltaj);
            else % This shouldn't ever happen!
                warning(['Indexing error in calculation of dH/ddelta! ' ...
                         'Resulting Jacobian matrix may have errors.']);
            end
        end
        
        % Generate J31
        J31 = sparse(jRow, jCol, jElem, R, N-1);
    end
    
    
    % J32 = dH/dV    (R) x (M)
    if type == 1
        % Find indices of all nonzero branch-to-bus connections where
        %   The row    exists in sets.H
        %   The column exists in sets.V
        [jRow, jCol]  = find( nzL( sets.H, sets.V ) );
        
        % NOTE: The indexing is now relative to J32.
        
        % Storage for J32 elements
        jElem = zeros( size(jRow) );
        
        % Generate each nonzero element of J32
        for count = 1:length(jRow)
            % Find indices with respect to J32
            ii = jRow(count);       % Branch-related
            kk = jCol(count);       % Bus-related
            
            % Find indices with respect to YBus and branch
            ij = sets.H(ii);        % Branch id
            i = branch.from(ij);    % From bus
            j = branch.to(ij);      % To bus
            k = sets.V(kk);         % Bus corresponding to 'V'
            
            % Voltage magnitudes and angles
            % 'from' bus must be modified by tap ratio.
            Vi = V(i) / abs(branch.tap(ij));
            deltai = delta(i) - angle(branch.tap(ij));

            % 'to' bus is used directly
            Vj = V(j);
            deltaj = delta(j);

            % Calculate result
            if (k == i)
                jElem(count) = -2 * branch.R_therm(ij) * ...
                                branch.g(ij) * ...
                                ( Vi - Vj * cos(deltai - deltaj) );
            elseif (k == j)
                jElem(count) = -2 * branch.R_therm(ij) * ...
                                branch.g(ij) * ...
                                ( Vj - Vi * cos(deltai - deltaj) );
            else % This shouldn't ever happen!
                warning(['Indexing error in calculation of dH/dV! ' ...
                         'Resulting Jacobian matrix may have errors.']);
            end
        end
        
        % Generate J32
        J32 = sparse(jRow, jCol, jElem, R, M);
    end
    
    % J13 = dP/dT    (N-1) x (R)
    if type == 1
        % Find indices of all nonzero bus-to-branch connections where
        %   The row    exists in sets.P
        %   The column exists in sets.T
        % Note that this requires using the transpose of nzL
        [jRow, jCol]  = find( nzL( sets.T, sets.P )' );
        
        % NOTE: The indexing is now relative to J13.
        
        % Storage for J13 elements
        jElem = zeros( size(jRow) );
        
        % Generate each nonzero element of J13
        for count = 1:length(jRow)
            % Find indices with respect to J13
            ii = jRow(count);       % Bus-related
            kk = jCol(count);       % Branch-related
            
            % Find indices with respect to YBus and branch
            i = sets.P(ii);         % Bus id for 'P'
            kn = sets.T(kk);        % Branch id
            k = branch.from(kn);    % 'from' bus for branch
            n = branch.to(kn);      % 'to' bus for branch
            
            % Compute dg/dT, db/dT for this branch
            dgdT = ( branch.X(kn)^2 - branch.R(kn)^2 ) / ...
                   ( branch.R(kn)^2 + branch.X(kn)^2 )^2 * ...
                     branch.R_ref(kn) / ...
                   ( branch.T_ref(kn) + branch.T_f(kn) );
            dbdT = ( 2 * branch.X(kn) * branch.R(kn) ) / ...
                   ( branch.R(kn)^2 + branch.X(kn)^2 )^2 * ...
                     branch.R_ref(kn) / ...
                   ( branch.T_ref(kn) + branch.T_f(kn) );
            
            % Calculate result
            if (k == i)
                % Voltage magnitudes and angles
                % 'from' bus must be modified by tap ratio.
                Vi = V(i) / abs(branch.tap(kn));
                deltai = delta(i) - angle(branch.tap(kn));
                
                % 'to' bus is used directly
                Vn = V(n);
                deltan = delta(n);

                % Compute dP/dT
                jElem(count) = ( Vi^2 - Vi * Vn * ...
                                 cos(deltai - deltan) ) * dgdT - ...
                               Vi * Vn * ...
                                 sin(deltai - deltan) * dbdT;
            elseif (n == i)
                % Voltage magnitudes and angles
                % 'from' bus must be modified by tap ratio.
                Vk = V(k) / abs(branch.tap(kn));
                deltak = delta(k) - angle(branch.tap(kn));

                % 'to' bus is used directly
                Vi = V(i);
                deltai = delta(i);

                % Compute dP/dT
                jElem(count) = ( Vi^2 - Vi * Vk * ...
                                 cos(deltai - deltak) ) * dgdT - ...
                               Vi * Vk * ...
                                 sin(deltai - deltak) * dbdT;
            else % This shouldn't ever happen!
                warning(['Indexing error in calculation of dP/dT! ' ...
                         'Resulting Jacobian matrix may have errors.']);
            end
        end
        
        % Generate J13
        J13 = sparse(jRow, jCol, jElem, N-1, R);
    end
    
    % J23 = dQ/dT    (M) x (R)
    if type == 1
        % Find indices of all nonzero bus-to-branch connections where
        %   The row    exists in sets.Q
        %   The column exists in sets.T
        % Note that this requires using the transpose of nzL
        [jRow, jCol]  = find( nzL( sets.T, sets.Q )' );
        
        % NOTE: The indexing is now relative to J23.
        
        % Storage for J23 elements
        jElem = zeros( size(jRow) );
        
        % Generate each nonzero element of J23
        for count = 1:length(jRow)
            % Find indices with respect to J23
            ii = jRow(count);       % Bus-related
            kk = jCol(count);       % Branch-related
            
            % Find indices with respect to YBus and branch
            i = sets.Q(ii);         % Bus id for 'Q'
            kn = sets.T(kk);        % Branch id
            k = branch.from(kn);    % 'from' bus for branch
            n = branch.to(kn);      % 'to' bus for branch
            
            % Compute dg/dT, db/dT for this branch
            dgdT = ( branch.X(kn)^2 - branch.R(kn)^2 ) / ...
                   ( branch.R(kn)^2 + branch.X(kn)^2 )^2 * ...
                     branch.R_ref(kn) / ...
                   ( branch.T_ref(kn) + branch.T_f(kn) );
            dbdT = ( 2 * branch.X(kn) * branch.R(kn) ) / ...
                   ( branch.R(kn)^2 + branch.X(kn)^2 )^2 * ...
                     branch.R_ref(kn) / ...
                   ( branch.T_ref(kn) + branch.T_f(kn) );
            
            % Calculate result
            if (k == i)
                % Voltage magnitudes and angles
                % 'from' bus must be modified by tap ratio.
                Vi = V(i) / abs(branch.tap(kn));
                deltai = delta(i) - angle(branch.tap(kn));

                % 'to' bus is used directly
                Vn = V(n);
                deltan = delta(n);
                
                % Compute dQ/dT
                jElem(count) = -Vi * Vn * ...
                                 sin(deltai - deltan) * dgdT + ...
                               ( -Vi^2 + Vi * Vn * ...
                                 cos(deltai - deltan) ) * dbdT;
            elseif (n == i)
                % Voltage magnitudes and angles
                % 'from' bus must be modified by tap ratio.
                Vk = V(k) / abs(branch.tap(kn));
                deltak = delta(k) - angle(branch.tap(kn));

                % 'to' bus is used directly
                Vi = V(i);
                deltai = delta(i);

                % Compute dQ/dT
                jElem(count) = -Vi * Vk * ...
                                 sin(deltai - deltak) * dgdT + ...
                               ( -Vi^2 + Vi * Vk * ...
                                 cos(deltai - deltak) ) * dbdT;
            else % This shouldn't ever happen!
                warning(['Indexing error in calculation of dQ/dT! ' ...
                         'Resulting Jacobian matrix may have errors.']);
            end
        end
        
        % Generate J23
        J23 = sparse(jRow, jCol, jElem, M, R);
    end
    
    
    % J33 = dH/dT    (R) x (R)
    if any( type == [1 2] )
        % Row, column, and element vectors for J33
        % NOTE: J33 has only diagonal entries
        jRow = 1:R; jCol = 1:R;
        jElem = zeros( size(jRow) );
        
        % Generate elements of J33
        for kk = 1:R
            % Indices
            kn = sets.T(kk);
            k = branch.from(kn);
            n = branch.to(kn);
            
            % Voltage magnitudes and angles
            % 'from' bus must be modified by tap ratio.
            Vk = V(k) / abs(branch.tap(kn));
            deltak = delta(k) - angle(branch.tap(kn));

            % 'to' bus is used directly
            Vn = V(n);
            deltan = delta(n);
            
            % Compute dg/dT, db/dT
            dgdT = ( branch.X(kn)^2 - branch.R(kn)^2 ) / ...
                   ( branch.R(kn)^2 + branch.X(kn)^2 )^2 * ...
                     branch.R_ref(kn) / ...
                   ( branch.T_ref(kn) + branch.T_f(kn) );
            
            % Compute dH/dT
        	jElem(kk) = 1 - branch.R_therm(kn) * ...
                        ( Vk^2 + Vn^2 - ...
                          2 * Vk * Vn * cos(deltak - deltan) ) * dgdT;
        end
        
        % Generate J33
        J33 = sparse(jRow, jCol, jElem, R, R);
    elseif type == 3
        % For FD-TDPF, J33 is approximated by the identity matrix
        J33 = speye(R);
    end
        
    %% Return Type-Dependent Output
    switch type
        case 1  % FC-TDPF
            J = [J11 J12 J13; J21 J22 J23; J31 J32 J33];
        case 2  % PD-TDPF
            J = {[J11 J12; J21 J22],J33};
        case 3  % FD-TDPF
            J = {J11,J22,J33};
        case 4  % Conventional PF or SD-TDPF
            J = [J11 J12; J21 J22];
        case 5  % Fast-decoupled conventional PF
            J = {J11,J22};
        otherwise
            error('Unrecognized Jacobian type.');
    end
end


