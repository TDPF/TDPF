%% function evalJacobian: Evaluate Jacobian Matrix
% Evaluates the Jacobian matrix for temperature-dependent power flow. The
% structure of the Jacobian(s) returned depends on the type of TDPF being
% performed, as specified by the user when the function is called.
% 
% SYNTAX:
%   J = evalJacobian(type,Sets,V,delta,T,G,B,branches)
%
% INPUTS:
%   type =      Indicates type of Jacobian(s) to return, based on type of
%               TDPF being performed. Options are:
%               1 -- Jacobian matrix for fully-coupled TDPF (FC-TDPF)
%               2 -- Jacobian submatrices J1-J4 and J9 for
%                    partially-decoupled TDPF (PD-TDPF)
%               3 -- Jacobian submatrices J1, J4, and J9 for
%                    fast-decoupled TDPF (FD-TDPF)
%               4 -- Jacobian matrix for conventional power flow (J1-J4);
%                    needed for conventional power flow or for
%                    sequentially decoupled TDPF (SD-TDPF)
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
%               these requirements. For 'type' == 4, this argument is
%               ignored (and may safely be replaced by an empty matrix).
%
% OUTPUTS:
%   J =         Jacobian matrix (or matrices) with structure dependent on
%               the input variable 'type' as follows:
%               1 -- J is a single matrix representing the entire Jacobian
%                    J1-J9 for use with FC-TDPF
%               2 -- J is a cell structure containing two Jacobian
%                    submatrices: J1-J4 and J9 (in that order).
%                    For use with PD-TDPF.
%               3 -- J is a cell structure containing three Jacobian
%                    submatrices: J1, J4 and J9 (in that order).
%                    For use with FD-TDPF.
%               4 -- J is single matrix representing the Jacobian J1-J4
%                    for conventional power flow. For use with conventional
%                    power flow or SD-TDPF.
%
% COMMENTS:
%
function J = evalJacobian(type,sets,V,delta,T,G,B,branch)
    %% Compute required vector/matrix dimensions
    N = length(V);          % Total number of buses
    M = length(sets.V);     % Number of PQ buses
    L = length(branch.id);  % Total number of branches
    R = length(sets.T);     % Number of temperature-dependent branches

    %% Compute some necessary quantities for branch elements
    % (Required for all cases except type==4)
    if type ~= 4
       % Compute series conductance of branch elements
       branch.g = branch.R ./ (branch.R.^2 + branch.X.^2);
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
    
    % The full Jacobian matrix is as follows:
    %
    % / J1 = dP/dd      J2 = dP/dV      J7 = dP/dT  \
    % | (N-1)x(N-1)     (N-1)x(M)       (N-1)x(R)   |
    % |                                             |
    % | J3 = dQ/dd      J4 = dQ/dV      J8 = dQ/dT  |
    % | (M)x(N-1)       (M)x(M)         (M)x(R)     |
    % |                                             |
    % | J5 = dH/dd      J6 = dH/dV      J9 = dH/dT  |
    % \ (R)x(N-1)       (R)x(M)         (R)x(R)     /
    %
    % (We may want to renumber later...)
    
    % J1 = dP/dd    (N-1) x (N-1)
    % Needed for all types
        % Find indices of all nonzero YBus elements + diagonals where:
        %   The row    exists in sets.P
        %   The column exists in sets.delta
        [jRow, jCol]  = find( nzY( sets.P, sets.delta ) );
        
        % NOTE: The indexing is now relative to J1 rather than YBus. This
        % is important when translating back to V and delta values!
        
        % Storage for J1 elements
        jElem = zeros( size(jRow) );
        
        % Generate each nonzero element of J1
        for count = 1:length(jRow)
            % Find indices with respect to J1
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
        
        % Generate J1
        J1 = sparse(jRow, jCol, jElem, N-1, N-1);
    
    
    % J2 = dP/dV    (N-1) x (M)
    if any( type == [1 2 4] )
        % Find indices of all nonzero YBus elements + diagonals where:
        %   The row    exists in sets.P
        %   The column exists in sets.V
        [jRow, jCol]  = find( nzY( sets.P, sets.V ) );
        
        % NOTE: The indexing is now relative to J2 rather than YBus. This
        % is important when translating back to V and delta values!
        
        % Storage for J2 elements
        jElem = zeros( size(jRow) );
        
        % Generate each nonzero element of J2
        for count = 1:length(jRow)
            % Find indices with respect to J2
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
        
        % Generate J1
        J2 = sparse(jRow, jCol, jElem, N-1, M);
    end
    
    
    % J3 = dQ/dd    (M) x (N-1)
    if any( type == [1 2 4] )
        % Find indices of all nonzero YBus elements + diagonals where:
        %   The row    exists in sets.Q
        %   The column exists in sets.delta
        [jRow, jCol]  = find( nzY( sets.Q, sets.delta ) );
        
        % NOTE: The indexing is now relative to J3 rather than YBus. This
        % is important when translating back to V and delta values!
        
        % Storage for J3 elements
        jElem = zeros( size(jRow) );
        
        % Generate each nonzero element of J3
        for count = 1:length(jRow)
            % Find indices with respect to J3
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
        
        % Generate J3
        J3 = sparse(jRow, jCol, jElem, M, N-1);
    end
    
    
    % J4 = dQ/dV	(M) x (M)
    % Needed for all types
        % Find indices of all nonzero YBus elements + diagonals where:
        %   The row    exists in sets.Q
        %   The column exists in sets.V
        [jRow, jCol]  = find( nzY( sets.Q, sets.V ) );
        
        % NOTE: The indexing is now relative to J4 rather than YBus. This
        % is important when translating back to V and delta values!
        
        % Storage for J4 elements
        jElem = zeros( size(jRow) );
        
        % Generate each nonzero element of J4
        for count = 1:length(jRow)
            % Find indices with respect to J4
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
        
        % Generate J4
        J4 = sparse(jRow, jCol, jElem, M, M);
    
    
    % J5 = dH/dd    (R) x (N-1)
    if type == 1
        % Find indices of all nonzero branch-to-bus connections where
        %   The row    exists in sets.H
        %   The column exists in sets.delta
        [jRow, jCol]  = find( nzL( sets.H, sets.delta ) );
        
        % NOTE: The indexing is now relative to J5.
        
        % Storage for J5 elements
        jElem = zeros( size(jRow) );
        
        % Generate each nonzero element of J5
        for count = 1:length(jRow)
            % Find indices with respect to J5
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
        
        % Generate J5
        J5 = sparse(jRow, jCol, jElem, R, N-1);
    end
    
    
    % J6 = dH/dV    (R) x (M)
    if type == 1
        % Find indices of all nonzero branch-to-bus connections where
        %   The row    exists in sets.H
        %   The column exists in sets.V
        [jRow, jCol]  = find( nzL( sets.H, sets.V ) );
        
        % NOTE: The indexing is now relative to J6.
        
        % Storage for J6 elements
        jElem = zeros( size(jRow) );
        
        % Generate each nonzero element of J6
        for count = 1:length(jRow)
            % Find indices with respect to J6
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
        
        % Generate J6
        J6 = sparse(jRow, jCol, jElem, R, M);
    end
    
    % J7 = dP/dT    (N-1) x (R)
    if type == 1
        % Find indices of all nonzero bus-to-branch connections where
        %   The row    exists in sets.P
        %   The column exists in sets.T
        % Note that this requires using the transpose of nzL
        [jRow, jCol]  = find( nzL( sets.T, sets.P )' );
        
        % NOTE: The indexing is now relative to J7.
        
        % Storage for J7 elements
        jElem = zeros( size(jRow) );
        
        % Generate each nonzero element of J7
        for count = 1:length(jRow)
            % Find indices with respect to J7
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
        
        % Generate J7
        J7 = sparse(jRow, jCol, jElem, N-1, R);
    end
    
    % J8 = dQ/dT    (M) x (R)
    if type == 1
        % Find indices of all nonzero bus-to-branch connections where
        %   The row    exists in sets.Q
        %   The column exists in sets.T
        % Note that this requires using the transpose of nzL
        [jRow, jCol]  = find( nzL( sets.T, sets.Q )' );
        
        % NOTE: The indexing is now relative to J8.
        
        % Storage for J8 elements
        jElem = zeros( size(jRow) );
        
        % Generate each nonzero element of J8
        for count = 1:length(jRow)
            % Find indices with respect to J8
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
        
        % Generate J8
        J8 = sparse(jRow, jCol, jElem, M, R);
    end
    
    
    % J9 = dH/dT    (R) x (R)
    if any( type == [1 2] )
        % Row, column, and element vectors for J9
        % NOTE: J9 has only diagonal entries
        jRow = 1:R; jCol = 1:R;
        jElem = zeros( size(jRow) );
        
        % Generate elements of J9
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
        
        % Generate J9
        J9 = sparse(jRow, jCol, jElem, R, R);
    elseif type == 3
        % For FD-TDPF, J9 is approximated by the identity matrix
        J9 = speye(R);
    end
        
    %% Return Type-Dependent Output
    switch type
        case 1  % FC-TDPF
            J = [J1 J2 J7; J3 J4 J8; J5 J6 J9];
        case 2  % PD-TDPF
            J = {[J1 J2; J3 J4],J9};
        case 3  % FD-TDPF
            J = {J1,J4,J9};
        case 4  % Conventional PF or SD-TDPF
            J = [J1 J2; J3 J4];
        otherwise
            error('Unrecognized Jacobian type.');
    end
end




