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
% TO DO:
%   * Implement the entries of the Jacobian!
%

function J = evalJacobian(type,sets,V,delta,T,G,B,branch)
    %% Compute required vector/matrix dimensions
    N = length(V);          % Total number of buses
    M = length(sets.V);     % Number of PQ buses
    L = length(T);          % Total number of branches
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
        J1 = zeros(N-1,N-1);
        for ii = 1:(N-1)
            for kk = 1:(N-1)
                i = sets.P(ii);
                k = sets.delta(kk);
                % Same bus
                if (k == i)
                    J1(ii,kk) = ...
                        V(i) * sum( V .* ( ...
                        -G(i,:) .* sin(delta(i) - delta) ...
                        +B(i,:) .* cos(delta(i) - delta) ...
                        )) - V(i)^2 * B(i,i);
                % Different buses
                else
                    J1(ii,kk) = ...
                        V(i)*V(k)*( ...
                        G(i,k)*sin(delta(i)-delta(k)) - ...
                        B(i,k)*cos(delta(i)-delta(k)) );
                end
            end
        end
        
    % J2 = dP/dV    (N-1) x (M)
    if any( type == [1 2 4] )
        J2 = zeros(N-1,M);
        for ii = 1:(N-1)
            for kk = 1:M
                i = sets.P(ii);
                k = sets.V(kk);
                % Same bus
                if (k == i)
                    J2(ii,kk) = V(i) * G(i,i) + ...
                        sum( V .* ( ...
                        G(i,:) .* cos(delta(i) - delta) + ...
                        B(i,:) .* sin(delta(i) - delta) ...
                        ));
                % Different buses
                else
                    J2(ii,kk) = V(i) * ( ...
                        G(i,k) * cos(delta(i) - delta(k)) + ...
                        B(i,k) * sin(delta(i) - delta(k)) );
                end
            end
        end
    end
        
    % J3 = dQ/dd    (M) x (N-1)
    if any( type == [1 2 4] )
        J3 = zeros(M,N-1);
        for ii = 1:M
            for kk = 1:(N-1)
                i = sets.Q(ii);
                k = sets.delta(kk);
                % Same bus
                if (k == i)
                    J3(ii,kk) = ...
                        V(i) * sum( V .* ( ...
                        G(i,:) .* cos(delta(i) - delta) + ...
                        B(i,:) .* sin(delta(i) - delta) ...
                        )) - V(i)^2 * G(i,i);
                % Different buses
                else
                    J3(ii,kk) = ...
                        V(i)*V(k)*( ...
                        -G(i,k)*cos(delta(i)-delta(k)) - ...
                        B(i,k)*sin(delta(i)-delta(k)) );
                end
            end
        end
    end
    
    % J4 = dQ/dV	(M) x (M)
    % Needed for all types
        J4 = zeros(M,M);
        for ii = 1:M
            for kk = 1:M
                i = sets.Q(ii);
                k = sets.V(kk);
                % Same bus
                if (k == i)
                    J4(ii,kk) = -V(i) * B(i,i) + ...
                        sum( V .* ( ...
                        G(i,:) .* sin(delta(i) - delta) - ...
                        B(i,:) .* cos(delta(i) - delta) ...
                        ));
                % Different buses
                else
                    J4(ii,kk) = V(i) * ( ...
                        G(i,k) * sin(delta(i) - delta(k)) - ...
                        B(i,k) * cos(delta(i) - delta(k)) );
                end
            end
        end
    
    % J5 = dH/dd    (R) x (N-1)
    if type == 1
        J5 = zeros(R,N-1);
        for ii = 1:R
            for kk = 1:(N-1)
                % Indices
                ij = sets.H(ii);
                i = branch.from(ij);
                j = branch.to(ij);
                k = sets.delta(kk);
                        
                % Voltage magnitudes and angles
                % 'from' bus must be modified by tap ratio.
                Vi = V(i) / abs(branch.tap(ij));
                deltai = delta(i) - angle(branch.tap(ij));

                % 'to' bus is used directly
                Vj = V(j);
                deltaj = delta(j);
                
                if (k == i)
                    J5(ii,kk) = -2 * branch.R_therm(ij) * ...
                                branch.g(ij) * ...
                                Vi * Vj * sin(deltai - deltaj);
                elseif (k == j)
                    J5(ii,kk) = 2 * branch.R_therm(ij) * ...
                                branch.g(ij) * ...
                                Vi * Vj * sin(deltai - deltaj);
                else
                    J5(ii,kk) = 0;
                end
            end
        end
    end
    
    % J6 = dH/dV    (R) x (M)
    if type == 1
        J6 = zeros(R,M);
        for ii = 1:R
            for kk = 1:M
                % Indices
                ij = sets.H(ii);
                i = branch.from(ij);
                j = branch.to(ij);
                k = sets.V(kk);
                
                % Voltage magnitudes and angles
                % 'from' bus must be modified by tap ratio.
                Vi = V(i) / abs(branch.tap(ij));
                deltai = delta(i) - angle(branch.tap(ij));

                % 'to' bus is used directly
                Vj = V(j);
                deltaj = delta(j);
                
                if (k == i)
                    J6(ii,kk) = -2 * branch.R_therm(ij) * ...
                                branch.g(ij) * ...
                                ( Vi - Vj * cos(deltai - deltaj) );
                elseif (k == j)
                    J6(ii,kk) = -2 * branch.R_therm(ij) * ...
                                branch.g(ij) * ...
                                ( Vj - Vi * cos(deltai - deltaj) );
                else
                    J6(ii,kk) = 0;
                end
            end
        end
    end
    
    % J7 = dP/dT    (N-1) x (R)
    if type == 1
        J7 = zeros(N-1,R);
        for ii = 1:(N-1)
            for kk = 1:R
                % Indices
                i = sets.P(ii);
                kn = sets.T(kk);
                k = branch.from(kn);
                n = branch.to(kn);
                
                if (k == i || n == i)
                    % Compute dg/dT, db/dT
                    dgdT = ( branch.X(kn)^2 - branch.R(kn)^2 ) / ...
                           ( branch.R(kn)^2 + branch.X(kn)^2 )^2 * ...
                             branch.R_ref(kn) / ...
                           ( branch.T_ref(kn) + branch.T_f(kn) );
                    dbdT = ( 2 * branch.X(kn) * branch.R(kn) ) / ...
                           ( branch.R(kn)^2 + branch.X(kn)^2 )^2 * ...
                             branch.R_ref(kn) / ...
                           ( branch.T_ref(kn) + branch.T_f(kn) );
                end
                if (k == i)
                    % Voltage magnitudes and angles
                    % 'from' bus must be modified by tap ratio.
                    Vi = V(i) / abs(branch.tap(kn));
                    deltai = delta(i) - angle(branch.tap(kn));

                    % 'to' bus is used directly
                    Vn = V(n);
                    deltan = delta(n);
                    
                    % Compute dP/dT
                    J7(ii,kk) = ( Vi^2 - Vi * Vn * ...
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
                    J7(ii,kk) = ( Vi^2 - Vi * Vk * ...
                                  cos(deltai - deltak) ) * dgdT - ...
                                Vi * Vk * ...
                                  sin(deltai - deltak) * dbdT;
                else
                	J7(ii,kk) = 0;
                end
            end
        end
    end
    
    % J8 = dQ/dT    (M) x (R)
    if type == 1
        J8 = zeros(M,R);
        for ii = 1:M
            for kk = 1:R
                % Indices
                i = sets.Q(ii);
                kn = sets.T(kk);
                k = branch.from(kn);
                n = branch.to(kn);
                
                if (k == i || n == i)
                    % Compute dg/dT, db/dT
                    dgdT = ( branch.X(kn)^2 - branch.R(kn)^2 ) / ...
                           ( branch.R(kn)^2 + branch.X(kn)^2 )^2 * ...
                             branch.R_ref(kn) / ...
                           ( branch.T_ref(kn) + branch.T_f(kn) );
                    dbdT = ( 2 * branch.X(kn) * branch.R(kn) ) / ...
                           ( branch.R(kn)^2 + branch.X(kn)^2 )^2 * ...
                             branch.R_ref(kn) / ...
                           ( branch.T_ref(kn) + branch.T_f(kn) );
                end
                if (k == i)
                    % Voltage magnitudes and angles
                    % 'from' bus must be modified by tap ratio.
                    Vi = V(i) / abs(branch.tap(kn));
                    deltai = delta(i) - angle(branch.tap(kn));

                    % 'to' bus is used directly
                    Vn = V(n);
                    deltan = delta(n);
                    
                    % Compute dQ/dT
                    J8(ii,kk) = -Vi * Vn * ...
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
                    J8(ii,kk) = -Vi * Vk * ...
                                  sin(deltai - deltak) * dgdT + ...
                                ( -Vi^2 + Vi * Vk * ...
                                  cos(deltai - deltak) ) * dbdT;
                else
                	J8(ii,kk) = 0;
                end
            end
        end
    end
    
    % J9 = dH/dT    (R) x (R)
    if any( type == [1 2] )
        J9 = zeros(R,R);
        % J9 has only diagonal entries
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
        	J9(kk,kk) = 1 - branch.R_therm(kn) * ...
                        ( Vk^2 + Vn^2 - ...
                          2 * Vk * Vn * cos(deltak - deltan) ) * dgdT;
        end
    elseif type == 3
        % For FD-TDPF, J9 is approximated by the identity matrix
        J9 = eye(R);
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




