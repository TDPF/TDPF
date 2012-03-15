%% Test the new indexing method for the Jacobian matrix
%% Data Import
branch = makeBranchStruct('IEEE_30_BRANCH.csv');
[bus,N,PQ,PV,Slack] = makeBusStruct('IEEE_30_BUS.csv');

%% Compute Index Sets
% Let:
%   N = Number of buses
%   M = Number of PQ buses
%   R = Number of temperature dependent lines

% Num. mismatches...
%   N-1     Power (P)
%   M       Reactive Power (Q)
%   R       Temperature difference (H)

% Num. variables...
%   N-1     Voltage angles (d)
%   M       Voltage magnitudes (V)
%   R       Temperatures (T)

% The following are the sets of bus and branch indices representing
% the variable and mismatch equation sets:
Sets.P = [PV PQ];
Sets.Q = PQ;
Sets.H = branch.id( branch.R > 0 & branch.T_rrise > 0 & ...
    branch.P_rloss > 0);
Sets.V = PQ;
Sets.d = [PV PQ];
Sets.T = Sets.H;

% Get appropriate lengths
% N already exists
M = length(Sets.V);
R = length(Sets.T);

%% Create 'Analytical' Jacobian Matrix
% Derivatives are...
%
% / dP/dd           dP/dV           dP/dT       \
% | (N-1)x(N-1)     (N-1)x(M)       (N-1)x(R)   |
% |                                             |
% | dQ/dd           dQ/dV           dQ/dT       |
% | (M)x(N-1)       (M)x(M)         (M)x(R)     |
% |                                             |
% | dH/dd           dH/dV           dH/dT       |
% \ (R)x(N-1)       (R)x(M)         (R)x(R)     /

% The following code uses a cell array to create a hypothetical listing
% of Jacobian elements for this system
J = cell( (N-1+M+R), (N-1+M+R) );

% J1 = dP/dd
for k = 1:(N-1)
    for n = 1:(N-1)
        % Offset of matrix J1 is (0,0)
        J{k,n} = ['dP' int2str(Sets.P(k)) '/dd' int2str(Sets.d(n))];
    end
end

% J2 = dP/dV
for k = 1:(N-1)
    for n = 1:M
        % Offset of matrix J2 is (0,N-1)
        J{k,n+(N-1)} = ['dP' int2str(Sets.P(k)) '/dV' int2str(Sets.V(n))];
    end
end

% J3 = dQ/dD
for k = 1:M
    for n = 1:(N-1)
        % Offset of matrix J3 is (N-1,0)
        J{k+(N-1),n} = ['dQ' int2str(Sets.Q(k)) '/dd' int2str(Sets.d(n))];
    end
end

% J4 = dQ/dV
for k = 1:M
    for n = 1:M
        % Offset of matrix J4 is (N-1,N-1)
        J{k+(N-1),n+(N-1)} = ...
            ['dQ' int2str(Sets.Q(k)) '/dV' int2str(Sets.V(n))];
    end
end

% J5 = dH/dd
for k = 1:R
    for n = 1:(N-1)
        % Offset of matrix J4 is (N-1+M,0)
        J{k+(N-1+M),n} = ...
            ['dH' int2str(Sets.H(k)) '/dd' int2str(Sets.d(n))];
    end
end

% J6 = dH/dV
for k = 1:R
    for n = 1:M
        % Offset of matrix J4 is (N-1+M,N-1)
        J{k+(N-1+M),n+(N-1)} = ...
            ['dH' int2str(Sets.H(k)) '/dV' int2str(Sets.V(n))];
    end
end

% J7 = dP/dT
for k = 1:(N-1)
    for n = 1:R
        % Offset of matrix J4 is (0,N-1+M)
        J{k,n+(N-1+M)} = ...
            ['dP' int2str(Sets.P(k)) '/dT' int2str(Sets.T(n))];
    end
end

% J8 = dQ/dT
for k = 1:M
    for n = 1:R
        % Offset of matrix J4 is (N-1,N-1+M)
        J{k+(N-1),n+(N-1+M)} = ...
            ['dQ' int2str(Sets.Q(k)) '/dT' int2str(Sets.T(n))];
    end
end

% J9 = dH/dT
for k = 1:R
    for n = 1:R
        % Offset of matrix J4 is (N-1+M,N-1+M)
        J{k+(N-1+M),n+(N-1+M)} = ...
            ['dH' int2str(Sets.H(k)) '/dT' int2str(Sets.T(n))];
    end
end

% Using this arrangement of Jacobian matrix entries, the vector of state
% variables [with bus domain order in brackets] becomes:
%	[ d[PQ], V[PV,PQ], T ]
%
% Thus, the various states can be updated at each iteration by extracting
% the variables in order from the state vector and plugging them back into
% the original vectors (e.g. voltage magnitude vector) at the appropriate
% indices (also in order).
