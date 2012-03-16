%% function makeYBus: Create bus adittance matrix
% Generates the bus admittance matrix YBus given a set of bus and branch
% data in structure format. Uses the current R, X, G, and B values
% specified in the input data (not the reference values).
% 
% SYNTAX:
%   [Y,G,B,YMag,YAng] = makeYBus(bus,branch)
%
% INPUTS:
%   bus =       Set of bus data from makeBusStruct()
%   branch =    Set of branch data from makeBranchStruct()
%
% OUTPUTS:
%   Y =         Bus admittance matrix (complex)
%   G =         Bus conductance matrix (real part of 'Y')
%   B =         Bus susceptance matrix (imaginary part of 'Y')
%   YMag =      Element-wise magnitudes of bus admittance matrix
%   YAng =      Element-wise angles of bus admittance matrix 
%
% COMMENTS:
%   1. Note that 'YAng' is in radians, not degrees.
%
function [Y,G,B,YMag,YAng] = makeYBus(bus,branch)
	% Get YBus dimensions
	N = length(bus.id);         % Number of buses
    L = length(branch.id);      % Number of branches
   
    % For off-nominal turns ratios, we define:
    %   Vt = voltage at 'to' bus (Z bus)
    %   It = current injected into branch at 'to' bus (Z bus)
    %   Vf = voltage at 'from' bus (tap bus)
    %   If = current injected into branch at 'from' bus (tap bus)
    %   y  = series admittance of branch 'from'-'to'
    %   b  = shunt susceptance (line charging) of branch 'from'-'to'
    %   a  = complex turns ratio for branch 'from'-'to' 
    %
    % For any given branch we then need the ABCD parameter relationships:   
    %   / If \   / A  B \   / Vf \   / yff  yft \   / Vf \
    %	|    | = |      | * |    | = |          | * |    | 
    %	\ It /   \ C  D /   \ Vt /   \ ytf  ytt /   \ Vt /
    %
    % Which are given by:
    %   A = yff = (y + jb/2) / (aa*) = ytt / (aa*)
    %   B = yft = -y / a*
    %   C = ytf = -y / a
    %   D = ytt + jb/2
    % 
    % Therefore, for I = YBus * V to work...
    %   YBus(k,k) = sum( ytt for k == t ) + sum( yff for k == f )
    %   YBus(k,n) = yft if kn == ft, else ytf if kn == tf
    
    % Compute line admittances
    branch.y = 1 ./ (branch.R + 1j * branch.X);
    
    % Compute diagonal elements
    diagIdx = 1:N; 
    diagElem = zeros(1,N);
    for kk = diagIdx
        % Get bus IDs
        k = bus.id(kk);
        
        % Find all branches with 'to' terminations at bus k
        tt = branch.id( branch.to == k );

        % Find all branches with 'from' terminations at bus k
        ff = branch.id( branch.from == k);

        % Compute Y(k,k)
        ytt = branch.y(tt) + 1j .* branch.B(tt) ./ 2;
        yff = (branch.y(ff) + 1j .* branch.B(ff) ./ 2 ) ...
              ./ (branch.tap(ff) .* conj(branch.tap(ff)));
        diagElem(kk) = sum( ytt ) + sum( yff ) + ...
                       bus.G(kk) + 1j * bus.B(kk);
    end
    
    % Compute off-diagonal elements
    % (Only need to compute these at nonzero locations in YBus)
    rowIdx = [branch.from, branch.to];
    colIdx = [branch.to, branch.from];
    offdiagElem = zeros(1,2*L);
    
    % NOTE: k = 'from' bus, n = 'to' bus
    
    % First, compute (row k, col n) entry 
    %   If [k] = ... + yft [k,n] * Vt [n]
    offdiagElem(1:L) = -branch.y ./ conj(branch.tap);
    
    % Next, compute (row n, col k) entry
    %   It [n] = ... + ytf [n,k] * Vf [k]
    offdiagElem((1:L)+L) = -branch.y ./ branch.tap;
    
    % Construct YBus
    Y = sparse([diagIdx, rowIdx], [diagIdx, colIdx], ...
               [diagElem, offdiagElem],N,N);
        
        
        % Off-diagonal elements
%         else
%             % Find all branchs from k to n (i.e. kn = ft)
%             ft = branch.id( (branch.from == k) & (branch.to == n) );
% 
%             % Find all branchs from n to k (i.e. kn = tf)
%             tf = branch.id( (branch.to == k) & (branch.from == n) );
% 
%             % Compute Y(k,n)
%             Y(kk,nn) = sum( -branch.y(ft) ./ conj(branch.tap(ft)) ) ...
%                      + sum( -branch.y(tf) ./ branch.tap(tf) );
%         end
    
    % Compute Real/Imaginary parts, Magnitude, Angle
    G = real(Y);
    B = imag(Y);
    YMag = abs(Y);
    YAng = angle(Y);
end
