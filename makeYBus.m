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

        % Find all branches with 'from' terminations at bus k
        ff = branch.id( branch.from == k);
        
        % Find all branches with 'to' terminations at bus k
        tt = branch.id( branch.to == k );
        
        % Compute Y(k,k)
        ytt = branch.y(tt) + 1j .* branch.B(tt) ./ 2;
        yff = (branch.y(ff) + 1j .* branch.B(ff) ./ 2 ) ...
              ./ (branch.tap(ff) .* conj(branch.tap(ff)));
        diagElem(kk) = sum( ytt ) + sum( yff ) + ...
                       bus.G(kk) + 1j * bus.B(kk);  % Shunt admittance
    end
    
    % Compute off-diagonal elements
    % (Only need to compute these at nonzero locations in YBus)
    offdiagRow = [branch.from, branch.to];
    offdiagCol = [branch.to, branch.from];
    
    % k = 'from' bus, n = 'to' bus
    % 'from' - 'to' -> If [k] = ... + yft [k,n] * Vt [n]
    yft = -branch.y ./ conj(branch.tap);
    
    % 'to' - 'from' -> It [n] = ... + ytf [n,k] * Vf [k]
    ytf = -branch.y ./ branch.tap;
    
    % Assemble off-diagonal elements
    offdiagElem = [yft, ytf];
    
    % Construct YBus
    Y = sparse([diagIdx, offdiagRow], [diagIdx, offdiagCol], ...
               [diagElem, offdiagElem],N,N);
    
    % Compute Real/Imaginary parts, Magnitude, Angle
    G = real(Y);
    B = imag(Y);
    YMag = abs(Y);
    YAng = angle(Y);
end
