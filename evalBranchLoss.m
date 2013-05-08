%% function evalBranchLoss: Evaluate Branch Losses
% Evaluates power transfer, loss, and (optionally) temperature for each
% system branch given system bus voltage data, branch connections, and
% branch impedances.
% 
% SYNTAX:
%   branch = evalBranchLoss(bus, branch, varargin)
%
% INPUTS:
%	bus =       Bus structure (as returned from importCaseData())
%               containing current state data (voltage magnitudes and
%               angles)
%	branch =    Branch structure (as returned from importCaseData())
%   varargin =  (Optional) Additional arguments passed as name-value pairs
%               (see OPTIONAL INPUTS below)
%
% OPTIONAL INPUTS:
%   The following Optional inputs may be passed as name-value pairs
%   following 'branch':
%
%   'updateTemps', [val]    Specify whether branch temperatures and
%                           resistances should be updated to match the
%                           given voltage states and computed losses.
%                           Default = FALSE
%   'tempTol', [val]        Specify temperature mismatch tolerance, in
%                           per-unit. (Used if updateTemps == TRUE.)
%                           Default = 1e-08
%
% OUTPUTS:
%   branch =    Modified branch structure with added fields:
%               .S_from     Complex power into the branch at the 'from' bus
%               .S_to       Complex power into the branch at the 'to' bus
%               .S_loss     Complex power loss within the branch
%                           (S_loss = S_from + S_to)
%               .loading    Per-unit branch loading (measured at 'from'
%                           bus)
%               If updateTemps == TRUE, the the branch fields 'T' and 'R'
%               may also be modified.
%
% COMMENTS:
%   1. All outputs are in per-unit consistent with the per-unit bases for
%      the bus and branch data.
%   2. If updateTemps == TRUE, then the function will iteratively update
%      the temperatures and resistances of temperature-dependent branches 
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
function branch = evalBranchLoss(bus, branch, varargin)
    %% Setup
    % Default control parameters
	updateTemps = false;	% Update temperatures/resistances? (T/F)
    tempTol = 1e-08;        % Tolerance for temperature update
    
	% Parse optional arguments
	while ~isempty(varargin)
		name = lower( varargin{1} );
		switch name
            % Control parameters
			case {'updatetemps'}
				updateTemps = varargin{2};
			case {'temptol'}
				tempTol = varargin{2};
		end
		
		% Clear these two entries
		varargin(1:2) = [];
    end
    
    %% Get Voltages
    % Bus indices for branch terminals
    [~, ff] = ismember(branch.from, bus.id);
    [~, tt] = ismember(branch.to,   bus.id);
    
    % Bus voltages for branch terminals
    % (Corrected at from bus for off-nominal turns ratios)
    Vf = bus.V_mag(ff) .* ...
        ( cosd(bus.V_angle(ff)) + 1j*sind(bus.V_angle(ff)) ) ./ ...
        branch.tap;
    Vt = bus.V_mag(tt) .* ...
        ( cosd(bus.V_angle(tt)) + 1j*sind(bus.V_angle(tt)) );
    
    %% Compute Powers/Losses
    % Series impedances
    y = 1 ./ (branch.R + 1j*branch.X);
    
    % Complex powers
    Sf = Vf .* conj(y .* (Vf - Vt));
    St = Vt .* conj(y .* (Vt - Vf));
    
    % Complex loss
    Sloss = Sf + St;
    
    %% Update Temperatures/Resistances
    % If applicable...
    if updateTemps
        % Get indices of temperature dependent branches
        td = find( branch.type > 0 );
        
        % Initial temperature mismatch
        mm = branch.T(td) - ...
            (branch.T_amb(td) + real(Sloss(td)) .* branch.R_therm(td));
        
        % Loop to temperature convergence
        while norm(mm,Inf) > tempTol
            % Maximum existing mismatch
            maxmm = norm(mm,Inf);
            
            % Track existing temperature estimate
            oldTemp = branch.T(td);
            
            % Update temperature estimate
            branch.T(td) = ...
                branch.T_amb(td) + real(Sloss(td)) .* branch.R_therm(td);
            
            % Update branch resistances
            branch.R(td) = branch.R_ref(td) .* ...
                ( (branch.T(td) + branch.T_f(td)) ./ ...
                (branch.T_ref(td) + branch.T_f(td)) );
            
            % Recompute series impedances
            y(td) = 1 ./ (branch.R(td) + 1j*branch.X(td));
            
            % Recompute complex powers and losses
            Sf(td) = Vf(td) .* conj(y(td) .* (Vf(td) - Vt(td)));
            St(td) = Vt(td) .* conj(y(td) .* (Vt(td) - Vf(td)));
            Sloss(td) = Sf(td) + St(td);
            
            % Compute new mismatch
            mm = oldTemp - branch.T(td);
            
            % Emergency hatch - for if mismatches are not converging
            if norm(mm,Inf) > maxmm
                warning([ ...
                    'Temperature data is not converging; temperature ' ...
                    'update to branch data aborted.']);
                break;
            end
        end
    end
    
    %% Compute Loading
    % Measured at 'from' bus
    Loading = abs(Sf) ./ branch.rating;
    
    %% Return
    % Assign data back to branch structure
    branch.S_from = Sf;
    branch.S_to = St;
    branch.S_loss = Sloss;
    branch.loading = Loading;
end


