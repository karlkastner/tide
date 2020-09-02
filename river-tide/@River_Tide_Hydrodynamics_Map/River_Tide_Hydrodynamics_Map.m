% Sat 14 Oct 16:39:01 CEST 2017
%% container class to store multiple river-tide scenarios
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
classdef River_Tide_Hydrodynamics_Map < Compute_Map
	properties
	end % properties
	methods
		function obj = River_Tide_Hydrodynamics_Map(varargin)
		         obj = obj@Compute_Map(varargin{:});
		end % constructor
	end % methods
end % class River_Tide_Map

