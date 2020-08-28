% Mon 24 Aug 10:24:24 +08 2020
%
%% container class to store multiple river-tide morphodyanics scenarios
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
classdef River_Tide_Morphodynamics_Map < Compute_Map
	properties
	end % properties
	methods
		function obj = River_Tide_Morphodynamics_Map(varargin)
		         obj = obj@Compute_Map(varargin{:});
		end % River_Tide_Morphodynamics_Map
	end % methods
end % class River_Tide_Morphodynamics_Map

