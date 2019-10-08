% Thu  8 Jun 13:03:02 CEST 2017
%% class for fitting models to at-a-station time series of tidal elevation
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
classdef River_Tide_Empirical < handle
	properties
		c     = struct();
		cs    = struct();
		c0    = [];
		model = 1;
		nlflag = true;
		opt   = struct('MaxFunEvals',1e3,'Display','off');
		fr
		fD
		fk
		fz
	end % properties
	methods
		function obj = River_Tide_Empirical(model)
			obj.model = model;
			% select model
			obj.rt_model();
		end
	end % methods
end % River_Tide_Empirical

