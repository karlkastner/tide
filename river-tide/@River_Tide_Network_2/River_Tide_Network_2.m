% Fri 22 Feb 14:41:57 CET 2019
%
%% tide in a fluvial delta channel network, extension of 1D river tide
%% the network is a directed graph
%% TODO convert from trig-to exponential form
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
classdef River_Tide_Network_2 < handle
	properties
		% array of River_Tide objects, one for each channel in the network
		rt

		% coupling conditions at channel junctions
		junction_condition = {};

		division_rule
		junction_Qs
		opt = struct('ignorertfordivision',false);

		% change over time
		evolution = struct('t', [], 'zb', []);
		
		% seasonal hydrograph
		season = struct('iorder',1,'Qmin',[],'Qmax',[]);
	
		% temporary storage, for reuse of initial conditions
		tmp	
	end % properties 
	methods
		function obj = River_Tide_Network_2(rt)
			if (nargin()>0)
				obj.rt = rt;
			end
			obj.division_rule = @obj.sediment_division_geometric;
		end

		% pseudo members
		function Q = Q(obj,id)
			Q = [];
			for idx=1:length(obj.rt)
				Q = [Q; obj.rt(idx).Q(id)];
			end	
		end
		
		function z = z(obj,id)
			z = [];
			for idx=1:length(obj.rt)
				z = [z; obj.rt(idx).z(id)];
			end	
		end	

		function init(obj)
			for idx=1:length(obj.rt)
				obj.rt(idx).init();
			end
		end
	end % methods
end % class River_Tide_Network_2

