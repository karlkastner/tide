% Sat 14 Oct 16:39:01 CEST 2017
%% container class to store individual river tide scenarios
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
classdef River_Tide_Map < handle
	properties
%		model_str = 'waveq';
		recompute = false;
		recflag   = false;
		filename
		rt
	end % properties
	methods
		function obj = River_Tide_Map(filename)
			obj.filename = filename;
		end % constructor

		function obj = init(obj)
			if (0 ~= exist(obj.filename,'file'))
				load(obj.filename,'rt');
				obj.rt = rt;
			else
				disp('file not existing, creating new container')
				obj.rt = containers.Map();
			end
		end % filename

		function obj = save(obj)
			if (obj.recflag)
				rt = obj.rt;
				save(obj.filename,'rt');
				% reset recflag
				obj.recflag = false;
			end
		end % filename

		% Sat 14 Oct 16:41:04 CEST 2017
		function [obj] = run(obj,val_C,pflag)
				val_C = {obj,val_C{:}};
				if (nargin()>2 && pflag)
					iterate_cell(@plot,val_C);
				else
					iterate_cell(@fun,val_C);
				end
		end % run
	end % methods
end % River_Tide_Map

