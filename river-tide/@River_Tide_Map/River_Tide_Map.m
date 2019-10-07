% Sat 14 Oct 16:39:01 CEST 2017
%% container class to store individual river tide scenarios
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

