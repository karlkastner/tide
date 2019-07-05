% Thu  8 Jun 13:03:02 CEST 2017
%% class for fitting models to at-a-station time series of tidal elevation
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

