% Sun  9 Aug 12:34:57 +08 2020

classdef River_Tide_Morphodynamics < River_Tide
	properties
		sediment = struct( 'd_mm', 0.2, ... % grain diameter
				   'p',   0.6, ...  % packing density
				   'rho', 2650 ...  % material density
				 );
		% order of accuracy for time stepping
		norder = 1;

		% interval of time steps to write values at time step to output
		ks = 1;
	end % properties
	methods
		% default constructor
		function obj =	River_Tide_Morphodynamics(varargin)
			obj = obj@River_Tide(varargin{:});
		end % constructor

	end % methods
end % class River_Tide_Morphodynamics
