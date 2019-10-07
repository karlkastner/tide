% Tue  4 Apr 09:21:39 CEST 2017
%% prediction of river tide by the method of Cai (2014)
classdef River_Tide_Cai < handle
	properties
		mode = 'lorentz';

		abstol = 1e-3;
		reltol = 1e-3;
		maxiter = 100;
		% relaxation factor
		relax = 0.1;
	end
	methods (Static)
	
		%Gamma = Gamma(k,zeta,lambda,mu,phi);
		%[delta lambda mu e] = rt_quantities(zeta,gamma,chi,rs,phi);
	end % methods (Static)
	methods
		% constructor
		function obj = River_Tide_Cai(mode)
			if (nargin()>0)
				obj.mode = mode;
			end
		end
	end % methods 
end% class River_Tide_Cai

