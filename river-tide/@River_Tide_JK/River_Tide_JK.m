% Thu  7 Jul 16:34:56 CEST 2016
% Karl Kastner, Berlin
% c.f. Jay & Jukulka 2003a,b
classdef < River_Tide_JK < handle
	properties
		% coefficients fit to the tide
		d
		% drag coefficient (constant)
		cD
		% period of tidal cycle (constant)
		T
		% 
		g = Constant.gravity;
	end % properties
	methods
		function obj = River_Tide_JK()
		end % constructor
		% pseudo members
		function omega_ = omega(obj)
			omega_ = 2*pi./obj.T;
		end % omega
	end % methods
end % classdef River_Tide_JK

