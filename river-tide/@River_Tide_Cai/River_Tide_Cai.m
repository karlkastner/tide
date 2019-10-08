% Tue  4 Apr 09:21:39 CEST 2017
%% prediction of river tide by the method of Cai (2014)
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

