% Thu  7 Jul 16:34:56 CEST 2016
% Karl Kastner, Berlin
% c.f. Jay & Jukulka 2003a,b
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

