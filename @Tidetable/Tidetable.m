% Fri 15 Jul 17:35:03 CEST 2016
% Karl Kastner, Berlin
%
%% Tide table
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
classdef Tidetable < handle
	properties
	% primary TPXO values
		time
		level
		ux
		uy
		umag
		udir
		dt

	% analysis values
		dmax
		dmin
		drange
		hdx
		ldx
		neap_dx
		range24
		spring_dx
		t24
		th
		tl
		tneap
		tspring
		vh
		vl
		x
		y
		placename
	end
	methods
		function obj = Tidetable()
		end
	end
end % class Tidetable

