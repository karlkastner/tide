% Wed 27 Sep 10:34:52 CEST 2017
%% wavelet transform of tidal time series
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
classdef Tide_wft < handle
	properties
		t0
		dt
		nt
		h
		h_tide
		h_low
		h_res
		w_coeff
		amplitude
		phase
		range
		midrange

		location

		F
		n
		F_low
		n_low
		%fc
		%nc

		% wavelet window type
		winstr = 'kaiserwin';

		% default gap length
		dt_max = 3/24;

		% ???
		pmin = 1;

		% interpolation method
		imethod = 'spline';

		nanflag  = true;
	end % properties
	methods
		function obj = Tide_wft(varargin)
			for idx=1:2:length(varargin)-1
				obj.(varargin{idx}) = varargin{idx+1};
			end
		end
		function [time obj] = time(obj)
			time = obj.t0 + (0:obj.nt-1)'*obj.dt;
		end
	end % methods
end % class Tide_wft

