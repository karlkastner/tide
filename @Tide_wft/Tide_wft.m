% Wed 27 Sep 10:34:52 CEST 2017
%% wavelet transform of tidal time series
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

