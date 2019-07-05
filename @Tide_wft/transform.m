% Sat  9 Jul 12:55:49 CEST 2016
% Karl Kastner, Berlin
%
%% wavelet transform tidal time series
%% input:
%% time   : [1xn] abszissa of input vector, for example time, must be equally spaced
%% val    : [1xn] signal, input data series (e.g water level or velocity)
%% F      : [1xm] base frequencies, 1, 1, 2, ... for mean level, diurnal, semidirunal ...
%% 		base periods from base frequencies T=1/F
%% n      : [1xm] wavelet window length in multiple of periods
%% fc, nc : [scalar] low frequency cutoff and window length in periods
%% winstr : [char] fourier windows (kaiser (recommended), hanning, box, etc)
%% dt_max : [scalar] maximum time to fill gaps in input data series (recommended 3/24 for tide)
%% output:
%% tide   : struct with fields
%%          w_coeff   : [1xn] wavelet coefficients (complex)
%%          amplitude : amplitude
%%          phase     : phase
%%          range     :
%%          h_tide    :
%%          h_low     :
%%          h
function obj = transform(obj,t0,dt,val)
	val = cvec(val);

	obj.dt   = dt;
	obj.t0   = t0;
	obj.nt   = length(val);

	time     = obj.time();

	% fill short gaps in input data series
	val      = fixnan(time,val,obj.dt_max,obj.imethod,NaN);

	% tidal species by wavelet transform
	obj.w_coeff    = wavelet_transform(val,dt,obj.F,obj.n,obj.winstr,obj.pmin);
	obj.amplitude  = abs(obj.w_coeff);
	obj.phase      = angle(obj.w_coeff);

	% reconstruct the signal
	obj.h_tide     = wavelet_reconstruct( obj.w_coeff,dt,obj.F,obj.n,obj.winstr );

	% floating mean depth (low frequency modulation by river flow)
	n_       = round(obj.n_low/dt);
	f_sample = 1./dt;
	fc       = obj.F_low;

	% TODO padd ends
	%h_low           = lowpass2(val,n_,f_sample,f_cutoff,nanflag);
	h_low           = lowpass2(val-obj.h_tide,n_,f_sample,fc,obj.nanflag);
	obj.h_low       = h_low;

	% tidal range
	% TODO timing of hhw with new function
	% TODO inflection points with new functions

	[timei, lmin, lmax, range, midrange] = tidal_envelope2(time,val);
	obj.range      = interp1_limited(timei,range,time,obj.dt_max,obj.imethod,NaN);
	obj.midrange   = interp1_limited(timei,midrange,time,obj.dt_max,obj.imethod,NaN);

	obj.h     = obj.h_tide + obj.h_low;

	obj.h_res = obj.h-val;

end % transform

