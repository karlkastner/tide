% Sa 20. Feb 10:28:32 CET 2016
% Karl Kastner, Berlin
%
%% surface levelation envelope of the tide
%% low water, high water and tidal range for lunar each day
%% 
%  function [timei, lmini, lmaxi, rangei, midrangei, phii] = tidal_envelope2(time,L,order)
%%
%% input:
%%        time  : 
%%        L     : surface elevation
%%        order : interpolation order (default 2)
%% ouput:
%%        timei : vector eqispaced
%%        lmini : minimum level
%%        lmaxi : maximum level
%%        rangei : range
%%        midrangei : (min + max)/2, usually different from mean
%%        phii  : pseudo phase
%%
%% Note: the pseudo phase phi jumps, this is because if the tide is semidiurnal,
%%       sometimes the lower hw becomes the next day higher then than the 
%%	current high water, e.g. there is no smooth transition by
%%       51min but a jump by 12h
% TODO resample to a multiple of 25/24
function [timei, lmini, lmaxi, rangei, midrangei, phii] = tidal_envelope2(time,L,order)
	if (nargin < 3)
		order = 2;
	end
	if (isvector(L))
		L = cvec(L);
	end

	% length of tidal day
	lunarday = 25/24;

	% maximum gap to interpolate over
	% TODO no magic numbers
	% must be 1.1, as maxima are at most 1.1 days apart
	dt_max   = 1.1;

	% resample to make steps uniform
	dt     = median(diff(time));
	t0     = time(1);
	timei  = (t0:dt:time(end))';
	n      = length(timei);
	n2     = round(lunarday/dt);
	n1     = floor(n/n2);
	n      = n1*n2;
	t2     = reshape(timei(1:n),n2,n1);

	nti   = length(timei);
	lmaxi = NaN(nti,size(L,2));
	lmini = NaN(nti,size(L,2));
	phii  = NaN(nti,size(L,2));

	% for each colum (time series)
	for ldx=1:size(L,2)
		level = L(:,ldx);
	
		% TODO do not resample here
		li     = interp1(time,level,timei);
	
		[lmax tmax lmin tmin] = interval_extrema(t0,dt,li,n2,n1);
	
		% pseudo phase
		%sunday = 1;
		%phi = 2*pi*(rvec(tmax)-floor(rvec(tmax)))./sunday;
		phi = 2*pi*(rvec(tmax)-floor(rvec(tmax)))./lunarday;
		%phi = 2*pi*(rvec(tmax)-t2(1,:))./lunarday;
		%im = 'nearest';
		im = 'linear';
		fdx = isfinite(tmax);
		if (sum(fdx) > 1)
			phii(:,ldx) = interp_angle(tmax(fdx),phi(fdx),timei,im);
		
			% make max and min continuous
			% limit interpolation to dt_max
			im='spline';
			lmaxi(:,ldx) = interp1_limited(tmax(fdx),lmax(fdx),timei,dt_max,im);
			fdx          = isfinite(tmin);
			lmini(:,ldx) = interp1_limited(tmin(fdx),lmin(fdx),timei,dt_max,im);
		end
	end % for ldx

	% range
	if (nargout() > 2)
		rangei     = lmaxi-lmini;
		midrangei  = 1/2*(lmaxi+lmini);
		% interpolation can make values negative in rare cases, avoid this
		% TODO, this is a quick fix, could be avoidbed by using pchip instead of spline
		fdx         = isfinite(rangei);
		rangei(fdx) = max(0,rangei(fdx));
	end
end % tidal_envelope2


