% 2017-05-02 13:48:53.721763802 +0200
%% mimimum and maximum within intervals of constant length,
%% intended for periodic functions
% TODO merge with tidal envelope2
function [lmax, tmax, lmin, tmin] = interval_extrema2(t0, dt, L, T)
	if (isvector(L))
		L = cvec(L);
	end

	% length of tidal day
	%lunarday = 25/24;

	% maximum gap to interpolate over
	% TODO no magic numbers
	% must be 1.1, as maxima are at most 1.1 days apart
	% dt_max   = 1.1;

	% resample to make steps uniform
	%dt     = median(diff(time));
	%timei  = (time(1):dt:time(end))';

	n      = size(L,1);
	%n      = length(timei);
	n2     = round(T/dt);
	n1     = floor(n/n2);
	% n      = n1*n2;
	%t2     = reshape(timei(1:n),n2,n1);

	% for each colum (time series)
	for ldx=1:size(L,2)
		level = L(:,ldx);

		% TODO do not resample here
		%li     = interp1(time,level,timei);
		li = level;

		%[lmax(:,ldx) tmax(:,ldx) lmin(:,ldx) tmin(:,ldx)] = interval_extrema(timei(2)-timei(1),li,n2,n1);
		[lmax(:,ldx), tmax(:,ldx), lmin(:,ldx), tmin(:,ldx)] = interval_extrema(t0,dt,li,n2,n1);
	end
	% shift
%	tmax = tmax + t0;
%	tmin = tmin + t0;
end

