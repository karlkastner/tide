% Sa 20. Feb 10:28:32 CET 2016
% Karl Kastner, Berln
%% times and evelations for high and low water
%
%function [lmax tmax lmin tmin] = interval_extrema(dt,l,n1,n2,order)
function [lmax tmax lmin tmin] = interval_extrema(t0,dt,l,n1,n2,order)
	if (~isscalar(dt))
		error('dt not scalar');
	end
	if (nargin()<6)
		order = 2;
	end
	time = t0 + dt*(0:n1*n2-1);

	% true extrema are preselected here,
	% as lmiting values for a day can occur at the boundary between days
	% and may not coincide with true extrema (f'' = 0)
	ismax = [false; (l(2:end-1) > l(1:end-2)) & (l(2:end-1)>=l(3:end)); false];
	ismax = reshape(ismax(1:n2*n1),n1,n2);
	ismin = [false; (l(2:end-1) < l(1:end-2)) & (l(2:end-1)<=l(3:end)); false];
	ismin = reshape(ismin(1:n2*n1),n1,n2);

	l2     = reshape(l(1:n1*n2),n1,n2);

	% maxima/minima must not be exactly zero, set zeros to machine precission
	l2(0 == l2) = eps(class(l2));

	% minimum and maximum water level in each interval
	l2_ = l2;
	l2_(~ismax) = -inf;
	[lmax imax] = max(l2_);
	l2_ = l2;
	l2_(~ismin) = inf;
	[lmin imin] = min(l2_);

	% translate indices from local intra-day to global entire time series
	imax        = (0:n2-1)*n1 + imax;
	imin        = (0:n2-1)*n1 + imin;

	% mark periods without extrema
	fmax = (ismax(imax) == 0);
	fmin = (ismin(imin) == 0);
	
	% higher order accurate determination of extrema by quadratic polynomials
	switch (order)
	case {0,1}
		tmax        = time(imax);
		tmin        = time(imin);
	otherwise
if (1)
		% use quadratic approximation to determine time and height/depth of hhw
		[lmax tmax] = extreme3(time,l,imax);
		[lmin tmin] = extreme3(time,l,imin);
else
		imax = max(2,imax);
		imax = min(n2*n1-1,imax);
		[tmax lmax] = extreme_quadratic([time(imax-1) time(imax) time(imax+1)],[cvec(l(imax-1)) cvec(l(imax))  cvec(l(imax+1))]);
		tmax = tmax';
		lmax   = lmax';

		imin = max(2,imin);
		imin = min(n2*n1-1,imin);
		[tmin lmin] = extreme_quadratic([time(imin-1) time(imin) time(imin+1)],[cvec(l(imin-1)) cvec(l(imin))  cvec(l(imin+1))]);
		tmin = tmin';
		lmin   = lmin';
end
	end
	% invaldate false extrema
	tmax(fmax) = NaN;
	lmax(fmax) = NaN;
	tmin(fmin) = NaN;
	lmin(fmin) = NaN;

	tmax = tmax';
	lmax = lmax';
	tmin = tmin';
	lmin = lmin';
end % interval_extrema

