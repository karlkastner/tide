% Sat Jul 13 11:16:40 UTC 2013
% Karl Kästner, Berlin
%
%% envelope of the tide
%%
%% input : t time in days
%%         f surface elevation
%% ouput : tl time of low water
%%         vl surface elevation at low water
%%         ldx index of low water
%%         th time of high water
%%         vh surface elevation at high water
%%         hdx index of high water
%%         ndx neap index
%%         sdx spring index
%%         dmax:
%%         drange: range per day
function [tl, vl, ldx, th, vh, hdx, ndx, sdx, dmax, dmin, drange] = tidal_envelope(t, f)
	% diurnal envelope
	n = round(86400*1/(t(2)-t(1)));
	fmin = ordfilt2(f, 1, ones(n,1));
	fmax = ordfilt2(f, n, ones(n,1));

	% low water
	ldx     = find(f(2:end-1) < f(1:end-2) & f(2:end-1) < f(3:end)) + 1;
	[tl vl] = inter(t,f,ldx);
	% high water
	hdx     = find(f(2:end-1) > f(1:end-2) & f(2:end-1) > f(3:end)) + 1;
	[th vh] = inter(t,f,hdx);

	% envelope for spring-neap cycle
	% ignores intradiurnal minor high and low tides
	fdx    = find(f(2:end-1) == fmin(2:end-1) & f(2:end-1) < f(1:end-2) & f(2:end-1) < f(3:end)) + 1;
	dmin   = interp1(fdx,f(fdx), (1:length(f))','spline');

	fdx    = find(f(2:end-1) == fmax(2:end-1) & f(2:end-1) > f(1:end-2) & f(2:end-1) > f(3:end)) + 1;
	dmax   = interp1(fdx,f(fdx), (1:length(f))','spline');

	% range
	drange = dmax - dmin;

	% find spring-neap
	n = round(5*86400*1/(t(2)-t(1)));
%	n = round(7*24*3600*1/(t(2)-t(1)));
	fmin = ordfilt2(drange, 1, ones(2*n,1));
	fmax = ordfilt2(drange, 2*n, ones(2*n,1));
	ndx  = find( drange(2:end-1) == fmin(2:end-1) & ...
		 drange(2:end-1) < drange(1:end-2) & drange(2:end-1) < drange(3:end)) + 1;
	sdx  = find( drange(2:end-1) == fmax(2:end-1) & ...
		 drange(2:end-1) > drange(1:end-2) & drange(2:end-1) > drange(3:end)) + 1;

	% interpolate the events from inbetween sampling intervals
	% syms f1 f2 f3 dt; y = [dt*dt -dt 1; 0 0 1; dt*dt dt 1] \ [f1; f2; f3]; t0=(-0.5*y(2)/y(1)), f0 = [t0^2 t0 1]*y
	function [ti fi] = inter(t,f,fdx)
		ti = zeros(size(fdx));
		fi = zeros(size(fdx));
		for idx=1:length(fdx)
			t2 = t(fdx(idx));
			dt = t(fdx(idx))-t(fdx(idx)-1);
			f1 = f(fdx(idx)-1);
			f2 = f(fdx(idx)  );
			f3 = f(fdx(idx)+1);
			ti(idx) = t2 + 0.5*dt*(f1 - f3)/(f1 - 2*f2 + f3);
			fi(idx) = f2 - 1/8*(f1 - f3)^2/(f1 - 2*f2 + f3);
			if (ti(idx) < t2 - dt || ti(idx) > t2 + dt)
				dt
				ti(idx)-t2
				f1
				f2
				f3
				error()
			end
		end % for idx
	end % function inter
end % envelope

