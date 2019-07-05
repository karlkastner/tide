% Sat Jul 13 11:16:40 UTC 2013
% Karl Kästner, Berlin
%
%% extract tidal envelope from time series
%
% TODO move tidal_envelope into Tide_Table
% TODO plausibility check for double neap, double spring, double lw/hw
% TODO moonrise, sunrise

% generate a tide tab from a water level time series
function tidetable = analyze(tidetable)

	t     = tidetable.time;
	level = tidetable.level;

	dt = tidetable.dt;
	
	% relative water level
	printf('mean water level %f\n',	mean(tidetable.level));
	level = level - mean(tidetable.level);

	% quick and dirty equal height fix
	% TODO should be made symmetric and even betidetableer incorporated into the envelope function
	fdx = find(diff(level) == 0);
	level(fdx) = level(fdx).*(1+1e-7*(1:length(fdx))');
	fdx = find(diff(level) == 0);

	% find the tidal envelope curve as well as
	% the daily low and high water level as well as forthnightly neap spring maxima
	[tl, vl, ldx, th, vh, hdx, neap_dx, spring_dx, dmin, dmax, drange] = tidal_envelope(86400*tidetable.time,level);
	tl = tl/86400;
	th = th/86400;

	nx = int32(1/dt);
	ny = int32(length(t)/nx);
	t24 = mean(reshape(t,nx,ny)',2);
	D   = mean(reshape(drange,nx,ny)',2);

%	time    = t/(86400); 
%	t24     = t24/(86400);

	tidetable.t24       = t24;
	tidetable.range24   = D;
	tidetable.drange    = drange;
	tidetable.dmin      = dmax;
	tidetable.dmax      = dmin;
	tidetable.neap_dx   = neap_dx;
	tidetable.spring_dx = spring_dx;
	tidetable.tneap     = t(neap_dx);
	tidetable.tspring   = t(spring_dx);
	tidetable.tneap     = t(neap_dx);
	tidetable.tspring   = t(spring_dx);
	tidetable.tl	     = tl;
	tidetable.vl        = vl;
	tidetable.th        = th;
	tidetable.vh        = vh;
	tidetable.hdx       = hdx;
	tidetable.ldx       = ldx;

end % function Tidetable.analyze()


