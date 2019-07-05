% Thu 21 Jul 12:26:11 CEST 2016
% Karl Kastner, Berlin
%% compute envelopes of hw and low water
function [hw, lw, range] = envelope_amplitude(T0,t0,dt,H)
	% number of samples per day
	m  = round(25/24*T0/dt);
	% number of days
	n  = floor(length(H)/m);
	nm = n*m;
	time = t0 + dt*(0:nm);
	
	% slice the time series into n periods of length T0 (m)
	Td  = reshape(time(1:nm),m,n);
	Hd  = reshape(H(1:nm),m,n);

	% high water
	[hw.l id] = max(Hd);
	id        = sub2ind([m n],id,1:n);
	hw.t      = Td(id);
	
	% low water
	[lw.l id] = min(Hd);
	id        = sub2ind([m n],id,1:n);
	lw.t      = Td(id);

	% range
	range.t   = mean(Td);
	range.l   = hw.l - lw.l;
end

