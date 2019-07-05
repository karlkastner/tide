% Wed 20 Jul 13:46:11 CEST 2016
% Karl Kastner, Berlin
%
%% slack water envelope of the tide
%
% t0 : start of series
% dt : time step
% T0 : period of principle component
% Q  : either discharge or velocity, axis pointing in downstream direction

% TODO interpolate extact zero or even fit a series (roots_fourier)
function [hws, lws] = envelope_slack_water(T0,t0,dt,Q)
	% number of samples per day
	m  = round(25/24*T0/dt);
	% number of days
	n  = floor(length(Q)/m);
	nm = n*m;
	time = t0 + dt*(0:nm);
	
	% slice the time series into n periods of length T0 (m)
	Td  = reshape(time(1:nm),m,n);
	Qd  = reshape(Q(1:nm),m,n);

	% padd last value
	Qd = [NaN, Qd(end,1:end-1)
               Qd ];
	Td = [Td(1)-dt, Td(end,1:end-1)
              Td ];
	% times in-between slots
	Td  = 0.5*(Td(1:end-1,:)+Td(2:end,:));
	
	% flat high water slack as 1
	% transition from discharging to reverse flow
	fhws      = (Qd(1:end-1,:)<=0) & (Qd(2:end,:) >=0);
	[void id] = max(fhws);
	id        = sub2ind([m n],id,1:n);
	hws.t       = Td(id);
	% invalid days, where there is no slack
	hws.t(0 == fhws(id)) = NaN;
	
	% flag low water slack as 1
	% transition from reverse flow to discharging)
	flws      = (Qd(1:end-1,:)>=0) & (Qd(2:end,:) <=0);
	[void id] = max(flws);
	id        = sub2ind([m n],id,1:n);
	lws.t       = Td(id);
	lws.t(0 == flws(id)) = NaN;
end

