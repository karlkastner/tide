% Tue  6 Jun 11:23:16 CEST 2017
%% shift phase according to time zone
function obj = shift_time_zone(obj,dt)
	% days to hours
	dt = 24*dt;
	% time zone correction (phase is now w/r to local time)
	T  = 1./obj.freq;
	p  = obj.phase;
	p  = p + 360*dt./T;
	p  = rad2deg(wrapTo2Pi(deg2rad(p)));
	obj.t.tidecon(:,3) = p;
end

