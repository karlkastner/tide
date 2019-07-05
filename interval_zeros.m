% Wed 29 Mar 16:13:27 CEST 2017
% Karl Kastner, Berlin

%% times of slack water determined frim velocity u
function [t_hws, t_lws, dt_hws, dt_lws] = interval_zeros(t0,dt,u,t_hw,t_lw,n1,n2,order)
	t       = t0+reshape(dt*(0:n1*n2-1),n1,n2);
%	t       = reshape(time(1:n1*n2),n1,n2);
%	dt = time(2)-time(1);

	% high water slack : change from flood to ebb flow
	fhw = [((u(1:end-1) >= 0) & (u(2:end) <= 0));false];
	fhw = fhw(1:n1*n2);

	% low water slack : change from flood to ebb flow
	flw = [((u(1:end-1) <= 0) & (u(2:end) >= 0));false];
	flw = flw(1:n1*n2);

	% time of mid-intervals
	tmid    = t + 1/2*dt;

	% time difference between slack and high water
	% 3 cases, slack water falls into day before, the same day or the day after
	mtmid         = tmid;
	mtmid(~(fhw)) = NaN;
	dt_hws  = bsxfun(@minus,  [NaN(n1,1),mtmid(:,1:end-1);    % day before
                                   mtmid;		          % same day
                                   mtmid(:,2:end),NaN(n1,1)], ... % day after
				   t_hw');
	mtmid         = tmid;
	mtmid(~(flw)) = NaN;
	dt_lws  = bsxfun(@minus,  [NaN(n1,1),mtmid(:,1:end-1);    % day before
                                   mtmid;		          % same day
                                   mtmid(:,2:end),NaN(n1,1)], ... % day after
				   t_lw');

	% nearest slack water
	[void id_hws] = min(abs(dt_hws)); % was max(.)
	fhw = isfinite(void);
	[void id_lws] = max(abs(dt_lws));
	flw = isfinite(void);

	% translate indices
	id_hws = id_hws-n1 + (0:n2-1)*n1;
	id_hws = max(2,id_hws);
	id_lws = id_lws-n1 + (0:n2-1)*n1;
	id_lws = max(2,id_lws);

	% time of slack water
	switch (order)
	case {0}
		% zeros from constant fit
		t_hws = tmid(id_hws);
		t_lws = tmid(id_lws);
	otherwise
		% zeros from linear fit
		t_hws  = zeros_p2([time(id_hws),time(id_hws+1)],[cvec(u(id_hws)),cvec(u(id_hws+1))]);
		t_lws  = zeros_p2([time(id_lws),time(id_lws+1)],[cvec(u(id_lws)),cvec(u(id_lws+1))]);
	end
	t_hws(~fhw) = NaN;
	t_lws(~flw) = NaN;

	dt_hws = t_hws - t_hw;
	dt_lws = t_lws - t_lw;

	%dt_hws  = bsxfun(@minus,tmid,t_hw);
	%dt_lws  = bsxfun(@minus,tmid,t_lw);

	% nearest slack water
	% make sign unique by wrapping
%	dt_hws(dt_hws<0) = dt_hws(dt_hws<0)+T;
%	dt_lws(dt_lws<0) = dt_lws(dt_lws<0)+T;

%	dt_hws(flw)  = 0;
%	dt_lws(flag) = 0;
	% mark days without high or low water slack respectively
%	fhws = (0 == dt_hws);
%	flws = (0 == dt_lws);
end % interval_zeros

