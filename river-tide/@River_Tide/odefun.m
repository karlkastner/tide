% Sun  8 Oct 13:08:39 CEST 2017
%
%% coefficients of the backwater and wave equation for river-tides
%
% TODO make q0 an unknown for networks (stacked odes with multiple bcs)
% TODO precompute cD, h, w, dhdx, dwdx,
%      this requires the bvp solve to accept a predefined mesh as an argument
% TODO account for inhomogeneous part
function [f, obj] = odefun(obj,x,y)
	%switch (obj.opt.hmode)
	%case {'matrix'}
	%	k = sum(obj.opt.oflag)+1;
	%otherwise
	%	k = sum(obj.opt.oflag);
	%end
	k = obj.neq;

	g      = obj.g;
	omega1 = obj.omega;
	flag   = obj.flag;

	w      = obj.width(x);
	D1_dx  = obj.D1_dx(x);
	dw_dx  = D1_dx*w;

	Q0     = obj.Q0_*ones(size(x));
	zb     = obj.zb(x);
	nx     = length(x);

	if (isscalar(y))
	        % TODO quick fix	
		y_ = repmat(y,1,k);
	else
		y_ = reshape(y,nx,k);
	end

        switch (obj.opt.hmode)                                                  
	case {'matrix'}
		z0 = y_(:,1); % y(1:nx);
		k_ = 2;
	otherwise
		% z0 is computed below
		%z0 = obj.z0(x);
		k_ = 1;
	end
	Qt = y_(:,k_:end);

	% TODO properly determine range and midrange
	%Qhr    = sum(abs(Qt),2);
	%Qmid   = Q0;
	%[Qhr,Qmid] = tidal_range_exp([Q0,Qt]);
	Qhr = sum(abs(Qt),2);
	Qmid = Q0;

	cf = obj.friction_coefficient(Qmid,Qhr);

	% TODO rename, iterate is idiosynchratic
	switch (obj.opt.hmode)
	case {'no-tide','predefined'}
		if (~isempty(obj.z0))
			z0 = obj.tmp.z0;
		else
			z0 = obj.solve_backwater(x,Q0,@(x) 0);
		end
		% neglect tidal influence on the tidally averaged water level
		k = 1;
	case {'iterate'}
		% update tide
		z0 = obj.solve_backwater(x,Q0,Qt);
		k           = 1;
	case {'matrix'}
		% depth is reiterated together with Qt
		z1     = obj.discharge2level(Qt(:,1),x);
		h0     = z0 - zb;
		h0     = max(h0,obj.hmin);
		cD     = obj.cd(x,h0);
		% TODO make cd a function
		f  = obj.odefun0(x, h0, z1, zb, w, dw_dx, Q0, Qhr, Qt, cD, cf);
		k  = 2;
	otherwise
		error('odefun');
	end

	h0     = z0 - zb;
	h0     = max(h0,obj.hmin);
	if (min(h0)<=0)
		warning('negative water depth')
	end

	h0     = max(h0,obj.hmin);
	cD     = obj.cd(x,h0);
	dh0_dx = D1_dx*h0;
	dz0_dx = D1_dx*z0;

	Q  = [Q0,Qt];
	QQ = fourier_quadratic_interaction_coefficients(Q,size(Q,2),2);

	% frequency components of the tide
	for idx=1:length(obj.opt.oflag)
	    if (obj.opt.oflag(idx))
		f(:,:,k) = obj.odefunk(idx, Q, QQ, Qhr, h0, dh0_dx, dz0_dx, ...
						      w, dw_dx, cD, cf, D1_dx);
	    	k=k+1;
	    end % if
	end % for 
%f(:,:,1)
%f(:,:,2)
%pause
end % River_Tide/odefun

