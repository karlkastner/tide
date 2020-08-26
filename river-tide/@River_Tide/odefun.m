% Sun  8 Oct 13:08:39 CEST 2017
%
%% coefficients of the backwater and wave equation for river-tides
%
% TODO precompute Cd, h, w, dhdx, dwdx,
%      this requires the bvp solve to accept a predefined mesh as an argument
% TODO heed identity abs(Qi)^2 = abs(Qi^2) = conj(Qi)*Qi
%      to move coupling terms from rhs (b) to lhs (A)
function [f, obj] = odefun(obj,x,y)
	% without any input, order of coupled odes is returned
	if (nargin()<2 || isempty(x))
		f(1) = 1;
		%if (obj.opt.dischargeisvariable)
		%	f(2) = -1;
		%	k = 2;
		%else
			k = 1;
		%end
		for idx=1:length(obj.opt.oflag)
			if (obj.opt.oflag(idx))
				k    = k+1;
				f(k) = 2;
			end
		end
		return;
	end

	k = obj.neq;

	g      = obj.g;
	omega1 = obj.omega;
	flag   = obj.flag;

	w      = obj.width(x);
	D1_dx  = obj.D1_dx(x);
	dw_dx  = D1_dx*w;

	zb     = obj.zb(x);
	nx     = length(x);

	[z0,Q0,Qt] = obj.extract(x,y);
	Q0         = repmat(Q0,nx,1);

	% TODO properly determine range and midrange
	%Qhr    = sum(abs(Qt),2);
	%Qmid   = Q0;
	%[Qhr,Qmid] = tidal_range_exp([Q0,Qt]);
	Qhr = sum(abs(Qt),2);
	Qmid = Q0;

	cf = -obj.friction_coefficient(Qmid,Qhr);

%	if (obj.opt.dischargeisvariable)
%		f = zeros(nx,4,obj.neq+1);
%	else
		f = zeros(nx,4,obj.neq);
%	end

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
		if (~isempty(Qt))
			z1     = obj.discharge2level(Qt(:,1),x);
		else
			z1 = 0;
		end
		h0     = z0 - zb;
		h0     = max(h0,obj.hmin);
		Cd     = obj.cd(x,h0);
		% TODO make cd a function
		if (~obj.opt.dischargeisvariable)
		f(:,:,1)  = obj.odefunz0(x, h0, z1, zb, w, dw_dx, Q0, Qhr, Qt, Cd, cf);
		else
		f(:,:,1)  = obj.odefunQ0(x, h0, z1, zb, w, dw_dx, Q0, Qhr, Qt, Cd, cf);
		end
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
	Cd     = obj.cd(x,h0);
	dh0_dx = D1_dx*h0;
	dz0_dx = D1_dx*z0;

	Q  = [Q0,Qt];
	QQ = fourier_quadratic_interaction_coefficients(Q,size(Q,2),2);

%	if (obj.opt.dischargeisvariable)
%		f(:,:,k) = obj.odefunQ0(x,h0,z1,zb,w,dw_dx,Q0,Qhr,Qt,Cd,cf);
%		k = 3;
%	end

	% frequency components of the tide
	for idx=1:length(obj.opt.oflag)
	    if (obj.opt.oflag(idx))
		f(:,:,k) = obj.odefunk(idx, Q, QQ, Qhr, h0, dh0_dx, dz0_dx, ...
						      w, dw_dx, Cd, cf, D1_dx);
	    	k=k+1;
	    end % if
	end % for 
%f(:,:,1)
%f(:,:,2)
%pause
end % River_Tide/odefun

