% Sun  8 Oct 13:08:39 CEST 2017
%
%% coefficients of the backwater and wave equation for river-tides
%
% TODO make q0 an unknown for networks (stacked odes with multiple bcs)
% TODO precompute cD, h, w, dhdx, dwdx,
%      this requires the bvp solve to accept a predefined mesh as an argument
% TODO account for inhomogeneous part
function [f, obj] = odefun(obj,x,y)
	% TODO dirty hack to determin the number of equations
	if (nargin()<2)
		switch (obj.opt.hmode)
		case {'matrix'}
			k=2;	
		otherwise
			k=1;
		end
		if (obj.opt.o2)
			k=k+1;
		end
			f = zeros(0,0,k);
		return;
	end

	% TODO, move to init
	if (isempty(obj.tmp.D1))
		obj.tmp.D1 = derivative_matrix_1_1d(x,[],2);
	end

	g      = obj.g;
	omega1 = obj.omega;
	flag   = obj.flag;

	w      = obj.fun.width(x);
	%dw_dx  = obj.tmp.D1*w;
	dw_dx  = derivative1(x,w); %obj.tmp.D1*w;

	Q0     = obj.fun.Q0(x);
	zb     = obj.fun.zb(x);
	nx = length(x);
        switch (obj.opt.hmode)                                                  
	case {'matrix'}
		z0 = y(1:nx); 
		Q1 = y(nx+1:2*nx);
		if (obj.opt.o2)
			Q2 = y(2*nx+1:3*nx);
		else
			Q2 = 0;
		end
	otherwise
		z0 = obj.tmp.z0(x);
		Q1  = y(1:nx);
		if (obj.opt.o2)
			Q2 = y(nx+1:2*nx);
		else
			Q2 = 0;
		end
	end
	% TODO consider Q2 to determine Qhr and pass it to ODEfun1 and ODEfun2
	% TODO consider influence of Q2 on z0
	Qhr    = abs(Q1);
	Qmid   = Q0;

	switch (obj.opt.friction_method)
	case {'neglect-river'} % lorentz
		% identical to Dronkers for Q0=0
		cf = obj.friction_coefficient_lorentz(abs(Qmid)./Qhr);
	case {'dronkers'}
		% half range of tidal amplitude
		cf   = -obj.friction_coefficient_dronkers(cvec(abs(Qmid)./Qhr));
	case {'godin'}
		cf   =  obj.friction_coefficient_godin(abs(Qmid)./Qhr);
	case {'savenije'}
		error('TODO');
	case {'no-reverse'}
		% identical to Dronkers for Q0>Qt
		c = [0,0,-pi];
	otherwise 
		error('objfun');
	end

	switch (obj.opt.hmode)
	case {'no-tide','predefined'}
		% neglect tidal influence on the tidally averaged water level
		k = 1;
	case {'iterate'}
		% recompute the backwater curve
	if (nargin(obj.fun.cd) < 2)
	C = @(x_) drag2chezy(obj.fun.cd(x_));
	else
	C = @(x_,h_) drag2chezy(obj.fun.cd(x_,h_));
	end
		Q1fun = @(x_)  interp1(x,Q1,x_,obj.opt.imethod,'extrap');
		obj.backwater.sopt.InitialStep = x(2)-x(1);
		%obj.backwater.sopt.RelTol = 10*eps^0.25;
		[x_, h0_, z0_] = obj.backwater.solve(Q0(1), ...
					Q1fun, ...
					C,...
					obj.fun.width,obj.fun.zb, ...
					obj.z0_downstream(1), ...
					obj.Xi);
		obj.tmp.x   = x_;
		obj.tmp.h0  = @(x) interp1(x_,h0_,x,obj.opt.imethod);
		obj.tmp.z0  = @(x) interp1(x_,z0_,x,obj.opt.imethod);
		z0          = obj.tmp.z0(x);
		k           = 1;
	case {'matrix'}
		% depth is reiterated together with Q1
		z1 = obj.q2z(x,Q1./w,omega1);
		f  = obj.odefun0(x, z0, z1, zb, w, dw_dx, Q0, Qhr, Q1, cD, cf);
		k  = 2;
	otherwise
		error('odefun');
	end
	h0     = z0 - zb;
	cD     = obj.fun.cd(x,h0);
	if (min(h0)<=0)
		error('negative water depth')
	end
	% TODO h and z were mixed up here
	%dh0_dx = obj.tmp.D1*h0;
	%dz0_dx = obj.tmp.D1*z0;
	% TODO store interpolators of dh and dz
	dh0_dx = derivative1(x,h0);
	dz0_dx = derivative1(x,z0);

	[f(:,:,k)] = obj.odefun1(Q0, Qhr, Q1, Q2, h0, dh0_dx, dz0_dx, w, dw_dx, cD, g, cf, omega1, flag);
	if (obj.opt.o2 == true)
		f(:,:,k+1) = obj.odefun2(Q0, Qhr, Q1, Q2, h0, dh0_dx, dz0_dx, w, dw_dx, cD, g, cf, omega1, flag);
	end
end % River_Tide/odefun

