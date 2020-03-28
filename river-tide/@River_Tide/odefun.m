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
			k = sum(obj.opt.oflag)+1;
		otherwise
			k = sum(obj.opt.oflag);
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
	dw_dx  = derivative1(x,w);

	Q0     = obj.fun.Q0(x);
	zb     = obj.fun.zb(x);
	nx = length(x);
	y_ = reshape(y,nx,[]);
        switch (obj.opt.hmode)                                                  
	case {'matrix'}
		z0 = y_(:,1); % y(1:nx);
		k_ = 2;
	otherwise
		z0 = obj.tmp.z0(x);
		k_ = 1;
	end
	Qt = y_(:,k_:end);
	%Q1 = y(nx+1:2*nx);
	%if (obj.opt.o2)
	%		Q2 = y(2*nx+1:3*nx);
	%		Qt = [Q1,Q2];
	%	else
	%		Q2 = 0;
	%		Qt = Q1;
	%	end
	%otherwise
	%	Q1  = y(1:nx);
	%	if (obj.opt.o2)
	%		Q2 = y(nx+1:2*nx);
	%		Qt = [Q1,Q2];
	%	else
	%		Q2 = 0;
	%		Qt = Q1;
	%	end
	%end
	% TODO properly determine range and midrange
	%Qhr    = sum(abs(Qt),2);
	%Qmid   = Q0;
	%[Qhr,Qmid] = tidal_range_exp([Q0,Qt]);
	Qhr = sum(abs(Qt),2);
	Qmid = Q0;


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
		Qtfun = @(x_)  interp1(x,Qt,x_,obj.opt.imethod,'extrap');
		obj.backwater.sopt.InitialStep = x(2)-x(1);
		%obj.backwater.sopt.RelTol = 10*eps^0.25;

		if (obj.bc(1,1).var == 'z')
			z0_downstream = obj.bc(1,1).rhs;
			[x_, h0_, z0_] = obj.backwater.solve( ...
					Q0(1), ...
					Qtfun, ...
					C,...
					obj.fun.width, ...
					obj.fun.zb, ...
					z0_downstream, ...
					obj.Xi);
						%obj.Q0_, Qt, Cfun, ...
						%obj.fun.width, obj.fun.zb, ...
						%z0_downstream, obj.Xi);
		else
			z0_downstream = obj.bc(2,1).rhs;
			[x_, h0_, z0_] = obj.backwater.solve( ...
					Q0(1), ...
					Qtfun, ...
					C,...
					obj.fun.width, ...
					obj.fun.zb, ...
					z0_downstream, ...
					[obj.Xi(2),obj.Xi(1)]);
			x_ = flipud(x_);
			h0_ = flipud(h0_);
			z0_ = flipud(z0_);
		end

		obj.tmp.x   = x_;
		obj.tmp.h0  = @(x) interp1(x_,h0_,x,obj.opt.imethod);
		obj.tmp.z0  = @(x) interp1(x_,z0_,x,obj.opt.imethod);
		z0          = obj.tmp.z0(x);
		k           = 1;
	case {'matrix'}
		% depth is reiterated together with Qt
		z1 = obj.discharge2level(Qt(:,1),x);
		h0     = z0 - zb;
		cD     = obj.fun.cd(x,h0);
		% TODO make cd a function
		f  = obj.odefun0(x, z0, z1, zb, w, dw_dx, Q0, Qhr, Qt, cD, cf);
		k  = 2;
	otherwise
		error('odefun');
	end

	h0     = z0 - zb;
	cD     = obj.fun.cd(x,h0);
	if (min(h0)<=0)
		error('negative water depth')
	end
	% TODO store interpolators of dh and dz
	dh0_dx = derivative1(x,h0);
	dz0_dx = derivative1(x,z0);

	Q  = [Q0,Qt];
	QQ = fourier_quadratic_interaction_coefficients(Q,size(Q,2),2);

	% frequency components of the tide
	for idx=1:length(obj.opt.oflag)
	if (obj.opt.oflag(idx))
		f(:,:,k) = obj.odefunk(idx, Q, QQ, Qhr, h0, dh0_dx, dz0_dx, w, dw_dx, cD, cf);
		k=k+1;
	end
	end

end % River_Tide/odefun

