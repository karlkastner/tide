% Wed 11 Oct 10:18:54 CEST 2017
%% solve for the oscillatory (tidal) componets
%%
%% function obj   = solve_wave(obj)
function obj   = solve_wave(obj)
	% solve the system of ode's
	[x, y, cflag] = feval(obj.opt.solver,@obj.odefun,@obj.bcfun,obj.Xi,obj.opt);
	
	% extract unknowns
	nx = length(x);
	switch (obj.opt.hmode)
	case {'matrix'}
		z0 = y(1:nx);
		Q1 = y(nx+1:2*nx);
		if (obj.opt.o2)
			Q2 = y(2*nx+1:3*nx);
		else
			Q2 = [];
		end
	otherwise
		% TODO compute h0 by function as zs0 and zb are stored
		z0 = obj.tmp.z0(x);
		Q1 = y(1:nx);
		if (obj.opt.o2)
			Q2 = y(nx+1:2*nx);
		else
			Q2 = [];
		end
	end

	obj.x  = x;
	obj.w  = obj.wfun(x);
	obj.Q_ = [obj.Q0fun(x).*Q1.^0,Q1,Q2];

	% TODO sloppy, make this better a non-static function
	% TODO pass width
	z1     = obj.q2z(x,Q1./obj.w,obj.omega);
	if (~isempty(Q2))
		z2     = obj.q2z(x,Q2./obj.w,2*obj.omega);
	else
		z2     = [];
	end
	obj.z_  = [z0,z1,z2];
	
	% TODO awkward
	obj.zb = obj.fun.zb(x);

	% TODO this is awkward to be stored here
	obj.initial.z0 = obj.tmp.z0(x);

	obj.cflag = cflag;
end % River_Tide/solve_wave

