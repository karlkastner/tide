% Wed 11 Oct 10:18:54 CEST 2017
%% solve for the oscillatory (tidal) componets
%%
%% function obj   = solve_wave(obj)
function obj   = solve_wave(obj)

	% solve the system of ode's
	[x, y, cflag] = feval(obj.opt.solver,@obj.odefun,@obj.bcfun,obj.Xi,obj.opt);
	
	% extract unknowns
	nx = length(x);
	obj.x  = x;
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
	end % switch
	obj.Q_ = [obj.fun.Q0(x), Q1, Q2];

	obj.z_ = [z0, obj.discharge2level([Q1, Q2])];
	
	% TODO awkward
	obj.zb = obj.fun.zb(x);

	% TODO this is awkward to be stored here
	obj.initial.z0 = obj.tmp.z0(x);

	% check result
	if (any(obj.tmp.h0(x)<=0))
		cflag = -2;
		warning('negative water depth');
	end
	if (any(~isfinite(obj.z_)))
		cflag = -1;
		warning('solution is not finite');
	end

	obj.cflag = cflag;
end % River_Tide/solve_wave

