% Fri  7 Aug 19:04:20 +08 2020
%
% function dzb_dt = dzb_dt(obj,t,zb)
%
function [dzb_dt,dt_max] = dzb_dt(obj,t,zb,ddir)
	p          = obj.sediment.p;
	rho        = obj.sediment.rho;

	% TODO quick hack, set a switch into RT, if zb is variable or not
	obj.fun.zb = @(xi) interp1(mid(obj.x),zb,xi,'linear','extrap');


	% flow
	obj.init();
	y0 = obj.solve();

	if (obj.out.cflag ~= 1)
		error('no-convergence');
	end

	% set initial condition for reuse in next time step
	%obj.opt.ifun = @(x) y0;
	obj.opt.ifun = @(x) obj.out.ypm;

	dx  = diff(obj.x);

	% transport
	% TODO store
	x = obj.x;
	xc     = mid(x);
	h      = obj.h0(xc);
	w      = obj.width(xc);
	Qs     = obj.sediment_transport(t,ddir);
%	Qs(end) = obj.bc_Qs;


		% differences
		switch (ddir)
		case {1}
			dQs = upwind_difference(x,Qs);
		case {-1}
			dQs = downwind_difference(x,Qs);
		otherwise % central
			dQs = central_difference(x,Qs);
		end

		% note that even if differences are upwinded, dx stays the same
		x   = obj.x;
		xc  = mid(x);
		wc  = obj.width(xc);
		hc  = obj.h0(xc);
		dx = diff(x);

		% scale of mass to elevation
		scale = -1./(wc.*rho.*p);

		% derivative
		dzb_dt = scale.*(dQs./dx);
		
		% time step limit
		c = transport_sensitivity_engelund_hansen(Qs(2:end-1),hc,wc);
		dt_max = min(abs(dx./c));

end % dzb_dt

