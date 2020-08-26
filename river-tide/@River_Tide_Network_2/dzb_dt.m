% Thu 20 Aug 10:25:58 +08 2020
%
% function dzb_dt = dzb_dt(obj,t,zb)
%
% TODO, store subiteration count
function [dzb_dt,dt_max] = dzb_dt(obj,t,zb,ddir)
	p          = obj.rt(1).sediment.p;
	rho        = obj.rt(1).sediment.rho;

	% TODO quick hack, set a switch into RT, if zb is variable or not
	nci = 0;
	for idx=1:length(obj.rt)
		x   = obj.rt(idx).x;
		xc  = mid(x);
		nxc = length(xc);
		obj.rt(idx).fun.zb = @(xi) interp1(xc,zb(nci+(1:nxc)),xi,'linear','extrap');
		nci = nci+nxc;
	end % for idx

	% allocate storage
	dzb_dt = zeros(nci,1);

	% determin discharge values and weights for integration over the seasons
	if (isa(obj.season.Qlim,'function_handle'))
		Qlim = obj.season.Qlim(t);
	else
		Qlim = obj.season.Qlim;
	end
	if (obj.season.iorder > 1)
		% integrate over the seasons,
		% without updating the bed-level during the year
		[weight, b] = int_1d_gauss(obj.season.iorder);
		Q0          = inv_hydrograph(  b*[0;1] ...
			 		     , Qlim(1) ...
					     , Qlim(2) ...
					    );
	else
		% single formative discharge, for which Qs ~ <Qs> ~ <Q^2>
		weight = 1;
		Q0     = formative_discharge(  Qlim(1) ...
					     , Qlim(2) ...
					    );
	end
	h0 = normal_flow_depth(Q0,obj.season.w0,obj.season.Cd,obj.season.S0,'Cd');
	Qs = total_transport_engelund_hansen(drag2chezy(obj.season.Cd), ...
		obj.rt(1).sediment.d_mm,Q0./(h0*obj.season.w0),h0,obj.season.w0);

	% compute seasonally averaged rate of bed level change
	c = 0;
	dt_max = inf;
%figure(1);
%clf
	for sdx=1:length(Q0)

	% TODO, this is a quick fix
	obj.rt(1).bc(2,1).rhs  = Q0(sdx);
	obj.rt(1).bc_Qs(2).val = Qs(sdx);

	obj.init();

	% reuse solution of last time step as initial condition
	for idx=1:length(obj.rt)
		obj.rt(idx).opt.ifun = @(x) obj.tmp(idx).ypm(:,sdx);
	end

	obj.solve();
	
	z0 = obj.rt(1).z(0);
%plot(z0)


	% TODO, the cflag belong to the network, not branches
	if (obj.rt(1).out.cflag ~= 1)
		error('no-convergence');
	end

	% save initial condition for reuse during next time step
	for idx=1:length(obj.rt)
		obj.tmp(idx).ypm(:,sdx) = obj.rt(idx).out.ypm;
	end

	% get transport at segment interfaces
	s = obj.sediment_transport(t,ddir);

	% derivative
	nci = 0;
%	dt_max_ = inf;
	for idx=1:length(obj.rt)
		% differences
		switch (ddir)
		case {1}
			dQsi = upwind_difference(x,s(idx).Qs);
		case {-1}
			dQsi = downwind_difference(x,s(idx).Qs);
		otherwise % central
			dQsi = central_difference(x,s(idx).Qs);
		end

		% note that even if differences are upwinded, dx stays the same
		x    = obj.rt(idx).x;
		xc   = mid(x);
		nxc  = length(xc);
		wc   = obj.rt(idx).width(xc);
		hc   = obj.rt(idx).h0(xc);
		dxi  = diff(x);

		% scale of mass to elevation
		scale = -1./(wc.*rho.*p);

		% derivative
		dzb_dt(nci+(1:nxc)) =   dzb_dt(nci+(1:nxc)) ...
				      + weight(sdx).*scale.*(dQsi./dxi);
		nci = nci+nxc;
		
		% time step limit
		c = transport_sensitivity_engelund_hansen(weight(sdx)*s(idx).Qs(2:end-1),hc,wc);
%		dt_max_ = min(dt_max_,min(abs(dxi./c)));
		dt_max = min(dt_max,min(abs(dxi./c)));
%		s(idx).c = s(idx)c 
	end % for idx
	% TODO, actually c should be added
%	dt_max = dt_max + weight(sdx)*dt_max_
	
%hold on
	end % for sdx
%pause
end % dzb_dt

