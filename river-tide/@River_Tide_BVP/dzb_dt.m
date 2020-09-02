% Thu 20 Aug 10:25:58 +08 2020
% Karl Kastner, Berlin
%
%% change of bed level over time, when width constant over time
%% dzb/dt + 1/(p rho w) dQs/dx = 0
%
% function dzb_dt = dzb_dt(obj,t,zb)
%
% TODO, store subiteration count
function [dzb_dt,dt_max] = dzb_dt(obj,t,zb,ddir)
	p          = obj.sediment.p;
	rho        = obj.sediment.rho;

	% TODO BVPSC should have global xc
	nci = 1;
	xc = [];
	for cdx=1:obj.nc
		x         = obj.x(cdx);
		xc        = [xc;mid(x)];
		nxc       = length(x)-1;
		ncii(cdx) = nci;
		nci       = nci+nxc;
		obj.tmp(cdx).zb   = [];
		obj.tmp(cdx).c.zb = [];
	end % for idx
	ncii(obj.nc+1) = nci;
	obj.fun.zb = @(cdx,xi) interp1(xc(ncii(cdx):ncii(cdx+1)-1),zb(ncii(cdx):ncii(cdx+1)-1),xi,'linear','extrap');

	% allocate storage
	dzb_dt = zeros(nci-1,1);

	% determin discharge values and weights for integration over the seasons
	% TODO evaluate in time, if given by function
	sbc = struct();

	% for each channel
	for cdx=1:obj.nc
	    % implicit assumption here, that mean-component has frequency index 1
	    % for each end point
	    for bdx=1:2
	        if (~isfield(obj.bc(bdx,1,cdx),'Qseason') || isempty(obj.bc(bdx,1,cdx).Qseason))
			weight = 1;
			sbc(bdx,cdx).rhs    = obj.bc(bdx,1,cdx).rhs*ones(1,obj.opt.iorder);
			sbc_Qs(bdx,cdx).rhs = obj.bc_Qs(bdx,cdx).rhs*ones(1,obj.opt.iorder);
		else
			Qlim = obj.bc(bdx,1,cdx).Qseason;

			if (obj.opt.iorder > 1)
				% discharges and weights for integrating over the seasons
				[weight, b] = int_1d_gauss(obj.opt.iorder);
				Q0          = inv_hydrograph(  b*[0;1] ...
					 		     , Qlim(1) ...
							     , Qlim(2) ...
							    );
			else
				% single formative discharge, for which Qs ~ <Qs> ~ <Q^2>
				% TODO pass manning/chezy flag
				weight = 1;
				Q0     = formative_discharge(  Qlim(1) ...
							     , Qlim(2) ...
							    );
			end
			%h0 = normal_flow_depth(Q0,obj.width(cdx,),obj.cd(cdx,),obj.season.S0,'Cd');
			%Qs = total_transport_engelund_hansen(drag2chezy(obj.season.Cd), ...
			%	obj.sediment.d_mm,Q0./(h0*obj.season.w0),h0,obj.season.w0);
			Qs = obj.bc_Qs(bdx,1,cdx).rhsfun(t,Q0);
			sbc(bdx,cdx).rhs    = Q0;
			sbc_Qs(bdx,cdx).rhs = Qs;
		end % else of if
	    end % for bdx
	end % for cdx

	% compute seasonally averaged rate of bed level change
	c      = 0;
	dt_max = inf;
	% for each integration point (aka season) over the annual hydrograph
	for sdx=1:obj.opt.iorder
		% for each channel
		for cdx=1:obj.nc
		 % for each boundary
		 for bdx=1:2
			% set to seasonal value
			if (~isempty(sbc(bdx,cdx).rhs))
			obj.bc(bdx,1,cdx).rhs    = sbc(bdx,cdx).rhs(sdx);
			obj.bc_Qs(bdx,cdx).rhs   = sbc_Qs(bdx,cdx).rhs(sdx);
			end
		 end % bdx
		end % cdx
	
		% reuse solution of last time step as initial condition
		if (isfield(obj.tmp(1),'ypm') && size(obj.tmp(1).ypm,2)>=sdx)
			obj.hydrosolver.inifun = @(cdx,x) obj.tmp(cdx).ypm(:,sdx);
		end % if

		obj.solve();
		
		% TODO, the cflag belong to the network, not branches
		if (obj.hydrosolver.out(1).cflag ~= 1)
			error('no-convergence');
		end % if
	
		% save initial condition for reuse during next time step
		for cdx=1:obj.nc
			obj.tmp(cdx).ypm(:,sdx) = obj.hydrosolver.out(cdx).ypm;
		end % for cdx
	
		% get transport at segment interfaces
		s = obj.sediment_transport(t,ddir);
	
		% derivative
		nci = 0;
	%	dt_max_ = inf;
		for cdx=1:obj.nc
			% differences
			switch (ddir)
			case {1}
				dQsi = upwind_difference(x,s(cdx).Qs);
			case {-1}
				dQsi = downwind_difference(x,s(cdx).Qs);
			otherwise % central
				dQsi = central_difference(x,s(cdx).Qs);
			end % switch ddir
	
			% note that even if differences are upwinded, dx stays the same
			x    = obj.x(cdx);
			xc   = mid(x);
			nxc  = length(xc);
			wc   = obj.width(cdx,xc);
			hc   = obj.h0(cdx,xc);
			dxi  = diff(x);
	
			% scale of mass to elevation
			scale = -1./(wc.*rho.*p);
	
			% derivative
			dzb_dt(nci+(1:nxc)) =   dzb_dt(nci+(1:nxc)) ...
					      + weight(sdx).*scale.*(dQsi./dxi);
			nci = nci+nxc;
			
			% time step limit
			c = transport_sensitivity_engelund_hansen(weight(sdx)*s(cdx).Qs(2:end-1),hc,wc);
	%		dt_max_ = min(dt_max_,min(abs(dxi./c)));
			dt_max = min(dt_max,min(abs(dxi./c)));
	%		s(idx).c = s(idx)c 
		end % for cdx
		% TODO, actually c should be added
	%	dt_max = dt_max + weight(sdx)*dt_max_
	
	end % for sdx
end % dzb_dt

