% Wed 18 Apr 13:24:05 CEST 2018
% TODO, this should go to the IVP class
%% offset of the tidally averaged surface elevation caused by tidal friction
%% Linear estimate of the mean water level offset (ignoring feed-back of tide)
function z0t = mwl_offset(obj)
	opt = struct();
	opt.MaxStep = 4*median(diff(obj.x));
	% TODO normalize domain lenght
	[x_,z0t_] = ode23s(@fun_,obj.Xi,0,opt);
	z0t = interp1(x_,z0t_,obj.x);

	function dz0t_dx = fun_(x_,z0t)
		dz0_dx   = obj.dz0_dx(cdx,x_);
		w        = obj.width(cdx,x_);
		%h0       = interp1(obj.x, obj.h0, x_);
		h0       = obj.h0(cdx,x_);
		%Q0       = obj.Q(0,x_);
		Q0       = obj.Q(0,cdx,x_);
		Q1       = obj.Q(1,cdx,x_);
		aQ1      = abs(Q1);
%interp1(obj.x, abs(obj.Q(1)), x_);
%		aQ1      = interp1(obj.x, abs(obj.Q(1)), x_);
		Qhr      = aQ1;
		cd       = obj.cd(cdx,x_);
		p        = obj.friction_coefficient_dronkers(abs(Q0)/Qhr);
		A0	 = w.*h0;
		g        = obj.g;

		% TODO, friction coefficient
%		dz0t_dx = (cd./g).*w./(pi*A0.^3).*(p(1) + p(2)*Q0*Qhr + p(3)*(Q0.*abs(Q0) + 1/2*abs(Q1).^2) ) ...
%                          - (z0t./h0).*dz0_dx;
%		dz0t_dx = -(cd./g).*w./(pi*A0.^3).*(p(2)*Q0*Qhr + p(3)*(0*Q0.*abs(Q0) + 1/2*abs(Q1).^2) ) ...
 %                         - (z0t./h0).*dz0_dx;
		% linearized:
%		dz0t_dx = -(cd*(6*Q0^2*z0t - Q1^2*h0 + 3*Q1^2*z0t))/(2*g*h0^4*w^2);
		% not linearized:
		dz0t_dx = -(cd*(  6*Q0^2*h0^2*z0t ...
			        + 6*Q0^2*h0*z0t^2 ...
				+ 2*Q0^2*z0t^3 ...
				- aQ1^2*h0^3)) ... 
			    / (2*g*h0^3*w^2*(h0 + z0t)^3);
	end % dz0t_dx
end % River_Tide/mwl_offset

