% 2018-04-25 14:14:53.330209405 +0200
% Karl Kastner, Berlin
%
function z0t = rt_mwl_simplified(x,Q0,Q1,z10,h0,w,cd,omega,dzb_dx)
	g = obj.g;
	X = [x(1),x(end)];
	Q1_  = @(xi) interp1(x,Q1,xi,'spline')
	h0_  = @(xi) interp1(x,h0,xi,'spline');
	ode = @(x,z0t) -cd*(   6*Q0^2*z0t ...
			   - abs(Q1_(x))^2*h0_(x) ...
			   + 0*3*abs(Q1_(x))^2*z0t ...
			) / (2*g*h0_(x)^4*w^2);
	ode = @(x,z0t) ( cd*( ...
					- 6*Q0^2*h0_(x)*z0t ...
					+ 12*Q0^2*z0t^2 ...
					+ abs(Q1_(x))^2*h0_(x)^2 ...
					- 3*abs(Q1_(x))^2*h0_(x)*z0t ...
					+ 6*abs(Q1_(x))^2*z0t^2) ...
				) /(2*g*h0_(x)^5*w^2)
		%ode = @(x,z0t) ((cd*(Q0^2 + abs(Q1(x))^2/2))/(w*(h0_(x) + z0t)^2) - (Q0^2*cd*(h0_(x) + z0t))/(h0_(x)^3*w))/(g*w*(h0_(x) + z0t));
		ode = @(x,z0t) -(cd*( ...
				  6*Q0^2*h0_(x)^2*z0t ...		% important
				+ 0*6*Q0^2*h0_(x)*z0t^2 ...		% unimportant
				+ 0*2*Q0^2*z0t^3 ...			% unimportant
				- abs(Q1_(x))^2*h0_(x)^3) ...		% important
				) / (2*g*h0_(x)^3*w^2*( ...
					h0_(x) ...
					+ 0*z0t ...			% unimportant
					)^3);
%	set(gca,'colororderindex',get(gca,'colororderindex')-1);
	sopt = struct();
	sopt.InitialStep = x(2)-x(1);
	sopt.MaxStep = 4*(x(2)-x(1));
	[x_, z0t] = ode23s(ode,[X(1),0.2*X(2)],0,sopt);
	z0t = interp1(x_,z0t,x,'spline');
end

