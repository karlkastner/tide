% 2017-04-08 15:16:22.323048338 +0200
%% linearly damped wave in rectangular channel
%% x_t = Ax + b
function [x, z] = damped_wave_ivp(L)
	% TODO no magic numbers
	g  = 9.81;
	T  = 86400;
	o  = 2*pi/T;
	cD = 5e-4;
	H  = 10;
	z0 = 1;
	u0 = z0*sqrt(g/H);

	% damping coefficient according to lorentz
	r  = -cD*8/(3*pi)*u0/H;
	% stabilisation
	s = 0; %-1e-6;

	a = 1/(H*g)*(-o^2 + r/H*1i*o)

	% blowing off as soon as argument(a) is > 0 or has non-zero imaginary part
	% a = -1e-4*(1-5e-3i)
	if (0)
		a = -1e-6;
		c0 = sqrt(g*H)
		a = -(o/c0)^2
		a = a*(1 - 0.1*i)
	end
	
	% initial slope
	dz0 = [     -4.078042362031606e-06 + 7.947489460229805e-06i]
	
	if (1)
		% solution in z and u
		A  = [0, 1/g*(-r+1i*o);
		      1i*o/H, 0 ];
		A = A+s*eye(2);
		%eig(A)
		y0 = [z0;u0];
	else
		% solution in z and dz_dx
		A = [0  1;  % dz = y
		     a, 0]; % dz = ...
		A = A - s*eye(2);
		y0 = [z0;dz0];
	end
	%	[V E] = eig(A);
	%	E
	%pause

	
	if (0)
		n = 1e2;
		x = linspace(0,L,n)';
		for idx=1:length(x)
			y(:,idx) = expm(A*x(idx))*y0;
		%	y(:,idx)       = (V*diag(exp(diag(E)*x(idx)))*V')*y0;
		end
	else
		% solve ivp
		% ode23 may fail sometimes
		[x y] = ode45(@dot,[0,L],y0);
	end
	z_ = y(1,:)';
	z = NaN*z0*1/2*(exp(sqrt(a)*x)-exp(-sqrt(a)*x));
	%z = z0/sqrt(a)*exp(1/2*a*x);
	%z(:,3) = 0*sin(sqrt(a)*x) + z0*cos(sqrt(-a)*x);

	
	
	function d = dot(x,y)
		d = A*y;
	%	d(1,1) = +1/g*(-r + 1i*o)*y(2);
	%	d(2,1) = s+1/H*1i*o*y(1);
	end
end % damped_wave_ivp

