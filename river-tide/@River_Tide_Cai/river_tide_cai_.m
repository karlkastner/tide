% Tue 21 Mar 12:09:23 CET 2017
% Karl Kastner, Berlin
%% determine the surface amplitude of the river-tide
%% c.f. Cai
function [X, Y] = river_tide_cai_(val0,omega,cd,Uscale,Xlim,zbfun,wfun,Qrfun)

	% acceleration by gravity
	g = 9.81;

	% TODO no magic numbers
	% storage ratio
	rs = 1;

	% solve system of ODEs
	[X Y] = ode23(@dval_dx, Xlim, val0);

function dval_dx = dval_dx(x,val)
	% tidal amplitude
	eta  = val(1);

	% mean surface elevation
	zs   = val(2);

	% get space dependent variables, that are independent of solution
	% bed level
	[zb  dzb_dx] = feval(zbfun, x);
	
	% width
	[w dw_dx]  = feval(wfun,x);

	% depth (average over tidal cycle)
	h = zs-zb;

	% normalised amplitude
	zeta = eta/h;

	% area
	A = h*w;

	% area convergence
	% TODO account for surface gradient (iteratively)
	dzs_dx = 0;

	% dA/dx = h dw/dx + w dh/dx = h dw/dx + w( dzs/dx + 
	dA_dx = h*dw_dx + w*( dzs_dx - dzb_dx );
	ainv = -1/A*dA_dx;
	a    = 1/ainv;

	% river dircharge, space dependent if there are confluences and distributaries
	Qr = feval(Qrfun,x);

	% river flow velocity
	Ur = Qr/A;

	% velocity scale
	phi = Ur/Uscale;
	if (phi > 1)
		error('Velocity scale too low');
	end

	% free celerity
	% eq. 8
	c0    = sqrt(g*h/rs);

	% friction number
	% tab 1
	chi0 = rs*g*c0/omega*cd/h;

	if ((4/3*zeta) >= 1)
		dval_dx = NaN*val;
		warning('depth too shallow');
		return;
	end

	chi  = chi0*zeta/(1 - (4/3*zeta)^2);

	% friction coefficient
%	p     = friction_coefficient_dronkers(phi);

	% estuarine shape number
	% c.f. cai 2016, tab 1 (here asymptotic area of 0 is assumed)
	gamma = c0/(omega*a);


	[delta lambda mu] = obj.rt_quantities(zeta,gamma,chi,rs,phi);

	% tidal velocity
	v = mu*rs*zeta*c0;

	% derivative of normalised tidal amplitude
	deta_dx = omega/c0*delta*eta;

	% derivative of the mean water level
	% eq 17 in cai 2016
	switch (lower(obj.mode))
	case {'dronkers'}
		p   = friction_coefficient_dronkers(phi); 
		Ft  = cd/(h*pi)*(1/2*p(3) + p(1))*v^2;
		Fr  = cd/(h*pi)*(p(3)-p(4)*phi)*Ur^2;
		Frt = cd/(h*pi)*(-p(2) - 3/2*p(4))*v*Ur;
		F = Ft + Fr + Frt;
	case {'lorentz'}
		L = friction_coefficient_lorentz(phi); 
		F = cd/h*(1/4*L(1)*v^2 + 1/2*L(2)*v*Ut);
	end

	dzs_dx = -F;

	dval_dx = [ deta_dx;
                    dzs_dx ];

end % deta_zs_dx_cai

end % amplitude_mwl_dx_cai

