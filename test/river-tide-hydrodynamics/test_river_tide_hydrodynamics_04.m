% Wed  9 Oct 15:23:10 PST 2019
function [out,rt] = river_tide_test_04(rt_map,pflag)
	if (nargin()<2)
		pflag = 1;
	end
	out.id = 4;
	out.name = 'wave along frictionless prismatic channel, reflection (standing wave)';

	% surface elevation at channel mouth
	z10       = 1;
	zs	  = [0,z10];

	% river discharge
	Q0        = 0;
	% width of channel
	w0        = 1;
	% drag/friction coefficient
	Cd	  = 0;
	% depth of channel at channel mouth
	h0        = 10;
	% bed level
	zb        = -h0;
	% base frequency
	T         = Constant.SECONDS_PER_DAY;
	omega     = 2*pi/T;

	% length of computational domain
	Lx        = 1e6;

	meta = river_tide_test_metadata();
	opt  = meta.opt;

	% reflection coefficient at right end of boundary
	qr = 1;

	rt = hydrodynamic_scenario(rt_map,zs,qr,zb,Q0,w0,Cd,omega,Lx,opt);

	Xi = rt.hydrosolver.xi;

	% check ode
	rmse = rt.check_continuity();

	% compare to analytical solution
	g = Constant.gravity;
	c0 = sqrt(g*h0);
	k0 = omega/c0;
	x = rt.x;
	z = z10*exp(-1i*k0*x) ...
	    - z10*exp(-1i*k0*x(end) - 1i*k0*(x(end)-x));

	rmse(2)  = rms(rt.z(1)-z);
	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
	nres_ = rms(cdiff(z,2));
	result = (rmse(1)/z10 > 0.05) || (rmse(2) > 10*nres_);
		
	zz      = NaN(size(rt.out.z));
	zz(:,2) = z;
	river_tide_test_plot(out.id,rt,zz,out.name,pflag);

	out.rmse   = rmse;
	out.result = result;
end % river_tide_test_04

