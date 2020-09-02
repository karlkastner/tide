% Wed  9 Oct 15:55:23 PST 2019
function [out,rt] = river_tide_test_10(rt_map,pflag)
	if (nargin()<2)
		pflag = 1;
	end
	out.id = 10;
	out.name = 'ideal estuary, width convergence cancel frictional damping';
	% tidal surface elevation
	z10 = 0.1;
	zs = [0,z10];
	% river discharge
	Q0        = 0;
	% width of channel
	w00        = 1e5;
	% depth of channel
	h0        = 10;
	% drag/friction coefficient
	Cd        = 2.5e-3;

	% base frequency
	T         = Constant.SECONDS_PER_DAY;
	omega     = 2*pi/T;

	g = Constant.gravity;
	c = sqrt(g*h0);
	k0 = omega/c;

	% damping by tide
	k = wave_number_tide(omega,Cd,h0,abs(z10));

	% damping modulus
	r =  -imag(k);

	% convergence length of width
	Lw = 0.5/r;
	w0      = @(x) w00*exp(-x/Lw);
	% bed level of channel
	zb        = -h0;

	Lx = 1e6;

	meta = river_tide_test_metadata();
	opt = meta.opt;
	% reflection coefficient at right end of boundary
	qr = 0;
	rt = hydrodynamic_scenario(rt_map,zs,qr,zb,Q0,w0,Cd,omega,Lx,opt);

	Xi = rt.hydrosolver.xi;

	% check ode
	rmse = rt.check_continuity();

	% compate to analytical solution
	x = rt.x;
	z = z10*exp(-1i*k*x + 0.5./Lw*x);

	rmse(2)  = rms(rt.z(1)-z);
	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
	nres_ = rms(cdiff(rt.z(1),1));
	result = (rmse(1) > 0.05/z10) || (rmse(2) > 10*nres_);

	out.rmse   = rmse;
	out.result = result;

	zz = NaN(size(rt.out.z));
	zz(:,2) = z;
	river_tide_test_plot(out.id,rt,zz,out.name,pflag);
	
end % river_tide_test_10

