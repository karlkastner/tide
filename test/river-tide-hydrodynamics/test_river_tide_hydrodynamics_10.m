% Wed  9 Oct 15:55:23 PST 2019
function [out, rt, d3d] = test_river_tide_hydrodynamics_10(rt_map,pflag)
	meta = test_river_tide_metadata();
	if (nargin()<1)
		rt_map      = River_Tide_Hydrodynamics_Map(meta.mapname_str);
	end
	if (nargin()<2)
		pflag = 1;
	end
	tab = readtable('test-river-tide.csv');
	out.id   = 10;
	fdx = find(tab.id == out.id)
	out.name = tab(fdx,:).name{1};

	% tidal surface elevation
	z10 = tab.z10(fdx);
	zs = [0,z10];

	% river discharge
	Q0 = tab.Q0(fdx);

	% width of channel at river mouth
	w00 = tab.w00(fdx);

	% width of channel
%	w0  = eval(tab(fdx,:).w0{1});

	% drag/friction coefficient
	Cd = tab.Cd(fdx);

	% depth of channel
	h0 = tab.h0(fdx);

	% bed level of channel
	zb        = eval(tab(fdx,:).zb{1});

	% base frequency
	T_d       = tab.T(fdx);
	T         = T_d*Constant.SECONDS_PER_DAY;
	omega     = 2*pi/T;

	g = Constant.gravity;
	c = sqrt(g*h0);
	k0 = omega/c;

	% damping by tide
	k = wave_number_tide(omega,Cd,h0,abs(z10));

	% damping modulus
	r1 =  -imag(k);

	% convergence length of width
	Lw = 0.5/r1;
	w0  = eval(tab(fdx,:).w0{1});

	% length of computational domain
	Lx = tab.Lx(fdx);

	% reflection coefficient at right end of boundary
	qr = tab.qr(fdx);

	opt = meta.opt;

	rt = hydrodynamic_scenario(rt_map,zs,qr,zb,Q0,w0,Cd,omega,Lx,opt);

	% generate d3d equivalent model for comparison
	rt.generate_delft3d(out.id,meta.param_silent,tab.Lc(fdx));
	[out.rmse_d3d, d3d] = test_rt_d3d_evaluate(rt,out.id,pflag);

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

