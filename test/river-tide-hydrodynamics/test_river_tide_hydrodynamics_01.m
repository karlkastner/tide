% Wed  9 Oct 15:23:10 PST 2019
function [out, rt, d3d] = test_river_tide_hydrodynamics_01(rt_map,pflag)
	meta = test_river_tide_metadata();
	if (nargin()<1)
		rt_map      = River_Tide_Hydrodynamics_Map(meta.mapname_str);
	end
	if (nargin()<2)
		pflag = 1;
	end
	tab = readtable('test-river-tide.csv');
	out.id   = 1;
	fdx = find(tab.id == out.id)
	out.name = tab(fdx,:).name{1};

	% surface elevation at channel mouth
	z10 = tab.z10(fdx);

	zs       = [0,z10];

	% river discharge
	Q0 = tab.Q0(fdx);

	% width at channel mouth
	w00 = tab.w00(fdx);

	% width of channel
	w0  = eval(tab(fdx,:).w0{1});

	% drag/friction coefficient
	Cd = tab.Cd(fdx);

	% depth at channel mouth
	h0 = tab.h0(fdx);

	% slope of channel bed
	S0         = -normal_flow_slope(Q0,h0,w00,drag2chezy(Cd));

	% bed level of channel
	zb        = eval(tab(fdx,:).zb{1});

	% base frequency
	T_d       = tab.T(fdx);
	T         = T_d*Constant.SECONDS_PER_DAY;
	omega     = 2*pi/T;

	% length of computational domain
	Lx = tab.Lx(fdx);

	% reflection coefficient at right end of boundary
	qr = tab.qr(fdx);

	opt  = meta.opt;

	rt = hydrodynamic_scenario(rt_map,zs,qr,zb,Q0,w0,Cd,omega,Lx,opt);

	% generate d3d equivalent model for comparison
	rt.generate_delft3d(out.id,meta.param_silent,tab.Lc(fdx));
	[out.rmse_d3d, d3d] = test_rt_d3d_evaluate(rt,out.id,pflag);

	Xi = rt.hydrosolver.xi;

	% check ode
	rmse = rt.check_continuity();

	% compare to analytical solution
	g = Constant.gravity;
	c0 = sqrt(g*h0);
	k0 = omega/c0;
	x = rt.x;
	z = z10*exp(-1i*k0*x);

	%rmse(2)  = rms(rt.out.z(:,2)-z);
	rmse(2)  = rms(rt.z(1)-z);
	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
	nres_ = 1/2*rms(cdiff(z,2))
	% 1/2*rms(cdiff(z,2)./cdiff(rt.x).^2)*dx^2
	% note : this is less than
	dx = (Xi(2)-Xi(1))/(opt.nx-1);
	nres__ = 1/2*z10*k0^2*dx^2

	result = (rmse(1)/z10 > 0.05) || (rmse(2) > 10*nres_);

	zz = NaN(size(rt.out.z));
	zz(:,2) = z;
	river_tide_test_plot(out.id,rt,zz,out.name,pflag);

	out.rmse   = rmse;
	out.result = result;
end % test_river_tide_hydrodynamics_01

