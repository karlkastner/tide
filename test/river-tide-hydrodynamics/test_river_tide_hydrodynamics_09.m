% Wed  9 Oct 15:23:10 PST 2019
function [out, rt, d3d] = test_river_tide_hydrodynamics_09(rt_map,pflag)
	meta = test_river_tide_metadata();
	if (nargin()<1)
		rt_map      = River_Tide_Hydrodynamics_Map(meta.mapname_str);
	end
	if (nargin()<2)
		pflag = 1;
	end
	tab = readtable('test-river-tide.csv');
	out.id   = 9;
	fdx = find(tab.id == out.id)
	out.name = tab(fdx,:).name{1};

	% mean surface elevation
	z00 = 0;

	% tidal surface elevation
	z10 = tab.z10(fdx);
	zs = [0,z10,0,0,0];

	% river discharge
	Q0 = tab.Q0(fdx);

	% width of channel at river mouth
	w00 = tab.w00(fdx);

	% width of channel
	w0  = eval(tab(fdx,:).w0{1});

	% drag/friction coefficient
	Cd = tab.Cd(fdx);

	% depth of channel
	h0 = tab.h0(fdx);

	% slope of channel bed
	S0        = -normal_flow_slope(Q0,h0,w00,drag2chezy(Cd));

	% bed level of channel
	zb        = eval(tab(fdx,:).zb{1});

	% base frequency
	T_d       = tab.T(fdx);
	T         = T_d*Constant.SECONDS_PER_DAY;
	omega     = 2*pi/T;

	% length of computational domain
	Lx = tab.Lx(fdx);

	% reflection coefficient at right end of boundary
	ql = tab.ql(fdx);
	qr = tab.qr(fdx);

	meta = test_river_tide_metadata();
	opt = meta.opt;

	rt = hydrodynamic_scenario(rt_map,zs,ql,qr,zb,Q0,w0,Cd,omega,Lx,opt);

	% generate d3d equivalent model for comparison
	d3dopt                = struct();
	d3dopt.Lc            = tab.Lc(fdx);
	d3dopt.bndisharmonic = true;
	folder = [meta.folder.d3d,num2str(out.id)];
	rt.generate_delft3d(folder,meta.param,meta.param_silent,d3dopt);
	[out.rmse_d3d, d3d] = test_rt_d3d_evaluate(rt,out.id,pflag);

	Xi = rt.hydrosolver.xi;

	% compare to analytical solution
	g = Constant.gravity;
	c = sqrt(g*h0);
	k = omega/c;
	x = rt.channel(1).x;

%	bw = Backwater1D();
%	nn = opt.nx;
%	[x_, h_, z0_] = bw.solve(-Q0,0,drag2chezy(Cd),wfun,zbfun,0,Xi);
%	[x_, h_, z0_] = bw.solve_analytic(-Q0,drag2chezy(Cd),w0,S0,h0,nn);
%	z0  = interp1(x_,z0_,rt.channel(1).x,'linear','extrap');
%	z0 = S0*x;
%	z0t = rt.channel(1).waterlevel(0) - z0;
	z2 = rt.channel(1).waterlevel(2);
	z2_ = rt.even_overtide_analytic(x,z10,h0,w0,abs(Q0),Cd,omega);
	%z0_ = interp1(x_,h_,rt.channel(1).x,'spline')+0*zbfun(rt.channel(1).x);
	% r = (1+1i)*sqrt(-Cd.*omega.*Q0/w0./(g*h0.^3));
	% z = z10*exp(-r*x);

	rmse(2)  = norm(z2_-z2);

	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
	nres_ = rms(cdiff(z2_,2))
	result = (rmse(1) > 0.05/z10) || (rmse(2) > 10*nres_);

	out.rmse   = rmse;
	out.result = result;

	zz = NaN(size(rt.channel(1).z));
	zz(:,3) = z2_;
	river_tide_test_plot(out.id,rt,zz,out.name,pflag);
	
end % river_tide_test_09

