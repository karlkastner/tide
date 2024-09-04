% Mon  7 Sep 16:02:26 +08 2020
function [out, rt, d3d] = test_river_tide_hydrodynamics_15(rt_map,pflag)
	meta = test_river_tide_metadata();
	if (nargin()<1 || isempty(rt_map))
		rt_map      = River_Tide_Hydrodynamics_Map(meta.mapname_str);
		rt_map.recompute = true;
	end
	if (nargin()<2)
		pflag = 1;
	end
	out.id   = 15;
	[rt, in] = hydrodynamic_scenario_from_table(rt_map, meta.rtspecfile_str, out.id, meta.opt); 
%	tab      = readtable(meta.testspec);
%	fdx      = find(tab.id == out.id)
%	out.name = tab(fdx,:).name{1};
%
%	% surface elevation at channel mouth
%	z10      = tab.z10(fdx);
%
%	zs = [0,z10,0,0];
%
%	% river discharge
%	Q0 = tab.Q0(fdx);
%
%	% width at channel mouth
%	w00 = tab.w00(fdx);
%
%	% width of channel
%	w0  = eval(tab(fdx,:).w0{1});
%
%	% drag/friction coefficient
%	Cd = tab.Cd(fdx);
%
%	% depth at channel mouth
%	h0 = tab.h0(fdx);
%
%	% slope of channel bed
%	S0         = -normal_flow_slope(Q0,h0,w00,drag2chezy(Cd));
%
%	% bed level of channel
%	zb        = eval(tab(fdx,:).zb{1});
%
%	% base frequency
%	T_d       = tab.T(fdx);
%	T         = T_d*Constant.SECONDS_PER_DAY;
%	omega     = 2*pi/T;
%
%	% length of computational domain
%	Lx = tab.Lx(fdx);
%
%	% reflection coefficient at right end of boundary
%	ql = tab.ql(fdx);
%	qr = tab.qr(fdx);
%
%	opt = meta.opt;
%	rt = hydrodynamic_scenario(rt_map,zs,ql,qr,zb,Q0,w0,Cd,omega,Lx,opt);

	% generate d3d equivalent model for comparison
	meta.param_silent.mdf.Sub2   = ' C ';   % activate transport
	d3dopt                = struct();
	d3dopt.Lc            = tab.Lc(fdx);
	d3dopt.bndisharmonic = true;
	folder = [meta.folder.d3d,num2str(out.id)];
	rt.generate_delft3d(folder,meta.param,meta.param_silent,d3dopt);
	rt.opt.stokes_order = 2;
	[out.rmse_d3d, d3d] = test_rt_d3d_evaluate(rt,out.id,pflag);

	% check sediment transport	


%	Xi = rt.hydrosolver.xi;
%
%	% check ode
%	rmse = rt.check_continuity();
%
%	% compare to analytical solution
%	g = Constant.gravity;
%	c0 = sqrt(g*h0);
%	k0 = omega/c0;
%	x = rt.channel(1).x;
%
%	r = (1+1i)*sqrt(-Cd.*omega.*Q0/w0./(g*h0.^3));
%	z = z10*exp(-r*x);
%
%	rmse(2)  = rms(rt.channel(1).waterlevel(1)-z)
%	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
%	nres_ = rms(cdiff(z,2));
%	result = (rmse(1)/z10 > 0.05) || (rmse(2) > 10*nres_);
%
	zz      = NaN(size(rt.channel(1).z));
%	zz(:,2) = z;
	river_tide_test_plot(out.id,rt,zz,out.name,pflag);
%
%	out.rmse   = rmse;
%	out.result = result;

end % test_river_tide_hydrodynamics_15
