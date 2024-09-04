% Wed  9 Oct 15:23:10 PST 2019
function [out, rt, d3d] = test_river_tide_hydrodynamics_01(rt_map,pflag)
	meta = test_river_tide_metadata();
	if (nargin()<1)
		rt_map      = River_Tide_Hydrodynamics_Map(meta.mapname_str);
	end
	if (nargin()<2)
		pflag = 1;
	end
	out.id   = 1;

	[rt, in] = hydrodynamic_scenario_from_table(rt_map, meta.rtspecfile_str, out.id, meta.opt); 

	% generate d3d equivalent model for comparison
	d3dopt               = struct();
	d3dopt.Lc            = tab.Lc(fdx);
	d3dopt.bndisharmonic = true;
	folder = [meta.folder.d3d,num2str(out.id)];
	rt.sediment = [];
	rt.generate_delft3d(folder,meta.param,meta.param_silent,d3dopt);
	[out.rmse_d3d, d3d] = test_rt_d3d_evaluate(rt,out.id,pflag);

	Xi = rt.hydrosolver.xi;

	% check ode
	rmse = rt.check_continuity();

	% compare to analytical solution
	g = Constant.gravity;
	c0 = sqrt(g*h0);
	k0 = omega/c0;
	x = rt.channel(1).x;
	z = z10*exp(-1i*k0*x);

	%rmse(2)  = rms(rt.channel(1).z(:,2)-z);
	rmse(2)  = rms(rt.channel(1).waterlevel(1)-z);
	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
	nres_ = 1/2*rms(cdiff(z,2))
	% 1/2*rms(cdiff(z,2)./cdiff(rt.x).^2)*dx^2
	% note : this is less than
	dx = (Xi(2)-Xi(1))/(opt.nx-1);
	nres__ = 1/2*z10*k0^2*dx^2

	result = (rmse(1)/z10 > 0.05) || (rmse(2) > 10*nres_);

	zz = NaN(size(rt.channel(1).z));
	zz(:,2) = z;
	river_tide_test_plot(out.id,rt,zz,out.name,pflag);

	out.rmse   = rmse;
	out.result = result;
end % test_river_tide_hydrodynamics_01

