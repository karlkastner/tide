% Fri  2 Sep 14:31:13 CEST 2022
function [rt, d3d] = river_tide_nonstationary_01(rt_map,pflag)
	meta = river_tide_nonstationary_metadata();

	if (nargin()<1)
		rt_map      = River_Tide_Hydrodynamics_Map(meta.mapname_str);
	end

	if (nargin()<2)
		pflag = 1;
	end

	for id=302 %[0:9] %,100]

	rt = hydrodynamic_scenario_from_table(rt_map, meta.rtspec_str, id, meta.opt); 

	% generate d3d equivalent model for comparison
	d3dopt               = struct();
	d3dopt.Lc            = rt.channel(1).x(end);
	folder = [meta.folder.d3d,num2str(id)];
	param = meta.param;
% TODO, warn if params are overwritten
	param_silent = meta.param_silent;
% with +1, falls dry for large amplitude (7)
	param.z0 = max(rt.channel(1).zb)+5;
	
	param_silent.mdf.Tstop = (20+10)*15*1440;

	[d3d,param,param_silent] = rt.generate_delft3d(folder,param,param_silent,d3dopt);
%	[out.rmse_d3d, d3d] = test_rt_d3d_evaluate(rt,out.id,pflag);

	% lake
	X        = d3d.mesh.X;
	X(3:4,:) = X(2,:) + 1e4*(1:2)';
	Y        = d3d.mesh.Y;
	Y(3:4,:) = ones(2,1)*Y(1,:);
	
	d3d.mesh.X = X;
	d3d.mesh.Y = Y;

	folder = [folder,'-lake'];
        d3d    = d3d.mesh.generate_delft3d(folder,param,param_silent);                 

	end

if (0)
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
end	
end % river_tide_test_09

