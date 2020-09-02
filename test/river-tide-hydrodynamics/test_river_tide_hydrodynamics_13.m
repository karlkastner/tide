% Wed  9 Oct 15:23:10 PST 2019
function [out,rt] = river_tide_test_13(rt_map,pflag)
	if (nargin()<2)
		pflag = 1;
	end
	out.id   = 13;
	out.name = 'quarter-diurnal tide, uniform flow';

	% frequency components surface elevation at channel mouth
	% (z0,z1,z2,z3,z4)
	z10 = 0.1;
	zs = [0,z10,0,0];
 
	% river discharge
	Q0        =  -10;
	% width of channel
	w0        = 1;
	% drag/friction coefficient
	Cd        = 2.5e-3;
	% depth of channel
	h0        = 5;
	% slope of channel bed
	S0         = -normal_flow_slope(Q0,h0,w0,drag2chezy(Cd));
	% bed level of channel
	zb     = @(x) -h0 + S0*x;


	% base frequency
	T         = Constant.SECONDS_PER_DAY;
	omega     = 2*pi/T;

	% domain length
	Lx	  = 2e5;

	meta = river_tide_test_metadata();
	opt = meta.opt;

	opt.oflag   = [true(1,4)];

	% reflection coefficient at right boundary
	qr = 0;
	rt = hydrodynamic_scenario(rt_map,zs,qr,zb,Q0,w0,Cd,omega,Lx,opt);

	% compare to analytical solution
	g = Constant.gravity;
	c = sqrt(g*h0);
	k = omega/c;
	x = rt.x;

%	bw = Backwater1D();
%	nn = opt.nx;
%	[x_, h_, z0_] = bw.solve(-Q0,0,drag2chezy(Cd),wfun,zbfun,0,Xi);
%	[x_, h_, z0_] = bw.solve_analytic(-Q0,drag2chezy(Cd),w0,S0,h0,nn);
%	z0  = interp1(x_,z0_,rt.x,'linear','extrap');
%	z0 = S0*x;
%	z0t = rt.z(0) - z0;
	z2 = rt.z(2);

	z2_ = rt.even_overtide_analytic(x,z10(1),h0,w0,abs(Q0),Cd,omega);
	%z0_ = interp1(x_,h_,rt.x,'spline')+0*zbfun(rt.x);
	% r = (1+1i)*sqrt(-Cd.*omega.*Q0/w0./(g*h0.^3));
	% z = z10*exp(-r*x);

	rmse(2)  = inf; %norm(z2_-z2);

	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
	nres_ = rms(cdiff(z2_,2))
	result = (rmse(1) > 0.05/z10) || (rmse(2) > 10*nres_);

	out.rmse   = rmse;
	out.result = result;

	zz = NaN(size(rt.out.z));
	zz(:,3) = z2_;
	river_tide_test_plot(out.id,rt,zz,out.name,pflag);
	
end % river_tide_test_09

