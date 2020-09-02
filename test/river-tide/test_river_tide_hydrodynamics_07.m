% Wed  9 Oct 15:23:10 PST 2019
function [out,rt] = river_tide_test_07(rt_map,pflag)
	if (nargin()<2)
		pflag = 1;
	end
	out.id       = 7;
	out.name      = 'backwater curve, no tide';
	% mean surface elevation
	z00 = 0;
	% tidal surface elevation
	% TODO deactivate?
	z10 = sqrt(eps);
	zs =[0,z10];
	% river discharge
	Q0        = -10;
	% width of channel
	w0        = 1;
	% drag/friction coefficient
	Cd        = 2.5e-3;
	% depth of channel at river mouth
	h0        = 10;
	% slope of bed
	S0         = -normal_flow_slope(Q0,h0,w0,drag2chezy(Cd));
	% bed level of channel
	zb         = @(x) -h0 + S0*x;
	% base frequency
	T         = Constant.SECONDS_PER_DAY;
	omega     = 2*pi/T;

	% length of computational domain
	Lx         = 1e6;

	meta = river_tide_test_metadata();
	opt = meta.opt;

	% reflection coefficient at right end of boundary
	qr = 0;
	rt = hydrodynamic_scenario(rt_map,zs,qr,zb,Q0,w0,Cd,omega,Lx,opt);

	Xi = rt.hydrosolver.xi;

	% TODO, check continuity of bw-curve
	rmse(1) = 0;
	% compare to analytical solution
	g = Constant.gravity;
	c0 = sqrt(g*h0);
	k0 = omega/c0;
	x = rt.x;

	bw = Backwater1D();
	nn = opt.nx;
%	[x_, h_, z0_] = bw.solve(-Q0,0,drag2chezy(Cd),wfun,zbfun,0,Xi);
	% TODO, why does it fail?
	dS = S0*sqrt(eps);
	[x_, h_, z0_] = bw.solve_analytic(Q0,drag2chezy(Cd),w0,S0+dS,h0,nn);
	z0  = interp1(x_,z0_,rt.x,'linear','extrap');
	%z0_ = interp1(x_,h_,rt.x,'spline')+0*zbfun(rt.x);
	% r = (1+1i)*sqrt(-Cd.*omega.*Q0/w0./(g*h0.^3));
	% z = z10*exp(-r*x);

	rmse(2)  = rms(rt.z(0)-z0);

	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
	nres_ = rms(cdiff(z0_,2))
	result = (rmse(1) > 0.05) || (rmse(2) > 10*nres_);

	out.rmse   = rmse;
	out.result = result;
	
	if (pflag)
		[Qlr,zlr] = rt.decompose();
		namedfigure(out.id,['Test: ',out.name]);
		clf();
		subplot(2,3,1);
		plot(rt.x,[zb(x),rt.z(0)]);
		hold on;
		plot(rt.x,z0,'--');
		legend('z_b','z_0');

		subplot(2,3,4);
		plot(rt.x,rt.width());
		legend('w_0');

		subplot(2,3,2);
		plot(rt.x,abs(rt.z(1)));
%		hold on;
%		plot(x,abs(z),'--');
		legend('|z_1|');

		subplot(2,3,5)
		plot(rt.x,angle(rt.z(1)));
%		hold on;
%		plot(x,angle(z),'--');
		legend('arg(z_1)');
		ylim(pi*[-1,1]);

		subplot(2,3,3)
		plot(rt.x,abs(rt.Q(1)));
		legend('|Q_1|');

		subplot(2,3,6)
		plot(rt.x,angle(rt.Q(1)));
		legend('arg(Q_1)');
		ylim(pi*[-1,1]);

		figure(100+tid);
		clf();
		subplot(2,2,1);
		plot(rt.x,abs(zlr));
		legend('z_1^-','z_1^+');
		subplot(2,2,2);
		plot(rt.x,abs(Qlr));
		legend('Q_1^-','Q_1^+');
		subplot(2,2,3);
		plot(rt.x,angle(zlr));
		ylim(pi*[-1,1]);
		subplot(2,2,4);
		plot(rt.x,angle(Qlr));
		ylim(pi*[-1,1]);

	end % if pflag
end % river_tide_test_07

