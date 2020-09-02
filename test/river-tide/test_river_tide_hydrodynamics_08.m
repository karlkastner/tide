% Wed  9 Oct 15:23:10 PST 2019
function [out,rt] = river_tide_test_08(rt_map,pflag)
	if (nargin()<2)
		pflag = 1;
	end
	out.id  = 8;
	out.name = 'subtidal water level offset, uniform flow';

	% mean surface elevation
	z00 = 0;
	% tidal surface elevation
	z10 = 0.1;
	zs = [0,z10];
	% river discharge
	Q0        = -10;
	% width of channel
	w0        = 1;
	% drag/friction coefficient
	Cd        = 2.5e-3;
	% depth of channel without mwl-offset
	h0        = 10;
	% slope of channel bed
	S0        = -normal_flow_slope(Q0,h0,w0,drag2chezy(Cd));
	% bed level of channel
	zb        = @(x) -h0 + S0*x;
	% base frequency
	T         = Constant.SECONDS_PER_DAY;
	omega     = 2*pi/T;

	% length of computational domain
	Lx         = 1e6;

	meta = river_tide_test_metadata();
	opt = meta.opt;
	opt.sopt.maxiter = 100;
	opt.sopt.reltol = 1e-8;

	% reflection coefficient at right end of boundary
	qr = 0;
	rt = hydrodynamic_scenario(rt_map,zs,qr,zb,Q0,w0,Cd,omega,Lx,opt);

	Xi = rt.hydrosolver.xi;

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
	z0 = S0*x;
	z0t = rt.z(0) - z0;
%mwl_offset_analytic(obj,x,z10,h0,w0,Cd,Q0) 
	z0t_ = rt.mwl_offset_analytic(x,z10,h0,w0,Cd,abs(Q0));
	%z0_ = interp1(x_,h_,rt.x,'spline')+0*zbfun(rt.x);
	% r = (1+1i)*sqrt(-Cd.*omega.*Q0/w0./(g*h0.^3));
	% z = z10*exp(-r*x);

	rmse(2)  = rms(z0t_-z0t);

	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
	nres_ = rms(cdiff(z0t_,2))
	result = (rmse(1) > 0.05/z10) || (rmse(2) > 10*nres_);

	out.rmse   = rmse;
	out.result = result;
	
	if (pflag)
		namedfigure(out.id,['Test: ',out.name]);
		clf();

		subplot(2,2,1);
		%[zbfun(x),rt.z(0)]);
		plot(rt.x,z0t);
		hold on;
		plot(rt.x,z0t_,'--');
		%plot(rt.x,z0t_./z0t,'-.');
		legend('z_0');

%		namedfigure(out.id,['Test: ',out.name]);
%		clf();
%		subplot(2,3,1);
%		plot(rt.x,[zbfun(x),rt.z(0)]);
%		hold on;
%		plot(x_,z0_,'--');
%		legend('z_b','z_0');
%
%		subplot(2,3,4);
%		plot(rt.x,[wfun(x)]);
%		legend('z_b','z_0');
%
%		subplot(2,3,2);
%		plot(rt.x,abs(rt.z(1))));
%%		hold on;
%%		plot(x,abs(z),'--');
%		legend('|z_1|');
%
%		subplot(2,3,5)
%		plot(rt.x,angle(rt.z(1))));
%%		hold on;
%%		plot(x,angle(z),'--');
%		legend('arg(z_1)');
%		ylim(pi*[-1,1]);
%
%		subplot(2,3,3)
%		plot(rt.x,abs(rt.Q(1)));
%		legend('|Q_1|');
%
%		subplot(2,3,6)
%		plot(rt.x,angle(rt.Q(1)));
%		legend('arg(Q_1)');
%
%		%dQ_dx = derivative1(rt.x,rt.Q(1));
%		%plot(rt.x,[real(dQ_dx),imag(dQ_dx)]);
%		%rt.Q(1)),imag(rt.Q(1))]);
%		%subplot(2,2,3)
%		%plot(diff(rt.x));
	end % if pflag
end % river_tide_test_1

