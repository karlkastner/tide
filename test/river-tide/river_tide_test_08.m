% Wed  9 Oct 15:23:10 PST 2019
function [fail,rmse,name,rt] = river_tide_test_08(rt_map,pflag)
	tid  = 8;
	s = 10;
	name = 'subtidal water level offset, uniform flow';
	% river discharge
	Qu        = -2*s;
	Q0        =  Qu;
	% width of channel
	w0        = 1;
	wfun      = @(x)   w0*ones(size(x));
	% drag/friction coefficient
	cD        = 2.5e-3;
	cdfun     = @(x)  cD*ones(size(x));
	% bed level of channel
	h0        = 1*s;
	S0         = normal_flow_slope(-Qu,w0,h0,drag2chezy(cD));
	zbfun     = @(x) -h0 + S0*x;
	bc        = struct();
	% mean sea level
	bc(1,1).var = 'z';
	bc(1,1).rhs = 0;
	% Dirichlet condition
	bc(1,1).p   = 1;
	% river discharge
	bc(2,1).var = 'Q';
	bc(2,1).rhs = -Q0;
	% Dirichlet condition
	bc(2,1).p   = 1;
	% wave entering from left
	bc(1,2).var = 'z';
	z10       = 0.01*s;
	bc(1,2).rhs = z10;
	bc(1,2).p   = [1,0];
	bc(1,2).q   = [1,0];
	% wave entering from right / reflected wave
	bc(2,2).var = 'z';
	bc(2,2).rhs =   0;
	bc(2,2).p   = [1,0];
	bc(2,2).q   = [0,1];
	% base frequency
	T         = Constant.SECONDS_PER_DAY;
	omega     = 2*pi/T;
	% domain size
	Xi        = [0,1e5*sqrt(s)];

	meta = river_tide_test_metadata();
	opt = meta.opt;
	opt.sopt.maxiter = 100;
	opt.sopt.reltol = 1e-8;
%	opt.hmode = 'iterate';

	% solve with model
	[rt]  = rt_map.fun({Xi} ... % Q0,
			, {wfun}, {cdfun}, {zbfun}, omega ...
			, {bc(1,1).var}, {bc(1,1).rhs}, {bc(1,1).p} ...
			, {bc(2,1).var}, {bc(2,1).rhs}, {bc(2,1).p} ...
			, {bc(1,2).var}, {bc(1,2).rhs}, {bc(1,2).p}, {bc(1,2).q} ...
			, {bc(2,2).var}, {bc(2,2).rhs}, {bc(2,2).p}, {bc(2,2).q} ...
			, {}, {}, {}, {} ...
			, {}, {}, {}, {} ...
			, opt);

	% compare to analytical solution
	g = Constant.gravity;
	c = sqrt(g*h0);
	k = omega/c;
	x = rt.x;

%	bw = Backwater1D();
%	nn = opt.nx;
%	[x_, h_, z0_] = bw.solve(-Q0,0,drag2chezy(cD),wfun,zbfun,0,Xi);
%	[x_, h_, z0_] = bw.solve_analytic(-Q0,drag2chezy(cD),w0,S0,h0,nn);
%	z0  = interp1(x_,z0_,rt.x,'linear','extrap');
	z0 = S0*x;
	z0t = rt.z_(:,1) - z0;

	z0t_ = rt.mwl_offset_analytic(z10);
	%z0_ = interp1(x_,h_,rt.x,'spline')+0*zbfun(rt.x);
	% r = (1+1i)*sqrt(-cD.*omega.*Q0/w0./(g*h0.^3));
	% z = z10*exp(-r*x);

	rmse(2)  = rms(z0t_-z0t);

	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
	nres_ = rms(cdiff(z0t_,2))
	fail = (rmse(1) > 0.05/z10) || (rmse(2) > 10*nres_);

	
	if (pflag)
		namedfigure(tid,['Test: ',name]);
		clf();

		subplot(2,2,1);
		%[zbfun(x),rt.z_(:,1)]);
		plot(rt.x,z0t);
		hold on;
		plot(rt.x,z0t_,'--');
		%plot(rt.x,z0t_./z0t,'-.');
		legend('z_0');

%		namedfigure(tid,['Test: ',name]);
%		clf();
%		subplot(2,3,1);
%		plot(rt.x,[zbfun(x),rt.z_(:,1)]);
%		hold on;
%		plot(x_,z0_,'--');
%		legend('z_b','z_0');
%
%		subplot(2,3,4);
%		plot(rt.x,[wfun(x)]);
%		legend('z_b','z_0');
%
%		subplot(2,3,2);
%		plot(rt.x,abs(rt.z_(:,2)));
%%		hold on;
%%		plot(x,abs(z),'--');
%		legend('|z_1|');
%
%		subplot(2,3,5)
%		plot(rt.x,angle(rt.z_(:,2)));
%%		hold on;
%%		plot(x,angle(z),'--');
%		legend('arg(z_1)');
%		ylim(pi*[-1,1]);
%
%		subplot(2,3,3)
%		plot(rt.x,abs(rt.Q_(:,2)));
%		legend('|Q_1|');
%
%		subplot(2,3,6)
%		plot(rt.x,angle(rt.Q_(:,2)));
%		legend('arg(Q_1)');
%
%		%dQ_dx = derivative1(rt.x,rt.Q_(:,2));
%		%plot(rt.x,[real(dQ_dx),imag(dQ_dx)]);
%		%rt.Q_(:,2)),imag(rt.Q_(:,2))]);
%		%subplot(2,2,3)
%		%plot(diff(rt.x));
	end % if pflag
end % river_tide_test_1

