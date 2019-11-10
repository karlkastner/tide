% Wed  9 Oct 15:23:10 PST 2019
function [fail,rmse,name,rt] = river_tide_test_11(rt_map,pflag)
	reltol = 1e-2;
	tid  = 12;
	name = 'swap left and right end';

	T   = Constant.SECONDS_PER_DAY;

	% river discharge
	Q0  = -10;
	z10 = [1,sqrt(eps)];
	z1L = [0,1];

	L = 5e5;
	for idx=1:2
	Qu = -Q0;
	% width of channel
	w0        = 1;
	wfun      = @(x)   w0*ones(size(x));
	% drag/friction coefficient
	cD        = 2.5e-3;
	cdfun     = @(x)  cD*ones(size(x));
	% bed level of channel
	h0        = 10;
	S0         = normal_flow_slope(-Qu,w0,h0,drag2chezy(cD))
%	figure(1)
	if (1==idx)
	clf
	zbfun     = @(x) -h0 + S0*x;
%	plot(zbfun(linspace(0,L)))
%	hold on
	else
		zbfun     = @(x) -h0 + S0*(L-x);
%	plot(zbfun(linspace(0,L)))
	end
	bc        = struct();

	% mean sea level
	if (1 == idx)
	bc(1,1).var = 'z';
	bc(1,1).rhs = 0;
	% river discharge
	bc(2,1).var = 'Q';
	bc(2,1).rhs = -Q0;
	else
	bc(1,1).var = 'Q';
	bc(1,1).rhs = Q0;
	% river discharge
	bc(2,1).var = 'z';
	bc(2,1).rhs = 0;
	end
	% Dirichlet condition
	bc(1,1).p   = 1;
	% Dirichlet condition
	bc(2,1).p   = 1;
q = 0;
	% wave entering from left
	bc(1,2).var = 'z';
	bc(1,2).rhs = z10(idx);
	bc(1,2).p   = [1,0];
	bc(1,2).q   = [1,q];

	bc(1,3).var = 'z';
	bc(1,3).rhs = 0;
	bc(1,3).p   = [1,0];
	bc(1,3).q   = [1,0];

	% wave entering from right / reflected wave
	bc(2,2).var = 'z';
	bc(2,2).rhs = z1L(idx);
	bc(2,2).p   = [1,0];
	bc(2,2).q   = [q,1];
	bc(2,3).var = 'z';
	bc(2,3).rhs = 0;
	bc(2,3).p   = [1,0];
	bc(2,3).q   = [0,1];

	% base frequency
	omega     = 2*pi/T;
	% domain size
	Xi        = [0,L];
	% model for river tide
	opt.model_str = 'wave';
	% solver of boundary value problem
	opt.solver = @bvp2c;
	%opt.solver = @bvp2fdm;
	% number of points along channel
	opt.nx     = 50; % 200
	% change of distance between points along channel 
	opt.xs     = 1; 
	opt.o2     = true;

	% solve with model
	[rt(idx)]  = rt_map.fun({Xi} ... % Q0,
			, {wfun}, {cdfun}, {zbfun}, omega ...
			, {bc(1,1).var}, {bc(1,1).rhs}, {bc(1,1).p} ...
			, {bc(2,1).var}, {bc(2,1).rhs}, {bc(2,1).p} ...
			, {bc(1,2).var}, {bc(1,2).rhs}, {bc(1,2).p}, {bc(1,2).q} ...
			, {bc(2,2).var}, {bc(2,2).rhs}, {bc(2,2).p}, {bc(2,2).q} ...
			, {bc(1,3).var}, {bc(1,3).rhs}, {bc(1,3).p}, {bc(1,3).q} ...
			, {bc(2,3).var}, {bc(2,3).rhs}, {bc(2,3).p}, {bc(2,3).q} ...
			, opt);

	end

	z = ([rt(1).z(1),flipud(rt(2).z(1))]);
	res = z(:,1)-z(:,2);
	rmse = rms(res);
	fail = (rmse > reltol*rms(z(:,1)));

	if (pflag)
		figure(1);
		clf();
		subplot(2,2,1)
		plot(rt(1).x,abs(z(:,1)));
		hold on
		plot(rt(1).x,abs(z(:,2)),'--');
		ylabel('|z|');
	
		subplot(2,2,2)
		plot(rt(1).x,angle(z(:,1)));
		hold on
		plot(rt(1).x,angle(z(:,2)),'--');
		ylabel('angle(z)');
	
		subplot(2,2,3)
		plot(rt(1).x,[rt(1).zb,rt(2).zb],'k');
		hold on
		plot(rt(1).x,[rt(1).z(0),rt(2).z(0)],'r--');
	end

	% compare to analytical solution
%	g = Constant.gravity;
%	c = sqrt(g*h0);
%	k = omega/c;
%	x = rt.x;

%	bw = Backwater1D();
%	nn = opt.nx;
%	[x_, h_, z0_] = bw.solve(-Q0,0,drag2chezy(cD),wfun,zbfun,0,Xi);
%	[x_, h_, z0_] = bw.solve_analytic(-Q0,drag2chezy(cD),w0,S0,h0,nn);
%	z0  = interp1(x_,z0_,rt.x,'linear','extrap');
%	z0 = S0*x;
%	z0t = rt.z_(:,1) - z0;
%	z2 = rt.z_(:,3);
%
%	z2_ = rt.even_overtide_analytic(z10);
%	%z0_ = interp1(x_,h_,rt.x,'spline')+0*zbfun(rt.x);
%	% r = (1+1i)*sqrt(-cD.*omega.*Q0/w0./(g*h0.^3));
%	% z = z10*exp(-r*x);
%
%	rmse(2)  = norm(z2_-z2);
%
%	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
%	nres_ = rms(cdiff(z2_,2))
%	fail = (rmse(1) > 0.05/z10) || (rmse(2) > 10*nres_);
%
%	zz = NaN(size(rt.z_));
%	zz(:,3) = z2_;
%	river_tide_test_plot(tid,rt,zz,name,pflag);
%	
end % river_tide_test_09

