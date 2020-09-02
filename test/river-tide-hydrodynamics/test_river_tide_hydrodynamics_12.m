% Wed  9 Oct 15:23:10 PST 2019
function [out,rt] = river_tide_test_11(rt_map,pflag)
	if (nargin()<2)
		pflag = 1;
	end
	reltol = 1e-2;
	out.id  = 12;
	out.name = 'swap left and right end';

	T   = Constant.SECONDS_PER_DAY;

	% river discharge
	Q0  = -10;
	z10 = [1,0];
	z1L = [0,1];

	L = 5e5;

	% width of channel
	w0        = 1;
	% drag/friction coefficient
	Cd        = 2.5e-3;
	% depth of channel
	h0        = 10;

	for idx=1:2
	if (2 == idx)
		Q0 = -Q0;
	end

	% slope of channel
	S0         = -normal_flow_slope(Q0,h0,w0,drag2chezy(Cd))

	% bed level of channel
	if (1==idx)
		zb  = @(x) -h0 + S0*x;
	else
		
		zb  = @(x) -h0 + S0*(x-L);
	end
	bc        = struct();

	% mean sea level
	if (1 == idx)
	bc(1,1).var = 'z';
	bc(1,1).rhs = 0;
	% river discharge
	bc(2,1).var = 'Q';
	bc(2,1).rhs = Q0;
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

	meta = river_tide_test_metadata();
	opt = meta.opt;
	opt.oflag   = [true(1,1)];

	% solve with model
	[rt(idx)]  = rt_map.fun({Xi} ... % Q0,
			, {w0}, {Cd}, {zb}, omega ...
			, {bc(1,1).var}, {bc(1,1).rhs}, {bc(1,1).p} ...
			, {bc(2,1).var}, {bc(2,1).rhs}, {bc(2,1).p} ...
			, {bc(1,2).var}, {bc(1,2).rhs}, {bc(1,2).p}, {bc(1,2).q} ...
			, {bc(2,2).var}, {bc(2,2).rhs}, {bc(2,2).p}, {bc(2,2).q} ...
			, {bc(1,3).var}, {bc(1,3).rhs}, {bc(1,3).p}, {bc(1,3).q} ...
			, {bc(2,3).var}, {bc(2,3).rhs}, {bc(2,3).p}, {bc(2,3).q} ...
			, opt);

	end


	zb = ([rt(1).zb,flipud(rt(2).zb)]);

	z0 = ([rt(1).z(0),flipud(rt(2).z(0))]);

	z1 = ([rt(1).z(1),flipud(rt(2).z(1))]);
	res = z1(:,1)-z1(:,2);
	rmse = rms(res);
	result = (rmse(1) > reltol*rms(z1(:,1)));

	% dummy value
	rmse(2) = NaN;
	out.rmse   = rmse;
	out.result = result;

	if (pflag)
		figure(1);
		clf();
		subplot(2,2,1)
		plot(rt(1).x,[z0(:,1) zb(:,1)]);
		hold on
		plot(rt(1).x,[z0(:,2) zb(:,2)],'--');
		ylabel('z_0');
	
		subplot(2,2,3)
		plot(rt(1).x,abs(z1(:,1)));
		hold on
		plot(rt(1).x,abs(z1(:,2)),'--');
		ylabel('|z_1|');
	
		subplot(2,2,4)
		plot(rt(1).x,angle(z1(:,1)));
		hold on
		plot(rt(1).x,angle(z1(:,2)),'--');
		ylabel('arg(z_1)');
	
		subplot(2,2,2)
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
%	[x_, h_, z0_] = bw.solve(-Q0,0,drag2chezy(Cd),wfun,zbfun,0,Xi);
%	[x_, h_, z0_] = bw.solve_analytic(-Q0,drag2chezy(Cd),w0,S0,h0,nn);
%	z0  = interp1(x_,z0_,rt.x,'linear','extrap');
%	z0 = S0*x;
%	z0t = rt.z(0) - z0;
%	z2 = rt.z_(:,3);
%
%	z2_ = rt.even_overtide_analytic(z10);
%	%z0_ = interp1(x_,h_,rt.x,'spline')+0*zbfun(rt.x);
%	% r = (1+1i)*sqrt(-Cd.*omega.*Q0/w0./(g*h0.^3));
%	% z = z10*exp(-r*x);
%
%	rmse(2)  = norm(z2_-z2);
%
%	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
%	nres_ = rms(cdiff(z2_,2))
%	fail = (rmse(1) > 0.05/z10) || (rmse(2) > 10*nres_);
%
%	zz = NaN(size(rt.out.z));
%	zz(:,3) = z2_;
%	river_tide_test_plot(out.id,rt,zz,out.name,pflag);
%	
end % river_tide_test_09

