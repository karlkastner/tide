% Tue 18 Aug 09:27:08 +08 2020
function [out] = river_tide_network_test_01(rt_map,pflag)
	tid  = 1;
	name = 'single channel'

	z10 = [1,2];
%,2];
%	L   = [5e3,

	for idx=1:length(z10)

	% river discharge
	Q0        = -10;
	% width of channel
	w0 = 1;
	wfun      = @(x)   w0*ones(size(x));
	% drag/friction coefficient
	cD        = 2.5e-3;
	cdfun     = @(x)  cD*ones(size(x));
	% bed level of channel
	h0        = 10;
	S0_ = -normal_flow_slope(-10,h0,w0,drag2chezy(cD))
	S0        = -normal_flow_slope(Q0,h0,w0,drag2chezy(cD));
	zbfun     = @(x) -h0 + S0*x;
	bc        = struct();
	% mean sea level
	bc(1,1).var = 'z';
	bc(1,1).rhs = 0;
	% Dirichlet condition
	bc(1,1).p   = 1;
	% river discharge
	bc(2,1).var = 'Q';
	bc(2,1).rhs = Q0;
	% Dirichlet condition
	bc(2,1).p   = 1;
	% wave entering from left
	bc(1,2).var = 'z';
	%z10         = 1; %sqrt(eps);
	bc(1,2).rhs = z10(idx);
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
	Xi        = [0,2*h0/S0_];


	meta = river_tide_test_metadata();
	opt = meta.opt;
%	opt.nx = 50
	opt.nx = 100
	opt.sopt.maxiter = 200;
	
	out(idx) = River_Tide( ...
				   'fun.zb',      zbfun ...
				 , 'fun.cd',      cdfun ...
				 , 'fun.width',   wfun ...
				 , 'omega',       omega ...
				 , 'opt',         opt ...
				 , 'Xi',          Xi ...
				);
	out(idx).opt.dischargeisvariable = true;
%	out(idx).opt.dischargeisvariable = false;
	out(idx).bc = bc;

	end
	rtn = River_Tide_Network_2(out);
	rtn.init();
	rtn.solve();
%	out(1).init;
%	out(1).solve;
%	rt = out(1);

	for idx=1:length(rtn.rt)
		rt = rtn.rt(idx)'
		x = rt.x;

		subplot(length(rtn.rt),4,4*(idx-1)+1);
		z0 = rt.z(0);
		zb = rt.zb(x);
		plot(x,[zb,z0]);
		ylabel('z0');
		
		subplot(length(rtn.rt),4,4*(idx-1)+2);
		Q0 = rt.Q(0);
		plot(x,Q0);
		ylabel('Q0')

		k = 1;
		subplot(length(rtn.rt),4,4*(idx-1)+3);
		z = rt.z(k);
		[Qlr,zlr]=rt.decompose();
		plot(x,[abs(z),abs(zlr)])
		%plot(x,[abs(z),real(z),imag(z)]);
		title('z1');	
	
		subplot(length(rtn.rt),4,4*(idx-1)+4);
		Q = rt.Q(k);
		plot(x,abs(Q));
		title('Q1');	

	end

%	% solve with model
%	[rt]  = rt_map.fun({Xi} ... % Q0,
%			, {wfun}, {cdfun}, {zbfun}, omega ...
%			, {bc(1,1).var}, {bc(1,1).rhs}, {bc(1,1).p} ...
%			, {bc(2,1).var}, {bc(2,1).rhs}, {bc(2,1).p} ...
%			, {bc(1,2).var}, {bc(1,2).rhs}, {bc(1,2).p}, {bc(1,2).q} ...
%			, {bc(2,2).var}, {bc(2,2).rhs}, {bc(2,2).p}, {bc(2,2).q} ...
%			, {}, {}, {}, {} ...
%			, {}, {}, {}, {} ...
%			, opt);
%
%	% check ode
%	rmse = rt.check_continuity();
%	% compare to analytical solution
%	g = Constant.gravity;
%	c0 = sqrt(g*h0);
%	k0 = omega/c0;
%	x = rt.x;
%
%	r = (1+1i)*sqrt(-cD.*omega.*Q0/w0./(g*h0.^3));
%	z = z10*exp(-r*x);
%
%	rmse(2)  = rms(rt.z_(:,2)-z)
%	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
%	nres_ = rms(cdiff(z,2));
%	fail = (rmse(1)/z10 > 0.05) || (rmse(2) > 10*nres_);
%
%	zz      = NaN(size(rt.z_));
%	zz(:,2) = z;
%	river_tide_test_plot(tid,rt,zz,name,pflag);
end % river_tide_test_06

