% Wed  9 Oct 15:23:10 PST 2019
function [fail,rmse,name,rt] = river_tide_test_06(rt_map,pflag)
	tid  = 6;
	name = 'infinitessimal wave along river with uniform flow';
	% river discharge
	Q0        = -100;
	% width of channel
	w0 = 10;
	wfun      = @(x)   w0*ones(size(x));
	% drag/friction coefficient
	cD        = 2.5e-3;
	cdfun     = @(x)  cD*ones(size(x));
	% bed level of channel
	h0        = 10;
	S0         = normal_flow_slope(-Q0,h0,w0,drag2chezy(cD));
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
	z10         = 1e-2; %sqrt(eps);
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
	Xi        = [0,1e6];

	meta = river_tide_test_metadata();
	opt = meta.opt;

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

	% check ode
	rmse = rt.check_continuity();
	% compare to analytical solution
	g = Constant.gravity;
	c0 = sqrt(g*h0);
	k0 = omega/c0;
	x = rt.x;

	r = (1+1i)*sqrt(-cD.*omega.*Q0/w0./(g*h0.^3));
	z = z10*exp(-r*x);

	rmse(2)  = rms(rt.z_(:,2)-z)
	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
	nres_ = rms(cdiff(z,2));
	fail = (rmse(1)/z10 > 0.05) || (rmse(2) > 10*nres_);

	zz      = NaN(size(rt.z_));
	zz(:,2) = z;
	river_tide_test_plot(tid,rt,zz,name,pflag);
end % river_tide_test_06

