% Wed  9 Oct 15:55:23 PST 2019
function [fail,rmse,name,rt] = river_tide_test_10(rt_map,pflag)
	tid = 10;
	name = 'ideal estuary, width convergence cancel frictional damping';
	% river discharge
	Q0        = 0;
	% width of channel
	w0        = 1e5;
	% drag/friction coefficient
	cD        = 2.5e-3;
	cdfun     = @(x)  cD*ones(size(x));
	% bed level of channel
	h0        = 10;
	zbfun     = @(x) -h0*ones(size(x));
	bc        = struct();
	% base frequency
	T         = Constant.SECONDS_PER_DAY;
	omega     = 2*pi/T;

	z10       = 0.1;
	g = Constant.gravity;
	c = sqrt(g*h0);
	k0 = omega/c;

	% damping by tide
	k = wave_number_tide(omega,cD,h0,abs(z10));
	% damping modulus
	r =  -imag(k);

	% convergence length of widh
	Lw = 0.5/r;
	%Lw = Inf;
	wfun      = @(x) w0*exp(-x/Lw);

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
	bc(1,2).rhs = z10;
	bc(1,2).p   = [1,0];
	bc(1,2).q   = [1,0];
	% wave entering from right / reflected wave
	bc(2,2).var = 'z';
	bc(2,2).rhs =   0;
	bc(2,2).p   = [1,0];
	bc(2,2).q   = [0,1];
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

	% compate to analytical solution
	x = rt.x;
	z = z10*exp(-1i*k*x + 0.5./Lw*x);

	rmse(2)  = rms(rt.z_(:,2)-z);
	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
	nres_ = rms(cdiff(rt.z_(:,2),1));
	fail = (rmse(1) > 0.05/z10) || (rmse(2) > 10*nres_);

	zz = NaN(size(rt.z_));
	zz(:,2) = z;
	river_tide_test_plot(tid,rt,zz,name,pflag);
	
	end % river_tide_test_10

