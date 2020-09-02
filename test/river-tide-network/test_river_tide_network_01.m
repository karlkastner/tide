% Tue 18 Aug 09:27:08 +08 2020
function [out] = test_river_tide_network_01(rt_map,pflag)
	if (nargin()<2)
		pflag = 1;
	end
	out.id   = 1;
	out.name = 'two independent (uncoupled) channels in one model';

	% tidal elevation
	z10 = [1,2];

	% river discharge
	Q0  = [-10,-5];

	% width of channel
	w0        = [1,1.5];
%	wfun      = @(cdx,x)   w0(cdx)*ones(size(x));

	% depth of channel
	h0        = [10,5];

	% drag/friction coefficient
	Cd        = 2.5e-3*[1,0.7];

	% slope of bed
	S0  = -normal_flow_slope(Q0,h0,w0,drag2chezy(Cd))

	% bed level of channel
	zb     = @(cdx,x) -h0(cdx) + S0(cdx)*x;


	% base frequency
	T         = Constant.SECONDS_PER_DAY;
	omega     = 2*pi/T;

	bc        = struct();

	% domain size
	Xi        = [0,2.0*h0(1)/S0(2);
                     0,1.5*h0(2)/S0(2)];

	meta             = river_tide_test_metadata();
	opt              = meta.opt;
	nx               = [100,200];
	opt.sopt.maxiter = 200;

	for idx=1:length(z10)

	% mean sea level
	bc(1,1,idx).var = 'z';
	bc(1,1,idx).rhs = 0;
	% Dirichlet condition
	bc(1,1,idx).p   = 1;
	% river discharge
	bc(2,1,idx).var = 'Q';
	bc(2,1,idx).rhs = Q0(idx);
	% Dirichlet condition
	bc(2,1,idx).p   = 1;
	% wave entering from left
	bc(1,2,idx).var = 'z';
	bc(1,2,idx).rhs = z10(idx);
	bc(1,2,idx).p   = [1,0];
	bc(1,2,idx).q   = [1,0];
	% wave entering from right / reflected wave
	bc(2,2,idx).var = 'z';
	bc(2,2,idx).rhs =   0;
	bc(2,2,idx).p   = [1,0];
	bc(2,2,idx).q   = [0,1];

	% solve individually
	hydrosolver    = BVPS_Characteristic();
	hydrosolver.xi = Xi(idx,:);
	hydrosolver.nx = nx(idx);
	opt.xs = opt.xs;
	opt.dischargeisvariable = true;
	
	zb_ = @(x) -h0(idx) + S0(idx)*x;

	% solve individually
	rt_a(idx) = River_Tide_BVP( ...
		   'zb',          zb_ ...
		 , 'cd',          Cd(idx) ...
		 , 'width',       w0(idx) ...
		 , 'omega',       omega ...
		 , 'opt',         opt ...
		 , 'hydrosolver', hydrosolver ...
		);
	rt_a(idx).bc = bc(:,:,idx);
	rt_a(idx).init();
	rt_a(idx).solve();
	end % for idx

	% solve as a "network"

	hydrosolver    = BVPS_Characteristic();
	hydrosolver.xi = Xi;
	hydrosolver.nx = nx;
	opt.xs = opt.xs;
	opt.dischargeisvariable = true;

	rtn = River_Tide_BVP( ...
		   'zb',          zb ...
		 , 'cd',          Cd ...
		 , 'width',       w0 ...
		 , 'omega',       omega ...
		 , 'opt',         opt ...
		 , 'hydrosolver', hydrosolver ...
		);
	rtn.bc = bc;
	rtn.init();
	rtn.solve();

	nc = rtn.nc;
	figure(1);
	clf();
	for idx=1:nc
		%rt = rtn.rt(idx)';
		x  = rtn.x(idx);

		subplot(nc,4,4*(idx-1)+1);
		z0 = [rt_a(idx).z(0),rtn.z(0,idx)];
%		zb = rt.zb(x);
		rmse(idx,1) = max(abs(diff(z0,[],2)));
		plot(x,z0); %[zb,z0]);
		ylabel('z0');
		
		subplot(nc,4,4*(idx-1)+2);
		%Q0 = rt.Q(0);
		Q0 = [rt_a(idx).Q(0),rtn.Q(0,idx)];
		rmse(idx,2) = max(abs(diff(Q0,[],2)));
		plot(x,Q0);
		ylabel('Q0');

		subplot(nc,4,4*(idx-1)+3);
		%Q0 = rt.Q(0);
		z1 = [rt_a(idx).z(1),rtn.z(1,idx)];
		rmse(idx,2) = max(abs(diff(z1,[],2)));
		plot(x,abs(z1));
		ylabel('z1');

%		k = 1;
%		subplot(nc,4,4*(idx-1)+3);
%		z = rt.z(k);
%		[Qlr,zlr]=rt.decompose();
%		plot(x,[abs(z),abs(zlr)])
%		%plot(x,[abs(z),real(z),imag(z)]);
%		title('z1');	
	
%		subplot(nc,4,4*(idx-1)+4);
%		Q = rt.Q(k);
%		plot(x,abs(Q));
%		title('Q1');	

	end % for idx

	out.rmse   = rmse
	out.result = any(flat(rmse>1e-2))

%	% solve with model
%	[rt]  = rt_map.fun({Xi} ... % Q0,
%			, {wfun}, {Cdfun}, {zbfun}, omega ...
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
%	r = (1+1i)*sqrt(-Cd.*omega.*Q0/w0./(g*h0.^3));
%	z = z10*exp(-r*x);
%
%	rmse(2)  = rms(rt.z_(:,2)-z)
%	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
%	nres_ = rms(Cdiff(z,2));
%	fail = (rmse(1)/z10 > 0.05) || (rmse(2) > 10*nres_);
%
%	zz      = NaN(size(rt.z_));
%	zz(:,2) = z;
%	river_tide_test_plot(tid,rt,zz,name,pflag);
end % river_tide_test_06

