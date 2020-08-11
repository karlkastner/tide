% Wed  9 Oct 15:23:10 PST 2019
function [t,zb,out] = river_tide_test_06(rt_map,pflag)
%function [t,zb,z1,fail,rmse,name,rt] = river_tide_test_06(rt_map,pflag)
	tid  = 6;
	name = 'infinitessimal wave along river with uniform flow';
	% river discharge
	Q0        = -10;
	% width of channel
	w0        = 1;
	wfun      = @(x)   w0*ones(size(x));
	% drag/friction coefficient
	cD        = 2.5e-3;
	cdfun     = @(x)  cD*ones(size(x));
	% initial bed level of channel
	h0        = 10;
	S0         = normal_flow_slope(-Q0,h0,w0,drag2chezy(cD));
	zbfun     = @(x) -h0 + S0*x + 0*1e-1*randn(size(x));
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
	z10         = 1;%sqrt(eps);
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
	Xi        = [0,12*h0/S0];


	meta       = river_tide_test_metadata();
	opt        = meta.opt;
	optmaxiter = 100;
	opt.nx     = 75;
	opt.sopt.maxiter = 100;
	opt

	% save struct, z1 Q0, L,nx, T,nt, uw/lf
	out = River_Tide_Morphodynamics( ...
				   'fun.zb',      zbfun ...
				 , 'fun.cd',      cdfun ...
				 , 'fun.width',   wfun ...
				 , 'omega',       omega ...
				 , 'opt',         opt ...
				 , 'Xi',          Xi ...
				);

	out.bc = bc;
	%out.set_bc_from_cell(bc0l_var,bc0l_val,bc0l_p,bc0r_var,bc0r_val,bc0r_p);

	secpyear =  86400*365.25;
	T        = [0,400*secpyear];
	nt       =  50;
	

	out.norder = 'upwind';
%	out.norder = 'leapfrog';
	[t,zb,z_] = out.evolve_bed_level(zbfun,T,nt);
	z1     = out.z(1);
	z0     = out.z(0);
	z0i = z_(:,1);
%z(0);

	figure(1)
	clf();
	subplot(2,2,1)
	plot([zb(:,1),zb(:,end)])
	ylabel('zb')
	
	subplot(2,2,2)
	plot([-(zb(:,end)-zb(:,1))])
	hold on;
	plot(abs(z1))
	hold on;
	plot([z0-z0i]);
	legend('zb(T)-zb(0)','|z_1|','z0(T)-z0(0)')
		

	subplot(2,2,3)
	d = [max(abs(diff(zb,[],2)))', max(abs(zb(:,1:end-1)-zb(:,end)))'];
	plot(d,'.-')
	legend('zb-zb(T-dT)','|zb-zb(T)|')

	subplot(2,2,4)
	plot([diff(z0)./diff(out.x)])
	hold on
	plot([diff(inner2outer(zb(:,end)))./diff(out.x)])
	hline(S0)
	ylabel('dzs/dx');		

	figure(2);
	subplot(2,1,1);
	[Qs,Qs0] = out.sediment_transport();
	plot([Qs,Qs0,Qs-Qs0]);
	abs(Qs([1:2 end-1:end])./Qs(end))
	Q = out.Q_;
	subplot(2,1,2)
	plot(abs(Q))
%
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

