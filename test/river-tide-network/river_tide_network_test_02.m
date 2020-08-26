% Tue 18 Aug 09:27:08 +08 2020
function [out] = river_tide_network_test_02(rt_map,pflag)
	tid  = 2;
	name = 'single channel, split in two'

z1flag = true;
%z1flag = false;

	z10 = [1];
%	L   = [5e3,

	%for idx=1:length(z10)

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
	S0        = -normal_flow_slope(Q0,h0,w0,drag2chezy(cD));
	zbfun     = @(x) -h0 + S0*x;

	% base frequency
	T         = Constant.SECONDS_PER_DAY;
	omega     = 2*pi/T;
	% domain size
	Xi        = [0,0.5;
                     0.5,1]*h0/S0;

	meta = river_tide_test_metadata();
	opt = meta.opt;
	optmaxiter =  100;
	opt.nx = 25;
	opt.sopt.maxiter = 100;

for idx=1:2
	bc        = struct();
	% mean sea level
	bc(1,1).var = 'z';
	bc(1,1).rhs = 0;
	% Dirichlet condition
	bc(1,1).p   = 1;

	% river discharge
	if (1 == idx)
		bc(2,1).var = '';
		bc(2,1).rhs = NaN;
	else
		bc(2,1).var = 'Q';
		bc(2,1).rhs = Q0;
	end
	%bc(2,1).var = '';
%	if (1==idx)
%		% upstream condition set by junction condition
%		bc(2,1).var = '';
%		% no discharge for downstream channel
%		bc(1,2).var = '';
%	else
%		% no upstream wl required
%		bc(2,1).var = '';
%
%		bc(1,2).var = 'Q';
%		bc(1,2).rhs = Q0;
%		% Dirichlet condition
%		bc(1,2).p   = 1;
%	end
%	% discharge has only one condition
%	bc(2,2).var = '';
	k = 2;
	% wave entering from left
	bc(1,k).var = 'z';
	%z10         = 1; %sqrt(eps);
	bc(1,k).rhs = z10;
	bc(1,k).p   = [1,0];
	bc(1,k).q   = [1,0];
	% wave entering from right / reflected wave
	bc(2,k).var = 'z';
	bc(2,k).rhs =   0;
	bc(2,k).p   = [1,0];
	bc(2,k).q   = [0,1];

	out(idx) = River_Tide( ...
				   'fun.zb',      zbfun ...
				 , 'fun.cd',      cdfun ...
				 , 'fun.width',   wfun ...
				 , 'omega',       omega ...
				 , 'opt',         opt ...
				 , 'Xi',          Xi(idx,:) ...
				);
if (~z1flag)
	out(idx).opt.oflag = false(1,4);
end
	out(idx).bc = bc;
end

function [cid,eid,p] = jfun()
	cid = [1,2];
	eid = [2,1];
	p   = [1./w0,1./w0];
end

	rtn = River_Tide_Network_2(out);
	rtn.junction_condition = {@jfun}

	rtn.init();

	rtn.solve();

	figure(1)
	clf
	for idx=1:length(rtn.rt)
		x = rtn.rt(idx).x;

		%subplot(length(rtn.rt),4,4*(idx-1)+1);
		subplot(1,4,1)
		z0 = rtn.rt(idx).z(0);
		plot(x,z0);
		hold on
		
		%subplot(length(rtn.rt),4,4*(idx-1)+2);
		subplot(1,4,2)
		Q0 = rtn.rt(idx).Q(0);
		plot(x,Q0);
		hold on
if (z1flag)
		k = 1;
		subplot(1,4,3)
		%subplot(length(rtn.rt),4,4*(idx-1)+3);
		z = rtn.rt(idx).z(k);
		plot(x,abs(z));
		hold on
		
		subplot(1,4,4)
		%subplot(length(rtn.rt),4,4*(idx-1)+4);
		Q = rtn.rt(idx).Q(k);
		%plot(x,abs(Q));
		plot(x,real(Q));
		%abs(Q));
		set(gca,'colororderindex',idx)
		plot(x,imag(Q),'--')
		hold on
end
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

