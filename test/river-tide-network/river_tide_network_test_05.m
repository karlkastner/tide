% Tue 18 Aug 09:27:08 +08 2020
function [out] = river_tide_network_test_05(rt_map,pflag)
	tid  = 2;
	name = 'two upstream channels, one downstream channel';

	% river discharge
	Q0   = -1;
	Q0_  = -10;
	z10 = 1;

	h0  = 10;
	w0  = [1,1/3,2/3];
	%w0  = [1,1/2,1/2];
	nx  = 20;
	Cd = 2.5e-3;
	Cd_ = 2.5e-3;

	S0  = -normal_flow_slope(Q0,h0,w0(1),drag2chezy(Cd));
	S0_  = -normal_flow_slope(Q0_,h0,w0(1),drag2chezy(Cd_));

	L = h0/S0_;

	% domain size
	Xi        = L*[0.5,1;
		             0,0.5;
		             0,0.5];

	meta = river_tide_test_metadata();
	opt=meta.opt;
	opt.maxiter =  100;
	opt.nx = nx;
	opt.sopt.maxiter = 100;
	
	for idx=1:length(w0)

	% width of channel
	wfun      = @(x)   w0(idx)*ones(size(x));

	% drag/friction coefficient
	cdfun     = @(x)  Cd*ones(size(x));

	% bed level of channel
	zbfun     = @(x) -h0 + S0*x;

	% base frequency
	T         = Constant.SECONDS_PER_DAY;
	omega     = 2*pi/T;


	bc        = struct();
	% mean sea level
	if (1~=idx)
		bc(1,1).var = 'z';
		bc(1,1).rhs = 0;
		% Dirichlet condition
		bc(1,1).p   = 1;
	else
		bc(1,1).var = '';
		bc(1,1).rhs = [];
	end

	% river discharge
	if (1==idx)
		bc(2,1).var = 'Q';
		bc(2,1).rhs = Q0;
	else
		bc(2,1).var = '';
		bc(2,1).rhs = [];
	end

	% wave entering from left
	if (1~=idx)
		bc(1,2).var = 'z';
		bc(1,2).rhs = z10;
		bc(1,2).p   = [1,0];
		bc(1,2).q   = [1,0];
	else
		bc(1,2).var = '';
		bc(1,2).rhs = [];
	end

	% wave entering from right
	if (1 == idx)
		bc(2,2).var = 'z';
		bc(2,2).rhs =   0;
		bc(2,2).p   = [1,0];
		bc(2,2).q   = [0,1];
	else
		bc(2,2).var = '';
		bc(2,2).rhs = [];
	end
	out(idx) = River_Tide( ...
				   'fun.zb',      zbfun ...
				 , 'fun.cd',      cdfun ...
				 , 'fun.width',   wfun ...
				 , 'omega',       omega ...
				 , 'opt',         opt ...
				 , 'Xi',          Xi(idx,:) ...
				);
	out(idx).bc = bc;
end

function [cid,eid,p] = jfun()
	cid = [1,2,3];
	eid = [1,2,2];
	p   = [1./w0];
end

	rtn = River_Tide_Network_2(out);
	rtn.junction_condition = {@jfun}

	rtn.init();

	rtn.solve();

	figure(1)
	clf();
e = 1e-2;
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
		%plot(x,Q0);
		plot(x,Q0*(1+e*idx));
		hold on
	
		[Qlr,zlr] = rtn.rt(idx).decompose();

		k = 1;
		subplot(1,4,3)
		%subplot(length(rtn.rt),4,4*(idx-1)+3);
		z = rtn.rt(idx).z(k);
		%set(gca,'colororderindex',1);
		plot(x,abs(z)*(1 + e*idx));
		hold on
		%plot(x,abs(zlr),'--');
		
		subplot(1,4,4)
		%subplot(length(rtn.rt),4,4*(idx-1)+4);
		Q = rtn.rt(idx).Q(k);
		%set(gca,'colororderindex',1);

		plot(x,abs(Q));
		%plot(x,real(Q));
		hold on
		%set(gca,'colororderindex',idx)
		%plot(x,imag(Q),'--')
if (1==idx)
	Q_(1) = Q(1);
else
	Q_(idx) = Q(end);
end
		%plot(x,abs(Q)*(1+0.02*idx));
		%hold on
		%plot(x,[real(Q),imag(Q)]);
		%plot(x,abs(Q));
		%hold on
		%plot(x,angle(Q))
		%plot(x,abs(Qlr),'--');
		hold on
	end
% test
w = [-1,1,1];
Q_
w.*abs(Q_)
sum(w.*Q_)
abs(sum(w.*Q_))./mean(abs(Q))

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
%	r = (1+1i)*sqrt(-Cd.*omega.*Q0/w0./(g*h0.^3));
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

