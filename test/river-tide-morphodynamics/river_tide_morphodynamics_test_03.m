% Tue 18 Aug 09:27:08 +08 2020
function [t,zb,rtn] = river_tide_morphodynamics_test_03(rt_map,pflag)
	tid  = 2;
	name = 'two upstream channels, one downstream channel';

	% river discharge
	Q0   = -10;
	Qlim = [-2,-20];
	iorder = 1;

	z10  =  1;

	% slope
	S0   = 2.5e-5;
	Cd   = 2.5e-3;

	w0   = [1,1/3,2/3];
	h0   = normal_flow_depth(Q0,w0(1),Cd,S0,'cd');

	L      = h0/S0
	z1     = 1;
	pz1    = 0;
	omega  = 2*pi/Physics.SECONDS_PER_DAY;
	pL     = 0.25;
	ruleQ  = NaN; % not yet used
	ruleQs = 'geometric';
	nx0    = 50;
	Xi     = 1.2*L;
	dx     = L*1/(nx0-1);
%	dx     = L/49;
	Ti     = Physics.SECONDS_PER_YEAR*2000;
	cfl    = 0.99;
	scheme = 'upwind';
	d_mm   = 0.2;

	rtn = River_Tide_Network_2();

	opt.xs     = 1; 

	map = River_Tide_Network_Map('mat/rt-mor-map.mat');
%	map.recompute = true;
	rtn = map.fun(   {z1} ...
		 , {pz1} ...	% [1]     reflected wave factor
		 , {omega} ...	% [rad/s] anguluar frequency of tide
		 , {Qlim(1)} ...	% [m^3/s] discharge at inflow bc
		 , {Qlim(2)} ...	% [m^3/s] discharge at inflow bc
		 , {iorder} ...
		 , {S0} ...	% [1]	  upstream slope
		 , {d_mm} ...	% 
		 , {Xi} ...	% [m]     doamin length
		 , {pL} ...	% [1]     relative location of bifurcation as seen
		     ...	%         from mouth with respect to domain length
		 , {w0(1)} ...	% [m]     width of upstream channel
		 , {w0(2)} ...	% [m]     width of upstream channel
		 , {w0(3)} ...	% [m]     width of upstream channel
		 , {Cd} ...
		 , {ruleQ} ...	%  -      rule for discharge division
		 , {ruleQs} ...	%  -      rule for sediment division
		 , {dx} ...	% [m]     spatial discretization step
		 , {Ti} ...	% [s]     simulated time span
		 , {cfl} ...	% [1]     cfl condition
		 , {scheme} ... 	%  -      numerical scheme
		 ... %, {opt} ...
	);
	map.save();

	t  = rtn.evolution.t;
	zb = rtn.evolution.zb;

	figure(20);
	x = [];
	for idx=1:length(rtn.rt)
		x = [x;mid(rtn.rt(idx).x)];
		plot(rtn.rt(idx).x/L, rtn.rt(idx).zb(rtn.rt(idx).x));
		hold on;
	end
%	set(gca,'colororderindex',1);
%	plot(x,zb(:,1),'.-');
%	plot(x,zb(:,1),'.-');

	figure(21);
	for idx=1:length(rtn.rt)
		%x = [x;mid(rtn.rt(idx).x)];
		plot(rtn.rt(idx).x/L, rtn.rt(idx).h0);
%(rtn.rt(idx).x));
		hold on
	end

	figure(3);
        clf();
	subplot(2,2,1);
	d(:,1) = max(abs(zb(:,1:end-1)-zb(:,end)));
	d(:,2) = max(abs(diff(zb,[],2)));
	plot(t(1:end-1),d);
	subplot(2,2,2);
	for idx=1:length(rtn.rt)
		Qs = rtn.rt(idx).sediment_transport(0,1);
		plot(mid(rtn.rt(idx).x),Qs(2:end-1));
		hold on;
	end

	figure(1)
	clf();
	e = 1e-2;
	for idx=1:length(rtn.rt)
		figure(1);
		x = rtn.rt(idx).x;

		%subplot(length(rtn.rt),4,4*(idx-1)+1);
		subplot(1,4,1);
		z0 = rtn.rt(idx).z(0);
		plot(x,z0);
		hold on;
		
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

