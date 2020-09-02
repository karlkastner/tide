% Wed  9 Oct 15:23:10 PST 2019
%function [t,zb,z1,fail,rmse,name,rt] = river_tide_test_06(rt_map,pflag)
function [out,t,ozb,rt] = test_river_tide_morphodyanmics_01(x0,zb0,pflag)
	if (nargin()<3)
		pflag = false;
	end
	dischargeisvariable = false;
	dischargeisvariable = true;

	out.id  = 1;
	out.name = 'infinitessimal wave along river with uniform flow';

	% tidal amplitude
	z10       = 1;

	% river discharge
	Q0        = -10;

	% width of channel
	w0        = 1;

	% drag/friction coefficient
	Cd        = 2.5e-3;

	% initial bed level of channel
	h0        = 10;
	h00 = h0;

	% inital slope of channel bed (asymmetotic slope at upstream boundary)
	S0         = -normal_flow_slope(Q0,h0,w0,drag2chezy(Cd));

	% inital bed level of channel 
	zb     = @(x) -h0 + S0*x;

	d_mm = 0.2;

	% base frequency
	T         = Constant.SECONDS_PER_DAY;
	omega     = 2*pi/T;

	bc        = struct();
	% mean sea level
	bc(1,1).var = 'z';
	bc(1,1).rhs = 0;
	% Dirichlet condition
	bc(1,1).p   = 1;
	% river discharge
	bc(2,1).rhs   = Q0;
	bc(2,1).var   = 'Q';
	% Dirichlet condition
	bc(2,1).p     = 1;
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

	bc_Qs        = struct();
	bc_Qs(1,1).p   = 0;
	bc_Qs(1,1).rhs = 0;
	bc_Qs(2,1).rhs = total_transport_engelund_hansen(drag2chezy(Cd),d_mm,Q0./(h0*w0),h0,w0)
	bc_Qs(2,1).p   = 1;

	% domain size
	L0 = h0/S0;
	Xi        = [0,2*h0/S0];

	nx               = 100;
	meta             = river_tide_test_metadata();
	opt              = meta.opt;
	opt.maxiter      = 100;
	opt.sopt.maxiter = 800;
	opt.dischargeisvariable = dischargeisvariable;

	hydrosolver    = BVPS_Characteristic();
	hydrosolver.xi = Xi;
	hydrosolver.nx = nx;

	morsolver        = Time_Stepper();
	morsolver.Ti     = [0,2000*Physics.SECONDS_PER_YEAR];
	% leapfrog-trapezoidal stable only for cfl <= 0.6
%	morsolver.cfl    = 0.6;
	morsolver.cfl    = 0.99;
	morsolver.scheme = 'upwind';
	morsolver.scheme = 'maccormack';
%	morsolver.scheme = 'leapfrog-trapezoidal';

	opt.stokes_order = 2;

	rt  = River_Tide_BVP( ...
				   'zb',          zb ...
				 , 'cd',          Cd ...
				 , 'width',       w0 ...
				 , 'omega',       omega ...
				 , 'opt',         opt ...
				 , 'hydrosolver', hydrosolver ...
				 , 'morsolver',   morsolver ...
				);
	rt.sediment.d_mm = d_mm;

	rt.bc    = bc;
	rt.bc_Qs      = bc_Qs;

tic();
	[t,zb] = rt.evolve_bed_level();
toc()
	%xc = mid(out.x);
	%zbfun = @(x) interp1(xc,zb(:,end),x,'linear');

	z1     = rt.z(1,1);
	z0     = rt.z(0,1);
%	z0i    = z_(:,1);
	z0i    = NaN;
	x      = rt.x(1);
	xc     = mid(x);
	z00    = S0*x;

	figure(1);
	clf();
	subplot(2,2,1)
	plot(mid(x),[zb(:,1),zb(:,end),rt.zb(1,mid(x)),rt.z(0,1,mid(x))])
%	ylabel('zb')
	legend('zb(0)','zb(T)');
	
	subplot(2,2,2)
	plot(mid(x),[-(zb(:,end)-zb(:,1))])
	hold on;
	plot(x,abs(z1))
	hold on;
	plot(x,[z0-z00,z0i-z00]);
	legend('-(zb(T)-zb(0))','|z_1|','z0(T)-S0 x','z0(0)-S0 x')

	subplot(2,2,3)
	d = [max(abs(diff(zb,[],2)))', max(abs(zb(:,1:end-1)-zb(:,end)))'];
	plot(mid(t)/(365*86400),d,'.-');
	legend('zb-zb(T-dT)','|zb-zb(T)|')
	xlabel('Years');
	title(['N = ' num2str(length(t))]);

	subplot(2,2,4)
	plot(mid(x),[diff(z0)./diff(x)]-S0)
	hold on;
	plot(mid(x),[abs(diff(z1))./diff(x)])
	plot(mid(x),[diff(inner2outer(zb(:,end)))./diff(x)]-S0)
%	hline(S0);
%	ylabel();		
	legend('dz0/dx-S0','dz1/dx','dzb/dx-S0');

	figure(2);
	subplot(2,1,1);
	t = NaN;
	iscentral = false;
	[Qs,trt,tas,tst] = rt.sediment_transport_(1,t,iscentral);
	Qs = Qs(2:end-1,:);
	Qs_rt     = Qs.*trt;
	Qs_stokes = Qs.*tst;
	Qs0 = Qs.*(1-trt-tst);
	plot(xc/L0,[Qs,Qs0,Qs_rt,Qs_stokes]);
%	abs(Qs([1:2 end-1:end])./Qs(end))
	ylabel('kg/s');
	legend('Total','River','RT','Stokes');

	Q = rt.out(1).Q;
	subplot(2,1,2)
	plot(x,abs(Q))
	ylabel('m^3/s');
	legend('Q0','|Qt|');

	h0=rt.h0(1)./h00;
	zs=rt.z(0);
	x=rt.x;
	S = cdiff(zs)./cdiff(x);
	cd=rt.cd(1);
	u0 = sqrt(9.81./cd.*h0.*S);
	Q0=rt.Q(0);
	u0(:,2) = Q0./h0;
	u1=abs(rt.Q(1))./h0;

	figure(3);
	subplot(2,2,1);
	clf();
	plot([u0,u1+1]);

	dzb_dt = rt.dzb_dt(0,zb(:,end),1);
	subplot(2,2,2);
	plot(dzb_dt);
	
	subplot(2,2,3);
	plot([h0/10,u1.^2+1]);
     
	%figure(4);
	subplot(2,2,4);
	plot(xc,Qs);

ozb = zb;
	x=rt.x;
 h=rt.h0(1);
 Q1=(rt.Q(1));
 Qs1=Qs-Qs0;
 figure(6);
 clf;
 p=-5;
 y = (h.^p-h(end).^p)./(h(1).^p-h(end).^p);
 x_ = rt.x*S0/h00;
 plot(x_,y,'linewidth',1,'color',[0,0,0]);
 hold on;
 plot(xc./L0,[Qs1./Qs1(1)],'--','linewidth',2.5);
 %legend('(h(x)-h_u)/(h(0)-h_u)',
 legend('(h(x)^{-5}-h_u^{-5})/(h(0)^{-5}-h_u^{-5})','(Qs_{rt}(x)/Qs_{rt}(0)');
%'(|Q_1(x)|/|Q_1(0)|)^2');
%  plot(mid(rt.x),[2-(zb(:,end)-zb(:,1))./(zb(end,end)-zb(end,1))]);
 xlim([0,0.5*sqrt(z10)]);
 xlabel('S_0 x/h_u');
 ylim([-0.001,1.001]);

	xl=xlim;
	xlim(xl);
	Q1=abs(rt.Q(1));
	xi=interp1(2*log(Q1./Q1(1)),rt.x,([-1:-1:-5]),'linear')*S0/h00;
	ax=addx();
	set(ax,'xlim',xl,'xtick',xi,'xticklabel',num2str((-1:-1:-5)'));
	xlabel({' ','ln(Qs_{rt}(x)/Qs_{rt}(0))'});
	set(gca,'ytick',[]);
	%pdfprint(11, 'img/tidal-river-schematic.pdf',2.5,[],'pdf',0.1);
if (pflag)
	pdfprint(6,'img/tidal-river-river-tide-transport-and-depth-along-river.pdf',2.5,[],'pdf',0.1);
	plot_;
end

	rtm_plot();

figure(12);
plot(diff(t(1:end-1)/86400/365.25))
ylabel('dt')
% plot(x,[abs((Q1./Q1(1)).^2),Qs1./Qs1(1)],'--','linewidth',1.5);
% legend('(h(x)-h_u)/(h(0)-h_u)','(|Q_1(x)|/|Q_1(0)|)^2');
%xlim([0,2e5])

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

