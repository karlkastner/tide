% Tue 25 Aug 18:21:53 +08 2020
function [out,t,ozb,rt] = river_tide_morphodyanmics_test_04_sesons(x0,zb0,pflag)
	if (nargin()<3)
		pflag = false;
	end
	out.id   = 5;
	out.name = 'integration over the hydrograph (seasons)';

	% number of points for integration over hydrograph
	iorder = [1,3,1,3];

	% tidal amplitude
	% TODO, why does it fail for e = 0?
	e = sqrt(eps);
	z10   = [e,e,1,1];

	% river discharge
	Qlim      = [-2,-20];
	Q0        = formative_discharge(Qlim(1),Qlim(2),'chezy'); 

	% width of channel
	w0        = 1;

	% drag/friction coefficient
	Cd        = 2.5e-3;

	% asymptotic depth in upstream reach
	h0        = 10;
	h00       = h0;

	% asymptotic bed slope in upstream reach
	S0         = -normal_flow_slope(Q0,h0,w0,drag2chezy(Cd));

	% initial bed level of channel
	zb         = @(x) -h0 + S0*x;

	d_mm = 0.2;

	% base frequency
	T         = Constant.SECONDS_PER_DAY;
	omega     = 2*pi/T;

	bc        = struct();
	% mean sea level
	bc(1,1).var = 'z';
	bc(1,1).rhs =  0;
	% Dirichlet condition
	bc(1,1).p   =  1;
	% river discharge
	%bc(2,1).rhs   = Q0;
	bc(2,1).Qseason  = Qlim;
	bc(2,1).var   = 'Q';
	% Dirichlet condition
	bc(2,1).p     = 1;
	% wave entering from left
	bc(1,2).var = 'z';
	bc(1,2).rhs = z10(1);
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

function Qs_ = Qsfun(t,Q0_)
	h0_ = normal_flow_depth(Q0_,w0,Cd,S0,'Cd');
	Qs_ = total_transport_engelund_hansen(drag2chezy(Cd),d_mm,Q0_./(h0_*w0),[],w0);
end

	bc_Qs(2,1).rhsfun = @Qsfun; 
	bc_Qs(2,1).p   = 1;

	% domain size
	L0        = h0/S0;
	Xi        = [0,1.2*L0];

	nx           = 25;

	meta             = river_tide_test_metadata();
	opt              = meta.opt;
	opt.maxiter      = 100;
	opt.sopt.maxiter = 800;
	opt.dischargeisvariable = true;

	hydrosolver = BVPS_Characteristic();
	hydrosolver.xi = Xi;
	hydrosolver.nx = nx;
	hydrosolver.opt.dischargeisvariable = true;

	morsolver        = Time_Stepper();
	morsolver.Ti     = [0,2000*Physics.SECONDS_PER_YEAR];
	morsolver.cfl    = 0.99;
	morsolver.scheme = 'upwind';

	rt = River_Tide_BVP( ...
				   'zb',      zb ...
				 , 'cd',      Cd ...
				 , 'width',   w0 ...
				 , 'omega',       omega ...
				 , 'opt',         opt ...
				 , 'hydrosolver',   hydrosolver ...
				 , 'morsolver',   morsolver ...
				);
	rt.sediment.d_mm = d_mm;

	rt.bc       = bc;
	rt.bc_Qs    = bc_Qs;

	figure(1);
	clf();

for idx=1:length(iorder)

tic();
	rt.opt.iorder = iorder(idx);
	rt.bc(1,2,1).var = 'z';
	rt.bc(1,2,1).rhs = z10(idx);
	rt.bc(1,2,1).p   = [1,0];
	rt.bc(1,2,1).q   = [1,0];
	[t,zb] = rt.evolve_bed_level();
toc()
	plot(zb(:,end));
	drawnow();
	hold on
end

if (0)
	%xc = mid(rt.x);
	%zbfun = @(x) interp1(xc,zb(:,end),x,'linear');

	z1     = rt.z(1);
	z0     = rt.z(0);
	%z0i    = z_(:,1);
	x      = rt.x;
	z00    = S0*x;

	figure(1);
	clf();
	subplot(2,2,1)
	plot(mid(x),[zb(:,1),zb(:,end)]);
%,zbfun(mid(x))])
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
	plot(mid(x),[diff(z0)./diff(rt.x)]-S0)
	hold on;
	plot(mid(x),[abs(diff(z1))./diff(rt.x)])
	plot(mid(x),[diff(inner2rter(zb(:,end)))./diff(rt.x)]-S0)
%	hline(S0);
%	ylabel();		
	legend('dz0/dx-S0','dz1/dx','dzb/dx-S0');

	figure(2);
	subplot(2,1,1);
	[Qs,Qs0] = rt.sediment_transport(0,0);
	Qs = Qs(2:end-1,:);
	plot(x/L,[Qs,Qs0,Qs-Qs0]);
%	abs(Qs([1:2 end-1:end])./Qs(end))
	Q = rt.Q_;
	ylabel('kg/s');
	legend('Qs total','Qs r','Qs rt');

	subplot(2,1,2)
	plot(x,[Q(:,1),abs(Q(:,2))])
	ylabel('m^3/s');
	legend('Q0','|Qt|');

	h0=rt.h0./h00;
	zs=rt.z(0);
	x=rt.x;
	S = cdiff(zs)./cdiff(x);
	cd=rt.cd;
	u0 = sqrt(9.81./cd.*h0.*S);
	Q0=rt.Q(0);
	u0(:,2) = Q0./h0;
	u1=abs(rt.Q(1))./h0;

	figure(3);
	subplot(2,2,1);
	clf
	plot([u0,u1+1]);

	dzb_dt = rt.dzb_dt(0,zb(:,end),1);
	subplot(2,2,2);
	plot(dzb_dt);
	
	subplot(2,2,3);
	plot([h0/10,u1.^2+1]);
     
	%figure(4);
	subplot(2,2,4);
	plot(Qs); 

 ozb = zb;
 x=rt.x;
 h=rt.h0;
 Q1=(rt.Q(1));
 Qs1=Qs-Qs0;
 figure(6);
 clf;
  p=-5; y = (h.^p-h(end).^p)./(h(1).^p-h(end).^p);
 x_ = rt.x*S0/h00;
 plot(x_,y,'linewidth',1,'color',[0,0,0]);
 hold on;
 plot(x_,[Qs1./Qs1(1)],'--','linewidth',2.5);
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
if (all(isfinite(xi)))
	set(ax,'xlim',xl,'xtick',xi,'xticklabel',num2str((-1:-1:-5)'));
end
	xlabel({' ','ln(Qs_{rt}(x)/Qs_{rt}(0))'});
	set(gca,'ytick',[]);
	%pdfprint(11, 'img/tidal-river-schematic.pdf',2.5,[],'pdf',0.1);
if (pflag)
	pdfprint(6,'img/tidal-river-river-tide-transport-and-depth-along-river.pdf',2.5,[],'pdf',0.1);
	plot_;
end

	rtm_plot();

figure(12);
plot(   mid(t(1:end-1))/Physics.SECONDS_PER_YEAR, ...
	diff(t(1:end-1)/Physics.SECONDS_PER_YEAR) );
ylabel('dt/y');

end

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

