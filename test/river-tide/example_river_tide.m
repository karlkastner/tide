% Thu 10 Oct 12:52:41 PST 2019
% exmaple river tide

	tid = 1000;
	pflag = 1;
	name = 'Example River\_Tide';

	% river discharge
	Qu        = -10;
	Q0        = 0.5*Qu;
	% width of channel
	w0        = 1;
	wfun      = @(x)   w0*ones(size(x));
	% drag/friction coefficient
	cD        = 2.5e-3;
	cdfun     = @(x)  cD*ones(size(x));
	% bed level of channel
	h0        = 10;
	S0         = normal_flow_slope(-Q0,h0,w0,drag2chezy(cD));
	zbfun     = @(x) -h0 + S0*x; %*ones(size(x));
	% mean component
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
	z10       = sqrt(eps);
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
	% model for river tide
	opt.model_str = 'wave';
	% solver of boundary value problem
	opt.solver = @bvp2c;
	% number of points along channel
	opt.nx     = 200;
	% change of distance between points along channel 
	opt.xs     = 1; 

	% solve with model

	rt = River_Tide( ...
		   'fun.zb',      zbfun ...
		 , 'fun.cd',      cdfun ...
		 , 'fun.width',   wfun ...
		 , 'omega',       omega ...
		 , 'Xi',          Xi ...
		 ... % , 'bc',          bc ...
		 , 'opt',         opt ...
		);
	% TODO there's something fishy with passing bc directly
	rt.bc = bc;

	rt.init();
	rt.solve();

	% check ode
	rmse = rt.check_continuity();
	% compare to analytical solution
	g = Constant.gravity;
	c0 = sqrt(g*h0);
	k0 = omega/c0;
	x = rt.x;

	r = (1+1i)*sqrt(-cD.*omega.*Q0/w0./(g*h0.^3));
	z = z10*exp(-r*x);

	rmse(2)  = rms(rt.z(1))-z)
	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
	nres_ = rms(cdiff(z,2));
	fail = (rmse(1)/z10 > 0.05) || (rmse(2) > 10*nres_);
	
	if (pflag)
		[Qlr,zlr] = rt.decompose();
		u_ = rt.velocity();
		namedfigure(tid,['Test: ',name]);
		clf();
		subplot(2,4,1);
		plot(rt.x,[zbfun(x),rt.z(0)]);
		legend('z_b','z_0');

		subplot(2,4,5);
		plot(rt.x,[wfun(x)]);
		legend('w_0');

		subplot(2,4,2);
		plot(rt.x,abs(rt.z(1))));
		hold on;
		plot(x,abs(z),'--');
		legend('|z_1|');
		ylim([0,max(ylim())]);

		subplot(2,4,6)
		plot(rt.x,angle(rt.z(1))));
		hold on;
		plot(x,angle(z),'--');
		legend('arg(z_1)');
		ylim(pi*[-1,1]);

		subplot(2,4,3)
		plot(rt.x,abs(u_(:,2)));
		legend('|u_1|');
		ylim([0,max(ylim())]);

		subplot(2,4,7)
		plot(rt.x,angle(u_(:,2)));
		legend('arg(u_1)');
		ylim(pi*[-1,1]);

		subplot(2,4,4)
		plot(rt.x,abs(rt.Q(1)));
		legend('|Q_1|');
		ylim([0,max(ylim())]);

		subplot(2,4,8)
		plot(rt.x,angle(rt.Q(1)));
		legend('arg(Q_1)');
		ylim(pi*[-1,1]);

		figure(100+tid);
		clf();
		subplot(2,2,1);
		plot(rt.x,abs(zlr));
		legend('z_1^-','z_1^+');
		subplot(2,2,2);
		plot(rt.x,abs(Qlr));
		legend('Q_1^-','Q_1^+');
		subplot(2,2,3);
		plot(rt.x,angle(zlr));
		ylim(pi*[-1,1]);
		subplot(2,2,4);
		plot(rt.x,angle(Qlr));
		ylim(pi*[-1,1]);

	end % if pflag

