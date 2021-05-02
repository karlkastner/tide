% Thu 10 Oct 12:52:41 PST 2019
% exmaple river tide

	pflag = 1;
	tid = 1;
	name  = 'Example River\_Tide';

	% seasonnaly averaged river discharge
	% not an input to river tide, only here to set the channel slope
	Q0_        = -10;

	% instantaneous river discharge
	Q0        = 0.75*Q0_;

	% width of channel
	w0        = 1;
	wfun      = @(x)   w0*ones(size(x));

	% drag/friction coefficient
	Cd        = 2.5e-3;
	cdfun     = @(x)  Cd*ones(size(x));

	% bed level of channel
	h0        = 8;
	S0         = normal_flow_slope(-Q0_,h0,w0,drag2chezy(Cd));
	zbfun     = @(x) -h0 + S0*x; %*ones(size(x));

	% boundary condition
	bc        = struct();

	% mean sea level
	bc(1,1).var = 'z';
	bc(1,1).rhs = 0;
	% Dirichlet condition
	bc(1,1).p   = 1;
	% river discharge
	bc(2,1).var = 'Q';
	bc(2,1).rhs = Q0;
	% Dirichlet condition
	bc(2,1).p   = 1;

	% wave entering from left
	bc(1,2).var = 'z';
	z10         = 1; %sqrt(eps);
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

	% length of river section to compute tide for
	L         = 2.5e5;

	% solver of boundary value problem
	%opt.solver = @bvp2c;

	% number of points along channel
	nx     = 200;

	% define a channel
	channel = River_Tide_Channel( ...
			 'zb', zbfun ...
			,'cd', cdfun ...
			,'width', wfun ...
		 	,'L',   L ...
			,'nx',  nx ...
			,'bc',  bc ...
	);

	% define a channel network
	rt = River_Tide_Network( ...
		   'channel',     channel ...
		  ,'omega',       omega ...
		);

	% initialize
	rt.init();

	% solve
	rt.solve();

	% check continuity
	rmse = rt.check_continuity();

	% compare to analytical solution
	g = Constant.gravity;
	c0 = sqrt(g*h0);
	k0 = omega/c0;
	x = rt.channel(1).x;

	r = (1+1i)*sqrt(-Cd.*omega.*Q0/w0./(g*h0.^3));
	z = z10*exp(-r*x);

	rmse(2)  = rms(rt.channel(1).waterlevel(1)-z);
	% err ~ C*df^2/dx^2*dx^2, where C sufficiently large constant
	nres_ = rms(cdiff(z,2));
	fail = (rmse(1)/z10 > 0.05) || (rmse(2) > 10*nres_);
	
	if (pflag)
		[Qlr,zlr] = rt.channel(1).decompose();
		namedfigure(tid,['Test: ',name]);
		clf();
		subplot(2,4,1);
		plot(rt.channel(1).x,[zbfun(x),rt.channel(1).waterlevel(0)]);
		legend('z_b','z_0');

		subplot(2,4,5);
		plot(rt.channel(1).x,[wfun(x)]);
		legend('w_0');

		subplot(2,4,2);
		plot(rt.channel(1).x,abs(rt.channel(1).waterlevel(1)));
		hold on;
		plot(x,abs(z),'--');
		legend('|z_1|');
		ylim([0,max(ylim())]);

		subplot(2,4,6)
		plot(rt.channel(1).x,angle(rt.channel(1).waterlevel(1)));
		hold on;
		plot(x,angle(z),'--');
		legend('arg(z_1)');
		ylim(pi*[-1,1]);

		subplot(2,4,3)
		u1 = rt.channel(1).velocity(1);
		plot(rt.channel(1).x,abs(u1));
		legend('|u_1|');
		ylim([0,max(ylim())]);

		subplot(2,4,7)
		plot(rt.channel(1).x,angle(u1));
		legend('arg(u_1)');
		ylim(pi*[-1,1]);

		subplot(2,4,4)
		plot(rt.channel(1).x,abs(rt.channel(1).discharge(1)));
		legend('|Q_1|');
		ylim([0,max(ylim())]);

		subplot(2,4,8)
		plot(rt.channel(1).x,angle(rt.channel(1).discharge(1)));
		legend('arg(Q_1)');
		ylim(pi*[-1,1]);

		figure(100+tid);
		clf();
		subplot(2,2,1);
		plot(rt.channel(1).x,abs(zlr));
		legend('z_1^-','z_1^+');
		subplot(2,2,2);
		plot(rt.channel(1).x,abs(Qlr));
		legend('Q_1^-','Q_1^+');
		subplot(2,2,3);
		plot(rt.channel(1).x,angle(zlr));
		ylim(pi*[-1,1]);
		subplot(2,2,4);
		plot(rt.channel(1).x,angle(Qlr));
		ylim(pi*[-1,1]);

	end % if pflag

