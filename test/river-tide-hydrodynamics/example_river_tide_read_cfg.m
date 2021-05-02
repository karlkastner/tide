% Thu 10 Oct 12:52:41 PST 2019
% exmaple river tide

	pflag = 1;
	tid   = 1;
	name  = 'Example River\_Tide';





	% define a channel network
	rt = River_Tide_Network();
	rt.read_cfg('river-tide-example.cfg');

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

