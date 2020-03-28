% Mon 14 Oct 18:36:06 PST 2019
function river_tide_test_plot(tid,rt,z,name,pflag)
	x = rt.x;
	zbfun = rt.fun.zb;
	wfun = rt.fun.width;
	%z = rt.z_;
	if (pflag)
		[Qlr,zlr] = rt.decompose();
		u_ = rt.velocity();

		namedfigure(tid,['Test: ',name]);
		clf();
		subplot(2,4,1);
		plot(rt.x,[zbfun(x),rt.z_(:,1)]);
		legend('z_b','z_0');

		subplot(2,4,5);
		plot(rt.x,[wfun(x)]);
		legend('w_0');

		subplot(2,4,2);
		plot(rt.x,abs(rt.z_(:,2)));
		hold on;
		plot(x,abs(z),'--');
		legend('|z_1|');
		ylim([0,max(ylim())]);

		subplot(2,4,6)
		plot(rt.x,angle(rt.z_(:,2)));
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
		plot(rt.x,abs(rt.Q_(:,2)));
		legend('|Q_1|');
		ylim([0,max(ylim())]);

		subplot(2,4,8)
		plot(rt.x,angle(rt.Q_(:,2)));
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

		if (size(z,2)>2)
			namedfigure(200+tid,['Test: ',name]);
			clf();
			for jdx=2:min(size(z,2),5)
				subplot(2,4,jdx-1);
				plot(rt.x,abs(rt.z_(:,jdx)));
				hold on;
				if (size(z,2)>=jdx)
					plot(rt.x,abs(z(:,jdx)),'--');
				end
				legend(sprintf('|z_%d|',jdx-1));
			
				subplot(2,4,4+jdx-1);
				plot(rt.x,angle(rt.z_(:,jdx)));
				hold on;
				if (size(z,2)>=jdx)
					plot(rt.x,angle(z(:,jdx)),'--');
				end
				legend(sprintf('arg(z_%d)',jdx-1));
			end
		end % if
	end % if pflag
end % river_tide_test_plot

