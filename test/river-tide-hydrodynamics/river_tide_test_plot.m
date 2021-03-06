% Mon 14 Oct 18:36:06 PST 2019
function river_tide_test_plot(tid,rt,z,name,pflag)
	if (pflag)
		x = rt.channel(1).x;

		[Qlr,zlr] = rt.channel(1).decompose();

		namedfigure(tid,['Test: ',name]);
		clf();
		subplot(2,4,1);
		plot(rt.channel(1).x,[rt.channel(1).zb(),rt.channel(1).waterlevel(0)]);
		legend('z_b','z_0');

		subplot(2,4,5);
		plot(rt.channel(1).x,[rt.channel(1).width()]);
		legend('w_0');

		subplot(2,4,2);
		plot(rt.channel(1).x,abs(rt.channel(1).waterlevel(1)));
		hold on;
		plot(rt.channel(1).x,abs(z),'--');
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
%		u2 = rt.velocity(2);
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

		if (size(z,2)>2)
			namedfigure(200+tid,['Test: ',name]);
			clf();
			for jdx=2:min(size(z,2),5)
				subplot(2,4,jdx-1);
				plot(rt.channel(1).x,abs(rt.channel(1).z(:,jdx)));
				hold on;
				if (size(z,2)>=jdx)
					plot(rt.channel(1).x,abs(z(:,jdx)),'--');
				end
				legend(sprintf('|z_%d|',jdx-1));
			
				subplot(2,4,4+jdx-1);
				plot(rt.channel(1).x,angle(rt.channel(1).z(:,jdx)));
				hold on;
				if (size(z,2)>=jdx)
					plot(rt.channel(1).x,angle(z(:,jdx)),'--');
				end
				legend(sprintf('arg(z_%d)',jdx-1));
			end
		end % if
	end % if pflag
end % river_tide_test_plot

