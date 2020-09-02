% Mon 14 Oct 18:36:06 PST 2019
function river_tide_test_plot(tid,rt,z,name,pflag)
	x = rt.x;
	if (pflag)
		[Qlr,zlr] = rt.decompose();

		namedfigure(tid,['Test: ',name]);
		clf();
		subplot(2,4,1);
		plot(rt.x,[rt.zb(),rt.z(0)]);
		legend('z_b','z_0');

		subplot(2,4,5);
		plot(rt.x,[rt.width()]);
		legend('w_0');

		subplot(2,4,2);
		plot(rt.x,abs(rt.z(1)));
		hold on;
		plot(x,abs(z),'--');
		legend('|z_1|');
		ylim([0,max(ylim())]);

		subplot(2,4,6)
		plot(rt.x,angle(rt.z(1)));
		hold on;
		plot(x,angle(z),'--');
		legend('arg(z_1)');
		ylim(pi*[-1,1]);

		subplot(2,4,3)
		u1 = rt.velocity(1);
		plot(rt.x,abs(u1));
		legend('|u_1|');
		ylim([0,max(ylim())]);

		subplot(2,4,7)
%		u2 = rt.velocity(2);
		plot(rt.x,angle(u1));
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

		if (size(z,2)>2)
			namedfigure(200+tid,['Test: ',name]);
			clf();
			for jdx=2:min(size(z,2),5)
				subplot(2,4,jdx-1);
				plot(rt.x,abs(rt.out.z(:,jdx)));
				hold on;
				if (size(z,2)>=jdx)
					plot(rt.x,abs(z(:,jdx)),'--');
				end
				legend(sprintf('|z_%d|',jdx-1));
			
				subplot(2,4,4+jdx-1);
				plot(rt.x,angle(rt.out.z(:,jdx)));
				hold on;
				if (size(z,2)>=jdx)
					plot(rt.x,angle(z(:,jdx)),'--');
				end
				legend(sprintf('arg(z_%d)',jdx-1));
			end
		end % if
	end % if pflag
end % river_tide_test_plot

