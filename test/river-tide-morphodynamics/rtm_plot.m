% 2020-08-24 14:32:03.916189156

function rtm_plot(rt,z10,h00,S0,pflag)

figure(11)
clf;

L = S0/h00;

x=rt.x;
 x_ = x/L;
 s=1;
 z0 = [  0+x*S0,rt.z(0)];
 zb =[-h00+x*S0,rt.zb(1)];
% amplify
%[z0(:,1)+s*(z0(:,2)-z0(:,1)),zb_(:,1)+s*(zb_(:,2)-zb_(:,1))]./h00
%rt2.zb(rt2.x)];
z0__ = z0(:,1) + s*(z0(end,2)-z0(end,1));
zb__ = zb(:,1) + s*(zb(end,2)-zb(end,1));
 plot(x_,[z0(:,1),zb(:,1)]./h00,'--','linewidth',1.5);
 hold on;
 set(gca,'colororderindex',1);
 plot(x_,[z0(:,2),zb(:,2)]./h00,'-','linewidth',1.5);
set(gca,'colororderindex',1);
plot(x_,[z0__,zb__]./h00,':','linewidth',1.5)
 ylabel('z/h_u');
 xlabel('S_0 x/h_u');
xlim([0,0.6]);
ylim([-1.501,0.651]);
dz = [diff(z0(end,:)), diff(zb(end,:))]
dz = mean(dz)
vline(-s*dz/S0 * (S0/h00))
hline(0);
e = z10/h00;
estr = num2str(e);
legend('z_s, z_1=0','z_b, z_1=0',['z_s, z_1= ' estr ' h_u'],['z_b, z_1=' estr ' h_u'],'location','southeast');

xl=xlim;
 Q1=abs(rt.Q(1));
 xi=interp1(2*log(Q1./Q1(1)),rt.x,([-1,-2,-3,-4,-5]),'linear')*S0/h00;
 ax=addx();
if (all(isfinite(xi)))
 set(ax,'xlim',xl,'xtick',xi,'xticklabel',num2str((-1:-1:-5)'));
end
 xlabel({' ','ln(Qs_{rt}(x)/Qs_{rt}(0))'}) 
set(gca,'ytick',[]); 
if (pflag) 
	pdfprint(11, 'img/tidal-river-schematic.pdf',2.5,[],'pdf',0.1);
end
