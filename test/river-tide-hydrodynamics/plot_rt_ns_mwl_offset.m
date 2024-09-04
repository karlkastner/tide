% Tue 13 Sep 17:57:36 CEST 2022
if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;
ps = 3.5;
col = {'k','r'}
lw  = 1;
% TODO plot at x0 where zs-hysteresis is largest
% TODO plot at time where hysteresis is largest (za at max hysteresis ztrms)
% TODO set dt-out to 30 min
%	ts : time of spring
%	tqmax : time of max current
%	then plot at tqmax and tqmin and when t=qbar
% -> for which range at the river mouth is the q_ta-range highest?
% series analysis : magnitude mwl-offset and hystersis vs Q and the 'magical' r_n/L parameter
% sin(t/14)*(sin(t/2)) = sin(t/(3/7)) + sin(t/(4/7))

id = 80;
nf = 2*24+1;

ls = {'-','--'};
folder = {'~/large/phd/src/lib/tide/test/mat/river-tide-nonstationary-302';
	  '~/large/phd/src/lib/tide/test/mat/river-tide-nonstationary-302-lake_' }
for ldx=1:2

m = Delft3D_Map('folder',folder{ldx});
m.init();

t  = m.time;
x=m.Y;
x=x(1,:);
nt = m.nt; %length(t);
zs = m.zs;
zs = squeeze(zs(:,2,:));
zb = m.zb;
zb = squeeze(zb(1,2,:));
hn = -zb(1);
dzb_dx = diff(zb(1:2))./diff(x(1:2));
L = hn./dzb_dx;
Cz = m.Chezy_u;
Cz=squeeze(Cz(1,2,:));
y =m.X;
w =abs(diff(y(1:2,1)));
Qn = -normal_flow_discharge(hn,w,Cz(1),dzb_dx)
q = m.discharge;
q = squeeze(q(:,1,:));

if (1)
	%td0=8*28*24;
	%nt_=round((nt-td0)/(28*24))*28*24;
	%td0=nt-nt_+1;
	k = 6;
	td0 = nt+1-k*24*14;
	tdx = td0:nt;
	% length(td0:nt)/(24*14)
	ft  = fourier_axis(t(tdx));
	fzs = fft(zs(tdx,:));
	fq  = fft(q(tdx,:));
	fdx = abs(ft) <= 0.5;
	fzs(fdx,:) = 0;
	fq(fdx,:)  = 0;
	zt = NaN(size(zs));
	qt  = NaN(size(q));
	zt(tdx,:) = ifft(fzs);
	qt(tdx,:) = ifft(fq);
	za        = zs - zt;
	qa	  = q  - qt;
	fzt_rms   = fft(zt(tdx,:).^2);
	fzt_rms(~fdx,:) = 0;
	zt_rms = NaN(size(zs));
	zt_rms(tdx,:) = sqrt(2*ifft(fzt_rms));	
	fqt_rms   = fft(qt(tdx,:).^2);
	fqt_rms(~fdx,:) = 0;
	qt_rms = NaN(size(q));
	qt_rms(tdx,:) = sqrt(2*ifft(fqt_rms));

	tdx = tdx((k-1)/2*14*24+(1:14*24))
%	flag = true(size(t));
%	flag(tdx) = false;
	t = t(tdx)-t(tdx(1));
	qt_rms = qt_rms(tdx,:);
	zt_rms = zt_rms(tdx,:);
	qt = qt(tdx,:);
	zt = zt(tdx,:);
	qa = qa(tdx,:);
	za = za(tdx,:);
	q = q(tdx,:);
	zs = zs(tdx,:);

	% shift neap to start of time series
	[zt_rms_min,mdx]=min(zt_rms(:,1));

	zs     = circshift(zs,-mdx+1);
	zt_rms = circshift(zt_rms,-mdx+1);
	qt_rms = circshift(qt_rms,-mdx+1);
	zt = circshift(zt,-mdx+1);
	qt = circshift(qt,-mdx+1);
	za = circshift(za,-mdx+1);
	qa = circshift(qa,-mdx+1);
	

%	qt_rms(flag,:) = NaN;
%	zt_rms(flag,:) = NaN;
%	qt(flag,:) = NaN;
%	zt(flag,:) = NaN;
%	qa(flag,:) = NaN;
%	za(flag,:) = NaN;
%	qt_rms(tdx(1:(k-1)/2*14*24)) = NaN;
%	zt_rms(tdx(1:(k-1)/2*14*24)) = NaN;
%	qt_rms(nt-(k-1)/2*14*24+1:end,:) = NaN;	
%	zt_rms(nt-(k-1)/2*14*24+1:end,:) = NaN;	
%	qa_rms(tdx(1:(k-1)/2*14*24)) = NaN;
%	za_rms(tdx(1:(k-1)/2*14*24)) = NaN;
%	qa_rms(nt-(k-1)/2*14*24+1:end,:) = NaN;	
%	za_rms(nt-(k-1)/2*14*24+1:end,:) = NaN;	
%	qa(tdx(1:(k-1)/2*14*24)) = NaN;
%	za(tdx(1:(k-1)/2*14*24)) = NaN;
%	qa(nt-(k-1)/2*14*24+1:end,:) = NaN;	
%	za(nt-(k-1)/2*14*24+1:end,:) = NaN;	
else
	q(1:2*24*14,:) = NaN;
	q(end-nf+1:end,:) = NaN;
	qa = trifilt1(q,nf);
	qt = q-qa;
	zs(1:2*24*14,:) = NaN;
	zs(end-nf+1:end,:) = NaN;
	za = trifilt1(zs,nf);
	zt = zs-za;
	zt_rms = sqrt(2*trifilt1(zt.^2,nf));
	qt_rms = sqrt(2*trifilt1(qt.^2,nf));
end
%q(1:24*14,:) = NaN;

% average level over th spring-neap cycle
za_     = mean(za,1);
% volume (prism) of the msf-tide
dx = x(2)-x(1);
dza = za-za_;
I_za = sum(dza,2)*dx; 

% over time
za_hyst_t = -(I_za(1:7*24) - I_za(end-1:-1:end-7*24));
[za_hyst_max_t,tdx_max] = max(za_hyst_t);

za_hyst_x = -(dza(1:7*24,:) - dza(end-1:-1:end-7*24,:));
za_hyst_x_max = max(za_hyst_x);
%[za_h_min,tdx_max] = max(za_hyst);
%[za_h_min,tdx_min] = min(za_hyst);

[za_h_x_mm, idx_max] = max(za_hyst_x_max)

zt0 = max(zt_rms(:,1));

A   = fourier_matrix_exp(14,t);
cq  = A \ q;
czt = A \ zt;
cz  = A \ zs;
[czt_max,idx_max] = max(abs(czt(2,:)));
if (1==ldx)
	id = idx_max;
	%id = 116;
	id = round(0.5*L/dx);
end

fdx = 1:24:length(t);

if (1 == ldx)
splitfigure([2,3],[1,1],fflag);;
cla();
plot(t,zs(:,[id]),'linewidth',lw);
hold on
plot(t,za(:,[id]),'linewidth',lw);
hold on
if (0)
	plot(m.time,zt_rms(:,[1]))
	plot(m.time,za(:,[1]));
	%plot(m.time,zt_rms(:,[id]))
end
%xlim(tlim-tlim(1))
legend('$z_s(x)$','$z_a$','interpreter','latex')
ylabel('Water level $z/\mathrm{m}$','interpreter','latex');
xlabel('day','interpreter','latex');
%ylim([0.6,1.6]);
vline(7,'color','k');
text(0.75,max(zs(:,[id])),sprintf('$x/L$ = %0.2f', round(x(id)/L,3)),'interpreter','latex');
[mv,mdx]=max(za(:,id));
vline(t(mdx),'color','r');
ylim([1.2,1.8]);
xlim(limits(t))

%fdx = (nt-24*18+1):24:nt;
%fdx = (nt-24*18+1):24:nt;
splitfigure([2,3],[1,2],fflag);;
cla();
y = zt_rms(:,:)./zt_rms(:,1);
plot(za(:,id),y(:,id));
hold on
quiver(za(fdx(1:end-1),id),y(fdx(1:end-1),id),diff(za(fdx,id)),diff(y(fdx,id)),0);

splitfigure([2,3],[1,3],fflag);;
cla();
y = zt_rms(:,:)./zt_rms(:,1);
plot(qa(:,1)/Qn,y(:,id));
hold on
quiver(qa(fdx(1:end-1),1)/Qn,y(fdx(1:end-1),id),diff(qa(fdx,1)/Qn),diff(y(fdx,id)),0);



splitfigure([2,3],[1,4],fflag);
cla();
plot(qa(:,[1,id])/Qn,zt_rms(:,id)./zt0)
hold on
legend('vs downstream q','vs at a station q')
quiver(qa(fdx(1:end-1),1)/Qn,zt_rms(fdx(1:end-1),id)./zt0,diff(qa(fdx,1))/Qn,diff(zt_rms(fdx,id)./zt0),0);
quiver(qa(fdx(1:end-1),id)/Qn,zt_rms(fdx(1:end-1),id)./zt0,diff(qa(fdx,id))/Qn,diff(zt_rms(fdx,id)./zt0),0);
xlabel('Tidally-averaged discharge $Q_a/Q_n$','interpreter','latex');
ylabel('$|z_t(x)|/\mathrm{max}(|z_t(0)|)$','interpreter','latex');
end
legend('without lake','with lake','location','southeast');

splitfigure([2,3],[1,5],fflag);
if (1 == ldx) cla(); end
id_ = id;
za_scale = x(id)/L*hn;
plot(za(:,id_)/za_scale,zt_rms(:,id_)./zt0,[col{ldx},'-'],'linestyle',ls{1},'linewidth',lw);
hold on
for idx=1:length(id_)
	qh = quiver(za(fdx(1:end-1),id_(idx))/za_scale,zt_rms(fdx(1:end-1),id_(idx)),diff(za(fdx,id_(idx)))/za_scale,diff(zt_rms(fdx,id_(idx))),0,'color',col{ldx});
	qh.HandleVisibility = 'off';
	if (1 == ldx)
		text(0.4,0.14,sprintf('x/L = %0.2f', round(x(id)/L,3)));
	end
end
xlabel('Tidally averaged water level $z_a(x) / z_n(x)$','interpreter','latex');
ylabel('Tidal range $|z_t(x)|/\mathrm{max}(|z_t(0)|)$','interpreter','latex');
legend('without lake','with lake','location','southeast');

if (1 == ldx)
splitfigure([2,3],[1,6],fflag);
cla();
%plot(t,I_za/L)
ax=plotyy(t,I_za/L,t,qa(:,1)/Qn)
xlabel('t/day','interpreter','latex');
ylabel('int (za-\bar za) dt/L','interpreter','latex');
ylabel(ax(2),'Qa(0)/Qn','interpreter','latex');
ylim(ax(1),0.151*[-1,1])
vline(7);
end

splitfigure([2,3],[2,1],fflag);
if (1 == ldx) cla(); end
%tdx_ = tdx(1) + [0,7*24];
tdx_ = [1,tdx_max, 7*24, 14*24-tdx_max];
if (0)
	plot(inner2outer(x)/L,za(tdx_(1:4),:)'/hn)
	legend('neap','rising','spring','falling')
else
	plot(inner2outer(x)/L,[min(za,[],1)',max(za,[],1)']/hn,'linestyle',ls{ldx},'linewidth',lw);
	hold on
	set(gca,'colororderindex',1)
end
xlim([0,1])
ylabel('Tidally-averaged water level $z_a/h_n$','interpreter','latex');
xlabel('Distance from mouth $x/L$','interpreter','latex');
legend('min','max','location','northwest');

splitfigure([2,3],[2,2],fflag);
if (1 == ldx)
	cla
	plot(inner2outer(x)/L,abs(cz(2,:))'/zt0,'linewidth',lw);
	yyaxis right
	cla
	h = plot(inner2outer(x)/L,[wrapToPi(angle(cz(2,:))-angle(cz(2,1)))'*14/(2*pi)],'linewidth',lw);
	h.HandleVisibility = 'off';
else
	yyaxis left
	hold on
	set(gca,'colororderindex',1);
	plot(inner2outer(x)/L,abs(cz(2,:))'/zt0,'linestyle',ls{ldx},'linewidth',lw);
	ylabel('Subtidal Water Level Amplitude $z_a(x)/|z_t(0)|$','interpreter','latex');
	yyaxis right
	hold on
	set(gca,'colororderindex',1);
	h = plot(inner2outer(x)/L,[wrapToPi(angle(cz(2,:))-angle(cz(2,1)))'*14/(2*pi)],'linestyle',ls{ldx},'linewidth',lw);
	h.HandleVisibility = 'off';
	hline(0);
	set(gca,'ytick',[-7:7])
	ylabel('Subtidal Water Level Phase lag $\mathrm{arg}(z_a) T_a/(2 pi) /$ day','interpreter','latex');
	ylim([-1,1]*14/2);
end
xlim([0,1]);
xlabel('Distance from mouth $x/L$','interpreter','latex');
legend('without lake','with lake');

splitfigure([2,3],[2,3],fflag);
if (1 == ldx)
	cla
	plot(x(2:end-1)/L,-abs(cq(2,:))'/Qn,'linewidth',lw);
	yyaxis right
	cla
	plot(x(2:end-1)/L,[wrapToPi(angle(cq(2,:)/Qn)-angle(cz(2,1)))'*14/(2*pi)],'linewidth',lw);
else
	yyaxis left
	hold on
	set(gca,'colororderindex',1);
	plot(x(2:end-1)/L,-abs(cq(2,:))'/Qn,'linestyle',ls{ldx},'linewidth',lw);
	ylabel('Subtidal Water Level Amplitude $|z_a(x)|/z_t(0)$','interpreter','latex');
	ylabel('Subtidal Discharge Amplitude $|Q_a(x)|/Q_n$','interpreter','latex');
	yyaxis right
	hold on
	set(gca,'colororderindex',1);
	plot(x(2:end-1)/L,[wrapToPi(angle(cq(2,:)/Qn)-angle(cz(2,1)))'*14/(2*pi)],'linestyle',ls{ldx},'linewidth',lw);
	hline(0);
	ylim([-1,1]*14/2);
	ylabel('Subtidal discharge phase lag $|\mathrm{arg}(Q_a) T_a/(2\pi) /$ day','interpreter','latex');
	set(gca,'ytick',[-7:7])
end
xlim([0,1]);
xlabel('Distance from mouth $x/L$','interpreter','latex');
legend('without lake','with lake');

splitfigure([2,3],[2,4],fflag);
if (1 == ldx) cla(); end
plot(inner2outer(x)/L,[max(zt_rms,[],1)/zt0],col{ldx},'linestyle',ls{ldx},'linewidth',lw)
%plot(inner2outer(x)/L,[min(zt_rms,[],2)/zt0,max(zt_rms,[],2)/zt0],'linestyle',ls{idx})
hold on
set(gca,'colororderindex',1);
xlabel('Distance from mouth $x/L$','interpreter','latex');
ylabel('Maximum Tidal range $|z_t(x)|/\mathrm{max}(|z_t(0)|)$','interpreter','latex');
xlim([0,1])
legend('without lake','with lake');

if (0)
plot(inner2outer(x)/L,zt_rms(tdx_(1:3),:)'/zt0)
hold on
set(gca,'colororderindex',4);
plot(inner2outer(x)/L,NaN*zt_rms(tdx_(4),:)'/zt0,'-')
set(gca,'colororderindex',4);
plot(inner2outer(x)/L,zt_rms(tdx_(4),:)'/zt0,'--','HandleVisibility','off')
legend('neap','rising','spring','falling')
xlabel('Distance from mouth $x/L$','interpreter','latex');
ylabel('Admittance $|z_t(x)|/\mathrm{max}(|z_t(0)|)$','interpreter','latex');
xlim([0,1])
end

splitfigure([2,3],[2,5],fflag);
if (1 == ldx) cla; end 
plot((x(2:end-1))/L,-qt_rms(tdx_(1:3),:)'/Qn)
hold on
set(gca,'colororderindex',4);
plot(x(2:end-1)/L,-NaN*qt_rms(tdx_(4),:)'/Qn,'-')
set(gca,'colororderindex',4);
plot(x(2:end-1)/L,-qt_rms(tdx_(4),:)'/Qn,'--','HandleVisibility','off')
legend('neap','rising','spring','falling')
xlabel('Distance from mouth $x/L$','interpreter','latex');
ylabel('Admittance $|Q_t(x)|/Q_n$','interpreter','latex');
xlim([0,1])

splitfigure([2,3],[2,6],fflag);
cla();
plotyy(t,zt_rms(:,1)/zt0,t,qa(:,1)/Qn);
xlabel('t/day','interpreter','latex');
% plot(inner2outer(x)/L,za_hyst_x_max/zt0);
xlabel('x/l','interpreter','latex');
% TODO normalize by zt0^2

if (0)
%q=q-nanmean(q);
plot(m.time,[q(:,1),trifilt1(q(:,1),2*24)])
end

end

if (pflag)
	pdfprint(11,'img/rt-ns-neap-spring-zs-vs-time.pdf',ps);
	pdfprint(15,'img/rt-ns-neap-spring-z_t-vs-z_ta.pdf',ps);
	pdfprint(16,'img/rt-ns-neap-spring-mfs-prism-vs-t.pdf',ps);
	pdfprint(21,'img/rt-ns-neap-spring-z_ta-along-channel.pdf',ps);
	pdfprint(22,'img/rt-ns-neap-spring-z-subtidal-along-channel.pdf',ps);
	pdfprint(23,'img/rt-ns-neap-spring-q-subtidal-along-channel.pdf',ps);
	pdfprint(24,'img/rt-ns-neap-spring-z_t-along-channel.pdf',ps);
end

