% Tue  6 Sep 12:01:15 CEST 2022
% -t(tdx(1))zs vs q
% range = b q^b(1 + a dq/dt)

if (~exist('pflag','var'))
	pflag = false;
end
fflag = pflag;
ps = 3.5;
nf = 2*24+1;
lw = 1;

q0 = 0.51;
T1 = 86400/2;
g  = Physics.gravity;

m     = Delft3D_Map();
m.folder = 'mat/river-tide-nonstationary-9';
m.folder = 'mat/river-tide-instationary-9';
m.init();

x=m.Y;
x=x(1,:)';

Cz = m.Chezy_u;
q  = m.discharge();
t  = m.time;
u  = m.u2;
zb = m.zb;
zs = m.zs;

Cz=squeeze(Cz(1,2,:));
q  = squeeze(q);
zs  = squeeze(zs(:,2,:));
u  = squeeze(u(:,2,:));
zb =squeeze(zb(1,2,:));
s = squeeze(zs(:,2,:));

dx = diff(x);
dzb_dx = diff(zb)./cdiff(x);
hn  = -zb(1);
y=m.X;
w=abs(diff(y(:,1)));
Qn  = normal_flow_discharge(-zb(1),w,Cz(1),dzb_dx(1));
An  = w*hn;
un  = Qn/An;
xbw = -zb(1)/dzb_dx(1);
ibw = find(x>=xbw,1,'first');


% maybe tidal admission here?
c  = sqrt(g*hn);
lambda0 = c*T1;
k0      = 2*pi/lambda0;

xmax = 1.2*xbw;
sx = 1./xbw;
L = 150;
%k= [2,round(25/3),round(50/3),round(100/3),round(150/3)];
xplot = [0.125,0.25,0.5,1]*xbw;
xplot_ = 0.33*xbw;
k= round(xplot/dx(1));
k_= round(xplot_/dx(1));
%[round(0.125*xbw/dx(1)),round(0.25*xbw/dx(1)),round(0.5*xbw/dx(1)),round(1.0*xbw/dx(1))];
n = length(k);

T0 = 15*20+2/24;
T0_ = -eps;
tdx = ((24*T0):m.nt)';
fx=fourier_axis(t(tdx));


% tide-river separtation
q = q/Qn;
%q0 = 3e3/Qn;
%q0 = 1e3/1e4;
%q0 = 1.2e3/1e4;
%q0 = 0.095
% TODO no magic numbers
 fsep = 1.0; % 0.5*f_tide
 fx   = fourier_axis(t(tdx));
 fdx  = (abs(fx) < fsep);
 fzs  = fft(zs(tdx,:));
 fq   = fft(q(tdx,:));
 fu   = fft(u(tdx,:));
 fzs(fdx,:) = 0;
 fq(fdx,:)  = 0;
 fu(fdx,:)  = 0;
 zst = NaN(size(zs));
 qt  = NaN(size(q));
 ut  = NaN(size(u));
 zst(tdx,:) = ifft(fzs);
 qt(tdx,:)  = ifft(fq);
 ut(tdx,:)  = ifft(fu);
% zst(end-nf:end,:) = NaN;
% qt(end-nf:end,:)  = NaN;
% ut(end-nf:end,:)  = NaN;
 zsr        = zs-zst;
 qr         = q-qt;
 ur         = u-ut;
 fzt2 = fft(2*zst(tdx,:).^2);
 fqt2 = fft(2*qt(tdx,:).^2);
 fut2 = fft(2*ut(tdx,:).^2);
 fzt2(~fdx,:) = 0;
 fqt2(~fdx,:)  = 0;
 fut2(~fdx,:)  = 0;
 zst_rms = NaN(size(zs));
 qt_rms  = NaN(size(q));
 ut_rms  = NaN(size(u));
 zst_rms(tdx,:) = sqrt(ifft(fzt2));
 qt_rms(tdx,:)  = sqrt(ifft(fqt2));
 ut_rms(tdx,:)  = sqrt(ifft(fut2));
if (0)
 nf = 2*24+1;
 zst_rms    = sqrt(2)*sqrt(trifilt1(zst.^2,nf));
 qt_rms     = sqrt(2)*sqrt(trifilt1(qt.^2,nf));
 ut_rms     = sqrt(2)*sqrt(trifilt1(ut.^2,nf));
end

% last cycle
nt = m.nt;
tdx = nt-24*15+1:nt;
t = t(tdx)-t(tdx(1));
zs = zs(tdx,:);
zsr = zsr(tdx,:);
zrt = zst(tdx,:);
zst_rms = zst_rms(tdx,:);
q = q(tdx,:);
qr = qr(tdx,:);
qt = qt(tdx,:);
qt_rms = qt_rms(tdx,:);
u = u(tdx,:);
ur = ur(tdx,:);
ut = ut(tdx,:);
ut_rms = ut_rms(tdx,:);

tdx=1:24*15;

% damping coefficient k = 1/z dz/dt
kz = -(diff(zst_rms,[],2)./cdiff(x)')./mid(zst_rms,2);
kq = -(diff(qt_rms,[],2)./diff(x(2:end-1))')./mid(qt_rms,2);
ku = -(diff(ut_rms,[],2)./cdiff(x)')./mid(ut_rms,2);

% tidal amplitude at river mouth
zst0 = nanmean(zst_rms(:,1));

% normalize
ur  = ur/un;
ut  = ut/un;
zst = zst / zst0;
zst_rms = zst_rms / zst0;


% q-zsr relation
%tdx_ = tdx(end-30*24+1:end-15*24);
tdx_ = tdx;
[qrmin,fdx]   = min(-qr(tdx_,ibw));
tdx_ = circshift(tdx_,-fdx);
[qrmax,maxdx] = max(-qr(tdx_,ibw));
qri = linspace(qrmin,qrmax)';

% interpolate
if (1)
zst_rms_i(:,:,1) = interp1(-qr(tdx_(1:maxdx),ibw),zst_rms(tdx_(1:maxdx),:),qri,'linear','extrap');
zst_rms_i(:,:,2) = interp1(-qr(tdx_(maxdx:end),ibw),zst_rms(tdx_(maxdx:end),:),qri,'linear','extrap');
else
zst_rms_i(:,:,1) = interp1(-qr(tdx_(1:maxdx),ibw)',zst_rms(tdx_(1:maxdx),:)',qri,'linear');
zst_rms_i(:,:,2) = interp1(-qr(tdx_(maxdx:end),ibw)',zst_rms(tdx_(maxdx:end),:)',qri,'linear');
end
if (0)
qri = min(-qr(tdx,:)) + linspace(0,1,100)'*range(-qr(tdx,:));
zsti = [];
for idx=2:size(qri,2)
	rising = find(cdiff(-qr(:,idx))>0 & t > T0);
	[qrs,fdx] = unique(-qr(rising,idx));
	zsti(:,idx,1) = interp1(qrs,zst(rising(fdx),idx),qri(:,idx));
	falling = find(cdiff(-qr(:,idx))<0 & t > T0);
	[qrs,fdx] = unique(-qr(falling,idx));
	zsti(:,idx,2) = interp1(qrs,zst(falling(fdx),idx),qri(:,idx));
end
end
zst_hysteresis = zst_rms_i(:,:,2) - zst_rms_i(:,:,1);
zst_mid = 0.5*(zst_rms_i(:,:,2) + zst_rms_i(:,:,1));
% integrate hysteresis
I_zsth  = sum(zst_hysteresis,2)'*(x(2)-x(1));
I_zstm  = sum(zst_mid,2)'*(x(2)-x(1));

% Figure : time series + spectra

splitfigure([2,3],[1,1],fflag);
cla();
%tlim = 415+7.5-4.5 + 15*[-1,1]/4;
tlim = limits(t);
plot(t(tdx)-tlim(1),-q(tdx,k_)/Qn,'linewidth',lw);
xlabel('t/d','interpreter','latex');
ylabel('Q/Qn','interpreter','latex');
xlim(tlim-0*mid(tlim));

splitfigure([2,3],[1,2],fflag);
cla();
plot(t(tdx),q(tdx,k)-circshift(q(tdx,k),24*15),'linewidth',lw);
ylim([-50,50])

if (0)
splitfigure([2,3],[1,3],fflag);
cla();
semilogy(fx,abs(fft(q(tdx,k))))
end

splitfigure([2,3],[1,4],fflag);
cla();
plot(t(tdx)-tlim(1),[zs(tdx,k_),zsr(tdx,k_)],'linewidth',lw);
xlabel('day','interpreter','latex');
ylabel('$z_s/\mathrm{m}$','interpreter','latex');
xlim(tlim-0*mid(tlim));
%axis auto
%xlim(limits(t))

if (0)
splitfigure([2,3],[1,5],fflag);
cla();
Tt = 0.5;
Tr = 15;
S=abs(fft(zs(tdx,:))).^2;
S(S<1e-7) = NaN;
stem((fftshift(fx)-1/Tt)*Tr,fftshift(S(:,k(n)),1));
xlim(10*[-1,1]);
ylim([0,1.1]*max(S(fx>0.5./Tt,k(n))));
set(gca,'yscale','lin');
xlabel('Frequency $$\frac{\omega-\omega_1}{\omega_r}$$','interpreter','latex');
title('Spectral Power');

splitfigure([2,3],[1,6],fflag);
cla();
semilogy(fftshift(fx),fftshift(abs(fft(zs(tdx,k))),1))

end

% figure : with respect to upstream discharge

splitfigure([2,2],[3,1],fflag);
cla();
plot(qri,zst_rms_i(:,k,1),'linewidth',lw);
hold on
set(gca,'colororderindex',1);
ph = plot(qri,zst_rms_i(:,k,2),'linewidth',lw);
for idx=1:length(ph) ph(idx).HandleVisibility = 'off'; end
xlabel('Upstream river discharge $Q_r(L)/Q_n$','interpreter','latex');
ylabel('Tidal admission $|z_t(x)|/|z_t(0)|$','interpreter','latex');
lh = legend(rats(cvec(xplot/xbw)));
title(lh,'$x/L$','interperter','latex');
s = 5;
c = colororder();
clear qh;
for idx=1:length(k)
	d = 10;
	dq = cdiff(qri);
	dzst = cdiff(zst_rms_i);
	qh(idx,1) = quiver(qri(1:d:end),zst_rms_i(1:d:end,k(idx)),s*dq(1:d:end),s*dzst(1:d:end,k(idx)),0, ...
			  '-','filled','linewidth',1,'Autoscale','off','color',c(idx,:),'MaxHeadSize',1);
	qh(idx,1).HandleVisibility = 'off';
	dq   = cdiff(flipud(qri));
	dzst = cdiff(flipud(zst_rms_i(:,:,2)));
	qh(idx,2) = quiver(flipud(qri(d:d:end)),flipud(zst_rms_i(d:d:end,k(idx),2)),s*dq(d:d:end),s*dzst(d:d:end,k(idx)),0, ...
			  '-','filled','linewidth',1,'Autoscale','off','color',c(idx,:),'MaxHeadSize',1);
	qh(idx,2).HandleVisibility = 'off';
end

splitfigure([2,2],[3,2],fflag);
cla();
plot(qri,-zst_hysteresis(:,k),'linewidth',lw);
ylabel('Hystersis |z_{t,up}-z_{t,down}|/|z_t(0)|','interpreter','latex');

splitfigure([2,2],[3,3],fflag);
cla();
plot(qri,-I_zsth./I_zstm,'linewidth',lw);
ylabel('Prism $$\frac{|P\uparrow-P\downarrow|}{|P\uparrow+P\downarrow|}\,\,\,\,\,\,\,\,$$','rot',90,'interpreter','latex');
xlabel('Upstream river discharge $Q_r(L)/Q_n$','interpreter','latex');

splitfigure([2,2],[3,4],fflag);
cla();
plot(sx*inner2outer(x),(-min(zst_hysteresis,[],1)),'linewidth',lw);
hold on
% note : mean hysteresis is not well defined
%plot(sx*inner2outer(x),(-mean(zst_hysteresis)),'linewidth',lw);
ylabel('max $(|z_{t}(x)\uparrow|-|z_{t}(x)\downarrow|) / |z_t(0)|$','interpreter','latex');
xlabel('Distance from mouth $x/L$','interpreter','latex');
xlim(sx*[0,xmax]);

% figure : with respect to at-a-station discharge

splitfigure([2,3],[2,1],fflag);
cla();
plot(-qr(tdx,k),zst_rms(tdx,k),'linewidth',lw);
hold on
tdx_ = tdx([1:24:24*15,1]);
dq = cdiff(qr);
dzst = cdiff(zst_rms);
s = 10;
c = colororder();
for idx=1:length(k)
	qh(idx) = quiver(-qr(tdx_,k(idx)),zst_rms(tdx_,k(idx)),-s*dq(tdx_,k(idx)),s*dzst(tdx_,k(idx)),0,'-','filled','linewidth',1,'Autoscale','off','color',c(idx,:),'MaxHeadSize',1);
end
xlabel('At-a-station river discharge $Q_r(x)/Q_{n}$','interpreter','latex');
ylabel('Tidal admission $z_t(x)/z_t(0)$','interpreter','latex');
if (pflag)
	pdfprint(21,'img/rt-instationary-zst-vs-qr-at-a-station.pdf',ps);
end

if (1)
splitfigure([2,3],[2,2],fflag);
cla();
plot(-qr(:,k),zsr(:,k),'linewidth',lw);
end

splitfigure([2,3],[2,3],fflag);
cla();
%plot(zm(:,k),zsr(:,k));
plot(zsr(tdx,k),zst_rms(tdx,k),'linewidth',lw);


if (1)
splitfigure([2,3],[2,4],fflag);
cla();
plot(sx*(x(2:end-1)),[min(qr)',nanmean(qr)',max(qr)'],'linewidth',lw);

splitfigure([2,3],[2,5],fflag);
cla();
plot(sx*inner2outer(x),[min(zsr)',nanmean(zsr)',max(zsr)'],'linewidth',lw);
end

splitfigure([2,3],[2,6],fflag);
cla();
plot(sx*inner2outer(x),[min(zst_rms)',nanmean(zst_rms)',max(zst_rms)'],'linewidth',lw);

tdx_  = find( t>T0_ & -qr(:,k(end)) < q0 & [-qr(2:end,k(end));NaN]>q0, 1, 'first')
tdx__ = find( t>T0_ & -qr(:,k(end)) > q0 & [-qr(2:end,k(end));NaN]<q0, 1, 'first')

% figure : along channel

splitfigure([3,3],[4,1],fflag);
cla
plot(sx*inner2outer(x(1:end)),zst_rms(tdx_,:),'linewidth',lw);
hold on
plot(sx*inner2outer(x(1:end)),zst_rms(tdx__,:),'linewidth',lw);
xlim(sx*[0,xmax]);
ylabel('Admission $|z_{t}(x)|/|z_t(0)|$','interpreter','latex');
xlabel('Distance from mouth $x/L$','interpreter','latex');
legend('rising','falling');

splitfigure([3,3],[4,2],fflag);
cla
plot(sx*inner2outer(x(1:end)),zsr(tdx_,:)/hn,'Handlevisibility','off','linewidth',lw);
hold on
plot(sx*inner2outer(x(1:end)),zsr(tdx__,:)/hn,'Handlevisibility','off','linewidth',lw);
plot(sx*inner2outer(x(1:end)),zb/hn,'--k','linewidth',lw);
xlim(sx*[0,xmax]);
legend('bed level','location','southeast');
ylabel('Tidally averaged water level $z_{r}/h_n$','interpreter','latex');
xlabel('Distance from mouth $x/L$','interpreter','latex');

splitfigure([3,3],[4,3],fflag);
cla
plot(sx*(x(1:end)),kz(tdx_,:)/k0,'linewidth',lw);
hold on
plot(sx*(x(1:end)),kz(tdx__,:)/k0,'linewidth',lw);
xlim(sx*[0,xmax]);
ylabel('Damping $Im(k_z)/k_0$','interpreter','latex');
xlabel('Distance from mouth $x/L$','interpreter','latex');
%ylim([0,0.12]);

splitfigure([3,3],[4,4],fflag);
cla();
plot((sx*x(2:end-1)),-qr(tdx_,:),'linewidth',lw);
hold on
plot((sx*x(2:end-1)),-qr(tdx__,:),'linewidth',lw);
ylabel('River discharge $Q_r/Q_{n}$','interpreter','latex');
xlabel('Distance from mouth $x/L$','interpreter','latex');
xlim(sx*[0,xmax]);

splitfigure([3,3],[4,5],fflag);
cla();
plot((sx*x(2:end-1)),qt_rms(tdx_,:),'linewidth',lw);
hold on
plot((sx*x(2:end-1)),qt_rms(tdx__,:),'linewidth',lw);
ylabel('Tidal discharge $||Q_t||/Q_{n}$','interpreter','latex');
xlabel('Distance from mouth $x/L$','interpreter','latex');
xlim(sx*[0,xmax]);

splitfigure([3,3],[4,6],fflag);
cla();
plot(sx*mid(x(2:end-1)),kq(tdx_,:)/k0,'linewidth',lw);
hold on
plot(sx*mid(x(2:end-1)),kq(tdx__,:)/k0,'linewidth',lw);
xlim(sx*[0,xmax]);
ylabel('Damping $Im(k_q)/k_0$','interpreter','latex');
xlabel('Distance from mouth $x/L$','interpreter','latex');

splitfigure([3,3],[4,7],fflag);
cla();
plot(sx*inner2outer(x),-ur(tdx_,:),'linewidth',lw);
hold on
plot(sx*inner2outer(x),-ur(tdx__,:),'linewidth',lw);
ylabel('Tidally averaged velocity $u_r/u_{n}$','interpreter','latex');
xlabel('Distance from mouth $x/L$','interpreter','latex');
xlim(sx*[0,xmax]);

splitfigure([3,3],[4,8],fflag);
cla();
plot(sx*inner2outer(x),ut_rms(tdx_,:),'linewidth',lw);
hold on
plot(sx*inner2outer(x),ut_rms(tdx__,:),'linewidth',lw);
ylabel('Tidal velocity $||u_t||/u_{n}$','interpreter','latex');
xlabel('Distance from mouth $x/L$','interpreter','latex');
xlim(sx*[0,xmax]);

splitfigure([3,3],[4,9],fflag);
cla();
plot(sx*x,ku(tdx_,:)/k0,'linewidth',lw);
hold on
plot(sx*x,ku(tdx__,:)/k0,'linewidth',lw);
xlim(sx*[0,xmax]);
ylabel('Damping $Im(k_u)/k_0$','interpreter','latex');
xlabel('Distance from mouth $x/L$','interpreter','latex');

if (pflag)
	pdfprint(14,'img/rt-instationary-zs-over-time.pdf',ps);
if(0)
	pdfprint(15,'img/rt-instationary-zs-spectral-power.pdf',ps);
end
	pdfprint(31,'img/rt-instationary-zst-vs-qr-upstream.pdf',ps);
	pdfprint(33,'img/rt-instationary-hysteresis-prism.pdf',ps);
	pdfprint(34,'img/rt-instationary-hysteresis-max.pdf',ps);
	pdfprint(41,'img/rt-hysteresis-zst-along-channel.pdf',ps);
	pdfprint(42,'img/rt-hysteresis-zsr-along-channel.pdf',ps);
	pdfprint(43,'img/rt-hysteresis-kz-along-channel.pdf',ps);
	pdfprint(44,'img/rt-hysteresis-qr-along-channel.pdf',ps);
	pdfprint(45,'img/rt-hysteresis-qt-along-channel.pdf',ps);
	pdfprint(46,'img/rt-hysteresis-kq-along-channel.pdf',ps);
	pdfprint(47,'img/rt-hysteresis-ur-along-channel.pdf',ps);
	pdfprint(48,'img/rt-hysteresis-ut-along-channel.pdf',ps);
	pdfprint(49,'img/rt-hysteresis-ku-along-channel.pdf',ps);
end

