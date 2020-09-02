% Tue 21 Mar 16:17:33 CET 2017

function test_amplitude_mwl_dx_cai

zbfun = @(x) deal(-10,0);
wfun  = @(x) deal(500,0); 
Qrfun = @(x) 3e3;

eta0  = 1;
zs0   = 0;
val0  = [eta0; zs0];
Xlim = [0 1e6];

% test
x = 0;

Qr = [1 2.5 5 7.5 10]*1e3; %[1 1e3 2e3 3e3 4e3];
Qr = linspace(1e3,1e4,5)';

T = 86400;
omega  = 2*pi/T;
Uscale = 2;
cd     = 1e-3;

figure(1)
clf
for idx=1:length(Qr)
	Qrfun = @(x) Qr(idx);


	val0(2) = 10;
	[X Y] = river_tide_jay(val0,omega,cd,Uscale,Xlim,zbfun,wfun,Qrfun)
figure(1);
subplot(2,2,1)
plot(X,Y(:,1),'.-');
title('\eta Jay');
hold on
subplot(2,2,2)
plot(X,Y(:,2),'.-');
title('MWL Jay');
hold on

	val0(2) = 0;
	[X Y] = river_tide_cai(val0,omega,cd,Uscale,Xlim,zbfun,wfun,Qrfun);

subplot(2,2,3)
plot(X,Y(:,1),'.-');
title('\eta Cai');
hold on
subplot(2,2,4)
plot(X,Y(:,2),'.-');
title('MWL Cai');
hold on

end% for

end

