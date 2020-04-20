% 2017-03-23 13:05:38.400718902 +0100
nt = 1e3;
t = (0:nt-1)'/nt;

yl = 2.1*[-1 1];

order = 5;

a     = [1 1];
phi   = [0 pi/3];
omega = 2*pi*[1 2];
dp    = (phi(2)-2*phi(1));

a = a(1); phi=phi(1); omega=omega(1);

Ur = [1 0.5 0.25 0];
l = (sum(abs(a))+max(abs(Ur)))^2;

clf();
for idx=1:length(Ur)
	
	u              = sum(bsxfun(@times,a,sin(bsxfun(@plus,t*omega,phi))),2) + Ur(idx);
	uau	       = -u.*abs(u);

	% decompose
	[tc uau_c] = stft(uau,t(2)-t(1),1,[1 1/2])
	uau_c      = uau_c;

	[uau_d void p] = friction_dronkers(u,[],[],order);
	[tc uau_t] = stft(uau_d,t(2)-t(1),1,[1 1/2]);
uau_t
	
%	[uau_t void void p_] = friction_trigonometric_dronkers([Ur(idx) a 0],dp,[],[],order);
	
	%uu = friction_dronkers(u,Ur(idx),1);
	%uu = friction_dronkers(u,1,midrange(u));
	
	subplot(4,3,3*idx-2)
	plot(t,[u uau]);
	hold on
	plot(t,uau_d,'--');
	ylim([-l l]);
	hline(0);
	
	subplot(4,3,3*idx-1)
	bar(1:5,[[uau_c], [uau_t(1:5)]]);
	title('Fourier decomposition u|u|');
	legend('numeric','analytic');
	ylim(yl);
	
	subplot(4,3,3*idx-0)
	bar(p)
	title('Coefficients');

end % for
