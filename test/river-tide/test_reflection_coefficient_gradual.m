g = 9.81;

h    = 10;
w0   = 1000;
winf = 1; % w0/4;
omega = 2*pi/86400;

c=sqrt(g*h);
lambda = 2*pi*c/omega;

nx = 1e6;
Xi = 4*lambda*[-1 1];
x = linspace(Xi(1),Xi(2),nx)';
xs = x/lambda;

% f' = 1/sqrt(2*pi*s^2)
df = (2.^(1:16))./lambda;
s = 1./(sqrt(2*pi)*df);

figure(1);
clf
r = [];
for idx=1:length(df)
	w  = w0 + (winf-w0)*normcdf(x,0,s(idx)*abs(winf-w0)); %*0.5*(winf+w0));
	w = cvec(w);
	dw_dx = cdiff(w)./cdiff(x);
	mdw_dx(idx,1) = max(abs(1./w.*dw_dx));
	r(idx,1) = reflection_coefficient_gradual(x,h,w,omega);

if (0)	
	subplot(2,2,1);
	plot(xs,w)
	hold on

	subplot(2,2,2);
	plot(xs,-(lambda./w).*dw_dx);
	ylabel('-\lambda/w dw/dx')
	hold on
end

%	ylim([1 1e5])

end
%p = sqrt(winf/w0);
%rmax = (1-p)/(1+p)
rmax = 1;
subplot(2,2,3)
semilogx(lambda*mdw_dx,abs(r)/rmax)
%semilogx(df,abs(r)/rmax)
ylabel('r/rmax');
xlabel('dw/dx');
xlim([1 100]);


h0   = 100;
hinf = 1;
w    = 1;

%depth convergence

df = (2.^(1:20))./lambda;
s  = 1./(sqrt(2*pi)*df);
figure(2);
clf
r  = [];
mdh_dx = [];
for idx=1:length(df)
	h  = h0 + (hinf-h0)*normcdf(x,0,s(idx)*abs(hinf-h0)); %*0.5*(winf+w0));
	h = cvec(h);
	dh_dx = cdiff(h)./cdiff(x);
	mdh_dx(idx,1) = max(abs(1./h.*dh_dx));

	r(idx,1) = reflection_coefficient_gradual(x,h,w,omega);


if (0)	
	subplot(2,2,1);
	plot(xs,h)
	hold on

	subplot(2,2,2);
	
	plot(xs,-(lambda./h).*dh_dx);
	ylabel('-\lambda/w dw/dx')
	hold on
end
%	ylim([1 1e5])
end
% TODO, there is still something wrong
%rmax = -(hinf.^0.5-h0.^0.5)/(hinf.^0.5+h0.^0.5)
%rmax = 1/3;
rmax = 1;
subplot(2,2,3);
semilogx(lambda*mdh_dx,abs(r)/rmax);
%semilogx(lambda*mdh_dx,abs(r)/rmax)
%semilogx(df,abs(r)/rmax)
ylabel('r/rmax');
xlabel('dh/dx');
xlim([1 100]);


