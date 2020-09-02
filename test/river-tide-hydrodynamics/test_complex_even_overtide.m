% Mon 23 Apr 20:16:26 CEST 2018

a1=2;b1=3;
a2=5;b2=7;
o  = 2*pi;
n1 = 2;
n2 = 3;

% a2=a1;b2=b1;n2=n1;

x = linspace(0,1)';

c1 = a1+1i*b1;
e1 = exp(i*n1*o*x);
c2 = a2+1i*b2;
e2 = exp(i*n2*o*x);

subplot(2,2,1)
plot(x,[real(c1*e1),1/2*(c1*e1+conj(c1*e1))])

subplot(2,2,2)
plot(x,[real(c2*e2),1/2*(c2*e2+conj(c2*e2))])

subplot(2,2,3)
y1 = real(c1*e1).*real(c2*e2);
y2 = 	1/4*( c1*c2*e1.*e2 ...
		+ c1*e1.*conj(c2*e2) ...
		+ conj(c1*e1).*(c2*e2) ...
		+ conj(c1*c2*e1.*e2) ...
	     );
y2 = 	1/4*(     c1*c2*e1.*e2 ...
		+ c1*conj(c2)*e1.*conj(e2) ...
		+ conj(c1)*c2*conj(e1).*e2 ...
		+ conj(c1*c2)*conj(e1.*e2) ...
	     );
y2 = 	1/4*(     c1*c2*(e1.*e2) ...
		+ c1*conj(c2)*e1.*conj(e2) ...
		+ conj(c1)*c2*conj(e1).*e2 ...
		+ conj(c1*c2)*(conj(e1.*e2)) ...
	     );
y2 = 	1/4*(     c1*c2*(e1.*e2) ...
		+ c1*conj(c2)*e1.*conj(e2) ...
		+ conj(c1)*c2*conj(e1).*e2 ...
		+ conj(c1*c2*e1.*e2) ...
	     );
y2 = 	1/2*(     real(c1*c2*e1.*e2) ...
		+ real(c1*conj(c2)*e1.*conj(e2)) ...
	     );
y2 = 	1/2*(     real(c1*c2*exp(i*(n1+n2)*o*x)) ...
		+ real(c1*conj(c2)*exp(i*(n1-n2)*o*x)) ...
	     );
y2 = 	1/2*(     real(c1*c2*exp(i*(n1+n2)*o*x)) ...
		+ real(c1*conj(c2)*exp(i*(n1-n2)*o*x)) );

y2 = 	1/2*(     real(c1*c2*exp(i*(n1+n2)*o*x)) ...
		+ real(conj(c1)*c2*exp(i*(n2-n1)*o*x)) );
if (0)
y2 = 	1/2*(     real(c1^2*exp(i*(n1+n2)*o*x)) ...
		+ abs(c1)^2 );
end
norm(y1-y2)
plot(x,[y1,y2])



if (0)
q1  = a1*cos(2*pi*x)+b1*sin(2*pi*x);
q2  = a2*cos(4*pi*x)+b2*sin(4*pi*x);
q1_ = sqrt(a1^2+b1^2)*exp(2i*pi*x - i*atan2(b1,a1));
q2_ = sqrt(a2^2+b2^2)*exp(4i*pi*x - i*atan2(b2,a2));

subplot(2,2,1);
plot(x,[q1,real(q1_)]);

subplot(2,2,2);
plot(x,[q2,real(q2_)]);

subplot(2,2,3);
plot(x,[q1.*q1 real(q1_.*q1_) real(q1_).*real(q1_) ])

subplot(2,2,4);
plot(x,[q1.*q2 real(q1_.*q2_)   real(q1_).*real(q2_) ])

 x=2*pi*linspace(0,1)'; a1=2; a2=5; p1=3; p2=2; e1=exp(i*p1*x);  e2=exp(p2*i*x); y1 = a1*e1; y2=a2*e2; plot([real(y1).*real(y2), 1/2*(a1*a2*real(exp(i*(p1-p2)*x) + exp(i*(p1+p2)*x)))])
end


