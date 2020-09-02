% Tue  1 May 11:09:29 CEST 2018

n = 1024^2;
x = 2*pi*(0:n-1).'/n;

omax = 4;
a = [ rand(1,omax+1)].';
b = [ 0, rand(1,omax)].';
m = 1; % only 0,1,2 as input, no terdiurnal tide
a(end-m:end) = 0;
b(end-m:end) = 0;

c = a+1i*b;


% mean diurnal, semidiurnal
cp = fourier_power_exp(c);

A = ones(n,1);
A_ = A;
c_ = c(1);
for odx=1:omax
	A(:,odx+1)    = exp(1i*odx*x);
	A_(:,2*odx)   = 1/2*exp(1i*odx*x);
	A_(:,2*odx+1) = 1/2*exp(-1i*odx*x);
	c_ = [c_; c(odx+1); conj(c(odx+1))];
end
%y1_ = A*c;
y1  = real(A*c);
y1_ = A_*c_;
'y1-y1_'
disp(norm(y1-y1_)/norm(y1))


for pdx=1:3

yp  = y1.^pdx;
%c
cpp  = A \ yp;
cpp_ = A_\ yp;
cpp__ = [cpp_(1);cpp_(2:2:end-1)];;
cpp__(2:end) = 1/2*(cpp__(2:end) + conj(cpp_(3:2:end)));
disp('yp - A*c, not correct, as real part was taken')
disp(norm(real(yp - A*cpp))/norm(yp))
disp('yp - A*c_, correct, re+im')
disp(norm(real(yp - A_*cpp_))/norm(yp))
disp('yp - yp__, correct, re reapplied)')
disp(norm(real(yp - A*cpp__))/norm(yp))
cpp = cpp__;
if (1==pdx)
disp('c-cpp__')
disp(norm(c-cpp__)/norm(c))
end

d = cp(1:omax+1,pdx)-cpp(1:omax+1);

disp(pdx)
disp(['analytic, numerical, a-n'])
disp([cp(1:omax+1,pdx), cpp, d])

if (norm(d) > sqrt(eps)*norm(cpp))
	warning('test failed')
end

end

% c0=1; c1=1/3+1/5*i; y=@(x) real(c0+c1*exp(1i*x)); quad( @(x) y(x).*cos(x),-pi,pi)/pi, quad( @(x) y(x).*sin(x),-pi,pi)/pi, quad( @(x) y(x),-pi,pi)/pi
if (0)
	a = 2; b = 3; c=a+1i*b; [cp,op]=complex_exp_product([c,c],[1,1]), fourier_power_exp([0,c,0])  
end

