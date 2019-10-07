syms r x dhdx h0 h;
 
 %r = sqrt(1/h^3);
 %R_ = -3/2*1/r*diff(r,h)*dhdx; %diff(h,x)

h0 = 1;
dhdx = 1;
h = 1;

h = h0+dhdx*x;
% h = h0*exp(-*x); 

 %r = %sqrt(1/h^3);
 r = 1i/sqrt(h);
 R = 1/r*diff(r,x);
 R = simplify(expand(R),'ignoreanalyticconstraints',true)

syms x1 % positive
Ir = int(r,x);
Ir = subs(Ir,x,x1) - subs(Ir,x,0)
Ir = subs(Ir,x1,x)

w = int(R*exp(-2*Ir),x)
pause


 f.r=matlabFunction(r,'vars',{'dhdx','h0','x'});
 f.R=matlabFunction(R,'vars',{'dhdx','h0','x'});
 %f.R_=matlabFunction(Dk_);
 f.h=matlabFunction(h,'vars',{'dhdx','h0','x'});

 h0=1;
 dhdx=-1;
 nx = 1e2;
 x=linspace(0,1,nx)';

 r=f.r(dhdx,h0,x);
 h=f.h(dhdx,h0,x);

 dk_ = [];
% dk_(:,1) =  f.R(dhdx,h0,x);
% dk_(:,2) =  f.R_(dhdx,h0,x);
 dk_(:,2) =  1./r.*cdiff(r)./cdiff(x);
% dk_(:,3) = -3/2*1./h*dhdx;

 plot(x,[dk_,-3/2*1./h*dhdx]);
 %ylim([0,100])     

 
