t=linspace(0,1)';
 u0 = 0;
 u = u0(1)+sin(2*pi*t);
 v=1+u0(1);
 m=2*mean(u)/range(u);
 p=rt.friction_coefficient_dronkers(0);
 U=range(u);
 u=u/range(u);
 U_=1;
p(5)=0;
 plot([u,u.*abs(u), -(p(1)+p(2)*u + p(3)*u.^2/U_ + p(4)*u.^3/U_.^2 )/pi])

