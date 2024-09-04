% Fri  2 Sep 12:53:52 CEST 2022
% TODO transform r into a time
function [t,Q,qq] = generate_river_discharge(Qmu,Qsd,Qs,Ts,rf,pf,T,dt)

n = round(T/dt);

t = (0:n-1)'*T/n;
Qmu_ = Qmu + Qs*sin(2*pi*t/Ts)';
Qmu_ = cvec(Qmu_);
%p = exprnd(mu,n,1);
%p = exprnd(Qmu_);
Qsd = Qsd.*(Qmu_./Qmu);
[a,b] = gam_moment2param(Qmu_,Qsd);
p = gamrnd(a,b);
qq=p;
for idx=1:pf
	qq(:,idx+1) = filter(1-rf,[1,-rf],qq(:,idx),Qmu_(1));
end

Q = qq(:,end);

end


