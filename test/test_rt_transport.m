% Sun  5 Jul 11:28:16 +08 2020
% determine magnitude
% there is reverse transport, as long as u0 < 0.2 u_t
p = [3,5];

a1 = 1;
a0 = linspace(-1,1,10)/2;
b1 = 0;
a2 = 0
b2 = linspace(-1,1,3)/2; 

f = @(t,a0,a1,a2,b2,p) (a0 + a1*sin(2*pi*t) + a2*sin(4*pi*t) + b2*cos(4*pi*t)).^p;

clf
for pdx=1:length(p)
	s = [];
	s_ = [];
for idx=1:length(a0)
	for jdx=1:length(b2)
		s(idx,jdx) = quad(@(t) f(t,a0(idx),a1,a2,b2(jdx),p(pdx)),0,1);
		s_(idx,jdx) = rt_transport(a0(idx),a1,a2,b2(jdx),p(pdx)); 
	end
end
subplot(2,2,pdx)
plot(a0,s);
hold on
set(gca,'colororderindex',1)
plot(a0,s_,'--')
grid on

end
	
%t = linspace(0,1);
%fut = @(t,u0,p) u0+(1-p)*cos(2*pi*t) + p*cos(4*pi*t); plot(t,ut); 
%p0 = fminsearch(@(p) -mean(fut(t,0,p).^ps), 0.5)
%u0=fzero(@(u0) mean(fut(t,u0,p0).^ps),0)
%ut_range = max(fut(t,0,p0)) - min(fut(t,0,p0))
%u0*ut_range
%end

