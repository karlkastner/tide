%TODO derive with fourier interaction coefficients

%Q1^(5-k) z1

l=linspace(-pi,pi,21);
g=10;
qs=[];
qs_=[];
qs__=[];

h0 = 10;
%Q0_ = h0*[0.125,0.25,0.5,1];
Q0_ = h0*[0.25,0.5,1]; %[0.125,0.25,0.5,1];
for jdx=1:length(Q0_)
 for idx=1:length(l);
e=1;
 t = innerspace(0,1,100)';
 Q0 = Q0_(jdx);
aQ1 = e*sqrt(g*h)
 Q1 = aQ1*cos(2*pi*t);
 z=e*cos(2*pi*t+l(idx));
e  = e*exp(1i*l(idx))
z_ = e*exp(2i*pi*t);
%z_ = e*cos(2*pi*t);
%abs([z,z_])
%[z-0.5*(z_+conj(z_))]
%pause
 Q  = Q0+Q1;
Qs0 = mean((Q0./h).^5);
qs(idx,jdx)=mean([(Q./(h+z)).^5])./Qs0;
%qs_(idx,jdx) = mean([(Q./h).^5.*(1-5*z./h)])./Qs0;
%qs_(idx,jdx) = mean([(Q0.^5 + 5*Q0.^4.*Q1 + 10*Q0.^3.*Q1.^2 + 10*Q0.^2*Q1.^3 + 5*Q0*Q1.^4 + Q1.^5)./h.^5])./Qs0;
n = 12;
qs_(idx,jdx)   = mean([(Q0.^5 + 5*Q0.^4.*Q1 + 10*Q0.^3.*Q1.^2 + 10*Q0.^2*Q1.^3 + 5*Q0*Q1.^4 + Q1.^5)./h.^5.*(1-5*z./h+0*15*(z./h).^2)])./Qs0;
n = 18;
qs__(idx,jdx)  = mean([(Q0.^5 + 5*Q0.^4.*Q1 + 10*Q0.^3.*Q1.^2 + 10*Q0.^2*Q1.^3 + 5*Q0*Q1.^4 + Q1.^5)./h.^5.*(1-5*z./h+15*(z./h).^2)])./Qs0;

%qr(idx,jdx) = stokes_transport_(aQ1./Q0,e./h0,0)	
%qr_(idx,jdx) = stokes_transport_(aQ1./Q0,e./h0,1)	
%qr__(idx,jdx) = stokes_transport_(aQ1./Q0,e./h0,2)	
%t = linspace(0,1,n);
%Q1 = e*sqrt(g*h)*cos(2*pi*t);
%z  = e*cos(2*pi*t+l(idx));
%Q  = Q0+Q1;
%qs__(idx,jdx)=mean([(Q./(h+z)).^5])./Qs0;

%qs__(idx,jdx) = mean([(Q0.^5 + 5*Q0.^4.*Q1 + 10*Q0.^3.*Q1.^2 + 10*Q0.^2*Q1.^3 + 5*Q0*Q1.^4 + 0*Q1.^5)./h.^5.*(1-5*z./h)])./Qs0;
%qs_(idx,jdx) = mean([(Q0.^5 + 10*Q0.^3.*Q1.^2 + 5*Q0.*Q1.^4)./h.^5.*(1-5*z./h)])./Qs0;
%qs__(idx,jdx) = mean([(Q./h).^5.*(1-5*z./h+15*(z./h).^2)])./Qs0;
 end;
end
 qs;
clf
 plot(l/(2*pi),qs)
hold on
set(gca,'colororderindex',1)
plot(l/(2*pi),qs_,'--')
set(gca,'colororderindex',1)
plot(l/(2*pi),qs__,'-.')
rms([qs-qs_])
rms([qs-qs__])
