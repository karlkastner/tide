
Q = Q0 + Q1*exp(2*pi*omega*t) + conj(Q1*exp(2*pi*omega*t))

% C0 = 1/(2*pi) int_-pi^pi |Q|Q dt
% C1 = 1/(2*pi) int_-pi^pi |Q|Q * exp(2piot) dt
% dronkers has then c1 minus as well !
% Lorents (c.f. Dronkers) just determines the timing of slack water
% and then intergrates the three parts with theire respective sign:
C0 = 1/(2*pi)*(   int(Q^2,-pi,omega*t2) ...
		- int(Q^2,omega*t1,omega*t2) ...
		+ int(Q^2,omega*t1,pi) )

% test
Q=linspace(0,1,10)';
 t=innerspace(0,1,10)';
 Q=sin(2*pi*t)+0.25;
 mQ=max(Q);
 [t.^0, sin(2*pi*t), cos(2*i*t) sin(4*pi*t) cos(4*pi*t)]'*([abs(Q).*Q, mQ.^2/pi*(c(2)*Q/mQ+c(4).*(Q/mQ).^3)])

