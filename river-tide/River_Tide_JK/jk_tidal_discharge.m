% Thu  7 Jul 16:37:13 CEST 2016
% Karl Kastner, Berlin
%
%% tidal discharge
%% c.f. Jay and Kukulka
%
function Qt = jk_tidal_discharge(x,R0,cD,h0,b,Qr,T)
	omega = 2*pi./T;
	r     = jk_damping_modulus(cD,h0,b,Qr,T);
	R     = jk_tidal_range(x,R0,cD,h0,b,Qr,T);  
	Qt    = omega*b./(r*sqrt(8)).*R;
end

