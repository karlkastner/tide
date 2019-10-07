% Thu  7 Jul 16:37:13 CEST 2016
% Karl Kastner, Berlin
%
%% tidal discharge
%% c.f. Jay and Kukulka
%
%% function Qt = tidal_discharge(obj,x,R0,h0,b,Qr)
% old: function Qt = tidal_discharge(obj,x,R0,cD,h0,b,Qr,T)
function [Qt,obj] = tidal_discharge(obj,x,R0,h0,b,Qr)
	omega = 2*pi./T;
	%r     = obj.damping_modulus(cD,h0,b,Qr,T);
	%R     = obj.tidal_range(x,R0,cD,h0,b,Qr,T);  
	r     = obj.damping_modulus(h0,b,Qr);
	R     = obj.tidal_range(x,R0,h0,b,Qr);  
	Qt    = omega*b./(r*sqrt(8)).*R;
end % River_Tide_JK/tidal_discharge

