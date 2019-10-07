% Thu  7 Jul 16:34:56 CEST 2016
%% predict tidal range
% function R = tidal_range(x,R0,cD,h0,b,Qr,T)
function [R, obj] = tidal_range(obj,x,R0,h0,b,Qr)
	%r = obj.damping_modulus(cD,h0,b,Qr,T);
	r = obj.damping_modulus(h0,b,Qr);
	R = R0.*exp(r.*x);
end % River_Tide_JK/tidal_range

