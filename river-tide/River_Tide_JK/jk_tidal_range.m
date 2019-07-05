% Thu  7 Jul 16:34:56 CEST 2016
%% predict tidal range
function R = jk_range(x,R0,cD,h0,b,Qr,T)
	r = jk_damping_modulus(cD,h0,b,Qr,T);
	R = R0.*exp(r.*x);
end

