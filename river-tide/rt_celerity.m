% Mon  4 Dec 10:56:55 CET 2017
%% celerity of the tidal wave
function [c, c0] = celerity(Q0,W,H,cd,omega,az1,S)
	if (~issym(Q0))
		g = Constant.gravity;
	else
		syms g
	end
	C = sqrt(g/cd);
	if (isempty(H))
		H = normal_flow_depth(Q0,W,C,S);
	end
	%r = damping_modulus(Q0,W,H,cd,omega,az1);
	r = damping_modulus_river(Q0,W,H,cd,omega);
	
	c0 = sqrt(g*H);
	k0 = omega./c0;
	c  = c0.*k0./sqrt(k0.^2 + r.^2);
end

