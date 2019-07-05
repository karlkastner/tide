% Sat 14 Oct 15:30:31 CEST 2017
%
%% damping modulus of the tidal wave for river flow only
%
% function [r, k, r_0, r_inf] = damping_modulus_river(Q0,W,H,cd,omega)
function [rr, ri, k, r_0, r_inf] = damping_modulus_river(Q0,W,H,cd,omega)
	if (issym(Q0) || issym(W) || issym(H) || issym(cd) || issym(omega))
		syms g
	else
		g = Constant.g;
	end

	if (0)
		az = 1;
		Qt = 0;
		[r, k, rk, r_lin] = damping_modulus(Q0,W,H,cd,omega,1,0);
	else

	q0 = Q0./W;

	% valid for small and large q0
%	r = real(sqrt(-omega.^2./(g.*H)     + 2i*omega.*cd./g.*q0./H.^3))

	re = -omega.^2./(g.*H);
	im = 2*omega.*cd./g.*q0./H.^3;
	[rr,ri]  = root_complex(re,im)

	% wave number
	k = rr+1i*ri;

	% small q0
	r_0   = (cd.*q0)./sqrt(g*H.^5);

	% large q0
	r_inf = sqrt(cd.*omega.*q0./(g*H.^3));
	end

%omega^2/(g*H)
%2*omega*cd/g*q0/H^3
%	r2 = real(r_)
%	(cd*q0*(1/(H^5*g))^(1/2))
%	r/r2
end

