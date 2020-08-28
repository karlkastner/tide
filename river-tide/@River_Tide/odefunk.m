% Sun  8 Oct 13:08:39 CEST 2017
% Karl Kastner, Berlin
%
%%
%% coefficients of the ordinary differential equation of the k-th freqeuncy
%% component of the tide
%%
%% f1 Q'' + f2 Q' + f3 Q + f4 = 0
%%
%% TODO rename f into c
%% TODO better pass dzb_dx instead of dz0_dx
%% TODO aa, oh and gh terms are not tested for width ~= 1
%%
% function [f, F3]  = odefunk(obj,k,Q, Qhr, h0, dh0_dx, dz0_dx, w, dw_dx, cD, c)
function [f, F3]  = odefunk(obj, k, Q, QQ, Qhr, h0, dh0_dx, dz0_dx, w0, dw0_dx, cD, c, D1_dx)

	n = size(Q,1);
	f = zeros(n,4);
	if (issym(Q))
		f = sym(f);
	end

	% jacobian without advective acceleration
	f = obj.odefun_swe_jacobian(f,k,h0,w0);

	% advective acceleration
	if (obj.flag.aa)
		f = obj.odefun_advective_acceleration(f,k,Q,QQ,h0,dh0_dx,w0,dw0_dx,D1_dx);
	end

	% forcing by friction
	f = obj.odefun_friction(f,k,Q,QQ,Qhr,h0,w0,cD,c,D1_dx);

	% forcing by width variation
	f = obj.odefun_width(f,k,w0,dw0_dx);

	% note : forcing by dzb/dx term is only non-zero for the mean (z0)

	% normalize
	f = f./(-1i.*k.*obj.omega.*obj.g.^-1.*w0.^-1);
end

