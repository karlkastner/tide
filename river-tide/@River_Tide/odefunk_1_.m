% Sun  8 Oct 13:08:39 CEST 2017
% Karl Kastner, Berlin
%
%%
%% coefficients of the ordinary differential equation of the k-th freqeuncy
%% component of the tide
%%
%% f1 Q'' + f2 Q' + f3 Q + f4 = 0
%%
%% function [f, F3]  = odefunk(obj, k, Q, QQ, Qhr, h0, dh0_dx, dz0_dx, w0, dw0_dx, Cd, c, D1_dx)
%  TODO : precomputed derivatives in the rhs should be avoided,
%	  instead the block-diagonal entries should be set accordingly
function [f, F3]  = odefunk_1_(obj, k, Q, QQ, Qhr, h0, dh0_dx, dz0_dx, w0, dw0_dx, Cd, c, D1_dx)

	n = size(Q,1);
	f = zeros(n,4);
	if (obj.issym)
		f = sym(f);
	end

	% jacobian without advective acceleration
	% homogeneous (forcing-free) part of the SWE
	f = obj.odefun_swe_jacobian(f,k,h0,w0);

	% advective acceleration
	% TODO, should  be part of the SWE-jacobian
	if (obj.opt.ode.advective_acceleration)
		f = obj.odefun_advective_acceleration(f,k,Q,QQ,h0,dh0_dx,w0,dw0_dx,D1_dx);
	end

	% forcing by friction
	f = obj.odefun_friction(f,k,Q,QQ,Qhr,h0,w0,Cd,c,D1_dx);

	% forcing by along-channel change of width
	f = obj.odefun_width(f,k,w0,dw0_dx);

	% note : forcing by dzb/dx term is only non-zero for the mean (z0)

	% normalize
	f = f./(-1i.*obj.omega(k).*obj.g.^-1.*w0.^-1);
end % odefunk

