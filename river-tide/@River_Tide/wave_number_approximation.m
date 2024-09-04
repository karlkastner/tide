% Fri  5 Jan 12:21:09 CET 2018
% Karl Kastner, Berlin
%
%% approximate wave number of the left and right traveling wave for variable coefficients
%%
%% TODO merge with wave_number_analytic
%%
%% function [k, k0, dk0_dx_rel, obj] = wave_numer_aproximation(obj)
%
function [k10, kz1, kq1, dk, obj] = wave_numer_aproximation(obj,x,Q,h0,z0,w0,Cd,dw_dx,D1_dx)
%	Q1    = obj.Q(1);
%	x     = obj.x;
	L     = x(end);

	% coefficients of the characteristic polynomial of the wave equation
	% TODO pass z0 when h is not precomputed
	c           = obj.odefun(x,Q);
	c           = c(:,:,1);
	% roots of the characteristic polynomial
	r           =  roots2(c(:,1:3));
	% wave number in case of constant coefficients
	% ode:                        c2 q'' + c1 q' + c0 q = 0 
	% with q = hat q e(i(ot-kx)): c2 (-ik)^2 + c1 ik + c0 = 0
	% characteristic poly:        c2 r^2 + c1 r + c0   = 0
	%                          => r = -ik <=> k = ir
	k10         = 1i*r;

	% relative change of coefficients along channel
	dr_dx     = derivative1(x,r);
	dr_dx_rel = bsxfun(@times,1./((r(:,1)-r(:,2))),dr_dx);

	dk = 1i*dr_dx_rel;

	% wave number of upstream and downstream traveling waves
	kq1 = [ k10(:,1) - dk(:,1), ...
                k10(:,2) - dk(:,2)];
	kz1 = [ k10(:,1) + dk(:,1), ...
                k10(:,2) + dk(:,2)];

end % River_Tide/wave_number_approximation

