% Sun  8 Oct 13:08:39 CEST 2017
% Karl Kastner, Berlin
%
%
%% coefficients of the wave equation for river-tides decomposed in frequency components
%% zero frequency component corresponds to backwater equation with tidal influence
%
% TODO heed identity abs(Qi)^2 = abs(Qi^2) = conj(Qi)*Qi
%      to move coupling terms from rhs (b) to lhs (A)
function [f, obj] = odefun(obj, x, Q, h0, zs, z0, zb, w0, Cd, dw_dx, D1_dx, D2_dx)
	if (min(h0)<=0)
		warning('negative water depth')
	end

	nx = length(x);

	Qmid   = Q(:,1);
	Qhr    = sum(abs(Q(:,2:end)),2);

	cf     = -2*obj.friction_coefficient(Qmid,Qhr);
	
	dw_dx  = D1_dx*w0;
	dh0_dx = D1_dx*h0;
	dz0_dx = D1_dx*z0;

	f = zeros(nx,4,obj.neq);

	z1 = obj.discharge2level(x,Q(:,2),w0);

	QQ = fourier_quadratic_interaction_coefficients(Q,size(Q,2),2);

	switch (obj.opt.odefunk)
	case {'odefunk_3'}
	f = obj.odefunk_3( ...
				Q, QQ, Qhr, zs, zb, h0, dh0_dx, dz0_dx, ...
					     w0, dw_dx, Cd, cf, D1_dx, D2_dx);
	otherwise
	f = obj.odefunk_1( ...
				Q, QQ, Qhr, zs, zb, h0, dh0_dx, dz0_dx, ...
					     w0, dw_dx, Cd, cf, D1_dx, D2_dx);
	end

	if (~obj.opt.dischargeisvariable)
		f(:,:,1)  = obj.odefunz0(x, h0, z1, zb, w0, dw_dx, Q(:,1), ...
						      Qhr, Q(:,2:end), Cd, cf);
	else
		f(:,:,1)  = obj.odefunQ0(x, h0, z1, zb, w0, dw_dx, Q(:,1), ...
						      Qhr, Q(:,2:end), Cd, cf);
	end
end % River_Tide/odefun

