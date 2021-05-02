% Sun  8 Oct 13:08:39 CEST 2017
% Karl Kastner, Berlin
%
%
%% coefficients of the wave equation for river-tides decomposed in frequency components
%% zero frequency component corresponds to backwater equation with tidal influence
%
% TODO heed identity abs(Qi)^2 = abs(Qi^2) = conj(Qi)*Qi
%      to move coupling terms from rhs (b) to lhs (A)

%        f = obj.odefun(x, [Q0, Qt], zs, zb, w0, Cd, dw_dx, D1_dx, D2_dx2);
function [f, obj] = odefun(obj, x, Q, zs, zb, w0, Cd, dw_dx, D1_dx, D2_dx2)
	nx  = length(x);

	neq = size(Q,2);

	% tidally averaged depth
	h0 = (zs(:,1)-zb);
	h0 = max(h0,obj.opt.hmin);

	Qmid   = Q(:,1);
	Qhr    = sum(abs(Q(:,2:end)),2);

	cf     = -2*obj.friction_coefficient(Qmid,Qhr);
	
	if (nargin()<8)
		dw_dx = derivative1(x,w0);
	end

	if (nargin()>8)	
%		dw_dx  = D1_dx*w0;
		dh0_dx = D1_dx*h0;
		dz0_dx = D1_dx*zs(:,1);
	else
		dz0_dx = derivative1(x,zs(:,1));

		D1_dx  = derivative_matrix_1_1d(x,[]);
		%dh0_dx = derivative1(x,h0);
		dh0_dx = D1_dx*h0;
%		D1_dx = derivative_matrix[];
		D2_dx2 = derivative_matrix_2_1d(x,[]);
	end

	f = zeros(nx,4,neq);

	z1 = obj.discharge2level(x,Q(:,2),w0);

	QQ = fourier_quadratic_interaction_coefficients(Q,size(Q,2),2);

	switch (obj.opt.odefunk)
	case {'odefunk_3'}
	f = obj.odefunk_3( ...
				Q, QQ, Qhr, zs, zb, h0, dh0_dx, dz0_dx, ...
					     w0, dw_dx, Cd, cf, D1_dx, D2_dx2);
	otherwise
	f = obj.odefunk_1( ...
				Q, QQ, Qhr, zs, zb, h0, dh0_dx, dz0_dx, ...
					     w0, dw_dx, Cd, cf, D1_dx, D2_dx2);
	end

	if (~obj.opt.dischargeisvariable)
		f(:,:,1)  = obj.odefunz0(x, h0, z1, zb, w0, dw_dx, Q(:,1), ...
						      Qhr, Q(:,2:end), Cd, cf);
	else
		f(:,:,1)  = obj.odefunQ0(x, h0, z1, zb, w0, dw_dx, Q(:,1), ...
						      Qhr, Q(:,2:end), Cd, cf);
	end
end % River_Tide/odefun

