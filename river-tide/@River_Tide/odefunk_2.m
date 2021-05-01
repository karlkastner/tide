% Fri 16 Oct 17:29:15 +08 2020

function [f] = odefunk_2(obj, x, Q, h0, z0, zb, w0, Cd, dw_dx, D1_dx)

	% allocate memory
	f = zeros(nx,4,nf+1);

	Qmid = Q(:,1);
	Qhr  = sum(abs(Q(:,2:end)),2);

	cf = -obj.friction_coefficient(Qmid,Qhr);
	
	% powers
	A2 = fourier_multiplicative_interaction_coefficients(A,A,nf+1);
	A3 = fourier_multiplicative_interaction_coefficients(A,A2,nf+1);
	Q2 = fourier_multiplicative_interaction_coefficients(Q,Q,nf+1);
	Q3 = fourier_multiplicative_interaction_coefficients(Q,Q2,nf+1);

	% zero derivative
	%dQ = Q;

	% first derivative in space
	% TODO get dQ/dx directly from
	%dQ(:,:,2) = D1*Q;

	% second derivative in space
	%dQ(:,:,3) = D2*Q;

	zs_x = D1*zs;
	w_x  = D1*w;

	% nope, not for A0!
	%A_x = 1./(1i*o).*dQ(:,:,3);
	% non-zero frequency terms
	%zs_x = (1./w).*A_x - (1./w.^2).*w_x.*A + zb_x;

	% antiderivative in time
	%zs_tx = -1/w.^2*w_x*dQ(:,:,2) + 1./w.*dQ(:,:,3);
	%zs_x  = 
	% zero frequency term
	%zs_x(:,1) = z0_x;

	Q_t    = (1i*omega*(0:nf)).*Q;
	AQ_t   = fourier_multiplicative_interaction_coefficients(A,Q_t);
	QQ_t   = fourier_multiplicative_interaction_coefficients(Q,Q_t);
	A2zs_x = fourier_multiplicative_interaction_coefficients(A2,zs_x);

	% frequency components of ode-coefficients w/o advective acceleration term
	% Q_xx
	oc(:,:,3) =  -g*w.*A3;
	% Q_x
	oc(:,:,2) =   g*w_x.*A3 - 3*g*w.^2.*A2zs_x - 2*w.^3.*AQ_t;
	% Q
	oc(:,:,1) =  (   w.^3.*A2.*d2_dt2 ...
		   + Cd.*w.^4.*(    f1.*d_dt ...
		                + 2*f2.*Q_t ...
		                + 3*f3.*QQ_t ...
                               ) ...
		 );

	% coefficients of advective acceleration term
 	if (obj.flag.advective_acceleration)
		Q_x  = D1*Q;
		AQ   = fourier_multiplicative_interaction_coefficients(A,Q);
		Q_Qx = fourier_multiplicative_interaction_coefficients(Q,Qx);

		% Q_xx
		oc{3} = oc{3} + w.^3.*Q2;
		% Q_x
		oc{2} = oc{2} + 2*w.^3.*(-Q_Qx + AQ_t + AQ*d_dt );
		% Q
		oc{1} = oc{1} - 2*w.^3*A_x.*Q_t;
	end % if advective acceleration

	f = coefficient_frequency_components(oc,dQ);

	f(:,:,1) = obj.odefunQ0(x, h0, z1, zb, w0, dw_dx, Q(:,1), Qhr, Q(:,2:end), Cd, cf);
end % odefunk_2

