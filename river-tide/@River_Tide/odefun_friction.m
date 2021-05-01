% Sun 15 Mar 15:03:27 +08 2020
% Karl Kastner, Berlin
%
% TODO : change h0.^3 to h0.^p for arbitrary power (manning)
% TODO : 1/h nonlinearity !
function f = odefun_friction(obj,f,k,Q,QQ,Qhr,h0,w0,Cd,c,D1_dx)
	g     = obj.g;
	omega = obj.omega;
	pi_   = obj.pi;
	fl = obj.opt.ode.oh || obj.opt.ode.gh;

	% - d/dt Cd/(g h^3 w^2) 1/pi(f0 Q_hr^2 + f1 Q + f2 Q^2)
	% d/dt Q_hr^2 = 0, drops
	%s  = - (1i*k*omega).*Cd./(pi_.*g.*h0.^3.*w0.^2);
	s  = - (1i*k*omega).*Cd./(pi*g.*h0.^3.*w0.^2);

	% f1 : Q''
	% f2 : Q'
	f(:,2) = f(:,2) + fl.*s.*(   -3./(1i.*omega.*w0).*c(:,1).*Qhr  ...
				   + -3./(1i.*omega.*w0).*c(:,2).*Q(:,1) );

	% quadratic interaction of Q and 1./k Q (!)
	Qk      = Q./(0:(size(Q,2)-1));
	Qk(:,1) = 0;

	% TODO the cubic interaction is more complicated,
	% because Q*Q^2 terms results in the same frequency in addition to Q*Q0^2
	% for the time being, it ends all up in the rhs,
	% which converges, as the terms are small
	DQ2k1 = fourier_multiplicative_interaction_coefficients(D1_dx*Qk,Q,size(Q,2),1);

	%DQ2k1 = fourier_multiplicative_interaction_coefficients(Q,D1_dx*Qk,size(Q,2),1);
	DQ3k1 = fourier_multiplicative_interaction_coefficients(DQ2k1,Q,size(Q,2),1);

	% f3 : Q : self damping
	f(:,3) = f(:,3) + s.*(-c(:,2).*Qhr - c(:,3).*abs(Q(:,1)) );

	% f4 : constant term, damping by interaction
	% sign changed according to test case 8
	% TODO the factor in front of c(3) seems wrong,
	% for z1=0.01/h0, it is too small, for z1=0.1/h0 it is just right,
	% maybe the 1/h non-linearity cannot be ignored
	f(:,4) = f(:,4) + s.*(   c(:,3).*QQ(:,k+1) ...
			       + fl.*(  -3./(1i.*omega.*w0).*c(:,2).*DQ2k1(:,k+1) ...
                                      + -3./(1i.*omega.*w0).*c(:,3).*DQ3k1(:,k+1)));
end % odefun_friction

