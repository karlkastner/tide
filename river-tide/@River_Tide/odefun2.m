% Sun 22 Apr 11:36:45 CEST 2018
%% coefficients of the ordinary differential quation of the even overtide
function f  = odefun2(Q0, Qhr, Q1, Q2, h0, dh0_dx, dz0_dx, w, dw_dx, cd, g, c, o1, flag)
	% specific discharges
	qhr = Qhr./w;
	q0  = Q0./w;
	q1  = Q1./w;
	q2  = Q2./w;
		
	c0  = c(:,1);
	c1  = c(:,2);
	c2  = c(:,3);

	% allocate memory
	n = max(length(Q0),length(Q1));
	if (~issym(Q1))
		f = zeros(n,4);
		Pi = pi;
	else
		Pi = sym('pi');
	end

	% q''
	f(:,1) =   g*(-1/(1i*o1))*ones(n,1);	     % change of surface elevation (z1' ~ Q1'')
	% q'
	f(:,2) =   g./(2i*o1*w).*dw_dx;	             % change of width
	% q
	f(:,3) = ( (4i*o1)./h0                   ... % wave travelling
	 	 + (2*cd.*c1.*qhr)./(Pi*h0.^3)   ... % self damping
		 + (4*cd.*c2.*q0 )./(Pi*h0.^3)   ... % damping by river flow
		 );
	% constant (generating term)
	% the lsh is mulitplied by Q2, not q2, so rhs has to be scaled by w
	f(:,4)  = -w.*(c2.*cd.*q1.^2)./(Pi*h0.^3); ... % generation by main component
end % River_Tide/odefun2

