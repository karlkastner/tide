% Sun  8 Oct 13:08:39 CEST 2017
%%
%% coefficients of the differential equation of the main tidal species
%%
%% f1 Q'' + f2 Q' + f3 Q + f4 = 0
%%
%% TODO rename f into c
%% TODO better pass dzb_dx instead of dz0_dx
%% TODO aa, oh and gh terms are not tested for width ~= 1
%%
% function [f F3]  = odefun1(Q0, Qhr, Q1, h0, dh0_dx, dz0_dx, w, dw_dx, cD, g, c, o1, flag)
function [f, F3]  = odefun1(Q0, Qhr, Q1, Q2, h0, dh0_dx, dz0_dx, w, dw_dx, cD, g, c, o1, flag)
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
	f(:,1) =   (1i*g)/o1; %*ones(n,1) ...		% change of surface elevation (z1' ~ Q1'')
		 - flag.aa*(1i*q0.^2)./(h0.^3*o1);
	% q'
	f(:,2) = (-1i*g./(o1*w).*dw_dx ...		% change of width
		 + flag.gh*(1i*g*dz0_dx)./(h0*o1) ...
		 + flag.aa*(2*q0)./h0.^2 ...
		 + flag.oh*( (3i*c2.*cD.*q0.^2) ...
			    +(3i*cD.*q0.*c1.*qhr) ...
			    +(3i*cD.*c0.*qhr.^2))./(h0.^4*o1*Pi) ...
		 + flag.oh*flag.aa*(2*1i*q0.^2.*dh0_dx)./(h0.^4*o1) ...
		 + flag.oh*flag.gh*(+( (3*1i*c0.*cD.*qhr.^2)...
			              +(3i*c1.*cD.*q0.*qhr) ...
                                      +(3i*c2.*cD.*q0.^2))./(h0.^3*o1*Pi) ...
                                    -( (2i*c0.*cD.*qhr.^2) ...
				      +(2i*c1.*cD.*q0.*qhr) ...
                                      +(2i*c2.*cD.*q0.^2) )./(h0.^4*o1*Pi)) ...
		 );
	% q
	f(:,3) = ( (1i*o1)./h0 ...			% wave propagation
		 + (  cD.*c1.*qhr)./(Pi*h0.^3) ...	% self damping
		 + (2*cD.*c2.*q0 )./(Pi*h0.^3) ...	% damping by river flow
		 - flag.aa.*(2*q0.*dh0_dx)./h0.^3 ...   % advective accleration
		 );
	if (nargout() > 1)
		F3 = [(1i*o1)./h0, (cD.*c1.*qhr)./(h0.^3*Pi), (2*cD.*q0.*c2)./(h0.^3*Pi)];
	end

	% inhomogeneous part
	% lhs is multiplied by Q1 not q1, so rhs must be scaled up by w
	% interaction with D2 : e1*e2 = a e1 + b e3
	f(:,4) = -w.*c2.*cD.*conj(q1).*q2./(Pi*h0.^3);
end % River_Tide/odefun1

