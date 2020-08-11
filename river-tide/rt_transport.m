% Sun  5 Jul 11:14:23 +08 2020
% 
% transport by river tide with even overtide,
% ignoring stokes transport
%
% s = int_0^1 (a_s w u.^p_s) dt
%
% approximate zeros:
% transport direction, for a1=0, b1 =0, a2 = 0, b2 = 0.5, p = 3
%	-a0 < 2.5 b2 : upstream (low flow)
%	-a0 = 2.5 b2 : no transport
%	-a0 > 2.5 b2 : downstream (high flow)
%
% a0 : amplitude of mean flow
% a1 : amplitude of main species (drops out, eliminated by shifting t0
% a1 : sin-amplitude of even overtide 
% b2 : cos amplitude of even overtide

% TODO shift time to make b2 zero
%
function s = rt_transport(a0,a1,a2,b2,p)
	switch (p)
	case {3}
		%s = a0^3 + 3*a0^2*b2 + 3*a0*b2^2 + b2^3;
		s = a0^3 + (3*a0*a1^2)/2 + (3*a0*a2^2)/2 + (3*a0*b2^2)/2 - (3*a1^2*b2)/4;
	case {5}
%		s = a0^5 + 5*a0^4*b2 + 10*a0^3*b2^2 + 10*a0^2*b2^3 + 5*a0*b2^4 + b2^5;
		s = a0^5 + 5*a0^3*a1^2 + 5*a0^3*a2^2 + 5*a0^3*b2^2 - (15*a0^2*a1^2*b2)/2 ...
		    + (15*a0*a1^4)/8 + (15*a0*a1^2*a2^2)/2 + (15*a0*a1^2*b2^2)/2 ...
		    + (15*a0*a2^4)/8 + (15*a0*a2^2*b2^2)/4 + (15*a0*b2^4)/8 ...
		    - (5*a1^4*b2)/4 - (15*a1^2*a2^2*b2)/8 - (15*a1^2*b2^3)/8;
	end
end

