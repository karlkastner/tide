% Fri 17 Nov 09:32:17 CET 2017
%% converts damping rate to drag coefficient
%% c.f. friedrichs, ippen harleman
function cd = rdamping_to_cdrag_tide(omega,r,h0,az1)
	full = true;
	g  = Constant.gravity; % is eliminated
	c  = sqrt(g*h0);
	k  = omega/c;
	u0 = az1*sqrt(g/h0);
	if (full)
		r_ = -(2*omega*r)/(k*(r^2/k^2 - 1));
		cd = -(3*h0*pi*r_)/(8*u0);
	else
		% linearised
		cd = (3*pi*h0^2*r)/(4*a);
	end
end % rdamping_to_cdrag_tide

