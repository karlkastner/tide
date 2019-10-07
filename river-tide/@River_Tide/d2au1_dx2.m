% Tue 17 Apr 14:25:43 CEST 2018
%% second derivative of the tidal velocity magnitude
%%
%% note: this is for finding zeros,
%%       the true derivative has to be scaled  up by z
function d2au1_dx2 = d2au1_dx2(obj,Q0,W,h,cd,omega,az1,Qt,dh_dx,dw_dx);
	% TODO
	d2h_dx2 = 0;
	% d2z0_dx2 = -3*cd/g*Q0^2./(W^2.*h.^4);

	[k, kq, kz] = obj.wave_number_analytic(
	dkq_dx    = obj.dkq_dx(Q0,W,h,cd,omega,az1,Qt,dh_dx,dw_dx);
	%2az1_dx2 = imag(dkz_dz) + imag(kz).^2;
	d2au1_dx2 = ( ...
			- 1./h.^3.*d2h_dx2 ...
			- 2./h.^2.*imag(kq) ...
			+ 1./h.*(imag(dkq_dx) + imag(kq).^2) ...
		    );
end % River_Tide/d2au1_dx2

