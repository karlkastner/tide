% Tue 17 Apr 14:28:26 CEST 2018
%% second derivative of the tidal surface elevation
%%
%% note: this is for finding zeros,
%%       the true derivative has to be scaled  up by z
function d2az1_dx2  = d2az1_dx2(obj,Q0,W,h,cd,omega,az1,Qt,dh_dx,dw_dx)
	[k, kq, kz] = out.wave_number_analytic(Q0,W,h,cd,omega,az1,Qt,dh_dx,dw_dx);
	dkz_dx = out.dkz_dx(Q0,W,h,cd,omega,az1,Qt,dh_dx,dw_dx);

	d2az1_dx2 = imag(dkz_dz) + imag(kz).^2;
end

