% Thu  7 Jul 19:19:57 CEST 2016
%% Regression of tidal coefficients according to Jay & Kulkulka
%%
%% coefficients of the r-regression factor 2 apart for specis (jay C7)
%% this can be repeated for each tidal species (diurnal, semidiurnal)
function [d,obj] = rivertide_regress(obj,z,z0,phi,phi0,h,Ur)
	m   = size(z,2);
	Ur  = abs(Ur);

	% log magnitude ratio
	Z = log(abs(z./z0));

	% phase difference
	delta_phi = phi - phi0;
	delta_phi = phasewrap(delta_phi);

	d = [];

	% for each species
	%d   = jk_rivertide_regress(z(fdx),z0(fdx),phi(fdx),phi0(fdx),h(fdx),Ur(fdx));
	for idx=1:m
		fdx = isfinite(Ur) & isfinite(z(:,idx)) & isfinite(z0(:,idx));
		n   = sum(fdx);

		% regression matrix
		% h(t,Qr) is dependent on Qr
		X = [ones(n,1),Ur(fdx),1./sqrt(Ur(fdx)).*abs(z0(fdx,idx)./h(fdx)).^2];


		% coefficients
		% solve by QR factorisation (this system can be ill conditioned)
		[Q R] = qr(X,0);
		d(:,:,idx) = R \ (Q'*[Z(fdx,idx), delta_phi(fdx,idx)]);
	end

	obj.d = d;
end % River_Tide_JK/rivertide_regress

