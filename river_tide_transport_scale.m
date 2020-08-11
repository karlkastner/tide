% Thu 30 Apr 12:13:24 +08 2020
% for p = 3, rt-transport exceeds r-transport when ut ~ 0.82 ur 
% for p = 5, rt-transport exceeds r-transport when ut ~ 0.43 ur 
function scale = transport_scale_river_tide(ut_div_ur,p)
	% Qs = Qsr*scale
	switch (p)
	case {1}
		scale = 1;
	case {2}
		scale = 1+1/2*(ut_div_ur).^2;
	case {3}
		% scale ~ 2.5 when ut/ur = 1
		% scale ~ 2   when ut/ur ~ 0.81
		% scale ~ 1 + 3/2 (ut/ur)^2 when ut/ur small (exact)
		scale = 1+3/2*(ut_div_ur).^2;
	case {4}
		scale = 1 + 3*(ut_div_ur).^2 + 3/8*(ut_div_ur).^4;
	case {5}
		% scale ~ 7.9 when ut_div_ur = 1
		% scale = 2   when ut_div_ur ~ 0.43
		% scale ~ 1 + 5*(ut/ur)^2 when ut/ur small
		scale = (1  + 5*(ut_div_ur).^2 + 15/8*(ut_div_ur).^4 );
	end
end

% constant term of sin^2n(x)
function c = sin2n_constant(n)
	c = 1./2.^(2*n)*nchoosek(2*n,n);
end

function derive()
	syms t ur ut; q=expand((ur+ut*sin(t)).^5), combine(expand(combine(q,'sincos')),'sincos')
end

