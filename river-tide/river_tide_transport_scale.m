% Thu 30 Apr 12:13:24 +08 2020
%
% scale of river transport to total transport that accounts for river-tide
% interaction and stokes transport
%
% input : ut_div_ur = Ut/U0 = Qt/Q0	ratio of tidal to river flow
%         Zr = zt/h0 		ratio or tidal surface amplitude to tidally averaged water depth
%	  p : degree of nonlinearity of transport
%	  (typically 1 for wash-load, 3 for bed-load, 5 for suspended bed-material load)
%	  
%
% output :
%	scale  : total scale
%	rt     : contribution by river-tide interaction (Q_0^(p-2k)*Q_t^(2k))
%	asym   : contribution by tidal asymmetry        (Q0^(p-3k)*Q1^(2k)*Q2^k)
%	stokes : contribution by Stokes transport       (1/h-nonlinearity)
%
% The leading term of rt is always of order (Ut/U0)^2
% Ratio ut/u0 required for rt-transport to equal river-transport (column 2)
% and rt-transport when ut = u0 (column 3)
%     | Qs_rt = Qs_0 | ut/u0 = 1
%  p  | ut/u0        | Qs_tot/Qs_0
%  3  | 0.82         | 2.5
%  5  | 0.43         | 7.9
%
function [scale, rt, asym, stokes] = river_tide_transport_scale(ut_div_ur,Zr,p,order)
	if (nargin()<2)
		Zr = zeros(1,size(ut_div_ur,2));
	end
	if (nargin()<3)
		p = 5;
	end
	if (nargin()<4)
		order = 2;
	end
	switch (p)
	case {1}
		scale = 1;
	case {2}
		scale = 1 + 1/2*(ut_div_ur).^2;
	case {3}
		% scale ~ 2.5 when ut/ur = 1
		% scale ~ 2   when ut/ur ~ 0.81
		% scale ~ 1 + 3/2 (ut/ur)^2 when ut/ur small (exact)
		scale = 1 + 3/2*(ut_div_ur).^2;
	case {4}
		scale = 1 + 3*(ut_div_ur).^2 + 3/8*(ut_div_ur).^4;
	case {5}
		% scale ~ 7.9 when ut_div_ur = 1
		% scale = 2   when ut_div_ur ~ 0.43
		% scale ~ 1 + 5*(ut/ur)^2 when ut/ur small
		% scale = (1  + 5*(ut_div_ur).^2 + 15/8*(ut_div_ur).^4 );
		[scale, rt, asym, stokes] = river_tide_transport_scale_5(ut_div_ur,Zr,order);
	end
end % river_tide_transport_scale

