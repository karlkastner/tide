% Fri 28 Aug 20:37:21 +08 2020
%
% scale of transport by mean flow only to total transport
% for the case, that the transport scales with u^5
% accounts for river tide interaction and stokes transport
% Stokes transport is approximated by a Taylor expansion
%
% input :
%   Qr : Q1/Q0
%   zr : z1/h0
%
% output:
% 	scale : scale of river to total transport
%	rt    : contribution by river-tide interaction
%	asym  : contribution by tidal-asymmetry
%	sto   : contribution by stokes transpot
%
% the leading omitted term of the series are -35*(z1/h0)^3 + 70*(z1/h0^4),
% error for order of terms:
%    e \ o	0	 1	 2	    3	      4 
%    0.1000   -0.3791    0.1209   -0.0291    0.0059   -0.0011
%    0.2000   -0.5981    0.4019   -0.1981    0.0819   -0.0301
%    0.3000   -0.7307    0.7693   -0.5807    0.3643   -0.2027
%    0.4000   -0.8141    1.1859   -1.2141    1.0259   -0.7661
function [scale,rt,asym,stokes] = river_tide_transport_scale_5(Qtr,ztr,order)
	if (nargin()<3)
		order = 2;
	end
	scale = 1;
	for idx=1:size(Qtr,2)
		scale = scale + rtts5(Qtr(:,idx),ztr(:,idx));
	end
	% TODO account for interaction between species
	asym = 0;

function scale = rtts5(Qr,zr)
	aQr = abs(Qr);

	% contribution by river tides (Q0*Q1)
	rt    =     (5*aQr.^2 + 15/8*aQr.^4);

	% contribution by stokes transport (1/h)
	fr = real(zr.*Qr./aQr);
	fi = imag(zr.*Qr./aQr);
	stokes = 0;
	if (order > 0)
		% pairings : 1, 3, 5
		%scale = scale - 5*( 5/16*Qr^1 + 3/8*10*Qr^3 + 1/2*Qr.^5)*f;
		% 0  1  2  3 4  5
		% 5  4  3  2 1  0 
		% 1  5 10 10 5  1
		stokes = stokes - 5/16*( 5*8*aQr.^1 + 10*6*aQr.^3 + 1*5*aQr.^5 ).*fr;
	    if (order > 1)
		% pairings : 0, 2, 4
		% This is more complicated, has also contribution, when out of phase
		stokes = stokes + 15/16*(8*1 + 10*6*aQr.^2 + 5*5*aQr.^4).*fr.^2;
		stokes = stokes + 15/16*(8*1 + 10*2*aQr.^2 + 5*1*aQr.^4).*fi.^2;
	    end
	end
	% "1" is for the transport by river flow
	scale = rt + stokes;
end

end % river_tide_transport_scale_5

