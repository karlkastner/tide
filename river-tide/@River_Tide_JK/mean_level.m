% Thu  7 Jul 16:37:34 CEST 2016
% Karl Kastner, Berlin
%
%% tidally averaged surface elevation
%% c.f. Jay and Kukulka
%
% time invariant:
% x  : distance from gauge to river mouth
% b  : width
% h0 : mean depth
% T  : major period
%
% time variant input:
% R0 : range at river mouth
% Qr : river dishcarge
%
% time invariant parameter:
% cD : drag coefficient
%
function h = mean_level(obj,x,alpha,R0,h0,b,Qr)
%function h = mean_level(obj,x,alpha,R0,cD,h0,b,Qr,T)
	if (1 == length(alpha))
		alpha(2) = alpha(1);
	end
	% Qt = obj.tidal_discharge(x,R0,h0,b,Qr,T);
	Qt = obj.tidal_discharge(x,R0,cD,h0,b,Q);
	% eq. 12 in J&K 2003b
	h = alpha(1) .* Qr.^(2/3) + 1/6*alpha(2) .*  abs(Qt).^2 ./ abs(Qr).^(4/3);
end % mean_level

