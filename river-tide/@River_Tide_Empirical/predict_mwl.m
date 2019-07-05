% Mon 11 Jul 15:43:41 CEST 2016
% Karl Kastner, Berlin
%
%% predict the mean water level
function [h, level, obj] = predict_mwl(obj,R0,Qr)
	Qr = abs(Qr);
	c = obj.c.mwl;
	if (~obj.nlflag)
		n  = length(Qr);
		X  = [ones(n,1),Qr.^(2/3), abs(R0).^2./(Qr.^(4/3))];	
		h  = X*c;
	else
		h = c(1) + c(2)*(Qr.^2 + c(3)*abs(R0).^2).^(1/3);
	end
	level = h-c(1);
end

