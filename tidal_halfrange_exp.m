% Mon  4 Nov 19:13:02 +08 2019
%
% half and mid range for a fourier series with mean and 2 frequency components
%
% z(t) = c(1) + Re(c(2)*exp(2i*pi*t) + c(3)*exp(4i*pi*t))
%
% function [hrange,mid,y_lim,r] = tidal_range_exp(c)
function [hrange,mid,y_lim,r] = tidal_halfrange_exp(c)
	[r,y]  = tide_low_high_exp(c);
	y      = real(y);
	y_lim  = [min(y,[],2),max(y,[],2)];
	hrange = 0.5*(y_lim(:,2)-y_lim(:,1));
	mid    = 0.5*(y_lim(:,1)+y_lim(:,2));
	hrange(isnan(hrange)) = 0;
	mid(isnan(mid)) = 0;
end

