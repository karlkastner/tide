% Mon  4 Nov 19:13:02 +08 2019
% tidal half range and mid-water for 2 frequency components and a non-zero mean
function [hrange,mid,y_lim,r] = tidal_range_2(c)
	[r,y] = tide_how_high_2(c);
	y_lim  = [min(y,[],2),max(y,[],2)];
	hrange = 0.5*(y_lim(:,2)-y_lim(:,1));
	mid    = 0.5*(y_lim(:,1)+y_lim(:,2));
end

