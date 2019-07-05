% Thu  1 Jun 17:44:36 CEST 2017
%% select a subsect of constituents
function obj = select(obj,fdx)
	obj.t.tidecon = obj.t.tidecon(fdx,:);
	obj.t.name    = obj.t.name(fdx,:);
	obj.t.freq    = obj.t.freq(fdx,:);
end

