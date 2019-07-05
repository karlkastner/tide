% Thu  1 Jun 12:26:36 CEST 2017
%% order constituents as specified by "name"
function obj = reorder(obj,name)
	id  = obj.build_index();
	n   = length(obj.name);
	m   = length(name);
	sdx = zeros(n,1);
	for idx=1:m %length(name)
		if (isfield(id,lower(name{idx})))
			sdx(id.(lower(name{idx}))) = idx;
		end % field id
	end % for idx
	flag      = (0 == sdx);
	sdx_      = m + cumsum(flag);
	sdx(flag) = sdx_(flag);
	mm = max(max(sdx),m);
	% TODO warn unused
	freq            = obj.t.freq;
	obj.t.freq      = NaN(mm,1);
	obj.t.freq(sdx) = freq;
	name            = obj.t.name;
	obj.t.name      = repmat(' ',mm,4);
	obj.t.name(sdx,:) = name;
	tc                   = obj.t.tidecon;
	obj.t.tidecon        = NaN(mm,4);
	obj.t.tidecon(sdx,:) = tc;

end % reorder

