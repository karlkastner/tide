% Fri 17 Feb 21:30:18 CET 2017
%% extract constituents of tpxo object
function [obj2, t, obj] = get_constituents(obj,name)
	id        = obj.build_index();
	n         = length(name);
	t         = struct();
	t.tidecon = NaN(n,4);
	t.freq    = NaN(n,1);
	t.name    = repmat(' ',n,4);
	for jdx=1:n
		% look up if field exists
		if (isfield(id,name{jdx}))
			% get index of field
			id_       = id.(name{jdx});
			% copy field values
			t.tidecon(jdx,:) = obj.t.tidecon(id_,:);
			t.name(jdx)    = obj.t.name(id_);
			t.freq(jdx)    = obj.t.freq(id_);
		else
			warning(['Const ',name{jdx},'not contained in input']);
		end % if
	end % for jdx
	obj2 = T_Tide(t);
end % get_constituent

