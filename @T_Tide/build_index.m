% 2017-02-17 20:57:10.318780258 +0100
%% build a structure whose field names contain the index
function [id, obj] = build_index(obj)
	name = obj.name;
	id   = struct();
	for idx=1:size(name,1)
		str = lower(chomp(name{idx})); %(idx,:)));
		if (~isempty(str) && str(1) ~= '2' && str(1) ~= '3')
			id.(str) = idx;
		end % if
	end % for idx
	obj.id = id;
end % build index
