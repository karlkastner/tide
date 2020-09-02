% Mon 10 Aug 12:39:35 +08 2020
% Karl Kastner, Berlin
%
%% determine hydrodynamics
%
function [y] = solve(obj)

	obj.hydrosolver.solve();

	% unstack solution for channels
	y = [];
	for cdx=1:obj.nc
		obj.postprocess(cdx, ...
				obj.hydrosolver.out(cdx).x, ...
				obj.hydrosolver.out(cdx).y, ...
				obj.hydrosolver.out(cdx).yc);
		%y = [y; obj.hydrosolver(cdx).out.y];
	end % for cdx
end % solve

