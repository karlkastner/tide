% Mon 10 Aug 12:39:35 +08 2020
function [y] = solve(obj)

	obj.rt(1).odesolver.solve();

	% unstack solution for channels
	y = [];
	for cdx=1:length(obj.rt)
		obj.rt(cdx).out = obj.rt(1).odesolver.out(cdx);
		obj.rt(cdx).postprocess(obj.rt(cdx).out.x, ...
					obj.rt(cdx).out.y, ...
					obj.rt(cdx).out.yc);
		y = [y; obj.rt(cdx).out.y];
	end % for cdx
end % solve

