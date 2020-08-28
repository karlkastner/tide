% Mon 10 Aug 12:39:35 +08 2020
function [y] = solve(obj)

%	obj.rt(1).odesolver.solve();
	obj.odesolver.solve();

	% unstack solution for channels
	y = [];
	for cdx=1:length(obj.rt)
		%obj.rt(cdx).out =  obj.rt(1).odesolver.out(cdx);
		%out = obj.rt(cdx).out;
		%obj.rt(cdx).postprocess(1, ...
		obj.postprocess(cdx, ...
				obj.odesolver(cdx).out.x, ...
				obj.odesolver(cdx).out.y, ...
				obj.odesolver(cdx).out.yc);
		y = [y; obj.odesolver(cdx).out.y];
	end % for cdx
end % solve

