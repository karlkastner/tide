% Mon 10 Aug 12:39:35 +08 2020
% Karl Kastner, Berlin
%
%% determine hydrodynamics
%
function y = solve(obj)

	obj.hydrosolver.solve();

	% unstack solution for channels
	y = [];
	cflag = zeros(obj.nc,1);
	for cdx=1:obj.nc
		cflag(cdx) = obj.channel(cdx).postprocess( ...
				obj.hydrosolver.out(cdx).x, ...
				obj.hydrosolver.out(cdx).y, ...
				obj.hydrosolver.out(cdx).yc);
	end % for cdx

	cflag = min(cflag);
end % solve

