% Wed 31 May 16:13:05 CEST 2017
%% raleigh criterion for resolving tidal constituents
%% T > 1/|f1-f2|
function [Tmin, T_, dT] = raleigh_criterion(varargin)
	n = length(varargin);
	T = zeros(n,1);
	tc = tidal_constituents();
	% fetch periods
	for idx=1:length(varargin)
		if (~isstr(varargin{idx}))
			T(idx) = varargin{idx};
		else
			try
				T(idx) = tc.(lower(varargin{idx}));
			catch e
				disp(e);
				error('here');
			end
		end % if
	end % for
	
	F    = 1./T;
	% get minimum difference in frequency
	dF = abs(bsxfun(@minus,F,F'));
	dT = 1./dF;
	% all unique pairs
	[dF sdx] = sort(dF(triu(true(size(dF)),+1)));
	% time required to get exactly one period difference
	T_   = 1./dF;
	Tmin = T_(1);
end

