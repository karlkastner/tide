	% estimate qt
	% TODO this is a poor estimate, ignoring frictional damping and bed slope
	if (nargin()<7 || isempty(Qt))
		qt = az1.*sqrt(g.*h);
		Qt = qt*w;
	else
		qt = Qt/w;
	end

