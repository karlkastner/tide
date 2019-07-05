% Wed 11 Oct 10:18:54 CEST 2017
%% call stationary or non-stationary solver respectively
function obj   = solve(obj)
	switch (obj.opt.model_str)
	case {'wave'}
		obj.solve_wave();
	case {'swe'}
		obj.solve_swe();
	end
end

