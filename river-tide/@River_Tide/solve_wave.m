% Wed 11 Oct 10:18:54 CEST 2017
%% solve for the oscillatory (tidal) componets
%%
%% function obj   = solve_wave(obj)
function obj   = solve_wave(obj)

	% solve the system of ode's
	[x, y, obj.out] = feval( obj.opt.solver, ...
			         @obj.odefun, ...
			         @obj.bcfun, ...
				 obj.Xi, ...
				 obj.opt);

	obj.postprocess(x,y);

end % River_Tide/solve_wave

