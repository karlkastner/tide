% Sun  8 Oct 13:08:39 CEST 2017
%% provide initial condition by solving the backwater equation for surface level
%% TODO this should not be solved as a ivp but included in the bvp iteration
%% TODO generate the mesh here and precompute fixed values instead of passing functions
%% TODO Q0 should not be a function
%% function obj = init(obj, Xi)
function obj = init(obj)
	obj.clear();

	switch (obj.opt.hmode)
	case {'matrix'}
		obj.neq = sum(obj.opt.oflag)+1;
	otherwise
		obj.neq = sum(obj.opt.oflag);
	end

	obj.bc_transformation();

	% initial condition function
	if (isempty(obj.opt.ifun))
	switch (obj.opt.hmode)
	case {'matrix'}
		obj.opt.ifun = @obj.initial_value; 
	otherwise
		% use default setting of solver
	end
	end

end % River_Tide/init

