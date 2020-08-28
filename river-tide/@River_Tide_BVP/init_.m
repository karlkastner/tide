% Sun  8 Oct 13:08:39 CEST 2017
%% provide initial condition by solving the backwater equation for surface level
%% TODO this should not be solved as a ivp but included in the bvp iteration
%% TODO generate the mesh here and precompute fixed values instead of passing functions
%% TODO Q0 should not be a function
%% function obj = init(obj, Xi)
function init(obj)
%	for idx=1:obj.nc
%		obj.init_(cdx);
%	end
	obj.clear();

	obj.neq = sum(obj.opt.oflag)+1;
	obj.nc = 

	%for idx=1:obj.nc
		%obj.rt(idx).opt.dischargeisvariable = true;
		%xi(idx,:)   = obj.rt(idx).opt.Xi;
		%nx(idx,1)   = obj.rt(idx).opt.nx;
	%end % for idx

	obj.bc_transformation();

	obj.odesolver.jfun    = obj.junction_condition;
	obj.odesolver.odefun  = @obj.odefun;
	obj.odesolver.bcfun   = @obj.bcfun;
	obj.odesolver.inifun  = @obj.fun.initial_value;
	% TODO make a field odeopt
	obj.odesolver.opt     = obj.opt;
%	obj.odesolver.xi      = xi;
%	obj.odesolver.nx      = nx;
	if (obj.nc > 1)
		obj.odesolver.opt.dischargeisvariable = true;
	end

	obj.odesolver.init();
%end
%function obj = init(obj)
%	for cdx=1:obj.nc
%	switch (obj.opt.hmode)
%	case {'matrix'}
%	otherwise
%		obj.neq = sum(obj.opt.oflag);
%	end
	% initial condition function
	%if (isfield(obj.opt,'initial_value'))
	%	obj.fun.init = @obj.opt.initial_value; 
	%else
	%	obj.fun.init = @obj.initial_value;
	%end
%	end % for cdx
end % River_Tide_BVP/init

