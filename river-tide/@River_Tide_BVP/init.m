% Sun  8 Oct 13:08:39 CEST 2017
% Karl Kastner, Berlin
%
%% initial condition
%
%% function obj = init(obj)
function init(obj)
	obj.clear();

	obj.neq = sum(obj.opt.oflag)+1;

	obj.bc_transformation();

	obj.hydrosolver.jfun    = obj.junction_condition;
	obj.hydrosolver.odefun  = @obj.odefun;
	obj.hydrosolver.bcfun   = @obj.bcfun;
	obj.hydrosolver.inifun  = @obj.initial_value;
	% TODO make a field odeopt
	obj.hydrosolver.opt     = obj.opt;

	if (obj.nc > 1)
		obj.hydrosolver.opt.dischargeisvariable = true;
	end

	obj.hydrosolver.init();

end % River_Tide_BVP/init

