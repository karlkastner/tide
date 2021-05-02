% Sun  8 Oct 13:08:39 CEST 2017
% Karl Kastner, Berlin
%
%% initial condition
%
%% function obj = init(obj)
function init(obj)
	obj.clear();
	if (isempty(obj.channel))
		error('network has no channels');
	end

	obj.hydrosolver.jfun    = obj.junction_condition;
	obj.hydrosolver.odefun  = @obj.odefun;
	obj.hydrosolver.bcfun   = @obj.bcfun;
	obj.hydrosolver.inifun  = @obj.initial_value;
	obj.hydrosolver.nx      = obj.channel.nx;
	hydrosolver.opt.xs      = obj.rt.opt.xs;
	
	Xi = [];
	for cdx=1:obj.nc
		Xi(cdx,:) = [0,obj.channel(cdx).L];
	end
	obj.hydrosolver.xi = Xi;

	% TODO make a field odeopt
	obj.hydrosolver.opt     = obj.rt.opt;
	if (obj.nc > 1)
		obj.hydrosolver.opt.dischargeisvariable = true;
	end

	for idx=1:length(obj.channel)
		obj.channel(idx).rt          = obj.rt;
	end

	obj.transform_bc();

	obj.hydrosolver.init();
end % River_Tide_BVP/init

