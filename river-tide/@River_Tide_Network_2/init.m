% 2020-08-20 16:47:03.708634914 +0800

function init(obj)
	for idx=1:length(obj.rt)
		obj.rt(idx).init();
	end

	for idx=1:length(obj.rt)
		obj.rt(idx).opt.dischargeisvariable = true;
		xi(idx,:)   = obj.rt(idx).Xi;
		nx(idx,1)   = obj.rt(idx).opt.nx;
	end % for idx

	% TODO, awkaward to choose solver of first rt-object
	obj.rt(1).odesolver.jfun    = obj.junction_condition;
	obj.rt(1).odesolver.odefun  = @obj.odefun;
	obj.rt(1).odesolver.bcfun   = @obj.bcfun;
	obj.rt(1).odesolver.inifun  = @obj.inifun;
	obj.rt(1).odesolver.opt     = obj.rt(1).opt;
	obj.rt(1).odesolver.xi	    = xi;
	obj.rt(1).odesolver.nx	    = nx;
	obj.rt(1).odesolver.opt.dischargeisvariable = true;

	obj.rt(1).odesolver.init();
end

