% Mon 10 Aug 12:39:35 +08 2020
function [y] = solve(obj)
	for idx=1:length(obj.rt)
		obj.rt(idx).opt.dischargeisvariable = true;
		rt = obj.rt(idx);
		opt= rt.opt;
		odefun{idx} = @rt.odefun;
		bcfun{idx}  = @rt.bcfun;
		inifun{idx} = @opt.ifun;
		xi(idx,:)   = rt.Xi;
		nx(idx,1)   = opt.nx;
		%100; %round(abs(diff(xi(idx,:)))/dx)+1;
	end % for idx

	% TODO, awkaward to choose options of first rt-object
	[out] = bvp2c(    odefun ...
			, bcfun ...
		        , inifun ...
			, xi ...
			, nx ...
			, obj.junction_condition ...
			, true ...
			, obj.rt(1).opt ...
		     );

	% unstack solution
	y = [];
	for idx=1:length(nx)
		obj.rt(idx).out = out(idx);
		obj.rt(idx).postprocess(out(idx).x,out(idx).y,out(idx).yc);
		y = [y;out(idx).y];
	end % for idx
end % solve

