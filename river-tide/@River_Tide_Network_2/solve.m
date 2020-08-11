% Mon 10 Aug 12:39:35 +08 2020
function [y_C,y] = solve(obj)
%	nx = obj.nx;
	for idx=1:length(obj.rt)
		odefun{idx} = @(x,y) obj.rt(idx).odefun(x,y);
		bcfun{idx}  = @(x,void,ccdx) obj.rt(idx).bcfun(x,void,ccdx);
		xi(idx,:)   = obj.rt(idx).Xi;
		nx(idx)     = 100; %round(abs(diff(xi(idx,:)))/dx)+1;
	end
	[x_C,y_C,out] = bvp2c_stacked(odefun,bcfun,xi,nx); %,varagin)

	% unstack solution
	for idx=1:length(nx)
		obj.rt(idx).postprocess(x_C{idx},y_C{idx});
	end
end

