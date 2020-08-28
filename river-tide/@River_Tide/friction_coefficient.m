% Sun  8 Oct 13:08:39 CEST 2017
% Karl Kastner, Berlin
%
function cf = friction_coefficient(obj,Qmid,Qhr)
	switch (obj.opt.friction_method)
	case {'neglect-river'} % lorentz
		% identical to Dronkers for Q0=0
		cf = obj.friction_coefficient_lorentz(abs(Qmid)./Qhr);
	case {'dronkers'}
		% half range of tidal amplitude
		cf   = -obj.friction_coefficient_dronkers(cvec(abs(Qmid)./Qhr));
	case {'godin'}
		cf   =  obj.friction_coefficient_godin(abs(Qmid)./Qhr);
	case {'savenije'}
		error('TODO');
	case {'no-reverse'}
		% identical to Dronkers for Q0>Qt
		c = [0,0,-pi];
	otherwise 
		error('objfun');
	end
end

