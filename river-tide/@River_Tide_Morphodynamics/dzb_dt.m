% Fri  7 Aug 19:04:20 +08 2020
function dzb_dt = dzb_dt(obj,t,zb)
	p      = obj.sediment.p;
	rho    = obj.sediment.rho;
	obj.fun.zb = @zbfun;
%(x) interp1(obj.x,zb,x);

	% flow
	obj.init();
	obj.solve();
	dx     = diff(obj.x);

	% transport
	w      = obj.width(mid(obj.x));
	Qs     = obj.sediment_transport();

	% derivative of bed-level with respect to time
	if (1)
		% least stable (e~1900)
		if (1)
		dzb_dt = 1./(p.*w.*rho).*(diff(Qs)./dx);
		% upwinding
		switch (obj.norder)
		case {'leapfrog'}
			% nothing to do
		case {'upwind'}
			dzb_dt = [dzb_dt(2:end);0];
		end
		% downwinding, super unstable
		% dzb_dt = [0;dzb_dt(1:end-1)];
		else
			% beam warming flux
			% TODO, here the dx^2 term is missing
			dQs = [0.5*(3*Qs(1:end-2) - 4*Qs(2:end-1) + Qs(3:end));
                              Qs(end-1)-Qs(end)];
			dzb_dt = 1./(p*w*rho).*dQs./dx;
		end
	else
		Qsc = obj.sediment_transport(true);
		% e ~ 150
		dzb_dt = 2./(p.*w.*rho).*((Qs(2:end)-Qsc)./dx);
		% e ~ 15, why is this closest to stable? should be more unstable
		%dzb_dt = 2./(p.*w.*rho).*((Qsc-Qs(1:end-1))./dx);
	end

	%dzb_dt = csmooth(dzb_dt);
	%dzb_dt = inner2outer(dzb_dt);
	dzb_dt(end) = 0;

function zbi = zbfun(x)
	zbi = interp1(mid(obj.x),zb,x,'linear','extrap');
end

end
