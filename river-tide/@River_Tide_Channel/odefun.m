% Sun  8 Oct 13:08:39 CEST 2017
% Karl Kastner, Berlin
%
%% coefficients of the backwater and wave equation for river-tides
%
% TODO precompute Cd, h, w, dhdx, dwdx,
%      this requires the bvp solve to accept a predefined mesh as an argument
% TODO heed identity abs(Qi)^2 = abs(Qi^2) = conj(Qi)*Qi
%      to move coupling terms from rhs (b) to lhs (A)
function [f, obj] = odefun(obj,x,y)
	% without any input, order of coupled odes is returned
	if ( nargin()<3 || isempty(x) )
		f(1) = 1;
		%if (obj.opt.dischargeisvariable)
		%	f(2) = -1;
		%	k = 2;
		%else
			k = 1;
		%end
		for idx=1:length(obj.rt.opt.oflag)
			if (obj.rt.opt.oflag(idx))
				k    = k+1;
				f(k) = 2;
			end
		end
		return;
	end

%	k = obj.neq;

	g      = obj.rt.g;
	omega1 = obj.rt.omega;
	%flag   = obj.flag;

	w0     = obj.width(x);
	D1_dx  = obj.D1_dx(x);
	D2_dx  = obj.D2_dx2(x);
	dw_dx  = D1_dx*w0;

	zb     = obj.zb(x);
	nx     = length(x);

	[z0,Q0,Qt] = obj.extract(x,y);
	Q0         = repmat(Q0,nx,1);
	zs         = [z0, obj.rt.discharge2level(x,Qt,w0)];

	% TODO properly determine range and midrange
	%Qhr    = sum(abs(Qt),2);
	%Qmid   = Q0;
	%[Qhr,Qmid] = tidal_range_exp([Q0,Qt]);

%	h0     = z0 - zb;
%	TODO hmin limitation?
	Cd     = obj.cd(x,z0-zb);

        %f = odefun@River_Tide(obj, x, [Q0, Qt], h0, zs, z0, zb, w0, Cd, dw_dx, D1_dx, D2_dx);
        %f = odefun@River_Tide(obj.rt, x, [Q0, Qt], h0, zs, z0, zb, w0, Cd, dw_dx, D1_dx, D2_dx);
        %f = obj.rt.odefun(x, [Q0, Qt], h0, zs, z0, zb, w0, Cd, dw_dx, D1_dx, D2_dx);
        f = obj.rt.odefun(x, [Q0, Qt], zs, zb, w0, Cd, dw_dx, D1_dx, D2_dx);
end % River_Tide_Channel/odefun

