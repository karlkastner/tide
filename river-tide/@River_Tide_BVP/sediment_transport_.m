% Fri  7 Aug 19:07:20 +08 2020
% Karl Kastner, Berlin
%
%% compute sediment transport for a single channel
%
function [Qs, rt, asym, stokes] = sediment_transport_(obj,cdx,t,ddir)
	d_mm = obj.sediment.d_mm;

	% discharge and channel properties at segment centres
	w   = obj.width(cdx);
	cd  = obj.cd(cdx);
	Cz  = drag2chezy(cd);
	z1  = obj.z(1,cdx);
	h0  = obj.h0(cdx);

	% at end-points
	if (1) %nargin()>3 && ddir)
		Q   = obj.out(cdx).Qc;
		w   = mid(w);
		% TODO hc from zc, to avoid recursive inner2outer
		h0  = mid(h0);
		Cz  = mid(Cz);
		z1  = mid(z1);
	else
		Q   = obj.out(cdx).Q;
	end

	U   = Q./(h0.*w);
	Ur = abs(U(:,2))./abs(U(:,1));
	zr = z1./h0;

	% sediment transport
	Qs0 = total_transport_engelund_hansen(Cz,d_mm,U(:,1),[],w);
	[scale, rt, asym, stokes] = river_tide_transport_scale(Ur,zr,5,obj.opt.stokes_order);

	Qs = scale.*Qs0;

	% buffer for values on boundary
	Qs = [0;Qs;0];

	% left boundary condition
	p       = obj.bc_Qs(1,cdx).p;
	val     = obj.bc_Qs(1,cdx).rhs;
	if (~isempty(val))
	if (ddir ~= 0)
		Qs(1)   = p*val + (1-p)*(1.5*Qs(2)-0.5*Qs(3));
	else
		Qs(1)   = p*val + (1-p)*(2*Qs(2)-Qs(3));
	end
	end

	% right boundary condition
	p       = obj.bc_Qs(2,cdx).p;
	val     = obj.bc_Qs(2,cdx).rhs;
	if (~isempty(val))
	if (ddir ~= 0)
		Qs(end) = p*val + (1-p)*(1.5*Qs(end-1)-0.5*Qs(end-2));
	else
		Qs(end) = p*val + (1-p)*(2*Qs(end-1)-Qs(end-2));
	end
	end

end % River_Tide_BVP_/sediment_transport

