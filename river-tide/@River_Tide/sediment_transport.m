% Fri  7 Aug 19:07:20 +08 2020
% Karl Kastner, Berlin
function [Qs,Qs0] = sediment_transport(obj,t,iscentral)
	d_mm = obj.sediment.d_mm;

	% discharge and channel properties at segment centres
	w   = obj.width;
	h   = obj.h0;	
	cd  = obj.cd;
	Cz  = drag2chezy(cd);

	% at end-points
	if (nargin()>1 && iscentral)
		Q   = obj.Qc_;
		%Q   = mid(Q);
		w   = mid(w);
		% TODO hc from zc, to avoid recursive inner2outer
		h   = mid(h);
		Cz  = mid(Cz);
	else
		Q   = obj.Q_;
	end

	U   = Q./(h.*w);

	% sediment transport
	Qs0 = total_transport_engelund_hansen(Cz,d_mm,U(:,1),h,w);
	Qs  = river_tide_transport_scale(abs(U(:,2))./abs(U(:,1)),5).*Qs0;

	% buffer for values on boundary
	Qs = [0;Qs;0];

	% left boundary condition
	p       = obj.bc_Qs(1).p;
	val     = obj.bc_Qs(1).val;
	if (iscentral)
		Qs(1)   = p*val + (1-p)*(1.5*Qs(2)-0.5*Qs(3));
	else
		Qs(1)   = p*val + (1-p)*(2*Qs(2)-Qs(3));
	end

	% right boundary condition
	p       = obj.bc_Qs(2).p;
	val     = obj.bc_Qs(2).val;
	if (iscentral)
		Qs(end) = p*val + (1-p)*(1.5*Qs(end-1)-0.5*Qs(end-2));
	else
		Qs(end) = p*val + (1-p)*(2*Qs(end-1)-Qs(end-2));
	end

end % sediment_transport

