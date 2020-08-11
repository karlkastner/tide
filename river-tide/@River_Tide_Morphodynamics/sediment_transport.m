% Fri  7 Aug 19:07:20 +08 2020
% Karl Kastner, Berlin
function [Qs,Qs0] = sediment_transport(obj,flag)
	d_mm = obj.sediment.d_mm;

	% discharge and channel properties at segment centres
	Q   = obj.Q_;
	w   = obj.width;
	h   = obj.h0;	
	cd  = obj.cd;
	Cz  = cd2chezy(cd);

	% at end-points
	if (nargin()>1 && flag)
		Q   = mid(Q);
		w   = mid(w);
		h   = mid(h);
		Cz  = mid(Cz);
	end

	U   = Q./(h.*w);

	% sediment transport
	Qs0 = total_transport_engelund_hansen(Cz,d_mm,U(:,1),h,w);
	Qs  = river_tide_transport_scale(abs(U(:,2))./abs(U(:,1)),5).*Qs0;
end

