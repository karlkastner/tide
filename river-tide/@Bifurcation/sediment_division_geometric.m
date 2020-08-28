% Thu 20 Aug 18:37:27 +08 2020
function Qs = sediment_division_geometric(obj,Q0,Qt,Qs,w,h,Cz)
	d_mm = obj.rt(1).sediment.d_mm;

	Czc = 0.5*(Cz(1) + Cz(2:3));
	hc = 0.5*(h(1)  + h(2:3));
	wc = 0.5*(abs(Q0(2:3)./Q0(1)).*w(1) + w(2:3));
	Ac = hc.*wc;
	U0c = Q0(2:3)./Ac;

	Q0 = abs(Q0);

	U0     =  Q0./(w.*h);
	Qs0_in =  total_transport_engelund_hansen(Cz(1),d_mm,U0(1),h(1),w(1));
	Qs_in  =  river_tide_transport_scale(abs(Qt(1))./abs(Q0(1)),5).*Qs0_in;

	% mean transport
	Qs0 = total_transport_engelund_hansen(Czc,d_mm,U0c,hc,wc);
	if (~obj.opt.ignorertfordivision)
		Qs  = river_tide_transport_scale(abs(Qt(2:3))./abs(Q0(2:3)),5).*Qs0;
		Qs  = Qs./sum(Qs)*(-Qs_in);
	else	
		% without influence of river tide on the division
		%$Qs = Qs0;
		Qs = Qs0./sum(Qs0).*(-Qs_in);
	end
% without rt-interaction, split 
%Qs_in
	Qs = [-sum(Qs),Qs];
end


