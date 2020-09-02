% Thu 20 Aug 18:37:27 +08 2020
% TODO make transport scale a function of River_Tide
function Qs = sediment_division_geometric(obj,Q0,Qt,zt,Qs,w,h,Cz,d_mm)
%	d_mm = obj.rt(1).sediment.d_mm;

	% this has to be equal in all channels
	zt  = mean(zt);
	Czc = 0.5*(Cz(1) + Cz(2:3));
	hc  = 0.5*(h(1)  + h(2:3));
	wc  = 0.5*(abs(Q0(2:3)./Q0(1)).*w(1) + w(2:3));
	Ac  = hc.*wc;
	U0c = Q0(2:3)./Ac;

	Q0 = abs(Q0);

	U0     =  Q0./(w.*h);
	% TODO should be same as Qs(1), but currently 0 (due to bc)
	Qs0_in =  total_transport_engelund_hansen(Cz(1),d_mm,U0(1),h(1),w(1));
	Qs_in  =  river_tide_transport_scale(abs(Qt(1))./abs(Q0(1)),zt./hc,5,obj.opt.stokes_order).*Qs0_in;

	% mean transport
	Qs0 = total_transport_engelund_hansen(Czc,d_mm,U0c,hc,wc);
	if (~obj.opt.ignore_rt)
		Qs  = river_tide_transport_scale(abs(Qt(2:3))./abs(Q0(2:3)),zt./hc,5,obj.opt.stokes_order).*Qs0;
		Qs  = Qs./sum(Qs)*(-Qs_in);
	else	
		% without influence of river tide on the division
		Qs = Qs0./sum(Qs0).*(-Qs_in);
	end
	Qs = [-sum(Qs),Qs];
end


