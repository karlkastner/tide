% Fri  7 Aug 19:07:20 +08 2020
% Karl Kastner, Berlin
function [out] = sediment_transport(obj,t,iscentral)
%	ni   = obj.ni;
%	Qs  = zeros(ni(end-1),1);
%	Qs0 = zeros(ni(end-1),1);
	out = struct();
	for idx=1:length(obj.rt)
		[out(idx).Qs,out(idx).Qs0] = obj.rt(idx).sediment_transport(t,iscentral);
	end


	% apply junction conditions
	for idx=1:length(obj.junction_Qs)
		[cid,eid] = feval(obj.junction_Qs{idx});

		% fetch values
		for cdx=1:length(cid)
			if (1 == eid(cdx))
				id = 1;
				Qs_j(cdx) = out(cid(cdx)).Qs(1);
			else
				ncx = length(obj.rt(cid(cdx)).x)-1;
				id  = ncx;
				Qs_j(cdx) = out(cid(cdx)).Qs(end);
			end

if (0)
			Q0_j = obj.rt(cid(cdx)).c.Q(id,1); 
			% TODO higher frequencies
			Qt_j = obj.rt(cid(cdx)).c.Q(id,2);
			w_j  = obj.rt(cid(cdx)).c.w(id);
			h_j  = obj.rt(cid(cdx)).c.h0(id);
			Cz_j = obj.rt(cid(cdx)).c.Cz(id);
else
			x_j       = obj.rt(cid(cdx)).Xi(eid(cdx));
			Q0_j(cdx) = obj.rt(cid(cdx)).Qc_(id,1); 
			% TODO higher frequencies
			Qt_j(cdx) = obj.rt(cid(cdx)).Qc_(id,2);
			w_j(cdx)  = obj.rt(cid(cdx)).width(x_j);
			h_j(cdx)  = obj.rt(cid(cdx)).h0(x_j);
			Cd_j(cdx) = obj.rt(cid(cdx)).cd(x_j);
			Cz_j(cdx) = drag2chezy(Cd_j(cdx));
end
			xmid  = mean(obj.rt(cid(cdx)).Xi);
			dir(cdx)  = sign(x_j-xmid);
		end % for cdx
		
		% flow direction with respect to bifurcation
		Q0_j = dir.*Q0_j;
		Qs_j = dir.*Qs_j;

		% compute division
		Qs_j = obj.sediment_division(Q0_j,Qt_j,Qs_j,w_j,h_j,Cz_j);

%dir
%Qs_j
		% flow direction with respect to channel
		Qs_j = dir.*Qs_j;
%Qs_j
%pause
		% write transport at junctions
		for cdx=1:length(cid)
			if (1 == eid(cdx))
				out(cdx).Qs(1) = Qs_j(cdx);
			else
				out(cdx).Qs(end) = Qs_j(cdx);
			end
		end % for cdx
	end % for idx
end % sediment_transport

