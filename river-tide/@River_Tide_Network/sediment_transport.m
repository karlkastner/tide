% Fri  7 Aug 19:07:20 +08 2020
% Karl Kastner, Berlin
%
%% compute sediment transport for the channel network, including routing at
%% junctions
function sediment_transport(obj,t,ddir)
	%out = struct();
	for cdx=1:obj.nc
		obj.channel(idx).sediment_transport(t,ddir);
		%[out(cdx).Qs,out(cdx).Qs0] = obj.sediment_transport_(cdx,t,ddir);
	end

	% apply junction conditions
	for idx=1:length(obj.junction_Qs)
		[cid,eid] = feval(obj.junction_Qs{idx});

		% fetch values
		for cdx=1:length(cid)
			if (1 == eid(cdx))
				id = 1;
				Qs_j(cdx) = obj.channel(cid(cdx)).Qs(1);
				%Qs_j(cdx) = out(cid(cdx)).Qs(1);
			else
				%nxc = length(obj.rt(cid(cdx)).x)-1;
				nx = obj.hydrosolver.nx(cid(cdx));
				%length(obj.rt(cid(cdx)).x)-1;
				idc  = nx-1;
				id   = nx;
				%Qs_j(cdx) = out(cid(cdx)).Qs(end);
				Qs_j(cdx) = obj.channel(cid(cdx)).Qs(end);
			end

if (0)
			Q0_j = obj.rt(cid(cdx)).c.Q(id,1); 
			% TODO higher frequencies
			Qt_j = obj.rt(cid(cdx)).c.Q(id,2);
			w_j  = obj.rt(cid(cdx)).c.w(id);
			h_j  = obj.rt(cid(cdx)).c.h0(id);
			Cz_j = obj.rt(cid(cdx)).c.Cz(id);
else
			x_j       = obj.hydrosolver.xi(cid(cdx),eid(cdx));
			% TODO, why c, should be at end points
			Q0_j(cdx) = obj.channel(cid(cdx)).Q(id,1); 
			%Q0_j(cdx) = obj.rt(cid(cdx)).Qc(id,1); 
			% TODO higher frequencies
			Qt_j(cdx) = obj.channel(cid(cdx)).Q(id,2);
			zt_j(cdx) = obj.channel(cid(cdx)).z(id,2);
			w_j(cdx)  = obj.width(cid(cdx),x_j);
			h_j(cdx)  = obj.h0(cid(cdx),x_j);
			Cd_j(cdx) = obj.cd(cid(cdx),x_j,h_j(cdx));
			Cz_j(cdx) = drag2chezy(Cd_j(cdx));
end
			%xmid  = mean(obj.rt(cid(cdx)).Xi);
			xmid      = mean(obj.hydrosolver.xi(cid(cdx),:));
			dir(cdx)  = sign(x_j-xmid);
		end % for cdx
		
		% flow direction with respect to bifurcation
		Q0_j = dir.*Q0_j;
		Qs_j = dir.*Qs_j;

		% compute division
		Qs_j = obj.bifurcation.sediment_division(Q0_j,Qt_j,zt_j,Qs_j,w_j,h_j,Cz_j,obj.sediment.d_mm);

		% flow direction with respect to channel
		Qs_j = dir.*Qs_j;

		% write transport at junctions
		for cdx=1:length(cid)
			if (1 == eid(cdx))
				obj.channel(cdx).Qs(1) = Qs_j(cdx);
				%out(cdx).Qs(1) = Qs_j(cdx);
			else
				obj.channel(cdx).Qs(end) = Qs_j(cdx);
				%out(cdx).Qs(end) = Qs_j(cdx);
			end
		end % for cdx
	end % for idx
end % sediment_transport

