% Sat Jul 13 11:16:40 UTC 2013
% Karl Kästner, Berlin
%
%% export tide table to csv file
%
function export_csv(tidetable,csvname)

	t   = tidetable.time;
	t24 = tidetable.t24;
	drange  = tidetable.drange;
	D   = tidetable.range24;
	tl  = tidetable.tl;
	vl  = tidetable.vl;
	ldx = tidetable.ldx;
	hdx = tidetable.hdx;
	th  = tidetable.th;
	vh  = tidetable.vh;
	
	neap_dx   = tidetable.neap_dx;
	spring_dx = tidetable.spring_dx;

	TT = [t24 'R'*ones(size(D))  D;
              tl 'L'*ones(size(ldx)) vl;
	      th 'H'*ones(size(hdx)) vh;
	      floor(t(neap_dx)/86400)*86400 'N'*ones(size(neap_dx)) drange(neap_dx);
	      floor(t(spring_dx)/86400)*86400 'S'*ones(size(spring_dx)) drange(spring_dx)];
	[s idx] = sort(TT(:,1));
	TT = TT(idx,:);
	
	kdx = 1;
	ndx = 1;
	SS = zeros(size(TT,1),20);
	[y mo d h mi] = datevec(TT(1,1) / (86400));
	SS(1:1:7) = [y mo d h mi TT(1,2:3)];
	% year month day S/N H/L hour minute level H/L hour minute level
%	if (nargin() > 2 & ~isempty(pout))
	fid = fopen(csvname,'w');
%	else
%		fid = 1;
%	end

	for idx=2:size(TT,1)
		[y mo d h mi] = datevec(TT(idx,1) / (86400));
		TT(idx,1:7) = [y mo d h mi TT(idx,2:3)];
		if (      TT(idx,1) == SS(kdx,1) ...
			& TT(idx,2) == SS(kdx,2) ...
			& TT(idx,3) == SS(kdx,3) )
			% add another event to the current row
			if (TT(idx,6) ~= 'R')
				SS(kdx,6+4*ndx:6+4*ndx+3) = TT(idx,4:7); % hour minute type level
				ndx = ndx+1;
			else
				SS(kdx,5) = TT(idx,7); % tidal range
			end
		else
			% TODO print last row behind the loop
			s = sprintf('%4d %2d %2d %c %3d', SS(kdx, 1), SS(kdx, 2), SS(kdx, 3), SS(kdx,4),round(100*SS(kdx,5)));
			for jdx=1:ndx
				s = [s sprintf(' %2d:%02d %c %4d', ...
					SS(kdx,2+4*jdx), ...
					SS(kdx,2+4*jdx+1), ...
					SS(kdx,2+4*jdx+2), ...
					round(100*SS(kdx,2+4*jdx+3))) ];
			end
			fprintf(fid,'%s\n',s);	
			% proceed to the next row
			kdx = kdx+1;
			SS(kdx,1:3) = TT(idx,1:3); % y m d
			switch (TT(idx,6))
				case {'S', 'N'}
					SS(kdx,4  ) = TT(idx,6);   % S/N
					ndx = 0;
				case {'R'}
					SS(kdx,4  ) = 32;          % space (no spring or neap event that day)
					SS(kdx,5  ) = TT(idx,7);   % tidal range for that day
					nd = 0;
				otherwise % else
					SS(kdx,4  ) = 32; % space (no spring or neap event that day)
					SS(kdx,6:9) = TT(idx,4:7); % hour minute type level
					ndx = 1;
			end % switch
		end % if
	end % for idx
%	if (nargin() > 2 & ~isempty(pout))
	fclose(fid);
%	end
%	TT = SS(1:kdx,:);

end % function tidetable_export_csv

