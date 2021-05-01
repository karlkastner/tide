% Mon  5 Oct 13:17:34 +08 2020
% restrict the bed level to a mesh with half the spatial resolution
% this does not care about varying segment length, as this is negligible,
% when the mesh segment length varies smoothly along the channel
function zp = mg_restrict(obj,z)
	% for each channel
	for cdx=1:obj.nc 
		if (obj.nxc(cdx) > 3)
			zp(obj.lower.nci(cdx):obj.lower.nci(cdx+1)-1)) = mg_restrict(z(obj.nci(cdx):obj.nci(cdx+1)-1)); 
		else
			% end of recursion, no restriction
			zp(obj.lower.nci(cdx):obj.lower.nci(cdx+1)-1)) = z(obj.nci(cdx):obj.nci(cdx+1)-1); 
		end % if nxc > 3
	end % for cdx
end % function mg_restrict

