% Mon  5 Oct 13:38:55 +08 2020
% intepolate bed level from a coarse mesh to a fine mesh with half the segment length
% for multigrid
function z = mg_interpolate(obj,z)
	% for each channel
	for cdx=1:obj.nc 
		if (obj.nxc(cdx) > 3)
			z(obj.nci(cdx):obj.nci(cdx+1)-1) = mg_interpolate(zp(obj.lower.nci(cdx):obj.lower.nci(cdx+1)-1));
		else
			% end of recursion, no projection
			z(obj.nci(cdx):obj.nci(cdx+1)-1) = zp(obj.lower.nci(cdx):obj.lower.nci(cdx+1)-1);
		end % if nxc > 3
	end % for cdx
end % function mg_interpolate

