% Sun 20 Feb 17:30:13 CET 2022
% apply boundary conditions
function y  = bcfun(rt,t,y);
	% for each channel
	for cdx = 1:rt.nc
		% for each end
		for edx=1:2
			switch (rt.bc(cdx,edx).var)
			case {'z1'}
				% transform
				% dQ_dx = w dz_dt = sum_j 1/(i o j w) dQ/dx
				cz0 = rt.bc(cdx,edx).rhs;
				if (isa(cz0,'function_handle'))
					cz0 = cz0(t);
				end
				ybc = 0;
				for idx=1:length(cz0)-1
				ybc = (  ybc ...
				       - 1i*idx*rt.omega*rt.w(cdx)*cz0(idx+1)*exp(+1i*idx*rt.omega*t) ...
				       + 1i*idx*rt.omega*rt.w(cdx)*conj(cz0(idx+1))*exp(-1i*idx*rt.omega*t) ...
		                      );
				% TODO p to Q
				% transform z into dQ_dx
				end
			case {'Q'}
				ybc = rt.bc(cdx,edx).rhs;
				if (isa(ybc,'function_handle'))
					ybc = ybc(t);
				end
			case {'junction'}
				% do not set values here, continue
				continue;
			otherwise
				error('here');
			end
			if (1 == edx)
				% first end
				at = rt.first(cdx);
				in = at+1;
			else
				% second end
				at = rt.last(cdx);
				in = at-1;
			end
			% distance towards junction point
			dx  = rt.x(at) - rt.x(in);

			% value at junction point
			y(at) = rt.bc(cdx,edx).q*ybc + (1-rt.bc(cdx,edx).q)*(y(in) + dx*ybc);
		end % for edx
	end % for cdx

	% conditions at junctions
	y = rt.apply_junction_condition(t,y);
end % bcfun

