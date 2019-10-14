% Sun  8 Oct 13:08:39 CEST 2017
% transform arbitrary to cs-integrated discharge boundary condition
function obj = bc_transformation(obj)
	n_z0 = 0;
	n_Q0 = 0;			

	bc = obj.bc;

	% for left and right end
	for id=1:2
		% mean component
		switch (bc(id,1).var)
		case {'z','z0'}
			obj.z0_downstream = bc(id,1).rhs;
			n_z0 = n_z0+1;
		case {'Q','Q0'}
			obj.Q0_ = bc(id,1).rhs;
			n_Q0 = n_Q0+1;
		otherwise
			error('bcfun')
		end
		if (bc(id,1).p ~= 1)
			error('only dirichlet condition supported for mwl so far');
		end 

		% for each tidal frequency component
		for jd=2:2 %size(bc,2)
			switch (bc(id,2).var)
			case {'z'}
				% dQ/dx = -1i*o*z
				w0 = obj.fun.width(obj.Xi(id));
				omega_j = (jd-1)*obj.omega;
				bc(id,jd).rhs = -1i*omega_j*w0*bc(id,jd).rhs;
				if ( bc(id,jd).p(2) ~= 0)
					error('bc of type dz/dx not yet implemented');
				end
				p = [0, bc(id,jd).p(1), 0];
				bc(id,jd).p   = p;
			case {'Q'}
				% nothing to do, native format
				bc(id,jd).rhs = bc(id,jd).rhs;
			case {'q'}
				% Q = w*q
				bc(id,jd).rhs = bc(id,jd).rhs*W0;
			case {'u'}
				error('u condition not yet implemented');
				% TODO Q = (z0-zb)*w0*u
			otherwise
				error('bc');
			end % switch
		end % for jd, frequency components
	end % for id, left and right end of doamin
	if (~isreal(obj.z0_downstream) || ~isreal(obj.Q0_))
		error('mean water level and mean discharge must be real');
	end
	if (n_z0 ~= 1 && n_Q0 ~= 1)
		error('Mean water level and mean discharge must be specified on opposit ends of the domain');
	end

	obj.bc = bc;
end % bc_transform

