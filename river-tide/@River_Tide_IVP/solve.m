% Sun  8 Oct 13:08:39 CEST 2017

function z0 = solve_backwater(obj,x,Q0,Qt)
		% recompute the backwater curve
		czfun = @(x,h) drag2chezy(obj.cd(x,h));
		% Qt class member becomes only available after the solver has terminated
		Qtfun = @(xi) interp1(x,Qt,xi,'linear');

		obj.backwater.sopt.InitialStep = 1e-4*diff(obj.Xi);

		if (obj.bc(1,1).var == 'z')
			% left to right
			z0_downstream = obj.bc(1,1).rhs;
			[x_, h0_, z0_] = obj.backwater.solve( ...
					Q0(1), ...
					Qtfun, ...
					czfun,...
					@obj.width, ...
					@obj.zb, ...
					obj.z0_downstream, ...
					obj.Xi);
		else
			% right to left
			z0_downstream = obj.bc(2,1).rhs;
			[x_, h0_, z0_] = obj.backwater.solve( ...
					Q0(1), ...
					Qtfun, ...
					czfun,...
					@obj.width, ...
					@obj.zb, ...
					obj.z0_downstream, ...
					[obj.Xi(2),obj.Xi(1)]);
			% correct direction
			x_  = flipud(x_);
			% h0_ = flipud(h0_);
			z0_ = flipud(z0_);
		end

		% interpolate here to x
		z0  = @(x) interp1(x_,z0_,obj.x,obj.opt.imethod);

end

