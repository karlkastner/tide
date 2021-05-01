% Thu 29 Oct 16:04:02 +08 2020
function f  = odefunk_3(obj, Q, Q2, Qhr, zs, zb, h0, dh0_dx, dz0_dx, w, w_x, Cd, ...
							      cf, D1_dx, D2_dx)
	
	cf = cf/pi;

	g = obj.g;
	fmic = @fourier_multiplicative_interaction_coefficients;
%	nf   = obj.nf;
	nf = size(Q,2)-1;
	nx = size(Q,1);

	% allocate memory
	f     = zeros(nx,4,nf+1);

	A     = w.*([zs(:,1)-zb,zs(:,2:end)]);
	A2    = fmic(A,A,nf+1);
	A3    = fmic(A,A2,nf+1);
	A4    = fmic(A,A3,nf+1);

	Qt    = obj.fourier_derivative(Q);
	%Q_tt  = obj.fourier_derivative(Q_t);
	A2Qt  = fmic(A2,Qt,nf+1);

	% powers of Q, stored in third dimension
	Qp        = Q;
	Qp(:,:,2) = Q2;
	for idx=3:size(cf,2)-1
		%Q3    = fmic(Q,Q2,nf+1);
		Qp(:,:,idx) = fmic(Q,Qp(:,:,idx-1),nf+1);
	end
	
	% derivatives of Q
	% TODO directly from bvp2m
	Qx  = D1_dx*Q;
	Qxx = D2_dx*Q;

	% signed square of friction term : F*Q = Q|Q|, F_t*Q = Q|Q|_t
	Qhr_ = Qhr.^2;
	F = cf(:,2).*Qhr_;

	for idx=2:size(cf,2)
		Qhr_ = Qhr_./Qhr;
		Qhr_(Qhr == 0) = 0;
		F = F + cf(:,idx).*Qp(:,:,idx-1).*Qhr_;
		if (2 == idx)
			F_t =       (idx-1)*cf(:,idx).*Qhr_;
		else
			F_t = F_t + (idx-1)*cf(:,idx).*Qp(:,:,idx-2).*Qhr_;
		end
	end
                  
	% coefficients of wave equation, not separated into frequency components

%	c(:,:,1) = -g*A4./w;
	c = -g*A4./w;
	f = obj.coefficient_frequency_components(f,c,Qxx,1,0);
	% c2*Q_x = (- A^2 Q_t - Cd*w*3*|Q|Q + g A^4 1/w^2 w_x) * Q_x 

	% TODO allof for variable length cf
	c = (         +A2Qt ...
                     + g*(w_x./w.^2).*A4 ...
		     + 3*Cd.*w.*F ...
		   );
	f = obj.coefficient_frequency_components(f,c,Qx,2,0);

  	% c3*Q = A^3*Q_tt + Cd w A |Q|Q_t 
	c = A3;
	f = obj.coefficient_frequency_components(f,c,Q,3,2);

	c = Cd.*w.*A.*F_t;
	f = obj.coefficient_frequency_components(f,c,Q,3,1);

	% advective-acceleration-term
	if (obj.opt.ode.advective_acceleration)
		Ax    = D1_dx*A;
		AQ    = fmic(A,Q,nf+1);
		AQ2   = fmic(A,Q2,nf+1);
		AxQ2  = fmic(Ax,Q2,nf+1);
		AAx   = fmic(A,Ax,nf+1);
		AQQx  = fmic(AQ,Qx,nf+1);
		AAxQt = fmic(AAx,Qt,nf+1);
		Qxt   = obj.fourier_derivative(Qx);
		A2Qxt = fmic(A,Qxt,nf+1);
		% TODO
		%AAxQt = dt*AAxQt;
		% c Qxx
	 	c = AQ2;
		f = obj.coefficient_frequency_components(f,c,Qxx,1,0);
		% c Qx
                c = 2*A2Qt + 4*AQQx - 3*AxQ2;
		f = obj.coefficient_frequency_components(f,c,Qx,1,0);
		c = -2*AAxQt + 2*A2Qxt; 
		f = obj.coefficient_frequency_components(f,c,Q,1,0);
	end

end % odefunk_3

