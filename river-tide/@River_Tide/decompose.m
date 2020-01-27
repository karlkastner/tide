% Fri 15 Dec 16:58:52 CET 2017
%% decompose the tide into a right and left travelling wave,
%% i.e. into incoming and reflected wave
%% TODO subtract forcing term
function [Q1lr, z1lr, obj] = decompose(obj)
	Q1    = obj.Q(1);
	x     = obj.x;
	dx    = diff(x);
	omega = obj.omega;

	% element edges to centres
	xc = 0.5*(x(1:end-1)+x(2:end));
	wc = obj.fun.width(xc);
	if (0)
		Q1c = 0.5*(Q1(1:end-1)+Q1(2:end));
		c  = obj.odefun(xc,Q1c);
	else
		switch (obj.opt.hmode)
		case {'matrix'}
			y = [obj.z0; Q1];
		otherwise
			y = Q1;
		end
		if (obj.opt.o2)
			y = [y;obj.Q(2)];
		end

		c = obj.odefun(x,y);
		c = 0.5*(c(1:end-1,:)+c(2:end,:));	
	end
	r  = roots2(c);
	
	% match values at segment end points (continuity)
	% constant coefficients in each segment
	% Q(-1/2 dx) = Qm exp(-1/2 rm dx) + Qp exp(-1/2 rp dx)
	% Q(+1/2 dx) = Qm exp(+1/2 rm dx) + Qp exp(+1/2 rp dx)
	A   = [diag(sparse(exp(-0.5*dx.*r(:,1)))), diag(sparse(exp(-0.5*dx.*(r(:,2)))));
	       diag(sparse(exp(+0.5*dx.*r(:,1)))), diag(sparse(exp(+0.5*dx.*(r(:,2)))))];
	% subtract constant part
	Q1_   = [Q1(1:end-1)-c(1:end-1,4); Q1(2:end)-c(2:end,4)];
	Q1lrc = A \ Q1_;
	Q1lrc = reshape(Q1lrc,[],2);

	% z = r/(i omega w) Q, as coefficients are constant along segments
	z1lrc = [r(:,1)./(1i*omega*wc).*Q1lrc(:,1), r(:,2)./(1i*omega*wc).*Q1lrc(:,2)];
	
	% segment centres to end points
	n = length(x);

	for idx=1:2
		A = sparse(n,n-1,n);
		A(sub2ind(size(A),1:n-1,1:n-1)) = exp(-0.5*dx.*r(:,idx));
		A(end,end)                      = exp(+0.5*dx(end).*r(end,idx));
		Q1lr(:,idx) = A*Q1lrc(:,idx);
		z1lr(:,idx) = A*z1lrc(:,idx);
	end % for idx

%	qlr = bsxfun(@times,Qlr,1./wc);
%	zlr = obj.q1_to_z1(xc,qlr,obj.omega);
end

