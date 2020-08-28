% Fri 15 Dec 16:58:52 CET 2017
% Karl Kastner, Berlin
%
%% decompose the tide into a right and left travelling wave,
%% i.e. into incoming and reflected wave
%% TODO subtract forcing term
%
% function [Q1lr, z1lr, obj] = decompose(obj)
%
function decompose(x,w0,z1,Q0,Qt)

	dx    = diff(x);

	% element edges to centres
	xc = mid(x); %0.5*(x(1:end-1)+x(2:end));
	wc = mid(w0);
	if (0)
		Q1c = 0.5*(Q1(1:end-1)+Q1(2:end));
		c  = obj.odefun(cdx,xc,Q1c);
	else
		switch (obj.opt.hmode)
		case {'matrix'}
			if (obj.opt.dischargeisvariable)
				y = [obj.z(cdx,0); Q0(1); Q1];
				k=2;
			else
				y = [obj.z(cdx,0); Q1];
				k = 2;
			end
		otherwise
			y = obj.Q(cdx,1);
			k = 1;
		end
		%if (obj.opt.oflag(2))
		%	y = [y;obj.Q_];
		%end

		c = obj.odefun(cdx,x,y);
		c_ = 0.5*(c(1:end-1,:,k)+c(2:end,:,k));	
	end
	r  = roots2(c_);
	
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
end % decompose

