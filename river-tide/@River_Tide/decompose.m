% Fri 15 Dec 16:58:52 CET 2017
% Karl Kastner, Berlin
%
%% decompose the tide into a right and left travelling wave,
%% i.e. into incoming and reflected wave
%% TODO subtract forcing term
%
% function [Q1lr, z1lr, obj] = decompose(obj)
%
function [Q1lr,z1lr] = decompose(obj,x,Q,zs,zb,w0,Cd)
	% TODO, quick fix
	%cdx = 1;

	omega = obj.omega;
	dx    = diff(x);

	% element edges to centres
	% TODO, w and x should be passed
	xc = mid(x);
	wc = mid(w0);
	
%	% TODO quick fix, pass Q2 and Q3 as well
%	Qt = Q1;
%	for idx=2:length(obj.opt.oflag)
%	if (obj.opt.oflag(idx))
%		Qt = [Qt;zeros(size(Q1,1),1)];
%	end
%	end
%
%	if (obj.opt.dischargeisvariable)
%		y = [z0; Q0(1); Qt];
%		k=2;
%	else
%		y = [z0; Qt];
%		k = 2;
%	end
	k = 2;

	% c  = obj.odefun(cdx,x,y);
%        f = obj.odefun(x, [Q0, Qt], zs, zb, w0, Cd, dw_dx, D1_dx, D2_dx);
	c = obj.odefun(x, Q, zs, zb, w0, Cd); %, dw_dx, D1_dx, D2_dx)

	c_ = 0.5*(c(1:end-1,:,k)+c(2:end,:,k));	
	r  = roots2(c_);
	
	% match values at segment end points (continuity)
	% constant coefficients in each segment
	% Q(-1/2 dx) = Qm exp(-1/2 rm dx) + Qp exp(-1/2 rp dx)
	% Q(+1/2 dx) = Qm exp(+1/2 rm dx) + Qp exp(+1/2 rp dx)
	A   = [diag(sparse(exp(-0.5*dx.*r(:,1)))), diag(sparse(exp(-0.5*dx.*(r(:,2)))));
	       diag(sparse(exp(+0.5*dx.*r(:,1)))), diag(sparse(exp(+0.5*dx.*(r(:,2)))))];

	% subtract constant part
	Q1_   = [Q(1:end-1,2)-c(1:end-1,4); Q(2:end,2)-c(2:end,4)];
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

end % decompose

