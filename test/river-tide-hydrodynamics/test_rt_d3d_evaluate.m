% Fri  4 Sep 20:00:10 +08 2020

function [nres, d3d] = test_rt_d3d_evaluate(rt,tid,pflag)
	if (nargin()<3)
		pflag = 0;
	end

	nres = NaN(1,6);
	d3d = Delft3D_Map();
	d3d.init(['mat/test-river-tide-',num2str(tid)]); 
	if (isempty(d3d.map))
		nres = NaN(1,4);
		return;
	end

	t  = d3d.time;
	nt = length(t);
	% number of samples per day
	dt=(t(2)-t(1));
	n  = 1./dt;
	% number of days - 1
	k = size(nt,1)/n-1;
	dn = round(n/12);

	d3d_x = inner2outer(head(d3d.Y,1));
	zs    = d3d.zs;
	h     = d3d.depth; 
	q     = d3d.qu;
	
	if (pflag)

	% solution
	figure(200+tid);
	clf();
	subplot(2,1,1)
	plot(squeeze(zs(end-n:dn:end,2,:))');
	ylabel('z_s');
	%colorbar
	
	% convergence in time
	nd=[];
	for k=1:size(zs,1)/n-1;
		d=squeeze(zs((k-1)*n+1:k*n,2,:)-zs(k*n+1:(k+1)*n,2,:));
		nd(k)=max(abs(d(:)));
	end

	subplot(2,2,3);
	semilogy(nd,'.-');
	ylabel('max|zs(t)-zs(T)|');
	zs_=zs(k*n+1:(k+1)*n,2,:);
	
	end

	t  = t(k*n+1:(k+1)*n);
	zs_= zs(k*n+1:(k+1)*n,2,:);
	q_ = q(k*n+1:(k+1)*n,1,:);

	T = 1./[1,2,3,4];
	A   = fourier_matrix_exp(T,t);
	zs_ = squeeze(zs_);
	q_  = squeeze(q_);
	%f.z = (A \ zs_).';
	%f.q = (A \ q_).';
	f.z = dt*(A' * zs_).';
	f.q = dt*(A' * q_).';

	% unpack coefficients
	f.z = f.z(:,[1,2:2:end]);
	f.q = f.q(:,[1,2:2:end]);
	%f.z = [f.z(:,1),f.z(:,2:2:end)];
	%f.q = [f.q(:,1),f.q(:,2:2:end)];
	
	z   = rt.channel(1).z;
	res = z-f.z(1:size(z,1),1:size(z,2));
	nres(1:size(z,2)) = rms(res);

	% scale
	az10 = abs(f.z(1,2));
	h0 = mean(h(k*n+1:(k+1)*n,2,1))
	e    = az10./h0;
	% 1 should actually scale as 1/e^2 !
	nres = nres./[e,az10,e.^2,e.^3,e.^4,NaN];
	
	if (pflag)
	namedfigure(300+tid,'water level, RT vs d3d');
	clf();
	namedfigure(400+tid,'Discharge, RT vs d3d');
	clf();
	%nres = [];
	for idx=1:4
		figure(300+tid);
		subplot(2,2,idx);
		plot(d3d_x,real(f.z(:,idx)),'--');
		title(['k = ',num2str(idx-1)])
		hold on;
		if (1~=idx)
			plot(d3d_x,imag(f.z(:,idx)),'--');
			plot(d3d_x,abs(f.z(:,idx)),'--');
		else
			zb = d3d.zb;
			plot(d3d_x,squeeze(zb(1,2,:)),'--','linewidth',1);
		end


		if (idx<=size(z,2)) %rt.channel(1).z,2)
		zt=rt.channel(1).waterlevel(idx-1);
		set(gca,'colororderindex',1);
		plot(rt.channel(1).x,real(zt),'-');
		if (1~=idx)
			plot(rt.channel(1).x,imag(zt),'-');
			plot(rt.channel(1).x,abs(zt),'-');
		else
			zb = rt.channel(1).zb;
			plot(rt.channel(1).x,zb,'-','linewidth',1);
			%plot(NaN(size(zt,1),1));
		end
		end
		%xlim([0,rt.nx]);
		xlim([0,max(rt.channel(1).x)]);
		if (1==idx)
			legend('d3d','d3d','rt','rt')
		else
			legend('d3d','d3d','d3d','rt','rt','rt')
		end
		%if (idx == 1)
			%set(gca,'colororderindex',get(gca,'colororderindex')-1);
		%end

		figure(400+tid);
		subplot(2,2,idx);
		if (1==idx)
			plot(mid(d3d_x(2:end-1)),f.q(:,idx),'--');
			hold on;
		else
		plot(mid(d3d_x(2:end-1)),[real(f.q(:,idx)),imag(f.q(:,idx))],'--');
		hold on;
		plot(mid(d3d_x(2:end-1)),abs(f.q(:,idx)),'--');
		end

		if (idx<=size(z,2)) %rt.channel(1).z,2)
			ut=rt.channel(1).discharge(idx-1);
			%ut=rt.u(idx-1);
			set(gca,'colororderindex',1);
			if (1==idx)
				plot(rt.channel(1).x,ut,'-');
			else
			plot(rt.channel(1).x,[real(ut),imag(ut)],'-');
			plot(rt.channel(1).x,abs(ut),'-');
		end
		end
		if (1==idx)
			legend('d3d','d3d','rt','rt')
		else
			legend('d3d','d3d','d3d','rt','rt','rt')
		end
		%xlim([0,rt.nx]);
		xlim([0,max(rt.channel(1).x)]);
	end % for idx
	end % if pflag
	%nres./abs(zt(1))

	% sediment transport
	v  = d3d.ssu([]);
	try 
		qs = rt.channel(1).sediment_transport(0,0);
		qs = qs.Qs(2:end-1);
	catch	e
		e
		Qs = NaN;
	end
	if (isempty(v))
		try
		close(500+tid);
		catch e
		end
	else
		namedfigure(500+tid,'transport');
		clf();
		v_=squeeze(v(k*n+1:(k+1)*n,2,:,1));
		f.qs = (A \ v_)';
		for idx=1:4
			subplot(2,2,idx);
			plot([real(f.qs(:,idx)),imag(f.qs(:,idx))]);
			hold on;
			plot(abs(f.qs(:,idx)));
			if (1==idx)
				hold on;
				set(gca,'colororderindex',1);
				plot([real(qs),imag(qs)],'--');
				plot(abs(qs),'--');
				res = qs - f.qs(1:length(qs),idx);
				nres(6) = rms(res)./abs(f.qs(1));
			end % if 1 == idx
		    xlim([0,rt.channel(1).nx]);
		end % idx=1:4
	end % if ~isempty(v)

	% erify D3D

	nf = 3;
	x=d3d.Y; x=x(1,2:end-1);
	nt = n;
	nx=size(q,1);
	t=d3d.time;
	t=t(end-nt+1:end);
	g = rt.rt.g;
	Cd = d3d.Chezy_u;
	zb = d3d.zb;
	h  = d3d.depth;
	zb = zb(:,:,2:end-1);
	zb = squeeze(repmat(zb(:,2,:),nt,1,1));
	h = h(:,:,2:end-1);
	Q = squeeze(q(end-nt+1:end,1,:));
%	zs = squeeze(zs(end-nt+1:end,2,:));
	h  = squeeze(h(end-nt+1:end,2,:));
	Cd = squeeze(Cd(end-nt+1:end,2,2:end-1));
%	zs = flat(zs);
	h = mid(h,2);
	zb = mid(zb,2);
	Cd = mid(Cd,2);
	Q = flat(Q');
	h = flat(h');
	zb = flat(zb');
	Cd = flat(Cd');
	% TODO, exact
	w = ones(size(zb));
	A = h.*w;
	[resc,resm,rescf,resmf] = residual_swe(t,x,A,Q,g,w,zb,Cd,nf);
	if (pflag)
	figure(1000+tid);
	clf;
	for idx=1:size(resmf,1);
		subplot(2,4,idx);
		plot(abs(resmf(idx,:))');
		xlim([0.5,rt(1).hydrosolver.nx+0.5]);
	end
	end
	
	if (0)
	 z1 = zt(22); z = z1*exp(2i*pi*t); figure(3);clf; plot([real(z),v_(:,22)])
	end

	if (0)
		% check initial condition
		zs = d3d.zs; figure(10); clf; zs=squeeze(zs(:,2,:)); plot((zs(1:2*24:end,:)-zs(1,:))')
	end
%catch e
%	e
%end	
end

