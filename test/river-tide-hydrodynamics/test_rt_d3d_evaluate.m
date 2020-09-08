% Fri  4 Sep 20:00:10 +08 2020

function [nres, d3d] = test_river_tide_evaluate_d3d(rt,tid,pflag)
	if (nargin()<3)
		pflag = 0;
	end

	nres = NaN(1,6);
	d3d = Delft3D_Map();
	d3d.init(['mat/test-river-tide-',num2str(tid)]); 
	t  = d3d.time;
	n  = 1./(t(2)-t(1));
	dn = round(n/12);
	v=d3d.zs;
	u=d3d.u2;
	u=d3d.qu;
	
	if (pflag)
	% solution
	figure(200+tid);
	clf
	subplot(2,1,1)
	plot(squeeze(v(end-n:dn:end,2,:))');
	%colorbar
	
	% convergence in time
	nd=[];
	for k=1:size(v,1)/n-1;
		d=squeeze(v((k-1)*n+1:k*n,2,:)-v(k*n+1:(k+1)*n,2,:));
		nd(k)=max(abs(d(:)));
	end
	subplot(2,2,3)
	semilogy(nd,'.-')
	v_=v(k*n+1:(k+1)*n,2,:);
	%[max(max(abs(v_(:,1,1:100)))) min(max(abs(v_(:,1,1:100)))), max(abs(v_(:))),max(abs(v_(:,1,end)))]
	
	end

	k = size(v,1)/n-1;
	v_= v(k*n+1:(k+1)*n,2,:);
	u_= u(k*n+1:(k+1)*n,1,:);
	
	t = t(k*n+1:(k+1)*n);
	T = 1./[1,2,3,4];
	A = fourier_matrix_exp(T,t);
	v_ = squeeze(v_);
	u_ = squeeze(u_);
	f.z = (A \ v_).';
	f.u = (A \ u_).';
	
	z   = rt.out.z;
	res = z-f.z(1:size(z,1),1:size(z,2));
	nres(1:size(z,2)) = rms(res);

	% scale
	az10 = abs(f.z(1,2));
	h  = d3d.depth; 
	h0 = mean(h(k*n+1:(k+1)*n,2,1))
	e    = az10./h0;
	% 1 should actually scale as 1/e^2 !
	nres = nres./[e,az10,e.^2,e.^3,e.^4,NaN];
	
	if (pflag)
	namedfigure(300+tid,'water level');
	clf();
	namedfigure(400+tid,'discharge');
	clf();
	%nres = [];
	for idx=1:4 %size(z,2)
		figure(300+tid);
		subplot(2,2,idx);
		plot(real(f.z(:,idx)),'.-');
		hold on;
		if (1~=idx)
			plot(imag(f.z(:,idx)),'.-');
			plot(abs(f.z(:,idx)),'.-');
		else
			%plot(NaN(size(f.z,1),1));
			zb = d3d.zb;
			plot(squeeze(zb(1,2,:)),'.-','linewidth',1);
		end


		if (idx<=size(z,2)) %rt.out.z,2)
		zt=rt.z(idx-1);
		set(gca,'colororderindex',1);
		plot(real(zt),'--');
		if (1~=idx)
			plot(imag(zt),'--.');
			plot(abs(zt),'--.');
		else
			zb = rt.zb(1);
			plot(zb,'--.','linewidth',1);
			%plot(NaN(size(zt,1),1));
		end
		end
		xlim([0,rt.nx]);

		%if (idx == 1)
			%set(gca,'colororderindex',get(gca,'colororderindex')-1);
		%end

		figure(400+tid);
		subplot(2,2,idx);
		plot([real(f.u(:,idx)),imag(f.u(:,idx))]);
		hold on;
		plot(abs(f.u(:,idx)));

		if (idx<=size(z,2)) %rt.out.z,2)
		ut=rt.Q(idx-1);
		%ut=rt.u(idx-1);
		hold on;
		set(gca,'colororderindex',1);
		plot([real(ut),imag(ut)],'--');
		plot(abs(ut),'--');
		end
		xlim([0,rt.nx]);
	end % for idx
	end % if pflag
	%nres./abs(zt(1))

	% sediment transport
	v  = d3d.ssu2;
	qs = rt.sediment_transport(0,0);
	qs = qs.Qs(2:end-1);
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
		    xlim([0,rt.nx]);
		end % idx=1:4
	end % if ~isempty(v)

	
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

