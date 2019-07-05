% Sat Jul 13 11:16:40 UTC 2013
%
%% plot average neap and spring tide
%
function tidetable = plot_neap_spring(tidetable,pflag,prefix)
	if (nargin() < 2 || isempty(pflag))
		pflag = false;
	end
	if (nargin() < 3)
		prefix = [];
	end

	% number of groups to split the neap and tide in (for testing)
	ks = 1;
	kn = 1;


	t   = tidetable.time;
	dt  = tidetable.dt;
	level = tidetable.level;
	u   = tidetable.umag;
	dir = tidetable.udir;
	ux  = tidetable.ux;
	uy  = tidetable.uy;
	spring_dx = tidetable.spring_dx;
	neap_dx   = tidetable.neap_dx;

	% plot the water level and the envelope curves
	namedfigure(1,'Neap spring cycles over one year');
	clf();
	plot(tidetable.time,[tidetable.level tidetable.dmin tidetable.dmax tidetable.drange]);
	hold on;
	plot(tidetable.tl,tidetable.vl,'k.');
	plot(tidetable.th,tidetable.vh,'k.');
	plot(tidetable.time(tidetable.neap_dx),tidetable.level(tidetable.neap_dx),'r.');
	plot(tidetable.time(tidetable.spring_dx),tidetable.level(tidetable.spring_dx),'r.');
	q = quantile(tidetable.drange,[1/7 6/7]);
	fdx = find(tidetable.drange.*(tidetable.drange > q(1) & tidetable.drange < q(2)));
	d = tidetable.drange;
	d(fdx) = NaN;
	plot(tidetable.time,d,'r'); hold on;
	lt=[]; ls = [];
	for y=2013:2013+19;
		for m=1:12;
			%ls{end+1} = sprintf('01-%02d-%04d',m,y);
			ls        = [ls; sprintf('01-%02d',m)];
			lt(end+1) = datenum(sprintf('01-%02d-%04d',m,y),'dd-mm-yyyy');
			%lt(end+1) = datenum(ls(end,:),'dd-mm-yyyy');
		end
	end

	%lt = datenum(ls);
	set(gca,'xtick',lt);
	set(gca,'xticklabel',ls);
	n = 365.25/dt;
	n = round(n);
	xlim([t(1) t(n)] + 8*30);
	if (pflag)
		pdfprint(1,['img/',prefix,'tide-year-2013']);
	end
%	figure(5)
	plot(tidetable.t24,tidetable.range24,'.g');

%	preparePrint();
%	print('-depsc','tide-year-2013.eps');
%	system('epstopdf tide-year-2013.eps');


	fdx = tidetable.spring_dx;
	ns  = round(13/24 * 1/dt);
	n0  = round(13/24 * 1/dt);
	n   = n0 + ns;
	idx = ones(size(fdx(2:end-1)))*(-n:n) + fdx(2:end-1)*ones(1,2*n+1);
	spring_l   = level(idx');
	spring_l_  = spring_l; % unshifted copy
	spring_T   = t(idx');
	spring_ux  = ux(idx');
	spring_uy  = uy(idx');
	spring_u   = u(idx');
	spring_dir = dir(idx');
	% shift the maximum spring amplitude to the centre of the time axis

	% make the phase shift of all the spring tides equal
	t_ = dt*(1:2*n+1)';
	T  = 1; % diurnal
	omega = 2*pi*t_/T;
	A = [sin(omega) cos(omega)];
	c = round(0.5*size(spring_l,1));
	gdx = [];
	for idx=1:size(spring_l,2)
		Y = A \ spring_l(:,idx);
		phi = atan2(Y(1),Y(2));
		phi = max(min(phi,phi+2*pi),phi-2*pi);
		delta = -phi/(2*pi)*T/dt;
		% monitor the maximum shift and exclude those which are shifted more than the plot range
		if (abs(delta) < ns)
			% shift the maximum / minimum to the centre
			spring_l(:,idx)   = circshift(spring_l(:,idx), round(delta));
%			spring_T(:,idx)   = circshift(spring_T(:,idx), round(delta));
			spring_ux(:,idx)  = circshift(spring_ux(:,idx), round(delta));
			spring_uy(:,idx)  = circshift(spring_uy(:,idx), round(delta));
			spring_u(:,idx)   = circshift(spring_u(:,idx), round(delta));
			spring_dir(:,idx) = circshift(spring_dir(:,idx), round(delta));
			gdx(end+1) = idx;
		end
	end % for idx
	length(spring_l)/length(gdx)
	% select the good ones
	l0 = size(spring_l,2)
	spring_l = spring_l(ns+1:end-ns,gdx);
%	spring_T = spring_T(ns+1:end-ns,gdx);
	spring_u = spring_u(ns+1:end-ns,gdx);
	spring_ux = spring_ux(ns+1:end-ns,gdx);
	spring_uy = spring_uy(ns+1:end-ns,gdx);
	spring_dir = spring_dir(ns+1:end-ns,gdx);
	spring_t   = (-(n-ns):(n-ns))*dt+0.5;

	namedfigure(2,'Individual spring tides');
	clf();
	% plot the first 4 tides after september
	t0 = datenum('23-09-2013','dd-mm-yyyy');
	t__ = zeros(size(t));
	t__(spring_dx) = t(spring_dx);
	fdx = find(t__ > t0);
	l = {};
	colour = {'b', 'g', 'r', 'k'};
	% plot the first m spring tides
	%m = 10;
	m = size(spring_l,1);
	colour = 0.5*[sin(2*pi*(0:m-1)'/m) sin(2*pi*(0:m-1)'/m + pi/3) sin(2*pi*(0:m-1)'/m + 2*pi/3) ] + 0.5;
	for idx = 1:m
		plot(spring_t,spring_l(:,idx), 'color', colour(idx,:));
		%pdx = fdx(idx):fdx(idx)+1/dt;
		%plot(t(pdx)-t(fdx(idx)),level(pdx), 'color', colour(idx,:));
		%colour{idx});
		hold on;
		l{end+1} = datestr(t(fdx(idx)),'dd-mm-yyyy');
	end
	legend(l);
	xtick = get(gca,'xtick');
	set(gca,'xticklabel',datestr(xtick,'HH-MM'));

	% cluster tides
	c = kmeans(spring_l',ks);
%	c = kmeans(spring_u',ks)

	namedfigure(3,'Individual spring tides');
	spring_T = spring_T;
	for idx=1:30
		% 0:00
		t0 = floor(spring_T(round(size(spring_T,1)/2),18+idx));
		[t0_ fdx] = min((spring_T(:,18+idx)-t0).^2);
		subplot(5,6,idx);
		%plot(spring_T(ns+1:end-ns,18+idx), spring_(ns+1:end-ns,18+idx)); % TODO : time shift
		try
		plot(spring_T(fdx:fdx+2*ns,18+idx), spring_(fdx:fdx+2*ns,18+idx)); % TODO : time shift
		datetick(gca);
		xlim([spring_T(fdx,18+idx) spring_T(fdx+2*ns,18+idx)]);
		ylim([-1 1]);
		title(datestr(t0,'yyyy-mm-dd'));
		catch
		end
	end

	namedfigure(4,'Average spring tide');
	clf();
	for idx=1:ks
		subplot(1,ks,idx); cla();
		fdx_    = find(c == idx);
		spring_mean = mean(spring_l(:,fdx_),2);
		mean_u  = mean(spring_u(:,fdx_),2); %.*sign(sin(deg2rad(pi+spring_dir(:,fdx_)))),2);
		mean_ux = mean(spring_ux(:,fdx_),2);
		mean_uy = mean(spring_uy(:,fdx_),2);
		mean_u_ = sqrt(mean_ux.^2 + mean_uy.^2);
		mean_phi   = atan2(mean_uy,mean_ux);
		mdx        = max(mean_u_);
		spring_std = std(spring_l(:,fdx_),[],2);
		plot(spring_t,spring_mean,'k');
		hold on;
	%	plot(t_,mean_ux,'r');
	%	plot(t_,mean_uy,'g');
	%	plot(t_,mean_u_,'b');
%		plot(t_,mean_u,'r');
		plot(spring_t, spring_mean+spring_std,'color',[0.75 0.75 0.75]);
		plot(spring_t, spring_mean-spring_std,'color',[0.75 0.75 0.75]);
%		xlim([-(n-ns) n-ns]);
		ylim([-1 1]);
		%title(['Average Spring Tide']); % num2str(round(100*length(fdx_)/l0)) '%']);
		ylabel('Tidal amplitude (m)');
		xlabel('Hours since low water');
		xlim([spring_t(1) spring_t(end)])
		hourtick(6/24);
%		l = (-30:6:30)/24;
%		set(gca,'xtick',l);
%		datetick(gca);
%		t__ = l'* dt*24/60;
%		s = sprintf('%3d:%02d\n', [floor(t_/60) mod(t_,60)]')
%		s = {};
%		for idx=1:length(t__)%
%			s{idx} = sprintf('%3d:%02d\n', [floor(t__(idx)/60) mod(t__(idx),60)]');
%		end
%		set(gca,'xticklabel', s);
		grid on
		set(gca, 'xminorgrid','off');	
		set(gca, 'yminorgrid','off');	
	end % for idx

%	preparePrint();
%	print('-depsc','spring.eps');
%	system('epstopdf spring.eps');
	if (pflag)
		pdfprint(2,['img/',prefix,'spring']);
	end



	% Neap tide

	fdx = neap_dx;
%	neap_t  = t(fdx);
	n0  = round(13/24 * 1/dt);
	ns  = 2*n0; %13* 1/dt;
	ns  = (round(48/24*1/dt) - n0);
	n   = n0 + ns;
	idx = ones(size(fdx(2:end-1)))*(-n:n) + fdx(2:end-1)*ones(1,2*n+1);
	neap_l = level(idx');
	% shift the lowest amplitude to the centre
	c = round(0.5*size(neap_l,1));
	gdx = [];
	t_ = dt*(1:2*n+1)';
	T  = 0.5; % semidiurnal
	omega = 2*pi*t_/T;
	for idx=1:size(neap_l,2)
		% remove the phase shift so that all the tidal cycle overlap
		A = [sin(omega) cos(omega)];
		Y = A \ neap_l(:,idx);
		phi = atan2(Y(1),Y(2));
		phi = max(min(phi,phi+2*pi),phi-2*pi);
		delta = ((pi-phi)/(2*pi)*T)/dt;
		% monitor the maximum shift and exclude those which are shifted more than the plot range
		if (abs(delta) < ns)
			neap_l(:,idx) = circshift(neap_l(:,idx), round(delta));
			gdx(end+1) = idx;
		end
	end % for idx
	% only select good ones
	l0 = size(neap_l,2);
	neap_l = neap_l(ns+1:end-ns,:);
	l1 = size(neap_l,2);

	% cluster the tides into two groups (for verification)
	if (kn > 1)
		c = clusterdata(neap_l, kn);
	else
		c = ones(size(neap_l,2),1);
	end
	for kdx=1:kn
		nc(kdx) = sum(c==kdx);
	end
	[v sdx] = sort(nc);

	namedfigure(5,['Average Neap ' num2str(round(100*l1/l0)) '%']);
	clf();
	for idx=1:length(sdx)
		h(idx) = subplot(1,kn,idx);
		fdx_ = find(c == sdx(end-idx+1));
		l1 = length(fdx_);
		neap_mean = mean(neap_l(:,fdx_),2);
		neap_std  = std(neap_l(:,fdx_),[],2);
		neap_t = (-(n-ns):n-ns)*dt + 0.5;
		plot(neap_t, neap_mean,'k');
		hold on;
		plot(neap_t, neap_mean - neap_std, 'color', [0.75 0.75 0.75]);
		plot(neap_t, neap_mean + neap_std, 'color', [0.75 0.75 0.75]);
		xlim([neap_t(1),neap_t(end)]);
		ylim([-1 1]);
		ylabel('Tidal amplitude (m)');
		xlabel('Hours since low water');
		hourtick(6/24);

%		l = -30:6:30;
%		set(gca,'xtick',l);
%		%l = get(gca,'xtick');
%		t_ = l'* dt/(60);
%		%s = sprintf('%3d:%02d\n', [floor(t_/60) mod(t_,60)]')
%		s = {};
%		for idx=1:length(t_)
%			s{idx} = sprintf('%3d:%02d\n', [floor(t_(idx)/60) mod(t_(idx),60)]');
%		end
%		set(gca,'xticklabel', s);
		
		grid on;
		set(gca, 'xminorgrid','off');
		set(gca, 'yminorgrid','off');
	end

	if (pflag)
		pdfprint(5,['img/',prefix,'neap']);
	end

	namedfigure(6,'Average neap and spring tide');
	clf();
	ax = axes();
	pause(1);
	% actually, the time is hours before low water, we fix this by shifting back 51 min
	dt = (51/(1440));
	plot(ax,spring_t+dt,spring_mean);
	hold(ax,'on')
	plot(ax,neap_t+dt,neap_mean);
	xlim(ax,[neap_t(1)+dt,neap_t(end)]);
	hourtick(6/24); % ax
	x  = (0:4)/4;
	set(ax,'xtick',x,'xticklabel',{'0:00','6:00','12:00','18:00','24:00'});

	ylabel(ax,'Tidal amplitude (m)');
	xlabel(ax,'Hours since low water');
	grid(ax,'on');
	set(ax, 'xminorgrid','off');
	set(ax, 'yminorgrid','off');

	% plot spring range
	mu = mean(spring_mean);
	ma = max(spring_mean);
	mi = min(spring_mean);
	range =  ma-mi;
	errorbar(24.75/24,mu,mu-mi,ma-mu,'k','linewidth',2);
	text(25.125/24,2/3*ma,num2str(range,'%3.2f m'))
	text(25.125/24,ma+0.1,'Range');

	% plot neap range
	mu = mean(neap_mean);
	ma = max(neap_mean);
	mi = min(neap_mean);
	range =  ma-mi;
	errorbar(24.75/24,mu,mu-mi,ma-mu,'r','linewidth',2);
	text(25.125/24,1/4*ma,num2str(range,'%3.2f m'))
	set(gca,'box','off')
	ylim([-0.8 0.8])
	
%	preparePrint();
%	print('-depsc','neap.eps');
%	system('epstopdf neap.eps');
	if (pflag)
		pdfprint(6,['img/',prefix,'neap-and-spring']);
	end
	%legend(ax,'location','northwest','spring','neap');
	legend(ax,'spring','neap');

	if (pflag)
		pdfprint(6,['img/',prefix,'neap-and-spring-legend'],3,2/3);
	end
end

