% Sat Jul 13 10:29:28 UTC 2013
% Karl Kästner, Berlin

	name = [ROOTFOLDER,'current/tide/table-pontianak/'];

	% load tide predicted by TPXO
	dat = load([base,filesep,'tide_.out']);

	% load tidal constituents as given by TPXO
	hc = load([base,filesep,'hc_.out']);
	A_orig = hc(:,1);
	P_orig = hc(:,2);

	% datenum gives the number of days since year 0
	t = 24*3600*datenum(dat(:,5), dat(:,3), dat(:,4), dat(:,6), dat(:,7), dat(:,8));
	level_orig = dat(:,9);

	% choose a 14 day calibration period
	% for a longer calibration period, say a complete metonic cycle
	% the amplitude if calculated accurately, but the phase is nowhere
	% really correct, due to the missing nodal corrections
	n = round(14*(24*3600/(t(2)-t(1))));
%	n = size(t,1);

	% calculate the tidal constituents by linear regression
	[A P T legend Tw res] = tidal_harmonic_analysis(t(1:n),level_orig(1:n));
%	T = T*3600;
%	Tw = Tw*3600;

	% regression residual
	nres = norm(res)

	% predict the tide with the extracted constituents
	level = sum(sin(t*(2*pi./T)' + ones(size(t))*P')*diag(A),2);

	% calculate the relative error over calibration periode
	norm(level(1:n)-level_orig(1:n))/norm(level_orig(1:n))

	% calculate error for one metonic cycle (18.6 years)
	norm(level-level_orig)/norm(level_orig)

	tide = t_tide_interface(name);

	P_orig = P_orig.*(P_orig > 0) + (P_orig+360).*(P_orig < 0);
	%P_orig = P_orig.*(P_orig > 0) + (P_orig+2*pi).*(P_orig < 0);
	%P_orig = rad2deg(P_orig);
	P = P.*(P > 0) + (P+2*pi).*(P < 0);
	P      = rad2deg(P);
	[A_orig A tide.value(:,1) P_orig P tide.value(:,3)]

	% monitor the error over time
	% this looks like a phase error, but is presumably due to the missing
	% nodal corrections the linear regression
	err = (level-level_orig);

	figure(1);	
	subplot(2,2,1)
	plot(t,err);
	subplot(2,2,2);
	bar([A_orig A tide.value(:,1)]);
	subplot(2,2,3);
	bar([P_orig P tide.value(:,3)]);
	subplot(2,2,4);
	plot(t(1:n),[level(1:n) level_orig(1:n)])
	figure(2);
	% 14-day overlay
	T = 14.75;
	T = 13.958333333372138;
	T = 13.649447278911564;
	T = 27.312831982108044
	n = round(T*24*3600 * 1/(t(2)-t(1)));
	m = 6;

	l = reshape(level_orig(1:m*n),n,m);
	figure(2);
	plot((1:n)*T/n,l);

