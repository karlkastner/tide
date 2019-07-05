% Fri 22 Feb 14:31:01 CET 2019
% Karl Kastner, Berlin
%
%% plot surface elevation amplitude
%
function plot_water_level_amplitude(obj)
	namedfigure(1,'Amplitude');
	clf();
	namedfigure(2,'Phase');
	clf();
	figure(10);
	clf
	figure(20);
	clf
	for idx=1:obj.nc
		[x,zm,phi,s,c] = obj.water_level_amplitude(idx);
		phi = unwrap(phi);
		%plot3(x,idx*ones(size(x));	
		figure(1);
		subplot(2,2,1);
		scatter3(x,idx*ones(size(x)),zm,[],zm,'.');
		hold on;
		view([0,90]);
		colorbar();
		figure(2)
		subplot(2,2,1);
		plot(x,zm);
		hold on

		figure(1);
		subplot(2,2,2);
		scatter3(x,idx*ones(size(x)),phi,[],phi,'.');
		hold on;
		view([0,90]);
		figure(2);
		subplot(2,2,2);
		plot(x,phi);
		hold on

		figure(1);
		subplot(2,2,3);
		scatter3(x,idx*ones(size(x)),s,[],s,'.');
		hold on;
		view([0,90]);
		figure(2)
		subplot(2,2,3);
		plot(x,s);
		hold on

		figure(1);
		subplot(2,2,4);
		scatter3(x,idx*ones(size(x)),c,[],c,'.');
		hold on;
		view([0,90]);
		figure(2)
		subplot(2,2,4);
		plot(x,c);
		hold on
	end % for idx
end

