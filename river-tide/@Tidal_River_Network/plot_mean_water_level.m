% Mon  4 Mar 11:14:30 CET 2019
%
%% plot tidally averaged water level
function plot_mean_water_level(obj)
	namedfigure(1,'MWL');
	clf();
	namedfigure(2,'Mean Discharge');
	clf();
	figure(10);
	clf
	figure(20);
	clf
	for idx=1:obj.nc
		[x,z0,Q0] = obj.mean_water_level(idx);

		figure(1);
		subplot(2,2,1);
		scatter3(x,z0,z0,[],z0,'.');
		hold on;
		view([0,90]);
		colorbar();
		
		subplot(2,2,2);
		plot(x,z0);
		hold on;

		figure(2)
		subplot(2,2,1);
		plot(x,Q0);
		hold on
	end % for idx
end

