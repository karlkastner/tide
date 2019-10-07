% 2016-07-12 15:06:16.443700194 +0200
% Karl Kastner, Berlin
%
%% predict tidal phase 
function [phi, k, obj] = predict_phase(obj,time,phi0,Ur,dUr_dt,x,x0,r0)

	c  = obj.c.phase;

	fk  = obj.fk;

	nt  = length(time);
	nd  = size(phi0,2);
	k   = NaN(nt,nd);

	% for each constituent
	for idx=1:size(phi0,2)
		% wave number
		k(:,idx) = fk(c(:,idx),idx,x,x0,Ur,r0,zeros(size(Ur)));
	end % for idx
	% phase
	phi = wrapTo2Pi(phi0 - k*x);
end % River_Tide_Empirical/predict_phase

