% 2016-07-12 15:06:16.443700194 +0200
% Karl Kastner, Berlin
%
%% predict the oscillatory components
function [D, r, re, obj] = predict_amplitude(obj,time,D0,Ur,dUr_dt,x,x0,r0)

	%[fD fr] = rt_model(obj.model);

	c = obj.c.amplitude;

	nt = length(time);
	nd = size(D0,2);
	r  = NaN(nt,nd);
	re  = NaN(nt,nd);
	D  = NaN(nt,nd);
	% for each constituent (species)
	for idx=1:size(D0,2)
		% damping modulus
		% r(:,idx) = obj.fr(c(:,idx),idx,x,x0,Ur,r0,D(:,1));

		% amplitude and damping modulus
		% D(:,idx) = obj.fD(c(:,idx),idx,D0(:,idx),x,x0,Ur,r0,D(:,1));
		[D(:,idx) r(:,idx) re(:,idx)] = obj.fD(c(:,idx),idx,D0(:,idx),x,x0,Ur,r0,D(:,1));
	end % for idx
end % River_Tide_Empirical/predict_amplitude

