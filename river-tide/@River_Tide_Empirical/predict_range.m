% 2016-07-12 15:06:16.443700194 +0200
% Karl Kastner, Berlin
%
%% predict the tidal range
function R = predict_range(obj,time,R0,Ur,dUr_dt,x,x0)

	%[fD fr] = rt_model(obj.model);

	c = obj.c.range;

%	nt = length(time);
%	r  = NaN(nt,nd);
%	re  = NaN(nt,nd);
%	R  = NaN(nt,nd);
	% for each constituent (species)
	
	% damping modulus
	% r(:,idx) = obj.fr(c(:,idx),idx,x,x0,Ur,r0,D(:,1));

	% amplitude and damping modulus
	% D(:,idx) = obj.fD(c(:,idx),idx,D0(:,idx),x,x0,Ur,r0,D(:,1));
	R = obj.fD(c,1,R0,x,x0,Ur,R0,[]);

end % River_Tide_Empirical/predict_range

