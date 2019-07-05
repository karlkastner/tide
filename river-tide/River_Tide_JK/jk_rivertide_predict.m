% 2016-07-07 19:19:42.784764868 +0200
% Karl Kastner, Berlin
%% predict river tide by the method of jay and kukulka
function [mag phi val] = jk_rivertide_predict(z0,phi0,Ur,h,d)
	m 	  = size(z0,2);
	n         = size(z0,1);
	Ur        = abs(Ur);

	mag = [];
	phi = [];
	for idx=1:m
		% regression/prediction matrix

		% this somehow contradicts 14a in J&F 1997:
		%	D = D0 exp( sqrt(A0/A) - sqrt(2/T c_D/omega |ur|/H_0)x ), mu = 0 (mean), i (species)
		% also c.f. godin 
		%	real(k) = 27/C^2 u_0^2 / (2 g H^2)
		X              = [ones(n,1),Ur,1./sqrt(Ur).*abs(z0(:,idx)./h).^2];

		% magnitude
		mag(:,idx)     = z0(:,idx).*exp(X*d(:,1,idx));

		% phase
		phi(:,idx)     = phi0(:,idx) + (X*d(:,2,idx));
	end

	% complex representation
	val = mag.*exp(1i*phi);
end

