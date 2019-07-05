% 2014-09-17 14:10:25.650206633 +0200
% Karl Kastner, Berlin

% TODO could be extended to slicing (consecutive widow)
function [amplitude, phase] = constituents(t,f,p,n,trendflag)
	if (nargin() > 4 & trendflag)
		amplitude = NaN(length(f),p+2);
		A = ones(n,2*p+2);
	else
		amplitude = NaN(length(f),p+1);
		A = ones(n,2*p+1);
	end
	phase = NaN(length(f),p);
	if (nargin() > 3)
		% moving window
		for idx=1:length(f)-n
			t_ = t(idx:idx+n-1);
			if (nargin() > 4 & trendflag)
				% trend
				A(:,2) = t_-mean(t_);
				for pdx=1:p
					A(:,3+2*(pdx-1)) = sin(pdx*2*pi*t_);
					A(:,4+2*(pdx-1)) = cos(pdx*2*pi*t_);
				end
				c = A \ f(idx:idx+n-1);
				c = c(:)';
				% TODO vectorise at the end from C
				amplitude(idx,:) = [c(1) c(2) sqrt(c(3:2:end).^2 + c(4:2:end).^2)];
				phase(idx,:)     = atan2(c(4:2:end),c(3:2:end));
			else
				for pdx=1:p
					A(:,2+2*(pdx-1)) = sin(pdx*2*pi*t_);
					A(:,3+2*(pdx-1)) = cos(pdx*2*pi*t_);
				end
				c = A \ f(idx:idx+n-1);
				c = c(:)';
				amplitude(idx,:) = [c(1) sqrt(c(2:2:end-1).^2 + c(3:2:end-1).^2)];
				phase(idx,:) = atan2(c(3:2:end),c(2:2:end));
			end
		end
	else
		% over all samples
		A = [ones(size(t)) sin(2*pi*t) cos(2*pi*t)];
		c = A \ f;
		amplitude = [c(1) sqrt(c(2).^2 + c(3).^2)];
		phase = atan2(c(3:2:end),c(2:2:end));
	end
end % constituents

