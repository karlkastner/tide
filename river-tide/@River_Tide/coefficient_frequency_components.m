% Fri 16 Oct 17:29:15 +08 2020
% (a exp(k) + a' exp(-k)) (b exp(l) + b' exp(-l))
% = a b exp(k-l) + a' b exp(l-k)) + a b' exp(k-l) + a'b' exp(-k-l)
% reorder : 
% = (a b exp(k-l) a'b' exp(-k-l) ) + (a' b exp(l-k)) + a b' exp(k-l)
function f = coefficient_frequency_components(obj,f,c,Q,ddx,dt)
	nf = size(Q,2)-1;
	% TODO this has to be modified for non-integer multiples
	omega = obj.omega(1);
	
	% for the zero, first, and second order derivative
	% (c_0 + ... + c_k) D^ddx (Q_0 + ... Q_k)
%	for ddx=1:3
		% ode-coefficients for frequency components
		% for freqeuency components of coefficients
		for kc = 1:nf+1
			% frequency compoenent of the Q-derivative
			% kq<kc
			%for kq = 1:kc
			for kq = 1:nf+1
				% plus : a*b*exp(io(kc+kl-1))
				k = kc + kq - 1;
				if (k <= nf+1)
				    if (k == kq)
					% lhs, factor Q indirectly from unknowns
					f(:,ddx,k) = f(:,ddx,k) + 0.5*(1i*(kq-1)*omega).^dt*c(:,kc);
				    else
					% forcing (rhs)
					% note that his could be moved to the lhs by extending c and the matrix
					f(:,4,k) = f(:,4,k) + 0.5*c(:,kc).*(1i*(kq-1)*omega).^dt.*Q(:,kq);
				    end
				end % k <= nf
	
				% minus : a*conj(b)*exp(k-l) <=> conj(a)*conj(b)*exp(l-k)
				% note kc = kq contributes to the mwl equation,
				% which is computed separately
				if (kc ~= kq)
					if (kc<kq)
						% forcing (rhs)
						k = kq - kc + 1;
						if (k == kq) % <=> kc == 1 (mean component)
							f(:,ddx,k) = f(:,ddx,k) + 0.5*conj(c(:,kc)).*(1i*(kq-1)*omega).^dt;
						else
							f(:,4,k) = f(:,4,k) + 0.5*conj(c(:,kc)).*(1i*(kq-1)*omega).^dt.*Q(:,kq);
						end
					else
						% forcing (rhs)
						k = kc - kq + 1;
						f(:,4,k) = f(:,4,k) + 0.5*c(:,kc).*conj((1i*(kq-1)*omega).^dt.*Q(:,kq));
					end
				end % if kc ~= kq
			end % for kq, kq<kc
		end % for kc
%	end % for ddx
end % coefficient_frequency_components

