% Fri  8 Dec 11:10:52 CET 2017
%% determine tidal discharge from water level for tidal wave
function q = z2q(x,z,hfun,omega)
	g = Constant.gravity;
	q = (-1i*omega)*cumintR(z,x);
	h = hfun(x);

	% determine integration constant (offset)
	% require E_pot = E_kin <=> d/dx Epot = d/dx Ekin
	% d/dx(Epot - Ekin) = 
	% TODO integrate if step length are not equal
	% delta_q = (cdiff(q_./h)./cdiff(x)) \ (cdiff(0.5*g*abs(z).^2 - 0.5*1./h.*abs(q_).^2 - 0.5*1./h*abs(delta_q).^2)./cdiff(x))
	q_      = q;
	% weights
	w       = (1./h.*real(q_));	% better than abs(w) and sqrt(g*abs(z).^2);
	a       = 0.5*sum(w.*1./h);
	b       = sum(w.*1./h.*real(q_));
	c       = 0.5*sum(w.*(1./h.*real(q_).^2 - g*real(z).^2)); % nore: in this row there was abs, not real
	delta_q = roots2([a b c]);
	% choose smaller
	[mv mdx] = min(abs(delta_q));
	delta_q = delta_q(mdx);
	q = q_+delta_q;
	% delta_q_ = (1./h.*real(q_)) \ (0.5*g*abs(z).^2 - 0.5./h.*abs(q_).^2 - 0.5./h.*abs(delta_q).^2)
end

