% Thu 20 Aug 17:46:30 +08 2020
% excluding zero flow, there are 8 conditions, of which 6 are feasible
% 3 confluences (into a, b or c)
% 3 bifurcations (from a, b, c)
% infeasible (sink from all, source into all)
function Qs = sediment_division(obj,Q0,Qt,zt,Qs,w,h,Cz,d_mm);

	dir       = sign(Q0);
	indicator = rvec(dir)*[1;2;4];


	% TODO, recover cases with zero flow in any channel
	if (any(0 == dir))
		error('not yet implemented');
	end

	switch (indicator)
	case {-5}  %   1, -2, -4 (a feeds)
		       % natural order 1,2,3 no permutation
		       Qs = obj.division_rule(Q0,Qt,zt,Qs,w,h,Cz,d_mm);
	case {-3}  %  -1,  2, -4 (b feeds)
		       id = [2,1,3];
		       % permute, compute and permute back
		       Qs(id) = obj.division_rule(Q0(id),Qt(id),zt(id),Qs(id),w(id),h(id),Cz(id),d_mm);
	case {-1}  %   1,  2, -4 (c drains)
		       Qs(3) = obj.confluence_rule(Qs(1),Qs(2));
	case {+1}  %  -1, -2,  4 (c feeds)
		       id = [3,1,2];
		       % permute, compute and permute back
		       Qs(id) = obj.division_rule(Q0(id),Qt(id),zt(id),Qs(id),w(id),h(id),Cz(id),d_mm);
	case {+3}  %   1, -2,  4 (b drains)
		       Qs(2) = obj.confluence_rule(Qs(1),Qs(3));
	case {+5}  %  -1,  2,  4 (a drains)
		       Qs(1) = obj.confluence_rule(Qs(2),Qs(3));
	otherwise % such as +7, -7
		Qs = [0,0,0];
		warning('flow at junction is inconsistent');
		%error('flow at junction is inconsistent');
	end % switch sum dir
end % sediment_division

% [ps1,r1,ps2,r2] = sediment_division_wang(Q1,Q2,w1,w2,k)
%	% wang
%	if (bifurcation)
%		Qs_a = Qs_0*(Qs_i./Qs_)^k;
%		Qs_b = Qs_0*
%	else % confluence
%
%	end

