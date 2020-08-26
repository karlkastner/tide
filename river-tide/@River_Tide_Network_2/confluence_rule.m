function Qs_out = confluence_rule(Qs)
	% flow out of branches remains the same
	% flow into the downstream branch is the combined flow
	Qs_out = -(Qs(2)+Qs(3));
end

