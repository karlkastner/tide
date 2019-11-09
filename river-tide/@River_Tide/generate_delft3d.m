%Tue  5 Nov 20:05:15 +08 2019
function generate_delft3d(obj,name)
%	nn = 2;
	w = obj.width();
	x = obj.x;
%	nn = [200,8];
%	x  = linspace(0,1e6,200).';
%	dS = x(2)-x(1);
%	y = zeros(size(x));
%	w = bfun(x);
	s = StructuredMesh();
	s.X = [x,x];
	% grid orientation
	s.Y = w*[-0.5,0.5];
	%s.Y = -s.Y;
	s.Z = obj.zb*[1,1]; 
	%s.generate_from_centreline(x,y,w,dS,nn(2));
	s.plot();
	%X = s.X;
	%s.Z = Z;
	%w  =
	%zb =
%	s.generate_1D()
	if (nargin()>1)
		s.export_delft3d_grd('mat/yangtze-200-008.grd');
		s.export_delft3d_dep('mat/yangtze-200-008.dep');
	end
end
