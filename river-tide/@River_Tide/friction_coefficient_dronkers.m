% Tue 21 Mar 08:47:01 CET 2017
% Karl Kastner, Berlin
%% friction coefficient according to Dronkers
%%
%% the coefficients are semi-autogenerated
%%
%% c.f. dronkers 1964
%% c.f. Cai 2016
%%
%% p = [p0,p1,p2,p3];
%% alpha = Ur/Ut = river velocity / tidal velocity amplitude = (umax+umin)/(umax-umin)
%%
%% function p = friction_coefficient_dronkers(alpha,order)
%%
function p = friction_coefficient_dronkers(obj,alpha)
	if (obj.issym)
		syms p0 p1 p2 p3 real
		p = [p0, p1, p2, p3];
		return;
	end

	% note: in cai 2016: alpha = -phi
	if (~issym(alpha))
		alpha = min(1,max(-1,alpha));
		pi_ = pi;
	else
		syms pi_
	end
	
	a = acos(alpha);
	alpha = a;
	switch (obj.opt.friction_order)
	case {0}
		p = - ( (3*sin(2*alpha))/2 + (cos(2*alpha) + 2)*(pi_/2 - alpha) );
	case {1}
		p = -[ pi_*(((3*sin(2*alpha))/2 + (cos(2*alpha) + 2)*(pi_/2 - alpha))/pi_ - (cos(alpha)*(sin(3*alpha)/3 + 3*sin(alpha) + cos(alpha)*(2*pi_ - 4*alpha)))/pi_)
                                                                                          sin(3*alpha)/3 + 3*sin(alpha) + cos(alpha)*(2*pi_ - 4*alpha) ].';
	case {2}
%		p = [  5/6*sin(2*alpha) - (cos(2*alpha) + 1)*(alpha - pi_/2) + sin(4*alpha)/12, ...
%                                   sin(3*alpha)/3 + 3*sin(alpha) - cos(alpha)*(4*alpha - 2*pi_), ...
%                                             -2*alpha + pi_ + 4/3*sin(2*alpha) - sin(4*alpha)/6
%			];
%		p = [ ...
%		 pi_^2*(((3*sin(2*alpha))/2 - (cos(2*alpha) + 2)*(alpha - pi_/2))/pi_ + (alpha - pi_/2 - (2*sin(2*alpha))/3 + sin(4*alpha)/12)/pi_)
 %                                                                pi_*(sin(3*alpha)/3 + 3*sin(alpha) - cos(alpha)*(4*alpha - 2*pi_))
  %                                                                   -2*pi_*(alpha - pi_/2 - (2*sin(2*alpha))/3 + sin(4*alpha)/12) ].';
		p = -[ pi_*(   ((3*sin(2*alpha))/2 + (cos(2*alpha) + 2).*(pi_/2 - alpha))/pi_ ...
			    + (cos(2*alpha).*(pi_/2 - alpha + (2*sin(2*alpha))/3 - sin(4*alpha)/12))/pi_ ...
                            - (cos(alpha).*(sin(3*alpha)/3 + 3*sin(alpha) + cos(alpha).*(2*pi_ - 4*alpha)))/pi_), ...
                       pi_*(   (sin(3*alpha)/3 + 3*sin(alpha) + cos(alpha).*(2*pi_ - 4*alpha))/pi_ ...
                            - (4*cos(alpha).*(pi_/2 - alpha + (2*sin(2*alpha))/3 - sin(4*alpha)/12))/pi_) ...
		       pi_ - 2*alpha + (4*sin(2*alpha))/3 - sin(4*alpha)/6];
%		p = -p/pi_;
	case {3}
		% note, this change of sign is not documented in dronkers, but necessary
		p = -[ -7/120*sin(2*a) + 1/24*sin(6*a) - 1/60*sin(8*a), ...
		       7/6*sin(a) - 7/30*sin(3*a) - 7/30*sin(5*a) + 1/10*sin(7*a), ...
		       pi_ - 2*a + 1/3*sin(2*a) + 19/30*sin(4*a) - 1/5*sin(6*a), ...
		       4/3*sin(a) - 2/3*sin(3*a) + 2/15*sin(5*a) ];	% ~ 4/3 a - 6/3 a + 2/3 a = 0 for high ur

	case {4}
	p = -[pi_.*(((3.*sin(2.*alpha))./2 + (cos(2.*alpha) + 2).*(pi_./2 - alpha))./pi_ + (cos(2.*alpha).*(pi_./2 - alpha + (2.*sin(2.*alpha))./3 - sin(4.*alpha)./12))./pi_ + (2.*(cos(alpha) - 2.*cos(2.*alpha).*cos(alpha)).*(sin(5.*alpha)./20 - sin(3.*alpha)./4 + sin(alpha)./2))./(3.*pi_) + ((cos(2.*alpha) + cos(alpha).*(2.*cos(alpha) - 4.*cos(2.*alpha).*cos(alpha))).*(sin(2.*alpha)./6 - (2.*sin(4.*alpha))./15 + sin(6.*alpha)./30))./(2.*pi_) - (cos(alpha).*(sin(3.*alpha)./3 + 3.*sin(alpha) + cos(alpha).*(2.*pi_ - 4.*alpha)))./pi_), ...
                                                                                               pi_.*((sin(3.*alpha)./3 + 3.*sin(alpha) + cos(alpha).*(2.*pi_ - 4.*alpha))./pi_ + (2.*(6.*cos(2.*alpha) + 3).*(sin(5.*alpha)./20 - sin(3.*alpha)./4 + sin(alpha)./2))./(3.*pi_) - (4.*cos(alpha).*(pi_./2 - alpha + (2.*sin(2.*alpha))./3 - sin(4.*alpha)./12))./pi_ + ((sin(2.*alpha)./6 - (2.*sin(4.*alpha))./15 + sin(6.*alpha)./30).*(cos(alpha).*(12.*cos(2.*alpha) + 6) - 6.*cos(alpha) + 4.*cos(2.*alpha).*cos(alpha)))./(2.*pi_)), ...
                                                                                                                                                                                                                                                          -pi_.*(((24.*cos(2.*alpha) + 16).*(sin(2.*alpha)./6 - (2.*sin(4.*alpha))./15 + sin(6.*alpha)./30))./(2.*pi_) - (2.*(pi_./2 - alpha + (2.*sin(2.*alpha))./3 - sin(4.*alpha)./12))./pi_ + (8.*cos(alpha).*(sin(5.*alpha)./20 - sin(3.*alpha)./4 + sin(alpha)./2))./pi_), ...
                                                                                                                                                                                                                                                                                                                                                  pi_.*((8.*(sin(5.*alpha)./20 - sin(3.*alpha)./4 + sin(alpha)./2))./(3.*pi_) + (16.*cos(alpha).*(sin(2.*alpha)./6 - (2.*sin(4.*alpha))./15 + sin(6.*alpha)./30))./pi_), ...
                                                                                                                                                                                                                                                                                                                                                                                                                                       (8.*sin(4.*alpha))./15 - (2.*sin(2.*alpha))./3 - (2.*sin(6.*alpha))./15];


	case {5}
 p = -[pi_.*(((3.*sin(2.*alpha))./2 + (cos(2.*alpha) + 2).*(pi_./2 - alpha))./pi_ + (cos(2.*alpha).*(pi_./2 - alpha + (2.*sin(2.*alpha))./3 - sin(4.*alpha)./12))./pi_ + (2.*(cos(alpha) - 2.*cos(2.*alpha).*cos(alpha)).*(sin(5.*alpha)./20 - sin(3.*alpha)./4 + sin(alpha)./2))./(3.*pi_) + ((cos(2.*alpha) + cos(alpha).*(2.*cos(alpha) - 4.*cos(2.*alpha).*cos(alpha))).*(sin(2.*alpha)./6 - (2.*sin(4.*alpha))./15 + sin(6.*alpha)./30))./(2.*pi_) + (2.*(sin(3.*alpha)./12 - sin(5.*alpha)./12 + sin(7.*alpha)./42).*(2.*cos(2.*alpha).*cos(alpha) - cos(alpha) + cos(alpha).*(2.*cos(2.*alpha) + 2.*cos(alpha).*(2.*cos(alpha) - 4.*cos(2.*alpha).*cos(alpha)))))./(5.*pi_) - (cos(alpha).*(sin(3.*alpha)./3 + 3.*sin(alpha) + cos(alpha).*(2.*pi_ - 4.*alpha)))./pi_), ...
 pi_.*((sin(3.*alpha)./3 + 3.*sin(alpha) + cos(alpha).*(2.*pi_ - 4.*alpha))./pi_ + (2.*(6.*cos(2.*alpha) + 3).*(sin(5.*alpha)./20 - sin(3.*alpha)./4 + sin(alpha)./2))./(3.*pi_) - (4.*cos(alpha).*(pi_./2 - alpha + (2.*sin(2.*alpha))./3 - sin(4.*alpha)./12))./pi_ - (2.*(sin(3.*alpha)./12 - sin(5.*alpha)./12 + sin(7.*alpha)./42).*(8.*cos(2.*alpha) + cos(alpha).*(4.*cos(alpha) - 8.*cos(2.*alpha).*cos(alpha)) - cos(alpha).*(2.*cos(alpha).*(12.*cos(2.*alpha) + 6) - 12.*cos(alpha) + 8.*cos(2.*alpha).*cos(alpha)) + 3))./(5.*pi_) + ((sin(2.*alpha)./6 - (2.*sin(4.*alpha))./15 + sin(6.*alpha)./30).*(cos(alpha).*(12.*cos(2.*alpha) + 6) - 6.*cos(alpha) + 4.*cos(2.*alpha).*cos(alpha)))./(2.*pi_)), ...
 -pi_.*((2.*(sin(3.*alpha)./12 - sin(5.*alpha)./12 + sin(7.*alpha)./42).*(cos(alpha).*(24.*cos(2.*alpha) + 12) - 24.*cos(alpha) + cos(alpha).*(48.*cos(2.*alpha) + 32) + 8.*cos(2.*alpha).*cos(alpha)))./(5.*pi_) - (2.*(pi_./2 - alpha + (2.*sin(2.*alpha))./3 - sin(4.*alpha)./12))./pi_ + ((24.*cos(2.*alpha) + 16).*(sin(2.*alpha)./6 - (2.*sin(4.*alpha))./15 + sin(6.*alpha)./30))./(2.*pi_) + (8.*cos(alpha).*(sin(5.*alpha)./20 - sin(3.*alpha)./4 + sin(alpha)./2))./pi_), ...
 pi_.*((8.*(sin(5.*alpha)./20 - sin(3.*alpha)./4 + sin(alpha)./2))./(3.*pi_) + (2.*(80.*cos(2.*alpha) + 60).*(sin(3.*alpha)./12 - sin(5.*alpha)./12 + sin(7.*alpha)./42))./(5.*pi_) + (16.*cos(alpha).*(sin(2.*alpha)./6 - (2.*sin(4.*alpha))./15 + sin(6.*alpha)./30))./pi_), ...
 -pi_.*((4.*(sin(2.*alpha)./6 - (2.*sin(4.*alpha))./15 + sin(6.*alpha)./30))./pi_ + (32.*cos(alpha).*(sin(3.*alpha)./12 - sin(5.*alpha)./12 + sin(7.*alpha)./42))./pi_), ...
 (8.*sin(3.*alpha))./15 - (8.*sin(5.*alpha))./15 + (16.*sin(7.*alpha))./105];
 case {6}
 p = -[                                                                                      pi_*(((3*sin(2*a))/2 + (cos(2*a) + 2)*(pi_/2 - a))/pi_ + (cos(2*a)*(pi_/2 - a + (2*sin(2*a))/3 - sin(4*a)/12))/pi_ + (2*(cos(a) - 2*cos(2*a)*cos(a))*(sin(5*a)/20 - sin(3*a)/4 + sin(a)/2))/(3*pi_) - ((sin(4*a)/20 - (2*sin(6*a))/35 + sin(8*a)/56)*(cos(2*a) - cos(a)*(4*cos(2*a)*cos(a) - 2*cos(a) + 2*cos(a)*(2*cos(2*a) + 2*cos(a)*(2*cos(a) - 4*cos(2*a)*cos(a)))) + cos(a)*(2*cos(a) - 4*cos(2*a)*cos(a))))/(3*pi_) + ((cos(2*a) + cos(a)*(2*cos(a) - 4*cos(2*a)*cos(a)))*(sin(2*a)/6 - (2*sin(4*a))/15 + sin(6*a)/30))/(2*pi_) + (2*(sin(3*a)/12 - sin(5*a)/12 + sin(7*a)/42)*(2*cos(2*a)*cos(a) - cos(a) + cos(a)*(2*cos(2*a) + 2*cos(a)*(2*cos(a) - 4*cos(2*a)*cos(a)))))/(5*pi_) - (cos(a)*(sin(3*a)/3 + 3*sin(a) + cos(a)*(2*pi_ - 4*a)))/pi_)
 pi_*((sin(3*a)/3 + 3*sin(a) + cos(a)*(2*pi_ - 4*a))/pi_ + (2*(6*cos(2*a) + 3)*(sin(5*a)/20 - sin(3*a)/4 + sin(a)/2))/(3*pi_) - (4*cos(a)*(pi_/2 - a + (2*sin(2*a))/3 - sin(4*a)/12))/pi_ - ((sin(4*a)/20 - (2*sin(6*a))/35 + sin(8*a)/56)*(cos(a)*(12*cos(2*a) + 6) - 8*cos(a) + 8*cos(2*a)*cos(a) + cos(a)*(16*cos(2*a) + 2*cos(a)*(4*cos(a) - 8*cos(2*a)*cos(a)) - 2*cos(a)*(2*cos(a)*(12*cos(2*a) + 6) - 12*cos(a) + 8*cos(2*a)*cos(a)) + 6) + cos(a)*(4*cos(2*a) + 4*cos(a)*(2*cos(a) - 4*cos(2*a)*cos(a)))))/(3*pi_) - (2*(sin(3*a)/12 - sin(5*a)/12 + sin(7*a)/42)*(8*cos(2*a) + cos(a)*(4*cos(a) - 8*cos(2*a)*cos(a)) - cos(a)*(2*cos(a)*(12*cos(2*a) + 6) - 12*cos(a) + 8*cos(2*a)*cos(a)) + 3))/(5*pi_) + ((sin(2*a)/6 - (2*sin(4*a))/15 + sin(6*a)/30)*(cos(a)*(12*cos(2*a) + 6) - 6*cos(a) + 4*cos(2*a)*cos(a)))/(2*pi_))
 -pi_*((2*(sin(3*a)/12 - sin(5*a)/12 + sin(7*a)/42)*(cos(a)*(24*cos(2*a) + 12) - 24*cos(a) + cos(a)*(48*cos(2*a) + 32) + 8*cos(2*a)*cos(a)))/(5*pi_) - (2*(pi_/2 - a + (2*sin(2*a))/3 - sin(4*a)/12))/pi_ + ((24*cos(2*a) + 16)*(sin(2*a)/6 - (2*sin(4*a))/15 + sin(6*a)/30))/(2*pi_) + (8*cos(a)*(sin(5*a)/20 - sin(3*a)/4 + sin(a)/2))/pi_ - ((sin(4*a)/20 - (2*sin(6*a))/35 + sin(8*a)/56)*(40*cos(2*a) + cos(a)*(8*cos(a) - 16*cos(2*a)*cos(a)) - cos(a)*(4*cos(a)*(12*cos(2*a) + 6) - 24*cos(a) + 16*cos(2*a)*cos(a)) - cos(a)*(2*cos(a)*(24*cos(2*a) + 12) - 48*cos(a) + 2*cos(a)*(48*cos(2*a) + 32) + 16*cos(2*a)*cos(a)) + 22))/(3*pi_))
 pi_*((8*(sin(5*a)/20 - sin(3*a)/4 + sin(a)/2))/(3*pi_) + ((sin(4*a)/20 - (2*sin(6*a))/35 + sin(8*a)/56)*(cos(a)*(48*cos(2*a) + 24) - 80*cos(a) + cos(a)*(96*cos(2*a) + 64) + cos(a)*(160*cos(2*a) + 120) + 16*cos(2*a)*cos(a)))/(3*pi_) + (2*(80*cos(2*a) + 60)*(sin(3*a)/12 - sin(5*a)/12 + sin(7*a)/42))/(5*pi_) + (16*cos(a)*(sin(2*a)/6 - (2*sin(4*a))/15 + sin(6*a)/30))/pi_)
 -pi_*((4*(sin(2*a)/6 - (2*sin(4*a))/15 + sin(6*a)/30))/pi_ + ((240*cos(2*a) + 192)*(sin(4*a)/20 - (2*sin(6*a))/35 + sin(8*a)/56))/(3*pi_) + (32*cos(a)*(sin(3*a)/12 - sin(5*a)/12 + sin(7*a)/42))/pi_)
 pi_*((32*(sin(3*a)/12 - sin(5*a)/12 + sin(7*a)/42))/(5*pi_) + (64*cos(a)*(sin(4*a)/20 - (2*sin(6*a))/35 + sin(8*a)/56))/pi_)
 (64*sin(6*a))/105 - (8*sin(4*a))/15 - (4*sin(8*a))/21];

 case {7}
 p = -[ pi_*(((3*sin(2*a))/2 + (cos(2*a) + 2)*(pi_/2 - a))/pi_ + (cos(2*a)*(pi_/2 - a + (2*sin(2*a))/3 - sin(4*a)/12))/pi_ + (2*(cos(a) - 2*cos(2*a)*cos(a))*(sin(5*a)/20 - sin(3*a)/4 + sin(a)/2))/(3*pi_) - ((sin(4*a)/20 - (2*sin(6*a))/35 + sin(8*a)/56)*(cos(2*a) - cos(a)*(4*cos(2*a)*cos(a) - 2*cos(a) + 2*cos(a)*(2*cos(2*a) + 2*cos(a)*(2*cos(a) - 4*cos(2*a)*cos(a)))) + cos(a)*(2*cos(a) - 4*cos(2*a)*cos(a))))/(3*pi_) + ((cos(2*a) + cos(a)*(2*cos(a) - 4*cos(2*a)*cos(a)))*(sin(2*a)/6 - (2*sin(4*a))/15 + sin(6*a)/30))/(2*pi_) - (2*(sin(5*a)/30 - sin(7*a)/24 + sin(9*a)/72)*(2*cos(2*a)*cos(a) - cos(a) + cos(a)*(2*cos(2*a) - 2*cos(a)*(4*cos(2*a)*cos(a) - 2*cos(a) + 2*cos(a)*(2*cos(2*a) + 2*cos(a)*(2*cos(a) - 4*cos(2*a)*cos(a)))) + 2*cos(a)*(2*cos(a) - 4*cos(2*a)*cos(a))) + cos(a)*(2*cos(2*a) + 2*cos(a)*(2*cos(a) - 4*cos(2*a)*cos(a)))))/(7*pi_) + (2*(sin(3*a)/12 - sin(5*a)/12 + sin(7*a)/42)*(2*cos(2*a)*cos(a) - cos(a) + cos(a)*(2*cos(2*a) + 2*cos(a)*(2*cos(a) - 4*cos(2*a)*cos(a)))))/(5*pi_) - (cos(a)*(sin(3*a)/3 + 3*sin(a) + cos(a)*(2*pi_ - 4*a)))/pi_),
 pi_*((sin(3*a)/3 + 3*sin(a) + cos(a)*(2*pi_ - 4*a))/pi_ + (2*(6*cos(2*a) + 3)*(sin(5*a)/20 - sin(3*a)/4 + sin(a)/2))/(3*pi_) + (2*(sin(5*a)/30 - sin(7*a)/24 + sin(9*a)/72)*(10*cos(2*a) - cos(a)*(2*cos(a)*(12*cos(2*a) + 6) - 16*cos(a) + 16*cos(2*a)*cos(a) + 2*cos(a)*(16*cos(2*a) + 2*cos(a)*(4*cos(a) - 8*cos(2*a)*cos(a)) - 2*cos(a)*(2*cos(a)*(12*cos(2*a) + 6) - 12*cos(a) + 8*cos(2*a)*cos(a)) + 6) + 2*cos(a)*(4*cos(2*a) + 4*cos(a)*(2*cos(a) - 4*cos(2*a)*cos(a)))) - cos(a)*(8*cos(2*a)*cos(a) - 4*cos(a) + 4*cos(a)*(2*cos(2*a) + 2*cos(a)*(2*cos(a) - 4*cos(2*a)*cos(a)))) + cos(a)*(8*cos(a) - 16*cos(2*a)*cos(a)) - cos(a)*(2*cos(a)*(12*cos(2*a) + 6) - 12*cos(a) + 8*cos(2*a)*cos(a)) + 3))/(7*pi_) - (4*cos(a)*(pi_/2 - a + (2*sin(2*a))/3 - sin(4*a)/12))/pi_ - ((sin(4*a)/20 - (2*sin(6*a))/35 + sin(8*a)/56)*(cos(a)*(12*cos(2*a) + 6) - 8*cos(a) + 8*cos(2*a)*cos(a) + cos(a)*(16*cos(2*a) + 2*cos(a)*(4*cos(a) - 8*cos(2*a)*cos(a)) - 2*cos(a)*(2*cos(a)*(12*cos(2*a) + 6) - 12*cos(a) + 8*cos(2*a)*cos(a)) + 6) + cos(a)*(4*cos(2*a) + 4*cos(a)*(2*cos(a) - 4*cos(2*a)*cos(a)))))/(3*pi_) - (2*(sin(3*a)/12 - sin(5*a)/12 + sin(7*a)/42)*(8*cos(2*a) + cos(a)*(4*cos(a) - 8*cos(2*a)*cos(a)) - cos(a)*(2*cos(a)*(12*cos(2*a) + 6) - 12*cos(a) + 8*cos(2*a)*cos(a)) + 3))/(5*pi_) + ((sin(2*a)/6 - (2*sin(4*a))/15 + sin(6*a)/30)*(cos(a)*(12*cos(2*a) + 6) - 6*cos(a) + 4*cos(2*a)*cos(a)))/(2*pi_)),
 pi_*((2*(pi_/2 - a + (2*sin(2*a))/3 - sin(4*a)/12))/pi_ - (2*(sin(3*a)/12 - sin(5*a)/12 + sin(7*a)/42)*(cos(a)*(24*cos(2*a) + 12) - 24*cos(a) + cos(a)*(48*cos(2*a) + 32) + 8*cos(2*a)*cos(a)))/(5*pi_) - ((24*cos(2*a) + 16)*(sin(2*a)/6 - (2*sin(4*a))/15 + sin(6*a)/30))/(2*pi_) - (8*cos(a)*(sin(5*a)/20 - sin(3*a)/4 + sin(a)/2))/pi_ + (2*(sin(5*a)/30 - sin(7*a)/24 + sin(9*a)/72)*(cos(a)*(48*cos(2*a) + 24) - 40*cos(a) + cos(a)*(48*cos(2*a) + 32) + 24*cos(2*a)*cos(a) + cos(a)*(80*cos(2*a) + 2*cos(a)*(8*cos(a) - 16*cos(2*a)*cos(a)) - 2*cos(a)*(4*cos(a)*(12*cos(2*a) + 6) - 24*cos(a) + 16*cos(2*a)*cos(a)) - 2*cos(a)*(2*cos(a)*(24*cos(2*a) + 12) - 48*cos(a) + 2*cos(a)*(48*cos(2*a) + 32) + 16*cos(2*a)*cos(a)) + 44) + cos(a)*(32*cos(2*a) + 4*cos(a)*(4*cos(a) - 8*cos(2*a)*cos(a)) - 4*cos(a)*(2*cos(a)*(12*cos(2*a) + 6) - 12*cos(a) + 8*cos(2*a)*cos(a)) + 12) + cos(a)*(8*cos(2*a) + 8*cos(a)*(2*cos(a) - 4*cos(2*a)*cos(a)))))/(7*pi_) + ((sin(4*a)/20 - (2*sin(6*a))/35 + sin(8*a)/56)*(40*cos(2*a) + cos(a)*(8*cos(a) - 16*cos(2*a)*cos(a)) - cos(a)*(4*cos(a)*(12*cos(2*a) + 6) - 24*cos(a) + 16*cos(2*a)*cos(a)) - cos(a)*(2*cos(a)*(24*cos(2*a) + 12) - 48*cos(a) + 2*cos(a)*(48*cos(2*a) + 32) + 16*cos(2*a)*cos(a)) + 22))/(3*pi_)),
 pi_*((8*(sin(5*a)/20 - sin(3*a)/4 + sin(a)/2))/(3*pi_) + ((sin(4*a)/20 - (2*sin(6*a))/35 + sin(8*a)/56)*(cos(a)*(48*cos(2*a) + 24) - 80*cos(a) + cos(a)*(96*cos(2*a) + 64) + cos(a)*(160*cos(2*a) + 120) + 16*cos(2*a)*cos(a)))/(3*pi_) + (2*(80*cos(2*a) + 60)*(sin(3*a)/12 - sin(5*a)/12 + sin(7*a)/42))/(5*pi_) + (16*cos(a)*(sin(2*a)/6 - (2*sin(4*a))/15 + sin(6*a)/30))/pi_ - (2*(sin(5*a)/30 - sin(7*a)/24 + sin(9*a)/72)*(160*cos(2*a) - cos(a)*(2*cos(a)*(48*cos(2*a) + 24) - 160*cos(a) + 2*cos(a)*(96*cos(2*a) + 64) + 2*cos(a)*(160*cos(2*a) + 120) + 32*cos(2*a)*cos(a)) + cos(a)*(16*cos(a) - 32*cos(2*a)*cos(a)) - cos(a)*(8*cos(a)*(12*cos(2*a) + 6) - 48*cos(a) + 32*cos(2*a)*cos(a)) - cos(a)*(4*cos(a)*(24*cos(2*a) + 12) - 96*cos(a) + 4*cos(a)*(48*cos(2*a) + 32) + 32*cos(2*a)*cos(a)) + 104))/(7*pi_)),
 -pi_*((4*(sin(2*a)/6 - (2*sin(4*a))/15 + sin(6*a)/30))/pi_ + ((240*cos(2*a) + 192)*(sin(4*a)/20 - (2*sin(6*a))/35 + sin(8*a)/56))/(3*pi_) + (2*(sin(5*a)/30 - sin(7*a)/24 + sin(9*a)/72)*(cos(a)*(96*cos(2*a) + 48) - 240*cos(a) + cos(a)*(192*cos(2*a) + 128) + cos(a)*(320*cos(2*a) + 240) + cos(a)*(480*cos(2*a) + 384) + 32*cos(2*a)*cos(a)))/(7*pi_) + (32*cos(a)*(sin(3*a)/12 - sin(5*a)/12 + sin(7*a)/42))/pi_),
 pi_*((32*(sin(3*a)/12 - sin(5*a)/12 + sin(7*a)/42))/(5*pi_) + (2*(672*cos(2*a) + 560)*(sin(5*a)/30 - sin(7*a)/24 + sin(9*a)/72))/(7*pi_) + (64*cos(a)*(sin(4*a)/20 - (2*sin(6*a))/35 + sin(8*a)/56))/pi_),
 -pi_*((32*(sin(4*a)/20 - (2*sin(6*a))/35 + sin(8*a)/56))/(3*pi_) + (128*cos(a)*(sin(5*a)/30 - sin(7*a)/24 + sin(9*a)/72))/pi_),
 (64*sin(5*a))/105 - (16*sin(7*a))/21 + (16*sin(9*a))/63];
	otherwise
		error('Not yet implemented, use derivation script and copy and paste equation');
	end
end % River_Tide/friction_coefficient_dronkers

