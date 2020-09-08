% Fri  1 Dec 14:21:20 CET 2017
% Karl Kastner, Berlin
%
%% analytic expression of the wave number of river tides
%%
%% valid for both tidally, river dominated and low friction conditions
%% and converging channels
%%
%% k10    : complex wave number for k and z in a reach with constant width and bed slope
%% im(k)  : damping modulus (rate of amplitude change)
%% re(k)  : actual wave number (rate of phase change)
%%
%% kq     : wave number for Q for a reach with changing width and depth
%% kz     : wave number for z for a reach with changing width and depth
%%
%% c.f. derive_wave_number
%%
% function [k10, kq, kz, c] = wave_number_analytic(obj,Q0,w,h,cd,omega,az1,Qt,dh_dx,dw_dx)
function [k10, kq, kz, c] = wave_number_analytic(obj,Q0,Qt,h,dh_dx,w,dw_dx,cd)
	g = obj.g;
	pi_ = obj.pi;

	% specific discharge
	q0 = Q0./w;
	qt = Qt./w;
	Q2 = 0;

	% discharge range
	Qhr = abs(Qt);

	alpha = q0./qt;
	p     = -obj.friction_coefficient_dronkers(alpha);
	p1    = p(:,2);
	p2    = p(:,3);
	if (0)
		c  = [ 1;
		      -1/w.*dw_dx;
		       omega^2/(g*h) - 1i*(omega*cd)/(pi*w*g*h^3)*(p1*Qt + p2*Q0)];
		k0 = roots(c);
	end
	flag.aa=false;
	flag.gh=false;
	flag.oh=false;

	% characteristic polynomial
	% [f, F3]  = odefunk(obj, k, Q, QQ, Qhr, h0, dh0_dx, dz0_dx, w0, dw0_dx, cD, c, D1_dx)
	c = obj.odefunk(1, Q0, Qt, Qhr, h, dh_dx, 0, w, dw_dx, cd, [], []);

	% roots of the characteristic polynomial
	r = roots2(c(:,1:3));
	
	% wave number
	k10 = 1i*r;

	Q = 1/2*(p1.*Qt + 2*p2.*Q0);

	% simplified for q
	kq = [dw_dx.*((1.0./2.0)./w-sqrt(pi).*sqrt(g).*h.^(3.0./2.0).*w.*1.0./sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*(1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(3.0./2.0).*1.0./w.^2.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*5.0e-1i-(1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(3.0./2.0).*1.0./sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*(Q.*cd.*omega.*8.0i+pi.*h.^2.*omega.^2.*w.*8.0).*2.5e-1i)./w).*1i)-(1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(3.0./2.0).*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*5.0e-1i)./w-sqrt(pi).*dh_dx.*sqrt(g).*h.^(3.0./2.0).*w.*((1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(5.0./2.0).*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*7.5e-1i)./w-sqrt(pi).*1.0./sqrt(g).*1.0./sqrt(h).*omega.^2.*w.*1.0./sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*2.0i).*1.0./sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*1i, ...
	dw_dx.*((1.0./2.0)./w-sqrt(pi).*sqrt(g).*h.^(3.0./2.0).*w.*1.0./sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*(1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(3.0./2.0).*1.0./w.^2.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*5.0e-1i-(1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(3.0./2.0).*1.0./sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*(Q.*cd.*omega.*8.0i+pi.*h.^2.*omega.^2.*w.*8.0).*2.5e-1i)./w).*1i)+(1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(3.0./2.0).*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*5.0e-1i)./w-sqrt(pi).*dh_dx.*sqrt(g).*h.^(3.0./2.0).*w.*((1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(5.0./2.0).*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*7.5e-1i)./w-sqrt(pi).*1.0./sqrt(g).*1.0./sqrt(h).*omega.^2.*w.*1.0./sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*2.0i).*1.0./sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*1i];

	% not expanded
kz = [(sqrt(pi).*sqrt(g).*h.^(3.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3.*2.0+Q.*cd.*omega.*w.*8.0i).*(Q.^3.*cd.^3.*omega.^3.*w.^2.*2.56e2i-pi.^3.*h.^6.*omega.^6.*w.^5.*3.2e1-pi.^2.*Q.*cd.*h.^4.*omega.^5.*w.^4.*1.92e2i+pi.^(3.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(7.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*3.0i+pi.^(5.0./2.0).*dh_dx.*dw_dx.^4.*g.^(5.0./2.0).*h.^(1.3e1./2.0).*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*3.0i+pi.*Q.^2.*cd.^2.*h.^2.*omega.^4.*w.^3.*3.84e2-pi.^3.*dh_dx.^2.*g.*h.^5.*omega.^4.*w.^5.*6.0+pi.^3.*dw_dx.^2.*g.*h.^7.*omega.^4.*w.^3.*1.6e1+pi.^3.*dw_dx.^4.*g.^2.*h.^8.*omega.^2.*w.*8.0-pi.^(5.0./2.0).*dw_dx.^3.*g.^(3.0./2.0).*h.^(1.3e1./2.0).*omega.^2.*w.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0i-pi.^3.*dh_dx.*dw_dx.^3.*g.^2.*h.^7.*omega.^2.*w.^2.*1.0e1+pi.^2.*Q.*cd.*dw_dx.^4.*g.^2.*h.^6.*omega.*4.0i+pi.^3.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^6.*omega.^2.*w.^3.*4.0+pi.*Q.^2.*cd.^2.*dh_dx.^2.*g.*h.*omega.^2.*w.^3.*2.4e1-pi.^2.*Q.*cd.*dh_dx.^2.*g.*h.^3.*omega.^3.*w.^4.*1.04e2i-pi.*Q.^2.*cd.^2.*dw_dx.^2.*g.*h.^3.*omega.^2.*w.*4.0e1+pi.^2.*Q.*cd.*dw_dx.^2.*g.*h.^5.*omega.^3.*w.^2.*3.2e1i-pi.^(5.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(1.1e1./2.0).*omega.^2.*w.^2.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0i+pi.^(3.0./2.0).*Q.*cd.*dw_dx.^3.*g.^(3.0./2.0).*h.^(9.0./2.0).*omega.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0+pi.^2.*Q.*cd.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^4.*omega.*w.^2.*4.8e1i-pi.^2.*Q.*cd.*dh_dx.*dw_dx.*g.*h.^4.*omega.^3.*w.^3.*5.6e1i-pi.^2.*Q.*cd.*dh_dx.*dw_dx.^3.*g.^2.*h.^5.*omega.*w.*3.6e1i-pi.*Q.^2.*cd.^2.*dh_dx.*dw_dx.*g.*h.^2.*omega.^2.*w.^2.*4.8e1).*(-1.0./4.0)-pi.*g.*h.^2.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(pi.*dw_dx.^3.*g.*h.^4.*1i-pi.*dh_dx.*h.^2.*omega.^2.*w.^3.*1i-pi.*dw_dx.*h.^3.*omega.^2.*w.^2.*2.0i+Q.*cd.*dh_dx.*omega.*w.^2.*6.0+Q.*cd.*dw_dx.*h.*omega.*w.*6.0).*(Q.^3.*cd.^3.*omega.^3.*w.^2.*2.56e2i-pi.^3.*h.^6.*omega.^6.*w.^5.*3.2e1-pi.^2.*Q.*cd.*h.^4.*omega.^5.*w.^4.*1.92e2i+pi.^(3.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(7.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*3.0i+pi.^(5.0./2.0).*dh_dx.*dw_dx.^4.*g.^(5.0./2.0).*h.^(1.3e1./2.0).*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*3.0i+pi.*Q.^2.*cd.^2.*h.^2.*omega.^4.*w.^3.*3.84e2-pi.^3.*dh_dx.^2.*g.*h.^5.*omega.^4.*w.^5.*6.0+pi.^3.*dw_dx.^2.*g.*h.^7.*omega.^4.*w.^3.*1.6e1+pi.^3.*dw_dx.^4.*g.^2.*h.^8.*omega.^2.*w.*8.0-pi.^(5.0./2.0).*dw_dx.^3.*g.^(3.0./2.0).*h.^(1.3e1./2.0).*omega.^2.*w.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0i-pi.^3.*dh_dx.*dw_dx.^3.*g.^2.*h.^7.*omega.^2.*w.^2.*1.0e1+pi.^2.*Q.*cd.*dw_dx.^4.*g.^2.*h.^6.*omega.*4.0i+pi.^3.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^6.*omega.^2.*w.^3.*4.0+pi.*Q.^2.*cd.^2.*dh_dx.^2.*g.*h.*omega.^2.*w.^3.*2.4e1-pi.^2.*Q.*cd.*dh_dx.^2.*g.*h.^3.*omega.^3.*w.^4.*1.04e2i-pi.*Q.^2.*cd.^2.*dw_dx.^2.*g.*h.^3.*omega.^2.*w.*4.0e1+pi.^2.*Q.*cd.*dw_dx.^2.*g.*h.^5.*omega.^3.*w.^2.*3.2e1i-pi.^(5.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(1.1e1./2.0).*omega.^2.*w.^2.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0i+pi.^(3.0./2.0).*Q.*cd.*dw_dx.^3.*g.^(3.0./2.0).*h.^(9.0./2.0).*omega.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0+pi.^2.*Q.*cd.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^4.*omega.*w.^2.*4.8e1i-pi.^2.*Q.*cd.*dh_dx.*dw_dx.*g.*h.^4.*omega.^3.*w.^3.*5.6e1i-pi.^2.*Q.*cd.*dh_dx.*dw_dx.^3.*g.^2.*h.^5.*omega.*w.*3.6e1i-pi.*Q.^2.*cd.^2.*dh_dx.*dw_dx.*g.*h.^2.*omega.^2.*w.^2.*4.8e1).*(1.0./2.0))./(pi.*g.*h.^3.*(pi.*h.^2.*omega.^2.*w.^2.*2.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*4.0i).^2.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^3-pi.^2.*g.^2.*h.^4.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^2.*(pi.*dw_dx.^3.*g.*h.^4.*1i-pi.*dh_dx.*h.^2.*omega.^2.*w.^3.*1i-pi.*dw_dx.*h.^3.*omega.^2.*w.^2.*2.0i+Q.*cd.*dh_dx.*omega.*w.^2.*6.0+Q.*cd.*dw_dx.*h.*omega.*w.*6.0).^2),...
(sqrt(pi).*sqrt(g).*h.^(3.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3.*2.0+Q.*cd.*omega.*w.*8.0i).*(Q.^3.*cd.^3.*omega.^3.*w.^2.*-2.56e2i+pi.^3.*h.^6.*omega.^6.*w.^5.*3.2e1+pi.^2.*Q.*cd.*h.^4.*omega.^5.*w.^4.*1.92e2i+pi.^(3.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(7.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*3.0i+pi.^(5.0./2.0).*dh_dx.*dw_dx.^4.*g.^(5.0./2.0).*h.^(1.3e1./2.0).*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*3.0i-pi.*Q.^2.*cd.^2.*h.^2.*omega.^4.*w.^3.*3.84e2+pi.^3.*dh_dx.^2.*g.*h.^5.*omega.^4.*w.^5.*6.0-pi.^3.*dw_dx.^2.*g.*h.^7.*omega.^4.*w.^3.*1.6e1-pi.^3.*dw_dx.^4.*g.^2.*h.^8.*omega.^2.*w.*8.0-pi.^(5.0./2.0).*dw_dx.^3.*g.^(3.0./2.0).*h.^(1.3e1./2.0).*omega.^2.*w.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0i+pi.^3.*dh_dx.*dw_dx.^3.*g.^2.*h.^7.*omega.^2.*w.^2.*1.0e1-pi.^2.*Q.*cd.*dw_dx.^4.*g.^2.*h.^6.*omega.*4.0i-pi.^3.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^6.*omega.^2.*w.^3.*4.0-pi.*Q.^2.*cd.^2.*dh_dx.^2.*g.*h.*omega.^2.*w.^3.*2.4e1+pi.^2.*Q.*cd.*dh_dx.^2.*g.*h.^3.*omega.^3.*w.^4.*1.04e2i+pi.*Q.^2.*cd.^2.*dw_dx.^2.*g.*h.^3.*omega.^2.*w.*4.0e1-pi.^2.*Q.*cd.*dw_dx.^2.*g.*h.^5.*omega.^3.*w.^2.*3.2e1i-pi.^(5.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(1.1e1./2.0).*omega.^2.*w.^2.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0i+pi.^(3.0./2.0).*Q.*cd.*dw_dx.^3.*g.^(3.0./2.0).*h.^(9.0./2.0).*omega.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0-pi.^2.*Q.*cd.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^4.*omega.*w.^2.*4.8e1i+pi.^2.*Q.*cd.*dh_dx.*dw_dx.*g.*h.^4.*omega.^3.*w.^3.*5.6e1i+pi.^2.*Q.*cd.*dh_dx.*dw_dx.^3.*g.^2.*h.^5.*omega.*w.*3.6e1i+pi.*Q.^2.*cd.^2.*dh_dx.*dw_dx.*g.*h.^2.*omega.^2.*w.^2.*4.8e1).*(-1.0./4.0)+pi.*g.*h.^2.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(pi.*dw_dx.^3.*g.*h.^4.*1i-pi.*dh_dx.*h.^2.*omega.^2.*w.^3.*1i-pi.*dw_dx.*h.^3.*omega.^2.*w.^2.*2.0i+Q.*cd.*dh_dx.*omega.*w.^2.*6.0+Q.*cd.*dw_dx.*h.*omega.*w.*6.0).*(Q.^3.*cd.^3.*omega.^3.*w.^2.*-2.56e2i+pi.^3.*h.^6.*omega.^6.*w.^5.*3.2e1+pi.^2.*Q.*cd.*h.^4.*omega.^5.*w.^4.*1.92e2i+pi.^(3.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(7.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*3.0i+pi.^(5.0./2.0).*dh_dx.*dw_dx.^4.*g.^(5.0./2.0).*h.^(1.3e1./2.0).*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*3.0i-pi.*Q.^2.*cd.^2.*h.^2.*omega.^4.*w.^3.*3.84e2+pi.^3.*dh_dx.^2.*g.*h.^5.*omega.^4.*w.^5.*6.0-pi.^3.*dw_dx.^2.*g.*h.^7.*omega.^4.*w.^3.*1.6e1-pi.^3.*dw_dx.^4.*g.^2.*h.^8.*omega.^2.*w.*8.0-pi.^(5.0./2.0).*dw_dx.^3.*g.^(3.0./2.0).*h.^(1.3e1./2.0).*omega.^2.*w.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0i+pi.^3.*dh_dx.*dw_dx.^3.*g.^2.*h.^7.*omega.^2.*w.^2.*1.0e1-pi.^2.*Q.*cd.*dw_dx.^4.*g.^2.*h.^6.*omega.*4.0i-pi.^3.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^6.*omega.^2.*w.^3.*4.0-pi.*Q.^2.*cd.^2.*dh_dx.^2.*g.*h.*omega.^2.*w.^3.*2.4e1+pi.^2.*Q.*cd.*dh_dx.^2.*g.*h.^3.*omega.^3.*w.^4.*1.04e2i+pi.*Q.^2.*cd.^2.*dw_dx.^2.*g.*h.^3.*omega.^2.*w.*4.0e1-pi.^2.*Q.*cd.*dw_dx.^2.*g.*h.^5.*omega.^3.*w.^2.*3.2e1i-pi.^(5.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(1.1e1./2.0).*omega.^2.*w.^2.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0i+pi.^(3.0./2.0).*Q.*cd.*dw_dx.^3.*g.^(3.0./2.0).*h.^(9.0./2.0).*omega.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0-pi.^2.*Q.*cd.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^4.*omega.*w.^2.*4.8e1i+pi.^2.*Q.*cd.*dh_dx.*dw_dx.*g.*h.^4.*omega.^3.*w.^3.*5.6e1i+pi.^2.*Q.*cd.*dh_dx.*dw_dx.^3.*g.^2.*h.^5.*omega.*w.*3.6e1i+pi.*Q.^2.*cd.^2.*dh_dx.*dw_dx.*g.*h.^2.*omega.^2.*w.^2.*4.8e1).*(1.0./2.0))./(pi.*g.*h.^3.*(pi.*h.^2.*omega.^2.*w.^2.*2.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*4.0i).^2.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^3-pi.^2.*g.^2.*h.^4.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^2.*(pi.*dw_dx.^3.*g.*h.^4.*1i-pi.*dh_dx.*h.^2.*omega.^2.*w.^3.*1i-pi.*dw_dx.*h.^3.*omega.^2.*w.^2.*2.0i+Q.*cd.*dh_dx.*omega.*w.^2.*6.0+Q.*cd.*dw_dx.*h.*omega.*w.*6.0).^2)];

end % River_Tide/wave_number_analytic

