% Fri  1 Dec 14:21:20 CET 2017
%% analytic expression of the wave number
%%
%% valid for both tidally, river dominated and low friction conditions
%% and converging channels
%%
%% k      : complex wave number in a reach with constant width and bed slope
%% im(k)  : damping modulus (rate of amplitude change)
%% re(k)  : actual wave number (rate of phase change)
%%
%% c.f. derive_wave_number
%%
% TODO dh_dx was x
% TODO remove az1 from arguments
function [k10, kq, kz, c] = wave_number_analytic(obj,Q0,w,h,cd,omega,az1,Qt,dh_dx,dw_dx)
	g = obj.g;

	% specific discharge
	q0 = Q0./w;
	Q2 = 0;

	% estimate qt
	% TODO this is a poor estimate, ignoring frictional damping and bed slope
	if (nargin()<7 || isempty(Qt))
		qt = az1.*sqrt(g.*h);
		Qt = qt*w;
	else
		qt = Qt/w;
	end
	if (~obj.issym)
		Pi = pi;	
	else
		Pi = sym('pi');
	end
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
	c = obj.odefun1(Q0, Qhr, Qt, Q2, h, dh_dx, 0, w, dw_dx, cd, g, p, omega,flag);

	% roots of the characteristic polynomial
	r = roots2(c(:,1:3));
	
	% wave number
	k10 = 1i*r;

	Q = 1/2*(p1.*Qt + 2*p2.*Q0);
	% simplified for q
	kq = [dw_dx.*((1.0./2.0)./w-sqrt(pi).*sqrt(g).*h.^(3.0./2.0).*w.*1.0./sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*(1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(3.0./2.0).*1.0./w.^2.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*5.0e-1i-(1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(3.0./2.0).*1.0./sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*(Q.*cd.*omega.*8.0i+pi.*h.^2.*omega.^2.*w.*8.0).*2.5e-1i)./w).*1i)-(1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(3.0./2.0).*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*5.0e-1i)./w-sqrt(pi).*dh_dx.*sqrt(g).*h.^(3.0./2.0).*w.*((1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(5.0./2.0).*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*7.5e-1i)./w-sqrt(pi).*1.0./sqrt(g).*1.0./sqrt(h).*omega.^2.*w.*1.0./sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*2.0i).*1.0./sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*1i, ...
	dw_dx.*((1.0./2.0)./w-sqrt(pi).*sqrt(g).*h.^(3.0./2.0).*w.*1.0./sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*(1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(3.0./2.0).*1.0./w.^2.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*5.0e-1i-(1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(3.0./2.0).*1.0./sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*(Q.*cd.*omega.*8.0i+pi.*h.^2.*omega.^2.*w.*8.0).*2.5e-1i)./w).*1i)+(1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(3.0./2.0).*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*5.0e-1i)./w-sqrt(pi).*dh_dx.*sqrt(g).*h.^(3.0./2.0).*w.*((1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(5.0./2.0).*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*7.5e-1i)./w-sqrt(pi).*1.0./sqrt(g).*1.0./sqrt(h).*omega.^2.*w.*1.0./sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*2.0i).*1.0./sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*1i];

if (0)
	% simplified
	kz = [(dh_dx.*-7.5e-1i)./h-(dw_dx.*1i)./w+(1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(3.0./2.0).*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*(1.0./2.0))./w-(Q.*cd.*dw_dx.*omega.*2.0)./(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i)+(pi.*dh_dx.*h.*omega.^2.*w.^2.*2.0i)./(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i)+(pi.*dw_dx.*h.^2.*omega.^2.*w.*2.0i)./(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i), ...
	(dh_dx.*-7.5e-1i)./h-(dw_dx.*1i)./w-(1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(3.0./2.0).*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i).*(1.0./2.0))./w-(Q.*cd.*dw_dx.*omega.*2.0)./(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i)+(pi.*dh_dx.*h.*omega.^2.*w.^2.*2.0i)./(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i)+(pi.*dw_dx.*h.^2.*omega.^2.*w.*2.0i)./(pi.*h.^2.*omega.^2.*w.^2.*4.0+Q.*cd.*omega.*w.*8.0i)];

%	kz = [1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(3.0./2.0).*1.0./omega.^(5.0./2.0).*1.0./w.^(7.0./2.0).*1.0./(Q.*cd.*2.0i+pi.*h.^2.*omega.*w).^(7.0./2.0).*(Q.^4.*cd.^4.*omega.^3.*w.^3.*-6.4e1i-pi.^4.*h.^8.*omega.^7.*w.^7.*4.0i+pi.^3.*Q.*cd.*h.^6.*omega.^6.*w.^6.*3.2e1-pi.*Q.^3.*cd.^3.*h.^2.*omega.^4.*w.^4.*1.28e2+pi.^2.*Q.^2.*cd.^2.*h.^4.*omega.^5.*w.^5.*9.6e1i-pi.^4.*dh_dx.*dw_dx.*g.*h.^8.*omega.^5.*w.^6+pi.^(3.0./2.0).*dh_dx.*sqrt(g).*h.^(5.0./2.0).*omega.^(7.0./2.0).*w.^(9.0./2.0).*(Q.*cd.*2.0i+pi.*h.^2.*omega.*w).^(5.0./2.0).*(1.0-2.0i)+pi.^(3.0./2.0).*dw_dx.*sqrt(g).*h.^(7.0./2.0).*omega.^(7.0./2.0).*w.^(7.0./2.0).*(Q.*cd.*2.0i+pi.*h.^2.*omega.*w).^(5.0./2.0).*(2.0-4.0i)+pi.^2.*Q.^2.*cd.^2.*dh_dx.*dw_dx.*g.*h.^4.*omega.^3.*w.^4.*2.8e1-pi.^3.*Q.*cd.*dh_dx.*dw_dx.*g.*h.^6.*omega.^4.*w.^5.*4.0i+sqrt(pi).*Q.*cd.*dh_dx.*sqrt(g).*sqrt(h).*omega.^(5.0./2.0).*w.^(7.0./2.0).*(Q.*cd.*2.0i+pi.*h.^2.*omega.*w).^(5.0./2.0).*(1.2e1+6.0i)+sqrt(pi).*Q.*cd.*dw_dx.*sqrt(g).*h.^(3.0./2.0).*omega.^(5.0./2.0).*w.^(5.0./2.0).*(Q.*cd.*2.0i+pi.*h.^2.*omega.*w).^(5.0./2.0).*(1.2e1+6.0i)+pi.*Q.^3.*cd.^3.*dh_dx.*dw_dx.*g.*h.^2.*omega.^2.*w.^3.*4.8e1i).*(1.0./4.0), ...
%	1.0./sqrt(pi).*1.0./sqrt(g).*1.0./h.^(3.0./2.0).*1.0./omega.^(5.0./2.0).*1.0./w.^(7.0./2.0).*1.0./(Q.*cd.*2.0i+pi.*h.^2.*omega.*w).^(7.0./2.0).*(Q.^4.*cd.^4.*omega.^3.*w.^3.*6.4e1i+pi.^4.*h.^8.*omega.^7.*w.^7.*4.0i-pi.^3.*Q.*cd.*h.^6.*omega.^6.*w.^6.*3.2e1+pi.*Q.^3.*cd.^3.*h.^2.*omega.^4.*w.^4.*1.28e2-pi.^2.*Q.^2.*cd.^2.*h.^4.*omega.^5.*w.^5.*9.6e1i+pi.^4.*dh_dx.*dw_dx.*g.*h.^8.*omega.^5.*w.^6+pi.^(3.0./2.0).*dh_dx.*sqrt(g).*h.^(5.0./2.0).*omega.^(7.0./2.0).*w.^(9.0./2.0).*(Q.*cd.*2.0i+pi.*h.^2.*omega.*w).^(5.0./2.0).*(1.0-2.0i)+pi.^(3.0./2.0).*dw_dx.*sqrt(g).*h.^(7.0./2.0).*omega.^(7.0./2.0).*w.^(7.0./2.0).*(Q.*cd.*2.0i+pi.*h.^2.*omega.*w).^(5.0./2.0).*(2.0-4.0i)-pi.^2.*Q.^2.*cd.^2.*dh_dx.*dw_dx.*g.*h.^4.*omega.^3.*w.^4.*2.8e1+pi.^3.*Q.*cd.*dh_dx.*dw_dx.*g.*h.^6.*omega.^4.*w.^5.*4.0i+sqrt(pi).*Q.*cd.*dh_dx.*sqrt(g).*sqrt(h).*omega.^(5.0./2.0).*w.^(7.0./2.0).*(Q.*cd.*2.0i+pi.*h.^2.*omega.*w).^(5.0./2.0).*(1.2e1+6.0i)+sqrt(pi).*Q.*cd.*dw_dx.*sqrt(g).*h.^(3.0./2.0).*omega.^(5.0./2.0).*w.^(5.0./2.0).*(Q.*cd.*2.0i+pi.*h.^2.*omega.*w).^(5.0./2.0).*(1.2e1+6.0i)-pi.*Q.^3.*cd.^3.*dh_dx.*dw_dx.*g.*h.^2.*omega.^2.*w.^3.*4.8e1i).*(1.0./4.0)];
else
	% not expanded
%	kz = [(pi.*g.*h.^2.*w.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(-pi.*dw_dx.^3.*g.*h.^4+pi.*dh_dx.*h.^2.*omega.^2.*w.^3+pi.*dw_dx.*h.^3.*omega.^2.*w.^2.*2.0+Q.*cd.*dh_dx.*omega.*w.^2.*6.0i+Q.*cd.*dw_dx.*h.*omega.*w.*6.0i).*(pi.^3.*dw_dx.^6.*g.^3.*h.^9.*(8.0-8.0i)+Q.^3.*cd.^3.*omega.^3.*w.^3.*5.12e2i-pi.^3.*h.^6.*omega.^6.*w.^6.*6.4e1+sqrt(pi).*dw_dx.*sqrt(g).*h.^(3.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(5.0./2.0).*(-4.0-4.0i)+pi.^(3.0./2.0).*dw_dx.^3.*g.^(3.0./2.0).*h.^(9.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(4.0+4.0i)+pi.^3.*dw_dx.^4.*g.^2.*h.^8.*omega.^2.*w.^2.*(-4.8e1+6.4e1i)-pi.^2.*Q.*cd.*h.^4.*omega.^5.*w.^5.*3.84e2i+pi.*Q.^2.*cd.^2.*h.^2.*omega.^4.*w.^4.*7.68e2+pi.^3.*dh_dx.^2.*g.*h.^5.*omega.^4.*w.^6.*(4.0-1.6e1i)+pi.^3.*dw_dx.^2.*g.*h.^7.*omega.^4.*w.^4.*(9.6e1-6.4e1i)+sqrt(pi).*dh_dx.*sqrt(g).*sqrt(h).*w.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(5.0./2.0).*(-3.0-3.0i)+pi.^(3.0./2.0).*dh_dx.*sqrt(g).*h.^(5.0./2.0).*omega.^2.*w.^3.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(8.0+8.0i)+pi.^(3.0./2.0).*dw_dx.*sqrt(g).*h.^(7.0./2.0).*omega.^2.*w.^2.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(8.0+8.0i)+pi.^3.*dh_dx.*dw_dx.^3.*g.^2.*h.^7.*omega.^2.*w.^3.*(-8.0-1.2e1i)+pi.^(5.0./2.0).*dw_dx.^3.*g.^(3.0./2.0).*h.^(1.3e1./2.0).*omega.^2.*w.^2.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(8.0-8.0i)+pi.^3.*dh_dx.*dw_dx.*g.*h.^6.*omega.^4.*w.^5.*(1.6e1-1.6e1i)-pi.^(3.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(7.0./2.0).*w.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*6.0+pi.^(5.0./2.0).*dh_dx.*dw_dx.^4.*g.^(5.0./2.0).*h.^(1.3e1./2.0).*w.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(-3.0+3.0i)+pi.^3.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^6.*omega.^2.*w.^4.*8.0i+pi.*Q.^2.*cd.^2.*dh_dx.^2.*g.*h.*omega.^2.*w.^4.*(-1.44e2+1.92e2i)+pi.^2.*Q.*cd.*dh_dx.^2.*g.*h.^3.*omega.^3.*w.^5.*(2.56e2+4.8e1i)+pi.^2.*Q.*cd.*dw_dx.^2.*g.*h.^5.*omega.^3.*w.^3.*(3.52e2+4.16e2i)+pi.*Q.^2.*cd.^2.*dw_dx.^2.*g.*h.^3.*omega.^2.*w.^2.*(-4.64e2+3.84e2i)+pi.^(5.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(1.1e1./2.0).*omega.^2.*w.^3.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(8.0-8.0i)+pi.^2.*Q.*cd.*dw_dx.^4.*g.^2.*h.^6.*omega.*w.*(-1.2e2-1.12e2i)+sqrt(pi).*Q.*cd.*dw_dx.*sqrt(g).*h.^(3.0./2.0).*omega.*w.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(-8.0+8.0i)-pi.^2.*Q.*cd.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^4.*omega.*w.^3.*9.6e1+pi.^2.*Q.*cd.*dh_dx.*dw_dx.*g.*h.^4.*omega.^3.*w.^4.*(2.56e2+1.44e2i)+pi.^(3.0./2.0).*Q.*cd.*dw_dx.^3.*g.^(3.0./2.0).*h.^(9.0./2.0).*omega.*w.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(8.0+8.0i)+pi.*Q.^2.*cd.^2.*dh_dx.*dw_dx.*g.*h.^2.*omega.^2.*w.^3.*(-2.88e2+1.92e2i)+pi.^2.*Q.*cd.*dh_dx.*dw_dx.^3.*g.^2.*h.^5.*omega.*w.^2.*(2.4e1-4.8e1i)).*(-1.0./4.0)+sqrt(pi).*sqrt(g).*h.^(3.0./2.0).*w.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*-4.0i+pi.*dw_dx.^2.*g.*h.^3.*2.0i+Q.*cd.*omega.*w.*8.0).*(pi.^3.*dw_dx.^6.*g.^3.*h.^9.*(8.0-8.0i)+Q.^3.*cd.^3.*omega.^3.*w.^3.*5.12e2i-pi.^3.*h.^6.*omega.^6.*w.^6.*6.4e1+sqrt(pi).*dw_dx.*sqrt(g).*h.^(3.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(5.0./2.0).*(-4.0-4.0i)+pi.^(3.0./2.0).*dw_dx.^3.*g.^(3.0./2.0).*h.^(9.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(4.0+4.0i)+pi.^3.*dw_dx.^4.*g.^2.*h.^8.*omega.^2.*w.^2.*(-4.8e1+6.4e1i)-pi.^2.*Q.*cd.*h.^4.*omega.^5.*w.^5.*3.84e2i+pi.*Q.^2.*cd.^2.*h.^2.*omega.^4.*w.^4.*7.68e2+pi.^3.*dh_dx.^2.*g.*h.^5.*omega.^4.*w.^6.*(4.0-1.6e1i)+pi.^3.*dw_dx.^2.*g.*h.^7.*omega.^4.*w.^4.*(9.6e1-6.4e1i)+sqrt(pi).*dh_dx.*sqrt(g).*sqrt(h).*w.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(5.0./2.0).*(-3.0-3.0i)+pi.^(3.0./2.0).*dh_dx.*sqrt(g).*h.^(5.0./2.0).*omega.^2.*w.^3.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(8.0+8.0i)+pi.^(3.0./2.0).*dw_dx.*sqrt(g).*h.^(7.0./2.0).*omega.^2.*w.^2.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(8.0+8.0i)+pi.^3.*dh_dx.*dw_dx.^3.*g.^2.*h.^7.*omega.^2.*w.^3.*(-8.0-1.2e1i)+pi.^(5.0./2.0).*dw_dx.^3.*g.^(3.0./2.0).*h.^(1.3e1./2.0).*omega.^2.*w.^2.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(8.0-8.0i)+pi.^3.*dh_dx.*dw_dx.*g.*h.^6.*omega.^4.*w.^5.*(1.6e1-1.6e1i)-pi.^(3.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(7.0./2.0).*w.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*6.0+pi.^(5.0./2.0).*dh_dx.*dw_dx.^4.*g.^(5.0./2.0).*h.^(1.3e1./2.0).*w.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(-3.0+3.0i)+pi.^3.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^6.*omega.^2.*w.^4.*8.0i+pi.*Q.^2.*cd.^2.*dh_dx.^2.*g.*h.*omega.^2.*w.^4.*(-1.44e2+1.92e2i)+pi.^2.*Q.*cd.*dh_dx.^2.*g.*h.^3.*omega.^3.*w.^5.*(2.56e2+4.8e1i)+pi.^2.*Q.*cd.*dw_dx.^2.*g.*h.^5.*omega.^3.*w.^3.*(3.52e2+4.16e2i)+pi.*Q.^2.*cd.^2.*dw_dx.^2.*g.*h.^3.*omega.^2.*w.^2.*(-4.64e2+3.84e2i)+pi.^(5.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(1.1e1./2.0).*omega.^2.*w.^3.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(8.0-8.0i)+pi.^2.*Q.*cd.*dw_dx.^4.*g.^2.*h.^6.*omega.*w.*(-1.2e2-1.12e2i)+sqrt(pi).*Q.*cd.*dw_dx.*sqrt(g).*h.^(3.0./2.0).*omega.*w.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(-8.0+8.0i)-pi.^2.*Q.*cd.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^4.*omega.*w.^3.*9.6e1+pi.^2.*Q.*cd.*dh_dx.*dw_dx.*g.*h.^4.*omega.^3.*w.^4.*(2.56e2+1.44e2i)+pi.^(3.0./2.0).*Q.*cd.*dw_dx.^3.*g.^(3.0./2.0).*h.^(9.0./2.0).*omega.*w.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(8.0+8.0i)+pi.*Q.^2.*cd.^2.*dh_dx.*dw_dx.*g.*h.^2.*omega.^2.*w.^3.*(-2.88e2+1.92e2i)+pi.^2.*Q.*cd.*dh_dx.*dw_dx.^3.*g.^2.*h.^5.*omega.*w.^2.*(2.4e1-4.8e1i)).*(1.0./8.0))./(pi.*g.*h.^3.*w.^2.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^3.*(pi.*h.^2.*omega.^2.*w.^2.*-2.0i+pi.*dw_dx.^2.*g.*h.^3.*1i+Q.*cd.*omega.*w.*4.0).^2-pi.^2.*g.^2.*h.^4.*w.^2.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^2.*(-pi.*dw_dx.^3.*g.*h.^4+pi.*dh_dx.*h.^2.*omega.^2.*w.^3+pi.*dw_dx.*h.^3.*omega.^2.*w.^2.*2.0+Q.*cd.*dh_dx.*omega.*w.^2.*6.0i+Q.*cd.*dw_dx.*h.*omega.*w.*6.0i).^2), ...
%(pi.*g.*h.^2.*w.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(-pi.*dw_dx.^3.*g.*h.^4+pi.*dh_dx.*h.^2.*omega.^2.*w.^3+pi.*dw_dx.*h.^3.*omega.^2.*w.^2.*2.0+Q.*cd.*dh_dx.*omega.*w.^2.*6.0i+Q.*cd.*dw_dx.*h.*omega.*w.*6.0i).*(pi.^3.*dw_dx.^6.*g.^3.*h.^9.*(8.0-8.0i)+Q.^3.*cd.^3.*omega.^3.*w.^3.*5.12e2i-pi.^3.*h.^6.*omega.^6.*w.^6.*6.4e1+sqrt(pi).*dw_dx.*sqrt(g).*h.^(3.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(5.0./2.0).*(4.0+4.0i)+pi.^(3.0./2.0).*dw_dx.^3.*g.^(3.0./2.0).*h.^(9.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(-4.0-4.0i)+pi.^3.*dw_dx.^4.*g.^2.*h.^8.*omega.^2.*w.^2.*(-4.8e1+6.4e1i)-pi.^2.*Q.*cd.*h.^4.*omega.^5.*w.^5.*3.84e2i+pi.*Q.^2.*cd.^2.*h.^2.*omega.^4.*w.^4.*7.68e2+pi.^3.*dh_dx.^2.*g.*h.^5.*omega.^4.*w.^6.*(4.0-1.6e1i)+pi.^3.*dw_dx.^2.*g.*h.^7.*omega.^4.*w.^4.*(9.6e1-6.4e1i)+sqrt(pi).*dh_dx.*sqrt(g).*sqrt(h).*w.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(5.0./2.0).*(3.0+3.0i)+pi.^(3.0./2.0).*dh_dx.*sqrt(g).*h.^(5.0./2.0).*omega.^2.*w.^3.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(-8.0-8.0i)+pi.^(3.0./2.0).*dw_dx.*sqrt(g).*h.^(7.0./2.0).*omega.^2.*w.^2.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(-8.0-8.0i)+pi.^3.*dh_dx.*dw_dx.^3.*g.^2.*h.^7.*omega.^2.*w.^3.*(-8.0-1.2e1i)+pi.^(5.0./2.0).*dw_dx.^3.*g.^(3.0./2.0).*h.^(1.3e1./2.0).*omega.^2.*w.^2.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(-8.0+8.0i)+pi.^3.*dh_dx.*dw_dx.*g.*h.^6.*omega.^4.*w.^5.*(1.6e1-1.6e1i)+pi.^(3.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(7.0./2.0).*w.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*6.0+pi.^(5.0./2.0).*dh_dx.*dw_dx.^4.*g.^(5.0./2.0).*h.^(1.3e1./2.0).*w.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(3.0-3.0i)+pi.^3.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^6.*omega.^2.*w.^4.*8.0i+pi.*Q.^2.*cd.^2.*dh_dx.^2.*g.*h.*omega.^2.*w.^4.*(-1.44e2+1.92e2i)+pi.^2.*Q.*cd.*dh_dx.^2.*g.*h.^3.*omega.^3.*w.^5.*(2.56e2+4.8e1i)+pi.^2.*Q.*cd.*dw_dx.^2.*g.*h.^5.*omega.^3.*w.^3.*(3.52e2+4.16e2i)+pi.*Q.^2.*cd.^2.*dw_dx.^2.*g.*h.^3.*omega.^2.*w.^2.*(-4.64e2+3.84e2i)+pi.^(5.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(1.1e1./2.0).*omega.^2.*w.^3.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(-8.0+8.0i)+pi.^2.*Q.*cd.*dw_dx.^4.*g.^2.*h.^6.*omega.*w.*(-1.2e2-1.12e2i)+sqrt(pi).*Q.*cd.*dw_dx.*sqrt(g).*h.^(3.0./2.0).*omega.*w.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(8.0-8.0i)-pi.^2.*Q.*cd.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^4.*omega.*w.^3.*9.6e1+pi.^2.*Q.*cd.*dh_dx.*dw_dx.*g.*h.^4.*omega.^3.*w.^4.*(2.56e2+1.44e2i)+pi.^(3.0./2.0).*Q.*cd.*dw_dx.^3.*g.^(3.0./2.0).*h.^(9.0./2.0).*omega.*w.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(-8.0-8.0i)+pi.*Q.^2.*cd.^2.*dh_dx.*dw_dx.*g.*h.^2.*omega.^2.*w.^3.*(-2.88e2+1.92e2i)+pi.^2.*Q.*cd.*dh_dx.*dw_dx.^3.*g.^2.*h.^5.*omega.*w.^2.*(2.4e1-4.8e1i)).*(-1.0./4.0)-sqrt(pi).*sqrt(g).*h.^(3.0./2.0).*w.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*-4.0i+pi.*dw_dx.^2.*g.*h.^3.*2.0i+Q.*cd.*omega.*w.*8.0).*(pi.^3.*dw_dx.^6.*g.^3.*h.^9.*(8.0-8.0i)+Q.^3.*cd.^3.*omega.^3.*w.^3.*5.12e2i-pi.^3.*h.^6.*omega.^6.*w.^6.*6.4e1+sqrt(pi).*dw_dx.*sqrt(g).*h.^(3.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(5.0./2.0).*(4.0+4.0i)+pi.^(3.0./2.0).*dw_dx.^3.*g.^(3.0./2.0).*h.^(9.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(-4.0-4.0i)+pi.^3.*dw_dx.^4.*g.^2.*h.^8.*omega.^2.*w.^2.*(-4.8e1+6.4e1i)-pi.^2.*Q.*cd.*h.^4.*omega.^5.*w.^5.*3.84e2i+pi.*Q.^2.*cd.^2.*h.^2.*omega.^4.*w.^4.*7.68e2+pi.^3.*dh_dx.^2.*g.*h.^5.*omega.^4.*w.^6.*(4.0-1.6e1i)+pi.^3.*dw_dx.^2.*g.*h.^7.*omega.^4.*w.^4.*(9.6e1-6.4e1i)+sqrt(pi).*dh_dx.*sqrt(g).*sqrt(h).*w.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(5.0./2.0).*(3.0+3.0i)+pi.^(3.0./2.0).*dh_dx.*sqrt(g).*h.^(5.0./2.0).*omega.^2.*w.^3.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(-8.0-8.0i)+pi.^(3.0./2.0).*dw_dx.*sqrt(g).*h.^(7.0./2.0).*omega.^2.*w.^2.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(-8.0-8.0i)+pi.^3.*dh_dx.*dw_dx.^3.*g.^2.*h.^7.*omega.^2.*w.^3.*(-8.0-1.2e1i)+pi.^(5.0./2.0).*dw_dx.^3.*g.^(3.0./2.0).*h.^(1.3e1./2.0).*omega.^2.*w.^2.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(-8.0+8.0i)+pi.^3.*dh_dx.*dw_dx.*g.*h.^6.*omega.^4.*w.^5.*(1.6e1-1.6e1i)+pi.^(3.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(7.0./2.0).*w.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*6.0+pi.^(5.0./2.0).*dh_dx.*dw_dx.^4.*g.^(5.0./2.0).*h.^(1.3e1./2.0).*w.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(3.0-3.0i)+pi.^3.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^6.*omega.^2.*w.^4.*8.0i+pi.*Q.^2.*cd.^2.*dh_dx.^2.*g.*h.*omega.^2.*w.^4.*(-1.44e2+1.92e2i)+pi.^2.*Q.*cd.*dh_dx.^2.*g.*h.^3.*omega.^3.*w.^5.*(2.56e2+4.8e1i)+pi.^2.*Q.*cd.*dw_dx.^2.*g.*h.^5.*omega.^3.*w.^3.*(3.52e2+4.16e2i)+pi.*Q.^2.*cd.^2.*dw_dx.^2.*g.*h.^3.*omega.^2.*w.^2.*(-4.64e2+3.84e2i)+pi.^(5.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(1.1e1./2.0).*omega.^2.*w.^3.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(-8.0+8.0i)+pi.^2.*Q.*cd.*dw_dx.^4.*g.^2.*h.^6.*omega.*w.*(-1.2e2-1.12e2i)+sqrt(pi).*Q.*cd.*dw_dx.*sqrt(g).*h.^(3.0./2.0).*omega.*w.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(8.0-8.0i)-pi.^2.*Q.*cd.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^4.*omega.*w.^3.*9.6e1+pi.^2.*Q.*cd.*dh_dx.*dw_dx.*g.*h.^4.*omega.^3.*w.^4.*(2.56e2+1.44e2i)+pi.^(3.0./2.0).*Q.*cd.*dw_dx.^3.*g.^(3.0./2.0).*h.^(9.0./2.0).*omega.*w.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(-8.0-8.0i)+pi.*Q.^2.*cd.^2.*dh_dx.*dw_dx.*g.*h.^2.*omega.^2.*w.^3.*(-2.88e2+1.92e2i)+pi.^2.*Q.*cd.*dh_dx.*dw_dx.^3.*g.^2.*h.^5.*omega.*w.^2.*(2.4e1-4.8e1i)).*(1.0./8.0))./(pi.*g.*h.^3.*w.^2.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^3.*(pi.*h.^2.*omega.^2.*w.^2.*-2.0i+pi.*dw_dx.^2.*g.*h.^3.*1i+Q.*cd.*omega.*w.*4.0).^2-pi.^2.*g.^2.*h.^4.*w.^2.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^2.*(-pi.*dw_dx.^3.*g.*h.^4+pi.*dh_dx.*h.^2.*omega.^2.*w.^3+pi.*dw_dx.*h.^3.*omega.^2.*w.^2.*2.0+Q.*cd.*dh_dx.*omega.*w.^2.*6.0i+Q.*cd.*dw_dx.*h.*omega.*w.*6.0i).^2)];

kz = [(sqrt(pi).*sqrt(g).*h.^(3.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3.*2.0+Q.*cd.*omega.*w.*8.0i).*(Q.^3.*cd.^3.*omega.^3.*w.^2.*2.56e2i-pi.^3.*h.^6.*omega.^6.*w.^5.*3.2e1-pi.^2.*Q.*cd.*h.^4.*omega.^5.*w.^4.*1.92e2i+pi.^(3.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(7.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*3.0i+pi.^(5.0./2.0).*dh_dx.*dw_dx.^4.*g.^(5.0./2.0).*h.^(1.3e1./2.0).*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*3.0i+pi.*Q.^2.*cd.^2.*h.^2.*omega.^4.*w.^3.*3.84e2-pi.^3.*dh_dx.^2.*g.*h.^5.*omega.^4.*w.^5.*6.0+pi.^3.*dw_dx.^2.*g.*h.^7.*omega.^4.*w.^3.*1.6e1+pi.^3.*dw_dx.^4.*g.^2.*h.^8.*omega.^2.*w.*8.0-pi.^(5.0./2.0).*dw_dx.^3.*g.^(3.0./2.0).*h.^(1.3e1./2.0).*omega.^2.*w.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0i-pi.^3.*dh_dx.*dw_dx.^3.*g.^2.*h.^7.*omega.^2.*w.^2.*1.0e1+pi.^2.*Q.*cd.*dw_dx.^4.*g.^2.*h.^6.*omega.*4.0i+pi.^3.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^6.*omega.^2.*w.^3.*4.0+pi.*Q.^2.*cd.^2.*dh_dx.^2.*g.*h.*omega.^2.*w.^3.*2.4e1-pi.^2.*Q.*cd.*dh_dx.^2.*g.*h.^3.*omega.^3.*w.^4.*1.04e2i-pi.*Q.^2.*cd.^2.*dw_dx.^2.*g.*h.^3.*omega.^2.*w.*4.0e1+pi.^2.*Q.*cd.*dw_dx.^2.*g.*h.^5.*omega.^3.*w.^2.*3.2e1i-pi.^(5.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(1.1e1./2.0).*omega.^2.*w.^2.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0i+pi.^(3.0./2.0).*Q.*cd.*dw_dx.^3.*g.^(3.0./2.0).*h.^(9.0./2.0).*omega.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0+pi.^2.*Q.*cd.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^4.*omega.*w.^2.*4.8e1i-pi.^2.*Q.*cd.*dh_dx.*dw_dx.*g.*h.^4.*omega.^3.*w.^3.*5.6e1i-pi.^2.*Q.*cd.*dh_dx.*dw_dx.^3.*g.^2.*h.^5.*omega.*w.*3.6e1i-pi.*Q.^2.*cd.^2.*dh_dx.*dw_dx.*g.*h.^2.*omega.^2.*w.^2.*4.8e1).*(-1.0./4.0)-pi.*g.*h.^2.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(pi.*dw_dx.^3.*g.*h.^4.*1i-pi.*dh_dx.*h.^2.*omega.^2.*w.^3.*1i-pi.*dw_dx.*h.^3.*omega.^2.*w.^2.*2.0i+Q.*cd.*dh_dx.*omega.*w.^2.*6.0+Q.*cd.*dw_dx.*h.*omega.*w.*6.0).*(Q.^3.*cd.^3.*omega.^3.*w.^2.*2.56e2i-pi.^3.*h.^6.*omega.^6.*w.^5.*3.2e1-pi.^2.*Q.*cd.*h.^4.*omega.^5.*w.^4.*1.92e2i+pi.^(3.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(7.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*3.0i+pi.^(5.0./2.0).*dh_dx.*dw_dx.^4.*g.^(5.0./2.0).*h.^(1.3e1./2.0).*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*3.0i+pi.*Q.^2.*cd.^2.*h.^2.*omega.^4.*w.^3.*3.84e2-pi.^3.*dh_dx.^2.*g.*h.^5.*omega.^4.*w.^5.*6.0+pi.^3.*dw_dx.^2.*g.*h.^7.*omega.^4.*w.^3.*1.6e1+pi.^3.*dw_dx.^4.*g.^2.*h.^8.*omega.^2.*w.*8.0-pi.^(5.0./2.0).*dw_dx.^3.*g.^(3.0./2.0).*h.^(1.3e1./2.0).*omega.^2.*w.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0i-pi.^3.*dh_dx.*dw_dx.^3.*g.^2.*h.^7.*omega.^2.*w.^2.*1.0e1+pi.^2.*Q.*cd.*dw_dx.^4.*g.^2.*h.^6.*omega.*4.0i+pi.^3.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^6.*omega.^2.*w.^3.*4.0+pi.*Q.^2.*cd.^2.*dh_dx.^2.*g.*h.*omega.^2.*w.^3.*2.4e1-pi.^2.*Q.*cd.*dh_dx.^2.*g.*h.^3.*omega.^3.*w.^4.*1.04e2i-pi.*Q.^2.*cd.^2.*dw_dx.^2.*g.*h.^3.*omega.^2.*w.*4.0e1+pi.^2.*Q.*cd.*dw_dx.^2.*g.*h.^5.*omega.^3.*w.^2.*3.2e1i-pi.^(5.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(1.1e1./2.0).*omega.^2.*w.^2.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0i+pi.^(3.0./2.0).*Q.*cd.*dw_dx.^3.*g.^(3.0./2.0).*h.^(9.0./2.0).*omega.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0+pi.^2.*Q.*cd.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^4.*omega.*w.^2.*4.8e1i-pi.^2.*Q.*cd.*dh_dx.*dw_dx.*g.*h.^4.*omega.^3.*w.^3.*5.6e1i-pi.^2.*Q.*cd.*dh_dx.*dw_dx.^3.*g.^2.*h.^5.*omega.*w.*3.6e1i-pi.*Q.^2.*cd.^2.*dh_dx.*dw_dx.*g.*h.^2.*omega.^2.*w.^2.*4.8e1).*(1.0./2.0))./(pi.*g.*h.^3.*(pi.*h.^2.*omega.^2.*w.^2.*2.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*4.0i).^2.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^3-pi.^2.*g.^2.*h.^4.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^2.*(pi.*dw_dx.^3.*g.*h.^4.*1i-pi.*dh_dx.*h.^2.*omega.^2.*w.^3.*1i-pi.*dw_dx.*h.^3.*omega.^2.*w.^2.*2.0i+Q.*cd.*dh_dx.*omega.*w.^2.*6.0+Q.*cd.*dw_dx.*h.*omega.*w.*6.0).^2),...
(sqrt(pi).*sqrt(g).*h.^(3.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3.*2.0+Q.*cd.*omega.*w.*8.0i).*(Q.^3.*cd.^3.*omega.^3.*w.^2.*-2.56e2i+pi.^3.*h.^6.*omega.^6.*w.^5.*3.2e1+pi.^2.*Q.*cd.*h.^4.*omega.^5.*w.^4.*1.92e2i+pi.^(3.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(7.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*3.0i+pi.^(5.0./2.0).*dh_dx.*dw_dx.^4.*g.^(5.0./2.0).*h.^(1.3e1./2.0).*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*3.0i-pi.*Q.^2.*cd.^2.*h.^2.*omega.^4.*w.^3.*3.84e2+pi.^3.*dh_dx.^2.*g.*h.^5.*omega.^4.*w.^5.*6.0-pi.^3.*dw_dx.^2.*g.*h.^7.*omega.^4.*w.^3.*1.6e1-pi.^3.*dw_dx.^4.*g.^2.*h.^8.*omega.^2.*w.*8.0-pi.^(5.0./2.0).*dw_dx.^3.*g.^(3.0./2.0).*h.^(1.3e1./2.0).*omega.^2.*w.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0i+pi.^3.*dh_dx.*dw_dx.^3.*g.^2.*h.^7.*omega.^2.*w.^2.*1.0e1-pi.^2.*Q.*cd.*dw_dx.^4.*g.^2.*h.^6.*omega.*4.0i-pi.^3.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^6.*omega.^2.*w.^3.*4.0-pi.*Q.^2.*cd.^2.*dh_dx.^2.*g.*h.*omega.^2.*w.^3.*2.4e1+pi.^2.*Q.*cd.*dh_dx.^2.*g.*h.^3.*omega.^3.*w.^4.*1.04e2i+pi.*Q.^2.*cd.^2.*dw_dx.^2.*g.*h.^3.*omega.^2.*w.*4.0e1-pi.^2.*Q.*cd.*dw_dx.^2.*g.*h.^5.*omega.^3.*w.^2.*3.2e1i-pi.^(5.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(1.1e1./2.0).*omega.^2.*w.^2.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0i+pi.^(3.0./2.0).*Q.*cd.*dw_dx.^3.*g.^(3.0./2.0).*h.^(9.0./2.0).*omega.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0-pi.^2.*Q.*cd.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^4.*omega.*w.^2.*4.8e1i+pi.^2.*Q.*cd.*dh_dx.*dw_dx.*g.*h.^4.*omega.^3.*w.^3.*5.6e1i+pi.^2.*Q.*cd.*dh_dx.*dw_dx.^3.*g.^2.*h.^5.*omega.*w.*3.6e1i+pi.*Q.^2.*cd.^2.*dh_dx.*dw_dx.*g.*h.^2.*omega.^2.*w.^2.*4.8e1).*(-1.0./4.0)+pi.*g.*h.^2.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*(pi.*dw_dx.^3.*g.*h.^4.*1i-pi.*dh_dx.*h.^2.*omega.^2.*w.^3.*1i-pi.*dw_dx.*h.^3.*omega.^2.*w.^2.*2.0i+Q.*cd.*dh_dx.*omega.*w.^2.*6.0+Q.*cd.*dw_dx.*h.*omega.*w.*6.0).*(Q.^3.*cd.^3.*omega.^3.*w.^2.*-2.56e2i+pi.^3.*h.^6.*omega.^6.*w.^5.*3.2e1+pi.^2.*Q.*cd.*h.^4.*omega.^5.*w.^4.*1.92e2i+pi.^(3.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(7.0./2.0).*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^(3.0./2.0).*3.0i+pi.^(5.0./2.0).*dh_dx.*dw_dx.^4.*g.^(5.0./2.0).*h.^(1.3e1./2.0).*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*3.0i-pi.*Q.^2.*cd.^2.*h.^2.*omega.^4.*w.^3.*3.84e2+pi.^3.*dh_dx.^2.*g.*h.^5.*omega.^4.*w.^5.*6.0-pi.^3.*dw_dx.^2.*g.*h.^7.*omega.^4.*w.^3.*1.6e1-pi.^3.*dw_dx.^4.*g.^2.*h.^8.*omega.^2.*w.*8.0-pi.^(5.0./2.0).*dw_dx.^3.*g.^(3.0./2.0).*h.^(1.3e1./2.0).*omega.^2.*w.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0i+pi.^3.*dh_dx.*dw_dx.^3.*g.^2.*h.^7.*omega.^2.*w.^2.*1.0e1-pi.^2.*Q.*cd.*dw_dx.^4.*g.^2.*h.^6.*omega.*4.0i-pi.^3.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^6.*omega.^2.*w.^3.*4.0-pi.*Q.^2.*cd.^2.*dh_dx.^2.*g.*h.*omega.^2.*w.^3.*2.4e1+pi.^2.*Q.*cd.*dh_dx.^2.*g.*h.^3.*omega.^3.*w.^4.*1.04e2i+pi.*Q.^2.*cd.^2.*dw_dx.^2.*g.*h.^3.*omega.^2.*w.*4.0e1-pi.^2.*Q.*cd.*dw_dx.^2.*g.*h.^5.*omega.^3.*w.^2.*3.2e1i-pi.^(5.0./2.0).*dh_dx.*dw_dx.^2.*g.^(3.0./2.0).*h.^(1.1e1./2.0).*omega.^2.*w.^2.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0i+pi.^(3.0./2.0).*Q.*cd.*dw_dx.^3.*g.^(3.0./2.0).*h.^(9.0./2.0).*omega.*sqrt(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).*8.0-pi.^2.*Q.*cd.*dh_dx.^2.*dw_dx.^2.*g.^2.*h.^4.*omega.*w.^2.*4.8e1i+pi.^2.*Q.*cd.*dh_dx.*dw_dx.*g.*h.^4.*omega.^3.*w.^3.*5.6e1i+pi.^2.*Q.*cd.*dh_dx.*dw_dx.^3.*g.^2.*h.^5.*omega.*w.*3.6e1i+pi.*Q.^2.*cd.^2.*dh_dx.*dw_dx.*g.*h.^2.*omega.^2.*w.^2.*4.8e1).*(1.0./2.0))./(pi.*g.*h.^3.*(pi.*h.^2.*omega.^2.*w.^2.*2.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*4.0i).^2.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^3-pi.^2.*g.^2.*h.^4.*(pi.*h.^2.*omega.^2.*w.^2.*4.0-pi.*dw_dx.^2.*g.*h.^3+Q.*cd.*omega.*w.*8.0i).^2.*(pi.*dw_dx.^3.*g.*h.^4.*1i-pi.*dh_dx.*h.^2.*omega.^2.*w.^3.*1i-pi.*dw_dx.*h.^3.*omega.^2.*w.^2.*2.0i+Q.*cd.*dh_dx.*omega.*w.^2.*6.0+Q.*cd.*dw_dx.*h.*omega.*w.*6.0).^2)];

end

	
if (0)
	dk0 = [((Q.*cd.*2+pi.*h.^2.*omega.*w.*1i).*(Q.*cd.*dw_dx.*h+Q.*cd.*dh_dx.*w.*3-pi.*dh_dx.*h.^2.*omega.*w.^2.*5.0e-1i).*2)./(pi.^2.*h.^5.*omega.^2.*w.^3.*4+Q.^2.*cd.^2.*h.*w.*16), ...
	       ((Q.*cd.*2+pi.*h.^2.*omega.*w.*1i).*(Q.*cd.*dw_dx.*h+Q.*cd.*dh_dx.*w.*3-pi.*dh_dx.*h.^2.*omega.*w.^2.*5.0e-1i).*2)./(pi.^2.*h.^5.*omega.^2.*w.^3.*4+Q.^2.*cd.^2.*h.*w.*16)];
	kq = k0 + dk0;
	kz = k0 - dk0;
end

	if (0) 
	
	% TODO, this is only valid for order == 2
	re = -omega.^2./(g.*h);
	im = (   p1./Pi.*omega.*cd./g.*qt./h.^3 ... 
	       + p2./Pi.*2.*omega.*cd./g.*q0./h.^3 );

	if (~obj.issym)
		k0_ = sqrt(re + 1i*im);
	else
		% Note: if the imaginary part later becomes zero,
		%       this is invalid
		[k_re, k_im]  = root_complex(re,im);
		k0_            = k_re+1i*k_im;
	end

	% TODO, this ignores change of width along channel
	% assumptions:
	% dQt/dx neglibible (river dominance)
	% d^2h_dx^2 negligible
	%
	a          = 1/pi*cd./(omega*h.^2.*w).*(p2.*Q0 + 1/2*p1.*Qt);
	dk0_dx_rel = -(3*a.^2 + 1i*a + 1/4)./(4*a.^2 + 1).*(dh_dx./h);
	kz         = k0 - dk0_dx_rel;
	kq         = k0 + dk0_dx_rel;


	if (nargout() > 1)
	if (~obj.issym)
		%k0         = [k0,-k0];
      		%dk0_dx     = derivative1(x,k0);                                         
	        dk0_dx_rel = bsxfun(@times,1./((k0(:,1)-k0(:,2))),dk0_dx);              

		% -(dhdxval*(3*a^2 + a*i + 1/4))/(h*(4*a^2 + 1))

		% wave number of upstream downstream traveling waves
		k = [ k0(:,1) - dk0_dx_rel(:,1), ...
        	     -k0(:,2) - dk0_dx_rel(:,2)];
	else
		% TODO, this is only with respect to depth
		syms dh_dx p2 real
		% TODO this needs to be simplified before differentiation
		% when p2==pi
		dk0_dx_rel = diff(k0,h)/(2*k0)
%		k0 = simplify(subs(k0,p2,Pi))
%		dk0_dx_rel = diff(k0,h)/(2*k0)
		k = k0 - dk0_dx_rel*dh_dx;
	end
	end

	end
                                                                                
	% linearised
	% TODO, make own class
	if (0)
	r_lin = real( (   1i*omega ...
		        + cd*p2./pi.*q0./h.^2 ...
		        + 1/2*cd.*p1./pi.*qt./h.^2 )./sqrt(g.*h) );
	end
end % wave_number_analytic
