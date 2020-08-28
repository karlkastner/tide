% Tue 17 Apr 14:49:20 CEST 2018
% Karl Kastner, Berlin
%
%% along-channel derivative of the wave number of the discharge
%% neglects width variation
%%
%% TODO, rederive with g as variable
function dkq_dx = dkq_dx(obj,Q0,w,h,cd,omega,az1,Qt,dh_dx,dw_dx)
	
	p     = -obj.friction_coefficient_dronkers(alpha);
	p1    = p(:,2);
	p2    = p(:,3);

dkq_dx = [1.0./W.^3.*1.0./sqrt(g).*1.0./h.^(5.0./2.0).*1.0./omega.^2.*1.0./(Q0.*cd.*p2.*2.0i+Qt.*cd.*p1.*1i-pi.*W.*h.^2.*omega).^2.*(W.^2.*dh_dx.*omega.^2.*(Q0.*cd.*p2.*2.0i+Qt.*cd.*p1.*1i-pi.*W.*h.^2.*omega).^2.*sqrt(-W.*omega.*(Q0.*cd.*p2.*2.0i+Qt.*cd.*p1.*1i-pi.*W.*h.^2.*omega)).*6.865874536577549e31i+pi.*W.^5.*dh_dx.^2.*sqrt(g).*h.^(9.0./2.0).*omega.^4.*1.274381380158623e32-pi.^2.*W.^5.*dh_dx.^2.*sqrt(g).*h.^(9.0./2.0).*omega.^4.*6.084722881095501e31+pi.*W.^3.*dh_dx.*h.^2.*omega.^3.*(Q0.*cd.*p2.*2.0i+Qt.*cd.*p1.*1i-pi.*W.*h.^2.*omega).*sqrt(-W.*omega.*(Q0.*cd.*p2.*2.0i+Qt.*cd.*p1.*1i-pi.*W.*h.^2.*omega)).*4.577249691051699e31i+Q0.*W.^4.*cd.*dh_dx.^2.*sqrt(g).*h.^(5.0./2.0).*omega.^3.*p2.*2.548762760317246e32i+Qt.*W.^4.*cd.*dh_dx.^2.*sqrt(g).*h.^(5.0./2.0).*omega.^3.*p1.*1.274381380158623e32i+Q0.^2.*W.^3.*cd.^2.*dh_dx.^2.*sqrt(g).*sqrt(h).*omega.^2.*p2.^2.*2.4338891524382e32+Qt.^2.*W.^3.*cd.^2.*dh_dx.^2.*sqrt(g).*sqrt(h).*omega.^2.*p1.^2.*6.084722881095501e31+pi.*Q0.*W.^4.*cd.*dh_dx.^2.*sqrt(g).*h.^(5.0./2.0).*omega.^3.*p2.*2.4338891524382e32i+pi.*Qt.*W.^4.*cd.*dh_dx.^2.*sqrt(g).*h.^(5.0./2.0).*omega.^3.*p1.*1.2169445762191e32i+Q0.*Qt.*W.^3.*cd.^2.*dh_dx.^2.*sqrt(g).*sqrt(h).*omega.^2.*p1.*p2.*2.4338891524382e32+pi.*Q0.*W.^3.*cd.*dh_dx.*dw_dx.*sqrt(g).*h.^(7.0./2.0).*omega.^3.*p2.*8.112963841460667e31i+pi.*Qt.*W.^3.*cd.*dh_dx.*dw_dx.*sqrt(g).*h.^(7.0./2.0).*omega.^3.*p1.*4.056481920730334e31i).*1.232595164407831e-32, ...
1.0./W.^3.*1.0./sqrt(g).*1.0./h.^(5.0./2.0).*1.0./omega.^2.*1.0./(Q0.*cd.*p2.*2.0i+Qt.*cd.*p1.*1i-pi.*W.*h.^2.*omega).^2.*(W.^2.*dh_dx.*omega.^2.*(Q0.*cd.*p2.*2.0i+Qt.*cd.*p1.*1i-pi.*W.*h.^2.*omega).^2.*sqrt(-W.*omega.*(Q0.*cd.*p2.*2.0i+Qt.*cd.*p1.*1i-pi.*W.*h.^2.*omega)).*-6.865874536577549e31i+pi.*W.^5.*dh_dx.^2.*sqrt(g).*h.^(9.0./2.0).*omega.^4.*1.274381380158623e32-pi.^2.*W.^5.*dh_dx.^2.*sqrt(g).*h.^(9.0./2.0).*omega.^4.*6.084722881095501e31-pi.*W.^3.*dh_dx.*h.^2.*omega.^3.*(Q0.*cd.*p2.*2.0i+Qt.*cd.*p1.*1i-pi.*W.*h.^2.*omega).*sqrt(-W.*omega.*(Q0.*cd.*p2.*2.0i+Qt.*cd.*p1.*1i-pi.*W.*h.^2.*omega)).*4.577249691051699e31i+Q0.*W.^4.*cd.*dh_dx.^2.*sqrt(g).*h.^(5.0./2.0).*omega.^3.*p2.*2.548762760317246e32i+Qt.*W.^4.*cd.*dh_dx.^2.*sqrt(g).*h.^(5.0./2.0).*omega.^3.*p1.*1.274381380158623e32i+Q0.^2.*W.^3.*cd.^2.*dh_dx.^2.*sqrt(g).*sqrt(h).*omega.^2.*p2.^2.*2.4338891524382e32+Qt.^2.*W.^3.*cd.^2.*dh_dx.^2.*sqrt(g).*sqrt(h).*omega.^2.*p1.^2.*6.084722881095501e31+pi.*Q0.*W.^4.*cd.*dh_dx.^2.*sqrt(g).*h.^(5.0./2.0).*omega.^3.*p2.*2.4338891524382e32i+pi.*Qt.*W.^4.*cd.*dh_dx.^2.*sqrt(g).*h.^(5.0./2.0).*omega.^3.*p1.*1.2169445762191e32i+Q0.*Qt.*W.^3.*cd.^2.*dh_dx.^2.*sqrt(g).*sqrt(h).*omega.^2.*p1.*p2.*2.4338891524382e32+pi.*Q0.*W.^3.*cd.*dh_dx.*dw_dx.*sqrt(g).*h.^(7.0./2.0).*omega.^3.*p2.*8.112963841460667e31i+pi.*Qt.*W.^3.*cd.*dh_dx.*dw_dx.*sqrt(g).*h.^(7.0./2.0).*omega.^3.*p1.*4.056481920730334e31i).*1.232595164407831e-32];

end

