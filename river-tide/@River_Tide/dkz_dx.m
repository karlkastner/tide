% Tue 17 Apr 14:47:47 CEST 2018
% Karl Kastner, Berlin
%
%% along channel derivative of the wave number of the tidal surface elevation
%% ignores width variation dh/dx and second order depth variation (d^2h/dx^2)
%% TODO rederive with g symbolic
function dkz_dx = dkz_dx(obj,Q0,w,h,cd,omega,az1,Qt,dh_dx,dw_dx)
	p     = -obj.friction_coefficient_dronkers(alpha);
	p1    = p(:,2);
	p2    = p(:,3);

dkz_dx = [ dh_dx*(((11973633753704223*g^(1/2)*h^(1/2)*(4*pi*W^2*h^2*omega^2 + cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(-8*i))^(5/2)*(384*pi*W^3*cd^2*h^2*omega^4*(Q0*p2 + (Qt*p1)/2)^2 - W^2*cd^3*omega^3*(Q0*p2 + (Qt*p1)/2)^3*256*i + W^4*cd*h^4*omega^5*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*192*i - (26182473019414641*W^5*dh_dx^2*g*h^5*omega^4)/140737488355328 - (8727491006471547*W^5*h^6*omega^6)/8796093022208 + W^4*cd*dh_dx^2*g*h^3*omega^3*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*104*i + 24*pi*W^3*cd^2*dh_dx^2*g*h*omega^2*(Q0*p2 + (Qt*p1)/2)^2))/18014398509481984 + (3991211251234741*g^(1/2)*h^(3/2)*(4*pi*W^2*h^2*omega^2 + cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(-8*i))^(5/2)*(768*pi*W^3*cd^2*h*omega^4*(Q0*p2 + (Qt*p1)/2)^2 + W^4*cd*h^3*omega^5*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*768*i - (130912365097073205*W^5*dh_dx^2*g*h^4*omega^4)/140737488355328 - (26182473019414641*W^5*h^5*omega^6)/4398046511104 + 24*pi*W^3*cd^2*dh_dx^2*g*omega^2*(Q0*p2 + (Qt*p1)/2)^2 + W^4*cd*dh_dx^2*g*h^2*omega^3*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*312*i))/9007199254740992 + pi*g*h*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)*(pi*dh_dx*W^3*h^2*omega^2*i + 6*cd*dh_dx*(Q0*p2 + (Qt*p1)/2)*W^2*omega)*(384*pi*W^3*cd^2*h^2*omega^4*(Q0*p2 + (Qt*p1)/2)^2 - W^2*cd^3*omega^3*(Q0*p2 + (Qt*p1)/2)^3*256*i + W^4*cd*h^4*omega^5*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*192*i - (26182473019414641*W^5*dh_dx^2*g*h^5*omega^4)/140737488355328 - (8727491006471547*W^5*h^6*omega^6)/8796093022208 + W^4*cd*dh_dx^2*g*h^3*omega^3*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*104*i + 24*pi*W^3*cd^2*dh_dx^2*g*h*omega^2*(Q0*p2 + (Qt*p1)/2)^2) + (pi*g*h^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)*(pi*dh_dx*W^3*h^2*omega^2*i + 6*cd*dh_dx*(Q0*p2 + (Qt*p1)/2)*W^2*omega)*(768*pi*W^3*cd^2*h*omega^4*(Q0*p2 + (Qt*p1)/2)^2 + W^4*cd*h^3*omega^5*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*768*i - (130912365097073205*W^5*dh_dx^2*g*h^4*omega^4)/140737488355328 - (26182473019414641*W^5*h^5*omega^6)/4398046511104 + 24*pi*W^3*cd^2*dh_dx^2*g*omega^2*(Q0*p2 + (Qt*p1)/2)^2 + W^4*cd*dh_dx^2*g*h^2*omega^3*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*312*i))/2 + (19956056256173705*pi*W^2*g^(1/2)*h^(5/2)*omega^2*(4*pi*W^2*h^2*omega^2 + cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(-8*i))^(3/2)*(384*pi*W^3*cd^2*h^2*omega^4*(Q0*p2 + (Qt*p1)/2)^2 - W^2*cd^3*omega^3*(Q0*p2 + (Qt*p1)/2)^3*256*i + W^4*cd*h^4*omega^5*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*192*i - (26182473019414641*W^5*dh_dx^2*g*h^5*omega^4)/140737488355328 - (8727491006471547*W^5*h^6*omega^6)/8796093022208 + W^4*cd*dh_dx^2*g*h^3*omega^3*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*104*i + 24*pi*W^3*cd^2*dh_dx^2*g*h*omega^2*(Q0*p2 + (Qt*p1)/2)^2))/2251799813685248 - 4*pi^2*W^2*g*h^3*omega^2*(pi*dh_dx*W^3*h^2*omega^2*i + 6*cd*dh_dx*(Q0*p2 + (Qt*p1)/2)*W^2*omega)*(384*pi*W^3*cd^2*h^2*omega^4*(Q0*p2 + (Qt*p1)/2)^2 - W^2*cd^3*omega^3*(Q0*p2 + (Qt*p1)/2)^3*256*i + W^4*cd*h^4*omega^5*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*192*i - (26182473019414641*W^5*dh_dx^2*g*h^5*omega^4)/140737488355328 - (8727491006471547*W^5*h^6*omega^6)/8796093022208 + W^4*cd*dh_dx^2*g*h^3*omega^3*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*104*i + 24*pi*W^3*cd^2*dh_dx^2*g*h*omega^2*(Q0*p2 + (Qt*p1)/2)^2) + pi^2*W^3*dh_dx*g*h^3*omega^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)*(384*pi*W^3*cd^2*h^2*omega^4*(Q0*p2 + (Qt*p1)/2)^2 - W^2*cd^3*omega^3*(Q0*p2 + (Qt*p1)/2)^3*256*i + W^4*cd*h^4*omega^5*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*192*i - (26182473019414641*W^5*dh_dx^2*g*h^5*omega^4)/140737488355328 - (8727491006471547*W^5*h^6*omega^6)/8796093022208 + W^4*cd*dh_dx^2*g*h^3*omega^3*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*104*i + 24*pi*W^3*cd^2*dh_dx^2*g*h*omega^2*(Q0*p2 + (Qt*p1)/2)^2)*i)/((2778046668940015*g^2*h^4*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)^2*(pi*dh_dx*W^3*h^2*omega^2*i + 6*cd*dh_dx*(Q0*p2 + (Qt*p1)/2)*W^2*omega)^2)/281474976710656 + pi*g*h^3*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(4*i) - 2*pi*W^2*h^2*omega^2)^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)^3) - (((3991211251234741*g^(1/2)*h^(3/2)*(4*pi*W^2*h^2*omega^2 + cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(-8*i))^(5/2)*(384*pi*W^3*cd^2*h^2*omega^4*(Q0*p2 + (Qt*p1)/2)^2 - W^2*cd^3*omega^3*(Q0*p2 + (Qt*p1)/2)^3*256*i + W^4*cd*h^4*omega^5*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*192*i - (26182473019414641*W^5*dh_dx^2*g*h^5*omega^4)/140737488355328 - (8727491006471547*W^5*h^6*omega^6)/8796093022208 + W^4*cd*dh_dx^2*g*h^3*omega^3*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*104*i + 24*pi*W^3*cd^2*dh_dx^2*g*h*omega^2*(Q0*p2 + (Qt*p1)/2)^2))/9007199254740992 + (pi*g*h^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)*(pi*dh_dx*W^3*h^2*omega^2*i + 6*cd*dh_dx*(Q0*p2 + (Qt*p1)/2)*W^2*omega)*(384*pi*W^3*cd^2*h^2*omega^4*(Q0*p2 + (Qt*p1)/2)^2 - W^2*cd^3*omega^3*(Q0*p2 + (Qt*p1)/2)^3*256*i + W^4*cd*h^4*omega^5*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*192*i - (26182473019414641*W^5*dh_dx^2*g*h^5*omega^4)/140737488355328 - (8727491006471547*W^5*h^6*omega^6)/8796093022208 + W^4*cd*dh_dx^2*g*h^3*omega^3*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*104*i + 24*pi*W^3*cd^2*dh_dx^2*g*h*omega^2*(Q0*p2 + (Qt*p1)/2)^2))/2)*((2778046668940015*g^2*h^3*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)^2*(pi*dh_dx*W^3*h^2*omega^2*i + 6*cd*dh_dx*(Q0*p2 + (Qt*p1)/2)*W^2*omega)^2)/70368744177664 + 3*pi*g*h^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(4*i) - 2*pi*W^2*h^2*omega^2)^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)^3 - 8*pi^2*W^2*g*h^4*omega^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(4*i) - 2*pi*W^2*h^2*omega^2)*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)^3 - (2778046668940015*pi*W^2*g^2*h^5*omega^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)*(pi*dh_dx*W^3*h^2*omega^2*i + 6*cd*dh_dx*(Q0*p2 + (Qt*p1)/2)*W^2*omega)^2)/17592186044416 - 24*pi^2*W^2*g*h^4*omega^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(4*i) - 2*pi*W^2*h^2*omega^2)^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)^2 + (pi*W^3*dh_dx*g^2*h^5*omega^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)^2*(pi*dh_dx*W^3*h^2*omega^2*i + 6*cd*dh_dx*(Q0*p2 + (Qt*p1)/2)*W^2*omega)*2778046668940015*i)/70368744177664))/((2778046668940015*g^2*h^4*(W*cd*omega*(Q0*p2 + (Qt*p1)/2)*(8*i) - 4*pi*W^2*h^2*omega^2)^2*(pi*W^3*dh_dx*h^2*omega^2*i + 6*W^2*cd*dh_dx*omega*(Q0*p2 + (Qt*p1)/2))^2)/281474976710656 + pi*g*h^3*(W*cd*omega*(Q0*p2 + (Qt*p1)/2)*(4*i) - 2*pi*W^2*h^2*omega^2)^2*(W*cd*omega*(Q0*p2 + (Qt*p1)/2)*(8*i) - 4*pi*W^2*h^2*omega^2)^3)^2), ...
 -dh_dx*(((11973633753704223*g^(1/2)*h^(1/2)*(4*pi*W^2*h^2*omega^2 + cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(-8*i))^(5/2)*(384*pi*W^3*cd^2*h^2*omega^4*(Q0*p2 + (Qt*p1)/2)^2 - W^2*cd^3*omega^3*(Q0*p2 + (Qt*p1)/2)^3*256*i + W^4*cd*h^4*omega^5*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*192*i - (26182473019414641*W^5*dh_dx^2*g*h^5*omega^4)/140737488355328 - (8727491006471547*W^5*h^6*omega^6)/8796093022208 + W^4*cd*dh_dx^2*g*h^3*omega^3*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*104*i + 24*pi*W^3*cd^2*dh_dx^2*g*h*omega^2*(Q0*p2 + (Qt*p1)/2)^2))/18014398509481984 + (3991211251234741*g^(1/2)*h^(3/2)*(4*pi*W^2*h^2*omega^2 + cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(-8*i))^(5/2)*(768*pi*W^3*cd^2*h*omega^4*(Q0*p2 + (Qt*p1)/2)^2 + W^4*cd*h^3*omega^5*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*768*i - (130912365097073205*W^5*dh_dx^2*g*h^4*omega^4)/140737488355328 - (26182473019414641*W^5*h^5*omega^6)/4398046511104 + 24*pi*W^3*cd^2*dh_dx^2*g*omega^2*(Q0*p2 + (Qt*p1)/2)^2 + W^4*cd*dh_dx^2*g*h^2*omega^3*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*312*i))/9007199254740992 - pi*g*h*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)*(pi*dh_dx*W^3*h^2*omega^2*i + 6*cd*dh_dx*(Q0*p2 + (Qt*p1)/2)*W^2*omega)*(384*pi*W^3*cd^2*h^2*omega^4*(Q0*p2 + (Qt*p1)/2)^2 - W^2*cd^3*omega^3*(Q0*p2 + (Qt*p1)/2)^3*256*i + W^4*cd*h^4*omega^5*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*192*i - (26182473019414641*W^5*dh_dx^2*g*h^5*omega^4)/140737488355328 - (8727491006471547*W^5*h^6*omega^6)/8796093022208 + W^4*cd*dh_dx^2*g*h^3*omega^3*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*104*i + 24*pi*W^3*cd^2*dh_dx^2*g*h*omega^2*(Q0*p2 + (Qt*p1)/2)^2) - (pi*g*h^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)*(pi*dh_dx*W^3*h^2*omega^2*i + 6*cd*dh_dx*(Q0*p2 + (Qt*p1)/2)*W^2*omega)*(768*pi*W^3*cd^2*h*omega^4*(Q0*p2 + (Qt*p1)/2)^2 + W^4*cd*h^3*omega^5*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*768*i - (130912365097073205*W^5*dh_dx^2*g*h^4*omega^4)/140737488355328 - (26182473019414641*W^5*h^5*omega^6)/4398046511104 + 24*pi*W^3*cd^2*dh_dx^2*g*omega^2*(Q0*p2 + (Qt*p1)/2)^2 + W^4*cd*dh_dx^2*g*h^2*omega^3*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*312*i))/2 + (19956056256173705*pi*W^2*g^(1/2)*h^(5/2)*omega^2*(4*pi*W^2*h^2*omega^2 + cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(-8*i))^(3/2)*(384*pi*W^3*cd^2*h^2*omega^4*(Q0*p2 + (Qt*p1)/2)^2 - W^2*cd^3*omega^3*(Q0*p2 + (Qt*p1)/2)^3*256*i + W^4*cd*h^4*omega^5*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*192*i - (26182473019414641*W^5*dh_dx^2*g*h^5*omega^4)/140737488355328 - (8727491006471547*W^5*h^6*omega^6)/8796093022208 + W^4*cd*dh_dx^2*g*h^3*omega^3*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*104*i + 24*pi*W^3*cd^2*dh_dx^2*g*h*omega^2*(Q0*p2 + (Qt*p1)/2)^2))/2251799813685248 + 4*pi^2*W^2*g*h^3*omega^2*(pi*dh_dx*W^3*h^2*omega^2*i + 6*cd*dh_dx*(Q0*p2 + (Qt*p1)/2)*W^2*omega)*(384*pi*W^3*cd^2*h^2*omega^4*(Q0*p2 + (Qt*p1)/2)^2 - W^2*cd^3*omega^3*(Q0*p2 + (Qt*p1)/2)^3*256*i + W^4*cd*h^4*omega^5*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*192*i - (26182473019414641*W^5*dh_dx^2*g*h^5*omega^4)/140737488355328 - (8727491006471547*W^5*h^6*omega^6)/8796093022208 + W^4*cd*dh_dx^2*g*h^3*omega^3*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*104*i + 24*pi*W^3*cd^2*dh_dx^2*g*h*omega^2*(Q0*p2 + (Qt*p1)/2)^2) - pi^2*W^3*dh_dx*g*h^3*omega^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)*(384*pi*W^3*cd^2*h^2*omega^4*(Q0*p2 + (Qt*p1)/2)^2 - W^2*cd^3*omega^3*(Q0*p2 + (Qt*p1)/2)^3*256*i + W^4*cd*h^4*omega^5*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*192*i - (26182473019414641*W^5*dh_dx^2*g*h^5*omega^4)/140737488355328 - (8727491006471547*W^5*h^6*omega^6)/8796093022208 + W^4*cd*dh_dx^2*g*h^3*omega^3*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*104*i + 24*pi*W^3*cd^2*dh_dx^2*g*h*omega^2*(Q0*p2 + (Qt*p1)/2)^2)*i)/((2778046668940015*g^2*h^4*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)^2*(pi*dh_dx*W^3*h^2*omega^2*i + 6*cd*dh_dx*(Q0*p2 + (Qt*p1)/2)*W^2*omega)^2)/281474976710656 + pi*g*h^3*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(4*i) - 2*pi*W^2*h^2*omega^2)^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)^3) - (((3991211251234741*g^(1/2)*h^(3/2)*(4*pi*W^2*h^2*omega^2 + cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(-8*i))^(5/2)*(384*pi*W^3*cd^2*h^2*omega^4*(Q0*p2 + (Qt*p1)/2)^2 - W^2*cd^3*omega^3*(Q0*p2 + (Qt*p1)/2)^3*256*i + W^4*cd*h^4*omega^5*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*192*i - (26182473019414641*W^5*dh_dx^2*g*h^5*omega^4)/140737488355328 - (8727491006471547*W^5*h^6*omega^6)/8796093022208 + W^4*cd*dh_dx^2*g*h^3*omega^3*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*104*i + 24*pi*W^3*cd^2*dh_dx^2*g*h*omega^2*(Q0*p2 + (Qt*p1)/2)^2))/9007199254740992 - (pi*g*h^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)*(pi*dh_dx*W^3*h^2*omega^2*i + 6*cd*dh_dx*(Q0*p2 + (Qt*p1)/2)*W^2*omega)*(384*pi*W^3*cd^2*h^2*omega^4*(Q0*p2 + (Qt*p1)/2)^2 - W^2*cd^3*omega^3*(Q0*p2 + (Qt*p1)/2)^3*256*i + W^4*cd*h^4*omega^5*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*192*i - (26182473019414641*W^5*dh_dx^2*g*h^5*omega^4)/140737488355328 - (8727491006471547*W^5*h^6*omega^6)/8796093022208 + W^4*cd*dh_dx^2*g*h^3*omega^3*((2778046668940015*Q0*p2)/281474976710656 + (2778046668940015*Qt*p1)/562949953421312)*104*i + 24*pi*W^3*cd^2*dh_dx^2*g*h*omega^2*(Q0*p2 + (Qt*p1)/2)^2))/2)*((2778046668940015*g^2*h^3*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)^2*(pi*dh_dx*W^3*h^2*omega^2*i + 6*cd*dh_dx*(Q0*p2 + (Qt*p1)/2)*W^2*omega)^2)/70368744177664 + 3*pi*g*h^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(4*i) - 2*pi*W^2*h^2*omega^2)^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)^3 - 8*pi^2*W^2*g*h^4*omega^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(4*i) - 2*pi*W^2*h^2*omega^2)*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)^3 - (2778046668940015*pi*W^2*g^2*h^5*omega^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)*(pi*dh_dx*W^3*h^2*omega^2*i + 6*cd*dh_dx*(Q0*p2 + (Qt*p1)/2)*W^2*omega)^2)/17592186044416 - 24*pi^2*W^2*g*h^4*omega^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(4*i) - 2*pi*W^2*h^2*omega^2)^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)^2 + (pi*W^3*dh_dx*g^2*h^5*omega^2*(cd*(Q0*p2 + (Qt*p1)/2)*W*omega*(8*i) - 4*pi*W^2*h^2*omega^2)^2*(pi*dh_dx*W^3*h^2*omega^2*i + 6*cd*dh_dx*(Q0*p2 + (Qt*p1)/2)*W^2*omega)*2778046668940015*i)/70368744177664))/((2778046668940015*g^2*h^4*(W*cd*omega*(Q0*p2 + (Qt*p1)/2)*(8*i) - 4*pi*W^2*h^2*omega^2)^2*(pi*W^3*dh_dx*h^2*omega^2*i + 6*W^2*cd*dh_dx*omega*(Q0*p2 + (Qt*p1)/2))^2)/281474976710656 + pi*g*h^3*(W*cd*omega*(Q0*p2 + (Qt*p1)/2)*(4*i) - 2*pi*W^2*h^2*omega^2)^2*(W*cd*omega*(Q0*p2 + (Qt*p1)/2)*(8*i) - 4*pi*W^2*h^2*omega^2)^3)^2)];

end

