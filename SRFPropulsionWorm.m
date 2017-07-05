
clc; clear all; close all;
%% Bjerknes force Fb
format long;
rou = 998; %kg/m^3
rou2 = 5e-9/(1e-3*pi*((30e-6)^2));
R01 = 283e-6/2;
L = 20e-6:1e-6:1e-3; %worm length
D = L.*80e-6/1e-3; %worm diameter
V = L.*pi.*D.^2./4;  %worm volume
R02 = (V.*3./4./pi).^(1/3); %equivalent worm size
d = R01+D/2; %equal bubble radius + worm radius
omega = 32e3*2*pi; %excitation frequency
zeta1 = 0.06e-6; %bubble oscillation amplitude, fitting parameter
zeta2 = 0.16e-6;
zeta3 = 0.25e-6;

Fb1 = -(6*rou*((rou - rou2)/(rou + (2*rou2)))*(R01^4)*(omega^2).*(zeta1.^2)).*((2*L.*((L.^2) ...
+ 6.*(d.^2)))./(3.*(d.^4).*(((L.^2)+ 4.*(d.^2)).^(3/2)))).*(pi/4.*D.^2);
Fb2 = -(6*rou*((rou - rou2)/(rou + (2*rou2)))*(R01^4)*(omega^2).*(zeta2.^2)).*((2*L.*((L.^2) ... 
+ 6.*(d.^2)))./(3.*(d.^4).*(((L.^2)+ 4.*(d.^2)).^(3/2)))).*(pi/4.*D.^2);
Fb3 = -(6*rou*((rou - rou2)/(rou + (2*rou2)))*(R01^4)*(omega^2).*(zeta3.^2)).*((2*L.*((L.^2) ...
+ 6.*(d.^2)))./(3.*(d.^4).*(((L.^2)+ 4.*(d.^2)).^(3/2)))).*(pi/4.*D.^2);

%% drag force Fd
Q = 0.08; % mL/m
A = 3e-3*1e-3;
v = Q*1e-6/60./A; % flow velocity
muo = 1e-3;
lambda = 0.75.*L; % worm wavelengh
Ct = 2*pi./(log((0.18.*lambda)./(D/2))); %
Cn = 4*pi./ ((log((0.18.*lambda)./(D/2)))+0.5);
Fd_min = L.*muo.*Ct.*v;
Fd_max = L.*muo.*Cn.*v;

%% propulsive force Fp
f = 1.5; % worm swimming frequency Hz
b = L./2; % worm swimming amplitude

Fp = (muo*Ct.*(2*(pi^2)*f.*(b.^2)./lambda.*(Cn./Ct-1)./(1 + 2*(pi^2).*(Cn./Ct).*(b./lambda).^2))).*L;

figure (1)
plot(L, Fb1,L, Fb2,L, Fb3, L, Fp); xlabel('length (m)'); ylabel('Force (N)')
legend('Bjerknes1','Bjerknes2','Bjerknes3','Propulsion')

figure (2)
plot(L, Fb3,L, Fd_min, L, Fd_max,L, Fp); xlabel('length (m)'); ylabel('Force (N)')
legend('Fb3','Drag min','Drag max','propulsion')
