% Author: Ted Simmons
% Date: October 13
% T


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IMPORTANT MODEL Variables %%%%

% Physical Parameters:
M_pl = 0.23; % kg
M_c = 1.0731; % note this different from Mc2 of 0.57 kg %kg Mass of cart with 3 cable connectors
l_pl = 0.3302; % m
Ip = 7.88e-3;%1/12 * (M_pl * l_pl.^2);
g = 9.81; % m/s^2
I_pl = 7.88e-03; %kgm^2
Bp = 0.0024; % Nms/rad damping force on the pendulum
cap_Lp = 0.6413;
% Conversion Parameters ( From setup_ip01_2_configuration.m): 
K_R2D = 180 / pi; % from radians to degrees
K_D2R = 1 / K_R2D; % from degrees to radians
K_IN2M = 0.0254; % from Inch to Meter
K_M2IN = 1 / K_IN2M; % from Meter to Inch
K_RDPS2RPM = 60 / ( 2 * pi ); % from rad/s to RPM
K_RPM2RDPS = 1 / K_RDPS2RPM; % from RPM to rad/s
K_OZ2N = 0.2780139; % from oz-force to N
K_N2OZ = 1 / K_OZ2N; % from oz-force to N

% Motor Parameters:
Rm = 2.6; % Ohms
Kt =  1.088 * K_OZ2N * K_IN2M; % (Nm/A)
Km = 0.804e-3 * K_RDPS2RPM; % (Vs/rad)
Kg = 3.71; % (dim)
r_mp = 0.5 / 2 * K_IN2M; 
Eff_g = 1;
Eff_m = 1; 
IMAX_UPM = 5;

% Simplified parameters:
Mc = M_c;
Mp = M_pl;% kg;
lp = l_pl;
Beq = 5.4; % this other damping coefficient
IC_ALPHA0 = 0;
motor_v2f = Eff_g*Eff_m*Kg*Kt/(Rm*r_mp);
motor_s2f =  Eff_g*Eff_m*Kg*Kg*Kt*Km/(Rm*r_mp*r_mp);
motor_f2v = 1/motor_v2f;
motor_s2v =  -motor_s2f/motor_v2f;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% ABCD EQUATIONS: (VOLTAGE -> Position, Angle %%%%%%%%%%%%%%%%%%%%%%

Av( 1, 1 ) = 0;
Av( 1, 2 ) = 0;
Av( 1, 3 ) = 1;
Av( 1, 4 ) = 0;
Av( 2, 1 ) = 0;
Av( 2, 2 ) = 0;
Av( 2, 3 ) = 0;
Av( 2, 4 ) = 1;
Av( 3, 1 ) = 0;
Av( 3, 2 ) = Mp^2*lp^2*g/(Mc*Ip+Mc*Mp*lp^2+Mp*Ip);
Av( 3, 3 ) = -(Ip*Eff_g*Kg^2*Eff_m*Kt*Km+Ip*Beq*Rm*r_mp^2+Mp*lp^2*Eff_g*Kg^2*Eff_m*Kt*Km+Mp*lp^2*Beq*Rm*r_mp^2)/Rm/r_mp^2/(Mc*Ip+Mc*Mp*lp^2+Mp*Ip);
Av( 3, 4 ) = -Mp*lp*Bp/(Mc*Ip+Mc*Mp*lp^2+Mp*Ip);
Av( 4, 1 ) = 0;
Av( 4, 2 ) = Mp*g*lp*(Mc+Mp)/(Mc*Ip+Mc*Mp*lp^2+Mp*Ip);
Av( 4, 3 ) = -Mp*lp*(Eff_g*Kg^2*Eff_m*Kt*Km+Beq*Rm*r_mp^2)/Rm/r_mp^2/(Mc*Ip+Mc*Mp*lp^2+Mp*Ip);
Av( 4, 4 ) = -Bp*(Mc+Mp)/(Mc*Ip+Mc*Mp*lp^2+Mp*Ip);

Bv( 1, 1 ) = 0;
Bv( 2, 1 ) = 0;
Bv( 3, 1 ) = Eff_g*Kg*Eff_m*Kt/r_mp*(Ip+Mp*lp^2)/Rm/(Mc*Ip+Mc*Mp*lp^2+Mp*Ip);
Bv( 4, 1 ) = Eff_g*Kg*Eff_m*Kt/r_mp*Mp*lp/Rm/(Mc*Ip+Mc*Mp*lp^2+Mp*Ip);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%% ABCD EQUATIONS: (Force-> Position, Angle %%%%%%%%%%%%%%%%%
fden = (Mc + Mp) * Ip + Mc*Mp*lp^2;

Af( 1, 1 ) = 0;
Af( 1, 2 ) = 0;
Af( 1, 3 ) = 1;
Af( 1, 4 ) = 0;
Af( 2, 1 ) = 0;
Af( 2, 2 ) = 0;
Af( 2, 3 ) = 0;
Af( 2, 4 ) = 1;
Af( 3, 1 ) = 0;
Af( 3, 2 ) = Mp^2 * lp^2 * g /fden;
Af( 3, 3 ) = -Beq*(Ip + Mp*lp^2)/fden;
Af( 3, 4 ) = -Mp*lp*Bp/fden;
Af( 4, 1 ) = 0;
Af( 4, 2 ) = (Mc+Mp)*Mp*g*lp/fden;
Af( 4, 3 ) = -Mp*lp*Beq/fden;
Af( 4, 4 ) = -(Mc+Mp)*Bp/fden;

Bf( 1, 1 ) = 0;
Bf( 2, 1 ) = 0;
Bf( 3, 1 ) = (Ip+Mp*lp^2)/fden;
Bf( 4, 1 ) = Mp*lp/fden;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% CONVERTING FROM STATE-SPACE TO TF %%%%%%%%%%%%%%%%%%%%%%%%%%%
C_v2x = [1 0 0 0];
C_v2a = [0 1 0 0];
C_v2s = [0 0 1 0];
C_v2w = [0 0 0 1];

C_f2x = [1 0 0 0];
C_f2a = [0 1 0 0];
C_f2s = [0 0 1 0];
C_f2w = [0 0 0 1];

C_f2p = [1 cap_Lp 0 0];

[N1, D1] = ss2tf(Av,Bv,C_v2x, 0);
G_v2x = tf(N1, D1); 
[N2, D2] = ss2tf(Av,Bv,C_v2a, 0);
G_v2a = tf(N2, D2);
[N3, D3] = ss2tf(Av,Bv,C_f2s, 0);
G_v2s = tf(N3, D3);
[N4, D4] = ss2tf(Av,Bv,C_f2w, 0);
G_v2w = tf(N4, D4);
[N5, D5] = ss2tf(Af,Bf,C_f2x, 0);
G_f2x = tf(N5, D5);
[N6, D6] = ss2tf(Af,Bf,C_f2a, 0);
G_f2a = tf(N6, D6);
[N7, D7] = ss2tf(Af,Bf,C_f2s, 0);
G_f2v = tf(N7, D7);
[N8, D8] = ss2tf(Af,Bf,C_f2w, 0);
G_f2w = tf(N8, D8);

[N9, D9] = ss2tf(Af,Bf,C_f2p, 0);
G_f2p = tf(N9, D9);

zp_f2a = zpk(G_f2a);
zp_f2x = zpk(G_f2x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Designing a Controller with Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%


% Plotting Root Locus of Angle and Position
a_ref = 0;
x_ref = 0;
p_ref = 0;

% Angle Controller
Kap = 5000;
Kad = 20;
Kai = 3000;
Ka_G = 1;
zac1 = 10;
zac2 = 9;
pac1 = 0;
G_ac = Ka_G * tf([Kad Kap Kai], [1 0]); %tf([1 zac1], 1) * tf([1 zac2], 1) * tf(1, [1 pac1]); %Ka_G * tf([Kad Kap Kai], [1 0]);
G_f2a_comp = G_ac * G_f2a;
zp_Gac = zpk(G_ac);

% Position Controller
Kxp = -2000;
Kxd = -20;
Kxi = -100;
Kx_G = 1;
G_xc = Kx_G * tf([Kxd Kxp Kxi], [1 0]);
G_f2x_comp = G_xc * G_f2x;
zp_Gxc = zpk(G_xc);

% Tip Position Controller
Kpp = 50000000;
Kpd = -300000;
Kpi = 1000000;
Kp_G = 1;
G_pc = Kp_G * tf([Kpd Kpp Kpi], [1 0]);
G_f2p_comp = G_pc * G_f2p;
zp_Gpc = zpk(G_pc);

% Angle Voltage Controller


% Position Voltage Controller


% LQR Voltage Controller

GAIN_LQR = [-50 180.2 -50.3 28.49];




G_xLQR = tf([GAIN_LQR(3) GAIN_LQR(1)], 1);
G_aLQR =  tf([GAIN_LQR(4) GAIN_LQR(2)], 1);
lqr_cut = 100;
G_xLQR_ip = imperfectTF(G_xLQR, lqr_cut);
G_aLQR_ip = imperfectTF(G_aLQR, lqr_cut);

G_xLQR_comp = G_xLQR * G_v2x;
G_aLQR_comp = G_aLQR * G_v2a;

% PD tip controller
GAIN_tPD = [20 2000 0];

G_tPD = tf([ GAIN_tPD(1) GAIN_tPD(2) GAIN_tPD(3)], [1 0]);

G_tPD_ip = imperfectTF(G_tPD, lqr_cut);

GAIN_tPD2 = [20 1000 0];

G_tPD2 = tf([ GAIN_tPD2(1) GAIN_tPD2(2) GAIN_tPD2(3)], [1 0]);

G_tPD_ip2 = imperfectTF(G_tPD2, lqr_cut);

lengthCoupling = 0.66;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Calculating Proper-Imperfect Transfer Functions %%%%%%%%%%

fcut = 1000;
% Angle
G_ac_ip = imperfectTF(G_ac, fcut);
% Position 
G_xc_ip = imperfectTF(G_xc, fcut);


% Tip Position Only
G_pc_ip = imperfectTF(G_pc, fcut);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(1);
rlocus(G_f2a);
title("Force-Angle Root Locus");
figure(2);
rlocus(G_f2a_comp);
title("Compensated Force-Angle Root Locus");

figure(3);
rlocus(G_f2x);
title("Force-Position Root Locus");
figure(4);
rlocus(G_f2x_comp);
title("Compensated Force-Position Root Locus");



figure(5);
rlocus(G_f2p);
title("Force - Tip Position Root Locus");
figure(6);
rlocus(G_f2p_comp);
title("Compensated Force - Tip Position Root Locus");

%}
%{
figure(7);
rlocus(G_v2a);
title("Voltage-Angle Root Locus");

figure(8);
rlocus(G_aLQR_comp);
title("LQR Compensated Voltage-Angle Root Locus");



figure(9);
rlocus(G_v2x);
title("Voltage-Position Root Locus");

figure(10);
rlocus(G_xLQR_comp);
title("LQR Compensated Voltage-Position Root Locus");

%}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Other Design Tools: %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Creating Routh-Hurwitz Array:
%{
a1 = 11.74;
a2 = -25.26;
a3 = -263.4;
b1 = 3.527;
b2 = 6.298e-15;
b3 = -1.12e-15;
syms kd kp
[RouthArray, RouthColumn] = routh_hurwitz([1 (a1+kd*b1) (a2+kd*b2+kp*b1) (a3+kp*b2) (kp*b3) 0], 3);
vpaRC = vpa(RouthColumn, 3);
%}

% Using PID Tuning tool:
%pidTuner(G_ac, 'pid');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ka_i = 0;
Kx_i = 0;
my_Gx = G_xLQR + tf(Kx_i, [1 0]);
my_Ga = G_aLQR + tf(Ka_i, [1 0]);
again = 1;
xgain = 1;
%%%% Controller Variables for Simulink TF Implementation %%%%%%%%
[ua_num, ua_den] = tfdata(my_Ga, 'v');

ua_long_num = zeros(1, 11);

for i = 1:length(ua_num)
    ua_long_num(i) = ua_num(length(ua_num) - i + 1);
end

Avec1 = ua_long_num(1);
Avec2 = ua_long_num(2);
Avec3 = ua_long_num(3);
Avec4 = ua_long_num(4);
Avec5 = ua_long_num(5);
Avec6 = ua_long_num(6);
Avec7 = ua_long_num(7);
Avec8 = ua_long_num(8);
Avec9 = ua_long_num(9);
Avec10 = ua_long_num(10);
Avec11 = ua_long_num(11);

ua_long_den = zeros(1, 11);

for i = 1:length(ua_den)
    ua_long_den(i) = ua_den(length(ua_den) - i + 1);
end


Bvec1 = ua_long_den(1);
Bvec2 = ua_long_den(2);
Bvec3 = ua_long_den(3);
Bvec4 = ua_long_den(4);
Bvec5 = ua_long_den(5);
Bvec6 = ua_long_den(6);
Bvec7 = ua_long_den(7);
Bvec8 = ua_long_den(8);
Bvec9 = ua_long_den(9);
Bvec10 = ua_long_den(10);
Bvec11 = ua_long_den(11);

[ux_num, ux_den] = tfdata(my_Gx, 'v');

ux_long_num = zeros(1, 11);

for i = 1:length(ux_num)
    ux_long_num(i) = ux_num(length(ux_num) - i + 1);
end

Cvec1 = ux_long_num(1);
Cvec2 = ux_long_num(2);
Cvec3 = ux_long_num(3);
Cvec4 = ux_long_num(4);
Cvec5 = ux_long_num(5);
Cvec6 = ux_long_num(6);
Cvec7 = ux_long_num(7);
Cvec8 = ux_long_num(8);
Cvec9 = ux_long_num(9);
Cvec10 = ux_long_num(10);
Cvec11 = ux_long_num(11);

ux_long_den = zeros(1, 11);

for i = 1:length(ux_den)
    ux_long_den(i) = ux_den(length(ux_den) - i + 1);
end

Dvec1 = ux_long_den(1);
Dvec2 = ux_long_den(2);
Dvec3 = ux_long_den(3);
Dvec4 = ux_long_den(4);
Dvec5 = ux_long_den(5);
Dvec6 = ux_long_den(6);
Dvec7 = ux_long_den(7);
Dvec8 = ux_long_den(8);
Dvec9 = ux_long_den(9);
Dvec10 = ux_long_den(10);
Dvec11 = ux_long_den(11);

w_cut = 2*pi*5;
damping_ratio = 0.9;

derivative_filter = tf([w_cut * w_cut 0], [1 2*w_cut*damping_ratio w_cut*w_cut]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%