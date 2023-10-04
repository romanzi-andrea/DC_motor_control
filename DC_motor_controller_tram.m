close all
clear all
clc

% System parameters
s = tf('s');
g = 9.81;
t_sim = 800;

V_supply = 600;
eta = 0.9;
tau_a = 0.01;
tau_e = 1;
Ve_n = 120;
ie_n = 1;
Re = Ve_n/ie_n;
Le = tau_e*Re;
Va_n = 600;   % = max supply voltage
ts = 25;

v_n = 60;       % km/h
w_n = 314;    %rad/s
v_n_mps = v_n/3.6;
r = v_n_mps/w_n;
an = v_n_mps/ts;

%load data
m_tram = 10e3;
m_person = 80;
n_people = 200;
m_people = n_people*m_person;
m_tot = m_tram+m_people;
J = m_tot*r^2;

%Pm = Fn*v_n_mps
Fn = m_tot*an;
Ff = 1/3*Fn;
Ftot = Fn+Ff;
Tn = Fn*r;
Pn = Tn*w_n;
P_mecc = Ftot*v_n_mps;
Pel = P_mecc/eta;
ia_n = Pel/Va_n;
En = P_mecc/ia_n;

Ra = Va_n*(Va_n-En)/(Pel);
La = tau_a*Ra;
K = En/w_n/ie_n;

beta = 1/3*Fn*v_n_mps/w_n^2;   %% Patt = beta*wn^2 = 1/3*Fn*vn

Ge = 1/Re/(1+tau_e*s);
Ga = 1/Ra/(1+tau_a*s); %understand how to compute Ra
Gm = 1/(J*s+beta); 

ts = 25;
vect_x = 0:1:10000;
vel_ref1 = 35/3.6*ones(1,1001);   
vel_ref2 = 60/3.6*ones(1,2000);
vel_ref3 = 60/3.6*ones(1,1000);
vel_ref4 = 75/3.6*ones(1,2000);
vel_ref5 = 60/3.6*ones(1,2000);
vel_ref6 = 60/3.6*ones(1,1000);
vel_ref7 = 35/3.6*ones(1,1000);
vel_ref = [vel_ref1, vel_ref2, vel_ref3, vel_ref4, vel_ref5, vel_ref6,vel_ref7];
w = vel_ref./r;
w_ref = timeseries(w);

s_perc1 = 0;
s_perc2 = 0;
s_perc3 = 0.05;
s_perc4 = 0;
s_perc5 = 0;
s_perc6 = -0.05;
s_perc7 = 0;

slope1 = atan(s_perc1)*ones(1,1001);    % rad
slope2 = atan(s_perc2)*ones(1,2000);
slope3 = atan(s_perc3)*ones(1,1000);
slope4 = atan(s_perc4)*ones(1,2000);
slope5 = atan(s_perc5)*ones(1,2000);
slope6 = atan(s_perc6)*ones(1,1000);
slope7 = atan(s_perc7)*ones(1,1000);
slope_ref = [slope1, slope2, slope3, slope4, slope5, slope6,slope7];

T_load_ref = m_tot*g*r*sin(slope_ref);

T_load = timeseries(T_load_ref);

w_b_ref = (K*ie_n*V_supply/Ra-T_load_ref)*1/(beta+K^2*ie_n^2/Ra);  
w_b = timeseries(w_b_ref);
%% current regulators (zero pole cancellation)
wc_curr_a_des = 1000;
Rca = wc_curr_a_des*1*Ra/s*(1+tau_a*s);
kp_ca =  0.4206*2; 
ki_ca = 42.06*2;
kb_ca = ki_ca/kp_ca;     %anti wind-up back calculation
Lca = Rca*Ga;
bode(Lca);

wc_curr_e_des = 1000;
Rce = 1/s*(1+tau_e*s)*Re*wc_curr_e_des;
kp_ce = 120000;
ki_ce = 120000;
kb_ce = ki_ce/kp_ce;
Lce = Rce*Ge;
figure
bode(Lce)

Fca = Lca/(1+Lca);

wc_w_des = 5/ts;
Rw = 100*(J*s+beta)/s; %= wc_w_des*(beta+J*s)/s;
kp_w = 7325;
ki_w = 97.67;
kb_w = ki_w/kp_w;
Lw = Rw*1/K*Fca*K*Gm;
% feed forward compensator
C = 1/Fca/(1+0.01*s);      %+ a pole in high freq to make it realizable
%we can use this compensator but it will slow down a lot the simulation of
%the system and means another tf so it is better to use a simple gain(given
%that the Fca is characterized by a very fast dynamics we can use simply a
%static gain K = 1) getting very good results(they are pratically the same)

PF = 1/(1+5.3*s); %prefilter to slow the w_ref and get the right acc for the tram

Rv = 1/s*(1+.01*s);    %voltage reg for flux weakening
kp_v = 0.1;
ki_v = 10;
kb_v = ki_v/kp_v;
figure
bode(Lw)
%% simulation
% out = sim('DC_motor_controller_tram_sim.slx'); %use solver ode23t
plot(w_ref);
hold on
plot(out.w)
%% computing the error
figure
plot(out.error)
total_err = 0;
for ind = 1:length(out.error.data)
    total_err = total_err+abs(out.error.data(ind));
end
ave = total_err/length(out.error.data)