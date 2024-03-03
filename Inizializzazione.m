%% Script di inizializzazione progetto ISI Brigida-Calzaretta-Massara

% definiamo i parametri del nostro sistema
% consideriamo come AUV il V-FEDIS alto 150cm, largo 240cm e lungo 330cm, 
% con un peso di circa 580kg
% velocità di 5.5 km/h = 1.527779 m/s

clear all
close all
clc

m = 580;             % massa     [kg]
lx = 3.3;            % lunghezza [m]
lz = 1.5;            % altezza   [m]

% calcoliamo il momento di inerzia considerando il corpo come un elissoide
% e usando la relazione I_y = m*(lx^2 + lz^2)/5

I_y = m*(lx^2 + lz^2)/5;    % momento di inerzia lungo l'asse y

% Definiamo la matrice di Damping (i valori scelti sono stati presi da una
% tesi sulla stima di traiettoria per il V-FIDES)

Damping_matrix = [3 0 0; 0 3 0; 0 0 16];    

% inseriamo i valori che abbiamo calcolato all'interno di una variabile S
% che ci servirà per la comunicazione con simulink

S = [m , I_y, Damping_matrix(1,1), Damping_matrix(2,2), Damping_matrix(3,3)]; 

% Posizione delle boe 

x_boa_1 = 100;
z_boa_1 = 0;
x_boa_2 = 0;
z_boa_2 = 150;


% Condizioni iniziali 

x_0 = rand() * x_boa_1*3/4;
z_0 = rand() * z_boa_2*3/4;

theta_0 = wrapTo2Pi(rand() * (pi/2) - pi);
u_0 = 0;

% Assumiamo che velocita trasversale e velocità angolare al tempo iniziale
% siano nulle
w_0 = 0;
q_0 = 0;

% definiamo le equazioni del nostro sistema

csi_0 = [x_0, z_0, theta_0, u_0, w_0, q_0]';      % vettore di stato a t = 0
syms x z theta u w q tau_du tau_dw tau_dq v_y_alpha_m v_z_m v_theta_m tau_u tau_q  dt m d1 d2 d3 I_y
csi = [x,z,theta,u,w,q]';                         % vettore di stato a t > 0
tau_d = [tau_du, tau_dw, tau_dq];                 % vettore dei disturbi
tau = [tau_u, tau_q];                             % vettore di forze e coppie applicate
symb_state_noise_ing = [csi', tau_d, tau];
x_dot = cos(theta) * u + sin(theta) * w;
z_dot = -sin(theta)*u+ cos(theta)*w;
theta_dot = q;
u_dot = (-m*w*q - d1*u + tau_du + tau_u)/m;
w_dot = (m*u*q - d2*w + tau_dw)/m;
q_dot = (-d3*q + tau_dq + tau_q)/I_y;

csi_dot = [x_dot,z_dot,theta_dot, u_dot, w_dot, q_dot]';

% calcoliamo il vettore di stato in tempo discreto

csi_disc = csi + dt*csi_dot;

% calcoliamo le matrici Jacobiane associate al vettore di stato

F = jacobian(csi_disc,[x,z,theta, u, w, q]); 
D = jacobian(csi_disc,[tau_du, tau_dw, tau_dq]);

% Modello dei sensori
y_alpha_m = atan2(z_boa_2 - z,x) + v_y_alpha_m;            % angolo y alpha misurato
z_m = z + v_z_m;                                            % posizione z misurata
theta_m = theta + v_theta_m;                                % angolo theta misurato

% calcoliamo le matrici Jacobiane associate ai sensori

H = jacobian([y_alpha_m, z_m, theta_m],[x,z,theta, u, w, q]);
M = jacobian([y_alpha_m, z_m, theta_m],[v_y_alpha_m, v_z_m, v_theta_m]);

sensor_model = [y_alpha_m, z_m, theta_m]';                  % vettore del modello del sensore  
v = [v_y_alpha_m, v_z_m, v_theta_m];                        % vettore dei disturbi sui sensori

% Inizializzazione Filtro di Kalman Esteso (EKF)

P_00 = diag([0.1,0.1,0.1,0.1,0.1,0.1]);     % Matrice cov dell'errore di stima

% Definiamo varianze su forze e coppie

var_tau_u = 0.05;
var_tau_w = 0.005;           
var_tau_q = 0.0001;          

% Definiamo la matrice Q = matrice di varianza del rumore bianco

Q = diag([var_tau_u, var_tau_w, var_tau_q]);

% Definiamo le varianze sulle misure dei sensori

var_y_alpha_m= 0.00001;
var_z_m=0.01;
var_theta_m=0.01;


% Definiamo la matrice R = matrice di varianza del disturbo di misura

R =diag([var_y_alpha_m, var_z_m, var_theta_m]);

% Tempo di campionamento di discretizzazione

dt = double(subs(dt,dt, 0.01));

is_particle = true;    % variabile che usiamo per il selettore particle/ext kalman filter
