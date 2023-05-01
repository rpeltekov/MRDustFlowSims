addpath('~/Downloads/MRI research/Code/BSmag-master/BSmag Core');
close all; clear all; clc;

% %% add coil structure and compute field
% 
% 
% theta_x = 0;
% theta_y = 0;
% 
% I = 0.00025;
% BSmag = BSmag_init();
% r = 0.35e-3;
% res_x = 10e-5; 
% res_y = 10e-5;
% res_z = 10e-5;
% theta_t = linspace(0,4*pi,100);        
% Gamma = [r*cos(theta_t'),r*sin(theta_t'),zeros(size(theta_t,2),1)]; % x,y,z [m,m,m]
% Rx = [1 0 0; 0 cos(theta_x) sin(theta_x); 0 -sin(theta_x) cos(theta_x)];
% Ry = [cos(theta_y) 0 -sin(theta_y);0 1 0; sin(theta_y) 0 cos(theta_y)];
% Gamma_r = Ry*Rx*Gamma';
% dGamma = 1e9; % filament max discretization step [m]
% [BSmag] = BSmag_add_filament(BSmag,Gamma_r',I,dGamma);
% 
% % Field points (where we want to calculate the field)
% x_M = linspace(-1e-3,1e-3,round(2e-3/res_x)); % x [m]
% y_M = linspace(-1e-3,1e-3,round(2e-3/res_y)); % y [m]
% z_M = linspace(-1e-3,1e-3,round(2e-3/res_z)); % z [m]
% 
% % Biot-Savart Integration
% [X_M,Y_M,Z_M]=meshgrid(x_M,y_M,z_M);
% [BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M);
% 
% % Plot B/|B|
% figure(1);
% normB=sqrt(BX.^2+BY.^2+BZ.^2);
% quiver3(X,Y,Z,BX./normB,BY./normB,BZ./normB,'b')
% axis tight

%% add coil structure and compute field


theta_x = 0;
theta_y = 0;

I = 0.00005;
BSmag = BSmag_init();
r = 0.375e-3;
res_x = 10e-5; 
res_y = 10e-5;
res_z = 10e-5;
theta_t = linspace(0,4*pi,100);        
Gamma1 = [r*cos(theta_t'),r*sin(theta_t'),ones(size(theta_t,2),1)*0.4e-6]; % x,y,z [m,m,m]
Gamma2 = [r*cos(theta_t'),r*sin(theta_t'),ones(size(theta_t,2),1)*0.2e-6];
Gamma3 = [r*cos(theta_t'),r*sin(theta_t'),ones(size(theta_t,2),1)*0];
Gamma4 = [r*cos(theta_t'),r*sin(theta_t'),-ones(size(theta_t,2),1)*0.2e-6];
Gamma5 = [r*cos(theta_t'),r*sin(theta_t'),-ones(size(theta_t,2),1)*0.4e-6];
Rx = [1 0 0; 0 cos(theta_x) sin(theta_x); 0 -sin(theta_x) cos(theta_x)];
Ry = [cos(theta_y) 0 -sin(theta_y);0 1 0; sin(theta_y) 0 cos(theta_y)];
% Gamma_r = Ry*Rx*Gamma';
dGamma = 5e9; % filament max discretization step [m]


%[BSmag] = BSmag_add_filament(BSmag,Gamma_r',I,dGamma);
[BSmag] = BSmag_add_filament(BSmag,Gamma1,I,dGamma);
[BSmag] = BSmag_add_filament(BSmag,Gamma2,I,dGamma);
[BSmag] = BSmag_add_filament(BSmag,Gamma3,I,dGamma);
[BSmag] = BSmag_add_filament(BSmag,Gamma4,I,dGamma);
[BSmag] = BSmag_add_filament(BSmag,Gamma5,I,dGamma);
% Field points (where we want to calculate the field)
x_M = linspace(-1e-3,1e-3,round(2e-3/res_x)); % x [m]
y_M = linspace(-1e-3,1e-3,round(2e-3/res_y)); % y [m]
z_M = linspace(-1e-3,1e-3,round(2e-3/res_z)); % z [m]

% Biot-Savart Integration
[X_M,Y_M,Z_M]=meshgrid(x_M,y_M,z_M);
[BSmag,X,Y,Z,BX,BY,BZ] = BSmag_get_B(BSmag,X_M,Y_M,Z_M);

% Plot B/|B|
figure(1);
normB=sqrt(BX.^2+BY.^2+BZ.^2);
quiver3(X,Y,Z,BX./normB,BY./normB,BZ./normB,'b')
axis tight

[r,c,h] = size(BX);
s_total = r*c*h;
%% Steady state simulation

dt = 4e-6; % 4 us sampling rate
gamma_bar = 4258; % gamma_bar = gamma / 2pi

rf = 60/360/(gamma_bar*dt); 
taper_length = 14;
tapered_RF = linspace(0,rf,taper_length);

[r,c,h] = size(BX);
s_total = r*c*h;

% tapered RF pulse
rep_time = 0.1e-3; % s
TR = 2e-3;
RF_reptime = int32(TR/rep_time);
RF_dur = int32(TR/dt);
reps = 1000;
total_reps = reps*RF_reptime;
total_dur = total_reps*(rep_time/dt);

b1_taper = zeros(taper_length*RF_dur,1);
b1_taper(1:2*RF_dur:((taper_length-1)*RF_dur+1)) = tapered_RF(1:2:(taper_length-1));
b1_taper(RF_dur:2*RF_dur:(taper_length*RF_dur+1)) = -tapered_RF(2:2:taper_length);
b1_ss = repmat([rf;zeros(RF_dur-1,1);-rf;zeros(RF_dur-1,1)],(reps-taper_length)/2,1);
b1 = [b1_taper; b1_ss];
figure; subplot(2,1,1);
plot(double(1:total_dur)*dt,b1,'LineWidth',2); 
xlabel('Time','FontSize',20);
ylabel('RF magnitude','FontSize',20);
title('Tapered RF pulse','FontSize',20);


% Simulation
t1 = 0.832;
t2 = 0.08;

mx_0 = 0;
my_0 = 0;
mz_0 = 1;

df = 0;
dp = 0;
gr = zeros(size(b1,1),1);
[mx,my,mz] = bloch(b1,gr,dt,t1,t2,df,dp,2,mx_0,my_0,mz_0);


% Plot a spin in the voxel
subplot(2,1,2); plot(double(1:total_dur)*dt,abs(mx), 'LineWidth',2); hold on
plot(double(1:total_dur)*dt,abs(my),'LineWidth',2);
plot(double(1:total_dur)*dt,abs(mz),'LineWidth',2);
legend('Mx', 'My', 'Mz','FontSize',20);
xlabel('Time','FontSize',20);
ylabel('Spin magnitude','FontSize',20);
title('Spin behavior in steady state','FontSize',20);

%%
% % Capture spin behavior at middle and end of TR
% figure; 
% % middle of TR
% subplot(1,2,1);
% plot(TR/2:TR:(reps*TR-TR/2),abs(mx((RF_dur/2):RF_dur:(RF_dur*reps-RF_dur/2))),'LineWidth',2);
% hold on;
% plot(TR/2:TR:(reps*TR-TR/2),abs(my((RF_dur/2):RF_dur:(RF_dur*reps-RF_dur/2))),'LineWidth',2);
% plot(TR/2:TR:(reps*TR-TR/2),abs(mz((RF_dur/2):RF_dur:(RF_dur*reps-RF_dur/2))),'LineWidth',2);
% legend('Mx', 'My', 'Mz');
% xlabel('Time');
% ylabel('Spin magnitude');
% title('Spin behavior in steady state middle of TR');
% 
% % end of TR
% subplot(1,2,2);
% plot(TR:TR:(reps*TR),abs(mx(RF_dur:RF_dur:(RF_dur*reps))),'LineWidth',2);
% hold on;
% plot(TR:TR:(reps*TR),abs(my(RF_dur:RF_dur:(RF_dur*reps))),'LineWidth',2);
% plot(TR:TR:(reps*TR),abs(mz(RF_dur/2:RF_dur:(RF_dur*reps))),'LineWidth',2);
% legend('Mx', 'My', 'Mz');
% xlabel('Time');
% ylabel('Spin magnitude');
% title('Spin behavior in steady state end of TR');


% % Compute total signal
% mss_kspace = mx+1i*my;
% sig_ss = mss_kspace*s_total;
% figure; subplot(1,2,1);
% plot(double(1:total_dur)*dt,abs(sig_ss),'LineWidth',2);
% xlabel('Time');
% ylabel('Magnitude');
% title('Mxy of voxel');
% subplot(1,2,2); 
% plot(TR:TR:(reps*TR),abs(sig_ss(RF_dur/2:RF_dur:(RF_dur*reps))),'LineWidth',2);
% hold on; plot(TR:TR:(reps*TR),abs(sig_ss(RF_dur:RF_dur:(RF_dur*reps))),'LineWidth',2);
% legend('Middle','End');
% xlabel('Time');
% ylabel('Magnitude');
% title('Mxy of voxel');


%% Start from steady state and test alternating field
mx_0_p = mx(end)*ones(1,s_total);
my_0_p = my(end)*ones(1,s_total);
mz_0_p = mz(end)*ones(1,s_total);
mx_0_s = mx(end);
my_0_s = my(end);
mz_0_s = mz(end);

BZ_f = BZ(:);
df = BZ_f*1e4*gamma_bar;
sig_perts = [];
sig_ss = [];

pert_freq = 1;
pert_period = 1/pert_freq;
t_num = pert_period/TR;
for i = 1:t_num*4
    if (mod(i,2) == 1)
        field = -df;
        b = [rf; zeros(RF_dur-1,1)];
    else
        field = df;
        b = [-rf; zeros(RF_dur-1,1)];
    end
    
    if ((mod(floor(i/t_num),2))==1)
        field = zeros(s_total,1); 
    end
    
    gr = zeros(size(b,1),1);
    [mx_p,my_p,mz_p] = bloch(b,gr,dt,t1,t2,field,dp,2,mx_0_p,my_0_p,mz_0_p);
    sig_p = mx_p+1i*my_p;
    sig_p_sum = sum(sig_p,2);
    sig_perts = [sig_perts; sig_p_sum(RF_dur/2)];
    mx_0_p = mx_p(end,:);
    my_0_p = my_p(end,:);
    mz_0_p = mz_p(end,:);
    
    [mx_s,my_s,mz_s] = bloch(b,gr,dt,t1,t2,0,dp,2,mx_0_s,my_0_s,mz_0_s);
    sig_s = mx_s+1i*my_s;
    sig_s_sum = sig_s * s_total;
    sig_ss = [sig_ss; sig_s_sum(RF_dur/2)];
    mx_0_s = mx_s(end,:);
    my_0_s = my_s(end,:);
    mz_0_s = mz_s(end,:);    
end

figure;
plot((1:length(sig_perts))*(dt*double(RF_dur)),abs(sig_perts),'LineWidth',1.5); hold on;
plot((1:length(sig_ss))*(dt*double(RF_dur)),abs(sig_ss),'LineWidth',1.5);
xlabel('Time (s)', 'FontSize',20);
ylabel('Signal Intensity', 'FontSize',20);
legend('with perturbation','steady state','FontSize',20);

figure; plot((1:length(sig_perts))*(dt*double(RF_dur)),abs(sig_perts)./abs(sig_ss),'LineWidth',2);
xlabel('Time (s)','FontSize',20);
ylabel('Signal Ratio','FontSize',20);
set(gca,'FontSize',18,'LineWidth',2);

% RF_dur = 500;
% perts1hz = load('1Hz_perts.mat');
% perts1hz = perts1hz.sig_perts;
% perts2hz = load('2Hz_perts.mat');
% perts2hz = perts2hz.sig_perts;
% perts5hz = load('5Hz_perts.mat');
% perts5hz = perts5hz.sig_perts;
% perts10hz = load('10Hz_perts.mat');
% perts10hz = perts10hz.sig_perts;
% perts20hz = load('20Hz_perts.mat');
% perts20hz = perts20hz.sig_perts;
% 
% figure; 
% plot((1:length(perts1hz))*(dt*double(RF_dur)),abs(perts1hz)/s_total,'LineWidth',1.5); hold on;
% plot((1:length(perts1hz))*(dt*double(RF_dur)),abs(perts2hz)/s_total,'LineWidth',1.5); 
% plot((1:length(perts1hz))*(dt*double(RF_dur)),abs(perts5hz)/s_total,'LineWidth',1.5); 
% plot((1:length(perts1hz))*(dt*double(RF_dur)),abs(perts10hz)/s_total,'LineWidth',1.5); 
% plot((1:length(perts1hz))*(dt*double(RF_dur)),abs(perts20hz)/s_total,'LineWidth',1.5); 
% xlabel('Time (s)');
% ylabel('Normalized signal intensity');
% set(gca,'FontSize',18,'LineWidth',2);
% legend('1Hz','2Hz','5Hz','10Hz','20Hz');
% 
% 
% 
