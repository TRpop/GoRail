clc
clear all
close all

%% GPS Uncertainty Ellipsoid

gps_u = [-0.5299   0.7256   -0.4389; -0.2442   -0.6262    -0.7404; -0.8122    -0.2852    0.5090];
gps_s = [0.6291  0  0; 0  0.1974  0; 0  0  0.0570];

rand_gps_x = sqrt(gps_s(1,1)); 
rand_gps_y = sqrt(gps_s(2,2)); 
rand_gps_z = sqrt(gps_s(3,3)); 

%% Car model

wheel_radius = 0.1;
length = 1;

velocity = 10; % m/s
heading = 0;

%% kalman set

H = [1 0 0; 0 1 0];
Q = [0.01 0 0;0 0.01 0; 0 0 0.00001];
R = [rand_gps_x 0;0 rand_gps_y];

P = 1*eye(3);

%% Simulation

sample_time = 0.1;
test_time = 10;
step = test_time/sample_time;

save_real = zeros(step,3);
save_gps = zeros(step,3);
save_kal = zeros(step,3);

real_x = 0;   % m
real_y = 0; % m
real_z = 0;     % m
real_theta = pi/4; % rad

x_kal = [real_x; real_y; real_theta]; % kalman init state
kal_theta = real_theta;
pre_kal_theta = kal_theta;

for i=1:step
   
   %% car position
    real_velocity = velocity + velocity/10*randn();
    real_x = real_x + real_velocity*sample_time*cos(real_theta); % m
    real_y = real_y + real_velocity*sample_time*sin(real_theta); % m
    real_theta = real_theta + real_velocity/length*tan(heading)*sample_time; % rad

   %% gps
    gps_x = real_x + rand_gps_x*randn();
    gps_y = real_y + rand_gps_y*randn();
    gps_z = real_z + rand_gps_z*randn();
    gps_xyz = [gps_x gps_y gps_z]';
    gps_xy = [gps_x gps_y]';
    
%    %% steroe
%     [st_u, st_s, st_cx, st_cy, ct_cz, st_x, st_y, st_z] = sj_stereo(real_x,real_y,real_z);
%     st_xyz = [st_x st_y st_z]';
%     
%    %% geometry sensor fusion
%     gps_w = inv(inv(gps_u*gps_s*gps_u')+inv(st_u*st_s*st_u'))*inv(gps_u*gps_s*gps_u');
%     st_w = inv(inv(gps_u*gps_s*gps_u')+inv(st_u*st_s*st_u'))*inv(st_u*st_s*st_u');
%     
%     sf_xyz = gps_w*gps_xyz + st_w*st_xyz;
    
   %% kalman filter
    pre_kal_theta = kal_theta;
    
    A_13 = -velocity*sample_time*sin(pre_kal_theta);
    A_23 = velocity*sample_time*cos(pre_kal_theta);
    A = [1 0 A_13; 0 1 A_23; 0 0 1];
    
    x_kal_p = A*x_kal;
    Pp = A*P*A' + Q;
    K = Pp*H'/(H*Pp*H'+R);
    
    z = gps_xy;%[sf_xyz(1); sf_xyz(2)]; % sensor fusion data
    x_kal = x_kal_p + K*(z - H*x_kal_p);
    
    P = Pp - K*H*Pp;
    
   %% data save
    save_real(i,:) = [real_x real_y real_z];
    save_gps(i,:) = gps_xyz';
    save_kal(i,:) = x_kal;
end

% error_gps = zeros(step,3);
% error_st = zeros(step,3);
% error_sf = zeros(step,3);
% error_kal = zeros(step,3);

error_gps = abs(save_real-save_gps);
error_kal = abs(save_real-save_kal);

plot_x = 1:step;

figure(1)
plot(plot_x, error_gps(:, 1))
hold on
plot(plot_x, error_kal(:, 1))
title('error of x')
xlabel('step')
ylabel('distance x(m)')
legend('gps error of x', 'kalman error of x')
grid on
hold off

figure(2)
plot(plot_x, error_gps(:, 2))
hold on
plot(plot_x, error_kal(:, 2))
title('error of y')
xlabel('step')
ylabel('distance y(m)')
legend('gps error of y', 'kalman error of y')
grid on
hold off

figure(3)
plot(plot_x, error_gps(:, 3))
hold on
plot(plot_x, error_kal(:, 3))
title('error of z')
xlabel('step')
ylabel('distance z(m)')
legend('gps error of z', 'kalman error of z')
grid on
hold off

figure(4)
plot(save_real(:, 1), save_real(:, 2))
hold on
plot(save_gps(:, 1), save_gps(:, 2), '.')
plot(save_kal(:, 1), save_kal(:, 2))
title('Trajactory')
xlabel('distance x(m)')
ylabel('distance y(m)')
legend('real', 'gps', 'kalman')
grid on
hold off


% 
% 

% % figure(2)
% % plot3(save_real(:,1),save_real(:,2),save_real(:,3),'b')
% % hold on
% % plot3(save_gps(:,1),save_gps(:,2),save_gps(:,3),'.m')
% % plot3(save_st(:,1),save_st(:,2),save_st(:,3),'k')
% % plot3(save_sf(:,1),save_sf(:,2),save_sf(:,3),'r')
% % title('Trajectory')
% % xlabel('distance x(m)')
% % ylabel('distance y(m)')
% % zlabel('distance z(m)')
% % legend('real','gps','st','fusion')
% % grid on
% 
% figure(2)
% plot(save_real(:,1),save_real(:,2),'b')
% hold on
% plot(save_gps(:,1),save_gps(:,2),'.m')
% plot(save_st(:,1),save_st(:,2),'k')
% %plot(save_sf(:,1),save_sf(:,2),'r')
% plot(save_kal(:,1),save_kal(:,2),'g')
% title('Trajectory')
% xlabel('distance x(m)')
% ylabel('distance y(m)')
% legend('real','gps','st','kalman')
% grid on
% 
% figure(3)
% plot(plot_x/10,save_real(:,1),'b');
% hold on
% plot(plot_x/10,save_gps(:,1),'m');
% plot(plot_x/10,save_st(:,1),'k');
% plot(plot_x/10,save_sf(:,1),'r');
% plot(plot_x/10,save_kal(:,1),'g');
% title('distance x')
% xlabel('time(s)')
% ylabel('distance x(m)')
% legend('real','gps','st','fusion','kalman')
% 
% figure(4)
% plot(plot_x/10,save_real(:,2),'b');
% hold on
% plot(plot_x/10,save_gps(:,2),'m');
% plot(plot_x/10,save_st(:,2),'k');
% plot(plot_x/10,save_sf(:,2),'r');
% plot(plot_x/10,save_kal(:,2),'g');
% title('distance y')
% xlabel('time(s)')
% ylabel('distance y(m)')
% legend('real','gps','st','fusion','kalman')
% 
% % figure(5)
% % plot(plot_x/10,save_real(:,3),'b');
% % hold on
% % plot(plot_x/10,save_gps(:,3),'m');
% % plot(plot_x/10,save_st(:,3),'k');
% % plot(plot_x/10,save_sf(:,3),'r');
% % title('distance z')
% % xlabel('time(s)')
% % ylabel('distance z(m)')
% % legend('real','gps','st','fusion')
% 
% figure(6)
% plot(plot_x/10,abs(save_real(:,1)-save_gps(:,1)),'m');
% hold on
% plot(plot_x/10,abs(save_real(:,1)-save_st(:,1)),'k');
% plot(plot_x/10,abs(save_real(:,1)-save_sf(:,1)),'r');
% plot(plot_x/10,abs(save_real(:,1)-save_kal(:,1)),'g');
% title('error x')
% xlabel('time(s)')
% ylabel('error x(m)')
% legend('gps','st','fusion','kalman')
% 
% figure(7)
% plot(plot_x/10,abs(save_real(:,2)-save_gps(:,2)),'m');
% hold on
% plot(plot_x/10,abs(save_real(:,2)-save_st(:,2)),'k');
% %plot(plot_x/10,abs(save_real(:,2)-save_sf(:,2)),'r');
% plot(plot_x/10,abs(save_real(:,2)-save_kal(:,2)),'g');
% title('error y')
% xlabel('time(s)')
% ylabel('error y(m)')
% legend('gps','st','kalman')
% 
% % figure(8)
% % plot(plot_x/10,abs(save_real(:,3)-save_gps(:,3)),'m');
% % hold on
% % plot(plot_x/10,abs(save_real(:,3)-save_st(:,3)),'k');
% % plot(plot_x/10,abs(save_real(:,3)-save_sf(:,3)),'r');
% % title('error z')
% % xlabel('time(s)')
% % ylabel('error z(m)')
% % legend('gps','st','fusion')