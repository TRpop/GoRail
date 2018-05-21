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

wheel_radius = 0.07;
length = 1;

velocity = 10; % m/s
heading = 0;

%% kalman set

H = [1 0 0 0; 0 1 0 0; 0 0 0 1];
Q = [1 0 0;0 1 0; 0 0 0.00001];
R = [rand_gps_x 0;0 rand_gps_y];

P = 1*eye(3);

%% Simulation

sample_time = 0.1;
test_time = 10;
step = test_time/sample_time;

save_real = zeros(step,3);
save_dist = zeros(step,3);
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
    save_dist(i, 1) = sqrt(real_x^2 + real_y^2);
    save_dist(i, 2) = sqrt(gps_xyz(1)^2 + gps_xyz(2)^2);
    save_dist(i, 3) = sqrt(x_kal(1)^2 + x_kal(2)^2);
    save_gps(i,:) = gps_xyz';
    save_kal(i,:) = x_kal;
end
%%
error_gps = abs(save_real-save_gps);
error_kal = abs(save_real-save_kal);

error_gps_dist = abs(save_dist(:, 1)-save_dist(:, 2));
error_kal_dist = abs(save_dist(:, 1)-save_dist(:, 3));

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

figure(5)
plot(plot_x, save_dist(:, 1))
hold on
plot(plot_x, save_dist(:, 2), '.')
plot(plot_x, save_dist(:, 3))
title('Distance')
xlabel('step')
ylabel('distance(m)')
legend('real', 'gps', 'kalman')
grid on
hold off

figure(6)
plot(plot_x, error_gps_dist)
hold on
plot(plot_x, error_kal_dist)
title('Error of Distance')
xlabel('step')
ylabel('distance(m)')
legend('gps', 'kal')
grid on
hold off

