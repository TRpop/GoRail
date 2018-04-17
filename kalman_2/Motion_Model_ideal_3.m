clc
clear all
close all

wheel_radius = 0.165;
width = 1.3;

x = 0;
y = 0;
theta = 0;
x_s = 0;
y_s = 0;
theta_s = 0;

vari_3 = 0.6;

w_r = 21;
w_l = 21;
dt = 0.1;
t = 50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 0.98;
b = 0.05;

x_gps_low = 0;
y_gps_low = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = [1 0 0; 0 1 0];
Q = [0.001 0 0;0 0.001 0; 0 0 0.0000001];
R = [0.36 0;0 0.36];
x_kal = [0; 0; 0];
P = 1*eye(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 0:dt:t;
Nsamples = length(n);
Nsamples_1 = Nsamples;

x_save = zeros(Nsamples,1);
y_save = zeros(Nsamples,1);
theta_save = zeros(Nsamples,1);

x_s_save = zeros(Nsamples,1);
y_s_save = zeros(Nsamples,1);
theta_s_save = zeros(Nsamples,1);

x_gps_save = zeros(Nsamples,1);
y_gps_save = zeros(Nsamples,1);

x_gps_mov_save = zeros(Nsamples,1);
y_gps_mov_save = zeros(Nsamples,1);

x_comp_save = zeros(Nsamples,1);
y_comp_save = zeros(Nsamples,1);

x_kal_save = zeros(Nsamples,1);
y_kal_save = zeros(Nsamples,1);

x_v = 1:1:Nsamples;
y_v_r = zeros(Nsamples,1);
y_v_l = zeros(Nsamples,1);
y_ideal_v_r = zeros(Nsamples,1);
y_ideal_v_l = zeros(Nsamples,1);
y_v = zeros(Nsamples,1);
y_w = zeros(Nsamples,1);
y_ideal_v = zeros(Nsamples,1);
y_ideal_w = zeros(Nsamples,1);

for i = 1: 1: Nsamples
    pre_theta = theta;
    [d_x, d_y, d_w, ideal_v_r, ideal_v_l, ideal_v, ideal_w] = Motion_Model_dt(w_r, w_l, wheel_radius ,width, theta,dt);
    x = x + d_x.*dt;
    y = y + d_y.*dt;
    theta = theta + d_w*dt;
    
    x_save(i) = x;
    y_save(i) = y;
    theta_save(i) = theta;
      
    [d_x, d_y, d_w, v_r, v_l, v] = Motion_Model_dt_error(w_r, w_l, wheel_radius ,width, theta_s,dt);
    x_s = x_s + d_x*dt;
    y_s = y_s + d_y*dt;
    theta_s = theta_s + d_w*dt;
    
  
    x_s_save(i) = x_s;
    y_s_save(i) = y_s;
    theta_s_save(i) = theta_s;
    
    x_gps = x_s + normrnd(0,vari_3);
    y_gps = y_s + normrnd(0,vari_3);
   
    x_gps_save(i) = x_gps;
    y_gps_save(i) = y_gps;
    
    A_13 = -ideal_v*dt*sin(pre_theta+ideal_w*dt/2);
    A_23 = ideal_v*dt*cos(pre_theta+ideal_w*dt/2);
    A = [1 0 A_13; 0 1 A_23; 0 0 1];
    x_kal_p = A*x_kal;
    Pp = A*P*A' + Q;
    K = Pp*H'*inv(H*Pp*H'+R);
    
    z = [x_gps; y_gps];
    x_kal = x_kal_p + K*(z - H*x_kal_p);
    
    P = Pp - K*H*Pp;
    
    x_kal_save(i) =x_kal(1);
    y_kal_save(i) =x_kal(2);

    x_gps_low = a*x_gps_low + (1-a)*x_gps;
    y_gps_low = a*y_gps_low + (1-a)*y_gps;
        
    x_gps_mov_save(i) = x_gps_low;
    y_gps_mov_save(i) = y_gps_low;
    
    x_comp = b*x + (1-b)*x_gps_low;
    y_comp = b*y + (1-b)*y_gps_low;
    
    x_comp_save(i) = x_comp;
    y_comp_save(i) = y_comp;
    
    y_v_r(i) = v_r;
    y_v_l(i) = v_l;
    y_ideal_v_r(i) = ideal_v_r;
    y_ideal_v_l(i) = ideal_v_l;
    y_v(i) = v;
    y_w(i) = d_w;
    y_ideal_v(i) = ideal_v;
    y_ideal_w(i) = ideal_w;
    
end

figure(1)
%plot(x_save,y_save);
%hold on
plot(x_s_save,y_s_save,'r');
hold on
 plot(x_gps_save,y_gps_save,'k.');
 hold on
plot(x_kal_save,y_kal_save,'k');
hold on
% plot(x_comp_save,y_comp_save,'k');
% hold on
% plot(x_gps_mov_save,y_gps_mov_save,'g');
% hold on
    
figure(2)
subplot(3,1,1)
plot(x_v,x_save);
hold on
plot(x_v,x_s_save,'r');
hold on
% plot(x_v,x_gps_save,'k.');
% hold on
plot(x_v,x_kal_save,'k');
hold on

subplot(3,1,2)
plot(x_v,y_save);
hold on
plot(x_v,y_s_save,'r');
hold on
% plot(x_v,y_gps_save,'k.');
% hold on
plot(x_v,y_kal_save,'k');
hold on

subplot(3,1,3)
plot(x_v,theta_save);
hold on
plot(x_v,theta_s_save,'r');
hold on

% figure(3)
% subplot(4,1,1)
% plot(x_v,y_ideal_v_r);
% hold on
% plot(x_v,y_v_r,'k');
% hold on
% 
% subplot(4,1,2)
% plot(x_v,y_ideal_v_l);
% hold on
% plot(x_v,y_v_l,'k');
% hold on
% 
% subplot(4,1,3)
% plot(x_v,y_ideal_v);
% hold on
% plot(x_v,y_v,'k');
% hold on
% 
% subplot(4,1,4)
% plot(x_v,y_ideal_w);
% hold on
% plot(x_v,y_w,'k');
% hold on

%     
% for k = 0 : 1 : 0
%     
%     w_r = 20;
%     w_l = 20;
%     t = 10;
%     
%     n = 0:dt:t;
%     Nsamples = length(n);
%     Nsamples_2 = Nsamples;
%     
%     x_save = zeros(Nsamples,1);
%     y_save = zeros(Nsamples,1);
%     
% 
%     x_s_save = zeros(Nsamples,1);
%     y_s_save = zeros(Nsamples,1);
%     
%     x_gps_save = zeros(Nsamples,1);
%     y_gps_save = zeros(Nsamples,1);
%     
%     x_gps_mov_save = zeros(Nsamples,1);
%     y_gps_mov_save = zeros(Nsamples,1);
%     
%     for i = 1: 1: Nsamples
%         pre_theta = theta;
%         [d_x, d_y, d_w, ideal_v_r, ideal_v_l, ideal_v, ideal_w ] = Motion_Model_dt(w_r, w_l, wheel_radius ,width, theta,dt);
%         x = x + d_x*dt;
%         y = y + d_y*dt;
%         theta = theta + d_w*dt;
%         
%         x_save(i) = x;
%         y_save(i) = y;
%         theta_save(i) = theta;
%         
%         
%         [d_x, d_y, d_w, v_r, v_l, v] = Motion_Model_dt_error(w_r, w_l, wheel_radius ,width, theta_s,dt);
%         x_s = x_s + d_x*dt;
%         y_s = y_s + d_y*dt;
%         theta_s = theta_s + d_w*dt;
%         
%         
%         x_s_save(i) = x_s;
%         y_s_save(i) = y_s;
%         theta_s_save(i) = theta_s;
%         
%         x_gps = x_s + normrnd(0,vari_3);
%         y_gps = y_s + normrnd(0,vari_3);
%         
%         x_gps_save(i) = x_gps;
%         y_gps_save(i) = y_gps;
%         
%         A_13 = -ideal_v*dt*sin(pre_theta+ideal_w*dt/2);
%         A_23 = ideal_v*dt*cos(pre_theta+ideal_w*dt/2);
%         A = [1 0 A_13; 0 1 A_23; 0 0 1];
%         x_kal_p = A*x_kal;
%         Pp = A*P*A' + Q;
%         K = Pp*H'*inv(H*Pp*H'+R);
%         
%         z = [x_gps; y_gps];
%         x_kal = x_kal_p + K*(z - H*x_kal_p);
%         
%         P = Pp - K*H*Pp;
%         
%         x_kal_save(i) =x_kal(1);
%         y_kal_save(i) =x_kal(2);
%         
%         
%         x_gps_low = a*x_gps_low + (1-a)*x_gps;
%         y_gps_low = a*y_gps_low + (1-a)*y_gps;
%         
%         x_gps_mov_save(i) = x_gps_low;
%         y_gps_mov_save(i) = y_gps_low;
%         
%         y_v_r(i) = v_r;
%         y_v_l(i) = v_l;
%         y_ideal_v_r(i) = ideal_v_r;
%         y_ideal_v_l(i) = ideal_v_l;
%         y_v(i) = v;
%         y_w(i) = d_w;
%         y_ideal_v(i) = ideal_v;
%         y_ideal_w(i) = ideal_w;
%     end
%     
%     figure(1)
%     plot(x_save,y_save);
%     hold on
%     plot(x_s_save,y_s_save,'r');
%     hold on
% %     plot(x_gps_save,y_gps_save,'m.');
% %     hold on
%     plot(x_kal_save,y_kal_save,'k');
%     hold on
%     plot(x_gps_mov_save,y_gps_mov_save,'g');
%     hold on
%     
%     figure(2)
%     subplot(3,1,1)
%     plot(x_v+Nsamples_1 + Nsamples_2*2*k,x_save);
%     hold on
%     plot(x_v+Nsamples_1 + Nsamples_2*2*k,x_s_save,'r');
%     hold on
% %     plot(x_v+Nsamples_1 + Nsamples_2*2*k,x_gps_save,'m.');
% %     hold on
%     plot(x_v+Nsamples_1 + Nsamples_2*2*k,x_kal_save,'k');
%     hold on
%     
%     subplot(3,1,2)
%     plot(x_v+Nsamples_1 + Nsamples_2*2*k,y_save);
%     hold on
%     plot(x_v+Nsamples_1 + Nsamples_2*2*k,y_s_save,'r');
%     hold on
% %     plot(x_v+Nsamples_1 + Nsamples_2*2*k,y_gps_save,'m.');
% %     hold on
%     plot(x_v+Nsamples_1 + Nsamples_2*2*k,y_kal_save,'k');
%     hold on
%     
%     subplot(3,1,3)
%     plot(x_v+Nsamples_1 + Nsamples_2*2*k,theta_save);
%     hold on
%     plot(x_v+Nsamples_1 + Nsamples_2*2*k,theta_s_save,'r');
%     hold on
%     
%     figure(3)
%     subplot(4,1,1)
%     plot(x_v+Nsamples_1 + Nsamples_2*2*k,y_ideal_v_r);
%     hold on
%     plot(x_v+Nsamples_1 + Nsamples_2*2*k,y_v_r,'k');
%     hold on
%     
%     subplot(4,1,2)
%     plot(x_v+Nsamples_1 + Nsamples_2*2*k,y_ideal_v_l);
%     hold on
%     plot(x_v+Nsamples_1 + Nsamples_2*2*k,y_v_l,'k');
%     hold on
%     
%     subplot(4,1,3)
%     plot(x_v+Nsamples_1 + Nsamples_2*2*k,y_ideal_v);
%     hold on
%     plot(x_v+Nsamples_1 + Nsamples_2*2*k,y_v,'k');
%     hold on
%     
%     subplot(4,1,4)
%     plot(x_v+Nsamples_1 + Nsamples_2*2*k,y_ideal_w);
%     hold on
%     plot(x_v+Nsamples_1 + Nsamples_2*2*k,y_w,'k');
%     hold on
%     
%     x_save = zeros(Nsamples,1);
%     y_save = zeros(Nsamples,1);
%     
%     x_s_save = zeros(Nsamples,1);
%     y_s_save = zeros(Nsamples,1);
%     
%     x_gps_save = zeros(Nsamples,1);
%     y_gps_save = zeros(Nsamples,1);
%     
%     x_gps_mov_save = zeros(Nsamples,1);
%     y_gps_mov_save = zeros(Nsamples,1);
%     
%     w_r = 20;
%     w_l = 20;
%     
%     for i = 1: 1: Nsamples
%         pre_theta = theta;
%         [d_x, d_y, d_w, ideal_v_r, ideal_v_l, ideal_v, ideal_w] = Motion_Model_dt(w_r, w_l, wheel_radius ,width, theta,dt);
%         x = x + d_x*dt;
%         y = y + d_y*dt;
%         theta = theta + d_w*dt;
%         
%         x_save(i) = x;
%         y_save(i) = y;
%         theta_save(i) = theta;
%         
%         
%         [d_x, d_y, d_w, v_r, v_l, v] = Motion_Model_dt_error(w_r, w_l, wheel_radius ,width, theta_s,dt);
%         x_s = x_s + d_x*dt;
%         y_s = y_s + d_y*dt;
%         theta_s = theta_s + d_w*dt;
%         
%         x_s_save(i) = x_s;
%         y_s_save(i) = y_s;
%         theta_s_save(i) = theta_s;
%         
%         x_gps = x_s + normrnd(0,vari_3);
%         y_gps = y_s + normrnd(0,vari_3);
%         
%         x_gps_save(i) = x_gps;
%         y_gps_save(i) = y_gps;
%         
%         A_13 = -ideal_v*dt*sin(pre_theta+ideal_w*dt/2);
%         A_23 = ideal_v*dt*cos(pre_theta+ideal_w*dt/2);
%         A = [1 0 A_13; 0 1 A_23; 0 0 1];
%         x_kal_p = A*x_kal;
%         Pp = A*P*A' + Q;
%         K = Pp*H'*inv(H*Pp*H'+R);
%         
%         z = [x_gps; y_gps];
%         x_kal = x_kal_p + K*(z - H*x_kal_p);
%         
%         P = Pp - K*H*Pp;
%         
%         x_kal_save(i) =x_kal(1);
%         y_kal_save(i) =x_kal(2);
%         
%         x_gps_low = a*x_gps_low + (1-a)*x_gps;
%         y_gps_low = a*y_gps_low + (1-a)*y_gps;
%         
%         x_gps_mov_save(i) = x_gps_low;
%         
%         y_gps_mov_save(i) = y_gps_low;
%         
%         y_v_r(i) = v_r;
%         y_v_l(i) = v_l;
%         y_ideal_v_r(i) = ideal_v_r;
%         y_ideal_v_l(i) = ideal_v_l;
%         y_v(i) = v;
%         y_w(i) = d_w;
%         y_ideal_v(i) = ideal_v;
%         y_ideal_w(i) = ideal_w;
%     end
%     
%     figure(1)
%     plot(x_save,y_save);
%     hold on
%     plot(x_s_save,y_s_save,'r');
%     hold on
% %     plot(x_gps_save,y_gps_save,'m.');
% %     hold on
%     plot(x_kal_save,y_kal_save,'k');
%     hold on
%     plot(x_gps_mov_save,y_gps_mov_save,'g');
%     hold on
%     
%     figure(2)
%     subplot(3,1,1)
%     plot(x_v+Nsamples_1 + Nsamples_2 + Nsamples_2*2*k,x_save);
%     hold on
%     plot(x_v+Nsamples_1 + Nsamples_2 + Nsamples_2*2*k,x_s_save,'r');
%     hold on
% %     plot(x_v+Nsamples_1 + Nsamples_2 + Nsamples_2*2*k,x_gps_save,'m.');
% %     hold on
%     plot(x_v+Nsamples_1 + Nsamples_2 + Nsamples_2*2*k,x_kal_save,'k');
%     hold on
%     
%     subplot(3,1,2)
%     plot(x_v+Nsamples_1 + Nsamples_2 + Nsamples_2*2*k,y_save);
%     hold on
%     plot(x_v+Nsamples_1 + Nsamples_2 + Nsamples_2*2*k,y_s_save,'r');
%     hold on
% %     plot(x_v+Nsamples_1 + Nsamples_2 + Nsamples_2*2*k,y_gps_save,'m.');
% %     hold on
%      plot(x_v+Nsamples_1 + Nsamples_2 + Nsamples_2*2*k,y_kal_save,'k');
%     hold on
%     
%     subplot(3,1,3)
%     plot(x_v+Nsamples_1 + Nsamples_2 + Nsamples_2*2*k,theta_save);
%     hold on
%     plot(x_v+Nsamples_1 + Nsamples_2 + Nsamples_2*2*k,theta_s_save,'r');
%     hold on
% 
%     figure(3)
%     subplot(4,1,1)
%     plot(x_v+Nsamples_1 + Nsamples_2 + Nsamples_2*2*k ,y_ideal_v_r);
%     hold on
%     plot(x_v+Nsamples_1 + Nsamples_2 + Nsamples_2*2*k,y_v_r,'k');
%     hold on
%     
%     subplot(4,1,2)
%     plot(x_v+Nsamples_1 + Nsamples_2 + Nsamples_2*2*k,y_ideal_v_l);
%     hold on
%     plot(x_v+Nsamples_1 + Nsamples_2 + Nsamples_2*2*k,y_v_l,'k');
%     hold on
%     
%     subplot(4,1,3)
%     plot(x_v+Nsamples_1 + Nsamples_2 + Nsamples_2*2*k,y_ideal_v);
%     hold on
%     plot(x_v+Nsamples_1 + Nsamples_2 + Nsamples_2*2*k,y_v,'k');
%     hold on
%     
%     subplot(4,1,4)
%     plot(x_v+Nsamples_1 + Nsamples_2 + Nsamples_2*2*k,y_ideal_w);
%     hold on
%     plot(x_v+Nsamples_1 + Nsamples_2 + Nsamples_2*2*k,y_w,'k');
%     hold on
% end

figure(1)
grid on
xlabel('x displacement(m)')
ylabel('y displacement(m)')
legend('real', 'GPS','kalman')

figure(2)
subplot(3,1,1)
xlabel('time (0.1s)')
ylabel('x displacement(m)')
legend('ideal','real', 'kalman')

subplot(3,1,2)
xlabel('time (0.1s)')
ylabel('y displacement(m)')
legend('ideal','real', 'kalman')

subplot(3,1,3)
xlabel('time (0.1s)')
ylabel('theta displacement(rad)')
legend('ideal','real')

% figure(3)
% subplot(4,1,1)
% xlabel('time (0.1s)')
% ylabel('right wheel velocity(m/s)')
% legend('ideal','real')
% 
% subplot(4,1,2)
% xlabel('time (0.1s)')
% ylabel('left wheel velocity(m/s)')
% legend('ideal','real')
% 
% subplot(4,1,3)
% xlabel('time (0.1s)')
% ylabel('car velocity(m/s)')
% legend('ideal','real')
% 
% subplot(4,1,4)
% xlabel('time (0.1s)')
% ylabel('car angular velocity(m/s)')
% legend('ideal','real')
% 


