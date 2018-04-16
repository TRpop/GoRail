function [ d_x, d_y, d_w, v_r, v_l, v, w] = Motion_Model_dt( w_r, w_l, wheel_r ,width, theta,dt)
% w_r = right wheel speed
% w_l = left wheel speed
% width = between right, left wheel width
% wheel_r = wheel radius  
% theta = angle of Model

v_r =  wheel_r *w_r;
v_l =  wheel_r *w_l;

v = (v_r + v_l)/2;
w = (v_r - v_l)/width;

d_x = cos(theta+w*dt/2.0)*v; % d_x = cos(theta)*(v + e);
d_y = sin(theta+w*dt/2.0)*v; % d_y = sin(theta)*(v + e);
d_w = w;

end
