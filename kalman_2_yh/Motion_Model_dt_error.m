function [ d_x, d_y, d_w, v_r,v_l,v] = Motion_Model_dt_error( w_r, w_l, wheel_r ,width, theta,dt)
% w_r = right wheel speed
% w_l = left wheel speed
% width = between right, left wheel width
% wheel_r = wheel radius  
% theta = angle of Model

v_r =  wheel_r * w_r;
v_l =  wheel_r * w_l;

e_r = normrnd(0,0.02);
while(e_r > 0)
   e_r = normrnd(0,0.02);
end

e_l = normrnd(0,0.02);
while(e_l > 0)
   e_l = normrnd(0,0.02);
end

v_r = v_r + v_r*e_r;
v_l = v_l + v_l*e_l;

v = (v_r + v_l)/2;
w = (v_r - v_l)/width;

e_w = w*normrnd(0,0.00001);

w = w+e_w;

d_x = cos(theta+(w+e_w)*dt/2)*(v);% + sin(theta)*e_T; % d_x = cos(theta)*(v + e);
d_y = sin(theta+(w+e_w)*dt/2)*(v);% + cos(theta)*e_T; % d_y = sin(theta)*(v + e);
d_w = w + e_w;
%d_w = d_w + normrnd(0,0.1);

end
