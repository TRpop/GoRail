function [v_x, v_y, dx, dy, w, new_theta] = ideal_model( w_r, w_l, wheel_r ,width, theta,dt)

    v_r = wheel_r*w_r;
    v_l = wheel_r*w_l;
    
    v = (v_r + v_l)/2.0;
    w = (v_r - v_l)/(width/2.0);
    
    new_theta = theta + w*dt;
    
    v_x = v*cos(new_theta);
    v_y = v*sin(new_theta);
    
    dx = v_x*dt;
    dy = v_y*dt;
    
end