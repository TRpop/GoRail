function [x_new, y_new] = getGPS(x, y)
    
    x_new = x + normrnd(0, 4);
    y_new = y + normrnd(0, 4);
end