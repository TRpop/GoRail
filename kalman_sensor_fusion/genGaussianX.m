function X = genGaussianX(T,AC)
% Generating a stationary Gaussian process, X of length T with
% auto-covariances
% ========================================================================
% [1]Wood & Chan (1994) Simulation of stationary Gaussian processes in [0,1]d
% [2]Davies (2001) Simulation of a stationary Gaussian time-series
% [3]Craigmile (2003) Simulating a class of stationary Gaussian processes using the Davies?Harte algorithm, with application to long memory processes
% ========================================================================
% [Input]
% T : length of time-series
% AC : vector of auto-covarainces with lag 0, 1, 2 ....
% ========================================================================

Lac = length(AC);
C = [];     % (2T-2)

for j = 1:(T-1)
   if j <= Lac 
       C(j) = AC(j);
   else
       C(j) = 0;
   end
end
for j = 2:T
    if j <= Lac
       C(2*T-j) = AC(j); 
    else
       C(2*T-j) = 0;  
    end
end
C = C';

G2 = [];    % (2T-2)
G2 = sqrt(fft(C));

Z = [];     % (2T-2)
Z(1) = normrnd(0,2);
Z(T) = normrnd(0,2);
for j = 1:T-2
    Z1 = normrnd(0,1);
    Z2 = normrnd(0,1);
    Z(1+j) = Z1 + i*Z2; 
    Z(2*T-1-j) = Z1 - i*Z2;
end
Z = Z';

ZG2 = [];   % (2T-2)
ZG2 = Z .* G2;

XX=[];      % (2T-2)
XX = ifft(ZG2);
X = (XX(1:T) - i*imag(XX(1:T))) * (T-1)^0.5 ;     % (T)
end

