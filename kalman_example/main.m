clear all;

dt = 0.2;
t = 0:dt:100;

Nsamples = length(t);

Xsaved = zeros(Nsamples, 1);
Zsaved = zeros(Nsamples, 1);

for k = 1:Nsamples
    z = 14 + 4*randn(1,1);
    volt = kalman_example(z);
    
    Xsaved(k) = volt;
    Zsaved(k) = z;
end

plot(t, Xsaved);
hold on
plot(t, Zsaved);
hold off