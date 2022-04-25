clear all
close all
clc
format long; 

k = 21;
j =(-1)^0.5;

d_f = 0.0001;
num = 2*k+1;
p = zeros(num, 1);
n = (0:num-1) / num;

F = 0:d_f:1;
leng_F = length(F);
Hd = zeros(leng_F, 1);
dn = leng_F / num;

t_1 = zeros(num, 1);
t = zeros(num, 1);

sample_int = round((1/num)/d_f);
position = zeros(num, 1);


for i = 1:(1/d_f)/2
    Hd(i, 1) = -j;
end

for i = ((1/d_f)/2+1):(1/d_f)
    Hd(i, 1) = j;
end

Hd(1, 1) =0; 


for i = 1:num
    p(i ,1) = Hd(sample_int*(i-1)+1, 1);
    position(i, 1) = sample_int*d_f*(i-1);
end


t_1 = ifft(p);
t = transpose(t_1);
t = [t(k+1:num), t(1:k)]; 
T = fft(t_1);
imR = spline(n,imag(T), F);

figure;
plot(0:d_f:1, imag(Hd), 'b:');
hold on;
plot(position, imag(p), 'bd');
plot(0:d_f:1, imR, 'r');
hold off;
title(['Frequency Response']);
xlabel('Normalized Frequency (Hz)');
legend('Hd(F)', 'Sample Point', 'T(F)');

figure; 
stem(n, t);
title(['Impulse Response']);

