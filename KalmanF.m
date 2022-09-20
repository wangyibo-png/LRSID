function [K, e, y, SNR] = KalmanF(A, B, C, D, Q, R, u, w, v, n, p, N)
%KALMANF 此处显示有关此函数的摘要
%   此处显示详细说明
sys = ss(A, [B eye(n)], C, D, 1);
[~, K] = kalman(sys, Q, R, [0;0]);
x = zeros(n, 1);
y = zeros(p, N);
for i = 1:N
    y(:, i) = C*x+D*u(:, i)+v(:, i);
    x = A*x+B*u(:, i)+w(:, i);
end
e = zeros(p, N);
x = zeros(n, 1);
for i = 1:N
    haty = C*x+D*u(:, i);
    e(:, i) = y(:, i)-haty;
    x = A*x+B*u(:, i)+K*e(:, i);
end
SNR = 10*log10(var(y)/var(e));
end

