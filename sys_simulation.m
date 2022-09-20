function [y_out] = sys_simulation(u, A, B, C, D, K, n, p, E)
%SYS_SIMULATION 此处显示有关此函数的摘要
%   此处显示详细说明
N = size(u, 2);
y_out = zeros(p, N);
x = zeros(n, 1);
for i = 1:N
    y_out(:, i) = C*x+D*u(:, i)+E(:, i);
    x = A*x+B*u(:, i)+K*E(:, i);
end
end