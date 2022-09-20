function [Tr_x] = toplitz_fft(x, n, m, N, F, H, G)
%TOPLITZ_FFT 利用fft生成toplitz矩阵
%   x为(N*n*m)维向量,x中的元素为n*m维矩阵,x为列向量
Tr_x = real(1/N*H'*diag(F*x)*G);
% Tr_x(abs(Tr_x) < 1e-6) = 0;
end

