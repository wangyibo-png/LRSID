function [Hr_x] = hankelize_fft(x, n, m, N, H, F, G)
%HANKELIZE_FFT 利用fft生成hankel矩阵
%   H,F,G分别为n维向量情况下的fft矩阵,x=[x(1);...;x(N)],x(i)为n*m为矩阵
Hr_x = real(1/N*H'*diag(F*vec(x))*G);
end