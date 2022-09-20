function [Tradj_X] = toplitz_adj_fft(X, N, F, H, G)
%TOPLITZ_ADJ_FFT 此处显示有关此函数的摘要
%   此处显示详细说明
Tradj_X = real(1/N*F'*diag(H*X*G'));
end

