function [Hradj_X] = hankelize_adj_fft(X, N, F, H, G)
%HANKELIZE_ADJ_FFT 利用fft实现hankel矩阵的伴随算子
%   X为rn*m维矩阵,r为hankel矩阵行块数,n为基本向量维度,N=n+m-1为总可用数据量
Hradj_X = real(1/N*F'*diag(H*X*G'));
end

