function [Hr_x] = hankelize(x, n, N, r)
%HANKELIZE 基于x构成hankel矩阵
%   r为hankel矩阵的行块数,x=[x(1);x(2);...;x(N)],x(i)为n维列向量
X = mat2cell(x, ones(N, 1)*n, 1);
index = hankel(1:r, r:N);
Hr_x = cell2mat(X(index));
end