function [M] = matrix_HrHr(N, F, H, G)
%MATRIX_HRHR 产生辅助矩阵M
%   x'*M*x = <Hr(x), Hr(x)>; Hr*(Hr(x)) = Mx
M = real(1/N^2*F'*((H*H').*((G*G').'))*F);
end

