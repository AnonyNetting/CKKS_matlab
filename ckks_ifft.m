function x = ckks_ifft(X)
% CKKS 专用IFFT（基于2N-th单位根）
% 输入: X 是频域信号（基于e^(2πi/2N)）
% 输出: x 是时域信号

N = length(X);
% 1. 构造长度 2N 的复数向量 Z
Z = zeros(1, 2*N);               % 默认零（偶数位置保持 0）

% 将 X 放到 Z 的奇数位置 (0-based 偶数 MATLAB 索引)
for j = 0:(N-1)
    idx = 2*j + 2;               % 这是 MATLAB 的 1-based 下标：m+1 = 2*j+2
    Z(idx) = X(j+1);
end

% 2. 对 Z 做一次标准 2N 点 IFFT
F = ifft(Z);                     % MATLAB 默认使用 exp(2πi/(2N)) 等价单位根

% 3. 取前 N 项并乘 2
x = reshape(2 * F(1:N),[N, 1]);
end