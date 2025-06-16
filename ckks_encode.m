function m_int = ckks_encode(z_slots, Delta)
% 使用 CKKS 专用的 IFFT
N = length(z_slots) * 2;

% 1. 共轭对称扩展，包含为支持旋转而进行的重排序
v = zeros(N, 1);
g = 5;
for i = 1 : N/2
    v((powmod(g,i,2*N)+1)/2) = z_slots(i);
    v(N + 1 - (powmod(g,i,2*N)+1)/2) = conj(z_slots(i));
end

% 让我们拿N=8举个例子。
% 使用下面被注释的ckks_encode时，输入[x0 x1 x2 x3]被扩展成[x0 x1 x2 x3 x3_conj x2_conj x1_conj x0_conj]，
% 进一步被填充上0后被执行2N点ifft，即扩展后的v的内容与旋转因子有如下对应关系：
% [x0 x1 x2 x3 x3_conj x2_conj x1_conj x0_conj] <=>omega^[1 3 5 7 9 11 13 15];

% 我们重排序的预期目标是使得扩展后的v的内容与旋转因子有如下对应关系：
% [x0 x1 x2 x3 x3_conj x2_conj x1_conj x0_conj] <=>omega^[5 9 13 1 ? ? ? ?];
% 其中5 = g^1 mod 2N, 9 = g^2 mod 2N, 13 = g^3 mod 2N, 1 = g^4 mod 2N。
% 考虑到扩展的对称性，扩展后的v的内容与旋转因子应当有如下对应关系
% [x0 x1 x2 x3 x3_conj x2_conj x1_conj x0_conj] <=>omega^[5 9 13 1 16-1 16-13 16-9 16-5];
% 进一步考虑由于调用ifft时无法设置旋转因子的位置，
% 我们只能对v = [x0 x1 x2 x3 x3_conj x2_conj x1_conj x0_conj]调换内部顺序，
% 即调整为[x3 x2_conj x0 x1_conj x1 x0_conj x2 x3_conj]，使其对应的旋转因子仍然是[1 3 5 7 9 11 13 15]。

% 2. 缩放
v_scaled = Delta * v;

% 3. CKKS 专用的逆FFT
m_float = reshape(ckks_ifft(v_scaled), [1, N]);

% 4. 取整并验证实数性
assert(max(abs(imag(m_float))) < 1e-10, 'IFFT 结果非实数！');
m_int = round(real(m_float));
end

%% 以下版本支持除旋转以外同态操作(换言之不能旋转）
% function m_int = ckks_encode(z_slots, Delta)
% % 使用 CKKS 专用的 IFFT
% N = length(z_slots) * 2;
% 
% % 1. 共轭对称扩展
% v = zeros(N, 1);
% v(1:N/2) = z_slots(:);
% v(N/2+1:end) = conj(z_slots(end:-1:1));
% 
% % 2. 缩放
% v_scaled = Delta * v;
% 
% % 3. CKKS 专用的逆FFT
% m_float = reshape(ckks_ifft(v_scaled), [1, N]);
% 
% % 4. 取整并验证实数性
% assert(max(abs(imag(m_float))) < 1e-10, 'IFFT 结果非实数！');
% m_int = round(real(m_float));
% end