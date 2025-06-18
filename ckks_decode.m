function z_decoded = ckks_decode(m_int, Delta)
% 使用 CKKS 专用的 FFT
N = length(m_int);
g = 5;
% 1. CKKS 专用的FFT
v_scaled = ckks_fft(m_int);

% 2. 缩放
v = v_scaled / Delta;

% 3. 提取明文
z_decoded = zeros(1, N/2);
for i = 1:N/2
    z_decoded(i) = v((powmod(g,i,2*N)+1)/2);
end
end

%% 以下版本支持除旋转以外同态操作(换言之不能旋转）
% function z_decoded = ckks_decode(m_int, Delta)
% % 使用 CKKS 专用的 FFT
% N = length(m_int);
% 
% % 1. CKKS 专用的FFT
% v_scaled = ckks_fft(m_int);
% 
% % 2. 缩放
% v = v_scaled / Delta;
% 
% % 3. 提取明文
% z_decoded = reshape(v(1:N/2), [1, N/2]);
% end