function z_decoded = ckks_decode(m_int, Delta)
% 使用 CKKS 专用的 FFT
N = length(m_int);

% 1. CKKS 专用的FFT
v_scaled = ckks_fft(m_int);

% 2. 缩放
v = v_scaled / Delta;

% 3. 提取明文
z_decoded = reshape(v(1:N/2), [1, N/2]);
end