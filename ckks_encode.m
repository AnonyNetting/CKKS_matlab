function m_int = ckks_encode(z_slots, Delta)
% 使用 CKKS 专用的 IFFT
N = length(z_slots) * 2;

% 1. 共轭对称扩展
v = zeros(N, 1);
v(1:N/2) = z_slots(:);
v(N/2+1:end) = conj(z_slots(end:-1:1));

% 2. 缩放
v_scaled = Delta * v;

% 3. CKKS 专用的逆FFT
m_float = reshape(ckks_ifft(v_scaled), [1, N]);

% 4. 取整并验证实数性
assert(max(abs(imag(m_float))) < 1e-10, 'IFFT 结果非实数！');
m_int = round(real(m_float));
end