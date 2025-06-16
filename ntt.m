%% ntt.m
% 计算 f 在 (ℤ_q[X]/(X^N+1)) 下的 negacyclic NTT 变换
% 输入：
%   f        — 1×N 整数向量，系数 ∈ [0, q-1]
%   q        — 素数模数，满足 q ≡ 1 mod (2N)
%   N        — 多项式阶，2 的幂
%   psi      — primitive 2N-th root mod q，使得 psi^N ≡ -1 mod q
% 输出：
%   F        — 1×N，在频域的 NTT 结果，满足 negacyclic 定义

function F = ntt(f, q, N, psi)
    % 参数检查：N 必须是 2 的幂
    assert(mod(log2(N),1) == 0, 'ntt_def: N must be a power of two.');

    % 1) 预计算 ψ^j (j=0..N-1)，然后 pre-multiply: f_tilde[j] = f[j] * ψ^j mod q
    psi_pows = zeros(1, N);
    psi_pows(1) = 1;
    for j = 2:N
        psi_pows(j) = modmul(psi_pows(j-1), psi, q);
    end
    f_tilde = modmul(f, psi_pows, q);  % 1×N 向量

    % 2) 对 f_tilde 做 in-place 位反转排列（Bit-reversal），为迭代 NTT 做准备
    bits = log2(N);
    for i = 0:(N-2)
        rev = bitrev(i, bits);
        if rev > i
            % MATLAB 下标是 1-based
            tmp = f_tilde(i+1);
            f_tilde(i+1) = f_tilde(rev+1);
            f_tilde(rev+1) = tmp;
        end
    end

    % 3) 迭代式 Cooley–Tuk radix‑2 NTT，使用 ω = ψ^2 mod q
    omega = modmul(psi, psi, q);
    len = 1;
    while len < N
        m = 2 * len;                         % 当前蝶形块大小
        w_m = powmod(omega, N / m, q);       % 本阶段的主旋转因子 ω_m = ω^(N/m)
        for start = 0 : m : (N-1)
            w = 1;
            for offset = 0:(len-1)
                u = f_tilde(start + offset + 1);            % 1-based
                v = modmul(w, f_tilde(start + offset + len + 1), q);
                % 蝶形（不含缩放）
                f_tilde(start + offset + 1)       = mod(u + v, q);
                f_tilde(start + offset + len + 1) = mod(u - v, q);
                w = modmul(w, w_m, q);
            end
        end
        len = m;
    end

    % 4) 返回 N 点 negacyclic NTT 结果 F(k) = ∑ f[n]*ψ^{n(2k+1)} mod q
    F = f_tilde;
end