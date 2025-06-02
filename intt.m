%% intt.m
% 计算频域向量 H 的逆 negacyclic NTT，得到时域多项式系数 (mod X^N+1) mod q
% 输入：
%   H        — 1×N 频域数据，系数 ∈ [0, q-1]
%   q, N     — 如上
%   psi      — primitive 2N-th root mod q
% 输出：
%   h        — 1×N，INTT 后得到的时域系数，范围 [0, q-1]

function h = intt(H, q, N, psi)
    % 参数检查：N 必须是 2 的幂
    assert(mod(log2(N),1) == 0, 'intt_def: N must be a power of two.');

    % 1) 先复制输入以免覆盖：Ht = H(1..N)
    Ht = H(:).';  % 保证是行向量

    % 2) 对 Ht 做位反转排列（Bit-reversal）
    bits = log2(N);
    for i = 0:(N-2)
        rev = bitrev(i, bits);
        if rev > i
            tmp = Ht(i+1);
            Ht(i+1) = Ht(rev+1);
            Ht(rev+1) = tmp;
        end
    end

    % 3) 迭代式 Cooley–Tuk 逆 NTT (radix‑2)
    %    使用 ω_inv = (ψ^2)^{-1} mod q
    omega = modmul(psi, psi, q);
    omega_inv = powmod(omega, q-2, q);  % 因 q 是素数，ω^{-1} = ω^{q-2} mod q
    len = 2;
    while len <= N
        w_m_inv = powmod(omega_inv, N / len, q);  % 本阶段主旋转因子 ω_m^{-1}
        half = len / 2;
        for start = 0 : len : (N-1)
            w = 1;
            for offset = 0:(half-1)
                u = Ht(start + offset + 1);                % 1-based
                v = modmul(Ht(start + offset + half + 1), w, q);
                % 正确的逆蝶形（先将右半边乘 w，再做加减）
                Ht(start + offset + 1)        = mod(u + v, q);
                Ht(start + offset + half + 1) = mod(u - v, q);
                w = modmul(w, w_m_inv, q);
            end
        end
        len = len * 2;
    end

    % 4) 整体除以 N（乘以 N^{-1} mod q）
    inv_N = modinv(N, q);
    for i = 1:N
        Ht(i) = modmul(Ht(i), inv_N, q);
    end

    % 5) 最后一步：对每个位置 j (0..N-1) 再乘以 ψ^{-j}
    %    因为正变换时做了 f_tilde[j] = f[j]*ψ^j，这里要消去它 ⇒ 乘 ψ^{-j}
    psi_inv = powmod(psi, q-2, q);  % ψ^{-1} mod q
    psinv_pows = zeros(1, N);
    psinv_pows(1) = 1;
    for j = 2:N
        psinv_pows(j) = modmul(psinv_pows(j-1), psi_inv, q);
    end

    % 6) 组合输出
    h = modmul(Ht, psinv_pows, q);
end