%% poly_mult_mod_ntt: 用 NTT/INTT 计算 (f⋆g mod (X^N+1)) mod q，
% 其中 f, g 为长度 N 的多项式系数向量，q 为当前模数，
% psi 为该模数下的 2N 次单位根。
function c = ckks_poly_mult_mod(f, g, q, N, psi)
    % 1. NTT(f), NTT(g)
    F = ntt(f, q, N, psi);
    G = ntt(g, q, N, psi);

    % 2. 频域逐点相乘 H = F .* G (mod q)
    H = modmul(F, G, q);

    % 3. INTT(H)
    c = intt(H, q, N, psi);
    % c 的系数已在 [0, q-1]
end