function poly_q = ckks_moddown(poly_pq, params)
%CKKS_MODDOWN  RNS-CKKS 的 ModDown 基还原（修正版）
%  poly_pq: cell, 长度 = k + M
%    {1..k}       是 p_i 分量 b^(p_i)
%    {k+1..k+M} 是 q_j 分量 b^(q_j)
%  返回 poly_q: cell, 长度 = l+1，对应新的 b'(q_0..q_l)

    % 读取参数
    P_primes         = params.P_primes;        % 1×k
    Q_primes         = params.Q_primes;        % 1×(L+1)
    hatpi_mod_qj     = params.hatpi_mod_qj;    % k×(L+1)
    inv_hatpi_mod_pi = params.inv_hatpi_mod_pi;% 1×k
    P_inv_mod_qj     = params.P_inv_mod_qj;    % 1×(L+1)

    k = params.k;
    M = numel(poly_pq) - k;

    poly_q = cell(1, M);
    N = params.N;

    for j = 1:M
        qj   = Q_primes(j);
        bkj  = poly_pq{k + j};         % 原始 b^(q_j)

        % === Step 1: Conv_B→C，计算 a(j) ===
        a_j = zeros(1, N);
        for i = 1:k
            p_i   = P_primes(i);
            % 1) v_i = b^(p_i) * inv_hatpi_mod_pi(i)   (mod p_i)
            v_i = modmul(poly_pq{i}, inv_hatpi_mod_pi(i), p_i);
            % 把 v_i 从 p_i 域提升到 q_j
            % 2) v_i_lift = v_i * hatpi_mod_qj(i,j)    (mod q_j)
            v_i_lift = modmul(mod(v_i, qj), hatpi_mod_qj(i, j), qj);
            % 3) 累加
            a_j = modadd(a_j, v_i_lift, qj);
        end

        % === Step 2: b'(j) = P^{-1} * (b(k+j) - a_j) mod q_j ===
        % 先做差
        diff = mod( bkj - a_j, qj );
        % 再乘以 P_inv_mod_qj(j)
        poly_q{j} = modmul(diff, P_inv_mod_qj(j), qj);
    end
end