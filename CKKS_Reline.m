function    ct_out = CKKS_Reline(params, keys, ct)
    M = numel(ct.d0);
    Q = params.Q_primes;
    P = params.P_primes;
    N = params.N;
    k = params.k;
    
    % 1. 重线性化：对 d2 使用 SWK
    if ~isKey(keys.map, 0)
        error('CKKS_Keyswitch: 未找到与 Relinearation 对应的 KeySwitchKey。');
    end
    idx_ks = keys.map(0);
    swk    = keys.swk{idx_ks};

    ct_modup_d2 = ckks_modup(ct.d2, params);
    
    % 预先分配 ct_ks.c0, ct_ks.c1（长度 = k+M）
    ct_ks.c0 = cell(1, numel(ct_modup_d2));
    ct_ks.c1 = cell(1, numel(ct_modup_d2));

    for i = 1:k
        p_i   = P(i);
        a_pi  = swk.a{i};  % 长度-N 多项式，系数 mod p_i
        b_pi  = swk.b{i};

        % 3.1 c0^{(p_i)} = (ct_rot1^{(p_i)} ⊛ b_pi) mod (X^N+1, p_i)
        ct_ks.c0{i} = ckks_poly_mult_mod(ct_modup_d2{i}, b_pi, p_i, N, params.psi{i});

        % 3.2 c1^{(p_i)} = (ct_rot1^{(p_i)} ⊛ a_pi) mod (X^N+1, p_i)
        ct_ks.c1{i} = ckks_poly_mult_mod(ct_modup_d2{i}, a_pi, p_i, N, params.psi{i});
    end

    %% 1.4. 在 {q_j} 下做多项式乘法
    for j = 1:M
        idx_q = k + j;        % MATLAB 索引 = k + j
        q_j   = Q(j);       % q_{j-1}
        a_qj  = swk.a{idx_q};      % 长度-N 多项式
        b_qj  = swk.b{idx_q};

        % 4.1 c0^{(q_j)} = (ct_rot1^{(q_j)} ⊛ b_qj) mod (X^N+1, q_j)
        ct_ks.c0{idx_q} = ckks_poly_mult_mod(ct_modup_d2{idx_q}, b_qj, q_j, N, params.psi{idx_q});

        % 4.2 c1^{(q_j)} = (ct_rot1^{(q_j)} ⊛ a_qj) mod (X^N+1, q_j)
        ct_ks.c1{idx_q} = ckks_poly_mult_mod(ct_modup_d2{idx_q}, a_qj, q_j, N, params.psi{idx_q});
    end

    %% 1.5. ModDown: 将 ct_ks.c0 和 ct_ks.c1 从 {p,q} 折回到 {q₀..q_{M-1}}
    c0_q = ckks_moddown(ct_ks.c0, params);  % 返回长度 = M 的 cell
    c1_q = ckks_moddown(ct_ks.c1, params);

    % 2. 按照 c0' = Δc0 + d0, c1' = Δc1 + d1 合并
    ct_out.c0 = cell(1, M);
    ct_out.c1 = cell(1, M);
    for j = 1:M
        qj = Q(j);
        % Δc0 + d0
        ct_out.c0{j} = modadd(c0_q{j}, ct.d0{j}, qj);
        % Δc1 + d1
        ct_out.c1{j} = modadd(c1_q{j}, ct.d1{j}, qj);
    end
end