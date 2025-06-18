function ct_out = CKKS_Keyswitch(params, keys, ct_rot, r)
% CKKS_Keyswitch  将“旋转后私钥”或任何旧私钥下的密文切换到原始私钥下
%
%   ct_out = CKKS_Keyswitch(params, keys, ct_rot, r)
%
% 说明：
%   本函数在 RNS‐CKKS 框架下，将在“旧私钥”（如旋转后的 s_rot）下加密的
%   密文 ct_rot 切换到“原始私钥” s 下。切换过程基于 Type II KeySwitch：
%     1. ModUp：将 ct_rot.c1 从 {q₀..q_{M-1}} 扩展到 {p₀..p_{k-1}} ∪ {q₀..q_{M-1}}。
%     2. 在 {p_i} 下，用 (a^{(p_i)}, b^{(p_i)}) 做多项式乘法 + 乘 inv_hatpi_mod_pi。
%     3. 在 {q_j} 下，用 (a^{(q_j)}, b^{(q_j)}) 做多项式乘法 + 乘 P^{-1} mod q_j。
%     4. ModDown：将结果从 {p,q} 折回到 {q₀..q_{M-1}}。
%     5. 合并原始 c0： c0' = c0_old + Δc0 (mod q_j)； c1' = Δc1 (mod q_j)。
%     简单来说:KeySwitch(ct_rot.c1, r) => (ct_ks.c0, ct_ks.c1),
%             ct_out = (ct_rot.c0, 0) + (ct_ks.c0, ct_ks.c1).
%
% 接口：
%   - params : CKKS_Setup 返回的参数结构体，需包含：
%       params.N                : 多项式环维度 N
%       params.k                : 特殊基素数个数 k
%       params.P_primes         : 长度 k 的行向量 [p₀,…,p_{k-1}]
%       params.Q_primes         : 长度 (L+1) 的行向量 [q₀,…,q_L]
%       params.hatpi_mod_qj     : k×(L+1) 矩阵，(P/p_i) mod q_j
%       params.inv_hatpi_mod_pi : 1×k 向量，(P/p_i)^{-1} mod p_i
%       params.P_inv_mod_qj     : 1×(L+1) 向量，P^{-1} mod q_j
%       params.q_lj_mod_pi      : (L+1)×(L+1)×k 张量，(Q_ℓ/q_j) mod p_i
%       params.inv_q_lj_mod_qj  : (L+1)×(L+1) 矩阵，(Q_ℓ/q_j)^{-1} mod q_j
%       params.psi              : 1×(k+L+1) cell, 每个元素为该模下的 NTT 原根
%   - keys   : 由 KSGen 生成的密钥结构，包含：
%       keys.map    : containers.Map, 键=旋转步长 r, 值=索引 idx
%       keys.swk   : cell 数组, swk{idx} = struct('power',r,'a',{a_all},'b',{b_all})
%           a_all, b_all: 1×(k+L+1) cell, 前 k 是 p_i 下的 (a',b')，后 L+1 是 q_j 下的 (a',b')
%   - ct_rot : 输入密文，仅在 M = numel(ct_rot.c0) 个近似基 {q₀..q_{M-1}} 下有分量
%       ct_rot.c0 : 1×M cell, 每个元素是长度-N “c0” 多项式系数 (mod q_{j-1})
%       ct_rot.c1 : 1×M cell, “c1” 多项式系数 (mod q_{j-1})
%   - r      : 旋转步长，用于从 keys.map 中查索引
%
% 输出：
%   ct_out : 切换后密文，仅保留 M 个近似基 {q₀..q_{M-1}} 分量
%       ct_out.c0 : 1×M cell, 新的 c0' mod q_{j-1}
%       ct_out.c1 : 1×M cell, 新的 c1' mod q_{j-1}
% 算法：
%   1.  被旋转结果调用时(0 < r < N)：
%       传入参数ct_rot包含两部分(ct_rot.c0, ct_rot.c1)
%       ct_out = (ct_rot.c0, 0) + KeySwitch(ct_rot.c1, swk)

    %% 1. 在 keys 中查找与 r 对应的 KeySwitchKey
    if r <= 0 || r >= params.N
        erorr('illegal rotation step r: r must be in the range (0, N)');
    end
    if ~isKey(keys.map, r)
        error('CKKS_Keyswitch: 未找到与 r 对应的 KeySwitchKey。');
    end
    idx_ks = keys.map(r);
    swk    = keys.swk{idx_ks};
    % swk.a, swk.b 都是长度 (k + L + 1) 的 cell 数组

    %% 2. ModUp: 把 ct_rot.c1 从 {q₀..q_{M-1}} 扩展到 {p₀..p_{k-1}} ∪ {q₀..q_{M-1}}
    % 注意：这里传入的是 “ct_rot.c1”，长度 = M
    ct_rot1 = ckks_modup(ct_rot.c1, params);  
    % ct_rot1 是长度 k+M 的 cell, 每个 cell{i} 是长度-N 多项式系数
    % 其中 i=1..k 对应 p₀..p_{k-1}； i=k+1..k+M 对应 q₀..q_{M-1}

    %% 3. 在 {p_i} 下做多项式乘法
    P_primes         = params.P_primes;         % 1×k
    k = params.k;
    N = params.N;

    % 预先分配 ct_ks.c0, ct_ks.c1（长度 = k+M）
    ct_ks.c0 = cell(1, k + numel(ct_rot.c0));
    ct_ks.c1 = cell(1, k + numel(ct_rot.c0));

    for i = 1:k
        p_i   = P_primes(i);
        a_pi  = swk.a{i};  % 长度-N 多项式，系数 mod p_i
        b_pi  = swk.b{i};

        % 3.1 c0^{(p_i)} = (ct_rot1^{(p_i)} ⊛ b_pi) mod (X^N+1, p_i)
        ct_ks.c0{i} = ckks_poly_mult_mod(ct_rot1{i}, b_pi, p_i, N, params.psi{i});

        % 3.2 c1^{(p_i)} = (ct_rot1^{(p_i)} ⊛ a_pi) mod (X^N+1, p_i)
        ct_ks.c1{i} = ckks_poly_mult_mod(ct_rot1{i}, a_pi, p_i, N, params.psi{i});
    end

    %% 4. 在 {q_j} 下做多项式乘法
    Q_primes     = params.Q_primes;      % 1×(L+1)
    offset = k;                          % ct_rot1{k+1..k+M} 对应 q₀..q_{M-1}
    M      = numel(ct_rot.c0);

    for j = 1:M
        idx_q = offset + j;        % MATLAB 索引 = k + j
        q_j   = Q_primes(j);       % q_{j-1}
        a_qj  = swk.a{idx_q};      % 长度-N 多项式
        b_qj  = swk.b{idx_q};

        % 4.1 c0^{(q_j)} = (ct_rot1^{(q_j)} ⊛ b_qj) mod (X^N+1, q_j)
        ct_ks.c0{idx_q} = ckks_poly_mult_mod(ct_rot1{idx_q}, b_qj, q_j, N, params.psi{idx_q});

        % 4.2 c1^{(q_j)} = (ct_rot1^{(q_j)} ⊛ a_qj) mod (X^N+1, q_j)
        ct_ks.c1{idx_q} = ckks_poly_mult_mod(ct_rot1{idx_q}, a_qj, q_j, N, params.psi{idx_q});
    end

    %% 5. ModDown: 将 ct_ks.c0 和 ct_ks.c1 从 {p,q} 折回到 {q₀..q_{M-1}}
    c0_q = ckks_moddown(ct_ks.c0, params);  % 返回长度 = M 的 cell
    c1_q = ckks_moddown(ct_ks.c1, params);

    %% 6. 合并原始 c0，构造最终 ct_out
    ct_out.c0 = cell(1, M);
    ct_out.c1 = cell(1, M); 
    for j = 1:M
        qj = Q_primes(j);  
        % 新的 c0_j = 原始 c0_j + delta c0_j  mod q_j
        ct_out.c0{j} = modadd(ct_rot.c0{j}, c0_q{j}, qj);
        % 新的 c1_j = delta c1_j
        ct_out.c1{j} = mod(c1_q{j}, qj);
    end
end