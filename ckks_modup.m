function poly_pq = ckks_modup(poly_q, params)
% ckks_modup_one  对单个多项式执行 “ModUp” 基扩展
%
%   poly_p = ckks_modup_one(poly_q, params) 将输入多项式 poly_q 
%   （仅在 {q₀,…,qₗ} 这 l+1 个模下有表示）扩展到包含特殊基 
%   {p₀,…,p_{k-1}} 以及原本的这 l+1 个模。返回 poly_p 是长度 
%   k + l + 1 的 cell 数组，分别表示在每个 p_i 和 q_j 下的系数向量。
%
%   输入参数：
%     poly_q:  cell 数组，长度 = l+1，每个元素 poly_q{j+1} 是多项式在 q_j 下的系数向量
%     params:  由 CKKS_Setup 返回的结构体，至少包含：
%            • params.P_primes         : 1×k 行向量，特殊基 {p₀,…,p_{k-1}}
%            • params.Q_primes         : 1×(L+1) 行向量，近似基 {q₀,…,q_L}
%            • params.q_lj_mod_pi      : (L+1)×(L+1)×k 张量，
%                                       q_lj_mod_pi(l+1,j+1,i) = (∏_{t=0}^l q_t / q_j) mod p_i
%            • params.inv_q_lj_mod_qj  : (L+1)×(L+1) 矩阵，
%                                       inv_q_lj_mod_qj(l+1,j+1) = (∏_{t=0}^l q_t / q_j)^{-1} mod q_j
%
%   输出参数：
%     poly_p:  cell 数组，长度 = k + l + 1，每个元素是多项式在对应 p_i 或 q_j 下的系数向量
%              • poly_p{1..k}         = coefficients mod p₀..p_{k-1}
%              • poly_p{k+1..k+l+1}   = 原始 poly_q{1..l+1}（直接复制）
%
%   准备：必须已有以下自定义函数在路径中：
%     • modmul(A, b, m) : 向量 A × 标量 b mod m
%     • modadd(A, B, m) : 向量 A, B 相加 mod m
%
    Q_primes = params.Q_primes;            % 1×(L+1)
    P_primes = params.P_primes;            % 1×k
    q_lj_mod_pi     = params.q_lj_mod_pi;  % (L+1)×(L+1)×k
    inv_q_lj_mod_qj = params.inv_q_lj_mod_qj; % (L+1)×(L+1)
    M = numel(poly_q);       % = l+1
    k = params.k;     % 特殊基大小
    N = params.N;    %多项式长度

    % 输出长度 = k + (l+1)
    poly_pq = cell(1, k + M);
    
    % 对 p_i (i=1..k) 做 CRT 重构
    for i = 1:k
        p_i = P_primes(i);        
        acc = zeros(1, N);
        for j = 1:M
            q_j = Q_primes(j);
            inv_factor = inv_q_lj_mod_qj(M, j); % (Q_l/q_j)^{-1} mod q_j
            coeffs_qj  = poly_q{j};
            vj = modmul(coeffs_qj, inv_factor, q_j);
            factor_mod_pi = q_lj_mod_pi(M, j, i); % (Q_l / q_j) mod p_i
            vj_lifted = modmul(mod(vj, p_i), factor_mod_pi, p_i);
            acc = modadd(acc, vj_lifted, p_i);
        end
        poly_pq{i} = acc;
    end

    % 直接复制原本在 q₀..qₗ 下的分量
    for j = 1:M
        poly_pq{k + j} = poly_q{j};
    end
end