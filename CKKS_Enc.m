function ct = CKKS_Enc(params, pk, m)
% CKKS_Enc_NTT  使用 NTT 加速实现 CKKS 第四章的 Enc 函数
%
% 输入：
%   params: 由 CKKS_SetUp 返回的结构体，包含以下字段：
%       .N            — 多项式环维度
%       .L            — 模数链深度
%       .Q_primes     — 长度 L+1 行向量 [q₀, …, q_L]
%       .encrypt_dist — 函数句柄，采样 v ← χ_enc（长度 N，系数 ∈ {0,±1}）
%       .error_dist   — 函数句柄，采样 e ← χ_err（长度 N，离散高斯）
%       .psi          — cell 数组，长度 = numel(params.D)，每个元素为该模数下的 ψ
%       .twiddles     — cell 数组，长度 = numel(params.D)，每个元素为长度 N 的 twiddle
%       因子(被取消)
%       .k            — special base 素数个数
%
%   pk: 结构体，公钥，包含：
%       .a{j}         — 对应 q_j 的 a_j (1×N 行向量)
%       .b{j}         — 对应 q_j 的 b_j (1×N 行向量)
%                    其中 j_matlab = 1..(L+1) 对应论文中的 j = 0..L
%
%   m:  1×N 整数行向量，表示要加密的明文多项式 m(X) 的系数
%       （在使用之前应当已近似编码并量化到整数域）
%
% 输出：
%   ct: 结构体，包含
%       .c0{j}  — 第一分量 c₀^{(j)} (1×N 多项式)，对应 b_j·v + m + e₀ mod q_j
%       .c1{j}  — 第二分量 c₁^{(j)} (1×N 多项式)，对应 a_j·v + e₁ mod q_j
%                j_matlab = 1..(L+1) 对应 j = 0..L
%
% 实现思路（参考论文第 4 章 Encpk 段落）：
%   1) 采样 v ← χ_enc；采样两个误差多项式 e0 ← χ_err, e1 ← χ_err
%   2) 对每个 j = 0..L:
%        • 计算 d0_j = b_j(X) ⋆ v(X) mod (X^N+1, q_j)  （用 NTT 加速）
%        • 计算 d1_j = a_j(X) ⋆ v(X) mod (X^N+1, q_j)  （用 NTT 加速）
%        • 将 e0, e1 对应取模 q_j；并将 m 对 q_j 取模
%        • c0_j = (d0_j + m + e0_j) mod q_j
%        • c1_j = (d1_j + e1_j)    mod q_j
%   3) 输出 ct.c0{j}, ct.c1{j}
%
% 注意：多项式乘法调用用户自行保存、并经 NTT 优化的 ckks_poly_mult_mod() 函数。

    %% 1. 准备
    N = params.N;
    L = params.L;
    Q_primes = params.Q_primes;         % [q₀, …, q_L]
    k = params.k;                       % special base 素数个数
    offset = k;                         % 在 params.psi / twiddles 中 q_j 的索引 = offset + j_matlab

    % 1.1 采样掩码多项式 v ← χ_enc
    v = params.encrypt_dist();          % 1×N，系数 ∈ {0,±1}

    % 1.2 采样误差多项式 e0, e1 ← χ_err
    e0 = params.error_dist();           % 1×N，整数
    e1 = params.error_dist();           % 1×N，整数

    % 确保 v, e0, e1 都是 1×N 行向量
    v  = reshape(v,  1, N);
    e0 = reshape(e0, 1, N);
    e1 = reshape(e1, 1, N);

    %% 2. 对每个 q_j 计算密文分量
    % 初始化输出 cell 数组
    ct.c0 = cell(1, L+1);
    ct.c1 = cell(1, L+1);

    for j_matlab = 1:(L+1)
        qj = Q_primes(j_matlab);
        psi_j = params.psi{offset + j_matlab};
        %twiddles_j = params.twiddles{offset + j_matlab};

        % 2.1 从公钥中取出 a_j, b_j（已经是 1×N 向量，系数 mod q_j）
        a_j = pk.a{j_matlab};   % 1×N
        b_j = pk.b{j_matlab};   % 1×N

        % 2.2 计算 b_j(X) ⋆ v(X) mod (X^N+1, q_j) → d0_j
        v_mod_j = mod(v, qj);
        d0_j = ckks_poly_mult_mod(b_j, v_mod_j, qj, N, psi_j);

        % 2.3 计算 a_j(X) ⋆ v(X) mod (X^N+1, q_j) → d1_j
        d1_j = ckks_poly_mult_mod(a_j, v_mod_j, qj, N, psi_j);

        % 2.5 将明文 m 对 q_j 取模
        m_j = mod(m, qj);       % 1×N

        % 2.6 计算密文两个分量
        c0_j = mod(d0_j + m_j + e0, qj);   % b_j·v + m + e0 mod q_j
        c1_j = mod(d1_j + e1,    qj);      % a_j·v + e1    mod q_j

        % 2.7 存储到输出
        ct.c0{j_matlab} = c0_j;
        ct.c1{j_matlab} = c1_j;
    end
end