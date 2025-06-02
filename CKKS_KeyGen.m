%% CKKS_KeyGen_NTT.m
% 使用 NTT/INTT 加速多项式乘法的 CKKS KeyGen 实现
%
% 假定 params = CKKS_SetUp(...) 已经运行，包含以下字段：
%   params.N            % 多项式环维度（2 的幂）
%   params.L            % 模数链深度
%   params.P_primes     % special base primes {p_i}
%   params.Q_primes     % approx base primes {q_j}
%   params.D            % 全部模数集合 [p_0,...,p_{k-1}, q_0,...,q_L]
%   params.psi          % cell 数组，长度 = numel(params.D)，每个元素为该模数下的 ψ（2N 次单位根）
%   params.twiddles     % cell 数组，长度 = numel(params.D)，每个元素为长度 N 的 twiddle
%   因子(暂时被取消)
%   params.key_dist     % 函数句柄：采样私钥多项式 s (长度 N)
%   params.error_dist   % 函数句柄：采样误差多项式 e (长度 N)
%
% 目标：生成 (sk, pk)：
%   sk = s(X) ∈ R₂ (隐含 mod 2)
%   pk = { (a_j(X), b_j(X)) ∈ R_{q_j}^2 : j = 0..L },
%   b_j(X) = a_j(X)·s(X) + e_j(X)  (mod X^N+1, mod q_j)
%
% 其中的多项式乘法都用 NTT/INTT 加速，用到 params.psi。

function keys = CKKS_KeyGen(params)
    %% 1. 准备工作
    N = params.N;
    L = params.L;
    P_primes = params.P_primes;
    Q_primes = params.Q_primes;        % 长度为 L+1
    D = params.D;                      % 长度为 k + (L+1)
    k = numel(P_primes);
    
    % offset：在 D 中，q_j 对应的索引 = k + (j+1)，其中 j=0..L 对应 MATLAB 索引 j_matlab=1..L+1
    offset = k;
    
    % 私钥 s ← χ_key
    s = params.key_dist();             % 1×N，元素 ∈ {0,±1}
    s = reshape(s, 1, N);              % 确保行向量

    % 准备输出结构
    pk.a = cell(1, L+1);
    pk.b = cell(1, L+1);

    % 2.2 采样误差 e
    e = params.error_dist();     % 1×N，可能正负
    
    %% 2. 对每个模 q_j = Q_primes(j) 生成 a_j, e_j, 并计算 b_j
    for j_matlab = 1:(L+1)
        % j_paper = j_matlab - 1
        qj = Q_primes(j_matlab);
        D_idx = offset + j_matlab;     % 在 D 中对应的索引，找 ψ

        % 2.1 采样 a_j: 均匀随机的多项式，系数在 [0, q_j-1]
        a_j = randi([0, qj-1], 1, N);

        % 2.2 用 NTT/INTT 加速多项式乘法 c = a_j * s (mod X^N+1, mod q_j)
        %     poly_mult_mod_ntt 会：
        %       (1) 将 a_j、s 做 NTT(·, q_j, ψ)
        %       (2) 在频域做点乘
        %       (3) 做 INTT(·, q_j, ψ^{-1})，并除以 N mod q_j
        c = ckks_poly_mult_mod(a_j, s, qj, N, params.psi{D_idx});

        % 2.3 计算 b_j = (c + e_j) mod q_j
        b_j = mod(- c + e, qj);

        % 2.4 存储到公钥
        pk.a{j_matlab} = a_j;
        pk.b{j_matlab} = b_j;
    end

    %% 3. 输出私钥与公钥
    keys.sk = s;
    keys.pk = pk;
end