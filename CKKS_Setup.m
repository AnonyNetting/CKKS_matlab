function params = CKKS_Setup(modulus_bit_length, L, eta, N)
    % CKKS_SetUp  根据 Cheon et al. (SAC 2018) 论文第 4 章预计算 CRT 和 NTT 常数
    %
    % 输入参数:
    %   modulus_bit_length: 模数比特长度（如 60、50 等）
    %   L:                  模数链深度（即 q₀, q₁, …, q_L 的个数减 1）
    %   N:                  多项式环维度（必须为 2 的幂，CKKS 中通常是 2^14, 2^15 等）
    %   eta:                近似基与特殊基选取时的精度范围（用以控制随机质数的位长偏差）
    %
    % 输出参数:
    %   params: 包含所有系统参数与预计算常数的结构体
    %
    % 参考：《A Full RNS Variant of Approximate Homomorphic Encryption》—
    % …… 第 4 章 “Setup” 中关于 CRT 重构常数与 NTT 预计算的描述，以及附录 A.1 “Modulus Extension”。
    %

    %% 1. 生成特殊素数基 B = { p₀, p₁, …, p_{k-1} }
    %   k 的值依安全性需求而定，此处示例设 k = 1，可根据实际可扩展。
    k = 1;                              % 特殊基中素数个数
    special_range = eta;                % 控制 p_i 位长时的搜索范围
    special_primes = ckks_generate_primes(k, N, modulus_bit_length + 10, special_range);
    % special_primes 是长度为 k 的行向量，包含 k 个满足 p_i ≡ 1 (mod 2N) 的大素数
    %
    % 论文中要求：每个 p_i 满足 p_i ≡ 1 mod (2N)，以保证在 R_{p_i}=Z_{p_i}[X]/(X^N+1) 中存在 N 次原根。

    %% 2. 生成近似基 C = { q₀, q₁, …, q_L }
    num_primes = L + 1;                 
    approx_range = eta;                 % 控制 q_j 位长时的搜索范围
    % 为了让第一个 q₀ 稍微“宽裕”一些，常可多给几比特：
    q0 = ckks_generate_primes(1, N, modulus_bit_length + 10, approx_range);
    % 然后其余 L 个 q₁…q_L 使用标准的 modulus_bit_length
    q_other = ckks_generate_primes(L, N, modulus_bit_length, approx_range);
    approx_primes = [q0, q_other];      % 长度为 (L+1)

    %% 3. 构建完整的模数集合 D = B ∪ C
    % 确认所有元素都两两不同(在p_i设置中心点与q_0相同的情况下，q_0会和p_i相同)
    % 若发现冲突就重生出冲突 p_i (或也可以重生 q_0)，直到两者都不同为止。
    if special_primes(1) == approx_primes(1)
        tmp = ckks_generate_primes(1 + k, N, modulus_bit_length + 10, approx_range);
        special_primes = tmp(2 : k+1);
    end
    D = [special_primes, approx_primes];   % 总共 k + (L+1) 个模数
    m = numel(D);                         % 模数总数 m = k + (L+1)

    %% 4. 预计算 CRT 及 NTT 常数
    % 4.1 将特殊基与近似基分离为更直观的符号
    P_primes = special_primes;         % {p₀, p₁, …, p_{k-1}}
    Q_primes = approx_primes;          % {q₀, q₁, …, q_L}

    % 4.2 计算 P_total = ∏_{i=0}^{k-1} p_i
    P_total = prod(P_primes);

    % 4.3 对于每个 i ∈ [1..k]，计算 \hat{p}_i = P_total / p_i
    hat_p = zeros(1, k);
    for i = 1:k
        hat_p(i) = P_total / P_primes(i);
    end

    % 4.4 计算常数: 
    %     (i)  \hat{p}_i mod q_j （论文记为 “p_i q_j”），对 i=0..k-1, j=0..L
    %    (ii)  (\hat{p}_i)^{-1} mod p_i 
    %   (iii)  P_total^{-1} mod q_j, 等同于 ∏_{i=0}^{k-1} (p_i)^{-1} mod q_j
    % 注意：Matlab 索引从 1 开始，但论文 i, j 从 0 开始。我们采用 MATLAB 索引时，
    %      i_matlab = i_paper + 1, j_matlab = j_paper + 1，以此类推。

    % 预分配
    pi_mod_q = zeros(k, num_primes);          % pi_mod_q(i, j) = \hat{p}_i mod q_j
    inv_hatp_mod_pi = zeros(1, k);            % inv_hatp_mod_pi(i) = (\hat{p}_i)^{-1} mod p_i
    P_inv_mod_q = zeros(1, num_primes);       % P_inv_mod_q(j) = P_total^{-1} mod q_j

    for i = 1:k
        % 4.4.1 计算 (\hat{p}_i)^{-1} mod p_i
        % 首先 \hat{p}_i mod p_i 恰好是 P_total/p_i mod p_i，因为 P_total = p_i * \hat{p}_i
        % 因此 \hat{p}_i mod p_i 必不为 0（p_i 与其他素数互素）。
        a_mod_pi = mod(hat_p(i), P_primes(i));
        inv_hatp_mod_pi(i) = modinv(a_mod_pi, P_primes(i));

        % 4.4.2 对于每个 q_j, 计算 \hat{p}_i mod q_j
        for j = 1:num_primes
            pi_mod_q(i, j) = mod(hat_p(i), Q_primes(j));
        end
    end

    % 4.4.3 计算 P_total^{-1} mod q_j → 等价于模 q_j 下 ∏_{i=0}^{k-1} p_i^{-1}
    for j = 1:num_primes
        % 先计算 P_mod_qj = P_total mod q_j
        P_mod_qj = mod(P_total, Q_primes(j));
        % 再求逆
        P_inv_mod_q(j) = modinv(P_mod_qj, Q_primes(j));
    end

    % 4.5 计算 q_{l, j} 及相关常数：  
    %     对于 0 ≤ j ≤ l ≤ L（MATLAB 中 j_matlab=1..(l+1), l_matlab=1..(L+1)）,
    %     定义 q_{l,j} = ∏_{0 ≤ j' ≤ l, j' ≠ j} q_{j'}.
    %     需要计算：
    %       (i)   q_{l,j} mod p_i   （即存储为 q_lj_mod_pi(l+1, j+1, i)）
    %      (ii)   (q_{l,j})^{-1} mod q_j   （即存储为 inv_q_lj_mod_qj(l+1, j+1)）
    %
    %   在 MATLAB 中，l_matlab = l_paper + 1, j_matlab = j_paper + 1.

    % 预分配张量
    q_lj_mod_pi = zeros(L+1, L+1, k);
    inv_q_lj_mod_qj = zeros(L+1, L+1);

    % 逐层 l = 0..L（在 MATLAB 对应 l_matlab = 1..L+1）
    for l_matlab = 1:(L+1)
        % l_matlab - 1 = l_paper
        for j_matlab = 1:l_matlab
            % j_matlab - 1 = j_paper
            % 构造 j' 的索引集合： 1..l_matlab, 去掉 j_matlab
            others = [1:(j_matlab-1), (j_matlab+1):l_matlab];
            % 计算 q_{l,j} = ∏_{j' ∈ others} Q_primes(j')
            qlj = 1;
            for idx = others
                qlj = qlj * Q_primes(idx);
            end
            % 然后对每个 p_i 取模
            for i = 1:k
                q_lj_mod_pi(l_matlab, j_matlab, i) = mod(qlj, P_primes(i));
            end
            % 最后计算 (q_{l,j})^{-1} mod q_j
            qlj_mod_qj = mod(qlj, Q_primes(j_matlab));
            inv_q_lj_mod_qj(l_matlab, j_matlab) = modinv(qlj_mod_qj, Q_primes(j_matlab));
        end
    end

    % 4.6 将所有 CRT 相关常数保存到 params
    params.P_primes            = P_primes;             % B 基
    params.Q_primes            = Q_primes;             % C 基
    params.P_total             = P_total;              % ∏ p_i
    params.hat_p               = hat_p;                % \hat{p}_i = P_total / p_i
    params.pi_mod_q            = pi_mod_q;             % \hat{p}_i mod q_j
    params.inv_hatp_mod_pi     = inv_hatp_mod_pi;      % (\hat{p}_i)^{-1} mod p_i
    params.P_inv_mod_q         = P_inv_mod_q;          % P_total^{-1} mod q_j
    params.q_lj_mod_pi         = q_lj_mod_pi;          % q_{l,j} mod p_i
    params.inv_q_lj_mod_qj     = inv_q_lj_mod_qj;      % (q_{l,j})^{-1} mod q_j

    %% 4.7 预计算 NTT 相关常数（同原函数）
    params.psi = cell(1, m);         % 存储每个模数 D(i) 下的 psi_i
    % params.twiddles = cell(1, m);    % 存储每个模数 D(i) 下的长度 N 的 twiddle 因子

    for i = 1:m
        qi = D(i);
        % 确保每个模数满足 qi ≡ 1 (mod 2N)，以便能找到 N 次原根
        if mod(qi - 1, 2 * N) ~= 0
            error('CKKS_SetUp: D(%d) = %d 不满足 ≡ 1 (mod 2N)。', i, qi);
        end
        % 找到 psi_i 使得 psi_i^N ≡ -1 (mod qi)
        psi_i = find_psi(qi, N);
        params.psi{i} = psi_i;

        % 计算 twiddle 因子：twiddle(j) = psi_i^(2*(j-1)+1) mod qi, j=1..N
        % twiddle = zeros(1, N);
        % for j = 1:N
        %     exponent = 2*(j-1) + 1;
        %     twiddle(j) = modpow(psi_i, exponent, qi);
        % end
        % params.twiddles{i} = twiddle;
    end

    %% 5. 选择概率分布（与原文附录 A.2 保持一致）
    params.h = 64;               % 私钥的汉明重量 h（论文示例 h=64）
    params.sigma = 3.2;          % 离散高斯分布的标准差 σ

    % 密钥分布 χ_key：均匀采样长度 N，汉明重量 = h，系数 ∈ {±1}
    params.key_dist = @() sample_hamming_weight(N, params.h);

    % 加密掩码分布 χ_enc：三元均匀分布 {0, ±1}
    params.encrypt_dist = @() sample_ternary_uniform(N);

    % 误差分布 χ_err：离散高斯分布（均值 0, 标准差 σ）
    params.error_dist = @() sample_discrete_gaussian(N, params.sigma);

    %% 6. 记录其他基本参数
    params.N = N;
    params.L = L;
    params.k = k;
    params.special_primes = special_primes;
    params.approx_primes = approx_primes;
    params.D = D;
    params.modulus_bit_length = modulus_bit_length;
    params.eta = eta;

    %% 7. 输出提示
    fprintf('CKKS 参数设置完成:\n');
    fprintf('  多项式维度 N = %d\n', N);
    fprintf('  模数链深度 L = %d\n', L);
    fprintf('  特殊基 p_i（k = %d 个）: %s\n', k, mat2str(special_primes));
    fprintf('  近似基 q_j（L+1 = %d 个）: %s\n', num_primes, mat2str(approx_primes));
    fprintf('  已预计算 CRT 重构常数和 NTT 预计算常数。\n');
end

%% =============================================================================
%% 辅助函数：寻找满足 psi^N ≡ -1 (mod q) 的原根 psi
function psi = find_psi(q, N)
    % 输入合法性检查并转换为 symbolic
    q_sym = sym(q);
    N_sym = sym(N);
    
    % 1. 检查 q ≡ 1 (mod 2N)
    if mod(q_sym - 1, 2 * N_sym) ~= 0
        error('find_psi: q 不满足 q ≡ 1 (mod 2N)。');
    end
    
    % 2. 计算 d = (q–1) / (2N)
    d_sym = (q_sym - 1) / (2 * N_sym);
    
    % 3. 对 d 做质因数分解，提取不重复质因子
    %    这里直接调用 sym/factor() 返回 d 的所有质因子（含重复）
    d_factors_all = factor(d_sym);       % 例如 [p1, p1, p2, p3, …]
    d_factors_unique = unique(d_factors_all);
    % d_factors_unique 现在是一个行向量，包含 d 的所有不同质因子
    
    % 4. 构造 “q-1 的全部不重复质因子集合” facs_qm1
    %    q-1 = 2 * N * d。其质因子必然包含 {2} 与 factor(N)、以及 d_factors_unique。
    N_factors = unique(factor(N_sym));               % N 本身的质因子
    facs_qm1 = unique([ sym(2), N_factors, d_factors_unique ]);
    % facs_qm1 为行向量，存储 q-1 的所有不同质因子
    
    % 5. 枚举 g，查找原根（generator）
    psi = [];   % 先置为空
    for g_val = uint64(2) : uint64(q_sym - 1)
        g_sym = sym(g_val);
        is_generator = true;
        % 对于 q-1 的每个质因子 p_i，测试 g^((q-1)/p_i) mod q 是否不等于 1
        for idx = 1 : numel(facs_qm1)
            pi = facs_qm1(idx);
            exp_i = (q_sym - 1) / pi;              % (q-1)/p_i
            % 使用 powermod（内含快速模幂算法）
            r = powermod(g_sym, exp_i, q_sym);     
            if r == 1
                is_generator = false;
                break;
            end
        end
        if is_generator
            % g_sym 现在是 “q 的一个原根”（阶 = q-1）
            % 6. 令 psi = g^d mod q；输出即可
            psi = powermod(g_sym, d_sym, q_sym);
            break;
        end
        % 否则继续尝试下一个 g_val
    end
    
    if isempty(psi)
        error('find_psi: 未找到原根；请检查 q 是否真的是素数 或 N/(2N) 设置是否正确。');
    end
    
    % （可选）验证：
    %   assert(powermod(psi, N_sym, q_sym) == q_sym - 1, 'psi^N ≠ -1');
    assert(powermod(psi, 2*N_sym, q_sym) == 1,        'psi^(2N) ≠ 1');
    psi = double(psi);
    % 最终输出 psi，保持 symbolic 类型。如果后续需要数值，可以调用 double(psi) 但注意位宽。
end

%% =============================================================================
%% 辅助函数：采样 - 汉明重量分布 χ_key
function vec = sample_hamming_weight(n, h)
    % sample_hamming_weight: 在 { -1, 0, +1 }^n 中均匀采样一个向量，
    %   恰好有 h 个非零位（±1，各 50% 概率），其余为 0。
    if h > n
        error('sample_hamming_weight: 要求的汉明重量 h 超过 n。');
    end
    vec = zeros(1, n);
    idx = randperm(n, h);
    signs = rand(1, h) > 0.5;  % 逻辑数组：true → -1, false → +1
    for t = 1:h
        if signs(t)
            vec(idx(t)) = -1;
        else
            vec(idx(t)) = +1;
        end
    end
end

%% 辅助函数：采样 - 离散高斯分布 χ_err
function vec = sample_discrete_gaussian(n, sigma)
    % sample_discrete_gaussian: 对每个系数独立采样 N(0, σ^2) 后四舍五入
    vec = round(randn(1, n) * sigma);
end

%% 辅助函数：采样 - 三元均匀分布 χ_enc
function vec = sample_ternary_uniform(n)
    % sample_ternary_uniform: 每个系数以 Prob=1/4 取 +1，1/4 取 -1，其余 1/2 取 0
    r = rand(1, n);
    vec = zeros(1, n);
    vec(r < 0.25)    =  1;
    vec(r >= 0.75)   = -1;
    % r ∈ [0.25, 0.75) 时保持为 0
end