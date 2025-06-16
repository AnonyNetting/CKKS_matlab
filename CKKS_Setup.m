function params = CKKS_Setup(q0_bit_length, p_bit_length, qother_bit_length, k, L, eta, N)
    % CKKS_SetUp  根据 Cheon et al. (SAC 2018) 论文第 4 章预计算 CRT 和 NTT 常数
    %
    % 输入参数:
    %   q0_bit_length:      q0模数比特长度（如 60、50 等）
    %   p_bit_length:       p_i模数比特长度
    %   qother_bit_length:  除q0以外模数比特长度
    %   k:                  P中包含p_i个数
    %   L:                  模数链深度（即 Q = q₀, q₁, …, q_L 的个数减 1）
    %   N:                  多项式环维度（必须为 2 的幂，CKKS 中通常是 2^14, 2^15 等）
    %   eta:                近似基与特殊基选取时的精度范围（用以控制随机质数的位长偏差）
    %
    % 输出参数:
    %   params: 包含所有系统参数与预计算常数的结构体
    %
    % 参考：《A Full RNS Variant of Approximate Homomorphic Encryption》— 
    % …… 第 4 章 “Setup” 中关于 CRT 重构常数与 NTT 预计算的描述，以及附录 A.1 “Modulus Extension”。

    %% 1. 生成特殊素数基 B = { p₀, p₁, …, p_{k-1} }
    special_range = eta;                % 控制 p_i 位长时的搜索范围
    special_primes = ckks_generate_primes(k, N, p_bit_length, special_range);
    % special_primes 是长度为 k 的行向量，包含 k 个满足 p_i ≡ 1 (mod 2N) 的大素数

    %% 2. 生成近似基 C = { q₀, q₁, …, q_L }          
    approx_range = eta;                 % 控制 q_j 位长时的搜索范围
    % 为了让第一个 q₀ 稍微“宽裕”一些，常可多给几比特：
    q0 = ckks_generate_primes(1, N, q0_bit_length, approx_range);
    % 然后其余 L 个 q₁…q_L 使用标准的 modulus_bit_length
    q_other = ckks_generate_primes(L, N, qother_bit_length, approx_range);
    approx_primes = [q0, q_other];      % 长度为 (L+1)

    %% 3. 构建完整的模数集合 D = B ∪ C
    % 确认所有元素都两两不同（此处可以多写一点用于在发生重复时生成新的素数，在此暂不实现）
    if numel(unique([special_primes, approx_primes])) ~= numel([special_primes, approx_primes])
        error("存在RNS分量重复，需要调整指定比特宽度");
    end
    D = [special_primes, approx_primes];   % 总共 k + (L+1) 个模数
    m = numel(D);

    %% 4. 预计算 CRT 及 NTT 常数
    % 4.1 将特殊基与近似基分离为更直观的符号
    P_primes = special_primes;         % {p₀, …, p_{k-1}}
    Q_primes = approx_primes;          % {q₀, …, q_L}

    % 4.2 删去

    %% 4.3 计算以下常数：
    %   (i)   \hat{p}_i mod q_j  （记为 hatpi_mod_qj(i,j)）
    %  (ii)   (\hat{p}_i)^{-1} mod p_i  （记为 inv_hatpi_mod_pi(i)）
    % (iii)   P_total^{-1} mod q_j    （等同于 ∏_{t=0}^{k-1} p_t^{-1} mod q_j）
    %
    %  以及：
    %  (iv)   对于 0 ≤ j ≤ l ≤ L，q_{l,j} mod p_i  （q_lj_mod_pi）
    %   (v)    (q_{l,j})^{-1} mod q_j      （inv_q_lj_mod_qj）

    % 预分配
    hatpi_mod_qj     = zeros(k, (L+1));     % \hat{p}_i mod q_j
    inv_hatpi_mod_pi = zeros(1, k);              % (\hat{p}_i)^{-1} mod p_i
    P_inv_mod_qj     = zeros(1, (L+1));     % P_total^{-1} mod q_j
    q_lj_mod_pi      = zeros(L+1, L+1, k);       % q_{l,j} mod p_i
    inv_q_lj_mod_qj  = zeros(L+1, L+1);          % (q_{l,j})^{-1} mod q_j

    % 为了使用 modmul, modinv 等自定义函数，先确保它们在路径中：
    %   modmul(a,b,q) 返回 (a*b) mod q，使用拆分法避免超过 2^52；
    %   modinv(a,m)  返回 a^{-1} mod m，使用整型扩展欧几里得。

    %% 4.4 计算 \hat{p}_i mod qj 以及其逆
    for i = 1:k
        % 4.4.1 计算 (\hat{p}_i) mod p_i 以及其逆
        %   \hat{p}_i = ∏_{t≠i} P_primes(t). 为避免大整数直接相乘，采用循环 modmul：
        pi_val = P_primes(i);               % p_i
        % 先计算 \hat{p}_i mod p_i：
        hatp_mod_pi = uint64(1);
        for t = 1:k
            if t == i, continue; end
            hatp_mod_pi = modmul(hatp_mod_pi, mod(P_primes(t), pi_val), pi_val);
        end
        % hatp_mod_pi ∈ [0, p_i); 由于 P_primes 中各素数两两互素，hatp_mod_pi ≠ 0
        inv_hatpi_mod_pi(i) = modinv(hatp_mod_pi, pi_val);

        % 4.4.2 计算 \hat{p}_i mod q_j, j = 1..(L+1)
        for j = 1:(L+1)
            qj = Q_primes(j);
            hatp_mod_qj = uint64(1);
            for t = 1:k
                if t == i, continue; end
                hatp_mod_qj = modmul(hatp_mod_qj, mod(P_primes(t), qj), qj);
            end
            hatpi_mod_qj(i, j) = double(hatp_mod_qj);
        end
    end

    %% 4.5 计算 P_total^{-1} mod q_j
    for j = 1:(L+1)
        qj = Q_primes(j);
        % 先用循环逐个乘 P_primes，再 mod qj：
        P_mod_qj = uint64(1);
        for t = 1:k
            P_mod_qj = modmul(P_mod_qj, mod(P_primes(t), qj), qj);
        end
        % P_mod_qj = (∏ p_i) mod q_j
        P_inv_mod_qj(j) = modinv(P_mod_qj, qj);
    end

    %% 4.6 计算 q_{l,j} mod p_i 以及其逆 mod q_j
    % 对于每个 l = 0..L，对应 MATLAB 索引 1..(L+1)
    for l_idx = 1:(L+1)
        % l_paper = l_idx - 1
        % 考虑 q_0, q_1, ..., q_l_paper → 对应 Q_primes(1:l_idx)
        for j_idx = 1:l_idx
            % j_paper = j_idx - 1
            % q_{l,j} = ∏_{t=0..l, t≠j} q_t  → 即 ∏ Q_primes(1..l_idx) 除去 Q_primes(j_idx)
            % 逐元素乘法并取模来计算三个对象：
            %   (a) q_{l,j} mod p_i  （对 i = 1..k）
            %   (b) q_{l,j} mod q_j  （其中 q_j = Q_primes(j_idx) → 用于求逆）

            % 先对 i = 1..k 计算 q_{l,j} mod p_i
            for i = 1:k
                pi_val = P_primes(i);
                qlj_mod_pi = uint64(1);
                for t = 1:l_idx
                    if t == j_idx, continue; end
                    qlj_mod_pi = modmul(qlj_mod_pi, mod(Q_primes(t), pi_val), pi_val);
                end
                q_lj_mod_pi(l_idx, j_idx, i) = double(qlj_mod_pi);
            end

            % 计算 (q_{l,j})^{-1} mod q_j
            qj = Q_primes(j_idx);
            qlj_mod_qj = uint64(1);
            for t = 1:l_idx
                if t == j_idx, continue; end
                qlj_mod_qj = modmul(qlj_mod_qj, mod(Q_primes(t), qj), qj);
            end
            inv_q_lj_mod_qj(l_idx, j_idx) = modinv(qlj_mod_qj, qj);
        end
    end

    %% 4.7 将所有 CRT 相关常数保存到 params
    params.P_primes            = P_primes;            
    params.Q_primes            = Q_primes;            
    params.hatpi_mod_qj        = hatpi_mod_qj;            
    params.inv_hatpi_mod_pi    = inv_hatpi_mod_pi;      
    params.P_inv_mod_qj        = P_inv_mod_qj;          
    params.q_lj_mod_pi         = q_lj_mod_pi;          
    params.inv_q_lj_mod_qj     = inv_q_lj_mod_qj;      

    %% 4.8 预计算 NTT 相关常数
    params.psi = cell(1, m);
    % params.twiddles = cell(1, m);

    for i = 1:m
        qi = D(i);
        % 确保每个模数满足 qi ≡ 1 (mod 2N)
        if mod(qi - 1, 2 * N) ~= 0
            error('CKKS_SetUp: D(%d) = %d 不满足 ≡ 1 (mod 2N)。', i, qi);
        end
        % 找到 psi_i 使得 psi_i^N ≡ -1 (mod qi)
        psi_i = find_psi(qi, N);
        params.psi{i} = psi_i;

        % 若需要 twiddle 表，可按需启用：
        % twiddle = zeros(1, N);
        % for j = 1:N
        %     exponent = 2*(j-1) + 1;
        %     % powmod 可以直接计算大指数的模幂
        %     twiddle(j) = powmod(sym(psi_i), sym(exponent), sym(qi));
        % end
        % params.twiddles{i} = double(twiddle);
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
    params.N                   = N;
    params.L                   = L;
    params.k                   = k;
    params.special_primes      = special_primes;
    params.approx_primes       = approx_primes;
    params.D                   = D;
    params.qother_bit_length  = qother_bit_length;
    params.eta                 = eta;

    %% 7. 输出提示
    fprintf('CKKS 参数设置完成:\n');
    fprintf('  多项式维度 N = %d\n', N);
    fprintf('  模数链深度 L = %d\n', L);
    fprintf('  特殊基 p_i（k = %d 个）: %s\n', k, mat2str(special_primes));
    fprintf('  近似基 q_j（L+1 = %d 个）: %s\n', (L+1), mat2str(approx_primes));
    fprintf('  已预计算 CRT 重构常数和 NTT 预计算常数。\n');
end

%% =============================================================================
%% 辅助函数：寻找满足 psi^N ≡ -1 (mod q) 的原根 psi
function psi = find_psi(q, N)
    % find_psi: 在 q ≡ 1 (mod 2N) 的前提下，找到 psi 使得 psi^N ≡ -1 (mod q)

    q_sym = sym(q);
    N_sym = sym(N);

    % 1. 检查 q ≡ 1 (mod 2N)
    if mod(q_sym - 1, 2 * N_sym) ~= 0
        error('find_psi: q 不满足 q ≡ 1 (mod 2N)。');
    end

    % 2. 计算 d = (q–1) / (2N)
    d_sym = (q_sym - 1) / (2 * N_sym);

    % 3. 构造 q−1 的不重复质因子集合
    %    q−1 = 2 * N * d，质因子来自 {2} ∪ factor(N) ∪ factor(d)
    d_factors_all    = factor(d_sym);
    d_factors_unique = unique(d_factors_all);
    N_factors        = unique(factor(N_sym));
    facs_qm1         = unique([ sym(2), N_factors, d_factors_unique ]);

    % 4. 枚举 g，直至找到原根
    psi = [];
    for g_val = uint64(2) : uint64(q_sym - 2)
        g_sym = sym(g_val);
        is_generator = true;
        for idx = 1:numel(facs_qm1)
            pi = facs_qm1(idx);
            exp_i = (q_sym - 1) / pi;
            if powermod(g_sym, exp_i, q_sym) == 1
                is_generator = false;
                break;
            end
        end
        if is_generator
            % g_val 为 q 的一个原根
            psi = powermod(g_sym, d_sym, q_sym);
            break;
        end
    end
    if isempty(psi)
        error('find_psi: 未找到原根；请检查 q 是否为素数。');
    end

    % 验证 psi^N ≡ -1 (mod q)：
    assert(powermod(psi, N_sym, q_sym) == q_sym - 1, 'psi^N ≠ -1');
    psi = double(psi);
end

%% =============================================================================
%% 辅助函数：采样 - 汉明重量分布 χ_key
function vec = sample_hamming_weight(n, h)
    if h > n
        error('sample_hamming_weight: 要求的汉明重量 h 超过 n。');
    end
    vec = zeros(1, n);
    idx = randperm(n, h);
    signs = rand(1, h) > 0.5;  % true → -1, false → +1
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
    vec = round(randn(1, n) * sigma);
end

%% 辅助函数：采样 - 三元均匀分布 χ_enc
function vec = sample_ternary_uniform(n)
    r = rand(1, n);
    vec = zeros(1, n);
    vec(r < 0.25)  = +1;
    vec(r >= 0.75) = -1;
end