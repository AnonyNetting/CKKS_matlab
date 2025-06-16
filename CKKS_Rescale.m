function ct_res = CKKS_Rescale(ct, params)
% CKKS_Rescale  对给定的 level‐l ciphertext 做 RNS 形式的 “Rescale” (RS(ct))
%
% 输入：
%   ct     — 一个结构体，包含两个字段：c0 和 c1，
%             每个字段都是一个 1×(l+1) 的 cell 数组：
%             ct.c0{j+1} ∈ ℤ_{q_j}[X]/(X^N+1)，ct.c1{j+1} ∈ ℤ_{q_j}[X]/(X^N+1)，
%             对应 j = 0,1,…,l，且都存储为长度 N 的行向量（多项式系数）。
%   params — CKKS_Setup 中生成的参数结构体，其中至少包含：
%             params.Q_primes （长度 = L+1），
%             以及可调用的 modmul(), modinv() 函数。
%
% 说明：
%   假设原 ciphertext ct 在“level‐l” 上，也就是说：
%     ct.c0{1}, …, ct.c0{l+1}  分别是 c₀^{(0)}, …, c₀^{(l)}（模 q₀,…,q_l）  
%     ct.c1{1}, …, ct.c1{l+1}  分别是 c₁^{(0)}, …, c₁^{(l)}（模 q₀,…,q_l）
%   Rescale 的目标：去掉最高层模 q_l，将 ciphertext 由 level‐l 下移到 level‐(l−1)，
%   令 q_l_inv(j) = q_l^{-1} mod q_j，j = 0,…,l−1，则对于 i = 0,1 以及 j = 0,…,l−1：
%
%     c′_i^{(j)} = ( q_l_inv(j) ⋅ c_i^{(j)}  −  c_i^{(l)} )  mod q_j
%
%   输出新的 ct_res.c0 和 ct_res.c1，都是长度为 l 的 cell 数组，对应 j = 0..(l−1)。
%
% 输出：
%   ct_res — 一个结构体，含字段 c0, c1，
%            ct_res.c0{j+1}, ct_res.c1{j+1} 分别对应 j=0..(l−1) 的 c′₀^{(j)}, c′₁^{(j)}。

    %—— 提取 Q_primes 以及层数 l
    Q_primes = params.Q_primes;       % 大小为 1×(L+1)，其中 l ≤ L
    l = numel(ct.c0) - 1;             % 这里假设 ct.c0, ct.c1 长度都 = l+1

    %—— 最高层模 q_l
    ql = Q_primes(l+1);

    %—— 预分配输出 ciphertext (level‐(l−1) 将有 l 个模量)
    ct_res.c0 = cell(1, l);
    ct_res.c1 = cell(1, l);

    %—— 先计算 q_l^{-1} mod q_j 对于 j = 0..(l−1)
    %   存到一个长度为 l 的向量 inv_ql_mod_qj
    inv_ql_mod_qj = zeros(1, l);
    for j = 1:l
        qj = Q_primes(j);
        % modinv 返回 ql mod qj 的逆元（uint64 或 double）
        inv_ql_mod_qj(j) = modinv(mod(uint64(ql), uint64(qj)), qj);
    end

    %—— 完成 “Rescale” 计算
    for j = 1:l
        qj = Q_primes(j);

        % 原 ciphertext 在层 j 的两个多项式（长度为 N 的行向量）
        c0_j = ct.c0{j};
        c1_j = ct.c1{j};

        % 原 ciphertext 在层 l 的两个多项式，注意下标 l+1
        c0_l = ct.c0{l+1};
        c1_l = ct.c1{l+1};

        % 取出归一化常数 inv(q_l) mod q_j
        invql = inv_ql_mod_qj(j);

        %—— 计算 c′₀^{(j)} = ( inv(q_l mod q_j) * (c₀^{(j)}  –  c₀^{(l)}) ) mod q_j
        % 注意：c0_j, c0_l 都是长度 N 的行向量
        term0 = mod(c0_j - c0_l, qj );         % c0_j - c0_l mod q_j
        tmp0  = modmul(invql, term0, qj); % 再乘上 inv(q_l) mod q_j，整体 mod q_j
        ct_res.c0{j} = tmp0;

        %—— 计算 c′₁^{(j)} = ( inv(q_l mod q_j) * (c₁^{(j)}  –  c₁^{(l)}) ) mod q_j
        term1 = mod(c1_j - c1_l, qj );         % c0_j - c0_l mod q_j
        tmp1  = modmul(invql, term1, qj); % 再乘上 inv(q_l) mod q_j，整体 mod q_j
        ct_res.c1{j} = tmp1;
    end

    %—— ct_res.c0 和 ct_res.c1 长度都为 l，分别对应 j=1..l (即 j=0..l-1)
end