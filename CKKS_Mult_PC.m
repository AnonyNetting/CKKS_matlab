function ct_out = CKKS_Mult_PC(params, ct, m)
% CKKS_Mult_PC  “明文×密文” 的同态乘法 (Plaintext × Ciphertext)，
% 利用 NTT/INTT 做多项式环 ℤ_{q_j}[X]/(X^N+1) 下的卷积。
%
% 输入：
%   params — CKKS_Setup(...) 返回的参数结构体，至少包含：
%            • params.Q_primes      — 行向量 [q₀, q₁, …, q_L]
%            • params.k             — 特殊模数量 k （默认为 CKKS_Setup 中设置的 k）
%            • params.psi           — cell 数组，长度 = k + (L+1)，
%                                   其中 ψ 对应 D = [p₀,…,p_{k-1}, q₀,…,q_L]
%            • params.N             — 多项式维度 N
%            • 并且路径中已有函数 ckks_poly_mult_mod(f,g,q,N,psi)、
%              modmul()、intt()、ntt()、modinv() 等。
%
%   ct     — 要乘的密文结构体，形如：
%            ct.c0 = { c0^{(0)}, c0^{(1)}, …, c0^{(L)} };
%            ct.c1 = { c1^{(0)}, c1^{(1)}, …, c1^{(L)} };
%            每个 c{i}^{(j)} ∈ ℤ_{q_j}[X]/(X^N+1) 都是一个 1×N 的行向量。
%            这里假设 ct 一共有 (L+1) 个 RNS 分量，对应 j = 0..L。
%
%   m      — 明文多项式（已经调用 ckks_encode），是一个长度为 N 的行向量。
%            在 RNS 体系下，我们需要对每个 q_j 先做模减至 [0, q_j)，
%            再与 c{i}^{(j)} 做 NTT 卷积。
%
% 输出：
%   ct_out — 乘完之后的新密文，结构与输入 ct 相同：
%            ct_out.c0 = { c0'^{(0)}, …, c0'^{(L)} };
%            ct_out.c1 = { c1'^{(0)}, …, c1'^{(L)} }。
%            每个 c{i}'^{(j)} = ( m · c{i}^{(j)} ) mod (X^N+1, q_j).
%
% 核心思路：
%   对每个 RNS 分量 j = 0..L：
%     1) f_j = mod(m, q_j)  → 将明文系数还原到 ℤ_{q_j}
%     2) c{i}_j = ct.c{i}{j+1}, i=0,1
%     3) ψ_j = params.psi{ k + j + 1 }  → 找到对应 q_j 的 2N 次原根 ψ
%     4) 用 ckks_poly_mult_mod(f_j, c{i}_j, q_j, N, ψ_j) 计算 negacyclic 卷积
%
    %—— 提取关键信息
    Q_primes = params.Q_primes;        % [q₀, q₁, …, q_L]
    L_plus1  = numel(Q_primes);        % 等于密文层数 L+1
    N        = params.N;               % 多项式维度
    k        = params.k;               % 特殊模个数
    psi_all  = params.psi;             % 总共 k + (L+1) 个 ψ 值:
                                       %   psi_all{1..k} 对应 p₀..p_{k-1},
                                       %   psi_all{k+1..k+L+1} 对应 q₀..q_L。

    %—— 预分配输出
    ct_out.c0 = cell(1, L_plus1);
    ct_out.c1 = cell(1, L_plus1);

    %—— 对每个 RNS 分量做 “m × c{i}^{(j)} mod (X^N+1, q_j)”：
    for j = 0:(L_plus1-1)
        qj   = Q_primes(j+1);                  % 第 j 层的模数 q_j
        psi_j = psi_all{k + (j+1)};            % 对应 q_j 的 ψ，存储在 params.psi{k+j+1}

        % 1) 把明文 m 在 ℤ_{q_j} 上取模，作为 f_j
        f_j = mod( double(m), double(qj) );    % 1×N 行向量，范围 [0, q_j)

        % 2) 从输入密文取出 c0^{(j)}, c1^{(j)}
        c0_j = ct.c0{j+1};  % 已知是 1×N 向量，整数或 double
        c1_j = ct.c1{j+1};

        %—— 3) 用 NTT+INTT 做多项式乘法
        %     ckks_poly_mult_mod 已封装 “negacyclic NTT 卷积 + mod”
        c0p_j = ckks_poly_mult_mod(f_j, c0_j, qj, N, psi_j);
        c1p_j = ckks_poly_mult_mod(f_j, c1_j, qj, N, psi_j);

        % 4) 存回输出
        ct_out.c0{j+1} = c0p_j;   % 还是 1×N 行向量，系数 ∈ [0, q_j)
        ct_out.c1{j+1} = c1p_j;
    end
end