function m_int = CKKS_Dec(params, sk, ct)
% CKKS_Dec  使用 NTT 加速实现 CKKS 解密，输出可直接交给 ckks_decode 的整数向量
%
% 输入：
%   params: 由 CKKS_SetUp 返回的结构体，包含：
%       .N            — 多项式维度
%       .Q_primes     — 长度 L+1 行向量 [q₀, …, q_L]
%       .psi          — cell 数组，长度 = numel(params.D)，存储每个模数下的 2N 次原根
%       .twiddles     — cell 数组，长度 = numel(params.D)，存储每个模数下的长度 N
%       twiddle(被取消)
%       .k            — special base 素数个数
%
%   sk:  1×N 私钥多项式 s（系数 ∈ {0,±1}，根据 key_dist 采样）
%
%   ct:  密文结构体，包含
%       .c0{j}  — 第一分量 c₀^{(j)} (1×N) ， j_matlab=1..(L+1) 对应 j=0..L
%       .c1{j}  — 第二分量 c₁^{(j)} (1×N)
%
% 输出：
%   m_int: 1×N 整数行向量，已从 c0,c1 恢复的明文“整数多项式”值，可直接传入 ckks_decode 进行浮点解码
%
% 解密思路（论文第 4 章 Dec 步骤）：
%   1) 只使用第 j=0 层（MATLAB 索引 j_matlab=1，对应 q₀）的 c₀, c₁
%   2) 计算 t(X) = c₀^{(0)}(X) – c₁^{(0)}(X) · s(X)  (mod X^N+1, q₀)
%   3) 对 t 中每个系数做“中心化”整数映射：若 > q₀/2，则减去 q₀
%   4) 返回长度 N 的整数系数向量 m_int
%
% `m_int` 就是“Δ·m”近似整数多项式，可直接交给 ckks_decode 反量化得到近似浮点明文。

    %% 1. 基本参数和第 0 层模数 q0
    N = params.N;
    q0 = params.Q_primes(1);            % j_matlab=1 对应论文 j=0
    k = params.k;
    offset = k;                         % c0/c1 对应 psi、twiddles 的起始索引
    
    % 1.1 取出 c0, c1 在第 0 层（j_matlab = 1）
    c0 = ct.c0{1};   % 1×N，系数 ∈ [0, q0-1]
    c1 = ct.c1{1};   % 1×N
    
    % 1.2 取出 sk：长度 N，系数 ∈ {0,±1}
    s = reshape(sk, 1, N);
    
    % 1.3 获取对应的 NTT 根与 twiddles
    psi_q0      = params.psi{offset + 1};      % 对应 q0
    %twiddles_q0 = params.twiddles{offset + 1}; % 对应 q0

    %% 2. 计算 c1(X) * s(X) mod (X^N + 1, q0) → cs
    %     使用 NTT 加速的多项式乘法
    s  = mod(s, q0);
    cs = ckks_poly_mult_mod(c1, s, q0, N, psi_q0);
    % cs 也是 1×N，系数 ∈ [0, q0-1]

    %% 3. 计算 t = c0 + cs (mod q0)
    t = mod(c0 + cs, q0);   % 1×N，仍然在 [0, q0-1]
    
    %% 4. 将 t 中的 “非中心化”系数映射回整数区间 [ -floor(q0/2) .. floor((q0-1)/2 ) ]
    %    若 t(i) > q0/2，则置 t(i) = t(i) - q0
    half_q0 = floor(q0/2);
    m_int = t;  % 复制一份
    idx = m_int > half_q0;
    m_int(idx) = m_int(idx) - q0;
    % m_int 中即为整数系数，可能为负

    %% 5. 输出 m_int（长度 N 的整数行向量）
end