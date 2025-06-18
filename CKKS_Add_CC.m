function ct_out = CKKS_Add_CC(params, ct1, ct2)
% CKKS_Add  同态加法：对两个同层级 ciphertext 做逐层逐系数相加
%
% 输入：
%   params — CKKS 参数结构体，包含至少 .Q_primes = [q0, q1, …, qL]
%   ct1    — 第一个 ciphertext，结构如下：
%            ct1.c0 = { c0^{(0)}, c0^{(1)}, …, c0^{(L)} }
%            ct1.c1 = { c1^{(0)}, c1^{(1)}, …, c1^{(L)} }
%            每个 c{i}^{(j)} 都是一个 1×N 的行向量，对应多项式系数 (mod q_j)。
%   ct2    — 第二个 ciphertext，结构完全同 ct1。两者必须是同一层级（同样的层数 L）。
%
% 输出：
%   ct_out — 输出 ciphertext，层级与输入相同，且在每个层级 j 将对应多项式系数取模 q_j 后相加：
%            ct_out.c0{j+1} = ( ct1.c0{j+1} + ct2.c0{j+1} ) mod q_j
%            ct_out.c1{j+1} = ( ct1.c1{j+1} + ct2.c1{j+1} ) mod q_j
%            其中 j = 0..L，MATLAB 下标 = j+1。
%
% 例外检查：
%   1) ct1 与 ct2 的层级（cell 数量）必须一致，否则报错。
%   2) 每个层级的多项式长度 N 必须一致（如果需要也可检查）。
%
    % 提取 Q_primes，以及层数 L
    Q_primes = params.Q_primes;            % 1×(L+1) 向量
    M = numel(ct1.c0);            % M = l + 1

    % 校验：ct1 和 ct2 必须有相同的层数
    if numel(ct1.c0) ~= M || numel(ct2.c0) ~= M ...
       || numel(ct1.c1) ~= M || numel(ct2.c1) ~= M
        error('CKKS_Add: 输入 ciphertext 层数不一致，无法相加。');
    end

    % 预分配输出结构
    ct_out.c0 = cell(1, M);
    ct_out.c1 = cell(1, M);

    % 遍历每个层级 j = 0..L 对应 MATLAB 下标 j+1
    for j = 1 : M
        qj = Q_primes(j);  

        % 取出第 j−1 层输入 ciphertext 的多项式系数
        a0 = ct1.c0{j};   % 1×N 向量 (mod qj)
        a1 = ct1.c1{j};   
        b0 = ct2.c0{j};   % 1×N 向量 (mod qj)
        b1 = ct2.c1{j};

        % 校验向量长度是否相同
        if numel(a0) ~= numel(b0) || numel(a1) ~= numel(b1)
            error('CKKS_Add: 第 %d 层多项式长度不一致。', j-1);
        end

        % 逐系数相加并做模 qj
        % 由于 a0, b0 可能已是 double 或 uint64，使用 double 再 mod 也行
        c0 = mod( double(a0) + double(b0), double(qj) );
        c1 = mod( double(a1) + double(b1), double(qj) );

        % 存回输出结构
        ct_out.c0{j} = c0;
        ct_out.c1{j} = c1;
    end
end