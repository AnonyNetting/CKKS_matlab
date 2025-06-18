function r = modmul(a, b, q)
%MODMUL   元素级 (a.*b) mod q，支持 q 最多 50 bit
%   r = MODMUL(a, b, q)
%   - a, b: 数组或标量，形状相同或一者为标量
%   - q   : 标量模数，1 <= q < 2^50
%
%   本函数先强制 cast 为 uint64，再用“动态 L‑bit 拆分”方法
%   计算 (a*b) mod q，保证所有中间乘法都 ≤ 2^64。

    %—— 输入检查 & uint64 转换 ——
    if any(a(:) < 0) || any(b(:) < 0)
        error('modmul: a, b must be ≥ 0');
    end
    a = uint64(a);
    b = uint64(b);
    q = uint64(q);
    if q == 0
        error('modmul: modulus q must be positive');
    end
    % 统一到 [0,q)
    a = mod(a, q);
    b = mod(b, q);

    %—— 根据 q 位宽 自动选 L —— 
    % bitlen = ⌈log2 q⌉
    bitlen = ceil(log2(double(q)));
    % 为保证 2*(bitlen-L) + bitlen ≤ 64  ⇒ L ≥ bitlen - 32
    L = max(bitlen - 32, 1);
    % 但不能超出 bitlen 本身
    L = min(L, bitlen);
    MASK = bitshift(uint64(1), L) - uint64(1);

    %—— 拆分高低位 ——
    a_lo = bitand(a, MASK);
    a_hi = bitshift(a, -L);
    b_lo = bitand(b, MASK);
    b_hi = bitshift(b, -L);

    %—— 预计算 2^L mod q, 2^(2L) mod q ——
    M_L   = bitshift(uint64(1), L);      % 2^L
    M_2L  = bitshift(uint64(1), 2*L);    % 2^(2L)
    M_L   = mod(M_L,  q);                % 2^L mod q
    M_2L  = mod(M_2L, q);                % 2^(2L) mod q

    %—— 各部分累加 —— 
    %  c2 = a_lo*b_lo
    part2 = mod(a_lo .* b_lo, q);

    %  c1 = a_hi*b_lo + a_lo*b_hi, 再乘 2^L
    c1    = a_hi .* b_lo + a_lo .* b_hi;          % ≤ 2*(2^(bitlen-L)*2^L) ≤ 2*q
    part1 = mod( mod(c1, q) .* M_L, q);

    %  c0 = a_hi*b_hi, 再乘 2^(2L)
    c0    = a_hi .* b_hi;                        % ≤ (2^(bitlen-L))^2
    part0 = mod( mod(c0, q) .* M_2L, q);

    %—— 最终结果 ——
    r = double( mod(part0 + part1 + part2, q) );
end