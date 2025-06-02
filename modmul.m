function r = modmul(a, b, q)
%==========================================================================
% 支持向量/矩阵输入的 (a .* b) mod q
%
%   r = modmul(a, b, q)
%
% 输入：
%   - a, b : 标量 或 数组（向量/矩阵），要求 a 和 b 形状相同，
%            或者其中一个是标量，另一个是数组。函数会对对应元素
%            逐元素计算 (a_i .* b_i) mod q。
%   - q    : 标量模数（假设 q < 2^40 且为正整数）
%
% 输出：
%   - r    : 与 a、b 相同形状的数组，每个元素等于 (a(i).*b(i)) mod q。
%
% 要点：
%   * 先把 a、b 直接 cast 为 uint64，再对 q 做一次 mod（两者同类型）。
%   * L = 26 位拆分：将 a、b 拆成高 26 位 + 低 26 位，分别计算分块乘积，
%     每一步都保持在 <2^52，最后再合并做一次 mod q。
%   * 最终输出默认为 double 类型（若需要 uint64，见备注）。
%==========================================================================
    %—— 1. 先把 a、b 都转为 uint64，再做 mod q
    a = uint64(a);                  % 无论 a 是标量还是向量，都先转 uint64
    b = uint64(b);
    a = mod(a, uint64(q));          % 此时 (uint64, uint64) 调用 mod，合法
    b = mod(b, uint64(q));

    %—— 2. 定义拆分长度 L = 26，掩码 MASK 用来取得“低 26 位”
    L    = 26;
    MASK = bitshift(uint64(1), L) - uint64(1);  % = 2^26 - 1

    %—— 3. 按元素拆分 a, b
    %    a_hi = floor(a/2^26)，a_lo = a mod 2^26；同理 b_hi, b_lo
    a_hi = bitshift(a, -L);     % 和 a 保持同形状
    a_lo = bitand(a, MASK);
    b_hi = bitshift(b, -L);
    b_lo = bitand(b, MASK);

    %—— 4. 预计算 2^26 mod q 和 2^52 mod q（都是 uint64 标量）
    M26_mod = uint64( mod( bitshift(uint64(1), L),   uint64(q) ) );   % = 2^26 mod q
    M52_mod = uint64( mod( bitshift(uint64(1), 2*L), uint64(q) ) );   % = 2^52 mod q

    %—— 5. 计算三大块（所有运算均为逐元素，且都在 uint64 范围内）
    % 5.1) tmp0 = a_hi .* b_hi mod q，再乘 (2^52 mod q)，再 mod q
    t0    = mod( a_hi .* b_hi, uint64(q) );                % < q
    part0 = mod( t0 .* M52_mod,      uint64(q) );          % 每个元素 < q

    % 5.2) c = a_hi.*b_lo + a_lo.*b_hi，先 mod q，再拆 26 位
    c     = mod( a_hi .* b_lo + a_lo .* b_hi, uint64(q) );  % < q
    c_hi  = bitshift(c, -L);     % < 2^(40-26) = 2^14
    c_lo  = bitand(c, MASK);     % < 2^26
    %     高 26 位乘 2^52 mod q，低 26 位乘 2^26 mod q，再相加 mod q
    part1 = mod( c_hi .* M52_mod + c_lo .* M26_mod, uint64(q) );

    % 5.3) t2 = a_lo .* b_lo mod q
    part2 = mod( a_lo .* b_lo, uint64(q) );

    %—— 6. 三部分累加，再做一次 mod q
    %    part0, part1, part2 每个元素都 < q，和 < 3q < 2^42，安全
    r_uint64 = mod( part0 + part1 + part2, uint64(q) );

    %—— 7. 转回 double 输出（若需要 uint64，可直接返回 r_uint64）
    r = double(r_uint64);
end
