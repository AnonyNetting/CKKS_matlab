function inv = modinv(a, m)
    % 扩展欧几里得求逆：计算 a^{-1} mod m （假设 gcd(a,m)=1）
    %
    % 输入：
    %   a, m — 可以是标量或相同形状的数组（但本例主要用在标量场景）。
    %          函数内部会先做 mod，确保 0 ≤ a < m，且 a、m 都为 uint64。
    % 输出：
    %   inv — 与输入 a 同形状的数组（或标量），等于 (a .^(-1)) mod m，
    %         如果 a 与 m 不互素则报错。
    %
    % 示例：
    %   q = 1099511480321;
    %   inv256 = modinv(256, q);   % 计算 256 在 mod q 意义下的逆元
    %

    %—— 1) 先把 a、m 转为 uint64 并做 mod
    a = mod(uint64(a), uint64(m));  % 确保 0 ≤ a < m，类型为 uint64
    m = uint64(m);                  % m 也转为 uint64
    if any(a(:) == 0)
        error('modinv: a and m must be coprime (a ≠ 0).');
    end

    %—— 2) 转为 int64，以支持中间出现负数
    ai = int64(a);
    mi = int64(m);

    % 如果输入是数组，逐元素循环（可改为 arrayfun，但这里示范用循环更直观）
    sz = size(ai);
    inv = zeros(sz);      % 预分配，最后再转成 double
    for idx = 1:numel(ai)
        x = ai(idx);
        M = mi(idx);

        % 扩展欧几里得初始化
        r0 = M;   r1 = x;
        t0 = int64(0);
        t1 = int64(1);

        %—— 3) 循环求 gcd 并记录 t0, t1
        while r1 ~= 0
            % 用 idivide 做 floor 整除，保证全程整数
            q = idivide(r0, r1, 'floor');  
            % 更新 (r0, r1)
            tmp_r = r0 - q * r1;
            r0 = r1;
            r1 = tmp_r;
            % 更新 (t0, t1)
            tmp_t = t0 - q * t1;
            t0 = t1;
            t1 = tmp_t;
        end

        if r0 ~= 1
            error('modinv: inverse does not exist for a = %d, m = %d', x, M);
        end

        %—— 4) t0 已经是 x 的 GCD=1 情况下的负/正系数，取 mod M
        inv_i = mod(t0, M);    % 这一步结果仍为 int64，且 ∈ [0, M)
        inv(idx) = double(inv_i);  % 最终返回 double。如需 uint64 可改为 uint64(inv_i)
    end
end
