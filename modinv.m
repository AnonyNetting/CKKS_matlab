function inv = modinv(a_in, m)
%MODINV   计算 a_in 在模 m 下的乘法逆元，内部用 modmul 避免溢出
%   inv = MODINV(a_in, m)
%   - a_in, m: 非负整数标量，0 < a_in < m < 2^50
%   返回 inv (uint64)，满足 modmul(a_in, inv, m) == 1。
%
%   如果 gcd(a_in, m) ≠ 1，则报错。

    %% 1. 输入检查
    if ~(isscalar(a_in) && isscalar(m) && a_in >= 0 && m > 0)
    error('modinv_safe: 输入必须是标量，且 0 ≤ a < m, m > 0');
    end
    % 转成 uint64 并规约
    m    = uint64(m);
    a    = mod(uint64(a_in), m);
    if a == 0
    error('modinv_safe: a 在模 m 下为 0，无逆元。');
    end
    
    %% 2. 初始化：我们跟踪
    %   aa, bb 做整数除余
    %   x0, x1 始终是对应系数 mod m
    aa = a;  
    bb = m;
    x0 = uint64(1);    % 对应 aa 初始
    x1 = uint64(0);    % 对应 bb 初始
    
    %% 3. 主循环：当 bb ≠ 0 时
    while bb ~= 0
    % 整数商与余数
    q = idivide(aa, bb, 'floor');      % floor(aa / bb), uint64
    r = mod(aa, bb);                   % aa - q*bb
    
    % 更新系数： x2 = x0 - q*x1  (mod m)，用 modmul 做乘法
    %    t = modmul(x1, q, m)   得到 (q*x1) mod m
    t  = modmul(x1, q, m);
    if x0 >= t
      x2 = x0 - t;
    else
      % 防止负数饱和：加上一个模再减
      x2 = x0 + m - t;
    end
    x2 = mod(x2, m);  % 确保落在 [0,m)
    
    % 向前推进 (aa,bb) 和 (x0,x1)
    aa = bb;
    bb = r;
    x0 = x1;
    x1 = x2;
    end
    
    %% 4. 检查 gcd
    if aa ~= 1
    error('modinv_safe: a 与 m 不互素，无逆元。');
    end
    
    %% 5. 最终 x0 就是逆元
    inv = x0;

    if modmul(inv, a_in, m) ~= 1
        error("modinv error");
    end
end