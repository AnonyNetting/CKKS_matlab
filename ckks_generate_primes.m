function primes_list = ckks_generate_primes(num_modulus, N, modulus_bit_length, range)
    % 计算中心点 (2^bit)
    center = pow2(modulus_bit_length);
    
    % 计算实际搜索范围 [min_value, max_value]
    min_value = max(1, round((1 - range) * center)); % 素数至少为2
    max_value = min(pow2(modulus_bit_length+1) - 1, round((1 + range) * center));
    
    % 确保范围有效
    if min_value >= max_value
        error('Invalid range: no numbers in [%d, %d]', min_value, max_value);
    end
    
    step = 2 * N;  % 同余条件步长
    primes_list = [];
    
    % 找到中心点附近满足同余条件的起始点
    start_point = center - mod(center - 1, step);
    if start_point < min_value
        start_point = start_point + step * ceil((min_value - start_point)/step);
    end
    
    % 初始化双向指针 (从中心对称位置开始)
    current_down = start_point;
    current_up = start_point + step;
    
    % 计数器确保双向交替搜索
    down_count = 0;
    up_count = 0;
    
    % 双向对称搜索素数
    while length(primes_list) < num_modulus
        % 对称向下搜索 (优先搜索接近中心的点)
        if down_count <= up_count && current_down >= min_value
            if current_down <= max_value && is_prime(current_down)
                primes_list(end+1) = current_down;
                if length(primes_list) == num_modulus
                    break;
                end
            end
            current_down = current_down - step;
            down_count = down_count + 1;
        end
        
        % 对称向上搜索
        if up_count <= down_count && current_up <= max_value
            if current_up >= min_value && is_prime(current_up)
                primes_list(end+1) = current_up;
                if length(primes_list) == num_modulus
                    break;
                end
            end
            current_up = current_up + step;
            up_count = up_count + 1;
        end
        
        % 边界检查
        if current_down < min_value && current_up > max_value
            break;
        end
    end
    
    % 验证是否找到足够素数
    if length(primes_list) < num_modulus
        error('Found only %d primes (requested %d)', length(primes_list), num_modulus);
    end
    
    % 按距离中心点的远近排序
    [~, idx] = sort(abs(primes_list - center));
    primes_list = primes_list(idx);
end

% 高效素数检测函数
function result = is_prime(n)
    % 处理小整数
    if n < 2
        result = false;
        return;
    end
    
    % 小整数使用内置函数
    if n < 2^50
        result = isprime(uint64(n));
        return;
    end
    
    % 大整数使用 Miller-Rabin 测试
    result = true;
    test_bases = [2, 325, 9375, 28178, 450775, 9780504, 1795265022];
    
    % 特殊情况处理
    if mod(n, 2) == 0
        result = (n == 2);
        return;
    end
    
    % 分解 n-1 = 2^s * d
    d = n - 1;
    s = 0;
    while mod(d, 2) == 0
        d = d / 2;
        s = s + 1;
    end
    
    % 对每个基数进行测试
    for a = test_bases
        if a >= n
            continue;
        end
        
        x = powermod(a, d, n);
        if x == 1 || x == n-1
            continue;
        end
        
        composite = true;
        for r = 1:s-1
            x = mod(x * x, n);
            if x == n-1
                composite = false;
                break;
            end
        end
        
        if composite
            result = false;
            return;
        end
    end
end

% 高效幂模运算
function result = powermod(base, exponent, modulus)
    result = 1;
    base = mod(base, modulus);
    
    while exponent > 0
        if bitand(exponent, 1) == 1
            result = mod(result * base, modulus);
        end
        exponent = bitshift(exponent, -1); % exponent = floor(exponent/2)
        base = mod(base * base, modulus);
    end
end