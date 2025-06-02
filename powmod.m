function r = powmod(a, e, q)
    % 二分幂：计算 a^e mod q
    r = uint64(1);
    base = mod(uint64(a), uint64(q));
    expv = uint64(e);
    while expv > 0
        if bitand(expv, 1) == 1
            r = modmul(r, base, q);
        end
        base = modmul(base, base, q);
        expv = bitshift(expv, -1);
    end
    r = double(r);
end