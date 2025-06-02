function rev = bitrev(i, bits)
    % 将 0 ≤ i < 2^bits 的整数 i 做 "bits" 位的二进制位反转
    rev = 0;
    for k = 1:bits
        rev = bitor(bitshift(rev, 1), bitand(i, 1));
        i = bitshift(i, -1);
    end
end