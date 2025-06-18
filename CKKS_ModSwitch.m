function ct_out = CKKS_ModSwitch(ct_in, params)
%CKKS_MODSWITCH  CKKS Modulus Switching（使深度掉落1层）
%   输入 ct_in: 结构体，fields .c0, .c1，均为 1×(L+1) cell  
%         params: 包含 Q_primes 的参数结构
%   输出 ct_out: 结构体，.c0, .c1 均为 1×L cell

    m = ckks_encode(ones(1, params.N/2), 2^params.qother_bit_length);

    ct_out = CKKS_Mult_PC(params, ct_in, m);

    ct_out = CKKS_Rescale(ct_out, params);
end