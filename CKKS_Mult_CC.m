%% CKKS_Mult_CC: 同态密文×密文乘法 + 重线性化
function ct_out = CKKS_Mult_CC(params, ct1, ct2)
    % 确保两个密文深度一致
    L1 = numel(ct1.c0);
    L2 = numel(ct2.c0);
    if L1 ~= L2
        error('CKKS_Mult_CC: 两密文深度不一致 (L1=%d, L2=%d)。', L1, L2);
    end
    M = L1;
    N = params.N;
    k = params.k;
    Q = params.Q_primes;

    % 1. 在深度 0..L 上分别计算 d0, d1, d2
    ct_out.d0 = cell(1, M);
    ct_out.d1 = cell(1, M);
    ct_out.d2 = cell(1, M);
    for j = 1:M
        qj  = Q(j);
        psi = params.psi{k + j};

        ct_out.d0{j} = ckks_poly_mult_mod(ct1.c0{j}, ct2.c0{j}, qj, N, psi);
        t1      = ckks_poly_mult_mod(ct1.c0{j}, ct2.c1{j}, qj, N, psi);
        t2      = ckks_poly_mult_mod(ct1.c1{j}, ct2.c0{j}, qj, N, psi);
        ct_out.d1{j} = modadd(t1, t2, qj);
        ct_out.d2{j} = ckks_poly_mult_mod(ct1.c1{j}, ct2.c1{j}, qj, N, psi);
    end
end