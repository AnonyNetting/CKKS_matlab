function ct_rot = CKKS_Rotation(params, ct, r)
% CKKS_Rotation: 对密文 ct 执行 Galois 旋转（坐标置换+符号翻转），不做 KeySwitch
%
% 输入：
%   params: 系统参数，包含 N, k, P_primes, Q_primes, psi 等
%   ct:     原始密文，ct.c0, ct.c1 是 1×M cell，每个 cell 是 1×N 向量
%   r:      旋转步长 (0<r<N)
%
% 输出：
%   ct_rot: 旋转后密文，结构同 ct
%
    N = params.N;
    %if r <= 0 || r >= N
    if r < 0 || r >= N/2
        error('CKKS_Rotation: r 必须满足 0<r<N/2');
    end
    
    % 1. 计算 g = 5^r mod 2N
    gal_gen = uint64(5);
    g = powmod(gal_gen, r, 2*N);
    
    % 2. 当前 RNS 深度
    M = numel(ct.c0);
    c0_rot = cell(1, M);
    c1_rot = cell(1, M);
    
    % 3. 对每个层做系数置换+符号翻转
    for j = 1 : M
        poly0 = ct.c0{j};  % 1×N
        poly1 = ct.c1{j};  % 1×N
        qj = params.Q_primes(j);
        tmp0 = zeros(1, N);
        tmp1 = zeros(1, N);
        for i = 1 : N
            % —— 检查 poly0(i), poly1(i) 是否真的在 [0,qj)
            if poly0(i) < 0 || poly0(i) >= qj || poly1(i) < 0 || poly1(i) >= qj
                warning('第 %d 层中，第 %d 个系数 poly0(i)=%d / poly1(i)=%d 不在 [0,%d)', ...
                        j-1, i, poly0(i), poly1(i), qj);
            end

            exp_i = uint64(i - 1);
            tgt = modmul(exp_i, g, 2*N);
            if tgt < N
                % 若 tgt < N, 直接搬移
                tmp0(double(tgt)+1) = poly0(i);
                tmp1(double(tgt)+1) = poly1(i);

                % —— 搬运后再检查 tmp0(tgt+1) 是否在 [0,qj)
                if tmp0(tgt+1) < 0 || tmp0(tgt+1) >= qj
                    warning('第 %d 层中，旋转后 tmp0(%d)=%d 不在 [0,%d)', ...
                            j-1, tgt+1, tmp0(tgt+1), qj);
                end
            else
                % tgt ∈ [N,2N-1], 搬移并乘 -1
                idx_new = double(tgt - N) + 1;
                tmp0(idx_new) = mod(-poly0(i), qj);
                tmp1(idx_new) = mod(-poly1(i), qj);

                % —— 搬运后再检查 tmp0(idx_new) 是否在 [0,qj)
                if tmp1(idx_new) < 0 || tmp1(idx_new) >= qj
                    warning('第 %d 层中，旋转后 tmp1(%d)=%d 不在 [0,%d)', ...
                            j-1, idx_new, tmp1(idx_new), qj);
                end
            end

        end
        c0_rot{j} = tmp0;
        c1_rot{j} = tmp1;
    end
    
    ct_rot.c0 = c0_rot;
    ct_rot.c1 = c1_rot;
end