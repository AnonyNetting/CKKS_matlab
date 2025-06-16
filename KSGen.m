%% KSGen: 生成 RNS-CKKS 的密钥切换密钥
function [ks_keys, ks_map] = KSGen(params, sk, rotation_steps)
    N = params.N;
    k = params.k;
    P_primes = params.P_primes;
    Q_primes = params.Q_primes;
    M = numel(Q_primes);

    rotation_steps = unique(rotation_steps);
    t = numel(rotation_steps);
    ks_keys = cell(1, t);
    ks_map  = containers.Map('KeyType','double','ValueType','double');

    % 统一噪声
    e_glob = params.error_dist();
    % RNS 模数列表 D
    D = [P_primes, Q_primes];

    for idx = 1:t
        r = rotation_steps(idx);
        % 计算 s2：r==0 时执行重线性化 (s^2)，否则执行旋转
        if r == 0
            % SWK: schoolbook 计算 s^2 mod (X^N+1)
            conv_len = 2*N - 1;
            tmp = fliplr( conv(fliplr(sk), fliplr(sk)) );  % 长度 2N-1
            % 模 X^N + 1：对于 i=N+1:2N-1，将系数的相反数累加到 tmp(i-N) 
            s_square = zeros(1, N);
            for i = 1:conv_len
                if i <= N
                    s_square(i) = tmp(i);
                else
                    j = i - N;
                    s_square(j) = s_square(j) - tmp(i);
                end
            end
            s2 = s_square;
        else
            % 旋转: g^r mod 2N
            gal = uint64(5);
            gr  = powmod(gal, r, 2*N);
            s2 = zeros(1, N);
            for i = 1:N
                tgt = mod(uint64(i-1)*gr, 2*N);
                if tgt < N
                    s2(tgt+1) = sk(i);
                else
                    s2(tgt-N+1) = -sk(i);
                end
            end
        end

        % 为每个模生成切换对 (a', b')
        m = numel(D);
        a_all = cell(1, m);
        b_all = cell(1, m);
        for i_pb = 1:k
            p_i = P_primes(i_pb);
            a = randi([0, p_i-1], 1, N);
            c = ckks_poly_mult_mod(a, mod(sk, p_i), p_i, N, params.psi{i_pb});
            b = mod(-c + mod(e_glob, p_i), p_i);
            a_all{i_pb} = a;
            b_all{i_pb} = b;
        end
        for j_q = 1:M
            q_j = Q_primes(j_q);
            D_idx = k + j_q;
            a = randi([0, q_j-1],1,N);
            c = ckks_poly_mult_mod(a, mod(sk, q_j), q_j, N, params.psi{D_idx});
            Pm = uint64(1);
            for ii = 1:k
                Pm = modmul(Pm, mod(P_primes(ii), q_j), q_j);
            end
            term = modmul(Pm, mod(s2, q_j), q_j);
            b = mod(-c + term + mod(e_glob, q_j), q_j);
            a_all{D_idx} = a;
            b_all{D_idx} = b;
        end

        ks_keys{idx} = struct('power', r, 'a', {a_all}, 'b', {b_all});
        ks_map(r)    = idx;
    end
end