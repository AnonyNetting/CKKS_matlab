%% CKKS_KeyGen.m
% 生成 RNS-CKKS 所需的公/私钥，并包含对密文乘法结果（c2）所需的 relinearization key (SWK) 的切换密钥
% 通过 rotation_steps == 0 识别 SWK（即密文乘法密钥切换）

function keys = CKKS_KeyGen(params, rotation_list)
    % params: CKKS_Setup 返回的系统参数
    % rotation_list: 用户指定的旋转步长向量（如 [1, -1, ...]），不含 0
    
    %% 1. 基础密钥生成
    N = params.N;
    L = params.L;
    k = params.k;
    Q_primes = params.Q_primes;
    
    % 1.1 私钥 s
    s = params.key_dist();         % 1xN, {0,±1}
    s = reshape(s,1,[]);

    % 1.2 公钥 (a_j, b_j)
    pk.a = cell(1, L+1);
    pk.b = cell(1, L+1);
    e = params.error_dist();       % 同一噪声可用于所有 q_j，也可重采样
    for j = 1:L+1
        qj = Q_primes(j);
        a_j = randi([0, qj-1], 1, N);
        s_q = mod(s, qj);
        c = ckks_poly_mult_mod(a_j, s_q, qj, N, params.psi{k+j});
        b_j = mod(-c + e, qj);
        pk.a{j} = a_j;
        pk.b{j} = b_j;
    end

    %% 2. 组合私/公钥
    keys.sk = s;
    keys.pk = pk;

    %% 3. 生成切换密钥：包括 SWK（r=0）和旋转切换（r∈rotation_list）
    % rotation_steps: 0 表示密文乘法 relinearization，其他为旋转
    rotation_steps = unique([0, rotation_list(:)']);
    [swk_and_rot, map] = KSGen(params, s, rotation_steps);
    keys.swk  = swk_and_rot;
    keys.map = map;
end