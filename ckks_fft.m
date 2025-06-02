function X = ckks_fft(x)
    N = length(x);
    N2 = 2 * N;

    % (1) 先把 x 填充到长度 2N 的向量 Z, 但只在前 N 位填数据，后 N 位填 0：
    Z = zeros(1, N2);
    Z(1:N) = x;        % Z(N+1:2N) 自动 = 0

    % (2) 直接调用 2N 点 fft：
    F = fft(Z);        % F[k] = Σ_{n=0..2N-1} Z[n] e^{−2πi n k/(2N)}

    % (3) 取出那些“奇数下标”对应的频域值：
    %      X[j] = F[2j+1] （0-based），但 MATLAB 里索引要 +1，所以写成 F(2*j+2)
    X = zeros(N,1);
    for j = 0:(N-1)
        X(j+1) = F(2*j + 2);
    end
end