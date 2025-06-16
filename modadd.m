function c = modadd(a, b, q)
%==========================================================================
% 支持向量/矩阵输入的 (a + b) mod q
%
%   r = modmul(a, b, q)
%
% 输入：
%   - a, b : 标量 或 数组（向量/矩阵），要求 a 和 b 形状相同，函数会对对应元素
%            逐元素计算 (a_i + b_i) mod q。
%   - q    : 标量模数
%
% 输出：
%   - r    : 与 a、b 相同形状的数组，每个元素等于 (a(i) + b(i)) mod q。
%
% 要点：
%   * 先把 a、b 直接 cast 为 uint64，再对 q 做一次 mod（两者同类型）。
%     a + b最后再合并做一次 mod q。
%   * 最终输出默认为 double 类型
%==========================================================================

    % 检查 a 和 b
    if ~isequal(size(a), size(b))
        error('输入 a 和 b 必须具有相同的尺寸。');
    end
    if min(a) < 0 || min(b) < 0 || max(a) >= q || max(b) >= q
        error("illegal input");
    end
    % 将 q 转为 uint64
    q_uint = uint64(q);

    % 将 a 和 b 转为 uint64
    a_uint = mod(uint64(a), q);
    b_uint = mod(uint64(b), q);

    % 先计算 uint64 类型下的加法
    sum_uint = a_uint + b_uint;

    % 最后对 (a + b) 做 mod q
    c_uint = mod(sum_uint, q_uint);

    % 转回 double 类型输出
    c = double(c_uint);
end