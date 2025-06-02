function ct_B = embedding_ct_into_slot(B, c)
    [d, n] = size(B);
    if mod(d, n) ~= 0
        error('d 必须是 n 的整数倍');
    end
    c_temp = c; % 保存c的值, 以防c < 1时数值丢失
    if c >= 1
        if mod(n, c) ~= 0
            error('当c>=1时，n必须是c的整数倍');
        end
    else
        if mod(1/c, 1) ~= 0 || mod(d, 1/c) ~= 0
            error('当c<1时，1/c必须是整数且d必须是1/c的整数倍');
        end
        c = 1;% c < 1时按照1执行计算
    end
    
    % 1. 分块
        num_blocks = d / n;
        blocks = cell(num_blocks, 1);
        for j = 1:num_blocks
            row_start = (j-1)*n + 1;
            row_end = j*n;
            blocks{j} = B(row_start:row_end, :);
        end
    
    % 2. 提取循环对角线
        L = cell(num_blocks, n);
        for j = 1:num_blocks
            block = blocks{j};
            for i = 0:n-1
                L{j, i+1} = get_circular_diag(block, i);
            end
        end
    % 3. 产生密文槽
        num_ct = n / c;
        % 使用与输入相同类型的零矩阵初始化
        temp_ct_B = zeros(num_ct, c*d, 'like', B);
        for m = 1:num_ct
            % 当前组的对角线范围
            diag_start = (m-1)*c + 1;
            diag_end = m*c;
            
            % 先按 L0, L1, ..., Lc-1 的顺序拼接所有位置
            concatenated = [];
            for k = diag_start:diag_end
                for pos = 1:n
                    for j = 1:num_blocks
                        concatenated = [concatenated, L{j, k}(pos)];
                    end
                end
            end
            temp_ct_B(m,:) = concatenated;
        end
    % 4. c < 1时的后处理
        c = c_temp;
        if c < 1
            num_chunks = 1/c; % 分成多少个大块
            chunk_size = c * d; % 每个大块的列数
            % 将 temp_ct_B 按列拆分成 num_chunks 个大块
            chunks = mat2cell(temp_ct_B, size(temp_ct_B, 1), repmat(chunk_size, 1, num_chunks));
            % 垂直堆叠这些大块
            ct_B = vertcat(chunks{:});
        else
            ct_B = temp_ct_B;
        end
end

function vec = get_circular_diag(block, i)
    % block: n×n 矩阵
    % i: 下对角线编号（0=主对角线，1=第一条循环下对角线）
    [n, ~] = size(block);
    % 使用与输入相同类型的零矩阵初始化
    vec = zeros(1, n, 'like', block);
    for k = 1:n
        row = mod(k + i - 1, n) + 1; % 循环行索引
        col = k;                     % 固定列索引
        vec(k) = block(row, col);
    end
end
