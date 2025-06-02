function B = inverse_embedding_from_slot(ct_B, n, c)
    % 处理 c < 1 的情况：先恢复垂直堆叠前的结构
    c_temp = c;
    if c < 1
        num_chunks = 1/c;
        chunk_rows = size(ct_B, 1) / num_chunks;
        chunk_cols = size(ct_B, 2);
        % 水平拼接分块
        chunks = mat2cell(ct_B, repmat(chunk_rows, 1, num_chunks), chunk_cols);
        ct_B = horzcat(chunks{:});
        c = 1; % 后续按 c=1 处理
    end
    
    [num_ct, slot_size] = size(ct_B);
    d = slot_size / c;
    if mod(d, n) ~= 0
        error('d 必须是 n 的整数倍');
    end
    
    num_blocks = d / n;
    B = zeros(d, n, 'like', ct_B);
    
    % 初始化存储所有 Li(Bj) 的 cell 数组
    L = cell(num_blocks, n);
    
    % 从 ct_B 中恢复 Li(Bj)
    for m = 1:num_ct
        diag_start = (m-1)*c + 1;
        diag_end = m*c;
        
        ct_Bm = ct_B(m, :);
        elements_per_diag = n * num_blocks;
        
        for k = diag_start:diag_end
            for pos = 1:n
                for j = 1:num_blocks
                    abs_pos = (k-diag_start)*elements_per_diag + (pos-1)*num_blocks + j;
                    if isempty(L{j,k})
                        L{j,k} = ct_Bm(abs_pos);
                    else
                        L{j,k} = [L{j,k}, ct_Bm(abs_pos)];
                    end
                end
            end
        end
    end
    
    % 从 Li(Bj) 重建原始块
    for j = 1:num_blocks
        block = zeros(n, n, 'like', ct_B);
        for i = 0:n-1
            vec = L{j, i+1};
            for k = 1:n
                row = mod(k + i - 1, n) + 1;
                col = k;
                block(row, col) = vec(k);
            end
        end
        row_start = (j-1)*n + 1;
        row_end = j*n;
        B(row_start:row_end, :) = block;
    end
    
    % 恢复原始 c 值（仅用于外部一致性，不影响计算）
    c = c_temp;
end
