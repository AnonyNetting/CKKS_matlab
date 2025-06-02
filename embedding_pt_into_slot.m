function pt_A = embedding_pt_into_slot(A, n, c)%先遍历第一个new_blocks中的U(0),直到U(n-1)；然后遍历第二个new_blocks...
    % 检查输入矩阵的维度是否合理
    [d, ~] = size(A);
    if mod(d, n) ~= 0
        error('d 必须是 n 的整数倍');
    end
    c_temp = c; % 保存c的值, 以防c < 1时数值丢失
    if c >= 1
        if mod(n, c) ~= 0
            error('当c>=1时，n必须是c的整数倍');
        end
    else
        if mod(1/c, 1) ~= 0 || mod(d, 1/c) ~= 0 || mod(c*d, d/n) ~= 0
            error('当c<1时，1/c必须是整数, d必须是1/c的整数倍, cd必须是d/n的整数倍(即nc >= 1');
        end
        c = 1;% c < 1时按照1执行计算
    end
    
    % 计算划分后的块数
    num_blocks = d / n;
    
    % 1. 将A划分成d/n×d/n个n×n矩阵，新的d/n×d/n矩阵是blocks
    blocks = cell(num_blocks, num_blocks);
    
    for i = 1:num_blocks
        for j = 1:num_blocks
            row_start = (i-1)*n + 1;
            row_end = i*n;
            col_start = (j-1)*n + 1;
            col_end = j*n;
            
            blocks{i, j} = A(row_start:row_end, col_start:col_end);
        end
    end

    % 2. 将blocks合并成为d/n×1维矩阵new_blocks，new_blocks的每个元素是n×d维矩阵
    new_blocks = cell(num_blocks, 1);
    
    % 遍历每个对角线
    for k = 0:num_blocks-1
        for j = 0:num_blocks-1
            new_blocks{k+1} = horzcat(new_blocks{k+1}, blocks{mod(k+j,d/n)+1, j+1});
        end
    end
    
    % 3. 修改后的提取new_blocks产生final_vectors
    % 初始化 final_vectors（d×d），使用与输入相同类型的零矩阵
    final_vectors = zeros(d, d, 'like', A);
    
    % 新的遍历顺序：先遍历所有block，再遍历每个block内的diag_idx
    k = 1;
    for block_idx = 1:num_blocks
        for diag_idx = 1:n
            % 提取 new_blocks{block_idx}
            current_block = new_blocks{block_idx}; % n×d 矩阵
            
            % 将 current_block 拆分为 (d/n) 个 n×n 块
            num_subblocks = d / n;
            subblocks = mat2cell(current_block, n, n*ones(1, num_subblocks));
            
            % 提取每个 n×n 块的第 diag_idx 条下对角线
            vec = [];
            for m = 1:n
                for i = 1:num_subblocks
                    vec = [vec, subblocks{i}(m, mod(diag_idx-1+m-1, n)+1)];
                end
            end
            
            % 存储到 final_vectors
            final_vectors(k, :) = vec;
            k = k + 1;
        end
    end
    
    % 4. 为每个final_vectors添加额外的旋转量,最终生成pt_A
    % 初始化 pt_A(dn/c×cd)，使用与输入相同类型的零矩阵
    temp_pt_A = zeros(d*n/c, c*d, 'like', A);
    
    % 计算每个n×d块的数量
    num_row_blocks = d / n;        % 中块数量 = d/n
    num_small_per_medium = n / c;  % 每个中块内的小块数量 = n/c
    
    for i = 0:c-1          % 外层循环：c次 (i ∈ [0,c-1])
        for j = 0:n/c-1     % 内层循环：n/c次 (j ∈ [0,n/c-1])
            % 创建临时矩阵存储旋转后的结果
            rotated_block = zeros(d, d, 'like', A);
            
            % ===== 第一步：行旋转 =====
            % 对每个n×d块单独进行行旋转（向下循环移动j*c行）
            for block = 0:num_row_blocks-1
                row_start = block*n + 1;
                row_end = (block+1)*n;
                rotated_block(row_start:row_end, :) = ...
                    circshift(final_vectors(row_start:row_end, :), [j*c, 0]);
            end
            
            % ===== 第二步：中块内小块行旋转 =====
            for medium_idx = 0:num_row_blocks-1  % 遍历每个中块
                medium_start = medium_idx*n + 1;
                medium_end = (medium_idx+1)*n;
                
                for p = 0:c-1  % 检查每个p行
                    if p + i >= c  % 满足条件时才旋转
                        % 收集当前中块内所有小块的p+1行
                        rows_to_rotate = zeros(num_small_per_medium, d, 'like', A);
                        for small_idx = 0:num_small_per_medium-1
                            row_in_medium = p + 1 + small_idx*c;
                            rows_to_rotate(small_idx+1, :) = ...
                                rotated_block(medium_start + row_in_medium - 1, :);
                        end
                        
                        % 执行向下循环旋转（最后一块移到第一块）
                        rotated_rows = [rows_to_rotate(end, :); rows_to_rotate(1:end-1, :)];
                        
                        % 将旋转后的行放回原位置
                        for small_idx = 0:num_small_per_medium-1
                            row_in_medium = p + 1 + small_idx*c;
                            rotated_block(medium_start + row_in_medium - 1, :) = ...
                                rotated_rows(small_idx+1, :);
                        end
                    end
                end
            end
            
            % ===== 第三步：列旋转 =====
            % 对整个矩阵进行列旋转(向左循环移动j*d/(n/c)+i*d/n列)
            rotated_block = circshift(rotated_block, [0, -j*d/(n/c)-i*d/n]);
            
            % ===== 第四步：填充pt_A =====
            % 将结果放入pt_A的对应位置
            temp_pt_A(j*d+1:(j+1)*d, i*d+1:(i+1)*d) = rotated_block;
        end
    end

    %5. c < 1时的后处理
        c = c_temp;
        if c < 1
            num_chunks = 1/c; % 分成多少个大块
            chunk_size = c * d; % 每个大块的列数
            % 将 temp_pt_A 按列拆分成 num_chunks 个大块
            chunks = mat2cell(temp_pt_A, size(temp_pt_A, 1), repmat(chunk_size, 1, num_chunks));
            % 垂直堆叠这些大块
            pt_A = vertcat(chunks{:});
        else
            pt_A = temp_pt_A;
        end
end