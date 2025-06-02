function ct_C = plaintext_ciphertext_matmul(pt_A, ct_B, d, n, c)
    % 明文-密文矩阵乘法
    % 输入参数：
    %   pt_A - 明文矩阵，维度为 R^{dn/c×cd}，由原始矩阵 A ∈ R^{d×d} 转换而来
    %   ct_B - 密文矩阵，维度为 R^{n/c×cd}，由原始矩阵 B ∈ R^{d×n} 转换而来
    %   d    - 矩阵A的维度
    %   n    - 矩阵B的维度 (需满足 n ≤ d)
    %   c    - 打包参数（每个密文包含的对角线数量）
    % 输出参数：
    %   ct_C - 结果密文矩阵，表示 C = AB 的加密结果，处理方式与 ct_B 相同
    
    % ==================== 初始化阶段 ====================
    % 计算对角线块数量：将d×d矩阵划分为(d/n)个n×n的对角线块
    d_tilde = d / n;  
    % 计算每个块内的密文数量：将n×n块划分为(n/c)个密文
    n_tilde = n / c;   
    
    % 预分配中间结果矩阵（维度为 (d/c) × (c*d)）
    ct_intermediate = zeros(d/c * max((d/n)/(c*d),1), c*d, 'like', pt_A);
    % 预分配最终结果矩阵（维度为 n_tilde × (c*d)）
    ct_C = zeros(n_tilde, c*d, 'like', pt_A);
    
    % ==================== 第一阶段：子矩阵并行乘法 ====================
    for m = 0:1/(min(c,1))-1
        % 遍历每个对角线块 (i ∈ [0, d_tilde-1])
        for i = 0:d_tilde-1
            % 遍历当前块内的每个密文位置 (r ∈ [0, n_tilde-1])
            for r = 0:(n/max(c,1))-1
                % 初始化累加器（存储当前密文位置的中间结果）
                temp_sum = zeros(1, c*d); 
                
                % 遍历B矩阵的密文索引 (j ∈ [0, n_tilde-1])
                for j = 0:(n/max(c,1))-1
                    % 遍历每个对角线位置 (ell ∈ [0, c-1])
                    for ell = 0:max(c,1)-1
                        % 步骤1：对密文进行循环移位（左移d*ell位）
                        rotated_ct = circshift(ct_B(m * n + j+1, :), -d*ell);
                        
                        % 步骤2：选择对应的明文行
                        % 计算明文行号 = 块偏移 + 密文偏移 + 对角线偏移
                        pt_row = m * (d*n) + i * (n^2/(max(c,1))) + r * n + j * max(c,1) + ell + 1;
                        pt_selected = pt_A(pt_row, :);
                        
                        % 调试打印
%                         disp('---------------------------------');
%                         fprintf('i=%d, r=%d, j=%d, ell=%d\n', i, r, j, ell);                    
%                         fprintf('pt_row = %d\n', pt_row);
%                         
%                         if isa(pt_selected, 'sym')
%                             disp('pt_selected (symbolic):');
%                             disp(char(pt_selected(1:end)));
%                         else
%                             disp('pt_selected (numeric):');
%                             disp(pt_selected(1:end));
%                         end
%                         
%                         if isa(rotated_ct, 'sym')
%                             disp('rotated_ct (symbolic):');
%                             disp(char(rotated_ct(1:end)));
%                         else
%                             disp('rotated_ct (numeric):');
%                             disp(rotated_ct(1:end));
%                         end
                        
                        % 步骤3：计算密文与明文的逐元素乘积
                        product = rotated_ct .* pt_selected;
                        
                        % 步骤4：累加到当前密文位置的中间结果
                        temp_sum = temp_sum + product;
                    end
                end
                
                % 将当前密文位置的中间结果存入矩阵
                ct_intermediate(m * d + i * n / (max(c,1)) + r + 1, :) = temp_sum;
                
                % 调试信息（显示中间结果）
%                 disp('-------------ct_intermediate--------------------');
%                 disp('ct_intermediate index:')
%                 disp(m * (d*n) + i * n / (max(c,1)) + r + 1);
%                 if isa(temp_sum, 'sym')
%                     disp('ct_intermediate (symbolic):');
%                     disp(char(temp_sum(1:end)));
%                 else
%                     disp('ct_intermediate (numeric):');
%                     disp(temp_sum(1:end));
%                 end
            end
        end
    end
    % 第一阶段，密文必要旋转次数(n/c) * (max(c,1)-1)，涉及密文的乘法次数dn/c
    
    % ==================== 第二阶段：带旋转的结果聚合 ====================
    % ct_intermediate的尺寸: d/c * max((d/n)/(c*d),1) × cd, 如果d/n(划分出来的小块的尺寸) 不与 cd 有整倍数关系,
    % 进一步, 若 d/n>cd 那么就需要跨行旋转; 若cd>d/n最后一个d/n块就需要跨行旋转。跨行在加密状态
    % 下就意味着不同的密文, 即使cd 和 d/n存在整数倍关系，若d/n>cd,
    % 跨行的旋转本身也会增加需要旋转的次数(下面的代码面向cd >= d/n 且整除, 若d/n>cd需要调整代码);
    % 例如cd = 4, d/n = 8时, 假定[a b c d][e f g h]这两行需要跨行一同右旋1, 则操作步骤:
    % 1. 提取[a b c 0][0 0 0 d] 和 [e f g 0][0 0 0 h], 2行各自需要1次乘法
    % 2. 旋转得到[0 a b c][d 0 0 0] 和 [0 e f g][h 0 0 0], 2行各自需要2次旋转
    % 3. 得到[h a b c] [d e f g]
    % 可见对于需要旋转的1×d/n，为了完成其内部旋转需要的乘法和旋转次数从这些数字在同一行时的1次乘法2次旋转
    % 变成了2次乘法4次旋转。
    % 进一步，当cd和d/n不存在整数倍关系时, 例如cd = 3, d/n = 4, [a b c][d e f][g h i][...]
    % 需要(a b c d)间, (e f g h)间内部旋转; 或者cd = 4, d/n = 3, [a b c d] [e f g h]
    % 需要(a b c)间, (d e f)间内部旋转, 这就会非常麻烦, 所以干脆避免一切复杂情况
    
    % 定义块大小参数：
    big_block_rows = (d / n)*(max((d/n)/(c*d),1));         % 每个大块的行数
    small_block_cols = min(c * d, d / n);       % 每个小块的列数
    num_small_blocks_per_row = c * n; % 每行包含的小块数量
    
    % 遍历每个输出密文位置 (i ∈ [0, n_tilde-1])
    for i = 0:n_tilde-1
        % 步骤1：提取当前大块数据（d/n行 × c*d列）
        big_block_start = i * big_block_rows + 1;
        big_block_end = (i+1) * big_block_rows;
        current_big_block = ct_intermediate(big_block_start:big_block_end, :);
        
        % 初始化旋转后的大块矩阵
        rotated_big_block = zeros(size(current_big_block), 'like', current_big_block);
        
        % 步骤2：处理大块中的每一行 (j ∈ [0, big_block_rows-1])
        for j = 0:big_block_rows-1
            % 获取当前行数据
            current_row = current_big_block(j+1, :);
            
            % 将行数据分割为多个小块（每个小块d/n列）
            small_blocks = mat2cell(current_row, 1, small_block_cols * ones(1, num_small_blocks_per_row));
            
            % 对每个小块进行循环右移j位
            rotated_small_blocks = cell(size(small_blocks));
            for k = 1:num_small_blocks_per_row
                rotated_small_blocks{k} = circshift(small_blocks{k}, [0, j]);
            end
            
            % 合并旋转后的小块为完整行
            rotated_row = [rotated_small_blocks{:}];
            rotated_big_block(j+1, :) = rotated_row;
        end
        
        % 步骤3：聚合旋转后大块的所有行，得到最终密文
        ct_C(i+1, :) = sum(rotated_big_block, 1);
    end
    % 第二阶段，密文必要旋转次数2*(d/n-1)*(n/c)，涉及密文的乘法次数(d/n-1)*n/c
    % 如果修改上面的代码使得可以处理d/n>cd且cd整除d/n，
    % 那时，c>=1或nc>=1的情况下旋转和乘法次数不变，nc<1时密文必要旋转次数2*(d/n-1)*(n/c)*(1/n)，涉及密文的乘法次数(d/n-1)*(n/c)*(1/n)
end
