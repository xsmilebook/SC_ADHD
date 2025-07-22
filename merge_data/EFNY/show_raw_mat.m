% --- 1. 设置路径和参数 ---
% 包含SC矩阵文件的文件夹路径
data_folder = '/ibmgpfs/cuizaixu_lab/congjing/brainproject/development/results/SCmat_Du15_all';

% 要匹配的文件后缀
file_pattern = '*_sift_invnodevol_count.csv';

% 图像保存路径和文件名
output_image_file = '/ibmgpfs/cuizaixu_lab/xuhaoshu/SC_ADHD/datasets/EFNY/figures/average_raw_SC_matrix_Du15.png';



original_matrix_size = 16;
final_matrix_size = 15;

brain_region_labels = {
    'VIS-P', 'CG-OP', 'DN-B', 'SMOT-B', 'AUD', ...
    'PM-PPr', 'dATN-B', 'SMOT-A', 'LANG', 'FPN-B', ...
    'FPN-A', 'dATN-A', 'VIS-C', 'SAL/PMN', 'DN-A'
};

if length(brain_region_labels) ~= final_matrix_size
    error('标签数量 (%d) 与矩阵维度 (%d) 不匹配。', ...
          length(brain_region_labels), final_matrix_size);
end

% --- 2. 查找文件 (代码不变) ---
search_path = fullfile(data_folder, file_pattern);
file_list = dir(search_path);
if isempty(file_list)
    error('未找到匹配文件。');
end
fprintf('找到了 %d 个匹配的文件。\n', length(file_list));

% --- 3. 读取和计算平均矩阵 (代码不变) ---
sum_matrix = zeros(final_matrix_size, final_matrix_size);
valid_file_count = 0;
fprintf('开始读取并处理矩阵...\n');
for i = 1:length(file_list)
    current_filepath = fullfile(data_folder, file_list(i).name);
    try
        original_matrix = readmatrix(current_filepath);
        if all(size(original_matrix) == [original_matrix_size, original_matrix_size])
            processed_matrix = original_matrix(1:final_matrix_size, 1:final_matrix_size);
            sum_matrix = sum_matrix + processed_matrix;
            valid_file_count = valid_file_count + 1;
        else
            fprintf('警告: 文件 %s 尺寸不正确，已跳过。\n', file_list(i).name);
        end
    catch ME
        fprintf('错误: 无法读取文件 %s。错误信息: %s\n', file_list(i).name, ME.message);
    end
end
if valid_file_count == 0
    error('没有文件被成功处理。');
end
average_matrix = sum_matrix / valid_file_count;
fprintf('共处理了 %d 个文件，平均矩阵已计算。\n', valid_file_count);

% --- 4. 可视化 (plotmatrix 风格) ---
fprintf('开始生成 plotmatrix 风格的图像...\n');
figure('Name', 'Average SC Matrix (Plotmatrix Style)', 'NumberTitle', 'off', 'Color', 'w');

% --- 4.1 创建发散型颜色图 ---
% 我们需要一个从蓝到白到红的颜色图
% cbrewer 函数可以提供高质量的颜色图，如果未安装，可以使用内置的
try
    cmap = cbrewer('div', 'RdBu', 256);
    cmap = flipud(cmap); % RdBu 默认是红-白-蓝，翻转得到蓝-白-红
catch
    fprintf('未找到 cbrewer 函数，使用自定义颜色图代替。\n');
    % 手动创建一个简单的蓝-白-红颜色图
    blue = [0, 0, 1];
    white = [1, 1, 1];
    red = [1, 0, 0];
    gradient1 = [linspace(blue(1), white(1), 128)', linspace(blue(2), white(2), 128)', linspace(blue(3), white(3), 128)'];
    gradient2 = [linspace(white(1), red(1), 128)', linspace(white(2), red(2), 128)', linspace(white(3), red(3), 128)'];
    cmap = [gradient1; gradient2];
end

% --- 4.2 确定对称的颜色范围 ---
max_abs_val = max(abs(average_matrix(:))); % 找到矩阵中绝对值的最大值
color_limits = [-max_abs_val, max_abs_val]; % 设置颜色范围

% --- 4.3 绘制热图 ---
imagesc(average_matrix, color_limits);
colormap(cmap);
colorbar;
axis square; % 确保图像是方形的

% --- 4.4 设置坐标轴和标签 ---
ax = gca;
ax.XTick = 1:final_matrix_size;
ax.YTick = 1:final_matrix_size;
ax.XTickLabel = brain_region_labels;
ax.YTickLabel = brain_region_labels;
ax.XTickLabelRotation = 90;
ax.TickLength = [0 0]; % 隐藏刻度线本身，只保留标签
set(ax, 'FontSize', 10, 'LineWidth', 1.5, 'FontName', 'Arial'); % 设置字体和边框线宽

% --- 4.5 添加边框和对角线 ---
hold on; % 保持当前图像，以便在其上添加新的图形元素

% 添加黑色边框 (画四条线)
plot([0.5, final_matrix_size + 0.5], [0.5, 0.5], 'k-', 'LineWidth', 1.5); % 上边
plot([0.5, final_matrix_size + 0.5], [final_matrix_size + 0.5, final_matrix_size + 0.5], 'k-', 'LineWidth', 1.5); % 下边
plot([0.5, 0.5], [0.5, final_matrix_size + 0.5], 'k-', 'LineWidth', 1.5); % 左边
plot([final_matrix_size + 0.5, final_matrix_size + 0.5], [0.5, final_matrix_size + 0.5], 'k-', 'LineWidth', 1.5); % 右边

% 添加黑色对角线
plot([0.5, final_matrix_size + 0.5], [0.5, final_matrix_size + 0.5], 'k-', 'LineWidth', 1.5);

hold off; % 结束在当前图像上绘图

% --- 4.6 添加标题 ---
title('Average Structural Connectivity Matrix (Du15 Atlas)');
xlabel('Brain Regions');
ylabel('Brain Regions');

% --- 5. 保存图像 ---
saveas(gcf, output_image_file);
fprintf('图像已成功保存到: %s\n', output_image_file);