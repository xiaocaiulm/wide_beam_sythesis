% function a = getw(w)
%     a = reshape([real(w)',imag(w)'], [], 1);
% end
function a = getw(w)
    % 获取矩阵的列数
    [~, t] = size(w);
    
    % 初始化输出向量
    a = [];
    
    % 遍历每一列，提取实部和虚部，然后加入到输出向量
    for i = 1:t
        real_part = real(w(:, i));
        imag_part = imag(w(:, i));
        
        % 将实部和虚部转置并添加到a中
        a = [a; real_part; imag_part];
    end
end
