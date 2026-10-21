clc
clear
%%
folders=dir('C:\Users\Xianpeng\Desktop\Testdata\nointerpolation\North\*.csv');
num=length(folders);

for i = 1:num
    fname=folders(i).name;
    fpath=fullfile(folders(i).folder,fname);
    data_temp=readmatrix(fpath);
    data{i}=data_temp(:,2);
end


max_length = 0; % 初始化
for i = 1:length(data)
    max_length = max(max_length, length(data{i}));
end

% 2. 将所有列填充到最大长度，不足的部分填充 NaN
for i = 1:length(data)
    current_length = length(data{i});
    
    % 如果当前列比最大长度短，填充 NaN
    if current_length < max_length
        data{i} = [data{i}; NaN(max_length - current_length, 1)];
    end
end

merged_data = cat(2, data{:});  
