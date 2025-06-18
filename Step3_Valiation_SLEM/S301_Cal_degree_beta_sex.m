clear; clc;
cd D:\Projects\Allo_Scaling\
% 0/1 男女
Roi_num = 400;
sub_num = 464;
%% SLEM
load Data\Networks\SLEM\subinfo_464.mat
gender = sub_info(:,2);
age = sub_info(:,3);


filelist.SC = dir('Data\Networks\SLEM\SCN\*.mat');
filelist.FC = dir('Data\Networks\SLEM\FCN\*.mat');
filelist.GMV = dir('Data\Networks\SLEM\MCN\*.mat');

net_names = fieldnames(filelist);

beta_all = {};

[~,temp,~] = gretna_read_image('Data\Template\Schaefer_400_7Net_MNI152_3mm.nii');

volume_effect = zeros(400,1);
for ind = 1 : 400
    region_map = (temp == ind);
    volume_effect(ind) = sum(region_map(:));
end
volume_effect = (volume_effect(:,1) + volume_effect(:,1)')/2;

for inet = 1 : length(net_names)
    tic
    result = {};

    net_name = net_names{inet};
    net_list =  filelist.(net_name);

    degree_Glo = zeros(1,sub_num);
    degree = zeros(Roi_num,sub_num);

    for isub = 1: sub_num
        file = load(fullfile(net_list(isub).folder, net_list(isub).name));
        tmp = fieldnames(file);

        network = file.(tmp{1});

        if strcmp(net_name,'FC')
            network(network<0) = 0;
        end

        if ~strcmp(net_name,'SC') %% 已经去过了
            network = network ./ volume_effect;
        end

        degree_Glo(isub) = sum(network(:));  degree(:,isub) = sum(network);
    end


    % 初始化结果存储
    result.beta = zeros(2, Roi_num);
    result.t_beta = zeros(2, Roi_num);
    result.p_beta = zeros(2, Roi_num);
    beta_SE = zeros(1, Roi_num);

    for igender = [0,1]
        sub_index = logical(gender-igender);
        for iroi = 1 : Roi_num

            X = log10(degree_Glo(sub_index)');
            Y = log10(degree(iroi,sub_index)');

            inf_ind = or(isinf(Y),isinf(X));
            sub_age = age(sub_index); 

            if sum(inf_ind/sum(sub_index)) > 0.3
                result.beta(igender+1, iroi) = nan;
                result.cov_eff(:,iroi) = nan;
                result.t_beta(igender+1, iroi) = nan;
                result.p_beta(igender+1, iroi) = nan;
                continue;
            end

            Y = Y(~inf_ind); X = X(~inf_ind); sub_age = sub_age(~inf_ind);
            [b,bint,~,~,~] = regress(Y,[X,sub_age,ones(length(sub_age),1)]);
            boot_beta = bootstrp(10000, @bootstrapFun, Y, X, sub_age);


            result.beta(igender+1, iroi) = b(1);

            beta_SE(iroi) = std(boot_beta);
            result.t_beta(igender+1, iroi) = (result.beta(igender+1, iroi) - 1)./ beta_SE(iroi);
            result.p_beta(igender+1, iroi) = 2 * (tcdf(-abs(result.t_beta(igender+1, iroi)), sum(sub_index) - 1));
        end
    end
    result.sig_thr(1) = gretna_FDR(result.p_beta(1,:),0.05);
    result.sig_thr(2) = gretna_FDR(result.p_beta(2,:),0.05);

    if ~isempty(result.sig_thr(1))
        sig_ind = result.p_beta(1,:) <= result.sig_thr(1);
        result.valid_beta(1,:) = result.beta(1,:) .* sig_ind;
    end

    if ~isempty(result.sig_thr(2))
        sig_ind = result.p_beta(2,:) <= result.sig_thr(2);
        result.valid_beta(2,:) = result.beta(2,:) .* sig_ind;
    end
    beta_all.(net_name) = result;
    toc
end

save Data\Beta\slem_beta_deg_sex.mat beta_all


function coefficients = bootstrapFun(Y, X, age)
% 在 Bootstrap 函数中重新拟合线性模型并返回回归系数

% 添加截距项到自变量矩阵 X
X = [X, age, ones(length(Y), 1)];

% 使用 regress 函数进行线性回归
b = regress(Y, X);

% 返回回归系数中的第一个系数（例如您需要的 beta）
coefficients = b(1);
end

