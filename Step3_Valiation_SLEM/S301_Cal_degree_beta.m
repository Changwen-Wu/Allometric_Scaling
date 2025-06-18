clear; clc;
cd D:\Projects\Allo_Scaling\

Roi_num = 400;
sub_num = 464;

load Data\Networks\SLEM\subinfo_464.mat
gender = sub_info(:,2);
age = sub_info(:,3); 

filelist.SC = dir('Data\Networks\SLEM\SCN\*.mat');
filelist.FC = dir('Data\Networks\SLEM\FCN\*.mat');
filelist.GMV = dir('Data\Networks\SLEM\MCN\*.mat');

net_names = fieldnames(filelist);

beta_all = {};

[~,temp,~] = gretna_read_image('Data\Template\Schaefer_400_7Net_MNI152_3mm.nii');




for inet = 1 : length(net_names)
    tic
    result = {};

    net_name = net_names{inet};
    net_list =  filelist.(net_name);

    degree_Glo = zeros(1,sub_num);
    degree = zeros(Roi_num,sub_num);

    volume_effect = zeros(400,1);
    for ind = 1 : 400
        region_map = (temp == ind);
        volume_effect(ind) = sum(region_map(:));
    end
    volume_effect = (volume_effect(:,1) + volume_effect(:,1)')/2;

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

        degree_Glo(isub) = sum(network(:)); degree(:,isub) = sum(network);
    end


    % 初始化结果存储
    result.beta = zeros(1, Roi_num);
    result.t_beta = zeros(1, Roi_num);
    result.p_beta = zeros(1, Roi_num);
    result.cov_eff = zeros(3, Roi_num);
    result.beta_SE = zeros(1, Roi_num);
    for iroi = 1 : Roi_num
        X = log10(degree_Glo');
        Y = log10(degree(iroi,:)');
        inf_ind = or(isinf(Y),isinf(X));

        if sum(inf_ind/sub_num) > 0.3
            result.beta(iroi) = nan;
            result.cov_eff(:,iroi) = nan;
            result.beta_SE(iroi) = nan;
            result.t_beta(iroi) = nan;
            result.p_beta(iroi) = nan;
            continue;
        end

        Y = Y(~inf_ind); X = X(~inf_ind);
        [b,bint,~,~,~] = regress(Y,[X,age(~inf_ind),gender(~inf_ind),age(~inf_ind).*gender(~inf_ind),ones(sum(~inf_ind),1)]);

        % 提取回归系数和标准误差
        b_SE = (bint(:,2) - bint(:,1)) / (2 * 1.96); % 计算标准误差（假设 95% 置信区间）
        t_values = b ./ b_SE; % 计算 t 统计量
        p_values = 2 * (1 - tcdf(abs(t_values), length(Y) - size(X, 2))); % 计算双尾 p 值

        result.beta(iroi) = b(1);
        result.cov_eff(:,iroi) = p_values(2:4);
        
        boot_beta = bootstrp(10000, @bootstrapFun, Y, X, age(~inf_ind), gender(~inf_ind));

        result.beta_SE(iroi) = std(boot_beta);
        result.t_beta(iroi) = (result.beta(iroi) - 1)./ result.beta_SE(iroi);
        result.p_beta(iroi) = 2 * (tcdf(-abs(result.t_beta(iroi)), sub_num - 1));
    end

    result.sig_thr = gretna_FDR(result.p_beta,0.05);

    if ~isempty(result.sig_thr)
        sig_ind = result.p_beta <= result.sig_thr;
        result.valid_beta = result.beta .* sig_ind;
    end

    beta_all.(net_name) = result;
    toc
end

save Data\Beta\beta_deg_slem.mat beta_all


function coefficients = bootstrapFun(Y, X, age, gender)
    % 在 Bootstrap 函数中重新拟合线性模型并返回回归系数

    % 添加截距项到自变量矩阵 X
    X = [X, age, gender, age .* gender, ones(length(Y), 1)];

    % 使用 regress 函数进行线性回归
    b = regress(Y, X);

    % 返回回归系数中的第一个系数（例如您需要的 beta）
    coefficients = b(1);
end

