clear; clc;
cd D:\Projects\Allo_Scaling

load Data\Beta\beta_deg_hcp.mat
mode = {'SC', 'FC', 'GMV'};

%% Enigma 数据库（weighted mean方法）
load Data\Disorder\DK_Sch_overlap.mat
load Data\Spin_rotation\DK40_rotation.mat
load Data\Disorder\summary_ct_dk68.mat

beta_dk = zeros(length(mode),70);


for imod = 1 : length(mode)
    beta = beta_all.(mode{imod}).beta;

    beta(isnan(beta)) = 0; % 在加权求和时忽略sc的nan
    for i = 1 :35
        if(sum(Dice_LH(:,i))) == 0
            continue
        end
        Dice = Dice_LH(:,i);
        beta_dk(imod,i) = beta(1:200)*Dice/ sum(Dice); %加权平均
        %         Dice = Dice_LH(:,i);
        %         beta_dk(imod,i) = beta(1:200)*Dice; %加权求和
        %         beta_dk(imod,i) = beta(1:200)*(Dice==max(Dice)); % winner takes all
    end

    for i = 1 :35
        if(sum(Dice_RH(:,i))) == 0
            continue
        end
        Dice = Dice_RH(:,i);
        beta_dk(imod,i+35) = beta(201:400)*Dice/ sum(Dice); %加权平均
        %         beta_dk(imod,i+35) = beta(201:400)*Dice; %加权求和
        %         beta_dk(imod,i+35) = beta(201:400)*(Dice==max(Dice)); % winner takes all
    end
end


beta_dk(:,[4,39]) = []; % 去掉皮层下

disorder_name = summary_ct_dk68(1,2:14);
summary_ct_dk68 = cell2mat(summary_ct_dk68(2:end,2:end));
summary_ct_dk68 = summary_ct_dk68(:,1:13);

% summary_ct_dk68 = abs([summary_ct_dk68, disorder_prop']);

disorder_result = {};
for imod = 1 : length(mode)
    result = {};
    R = zeros(1,size(summary_ct_dk68,2));
    P = zeros(1,size(summary_ct_dk68,2));
    map1 = beta_dk(imod,:);
    for i = 1 : size(summary_ct_dk68,2)
        map2 = summary_ct_dk68(:,i);
        if size(map1, 1) == 1; map1 = map1.'; end
        if size(map2, 1) == 1; map2 = map2.'; end
        [r,~] = corr(map1, map2, 'type', "Spearman");
        [p_spin, ~] = perm_sphere_p_wcw(map1, map2, perm_id, "Spearman");
        R(i) = r;
        P(i) = p_spin;
    end

    thr = gretna_FDR(P,0.05);
    result.R = R; result.P = P;
    result.thr = thr;

    disorder_result.(mode{imod}) = result;
end

save Result\S202_Corr_disorder\hcp_disorder_enigma.mat disorder_result
clear

%% Misic CBF 的文章
cd D:\Projects\Allo_Scaling\

load Data\Beta\beta_deg_hcp.mat
load Data\Spin_rotation\Schaefer400_rotation.mat
new_disorder = dir('Data\Disorder\CBF_misic\*.npy');

disorder_data = nan(400,length(new_disorder));
for i = 1 : length(new_disorder)
    data = readNPY(fullfile(new_disorder(i).folder, new_disorder(i).name));
    disorder_data(:,i) = data;
end

disorder_name = {'3Rtau' '4Rtau' 'DLB' 'EOAD' 'LOAD' 'PS1' 'TDP43A' 'TDP43C'};

[result] = spin_spearman(disorder_data, beta_all, perm_id);
result.label = disorder_name;

save Result\S202_Corr_disorder\hcp_disorder_misic.mat result






function [result] = spin_spearman(signatures, beta_all, perm_id)

if size(signatures, 2) == 400; signatures = signatures.'; end

mode = {'SC', 'FC', 'GMV'};
result = {};

P = zeros(3,size(signatures, 2));
R = zeros(3,size(signatures, 2));

thr = zeros(1,length(mode));

for imod = 1 : length(mode)
    beta = beta_all.(mode{imod}).beta;
    for ifield = 1 : size(signatures, 2)
        data = signatures(:,ifield);

        if size(beta, 1) == 1; beta = beta.'; end

        [P(imod,ifield), ~] = perm_sphere_p_wcw(beta, data, perm_id, 'spearman');
        [R(imod,ifield),~] = corr(beta,data,'type','Spearman','rows','pairwise');
    end
   thr(imod) = gretna_FDR(P(imod,:),0.05);
end

result.P = P; result.R = R; result.thr = thr;

end