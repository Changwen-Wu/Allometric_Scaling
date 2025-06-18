clear; clc
cd D:\Projects\Allo_Scaling

load('Data\Beta\beta_deg_hcp.mat'); HCP = beta_all;
load('Data\Beta\beta_deg_slem.mat'); SLEM = beta_all;

load('Data\Spin_rotation\Schaefer400_rotation.mat')

mod = {'SC', 'FC', 'GMV'};
P = zeros(1,length(mod)); R = P;


for imod = 1 : length(mod)
    ori = HCP.(mod{imod}).beta';
    val = SLEM.(mod{imod}).beta';

    [P(imod), ~] = perm_sphere_p_wcw(ori, val, perm_id, 'spearman');
    [R(imod),~] = corr(ori,val,'type','Spearman','rows','pairwise');
end

