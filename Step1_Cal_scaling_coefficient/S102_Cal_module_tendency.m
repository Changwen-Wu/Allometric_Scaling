clear;
cd E:\Projects\Allo_Scaling\

load Data\Beta\beta_deg_hcp.mat
load Data\Spin_rotation\Schaefer400_rotation.mat


class = table2array(readtable('Data\Module_division\voneconomo_Schaefer400.csv'));
% ['pm', 'assoc1', 'assoc2', 'pss', 'ps', 'lim', 'ins']
result = permutation_class(beta_all,class,perm_id);
result.FC.sig_percentage = result.FC.sig_deg_number./sum(result.FC.sig_deg_number,2);
result.SC.sig_percentage = result.SC.sig_deg_number./sum(result.SC.sig_deg_number,2);
result.GMV.sig_percentage = result.GMV.sig_deg_number./sum(result.GMV.sig_deg_number,2);
save Result\S102_Module_tendency\hcp_von_tendency_deg.mat result

load Data\Module_division\Yeo_7Net_400.mat
subnet(subnet>7) = subnet(subnet>7)-7;
class = subnet;
result = permutation_class(beta_all,class,perm_id);
result.FC.sig_percentage = result.FC.sig_deg_number./sum(result.FC.sig_deg_number,2);
result.SC.sig_percentage = result.SC.sig_deg_number./sum(result.SC.sig_deg_number,2);
result.GMV.sig_percentage = result.GMV.sig_deg_number./sum(result.GMV.sig_deg_number,2);
save Result\S102_Module_tendency\hcp_yeo_tendency_deg.mat result

class = table2array(readtable('Data\Module_division\mesulam_schaefer400.csv'));
% 1  paralimbic
% 2  hetermodal
% 3  unimodal
% 4  idiotypic
result = permutation_class(beta_all,class,perm_id);
result.FC.sig_percentage = result.FC.sig_deg_number./sum(result.FC.sig_deg_number,2);
result.SC.sig_percentage = result.SC.sig_deg_number./sum(result.SC.sig_deg_number,2);
result.GMV.sig_percentage = result.GMV.sig_deg_number./sum(result.GMV.sig_deg_number,2);
save Result\S102_Module_tendency\hcp_laminar_tendency_deg.mat result




function result = permutation_class(beta_all,class,perm_id)

Permutation_times = 10000;
mode = {'SC', 'FC', 'GMV'};
class_num = length(unique(class));

result = {};

for imod = 1 : length(mode)
    thr = beta_all.(mode{imod}).sig_thr;
    if isempty(thr)
        continue;
    end
    beta = beta_all.(mode{imod}).beta;
    beta_p = beta_all.(mode{imod}).p_beta;

    beta_pos_neg = zeros(2,length(beta));
    beta_pos_neg(1,:) = (beta>1).*(beta_p<=thr);
    beta_pos_neg(2,:) = (beta<1).*(beta_p<=thr);
    P_class = zeros(2,class_num);
    sig_deg_number = zeros(2,class_num);
    for i = 1:2 % 1:pos 2:neg
        beta = beta_pos_neg(i,:);
        sig_deg_number_perm = zeros(Permutation_times,class_num);

        for iclass = 1 : class_num
            class_deg = beta(class == iclass);
            sig_deg_number(i,iclass) = sum(class_deg(:));
        end

        perm_beta = beta(perm_id);
        for iclass = 1 : class_num
            class_deg = perm_beta(class == iclass,:);
            sig_deg_number_perm(:,iclass) = sum(class_deg);
        end

        P_class(i,:) = (sum((sig_deg_number(i,:) <= sig_deg_number_perm)) +1) / (Permutation_times + 1);
    end

    result.(mode{imod}).thr = gretna_FDR(P_class(:),0.05);
    result.(mode{imod}).P_class = P_class;
    result.(mode{imod}).sig_deg_number = sig_deg_number;
end
end