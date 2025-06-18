function Results = gretna_allometric_scaling(X_subcomponent, X_global, Beta_predicted, Covariates)

%==========================================================================
% This function is used to perform regression analysis to study allometric
% scaling law. Allometric scaling law is a standard approach to explore how
% biological characteristics change with size. Recently, it is also used to
% quantify how the size of individual brain regions changes relative to
% changes in overall brain size.
%
% Syntax: function Results = gretna_allometric_scaling(X_subcomponent, X_global, Beta_predicted, Covariates)
%
% Inputs:
%  X_subcomponent:
%                 M*N data array with M being observations and N variables
%                 (e.g., cortical thickness of N regions from M subjects).
%        X_global:
%                 M*1 data array with M being observations (e.g., total
%                 cortical thickness or brain size of the M subjects).
%  Beta_predicted:
%                 A scalar denoting the predicted scaling exponent under
%                 the null hypothesis that cortical structure is invariant
%                 under scaling. Specifically, the predicted scaling
%                 exponents of total intracranial volume for cortical volume,
%                 surface area, cortical thickness, sulcal depth, and mean
%                 curvature are 1, 2/3, 1/3, 1/3, and ¨C1/3, respectively. 
% Covariates (opt):
%                 M*L data array with M being observations and L variables
%                 (e.g., age and sex of the M subjects).
%
% Outputs:
%    Results:
%           .beta:
%                 1*N data array storing allometric scaling coefficients
%                 of all variables.
%        .se_beta:
%                 1*N data array storing the standard errors of all
%                 allometric scaling coefficients derived from
%                 bootstrapping (n = 10000).
%         .t_beta:
%                 1*N data array storing t-values of all allometric scaling
%                 coefficients. The t-tests examine whether the allometric
%                 scaling coefficients deviate from 1.
%         .p_beta:
%                 1*N data array storing p-values of the t-values.
%
% References:
%  1. de Jong et al., (2017) Allometric scaling of brain regions to
%     intra-cranial volume: An epidemiological MRI study. Human Brain
%     Mapping, 38:151-164.
%  2. Janssen et al., (2022) Longitudinal Allometry of Sulcal Morphology
%     in Health and Schizophrenia. Journal of Neuroscience, 42:3704-3715.
%  3. Im et al., (2008) Brain Size and Cortical Structure in the Adult
%     Human Brain. Cerebral Cortex, 18:2181--2191.
%
% Jiashun Zhang, IBRR, SCNU, Guangzhou, 2025/3/25, 2024024661@m.scnu.edu.cn
% Jinhui Wang, IBRR, SCNU, Guangzhou, 2025/3/25, jinhui.Wang.1982@gmail.com
%==========================================================================

Num_subs = size(X_global,1);

if size(X_subcomponent,1) ~= Num_subs
    error('The number of subjects are not equal between ''X_global'' and ''X_subcomponent''!');
end

if nargin == 4
    if size(Covariates,1) ~= Num_subs
        error('The number of subjects are not equal between ''X_global'' and ''Covariates''!');
    end
    
    if size(Covariates,1) ~=  size(X_subcomponent,1)
        error('The number of subjects are not equal between ''X_subcomponent'' and ''Covariates''!');
    end
end

if any(X_subcomponent(:) <= 0)
    X_subcomponent(X_subcomponent <= 0) = nan;
    warning('Negative values in ''X_subcomponent'' are ignored!');
end

if any(X_global <= 0)
    X_global(X_global <= 0) = nan;
    warning('Negative values in ''X_global'' are ignored!');
end

if nargin == 3;
    Covariates = [];
end

Num_subcomponents = size(X_subcomponent,2);
Results.beta      = zeros(1, Num_subcomponents);
Results.se_beta   = zeros(1, Num_subcomponents);
Results.t_beta    = zeros(1, Num_subcomponents);
Results.p_beta    = zeros(1, Num_subcomponents);

X_global   = log10(X_global);
Desnmatrix = [X_global Covariates ones(Num_subs, 1)];

for isubcom = 1:Num_subcomponents
    
    X_subcomponent_isubcom = log10(X_subcomponent(:,isubcom));
    
    beta      = regress(X_subcomponent_isubcom, Desnmatrix);
    boot_beta = bootstrp(10000, @regress, X_subcomponent_isubcom, Desnmatrix);
    
    Results.beta(isubcom)    = beta(1);
    Results.se_beta(isubcom) = std(boot_beta(:,1));
    Results.t_beta(isubcom)  = (Results.beta(isubcom) - Beta_predicted) / Results.se_beta(isubcom);
    Results.p_beta(isubcom)  = 2 * (tcdf(-abs(Results.t_beta(isubcom)), Num_subs - 1)); % two tail
end

return