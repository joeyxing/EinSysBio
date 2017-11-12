%% clear
clear
close all
%% initial state: no protein or mRNA
% only amino acid and nuleotide
% A_d = 1000;
% R_d = 1000; % amino acid for A and R protein
% 
% A_mRNA_d = 1000;
% R_mRNA_d = 1000; % Nuleotide for mRNA

A_gene = 1;
R_gene = 1; % number of DNA 0/1

A_promoter = 0;
R_promoter = 0; % number of promoter bound 0/1

A = 0;
R = 0;
AR_complex = 0;

A_mRNA = 0;
R_mRNA = 0;

%% biochem reaction list, p means propensity
% delta: reactant change amount
%[A_d, R_d, A_mRNA_d, R_mRNA_d, A_gene, R_gene, A_promoter, R_promoter, A,   R , AR_complex, A_mRNA, R_mRNA]
%[ 1,   2,      3,        4,||[     5,      6,         7,          8,     9,  10 ,     11,       12,     13  ]
%             [ 1, 2, 3, 4, 5, 6, 7, 8, 9]
delta(1,:)  = [ 0, 0, 0, 0, 1, 0, 0, 0, 0]; % A_d -> A, p = A_mRNA
delta(2,:)  = [ 0, 0, 0, 0, 0, 1, 0, 0, 0]; % R_d -> R, p = R_mRNA*0.1

delta(3,:)  = [ 0, 0, 0, 0, 0, 0, 0, 1, 0]; % A_mRNA_d -> A_mRNA, p = 4*A_gene + 40*A_promoter
delta(4,:)  = [ 0, 0, 0, 0, 0, 0, 0, 0, 1]; % R_mRNA_d -> R_mRNA, p = 0.001*R_gene + 2*R_promoter

delta(5,:)  = [-1, 0, 1, 0,-1, 0, 0, 0, 0]; % A_gene + A -> A_promoter, p = 10*A_gene*A
delta(6,:)  = [ 0,-1, 0, 1,-1, 0, 0, 0, 0]; % R_gene + A -> R_promoter, p = 10*R_gene*A 

delta(7,:)  = [ 1, 0,-1, 0, 1, 0, 0, 0, 0]; % A_promoter -> A + A_gene, p = 10*A_promoter
delta(8,:)  = [ 0, 1, 0,-1, 1, 0, 0, 0, 0]; % R_promoter -> A + R_gene, p = 10*R_promoter

delta(9,:)  = [ 0, 0, 0, 0,-1,-1, 1, 0, 0]; % A + R -> AR_complex, p = 100*A*R

delta(10,:) = [ 0, 0, 0, 0,-1, 0, 0, 0, 0]; % A -> A_d, p = 0.1*A
delta(11,:) = [ 0, 0, 0, 0, 0,-1, 0, 0, 0]; % R -> R_d, p = 0.01*A

delta(12,:) = [ 0, 0, 0, 0, 0, 1,-1, 0, 0]; % AR_complex -> A_d + R, p = 0.1*AR_complex
delta(13,:) = [ 0, 0, 0, 0, 1, 0,-1, 0, 0]; % AR_complex -> R_d + A, p = 0.01*AR_complex

delta(14,:) = [ 0, 0, 0, 0, 0, 0, 0,-1, 0]; % A_mRNA -> A_mRNA_d, p = 1*A_mRNA
delta(15,:) = [ 0, 0, 0, 0, 0, 0, 0, 0,-1]; % R_mRNA -> R_mRNA_d, p = 0.02*R_mRNA

%% Gillespie Simulation
t_end = 20e2;
t = 0;
idx = 1;
reactant = [A_gene, R_gene, A_promoter, R_promoter, A, R, AR_complex, A_mRNA, R_mRNA];
state(:,1) = [t, reactant];
while t < t_end
    % propensity list
%     p_list(1) = reactant(1)*reactant(12);              % A_d + A_mRNA -> A + A_mRNA
%     p_list(2) = 0.1*reactant(2)*reactant(13);          % R_d + R_mRNA -> R + R_mRNA
    p_list(1) = reactant(8);                          % A_d + A_mRNA -> A + A_mRNA    
    p_list(2) = 0.1*reactant(9);                      % R_d + R_mRNA -> R + R_mRNA
    
    p_list(3) = 4*reactant(1) + 40*reactant(3);       % A_mRNA_d -> A_mRNA
    p_list(4) = 0.001*reactant(2) + 2*reactant(4);    % R_mRNA_d -> R_mRNA
    
    p_list(5) = 10*reactant(1)*reactant(5);           % A_gene + A -> A_promoter
    p_list(6) = 10*reactant(2)*reactant(5);           % R_gene + A -> R_promoter
    
    p_list(7) = 10*reactant(3);                       % A_promoter -> A + A_gene
    p_list(8) = 10*reactant(4);                       % R_promoter -> A + R_gene
    
    p_list(9) = 100*reactant(5)*reactant(6);          % A + R -> AR_complex
    
    p_list(10) = 0.1*reactant(5);                     % A -> A_d
    p_list(11) = 0.01*reactant(6);                    % R -> R_d
    
    p_list(12) = 0.1*reactant(7);                     % AR_complex -> A_d + R
    p_list(13) = 0.01*reactant(7);                    % AR_complex -> R_d + A
    
    p_list(14) = reactant(8);                         % A_mRNA -> A_mRNA_d
    p_list(15) = 0.02*reactant(9);                    % R_mRNA -> R_mRNA_d
    
    lambda = sum(p_list);
    tau = exprnd(1/lambda);
    
    t = t + tau;
    
    prob = rand();
    c = cumu(p_list/lambda); % distribution of probability
    nr = find(c > prob, 1) - 1; % number of reaction
    reactant = reactant + delta(nr,:);
    
    idx = idx + 1;
    state(:, idx) = [t, reactant];
end
%% plot total number of A and R
figure;hold on
plot(state(1,:), state(6,:)+state(8,:))
plot(state(1,:), state(7,:)+state(8,:))
legend('A','R')
xlim([0 t_end])

%% plot moving averge of the state of two promotor
ma_Ap = zeros(1, size(state, 2));
ma_Rp = zeros(1, size(state, 2));
s1 = 0;
s2 = 0;
for idx = 1:size(state, 2)
    s1 = s1 + state(4,idx);
    s2 = s2 + state(5,idx);
    ma_Ap(idx) = s1 / idx;
    ma_Rp(idx) = s2 / idx;
end
figure; hold on
plot(state(1,:), ma_Ap)
plot(state(1,:), ma_Rp)
legend('A promoter', 'R promoter')
xlim([0 t_end])

%% auxiliary function to calculate cumulative value of probability
function c = cumu(m_l)
    % calculate the cumulative of a list
    c = zeros(1, size(m_l,2)+1);
    % m_l = m_l / sum(m_l);
    c(1) = 0;
    for i = 1:size(m_l,2)
        c(i+1) = sum(m_l(1:i));
    end
end
