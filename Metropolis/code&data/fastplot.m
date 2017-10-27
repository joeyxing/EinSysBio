%% Assignment 2: Metropolis algorithm
clear
close all

global stdev ReprData1718
stdev = 1.00;
%                t    m1     m2    m3
ReprData1718 = [4.0  38.3   3.9   8.6   
                8.0   3.6  52.1  10.8   
                12.0  29.5   4.3  23.6   
                16.0   4.8  64.5   9.3   
                20.0  26.6   4.4  31.4   
                24.0   6.7  66.0   7.7   
                28.0  23.7   3.6  34.1];
%% Full sampling marginal probability
% follow the procedure on page 16 :-)
N = 40;
min_alpha0 = 0.85; max_alpha0 = 1.15; N_alpha0 = N;
min_n = 1.9; max_n = 2.1; N_n = N;
min_beta = 4; max_beta = 6; N_beta = N;
min_alpha = 850; max_alpha = 1100; N_alpha = N;

N = 40;
load('FullSamplingResult.mat')
dim_names = {'\alpha_{0}','n','\beta','\alpha'};
h_mtrx = [[min_alpha0,max_alpha0,N_alpha0]
             [min_n,max_n,N_n]
             [min_beta,max_beta,N_beta]
             [min_alpha,max_alpha,N_alpha]];

PX(1,:) = reshape(sum(sum(sum(PP,2),3),4), [N 1]);
PX(2,:) = reshape(sum(sum(sum(PP,1),3),4), [N 1]);
PX(3,:) = reshape(sum(sum(sum(PP,1),2),4), [N 1]);
PX(4,:) = reshape(sum(sum(sum(PP,1),2),3), [N 1]);

i_figure = 0;
for ii = 1:4
    i_figure = i_figure + 1;
    figure(i_figure)
    subplot(1,2,1)
    plot(linspace(h_mtrx(ii,1),h_mtrx(ii,2),h_mtrx(ii,3)), PX(1,:));
    title(['Marginal probability on ' dim_names{ii} ' (Bayesian)'])
    xlabel(dim_names{ii})
    ylabel('Probability')
end

% 2-D marginal probability
PXX = zeros(N,N,6);
dim_cmb = nchoosek(1:4,2);

for ii = 1:size(dim_cmb,1)
    compensate = setdiff(1:4, dim_cmb(ii,:));
    PXX(:,:,ii) = reshape(sum(sum(PP,compensate(1)), compensate(2)), [N N]);
    i_figure = i_figure+1;
    figure(i_figure)
    subplot(1,2,1)
    [g1,g2]=meshgrid(linspace(h_mtrx(dim_cmb(ii,1),1),...
                              h_mtrx(dim_cmb(ii,1),2),...
                              h_mtrx(dim_cmb(ii,1),3)),...
                     linspace(h_mtrx(dim_cmb(ii,2),1),...
                              h_mtrx(dim_cmb(ii,2),2),...
                              h_mtrx(dim_cmb(ii,2),3)));
    surf(g1,g2,PXX(:,:,ii))
    title(['Marginal probability on ' dim_names{dim_cmb(ii,1)} ' and ' dim_names{dim_cmb(ii,2)} ' (Bayesian)']);
    xlabel(dim_names{dim_cmb(ii,1)})
    ylabel(dim_names{dim_cmb(ii,2)})
    zlabel('Probability')
end

%% Part 6 metropolis sampling visualize RR
load('MetropolisResult.mat')

n_steps = 3e5;
i_figure = 0;

% 1-D
for ii = 1:4
    i_figure = i_figure + 1;
    figure(i_figure)
    subplot(1,2,2)
    [v,b] = histcounts(RR(:,ii),N);
    b = (b(1:end-1)+b(2:end))/2;
    plot(b,v/n_steps);
    title(['Marginal probability on ' dim_names{ii} ' (Metropolis)'])
    xlabel(dim_names{ii})
    ylabel('Probability')
end

% 2-D
for ii = 1:size(dim_cmb,1)
    %compensate = setdiff(1:4, dim_cmb(ii,:));
    i_figure = i_figure+1;
    figure(i_figure)
    subplot(1,2,2)
    [v,bx,by] = histcounts2(RR(:,dim_cmb(ii,1)),RR(:,dim_cmb(ii,2)),N);
    bx = (bx(1:end-1)+bx(2:end))/2;
    by = (by(1:end-1)+by(2:end))/2;
    
    surf(bx,by,v/n_steps)
    title(['Marginal probability on ' dim_names{dim_cmb(ii,1)} ' and ' dim_names{dim_cmb(ii,2)} ' (Metropolis)']);
    xlabel(dim_names{dim_cmb(ii,1)})
    ylabel(dim_names{dim_cmb(ii,2)})
    zlabel('Probability')
end
