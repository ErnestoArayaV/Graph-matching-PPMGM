%
%PPM for graph matching
%PPMGM used to boost seedless graph matching methods
%code for generating Fig.2a
%
% We generate and match matrices A,B generated from a correlated random
% graph model

%coded by Ernesto Araya Valdivia. 
%matching_robust_spectral.m was taken from from
%https://github.com/Leron33/grampa authored by Jiaming Xu. 
%%
clc; clear;

%% initialization
maxiter=5;
num_run = 25; %number of Montecarlo runs 

vec_dim = 800; %500:500:1500;%50:200:850; %size of the A,B matrices (can be an array)
len_dim = length(vec_dim);

vec_noise =0:0.05:1; %vector of 'noise' parameters sigma
len_noise = length(vec_noise);

pi_corr = zeros(len_dim, len_noise, num_run); %to store the correlations with the ground truth matching 
pi2_corr = zeros(len_dim, len_noise, num_run);
pi3_corr = zeros(len_dim, len_noise, num_run);
pi4_corr = zeros(len_dim, len_noise, num_run);
pi5_corr = zeros(len_dim, len_noise, num_run);

tic
%% Iteration over independent samples 
for ind_run = 1:num_run
fprintf('Iteration %i',ind_run);
    %% Iteration over dimensions 
    for ind_dim = 1:len_dim
        n = vec_dim(ind_dim);
        %% Iteration over noise levels 
        for ind_noise = 1:len_noise
            %Generate random A,B according to a model
            sigma = vec_noise(ind_noise); 
            fprintf('Iteration %i , Matrix dimension %i , Sigma %f \n', ind_run, n, sigma);
            [A, B, P_rnd] = generate_wig(n,sigma); %CGW model
            %[A, B, A0, B0, P_rnd] = generate_er(n, p, sigma); %CER model
            %other models below (insert here new model if desired)
            %c1=n/2;
            %c2=n/2;
            %[A,B,A0,B0,P_rnd]=generate_SMB_2(c1,c2, p,q,sigma);%cor. SBM
            %[A, B, P_rnd] = generate_wish(n,d,1,sigma); %cor. Wishart mat.
           
           %different algorithms 
            %% Grampa
            
            P= matching_robust_spectral(A, B, 0.2);
            %% Grampa+PPMGM
            
            P3= matching_proj_it_in(A, B,maxiter,P);
            %% Umeyama
            
            P2=matching_umeyama(A, B);
            %% Umeyama+PPMGM
            
            P4= matching_proj_it_in(A,B,maxiter^2,P2);
            %% PPM random init
            %choose at random a number (in_ball param.) of columns to permute of the ground truth
            %perm. 
            frac_of_fps=1/10;% the fraction of fixed points(where P_rnd does not change) reatined in the final permutation
            in_ball=floor((1-frac_of_fps)*n); %number of permuted columns
            P_init =initial_perm(n,in_ball,P_rnd); %Initialization matrix with frac_of_fps*n fixed points
            P5 = matching_proj_it_in(A, B,maxiter,P_init);
           
            %% compute accuracy measures (correlation)
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            fix_pt_ratio2= sum(dot(P_rnd, P2)) / n;
            fix_pt_ratio3= sum(dot(P_rnd, P3)) / n;
            fix_pt_ratio4= sum(dot(P_rnd, P4)) / n;
            fix_pt_ratio5= sum(dot(P_rnd, P5)) / n;
           
            pi_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            pi2_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio2;
            pi3_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio3;
            pi4_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio4;
            pi5_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio5;
  
            toc;
    
        end
    end
end
%compute statistics for plotting
%mean
m_pi_corr=mean(pi_corr,3);
m_pi2_corr=mean(pi2_corr,3);
m_pi3_corr=mean(pi3_corr,3);
m_pi4_corr=mean(pi4_corr,3);
m_pi5_corr=mean(pi5_corr,3);
%std
std_pi_corr=std(pi_corr,0,3);
std_pi2_corr=std(pi2_corr,0,3);
std_pi3_corr=std(pi3_corr,0,3);
std_pi4_corr=std(pi4_corr,0,3);
std_pi5_corr=std(pi5_corr,0,3);
%quantiles
p1=0.95; %define the quantiles for 'confidence' regions
p2=0.05;

q1_pi_corr=quantile(pi_corr,p1,3);
q1_pi2_corr=quantile(pi2_corr,p1,3);
q1_pi3_corr=quantile(pi3_corr,p1,3);
q1_pi4_corr=quantile(pi4_corr,p1,3);
q1_pi5_corr=quantile(pi5_corr,p1,3);

q2_pi_corr=quantile(pi_corr,p2,3);
q2_pi2_corr=quantile(pi2_corr,p2,3);
q2_pi3_corr=quantile(pi3_corr,p2,3);
q2_pi4_corr=quantile(pi4_corr,p2,3);
q2_pi5_corr=quantile(pi5_corr,p2,3);
%mins
min_pi_corr=min(pi_corr,3);
min_pi2_corr=min(pi2_corr,3);
min_pi3_corr=min(pi3_corr,3);
min_pi4_corr=min(pi4_corr,3);
min_pi5_corr=min(pi5_corr,3);
%max
max_pi_corr=max(pi_corr,3);
max_pi2_corr=max(pi2_corr,3);
max_pi3_corr=max(pi3_corr,3);
max_pi4_corr=max(pi4_corr,3);
max_pi5_corr=max(pi5_corr,3);

%plotting loop
for j=1:len_dim
    figure;hold on;
    %data for the curves
    hdata1=line(vec_noise, m_pi_corr(j,:));
    hdata2=line(vec_noise, m_pi2_corr(j,:));
    hdata3=line(vec_noise, m_pi3_corr(j,:));
    hdata4=line(vec_noise, m_pi4_corr(j,:));
    hdata5=line(vec_noise, m_pi5_corr(j,:));
    %define shaded areas
    curve11 = q1_pi_corr (j,:);
    curve12 = q2_pi_corr (j,:);
    hdata1_1 = [vec_noise, fliplr(vec_noise)];
    inBetween = [curve11, fliplr(curve12)];
    fill(hdata1_1, inBetween, [.75 .75 1],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    curve21 = q1_pi2_corr (j,:);
    curve22 = q2_pi2_corr (j,:);
    hdata2_1 = [vec_noise, fliplr(vec_noise)];
    inBetween = [curve21, fliplr(curve22)];
    fill(hdata2_1, inBetween, [0.75 0 0],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    curve31 = q1_pi3_corr (j,:);
    curve32 = q2_pi3_corr (j,:);
    hdata3_1 = [vec_noise, fliplr(vec_noise)];
    inBetween = [curve31, fliplr(curve32)];
    fill(hdata3_1, inBetween, [.75 .75 1],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    curve41 = q1_pi4_corr (j,:);
    curve42 = q2_pi4_corr (j,:);
    hdata4_1 = [vec_noise, fliplr(vec_noise)];
    inBetween = [curve41, fliplr(curve42)];
    fill(hdata4_1, inBetween, [0.75 0 0],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    curve51 = q1_pi5_corr (j,:);
    curve52 = q2_pi5_corr (j,:);
    hdata5_1 = [vec_noise, fliplr(vec_noise)];
    inBetween = [curve51, fliplr(curve52)];
    fill(hdata5_1, inBetween, [0 0.5 0],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    %set the curves
    set(hdata1, 'LineStyle','--','Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.75 .75 1],'LineWidth', 2);
    set(hdata2, 'LineStyle','--','Color', 'r','Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.75 0 0],'LineWidth', 2);
    set(hdata3, 'LineStyle','-','Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.75 .75 1],'LineWidth', 2);
    set(hdata4, 'LineStyle','-','Color', 'r','Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', [0.75 0 0], 'MarkerFaceColor',  [0.75 0 0],'LineWidth', 2);
    set(hdata5, 'LineStyle','--','Color', [0 0.5 0],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', [0 0.5 0], 'MarkerFaceColor',[0 0.5 0],'LineWidth', 2);
    %title
    str = sprintf('PPMGM with n=%i, Wigner',vec_dim(j));
    hTitle=title(str);set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on','XColor', [.4 .4 .4], 'YColor', [.4 .4 .4],'LineWidth', 1,'Fontsize',14);
    %labels
    hXLabel = xlabel('Noise $\sigma$','interpreter','latex');
    hYLabel = ylabel('Recovery fraction','interpreter','latex');
    set([hXLabel, hYLabel], 'FontSize', 20);
    %legends
    hLegend = legend([hdata1, hdata2,hdata3,hdata4,hdata5], 'Grampa','Umeyama','Grampa+PPMGM', 'Umeyama+PPMGM','PPMGM rand. init.');  
    set(hLegend, 'FontSize', 11);hold off
    
end

