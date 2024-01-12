%
%PPM for graph matching
%PPMGM used to test sparsification strategies for graph matching methods
%code for generating Fig.5 (a and b)
%
% We generate and match matrices A,B generated from a correlated random
% graph model, then sparsify them according to three strategies:
% Spar1,Spar2, Spar3. We refer to Sec.5.2 of "Seeded graph matching for the correlated Gaussian Wigner
% model via the projected power method" for the details. 
% 

%coded by Ernesto Araya Valdivia. 
%%
clc; clear;

%% initialization
maxiter=5;
num_run = 25; %number of Montecarlo runs 
vec_dim = 1000;%500:500:1500;%50:200:850; %size of the matrix (can be an array)
len_dim = length(vec_dim);
frac_of_fps=0.1;
len_fps=length(frac_of_fps);
in_ball=floor((1-frac_of_fps)*vec_dim);

vec_noise =0:0.05:1; %'noise' parameter sigma
len_noise = length(vec_noise);

pi_corr = zeros(len_dim, len_noise, num_run);%, len_th);         %measures of accuracy 
pi2_corr = zeros(len_dim, len_noise, num_run);%,len_th);
pi3_corr = zeros(len_dim, len_noise, num_run);%,len_th);
pi4_corr = zeros(len_dim, len_noise, num_run);%,len_th);

tic
%% Iteration over independent samples 
for ind_run = 1:num_run

    %% Iteration over dimensions 
    for ind_dim = 1:len_dim
        n = vec_dim(ind_dim);
        ss=log(n)/(2*n);
        t=(sqrt(2)*erfcinv(2*ss))/sqrt(n);
        fprintf('Iteration %i',ind_run);
        %% Iteration over noise levels 
        for ind_noise = 1:len_noise
            
            %Generate random A,B according to a model
            sigma = vec_noise(ind_noise); %disp(sigma);
            fprintf('Iteration %i , Matrix dimension %i , Sigma %f \n', ind_run, n, sigma);
            [A, B, P_rnd] = generate_wig(n,sigma);
            %prepare the copies of 'A and B' for thresholding
            A2=A-diag(diag(A));B2=B-diag(diag(B));A3=A;B3=B;
            %first thresholding (t1), from the 'mother graphs' A,B
            A2(abs(A2)>t)  = 1;
            A2(A2~=1) = 0;
            B2(abs(B2)>t)  = 1;
            B2(B2~=1) = 0;
            %second thresholding (t2), from the 'mother graphs' A,B
            A3(abs(A)<t)  = 0;
            B3(abs(B)<t)  = 0;
            %top k per row (t3)
            k=ceil(2*log(n));
            A4=topk_perRow(A,k);
            B4=topk_perRow(B,k);
            
            %% PPMGM
            % without thresholding
            P_init =initial_perm(n,in_ball,P_rnd);
            P = matching_proj_it_in(A, B,maxiter,P_init);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            pi_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            % with sparsification strat. 1
            P = matching_proj_it_in(A2, B2,maxiter,P_init);
            fix_pt_ratio2 = sum(dot(P_rnd, P)) / n;
            pi2_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio2;
            % with sparsification strat. 2
            P = matching_proj_it_in(A3, B3,maxiter,P_init);
            fix_pt_ratio3 = sum(dot(P_rnd, P)) / n;
            pi3_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio3;
            % with sparsification strat. 3
            P = matching_proj_it_in(A4, B4,maxiter,P_init);
            fix_pt_ratio4 = sum(dot(P_rnd, P)) / n;
            pi4_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio4;

            toc;
           
        end
     %   end
    end
end
%mean
m_pi_corr=mean(pi_corr,3);
m_pi2_corr=mean(pi2_corr,3);
m_pi3_corr=mean(pi3_corr,3);
m_pi4_corr=mean(pi4_corr,3);

%quantiles
p1=0.95;
p2=0.05;
q1_pi_corr=quantile(pi_corr,p1,3);
q1_pi2_corr=quantile(pi2_corr,p1,3);
q1_pi3_corr=quantile(pi3_corr,p1,3);
q1_pi4_corr=quantile(pi4_corr,p1,3);
%
q2_pi_corr=quantile(pi_corr,p2,3);
q2_pi2_corr=quantile(pi2_corr,p2,3);
q2_pi3_corr=quantile(pi3_corr,p2,3);
q2_pi4_corr=quantile(pi4_corr,p2,3);
% 

for j=1:len_dim
    figure;hold on;
    %set data
    hdata1=line(vec_noise, m_pi_corr(j,:));
    hdata2=line(vec_noise, m_pi2_corr(j,:));
    hdata3=line(vec_noise, m_pi3_corr(j,:));
    hdata4=line(vec_noise, m_pi4_corr(j,:));
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
    fill(hdata3_1, inBetween, [0 0.5 0],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    curve41 = q1_pi4_corr (j,:);
    curve42 = q2_pi4_corr (j,:);
    hdata4_1 = [vec_noise, fliplr(vec_noise)];
    inBetween = [curve41, fliplr(curve42)];
    fill(hdata4_1, inBetween, [.5 0 .5],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    %set curves
    set(hdata1, 'LineStyle','--','Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.75 .75 1],'LineWidth', 2);
    set(hdata2, 'LineStyle','--','Color', [0.75 0 0],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.75 0 0],'LineWidth', 2);
    set(hdata3, 'LineStyle','--','Color',[0 0.5 0],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0 0.5 0],'LineWidth', 2);
    set(hdata4, 'LineStyle','--','Color', [.5 0 .5],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', [.5 0 .5], 'MarkerFaceColor',  [.5 0 .5],'LineWidth', 2);
    %title
    str = sprintf('PPMGM with n=%i, Wigner',vec_dim(j));
    hTitle=title(str);set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
    %labels
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on','XColor', [.4 .4 .4], 'YColor', [.4 .4 .4],'LineWidth', 1,'Fontsize',14);
    hXLabel = xlabel('Noise $\sigma$','interpreter','latex');
    hYLabel = ylabel('Recovery fraction','interpreter','latex');
    set([hXLabel, hYLabel], 'FontSize', 20);
    %legend
    hLegend = legend([hdata1, hdata2,hdata3,hdata4], 'PPMGM','PPMGM thr.1','PPMGM thr.2', 'PPMGM topk');
    set(hLegend, 'FontSize', 11);hold off
    
end

