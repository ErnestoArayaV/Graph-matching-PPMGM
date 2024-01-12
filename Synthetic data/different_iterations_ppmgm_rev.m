%
%PPM for graph matching
%PPMGM with different number of iterations (parameter N)
%code for generating Fig.4
%
% We generate and match matrices A,B generated from a correlated random
% graph model

%coded by Ernesto Araya Valdivia.  
%%
clc; clear;

%% initialization
num_run = 25; %number of Montecarlo runs 
vec_dim = 1000;%500:500:1500;%50:200:850; %size of the matrix (can be an array)

len_dim = length(vec_dim);
frac_of_fps=0.1;

in_ball=floor((1-frac_of_fps)*vec_dim);

vec_noise =0:0.05:1; %'noise' parameter sigma
len_noise = length(vec_noise);

pi_corr = zeros(len_dim, len_noise, num_run);
pi2_corr = zeros(len_dim, len_noise, num_run);
pi3_corr = zeros(len_dim, len_noise, num_run);
pi4_corr = zeros(len_dim, len_noise, num_run);
pi5_corr = zeros(len_dim, len_noise, num_run);

tic
%% Iteration over independent samples 
for ind_run = 1:num_run

    %% Iteration over dimensions 
    for ind_dim = 1:len_dim
        n = vec_dim(ind_dim);
        
        fprintf('Iteration %i',ind_run);
        %% Iteration over noise levels 
        for ind_noise = 1:len_noise
            
            %Generate random A,B according to a model
            sigma = vec_noise(ind_noise); %disp(sigma);
            fprintf('Iteration %i , Matrix dimension %i , Sigma %f \n', ind_run, n, sigma);
            [A, B, P_rnd] = generate_wig(n,sigma);
            
            %% PPMGM
            %it1
            P_init =initial_perm(n,in_ball,P_rnd);
            P = matching_proj_it_in(A, B,1,P_init);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            pi_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            %it2
            P_init =initial_perm(n,in_ball,P_rnd);
            P = matching_proj_it_in(A, B,2,P_init);
            fix_pt_ratio2 = sum(dot(P_rnd, P)) / n;
            pi2_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio2;
            %it3
            P_init =initial_perm(n,in_ball,P_rnd);
            P = matching_proj_it_in(A, B,4,P_init);
            fix_pt_ratio3 = sum(dot(P_rnd, P)) / n;
            pi3_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio3;
            %it4
            P_init =initial_perm(n,in_ball,P_rnd);
            P = matching_proj_it_in(A, B,8,P_init);
            fix_pt_ratio4 = sum(dot(P_rnd, P)) / n;
            pi4_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio4;
            %it5
            P_init =initial_perm(n,in_ball,P_rnd);
            P = matching_proj_it_in(A, B,30,P_init);
            fix_pt_ratio5 = sum(dot(P_rnd, P)) / n;
            pi5_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio5;
            
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
m_pi5_corr=mean(pi5_corr,3);
%quantiles
p1=0.95;
p2=0.05;
q1_pi_corr=quantile(pi_corr,p1,3);
q1_pi2_corr=quantile(pi2_corr,p1,3);
q1_pi3_corr=quantile(pi3_corr,p1,3);
q1_pi4_corr=quantile(pi4_corr,p1,3);
q1_pi5_corr=quantile(pi5_corr,p1,3);
%
q2_pi_corr=quantile(pi_corr,p2,3);
q2_pi2_corr=quantile(pi2_corr,p2,3);
q2_pi3_corr=quantile(pi3_corr,p2,3);
q2_pi4_corr=quantile(pi4_corr,p2,3);
q2_pi5_corr=quantile(pi5_corr,p2,3);
% 

for j=1:len_dim
    figure;hold on;
    %set data
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
    fill(hdata3_1, inBetween, [0 0.5 0],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    curve41 = q1_pi4_corr (j,:);
    curve42 = q2_pi4_corr (j,:);
    hdata4_1 = [vec_noise, fliplr(vec_noise)];
    inBetween = [curve41, fliplr(curve42)];
    fill(hdata4_1, inBetween, [.5 0 .5],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    curve51 = q1_pi5_corr (j,:);
    curve52 = q2_pi5_corr (j,:);
    hdata5_1 = [vec_noise, fliplr(vec_noise)];
    inBetween = [curve51, fliplr(curve52)];
    fill(hdata5_1, inBetween, [0.9290 0.6940 0.1250],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    %set curves
    set(hdata1, 'LineStyle','--','Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.75 .75 1],'LineWidth', 2);
    set(hdata2, 'LineStyle','--','Color', [0.75 0 0],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.75 0 0],'LineWidth', 2);
    set(hdata3, 'LineStyle','--','Color',[0 0.5 0],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0 0.5 0],'LineWidth', 2);
    set(hdata4, 'LineStyle','--','Color', [.5 0 .5],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', [.5 0 .5], 'MarkerFaceColor',  [.5 0 .5],'LineWidth', 2);
    set(hdata5, 'LineStyle','--','Color', [0.9290 0.6940 0.1250],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', [0.9290 0.6940 0.1250], 'MarkerFaceColor',[0.9290 0.6940 0.1250],'LineWidth', 2);
    %title
    str = sprintf('PPMGM with n=%i, Wigner',vec_dim(j));
    hTitle=title(str);set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
    %labels
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on','XColor', [.4 .4 .4], 'YColor', [.4 .4 .4],'LineWidth', 1,'Fontsize',14);
    hXLabel = xlabel('Noise $\sigma$','interpreter','latex');
    hYLabel = ylabel('Recovery fraction','interpreter','latex');
    set([hXLabel, hYLabel], 'FontSize', 20);
    %legend
    hLegend = legend([hdata1, hdata2,hdata3,hdata4,hdata5], 'PPMGM it.1','PPMGM it.2','PPMGM it.3', 'PPMGM it.4','PPMGM it.5');
    set(hLegend, 'FontSize', 11);hold off
    
end

