%
%PPMGM used to boost seedless graph matching methods: grampa and full
%quadratic programming 

%code for generating Fig.3 

%matching_robust_spectral.m and matching_full_qp.m were taken from from
%https://github.com/Leron33/grampa authored by Jiaming Xu. 
%for quadprog_admm.m refer to http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
%
%% initialization
clear;

num_run = 25;
vec_dim = 200;%500:10:500;
len_dim = length(vec_dim);
vec_noise = 0:0.05:0.8;
len_noise = length(vec_noise);

full_qp_corr = ones(len_dim, len_noise, num_run);
full_qp_run = zeros(len_dim, len_noise, num_run);
qm_PPM_corr=ones(len_dim, len_noise, num_run);
qm_PPM_run=zeros(len_dim, len_noise, num_run);
robust_corr = ones(len_dim, len_noise, num_run);
robust_run = zeros(len_dim, len_noise, num_run);
gp_PPM_corr=ones(len_dim, len_noise, num_run);
gp_PPM_run=zeros(len_dim, len_noise, num_run);

%% Iteration over independent samples 
for ind_run = 1:num_run
    fprintf('Iteration %i \n', ind_run);

    %% Iteration over dimensions 
    for ind_dim = 1:len_dim
        n = vec_dim(ind_dim);
        fprintf('Matrix dimension %i \n', n);
        
        %% Iteration over noise levels 
        for ind_noise = 3:len_noise 
            sigma = vec_noise(ind_noise); disp(sigma);
            [A, B, P_rnd] = generate_wig(n,sigma); %CWG model
            % one can put here different models
            %[A, B, A0, B0, P_rnd] = generate_er(n, p, sigma);%CER model
            
            %% Full QP 
            tic;
            P = matching_full_qp(A, B);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            full_qp_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            full_qp_run(ind_dim, ind_noise, ind_run) = toc;
            %% Full QP+PPMGM
            P=matching_proj_it_in(A,B,5,P);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            qm_PPM_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            qm_PPM_run(ind_dim, ind_noise, ind_run) = toc;

            %% Grampa
            tic;
            P = matching_robust_spectral(A, B, 0.2);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            robust_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            robust_run(ind_dim, ind_noise, ind_run) = toc;
            %% Grampa + PPMGM 
            P=matching_proj_it_in(A,B,5,P);
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            gp_PPM_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            gp_PPM_run(ind_dim, ind_noise, ind_run) = toc;
        end
    end
end
%mean correlation 
full_qp_corr_mean = mean(full_qp_corr, 3);
qm_PPM_corr_mean=mean(qm_PPM_corr,3);
robust_corr_mean = mean(robust_corr, 3);
gp_PPM_corr_mean=mean(gp_PPM_corr,3);

full_qp_run_mean = mean(mean(mean(full_qp_run(1, 3:end, :))));
qm_PPM_run_mean=mean(mean(mean(qm_PPM_run(1, 3:end, :))));
robust_run_mean = mean(mean(mean(robust_run(1, 3:end, :))));
gp_PPM_run_mean=mean(mean(mean(gp_PPM_run(1, 3:end, :))));

%quantiles for 'confidence regions'
p1=0.95;
p2=0.05;
q1_pi_corr=quantile(full_qp_corr,p1,3);
q1_pi2_corr=quantile(qm_PPM_corr,p1,3);
q1_pi3_corr=quantile(robust_corr,p1,3);
q1_pi4_corr=quantile(gp_PPM_corr,p1,3);

q2_pi_corr=quantile(full_qp_corr,p2,3);
q2_pi2_corr=quantile(qm_PPM_corr,p2,3);
q2_pi3_corr=quantile(robust_corr,p2,3);
q2_pi4_corr=quantile(gp_PPM_corr,p2,3);

disp([robust_run_mean,qm_PPM_corr_mean, full_qp_run_mean]);

for j=1:len_dim
    figure;hold on;
    %set data
    hdata3=line(vec_noise, full_qp_corr_mean(j,:));
    hdata4=line(vec_noise, qm_PPM_corr_mean(j,:));
    hdata1=line(vec_noise, robust_corr_mean(j,:));
    hdata2=line(vec_noise, gp_PPM_corr_mean(j,:));
    %define shaded areas
    curve11 = q1_pi_corr (j,:);
    curve12 = q2_pi_corr (j,:);
    hdata1_1 = [vec_noise, fliplr(vec_noise)];
    inBetween = [curve11, fliplr(curve12)];
    fill(hdata1_1, inBetween, [0.75 0 0],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
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
    fill(hdata4_1, inBetween, [.75 .75 1],'FaceAlpha',0.2,'LineStyle',':','EdgeColor','None');
    %set curves
    set(hdata1, 'LineStyle','--','Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.75 .75 1],'LineWidth', 2);
    set(hdata3, 'LineStyle','--','Color', [0.75 0 0],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', 'r', 'MarkerFaceColor', [0.75 0 0],'LineWidth', 2);
    set(hdata2, 'LineStyle','-','Color', [.75 .75 1],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', [.75 .75 1], 'MarkerFaceColor', [.75 .75 1],'LineWidth', 2);
    set(hdata4, 'LineStyle','-','Color', [0.75 0 0],'Marker', 'none', 'MarkerSize', 5,'MarkerEdgeColor', [0.75 0 0], 'MarkerFaceColor', [0.75 0 0],'LineWidth', 2);
    %title
    str = sprintf('PPMGM with n=%i, Wigner',vec_dim(j));
    hTitle=title(str);set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
    %labels
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on','XColor', [.4 .4 .4], 'YColor', [.4 .4 .4],'LineWidth', 1,'Fontsize',14);
    hXLabel = xlabel('Noise $\sigma$','interpreter','latex');
    hYLabel = ylabel('Recovery fraction','interpreter','latex');
    set([hXLabel, hYLabel], 'FontSize', 20);
    %legend
    hLegend = legend([hdata1, hdata2,hdata3,hdata4], 'Grampa', 'Grampa+PPMGM','QPADMM','QPADMM+PPMGM');%,'PPMGM in.5','PPMGM in.6');%,'PPMGM in.5','PPMGM in.6');%,'Top-eigen','Top-eigen+PPM');
    set(hLegend, 'FontSize', 11);hold off
end
%save('.\mat_files\comparison_sp_dp_qp.mat');