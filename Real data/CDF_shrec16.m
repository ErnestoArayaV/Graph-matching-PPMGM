%
%Code for testing ppmgm with the SHREC'16 computer vision dataset. cf https://cvg.cit.tum.de/_media/spezial/bib/shrec16-topology-slides.pdf
%Computes the cdf curves for Grampa and ppmgm (n=30) and ppmgm(N=100)
%imporvements over the grqmpa baseline. 

%It requires the datafile TOPKIDS. 
%It uses helpers: cal_tri_areas.m, calc_dist_matrix.m,
%geodesic_distance.m, getoptions.m,load_off.m,merge_ground_truth.m,read_corerspondance.m, triangulation2adjacency.m (copyrights in the specific files)

%Based on a code by Zorah Lähner (laehner@in.tum.de). Revised by Jiaming Xu
%Feb. 8 2020. Revised by Ernesto Araya Valdivia 10 Jan.2024

%% PREPARE
path_kids = './Real data/TOPKIDS';          % path to the complete TOPKIDS data set
track = '/low resolution/';  % low or high resolution
%% CALCULATE CURVES
thresholds = 0:0.01:1;
nb_examples = 90;
%initialize the CDF curves
curve1 = zeros(nb_examples,length(thresholds));
curve2 = zeros(nb_examples,length(thresholds));
curve_grampa = zeros(nb_examples,length(thresholds));
for k=1:nb_examples
    i=ceil(k/9);
    j=k-(i-1)*9;
    if j<i
        j=j+15;
    else
        j=j+16;
    end
    i=i+15;
    %% load data
    % We load to shapes M, N to be matched 
    strcat(path_kids,track, 'kid', num2str(i,'%02d'), '.off')
    M = load_off(strcat(path_kids,track, 'kid', num2str(i,'%02d'), '.off'));
    N = load_off(strcat(path_kids,track, 'kid', num2str(j,'%02d'), '.off'));
    
    V1=M.VERT;                  % 3-d coordinates of vertices
    F1=M.TRIV;                  % face for triangulation
    V2=N.VERT;
    F2=N.TRIV;
    
    %triangulation to adjacency
    adj1 = triangulation2adjacency(M.TRIV);     % adj. after triangulation
    adj2 = triangulation2adjacency(N.TRIV);
    dist=geodesic_distance(N.TRIV,N.VERT);      %Added by JX
    dist=sparse(dist);
    n1=size(adj1,1);
    n2=size(adj2,1);
    
%     EYE1=sparse(1:n1,1:n1,1,n1,n1);
%     EYE2=sparse(1:n2,1:n2,1,n2,n2);
%     W21=double((double((adj1*adj1)>0)-adj1-EYE1)>0);
%     W22=double((double((adj2*adj2)>0)-adj2-EYE2)>0);
    %% Read ground truth
    % recover the gt correct matching
    gt_M_null = read_correspondence(strcat(path_kids, track, 'kid', num2str(i,'%02d'), '_ref.txt'));
    gt_N_null = read_correspondence(strcat(path_kids, track, 'kid', num2str(j,'%02d'), '_ref.txt'));
    gt = merge_ground_truth(gt_M_null, gt_N_null);
    P_rnd=zeros(n2,n1);
    for ind=1:length(gt(:,1))
        P_rnd(gt(ind,2),gt(ind,1))=1;
    end
    P_rnd=sparse(P_rnd);

    %% Grampa baseline
    eta =0.2; %as recomended in Grampa paper
    P_sp=matching_robust_spectral(full(adj1), full(adj2), eta);%E.A:in my experiments the full version perfoms better than the sparse one
   
    if n1<=n2
        init_idx1 = [1:size(adj1,1)]';
        idx2=1:size(adj2,1);
        init_idx2 = P_sp*idx2';   %contains nodes in adj2 corresponding to node i in adj1
        %visualization
%            PlotResultAfterLocalMinimization(V1',F1',V2',F2',init_idx1,init_idx2,'source','target');
    else
        init_idx2=[1:size(adj2,1)]';
        idx1=1:size(adj1,1);
        init_idx1=P_sp'*idx1';
%            PlotResultAfterLocalMinimization(V1',F1',V2',F2',init_idx1,init_idx2,'source','target');
    end
    
    %compute the errors for Grampa according to the Princeton protocol
    corr_grampa=[init_idx1,init_idx2];
    errors_grampa = zeros(size(corr_init,1), 1);
    
    for m=1:size(corr_init,1)
        
        if (strcmp(track, 'low resolution/'))
            gt_match = gt(gt(:,1) == corr_init(m,1), 2);
            match = corr_grampa(m,2);
            
            if ~isempty(gt_match)
                % using the geodesic distance of the second graph
                errors_grampa(m) = dist(gt_match, match); % TODO include your geodesics here
            else
                errors_grampa(m) = 200;
            end
        else
            errors_grampa(m) = 100;
        end
        
    end
    
    diameters = sqrt(sum(calc_tri_areas(N)));
    errors_grampa = errors_grampa / diameters;
    for m=1:length(thresholds)
        curve_grampa(k,m) = 100*sum(errors_init <= thresholds(m)) / length(errors_init);
    end
    
    %compute the correlation of Grampa with the ground truth
    corr_sp=full(sum(dot(P_rnd,P_sp'))/length(errors_init));
    
    %% PPMGM (rectangular version)
    max_iter =100; 
    iter = 30;
    r_old=corr_sp;
    S_100=P_sp'; %take the result of Grampa as initial point, for ppmgm with N=100 iters.
    
    for iter_count=1:1:max_iter
        X=adj2*S_100*adj1;
        [S_100,] = greedy_match(X);
        r=full(sum(dot(P_rnd,S_100))/length(errors_init));
%         if abs(r-r_old) <1e-6  % uncomment for tol. convergence criterion 
%             break;
%         end
        r_old=r;
        if iter_count==30
            S_30=S_100; %define the ppmgm with N=30 iterations
        end
    end
    %compute the matching for S_100
    if n1<=n2
        final_idx1 = [1:n1]';
        idx2=1:n2;
        final_idx2 = S_100'*idx2';
        %        PlotResultAfterLocalMinimization(V1',F1',V2',F2',final_idx1,final_idx2,'source','target');
    else
        final_idx2=[1:n2]';
        idx1=1:n1;
        final_idx1=S_100*idx1';
        %        PlotResultAfterLocalMinimization(V1',F1',V2',F2',final_idx1,final_idx2,'source','target');
    end
    corr=[final_idx1,final_idx2];
    
    errors = zeros(size(corr,1), 1);
    
    for m=1:size(corr,1)
        
        if (strcmp(track, 'low resolution/'))
            gt_match = gt(gt(:,1) == corr(m,1), 2);
            match = corr(m,2);
            
            if ~isempty(gt_match)&& match>0
                % using the second graph
                errors(m) = dist(gt_match, match); 
                errors(m) = -2;
            end
        else
            errors(m) = -1;
        end
        
    end
    diameters = sqrt(sum(calc_tri_areas(N)));
    errors = errors / diameters;
    for m=1:length(thresholds)
        curve1(k,m) = 100*sum(errors <= thresholds(m)) / length(errors);
    end
     %compute the matching for S_30
    if n1<=n2
        final_idx1 = [1:n1]';
        idx2=1:n2;
        final_idx2 = S_30'*idx2';
        %        PlotResultAfterLocalMinimization(V1',F1',V2',F2',final_idx1,final_idx2,'source','target');
    else
        final_idx2=[1:n2]';
        idx1=1:n1;
        final_idx1=S_30*idx1';
        %        PlotResultAfterLocalMinimization(V1',F1',V2',F2',final_idx1,final_idx2,'source','target');
    end
    corr2=[final_idx1,final_idx2];
    
    errors2 = zeros(size(corr2,1), 1);
    
    for m=1:size(corr2,1)
        
        if (strcmp(track, 'low resolution/'))
            gt_match = gt(gt(:,1) == corr(m,1), 2);
            match = corr2(m,2);
            
            if ~isempty(gt_match)&& match>0
                % using the second graph
                errors2(m) = dist(gt_match, match); 
                errors2(m) = -2;
            end
        else
            errors2(m) = -1;
        end
        
    end
    diameters = sqrt(sum(calc_tri_areas(N)));
    errors2 = errors2 / diameters;
    for m=1:length(thresholds)
        curve2(k,m) = 100*sum(errors2 <= thresholds(m)) / length(errors2);
    end
    
end
%plot the curves 
figure;hold on;
plot(thresholds', mean(curve1/100, 1)','b'), ylim([0 1]);
plot(thresholds', mean(curve2/100, 1)','color',[0.4660 0.6740 0.1880]), ylim([0 1]);
plot(thresholds', mean(curve_grampa/100, 1)','red'), ylim([0 1]);
legend('PPMGM (N=100)','PPMGM (N=30)','GRAMPA');
line_width=1.5;
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',line_width);
xlabel('${\epsilon}$','interpreter','latex', 'FontSize',25)
ylabel('CDF', 'FontSize',20)