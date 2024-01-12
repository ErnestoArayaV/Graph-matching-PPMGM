function [S_100,S_30,P_sp,curve1,curve2,curve_grampa]=matching_shapes_grampa_ppmgm(path_i,path_j)

thresholds = 0:0.01:1;
%initialize the CDF curves
curve1 = zeros(nb_examples,length(thresholds));%CDF for ppmgm(N=100)
curve2 = zeros(nb_examples,length(thresholds));%CDF for ppmgm(N=30)
curve_grampa = zeros(nb_examples,length(thresholds));%CDF for grampa
%% load data
    % We load to shapes M, N to be matched 
    M = load_off(path_i);
    N = load_off(path_j);
    
    V1=M.VERT;                  % 3-d coordinates of vertices
    F1=M.TRIV;                  % face for triangulation
    V2=N.VERT;
    F2=N.TRIV;
    
    %triangulation to adjacency
    adj1 = triangulation2adjacency(F1,V1);     % adj. after triangulation
    adj2 = triangulation2adjacency(F2,V2);
    dist=geodesic_distance(N.TRIV,N.VERT);      %Added by JX
    dist=sparse(dist);
    n1=size(adj1,1);
    n2=size(adj2,1);
    
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
        curve_grampa(k,m) = 100*sum(errors_grampa <= thresholds(m)) / length(errors_init);
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