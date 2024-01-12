path_kids = '/Users/biernack/Documents/MATLAB/grampa-master-12Jul21/TOPKIDS/';%'../';          % path to the complete TOPKIDS data set
track = 'low resolution/';  % low or high resolution
% pair of graphs to be plotted
i=7;j=4; %or select a different pair...
%% prepare seed matching (grampa) and ppmgm matching 
%we load the precomputed results for fig.7. For a different pair, uncomment
%the code below
PPM_output=load('Data_fig7.mat').S;%ppmm result
Seed= transpose(load('Data_fig7.mat').P_sp);%grampa baseline
% %Code for using a different pair of shapes (different from i=7,j=4) below...
% path_to_i= strcat(path_kids,track, 'kid', num2str(i,'%02d'), '.off');
% path_to_j= strcat(path_kids,track, 'kid', num2str(j,'%02d'), '.off'));
% [PPM_output,~,Seed,~,~,~]=matching_shapes_grampa_ppmgm(path_to_i,path_to_j);

%% plot results
M = load_off(strcat(path_kids,track, 'kid', num2str(i,'%02d'), '.off'));
N = load_off(strcat(path_kids,track, 'kid', num2str(j,'%02d'), '.off'));

V1=M.VERT;                  % 3-d coordinates of vertices
F1=M.TRIV;                  % face for triangulation
V2=N.VERT;
F2=N.TRIV;

adj1 = triangulation2adjacency(F1,V1');     % adj after triangulation
adj2 = triangulation2adjacency(F2,V2');

n1=size(adj1,1);
n2=size(adj2,1);

if n1<=n2
    final_idx1 = [1:n1]';
    idx2=1:n2;
     final_idx2 = Seed*idx2';  
     final_idx3 =  PPM_output*idx2';
         
        PlotResultAfterLocalMinimization(V1',F1',V2',F2',final_idx1,final_idx2,final_idx3 ,'Image A','Image B/Grampa','Image B/Grampa+PPM','Comparison weighted');
else
    final_idx1= [1:n2]';
    idx1=1:n1;
    final_idx2 = Seed'*idx1'; 
    final_idx3 = PPM_output'*idx1'; 
        PlotResultAfterLocalMinimization(V2',F2',V1',F1',final_idx1,final_idx2,final_idx3 ,'Image A','Image B/Grampa','Image B/Grampa+PPM','Comparison weighted');
end