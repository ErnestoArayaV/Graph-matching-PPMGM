% code to obtain the heatmap in fig. 6
% we load data obtained with the function thresholding_experiments_rev.m

%% load data 
vec_noise = load('Data_fig6.mat').vec_noise;
vec_th = load('Data_fig6.mat').vec_th;
Mat = load('Data_fig6.mat').Mat;
%plot heatmap
selected_colormap = copper;%change the colormap if desired
figure;
h=heatmap(vec_noise,vec_th,Mat','Colormap',selected_colormap);
h.XLabel='Noise';
h.YLabel='Threshold';

