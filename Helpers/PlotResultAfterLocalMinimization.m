function PlotResultAfterLocalMinimization(V1,F1,V2,F2,idx1,idx2,idx3,title1,title2,title3,title4)
%figure
colorVerts = V1(:,idx1)';
scattColor = bsxfun(@rdivide,bsxfun(@minus,colorVerts,min(colorVerts,[],1)), (max(colorVerts,[],1)-min(colorVerts,[],1)));
%figure('Position', [100, 100, 1000, 1000]);
%figure;
f=figure('units','normalized','position',[0 0 0.5 0.5]);
f.Name = title4;
%[ha, pos]=tight_subplot(1,2,0,0.01,0.01);
subplot(1,3,1)
title(title1,'FontSize',20)
params.scattColor = scattColor;
params.verInd  = idx1;
plotMeshAndPoints( V1, F1, params )
%axes(ha(1));
subplot(1,3,2)
title(title2,'FontSize',20)
params.scattColor = scattColor;
params.verInd  = idx2;
plotMeshAndPoints( V2, F2, params )
subplot(1,3,3)
title(title3,'FontSize',20)
params.scattColor = scattColor;
params.verInd  = idx3;
plotMeshAndPoints( V2, F2, params )
%axes(ha(2));
%set(gcf,'Position',[100, 100, 800, 800]);
%pos = get(gca, 'Position');
   % pos(1) = 0.4;
   % set(gca, 'Position', pos)
   % set(gcf,'color','w');
    %a = axes;
%t1 = title(title3);
%t1.Visible='on';
end