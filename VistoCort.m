function [] = VistoCort(overlays,mapValues,fsData)
% attempt to go from the visual map to the cortical map. Coords are taken
% from the specific subject involved <<DA>>
sub = 1;
if ieNotDefined('fitVals'), fitVals = [1.5 4.5]; end

subData = mapValues(mapValues(:,1)==1,:); %take data from first subject for now (BM)
subData = subData(1:32,:); %only take one eye for now (this is BM's bad eye)
hemi = 'RH';
A = 0.4990;% actual values
B = 0.1845;
% q = overlays(2,sub).q; %first line is bm. [pa ecc]
%plot the visual map with the responses of left and right eye
% subplot(2,2,1:2)
% figure, polar(d2r(subData(:,3)),subData(:,2),'.r') %plots reference data
hold on
subRange = [2 3];
% calculate coords
% bensonECCCoords = @(ecc,q) (q(2) + log(ecc/90))/q(2); %found solution using solve()
% bensonPACoords = @(pa,q) -(pa/180).^(1/q(1));%if all values are real (found solution using solve('Real',1))
nStep = 8;
%coords are in el-coordinates

paRange = [90 270];
switch hemi
    case 'LH'
        inbounds = subData(:,3) < paRange(1) | subData(:,3) > paRange(2);
    case 'RH'
        inbounds = subData(:,3) < paRange(2) & subData(:,3) > paRange(1);
end
subData = subData(inbounds,:);

[~,idx] = findClosestVoxel(subData,overlays,4,5); %finds the voxels/coordinates that have pa/ecc values closest to vis pa/ecc values
[~,idxRef] = findClosestVoxel(subData,overlays,2,3); %reference points

%make a subplot with the same data, but different color overlays!
for t=1:6,
    figure
    set(gca,'FontSize',25,'FontName','Trebuchet')
    set(gca,'XTick',[],'LineWidth',2);
    set(gca,'YTick',[],'LineWidth',2);
    box('on');
    switch t
        case 1
            plotOver = overlays(2,1).values; %pa data
            select = overlays(2,1).restrict;
            scalVals = [-pi pi];
            colormap('hsv');
            title('polar angle data');
        case 2
            plotOver = overlays(2,1).model; %pa model
            select = overlays(2,1).restrict;
            scalVals = [-pi pi];
            colormap('hsv');
            title('polar angle model');
        case 3
            plotOver = overlays(2,1).residuals; %pa residual
            select = overlays(2,1).restrict;
            scalVals = [0 40];
            colormap('Summer');
            title('polar angle residual');
        case 4
            plotOver = overlays(3,1).values; %ecc data
            select = overlays(2,1).restrict; %use pa restriction for now, ecc restriction doesn't leave much there
            scalVals = [0 8];
            colormap('hsv');
            title('eccentricity data');
        case 5
            plotOver = overlays(3,1).model; %ecc model
            select = overlays(2,1).restrict;
            scalVals = [0 8];
            colormap('hsv');
            title('eccentricity model');
        case 6
            plotOver = overlays(3,1).residuals; %ecc residual
            select = overlays(2,1).restrict;
            scalVals = [0 40];
            colormap('Summer');
            title('eccentricity residual');
    end
    
    %     subplot(subRange(1),subRange(2),t),
    hold on
    scatter(fsData.patch.W(:,1),fsData.patch.W(:,2),60,[.6 .6 .6]); % gray color all voxels inbound
    scatter(fsData.patch.W(select,1),fsData.patch.W(select,2),60,plotOver(select),'filled'); % color plot the inbound+restrict voxels
    colorbar;
    caxis([scalVals(1) scalVals(2)]);
    
    %idx is an index based on coords that are inbound!
    wIB = [fsData.patch.W(overlays(2,1).inbounds,1) fsData.patch.W(overlays(2,1).inbounds,2)];
    wIBidx = [wIB(idx,1) wIB(idx,2)];
    wIBidxRef = [wIB(idxRef,1) wIB(idxRef,2)];
    scatter(wIBidx(:,1),wIBidx(:,2),40,[1 0 0],'filled');
    scatter(wIBidxRef(:,1),wIBidxRef(:,2),40,[0 0 0],'filled');
    colourPal = [1 0 0; 0 0 0];
    %draw nice lines 'Benson-style'
    
    for j=1:2, %ref points and data points
        
        if j==2,
            wIBidx = wIBidxRef; %take reference points now
        end
        for i = 1:8,
            switch i
                case 1
                    numID = [1 5 9 13];
                case 2
                    numID = [2 6 10 14];
                case 3
                    numID = [3 7 11 15];
                case 4
                    numID =[4 8 12 16];
                case 5
                    numID =[1 2 3 4];
                case 6
                    numID =[5 6 7 8];
                case 7
                    numID =[9 10 11 12];
                case 8
                    numID =[13 14 15 16];
            end
            p = polyfit(wIBidx(numID,1),wIBidx(numID,2),3); %increase polynomial value if needed
            f = polyval(p,wIBidx(numID,1));%(p,pafit.sortedValues);
            plot(wIBidx(numID,1),f,'Color',colourPal(j,:),'LineWidth',2);
        end
      
    end

 print(gcf, '-dpdf', '-r150', ['Sub_' num2str(sub) 'Hemi_' hemi '_Figure_' num2str(t) '.pdf']) %print figure
    
end
figure,
pol=polar(d2r(subData(:,3)),subData(:,2),'.k');
set(gca,'FontSize',40);
set(pol,'MarkerSize',25);
hold on
pol=polar(d2r(subData(:,5)),subData(:,4),'.r');
set(gca,'FontSize',40);
set(pol,'MarkerSize',25);

print(gcf, '-dpdf', '-r150', ['Sub_' num2str(sub) 'Hemi_' hemi  '_PolarPlot.pdf']) %print figure

end
function [VisCort,idx] = findClosestVoxel(subData,overlays,a,b)

xMap = subData(:,a); %ecc map
yMap = d2r(subData(:,b))-pi; %pa map in degrees
xCort = overlays(3,1).model(overlays(2,1).inbounds);%
yCort = overlays(2,1).model(overlays(3,1).inbounds);%

mat = nan(length(xMap),length(xCort));
for i = 1:length(xMap),
    for j = 1:length(xCort),
        mat(i,j) = sqrt( abs((xCort(j) - xMap(i)).^2) + abs((yCort(j) - yMap(i)).^2) );
        
    end
end

minMat = min(mat,[],2);
for ii=1:length(minMat),
    idx(ii) = find(mat(ii,:)==minMat(ii));
end
%now we know the model values, convert them from ellipticla to cardinal coordinates!
% [xCart,yCart] = elliptical2cart(xCort(idx),yCort(idx));

yMap = yMap + pi;


figure, scatter(xCort,yCort,30,[0 0 0]);
hold on
scatter(xCort(idx),yCort(idx),30,[1 0 0],'filled');
% now we know which voxels, find the right coordinates
VisCort = [xCort(idx) yCort(idx)];

end