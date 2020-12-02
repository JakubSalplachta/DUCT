%------------------------------------------
% Matlab pipeline for analysis of PV and BD systems, please refer to [1]
%   - input: binary masks for each system, binary masks of main branches
%   - output: numerical results of evaluated parameters
%   
%   - pipeline is based on skeletonization, followed by graph conversion
%   and subsequent analysis
%       - skeletozation based on: SKELETON3D by Philip Kollmannsberger (philipk@gmx.net) [2]
%       - graph coversion based on: SKEL2GRAPH3D by Philip Kollmannsberger (philipk@gmx.net) [2]
% 
% Jakub Salplachta (jakub.salplachta@ceitec.vutbr.cz)

% References: 
% [1] ...
% [2] Kerschnitzki, Kollmannsberger et al., 
% "Architecture of the osteocyte network correlates with bone material quality." 
% Journal of Bone and Mineral Research, 28(8):1837-1845, 2013. 
%------------------------------------------

clearvars
close all
clc
%% Input data
path1PV=''; % Path to whole system mask
path2PV=''; % Path to main branch mask
path1BD=''; % Path to whole system mask
path2BD=''; % Path to main branch mask

%% Output folder
output_folder=''; % Path to output folder
mkdir(output_folder); 

%% Other parameters
sample_id = ''; % Sample ID
voxel_size=0; % Voxel size 
vis=0; % True/false - visualization of branching analysis

% Histogram settings
bin_size = 0.015; % Size of a bin in mm
nbins = 200; % Number of bins - Used by all distance histograms

%% 1. Initial analysis of the systems without bounding box
tic0=tic;
bb=[];
% PV system analysis
[maskPV, ~, ~,~, ~,skelPV, bbPV, ~, T1PV, T2PV,size_orig_PV] = Liver_analysis_func(path1PV, path2PV, bb,voxel_size,vis); % Initial PV system analysis
t1=toc(tic0);

% PV analysis results saving
filename=fullfile(output_folder, strcat(sample_id, 'parameters_PV.xlsx'));
writetable(T1PV,filename,'WriteVariableNames',false,'Sheet','System parameters');
writetable(T2PV,filename,'Sheet','Main branch diamters');
disp(['PV whole system analysis time: ' num2str(t1) ' s'])

tic1=tic;
% BD system analysis
[maskBD, ~, ~, ~, ~, skelBD, bbBD ,~, T1BD, T2BD,size_orig_BD] = Liver_analysis_func(path1BD, path2BD,bb, voxel_size,vis); % Initial BD system analysis
t2=toc(tic1);

% BD analysis results saving
filename=fullfile(output_folder, strcat(sample_id, 'parameters_BD.xlsx'));
writetable(T1BD,filename,'WriteVariableNames',false,'Sheet','System parameters');
writetable(T2BD,filename,'Sheet','Main branch diamters');
disp(['BD whole system analysis time: ' num2str(t2) ' s'])

%% Analysis with bounding box

% bb = [max(bbPV(1), bbBD(1)), min(bbPV(2), bbBD(2)), max(bbPV(3), bbBD(3)), min(bbPV(4), bbBD(4)), max(bbPV(5), bbBD(5)), min(bbPV(6), bbBD(6))];
bb=bbBD; % Systems are cut based to BD system 

% Masks and skeletons trimming
pom=zeros(size_orig_PV(1),size_orig_PV(2),size_orig_PV(3),'single');
pom(bbPV(1):bbPV(2), bbPV(3):bbPV(4), bbPV(5):bbPV(6))=maskPV;
maskPV=pom(bb(1):bb(2), bb(3):bb(4), bb(5):bb(6));

pom=zeros(size_orig_BD(1),size_orig_BD(2),size_orig_BD(3),'single');
pom(bbBD(1):bbBD(2), bbBD(3):bbBD(4), bbBD(5):bbBD(6))=maskBD;
maskBD=pom(bb(1):bb(2), bb(3):bb(4), bb(5):bb(6));

pom=zeros(size_orig_PV(1),size_orig_PV(2),size_orig_PV(3),'single');
pom(bbPV(1):bbPV(2), bbPV(3):bbPV(4), bbPV(5):bbPV(6))=skelPV;
skelPV=pom(bb(1):bb(2), bb(3):bb(4), bb(5):bb(6));

pom=zeros(size_orig_BD(1),size_orig_BD(2),size_orig_BD(3),'single');
pom(bbBD(1):bbBD(2), bbBD(3):bbBD(4), bbBD(5):bbBD(6))=skelBD;
skelBD=pom(bb(1):bb(2), bb(3):bb(4), bb(5):bb(6));

clear pom

tic2=tic;
[bifuPV, trifuPV, morePV, node1PV,linkPV] = Liver_analysis_func_boxed(skelPV, voxel_size); % Analysis of bounded PV system
t3=toc(tic2);
disp(['PV analysis time: ' num2str(t3) ' s'])

tic3=tic;
[bifuBD, trifuBD, moreBD, node1BD,linkBD] = Liver_analysis_func_boxed(skelBD,voxel_size); % Analysis of bounded BD system
t4=toc(tic3);
disp(['BD analysis time: ' num2str(t4) ' s'])

%% Branching distances

%Create point clouds from bifu, trifu, more
tic4=tic;
nTotPV=length(node1PV);
nBifuPV=sum(bifuPV);
nTrifuPV=sum(trifuPV);
nMorePV=sum(morePV);
nPointsPV=nBifuPV+nTrifuPV+nMorePV;
pointsPV=zeros(nPointsPV,3);
j=1;
for i=1:nTotPV
    if bifuPV(i)==1 || trifuPV(i)==1 || morePV(i)==1
        pointsPV(j,:)=[node1PV(i).comx, node1PV(i).comy, node1PV(i).comz] * voxel_size;
        j = j+1;
    end
end 
ptCloudPV = pointCloud(pointsPV);
t5=toc(tic4);
disp(['PV create point cloud time: ' num2str(t5) ' s'])

tic5=tic;
nTotBD=length(node1BD);
nBifuBD=sum(bifuBD);
nTrifuBD=sum(trifuBD);
nMoreBD=sum(moreBD);
nPointsBD=nBifuBD+nTrifuBD+nMoreBD;
pointsBD=zeros(nPointsBD,3);
j=1;
for i=1:nTotBD
    if bifuBD(i)==1 || trifuBD(i)==1 || moreBD(i)==1
        pointsBD(j,:)=[node1BD(i).comx, node1BD(i).comy, node1BD(i).comz] * voxel_size;
        j = j+1;
    end
end
ptCloudBD = pointCloud(pointsBD);
t6=toc(tic5);
disp(['BD create point cloud time: ' num2str(t6) ' s'])

%Loop over all branching points in BD - calculation of distances to neasrest PV
%branching points
tic6=tic;
pointDistancesBD2PV=zeros(nPointsBD, 1);
for j=1:nPointsBD
    [indices,dists] = findNearestNeighbors(ptCloudPV,pointsBD(j,:),1);
    pointDistancesBD2PV(j)=dists(1);
end

avePointDistanceBD2PV=mean(pointDistancesBD2PV);
stdPointDistanceBD2PV=std(pointDistancesBD2PV);
t7=toc(tic6);
disp(['Point distance BD->PV mean: ' num2str(avePointDistanceBD2PV) ' std: ' num2str(stdPointDistanceBD2PV)])
disp(['BD distances to PV cloud time: ' num2str(t7) ' s'])
figure
% maxd = max(pointDistancesBD2PV);
% nbins = 1 + floor(maxd/bin_size);
lastedge = nbins * bin_size;
edges = 0:bin_size:lastedge;
h=histogram(pointDistancesBD2PV, edges);
title('Histogram of Point distances BD->PV');
xlabel('Distance (mm)');
ylabel('Count');
hm = zeros(h.NumBins, 2);
hm(:,1) = h.BinEdges(2:h.NumBins+1);
hm(:,2) = h.Values;
filepath = fullfile(output_folder, strcat(sample_id, 'histPointDistancesBD2PV.csv'));
FID = fopen(filepath, 'w');
fprintf(FID, 'Distance (mm), Count\n'); 
fclose(FID);
dlmwrite(filepath,hm,'delimiter',',','-append');
saveas(h, fullfile(output_folder, strcat(sample_id, 'pointDistancesBD2PV.png')));
filepath = fullfile(output_folder, strcat(sample_id, 'pointDistancesBD2PV.csv'));
FID = fopen(filepath, 'w');
fprintf(FID, 'Distance (mm)\n'); 
fclose(FID);
dlmwrite(filepath, pointDistancesBD2PV,'delimiter',',','-append');
%xlswrite(fullfile(output_folder, strcat(sample_id, 'pointDistancesBD2PV.xlsx'), pointDistancesBD2PV);


%% Surface distances analysis - from BD to PV

% Create point clouds of BD and PV skeleton points
tic7=tic;
[sizer,sizec,sizev]=size(maskPV);

pointsPV=[];
ind=1;
for i=1:length(linkPV)
    pom=linkPV(i).point;
    N = length(pom);
    for j=1:N
        [x,y,z]=ind2sub([sizer, sizec, sizev],pom(j));
        pointsPV(ind,:)=[x,y,z];
        ind=ind+1;
    end
end

ptCloudPV = pointCloud(pointsPV);

pointsBD=[];
ind=1;
for i=1:length(linkBD)
    pom=linkBD(i).point;
    N = length(pom);
    for j=1:N
        [x,y,z]=ind2sub([sizer, sizec, sizev],pom(j));
        pointsBD(ind,:)=[x,y,z];
        ind=ind+1;
    end
end

ptCloudBD=pointCloud(pointsBD);
nPointsBD=length(pointsBD);
t8=toc(tic7);
disp(['Create point clouds time: ' num2str(t8) ' s'])

%Loop over all skeleton points in BD
tic8=tic;
pointDistancesBD2PV_total=zeros(nPointsBD, 1);
pointDistancesBD2PV_surface=zeros(nPointsBD, 1);
maskTotal=maskPV+maskBD;

parfor j=1:nPointsBD
    [indices,dists] = findNearestNeighbors(ptCloudPV,pointsBD(j,:),1);
    pointDistancesBD2PV_total(j)=dists(1);
    r1=pointsBD(j,:);
    r2=pointsPV(indices,:);
    e1_2_2=r2-r1;

    dist1=max([abs(r1(1)-r2(1)),abs(r1(2)-r2(2)),abs(r1(3)-r1(3))]); % Chessboard distance calculation
    e1_2_2=e1_2_2/dist1;
    
    % Creation of line connection between BD skel point and the nearest PV
    % skel poing
    coordinates=zeros(dist1+1,3);
    coordinates(1,:)=r1;
    lineprofile=zeros(1,dist1+1); 
    lineprofile(1)=maskTotal(r1(1),r1(2),r1(3));
    for i=1:dist1
        coordinates(i+1,:)=round(r1+i*e1_2_2); % Coordinate calculation
        lineprofile(i+1)=maskTotal(coordinates(i,1),coordinates(i,2),coordinates(i,3));

    end
    
    % Line profile analysis
    ind=find(lineprofile==0); 
    if isempty(ind)
    else
        r1=coordinates(ind(1)-1,:);
        r2=coordinates(ind(end),:);                
        dist2=norm(r2-r1); % Euclidean distance from surface to surface
        pointDistancesBD2PV_surface(j)=dist2;
    end    
end

%Distances from BD to PV
surfaceDistancesBD2PV = pointDistancesBD2PV_surface * voxel_size; %Scale by voxel size
aveSurfaceDistanceBD2PV=mean(surfaceDistancesBD2PV);
stdSurfaceDistanceBD2PV=std(surfaceDistancesBD2PV);
t9=toc(tic8);
disp(['Surface distance BD->PV mean: ' num2str(aveSurfaceDistanceBD2PV) ' std: ' num2str(stdSurfaceDistanceBD2PV)])
disp(['BD distance to PV time: ' num2str(t9) ' s'])
figure
maxd = max(surfaceDistancesBD2PV);
% nbins = 1 + floor(maxd/bin_size);
lastedge = nbins * bin_size;
edges = 0:bin_size:lastedge;
h=histogram(surfaceDistancesBD2PV, edges);
title('Histogram of Surface distances BD->PV');
xlabel('Distance (mm)');
ylabel('Count');
hm = zeros(h.NumBins, 2);
hm(:,1) = h.BinEdges(2:h.NumBins+1);
hm(:,2) = h.Values;
filepath = fullfile(output_folder, strcat(sample_id, 'histSurfaceDistancesBD2PV.csv'));
FID = fopen(filepath, 'w');
fprintf(FID, 'Distance (mm), Count\n'); 
fclose(FID);
dlmwrite(filepath,hm,'delimiter',',','-append');
saveas(h, fullfile(output_folder, strcat(sample_id, 'surfaceDistancesBD2PV.png')));
%csvwrite(fullfile(output_folder, strcat(sample_id, 'surfaceDistancesBD2PV.csv'), surfaceDistancesBD2PV);
% xlswrite(fullfile(output_folder, strcat(sample_id, 'surfaceDistancesBD2PV.xlsx'), surfaceDistancesBD2PV));

filepath = fullfile(output_folder, strcat(sample_id, 'surfaceDistancesBD2PV.csv'));
FID = fopen(filepath, 'w');
fprintf(FID, 'Distance (mm)\n'); 
fclose(FID);
dlmwrite(filepath,surfaceDistancesBD2PV,'delimiter',',','-append');

%%
tend=toc(tic0);
disp(['Total time: ' num2str(tend) ' s'])
