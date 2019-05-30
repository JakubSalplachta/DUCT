function [mask, bifu, trifu, more, node1, skel, bb, link1, T1,T2,size_orig] = Liver_analysis_func(path1, path2,bb, voxel_size,vis)
% Liver_analysis_func - calculate PV/BD system morphological ana branching parameters
% - skeletozation based on: SKELETON3D by Philip Kollmannsberger (philipk@gmx.net)
% - graph coversion based on: SKEL2GRAPH3D by Philip Kollmannsberger (philipk@gmx.net)

%   - Inputs:
%       path1      - path to whole system binary mask
%       path2      - path to main branch binary mask
%       bb         - initial bounding box, if specified
%       voxel_size - size of a voxel [mm]
%       vis        - True/false - visualization of branching analysis
%
%   - Outputs:
%       mask       - loaded binary mask of whole system
%       bifu       - number of bifurcations
%       trifu      - number of trifurcations
%       more       - number of quadrifurcations and more
%       node1      - structure describing node properties
%       skel       - whole system skeleton
%       bb         - calculated bounding box coordinates
%       link1      - structure describing link properties
%       T1         - table containing results of the analysis
%       T2         - table containing results of diameter analysis of
%                    main branch
%       size_orig  - size of input binary mask
%
% Jakub Salplachta (jakub.salplachta@ceitec.vutbr.cz)
%-------------------------------
extension=300; % Number of extension voxels in each dimmension - needed for diamater calculation
th1=round(0.2/voxel_size);% Minimum distance of two nodes not to be merged
th2=round(1.5/voxel_size);% Distance measure for diameter calculation 
%% Whole system analysis

% 1. Data loading
files=dir([path1 '\' '*.tif']);
if isempty(bb)
    files=files;
else
    files=files(bb(5):bb(6),:);
end
pom=single(imread([path1 '\' files(1).name]));

mask=zeros(size(pom,1),size(pom,2),length(files),'single');
parfor i=1:length(files)
    mask(:,:,i)=single(imread([path1 '\' files(i).name]));
end
mask=logical(mask);
size_orig=[size(mask,1),size(mask,2),size(mask,3)];

if isempty(bb)
    bb = bb3(mask);
else
    bb(5)=1;
    bb(6)=size(mask,3);
end

mask = mask(bb(1):bb(2), bb(3):bb(4), bb(5):bb(6));
[a,b,c]=size(mask);

% 2. Skeletonization and Graph conversion
skel=Skeleton3D(mask);% Skeletonization
[A,node1,link1] = Skel2Graph3D(skel,th1);% Graph conversion 

% First graph analysis - false and true edpoints detection, false nodes detection
endpoints=[];
false_nodes=[];
false_endpoints=[];
for i=1:length(node1)
    if node1(i).ep==1
        if isempty(node1(i).conn)
            false_endpoints=[false_endpoints i];
        else
            endpoints=[endpoints i];
                   
        end
    else 
        if length(unique(node1(i).conn))<3 && length(node1(i).links)<3
            false_nodes=[false_nodes i];
        end
    end
end

% Detection of nodes to be merged - false endpoints, true endpoints and false nodes are omitted
A_triangular=triu(full(A));
A_triangular(A_triangular>th1)=0;
A_triangular=logical(A_triangular);
Node_matrix=ones(length(node1));
Node_matrix([endpoints false_endpoints false_nodes],:)=0;
Node_matrix(:,[endpoints false_endpoints false_nodes])=0;
merging=A_triangular.*Node_matrix; % Nodes to be merged

% Auxilary variables
pom=sum(merging,2)+1;
pom([endpoints false_endpoints false_nodes])=0;
pom3=ones(length(node1),1);
pom3([endpoints false_endpoints false_nodes])=0;

% Variables with merging and alysis results
merged=zeros(1,length(node1));
nodes_merging=[];
bifu=zeros(1,length(node1)); % Bifurcations
trifu=zeros(1,length(node1)); % Trifurcations
more=zeros(1,length(node1)); % Quadrifurcation and more
NA=zeros(1,length(node1)); % Not classified

while sum(pom)>0
    [M, I]=max(pom);
    for i=1:length(I)
            nodes=find(merging(I(i),:)>0);
            pom2=node1(I(i)).conn;
            conns=pom2;
            for j=1:length(nodes)
                pom2=node1(nodes(j)).conn;
                conns=[conns, pom2];
            end
            pos=find(conns==I(i));
            conns(pos)=[];          
            for j=1:length(nodes)
                pos=find(conns==nodes(j));
                conns(pos)=[]; 
            end
            [counts]=histcounts(conns,'BinMethod','integers');
            pom2=find(counts==1);
            val=length(pom2)-1;
        if val==2
            bifu(I(i))=1;
        elseif val==3
            trifu(I(i))=1;
        elseif val>3
            more(I(i))=1;
        else
            NA(I(i))=1;
        end
        pom3([I(i) nodes])=0;
         if M>1
             pom4=[I(i) nodes];
             for k=1:length(pom4)
                 nodes_merging(pom4(k)).merging=1;
                 nodes_merging(pom4(k)).nodes=sort([I(i) nodes]);
             end
             merged(pom4)=1;
         end
    end  
    merging([I nodes],:)=0;
    merging(:,[I nodes])=0;
    pom=(sum(merging,2)+1).*pom3;
    val=0;  
end

% 4. Length analysis

% Total length calculation
w_total_real=0;
for i=1:length(link1)
    pom=link1(i).point;
    for j=1:length(pom)-1
        [x1,y1,z1]=ind2sub([a,b,c],pom(j));
        [x2,y2,z2]=ind2sub([a,b,c],pom(j+1));
        distance=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2); % 3D Euclidean distance
        w_total_real=w_total_real+distance;
    end
end
% w_total_real= sum(cellfun('length',{link1.point}));% Total length of network

% Theoretical length of ideal system
w_total_ideal=0;
for i=1:length(link1)
    pom=link1(i).point;
    [x1,y1,z1]=ind2sub([a,b,c],pom(1));
    [x2,y2,z2]=ind2sub([a,b,c],pom(end));
    distance=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2); % 3D Euclidean distance
    w_total_ideal=w_total_ideal+distance;
end

%% Visualization of whole system branching analysis results
if vis==1
    x=[];
    xx=[];
    y=[];
    yy=[];
    z=[];
    zz=[];
    figure;
    hold on;
    pos=1:length(node1);
    pos(false_endpoints)=[];
    false_nodes_pom=zeros(1,length(node1));
    false_nodes_pom(false_nodes)=1;
    for i=pos
    x1 = node1(i).comx;
    y1 = node1(i).comy;
    z1 = node1(i).comz;

    for j=1:length(node1(i).links)    % Draw all connections of each node          
        for k=1:length(link1(node1(i).links(j)).point)-1   % Draw edges as lines using voxel positions          
            [x3,y3,z3]=ind2sub([a,b,c],link1(node1(i).links(j)).point(k));
            [x2,y2,z2]=ind2sub([a,b,c],link1(node1(i).links(j)).point(k+1));
            x(end+1)=x2;
            y(end+1)=y2;
            z(end+1)=z2;
            xx(end+1)=x3;
            yy(end+1)=y3;
            zz(end+1)=z3;
           line([y3 y2],[x3 x2],[z3 z2],'Color','k','LineWidth',2);
        end
    end

    % Draw all nodes - color selected based on node type 
    if(node1(i).ep==1)
       ncol = 'r';
           plot3(y1,x1,z1,'o','Markersize',9,...
        'MarkerFaceColor',ncol,...
        'Color','k');
    else
       if bifu(i)==1
           ncol='y';
           plot3(y1,x1,z1,'o','Markersize',9,...
          'MarkerFaceColor',ncol,...
          'Color','k');
       elseif trifu(i)==1
           ncol='blue';
           plot3(y1,x1,z1,'o','Markersize',9,...
          'MarkerFaceColor',ncol,...
          'Color','k');
       elseif more(i)==1
           ncol='g';
           plot3(y1,x1,z1,'o','Markersize',9,...
          'MarkerFaceColor',ncol,...
          'Color','k');
       end
    end
    end
    axis image;axis off;
    set(gcf,'Color','white');
    drawnow;
end
%% Main Branch analysis

% 1. Data loading
files=dir([path2 '\' '*.tif']);
files=files(bb(5):bb(6),:);

% pom=single(imread([path2 '\' files(1).name]));

mask2_orig=zeros(size(mask),'single');
parfor i=1:length(files) 
    pom=imread([path2 '\'  files(i).name]);
    mask2_orig(:,:,i)=pom(bb(1):bb(2), bb(3):bb(4));
end
% mask2_orig=logical(mask2_orig(bb(1):bb(2), bb(3):bb(4),:));

mask2_orig=logical(mask2_orig);
mask2_orig=permute(mask2_orig,[3 2 1]); % Dimensions permutation - for optimal diameter calculation
[a,b,c]=size(mask2_orig);

% 2. Skeletonization and Graph conversion

skel2=Skeleton3D(mask2_orig); % Skeletonization

% Extension - needed for diameter analysis
pom=zeros(a+extension,b+extension,c+extension);
pom(extension/2+1:end-extension/2,extension/2+1:end-extension/2,extension/2+1:end-extension/2)=mask2_orig;
mask2_orig=logical(pom);

pom=zeros(a+extension,b+extension,c+extension);
pom(extension/2+1:end-extension/2,extension/2+1:end-extension/2,extension/2+1:end-extension/2)=skel2;
skel2=logical(pom);

[a,b,c]=size(mask2_orig);

[A,node2,link2] = Skel2Graph3D(skel2,th1);% Graph conversion 


% 3. Length analysis

% Total length calculation
w_branch_real=0;
for i=1:length(link2)
    pom=link2(i).point;
    for j=1:length(pom)-1
        [x1,y1,z1]=ind2sub([a,b,c],pom(j));
        [x2,y2,z2]=ind2sub([a,b,c],pom(j+1));
        distance=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2); % 3D Euclidean distance
        w_branch_real=w_branch_real+distance;
    end
end

% Theoretical length of ideal system
w_branch_ideal=0;
for i=1:length(link2)
    pom=link2(i).point;
    [x1,y1,z1]=ind2sub([a,b,c],pom(1));
    [x2,y2,z2]=ind2sub([a,b,c],pom(end));
    distance=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2); % 3D Euclidean distance
    w_branch_ideal=w_branch_ideal+distance;
end

% 4. Diameter calculation

% Mask boundary detection
% SE=strel('sphere',1);
load('SE.mat')
pompom=imerode(mask2_orig,SE);
boundary=mask2_orig-pompom;


points=[];
for i=1:length(link2)
    points=[points link2(i).point];
end

points=sort(unique(points));
delta=2;
positions=1+delta:th2:length(points)-delta;
area_size=50; % Number of voxels of search area in each direction and dimmension
ind=-delta:delta;
pom=[];
diameters=zeros(length(positions),1);
for i=1:length(positions)
    for j=1:length(ind)
        [x,y,z]=ind2sub([a,b,c],points(positions(i)+ind(j)));
        boundary_pom=boundary(x-area_size:x+area_size,y-area_size:y+area_size,z-area_size:z+area_size);
        centroid=zeros(size(boundary_pom));
        centroid(area_size+1,area_size+1,area_size+1)=1;
        edt=bwdist(centroid);
        pom(j)=min(edt(boundary_pom==1));
    end
    diameters(i)=mean(pom);
end

ratio=[];
for i=1:length(positions)-1
    ratio(i)=diameters(i+1)/diameters(i);
end
%% Visualization of main branch skeletonization and graph coversion result
if vis
    x=[];
    xx=[];
    y=[];
    yy=[];
    z=[];
    zz=[];
    figure;
    hold on;
    pos=1:length(node2);
    for i=pos
       x1 = node2(i).comx;
       y1 = node2(i).comy;
       z1 = node2(i).comz;

        for j=1:length(node2(i).links)    % Draw all connections of each node
            for k=1:length(link2(node2(i).links(j)).point)-1 % Draw edges as lines using voxel positions           
                [x3,y3,z3]=ind2sub([a,b,c],link2(node2(i).links(j)).point(k));
                [x2,y2,z2]=ind2sub([a,b,c],link2(node2(i).links(j)).point(k+1));
                x(end+1)=x2;
                y(end+1)=y2;
                z(end+1)=z2;
                xx(end+1)=x3;
                yy(end+1)=y3;
                zz(end+1)=z3;
               line([y3 y2],[x3 x2],[z3 z2],'Color','k','LineWidth',2);
            end
        end

        % Draw all nodes - color selected based on node type 
       if(node2(i).ep==1)
           ncol = 'r';
               plot3(y1,x1,z1,'o','Markersize',9,...
            'MarkerFaceColor',ncol,...
            'Color','k');
       end
    end
    axis image;axis off;
    set(gcf,'Color','white');
    drawnow;
end

%% Summary
system_length_real=w_total_real*voxel_size;
system_length_ideal=w_total_ideal*voxel_size;
main_branch_length_real=w_branch_real*voxel_size;
main_branch_length_ideal=w_branch_ideal*voxel_size;
system_volume=sum(mask(:))*voxel_size^3;
bifurcations=sum(bifu);
trifurcations=sum(trifu);
more_nodes=sum(more);
average_diameters=diameters.*voxel_size*2;

disp(' ')
disp(['Total length of system: ' num2str(system_length_real) ' mm'])
disp(['Theoretical length of system: ' num2str(system_length_ideal) ' mm'])
disp(['Average length from trunk to tip (main branch): ' num2str(main_branch_length_real) ' mm'])
disp(['Theoretical length of the main branch: ' num2str(main_branch_length_ideal) ' mm'])
disp(['System volume: ' num2str(system_volume) ' mm^3'])
disp(['Number of bifurcations: ' num2str(bifurcations)])
disp(['Number of trifurcations: ' num2str(trifurcations)])
disp(['Number of quadri(and more)furcations: ' num2str(more_nodes)])
disp(['Average diameters of the main branch: ' num2str(mean(average_diameters))])


Parameters={'Total length of system [mm]'; 'Theoretical length of system [mm]';'Average length from trunk to tip (main branch) [mm]'; 'Theoretical length of the main branch [mm]';  'System volume [mm^3]'; 'Number of bifurcations'; 'Number of trifurcations';'Number of quadri(and more)furcations'; 'Average diameter of the main branch [mm]'};
Values=[system_length_real;system_length_ideal; main_branch_length_real; main_branch_length_ideal; system_volume;bifurcations;trifurcations; more_nodes; mean(average_diameters)];
T1=table(Parameters, Values);
T2=table(average_diameters);
T2.Properties.VariableNames={'Diameters'};
T2.Properties.VariableUnits={'mm'};

end
