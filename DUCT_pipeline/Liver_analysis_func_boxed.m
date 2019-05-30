function [ bifu, trifu, more, node1, link1] = Liver_analysis_func_boxed(skel, voxel_size)
% Liver_analysis_func_boxed - calculate PV/BD system branching parameters
%
% - graph coversion based on: SKEL2GRAPH3D by Philip Kollmannsberger (philipk@gmx.net)
%
%   - Inputs:
%       skel       - trimmed skeleton of whole system
%       voxel_size - size of a voxel [mm]
%
%   - Outputs:
%       bifu       - number of bifurcations
%       trifu      - number of trifurcations
%       more       - number of quadrifurcations and more
%       node1      - structure describing node properties
%       link1      - structure describing link properties
%
% Jakub Salplachta (jakub.salplachta@ceitec.vutbr.cz)
% ---------------------------------------------------

th1=round(0.2/voxel_size);% Minimum distance of two nodes not to be merged

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
end
