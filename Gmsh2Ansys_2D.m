%% Copyright

% Created by G Y SANDESH REDDY , MAYANK CHAUHAN and SHUBHAM SAMANT
% Part of the research work @ Dr Pardha and Dr Supradeepan Research Group, BITS Hydearbad

function Gmsh2Ansys_2D(gmshFileName,ansysFileName)

% Open gmsh file
fileName = [gmshFileName '.msh'];
fid = fopen(fileName);

GmshMeshFileFormatFlag = 0;
GmshPhysicalNameFlag = 0;
GmshNodesFlag = 0;
domain_dimension = 0;

while ~feof(fid)
    xp = fscanf(fid, '%s', 1);
    
    % Read in Gmsh File Format
    if strcmp(xp,'$MeshFormat')
        GmshMeshFileFormatFlag = 1;
        fscanf(fid,'%f',1);
        fscanf(fid,'%d',1);
        fscanf(fid,'%d',1);    
    end
    
    % Read in the deatils of Physical Names from Gmsh
    if strcmp(xp,'$PhysicalNames')
        if (GmshMeshFileFormatFlag ~= 1)
            error('Invalid Gmsh File. MeshFormat tag not found in the mesh file');
        end
        GmshPhysicalNameFlag = 1;
        no_of_names= fscanf(fid, '%d', 1);
       
        zone_ids = [];
        zone_id_names = [];
        
        for i=1:no_of_names
            current_dimension = fscanf(fid,'%d',1);
            domain_dimension=max(domain_dimension,current_dimension);
            zone_id=fscanf(fid,'%d',1);
            boundary_name= fscanf(fid, '%s', 1);
            
            if (domain_dimension == 3)
                error('Use Gmsh2Ansys_3D instead of Gmsh2Ansys_2D. Three dimensional elements found in the Gmsh file.');
            end
            
            if (current_dimension == 2)
                continue;
            end
           zone_ids = [zone_ids;zone_id];
           zone_id_names = [zone_id_names;{boundary_name}];           
        end
    end
    
    %nodes
    if strcmp(xp,'$Nodes')
        if (GmshPhysicalNameFlag ~= 1)
            error('Invalid Gmsh File. PhysicalName tag not found in the mesh file');
        end
        GmshNodesFlag = 1;
        
        % read in number of nodes
        no_nodes = fscanf(fid, '%d', 1);              
        
        previousNodes_num = 0;
        i=1;
        while (i<=no_nodes)
            text = fgetl(fid);
            if (strcmp(text,'') | strcmp(text,' ') | isempty(text))
                continue;
            end
            text_node = str2num(text);
            Nodes_num = length(text_node);
            if Nodes_num ==0                         %blank line to be skipped
                continue;
            end
            if (previousNodes_num == 0)
                previousNodes_num = Nodes_num;            %when 3rd column is not there i.e only 2 columns
            elseif (Nodes_num>previousNodes_num)                 
                nodes(:,(previousNodes_num+1):Nodes_num) = zeros(size(nodes,1),(Nodes_num-previousNodes_num)); %when 8th column is there where we have added zeros column matrix on 8th column and then iteraion for i is done
                previousNodes_num = Nodes_num;
            end
            nodes(i,1:Nodes_num) = text_node;
            i=i+1;
        end
        
        %elements
        elseif strcmp(xp,'$Elements')
            if (GmshNodesFlag ~= 1)
                error('Invalid Gmsh File. Nodes tag found before Elements tag in the mesh file');
            end
        no_elements = fscanf(fid, '%d', 1);             
        previousNumNodes = 0;
        i=1;
        while (i<=no_elements)
            text = fgetl(fid);
            if (strcmp(text,'') | strcmp(text,' ') | isempty(text))
                continue;
            end
            data = str2num(text);   %reads line by line
            noColumns = length(data);
            if noColumns ==0                         %blank line to be skipped
                continue;
            end
            startElemNodes = 4+data(3);
            noElemNodes = noColumns-startElemNodes+2;
            if (previousNumNodes == 0)
                previousNumNodes = noElemNodes;            %when 8th column is not there i.e only 7 columns
            elseif (noElemNodes>previousNumNodes)                 
                elements(:,(previousNumNodes+1):noElemNodes) = zeros(size(elements,1),(noElemNodes-previousNumNodes)); %when 8th column is there where we have added zeros column matrix on 8th column and then iteraion for i is done
                previousNumNodes = noElemNodes;
            end
            elements(i,1:noElemNodes) = data([4,startElemNodes:noColumns]);        %element number     elem-type      number of tags       <tags>    node number list
            i=i+1;
        end
    end 
end
           
fclose(fid);

node = nodes(:,2:end);%x y cordinates

elemF=elements(elements(:,1)==zone_id,2:end);                                  
nodeF = node(sort(unique(elemF(:))),:);
I_F=elemF(:);                      %all nodes in single column(fluid)
m_F=sort(unique(I_F));                             %sorted nodes
MAPf=zeros(size(nodeF,1),1);                   %mapped fluid nodes
for i=1:size(m_F)
    MAPf(m_F(i,1),1)=i;                  %mapped fluid nodes
end
elemF(:)= MAPf(elemF(:));
size_elemF=size(elemF,1);                  %size of fluid elements
edgeF=[elemF(1:size_elemF,1:2);elemF(1:size_elemF,2:3);elemF(1:size_elemF,3),elemF(1:size_elemF,1)];   %edge matrix reshaped  
sorted_edgeF=sort(edgeF,2);   %%edge matrix sorted i.e ascending order      
                                                
[~,L1_F]=unique(sorted_edgeF,'first','rows');   
[~,L2_F]=unique(sorted_edgeF,'last','rows');
internal_edgeF=(L1_F~=L2_F);%logic created whether equal or not

%shows elements sharing common edge
elements_sharing_edgeF = zeros(size(sorted_edgeF,1),2);

elements_sharing_edgeF(:,1)=1:(size(sorted_edgeF,1));
elements_sharing_edgeF(L1_F(internal_edgeF),:)=[L1_F(internal_edgeF>0) L2_F(internal_edgeF)];
elements_sharing_edgeF(L2_F(internal_edgeF),:)=[L1_F(internal_edgeF>0) L2_F(internal_edgeF)];

elements_sharing_edgeF(elements_sharing_edgeF>(2*size_elemF))=elements_sharing_edgeF(elements_sharing_edgeF>(2*size_elemF))-size_elemF;
elements_sharing_edgeF(elements_sharing_edgeF>(1*size_elemF))=elements_sharing_edgeF(elements_sharing_edgeF>(1*size_elemF))-size_elemF; 

% sorting the elements shared by a edge in the order required by Ansys
elem_1=elements_sharing_edgeF(:,1);
elem_2=elements_sharing_edgeF(:,2);
Faces=([  (sorted_edgeF)  ,  elem_1  ,elem_2]);

Faces=unique(Faces,'rows');

%boundary elements
boundary_elements = elements(elements(:,1)~=zone_id,2:end); 
total_num_Faces=size(Faces,1);
writeAnsysMsh_2D(ansysFileName,nodeF,elemF,total_num_Faces,Faces,boundary_elements,zone_ids,zone_id_names)
