%% Copyright

% Created by G Y SANDESH REDDY , MAYANK CHAUHAN and SHUBHAM SAMANT
% Part of the research work @ Dr Pardha and Dr Supradeepan Research Group, BITS Hydearbad

function writeAnsysMsh_2D(mshFileName,node,elem,total_num_faces,faces,boundary_elements,zone_ids,zone_id_names)
node_size = size(node,1);
fid = fopen([mshFileName '.msh'],'w');
%fprintf(fid,'(0 "Gmsh to Ansys converter")\n');
fprintf(fid,'(0 "Created by G Y SANDESH REDDY , MAYANK CHAUHAN and SHUBHAM SAMANT")\n');
fprintf(fid,'(0 "Part of the research work @ Dr Pardha and Dr Supradeepan Research Group, BITS Hydearbad")\n');
%
fprintf(fid,'(2 2)\n');
fprintf(fid,'(0 "Grid dimensions:")\n');
fprintf(fid,'(10 (0 1 %.2x 0 3))\n',node_size);
fprintf(fid,'(12 (0 1 %.2x 0 0))\n',size(elem,1)+size(boundary_elements,1));
fprintf(fid,'(13 (0 1 %.2x 0 0))\n\n\n',total_num_faces);

fprintf(fid,'(10 (1 1 %.2x  1 3)(\n', node_size);

fprintf(fid,'%f %f %f \n',node');
fprintf(fid,'))\n \n');

face_type = 2;
boundary_type_2 = 2;
boundary_type_1 = 5;

% internal faces
face_id = 2;
int_faces = faces(faces(:,end)~=0,:);
num_faces_start = 1;
num_faces_end = size(int_faces,1);
fprintf(fid,'(13 (%.2x %.2x %.2x %d %d)\n( \n',face_id,num_faces_start,num_faces_end,boundary_type_2,face_type);
fprintf(fid,'%.2x %.2x  %.2x %.2x  \n',int_faces');
fprintf(fid,'))\n \n');

% boundary faces
for i=1:length(zone_ids)
    face_id = face_id+1;
    %boundary_face = boundary_elements(boundary_elements(:,1)==zone_ids(i),3:4);
    boundary_face = boundary_elements(boundary_elements(:,1)==zone_ids(i),1:2);
    boundary_face = sort(boundary_face,2);
    [~,index]=ismember(boundary_face,faces(:,1:2),'rows');
     
    boundary_face = faces(index,:);
    num_faces_start = num_faces_end+1;
    num_faces_end = num_faces_end + size(boundary_face,1);
    fprintf(fid,'(13 (%.2x %.2x %.2x %d %d)\n( \n',face_id,num_faces_start,num_faces_end,boundary_type_1,face_type);
    fprintf(fid,'%.2x %.2x  %.2x %.2x  \n',boundary_face');
    fprintf(fid,'))\n \n');
end

% elements
fprintf(fid,'(12 (1 1 %.2x 1 1)(\n',size(elem,1));
NO_OF_ELEMENTS=repmat(size(elem,2),size(elem,1),1); 
fprintf(fid,'%d %d %d %d %d %d %d %d \n ',NO_OF_ELEMENTS);
fprintf(fid,')())\n \n');

% block 39
fprintf(fid,'(39 (1 fluid fluid)()) \n');
face_id = 2;
fprintf(fid,'(39 (%d interior interior)())\n',face_id);

for i=1:length(zone_ids)
    face_id = face_id+1;
    face_name = zone_id_names{i};
    fprintf(fid,'(39 (%d pressure-outlet %s)())\n',face_id,face_name(2:(end-1)));
end

fclose(fid);
%save('test.msh','-ascii')
