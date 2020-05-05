function [] = FOLD_Distorter(filename,a,b,m,n,deviation,rundocfolder)
%FOLD Distorter Function distorts a .FOLD with all vertices deviated from
%their nominal positions according to an input distribution.

%% Reading values from existing .FOLD file
[~,frame,vertices,~,faces,~,~] = FOLD_reader(filename,rundocfolder);

%

%% Deviating Coordinate Points
dev = (deviation/100)*(1/2)*(a+b);

VERTEX = cell2mat(vertices.coords);
% Vertices = zeros(length(vertices.coords),3);
% 
% for ii = 1:length(vertices.coords)/3
%     Vertices(ii,1)v = VERTEX((ii-1)*3+1);
%     Vertices(ii,2) = VERTEX((ii-1)*3+2);
%     Vertices(ii,3) = VERTEX((ii-1)*3+3);
% end

for ii = 1:length(VERTEX)
    VERTEX(ii,1) = VERTEX(ii,1) + (dev*randn);
    VERTEX(ii,2) = VERTEX(ii,2) + (dev*randn);
    VERTEX(ii,3) = VERTEX(ii,3) + (dev*randn);
end

vertices.coords = VERTEX;
%

%% Converting values into correct forms for writing
faces.vertices = cell2mat(faces.vertices);
file.spec = '"1.1"';
file.creator = '"Matlab Code"';
file.author = '"Steve Grey"';
file.classes = '"singleModel"';
frame.title = '"Miura-or Sheet"';
frame.classes = '"foldedForm"';
frame.atrributes = '"3D"';
frame.unit = '"mm"';
%

% %% Plot to check
% h = figure;
% for ii = 1:size(faces.vertices,1)
%     x = vertices.coords(faces.vertices(ii,:),1);
%     y = vertices.coords(faces.vertices(ii,:),2);
%     z = vertices.coords(faces.vertices(ii,:),3);
%     patch(x,y,z,'c')
%     hold on
%     alpha(1)
% end
% 
% for ii = 1:length(vertices.coords)
%     x = vertices.coords(ii,1);
%     y = vertices.coords(ii,2);
%     z = vertices.coords(ii,3);
%     hold on
%     scatter3(x,y,z,50,'m','filled')
% end
% % Make it look presentable
% view(0,90)
% axis equal
% ax = gca;
% ax.XAxis.Visible = 'off';
% ax.YAxis.Visible = 'off';
% ax.ZAxis.Visible = 'off';
% h.Color = [1 1 1];
% set(gcf, 'Position', [680, 49, 1037, 948])
% h.PaperSize = [25 25];

%% Writing to the file

% Open the file and write { to start
fid = fopen(fullfile(rundocfolder,filename),'w');
fprintf(fid,'{\n');

% Write the contents of the 'file' structure
file_names = fieldnames(file);
for ii = 1:length(file_names)
    command_str = sprintf('var_contents = num2str(file.%s);',file_names{ii});
    eval(command_str)
    fprintf(fid,'  "file_%s": %s,\n',file_names{ii},var_contents);
end

% Write the contents of the 'frame' structure
frame_names = fieldnames(frame);
for ii = 1:length(frame_names)
    command_str = sprintf('var_contents = num2str(frame.%s);',frame_names{ii});
    eval(command_str)
    fprintf(fid,'  "frame_%s": %s,\n',frame_names{ii},var_contents);
end

% Write the vertex locations
fprintf(fid,'  "vertices_coords": [\n');
for ii = 1:size(vertices.coords,1)-1
    fprintf(fid,'    [%g,%g,%g],\n',vertices.coords(ii,1),vertices.coords(ii,2),vertices.coords(ii,3));
end
fprintf(fid,'    [%g,%g,%g]\n',vertices.coords(ii+1,1),vertices.coords(ii+1,2),vertices.coords(ii+1,3));
fprintf(fid,'  ],\n');

% Write the faces
fprintf(fid,'  "faces_vertices": [\n');
for ii = 1:size(faces.vertices,1)-1
    fprintf(fid,'    [%g,%g,%g,%g],\n',faces.vertices(ii,1),faces.vertices(ii,2),faces.vertices(ii,3),faces.vertices(ii,4));
end
fprintf(fid,'    [%g,%g,%g,%g]\n',faces.vertices(ii+1,1),faces.vertices(ii+1,2),faces.vertices(ii+1,3),faces.vertices(ii+1,4));
fprintf(fid,'  ],\n');

fprintf(fid,'}');
[~] = fclose(fid);

end