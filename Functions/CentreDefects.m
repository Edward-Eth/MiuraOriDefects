function [filename,width] = CentreDefects(a,b,alphaangle,theta,m,n,deviation,repeat_number,repeats,rundocfolder)
%CentreDefects Function generates a .FOLD with unit cell centre points
%deviated according to an input distribution.


%% Start by defining some of the things about the file we are going to create

% .FOLD Version
file.spec = 1.1;

% Made in MATLAB - Keep the quotes
file.creator = '"MATLAB Code"';

% Author - Keep the quotes
file.author = '"Steve Grey"';

% Class:
% - "singleModel": A single origami model, possibly still in multiple frames to represent crease pattern, folded form, etc
file.classes = '"singleModel"';

frame.title = '"Miura-ori Sheet"';

% This code is geared towards only doing the folded form here
frame.classes = '"foldedForm"';

% Again this code is only aimed towards 3D here
frame.attributes = '"3D"';

% NOTE: all units are in mm throughout
frame.unit = '"mm"';

%% Next define some of the things about the file we are going to create
dev = (deviation/100)*(1/2)*(a+b);


% Define ideal unit cell parameters [Schenk & Guest (2013)] See image in
% folder for reference
H = a*sind(theta)*sind(alphaangle);
S = b*cosd(theta)*tand(alphaangle)/sqrt(1+cosd(theta)^2*tand(alphaangle)^2);
L = a*sqrt(1-sind(theta)^2*sind(alphaangle)^2);
V = b/sqrt(1+cosd(theta)^2*tand(alphaangle)^2);
width = n*2*S;

        unit_cell_vertices(1,:) = [0 0 0];
        unit_cell_vertices(2,:) = [0 L H];
        unit_cell_vertices(3,:) = [0 2*L 0];
        
        unit_cell_vertices(4,:) = [S V 0];
        %Unit cell vertex 5 defined in loop to allow for variation
        unit_cell_vertices(6,:) = [S 2*L+V 0];
        unit_cell_vertices(5,:) = [S+(dev*randn) L+V+(dev*randn) H];       
        unit_cell_vertices(7,:) = [2*S 0 0];
        unit_cell_vertices(8,:) = [2*S L H];
        unit_cell_vertices(9,:) = [2*S 2*L 0];

count = 0;
for jj = 1:m
    for ii = 1:n
        % Define the vertices in the current unit cell (From back left up to the top
        % down to the front and finally to the bottom. Then moving on to
        % the middle and finally the right)

        unit_cell_vertices(5,:) = [S+(dev*randn) L+V+(dev*randn) H];
        
        count = count + 1;
        vertices.coords((count-1)*9+1:(count-1)*9+9,1) = unit_cell_vertices(:,1) + 2*S*(ii-1);
        vertices.coords((count-1)*9+1:(count-1)*9+9,2) = unit_cell_vertices(:,2) + 2*L*(jj-1);
        vertices.coords((count-1)*9+1:(count-1)*9+9,3) = unit_cell_vertices(:,3);

    end
end

vertices.coords = remove_extra_vertices(vertices.coords);
%         vertices.coords = unique(vertices.coords,'rows','stable');

% Define the faces in a unit cell on the top
unit_cell_faces(1,:) = [1 2 5 4];
unit_cell_faces(2,:) = [4 5 8 7];
unit_cell_faces(3,:) = [2 3 6 5];
unit_cell_faces(4,:) = [6 5 8 9];

for ii = 1:n
    for jj = 1:m
        if jj == 1
            node_nums(ii,jj,:) = [1:9]+6*(ii-1); %#ok<NBRAK>
        elseif jj == 2
            node_nums(ii,jj,[2 3]) = node_nums(1,jj-1,[2 3]) + 6*n+2 + (ii-1)*4; %#ok<*SAGROW>
            node_nums(ii,jj,[5 6]) = node_nums(1,jj-1,[5 6]) + 6*n+1 + (ii-1)*4;
            node_nums(ii,jj,[8 9]) = node_nums(1,jj-1,[8 9]) + 6*n+0 + (ii-1)*4;
            node_nums(ii,jj,[1 4 7]) = node_nums(ii,jj-1,[3 6 9]);
        else
            node_nums(ii,jj,:) = node_nums(ii,jj-1,:) + 4*n+2;
            node_nums(ii,jj,[1 4 7]) = node_nums(ii,jj-1,[3 6 9]);
            
        end
    end
end

% Use the unit cell faces to define the vertices which surround
% every face in the tube
count = 0;
for ii = 1:n
    for jj = 1:m
        for kk = 1:4
            faces.vertices(4*count + kk,:) = node_nums(ii,jj,unit_cell_faces(kk,:));
        end
        count = count + 1;
    end
end

%% Plot to check
h = figure;
for ii = 1:size(faces.vertices,1)
    x = vertices.coords(faces.vertices(ii,:),1);
    y = vertices.coords(faces.vertices(ii,:),2);
    z = vertices.coords(faces.vertices(ii,:),3);
    patch(x,y,z,'c')
    hold on
    alpha(1)
end

for ii = 1:length(vertices.coords)
    x = vertices.coords(ii,1);
    y = vertices.coords(ii,2);
    z = vertices.coords(ii,3);
    hold on
    scatter3(x,y,z,50,'m','filled')
end
% Make it look presentable
view(0,90)
axis equal
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ax.ZAxis.Visible = 'off';
h.Color = [1 1 1];
set(gcf, 'Position', [680, 49, 1037, 948])
h.PaperSize = [25 25];

%% Writing to the file

% Define a descriptive filename
if repeats == 1
    filename = sprintf('n%i_m%i_a%i_b%i_centredevpercent%g.FOLD',n,m,a,b,deviation);
else
    filename = sprintf('n%i_m%i_a%i_b%i_centredevpercent%g_repeat%i.FOLD',n,m,a,b,deviation,repeat_number);
end

% Open the file and write { to start
fid = fopen((fullfile(rundocfolder,filename)),'w');
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