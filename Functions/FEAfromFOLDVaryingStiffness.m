function [] = FEAfromFOLDVaryingStiffness(filename,meshsize,displacement,m,n,stiffnessvary,rundocfolder)
%FEAfromFOLD Function generates an abaqus input file from a .FOLD file,
%with preset loading conditions with fold line stiffnesses varying randomly
%according to an input variance level.


%% Obtain the geometry from the file

% m = str2num(extractBefore((extractAfter(filename,'_m')),'_'));
% n = str2num(extractBefore((extractAfter(filename,'_n')),'_'));

% [file,frame,vertices,edges,faces,faceOrders,edgeOrders] = FOLD_reader(filename)
[~,frame,vertices,~,faces,~,~] = FOLD_reader(filename,rundocfolder);

% Check units and convert to m
if frame.unit == 'mm'
    for ii = 1:length(vertices.coords)
        vertices.coords{ii} = vertices.coords{ii}/1000;
    end
elseif frame.unit == 'm'
else
   error('Unrecognised units in .FOLD file') 
end

%% Define the material properties

E = sqrt(2300e6*6100e6);      % Young's Modulus (N/m2)
nu = sqrt(0.4*0.17); % Poisson's Ratio
rho = 0.297/370e-6;     % kg/m3

t = 0.37e-3;  % (m)

NominalFoldStiffnessPerLength = 0.1; % Nmm/rad/mm

%% Define the mesh seed

% Element Type
eltype = 'S4R'; % Shell elements

% Run through every face seeding the edges - Defines where the mesh goes
for jj = 1:length(faces.vertices)
    % For the first 3 edges
    for ii = 1:3
        % for x,y, and z
        for kk = 1:3
            meshseed(jj,ii,:,kk) = linspace(vertices.coords{faces.vertices{jj}(ii)}(kk),vertices.coords{faces.vertices{jj}(ii+1)}(kk),meshsize+1);
        end
    end
    
    % Edge 4 links back to edge 1
    for kk = 1:3
        meshseed(jj,4,:,kk) = linspace(vertices.coords{faces.vertices{jj}(4)}(kk),vertices.coords{faces.vertices{jj}(1)}(kk),meshsize+1);
    end
    
    % Plot to check
%     for kk = 1:4
%         pos(kk,:) = vertices.coords{faces.vertices{jj}(kk)};
%     end
%     patch(pos(:,1),pos(:,2),pos(:,3),'g')
%     hold on
%     
%     for ii = 1:4
%         x(:) = meshseed(jj,ii,:,1);
%         y(:) = meshseed(jj,ii,:,2);
%         z(:) = meshseed(jj,ii,:,3);
%         scatter3(x,y,z,[],'r')
%     end
    
end
% axis equal

%% Define the mesh

% Total nodes on a face
total_nodes = (meshsize+1)^2;

% Run through and define the positions of all the nodes
for jj = 1:length(faces.vertices)
    % Initialise a counter to point to 
    count = 0;
    
    % Each row consists of an even distribution of nodes from one side to
    % the other, defined by the meshseed
    for ii = 1:meshsize + 1
        % x, y, and z
        for kk = 1:3
            faces.mesh(jj,count + 1:count + meshsize + 1,kk) = linspace(meshseed(jj,1,ii,kk),meshseed(jj,3,meshsize+2-ii,kk),meshsize+1);
        end
        % Increment counter for the next row
        count = count + meshsize + 1;
    end
    % Plot to check
%     clear x y z
%     x(:) = faces.mesh(jj,:,1);
%     y(:) = faces.mesh(jj,:,2);
%     z(:) = faces.mesh(jj,:,3);
%     scatter3(x,y,z,[],'b')  
end

%% Plot to check and define elements

% Take note of where the edge nodes are. This is currently set-up for an
% even number of nodes on each edge and could be made more intelligent
edge_nodes(1,:) = 1:meshsize+1;
edge_nodes(2,:) = meshsize+1:meshsize+1:total_nodes;
edge_nodes(3,:) = total_nodes-meshsize:total_nodes;
edge_nodes(4,:) = 1:meshsize+1:total_nodes-meshsize;

% h = figure;
% %Plot the facets from the edge nodes
% for ii = 1:length(faces.vertices)
%     
%     % Plot the structure to check (run around all the edges nodes of the
%     % facet and plot)
%     xplot = [faces.mesh(ii,edge_nodes(1,:),1) faces.mesh(ii,edge_nodes(2,:),1) fliplr(faces.mesh(ii,edge_nodes(3,:),1)) fliplr(faces.mesh(ii,edge_nodes(4,:),1))];
%     yplot = [faces.mesh(ii,edge_nodes(1,:),2) faces.mesh(ii,edge_nodes(2,:),2) fliplr(faces.mesh(ii,edge_nodes(3,:),2)) fliplr(faces.mesh(ii,edge_nodes(4,:),2))];
%     zplot = [faces.mesh(ii,edge_nodes(1,:),3) faces.mesh(ii,edge_nodes(2,:),3) fliplr(faces.mesh(ii,edge_nodes(3,:),3)) fliplr(faces.mesh(ii,edge_nodes(4,:),3))];
%     patch(xplot,yplot,zplot,'b')
%     hold on
%     % Part numbers in red next to the structure
%     text(mean(xplot),mean(yplot),mean(zplot)+0.01,num2str(ii),'FontSize',14,'Color','red')
% end
% axis equal

% Join the nodes together into the elements
elementcounter = 0;
% create elements by going through the local nodes
for teller_xi = 1:meshsize % along xi-axis
    for teller_zeta = 1:meshsize % along zeta-axis
        elementcounter = elementcounter + 1;
        elFts(elementcounter,:) = [...
            (meshsize+1)*(teller_xi-1)+teller_zeta, ...
            (meshsize+1)*(teller_xi-1)+teller_zeta+1, ...
            (meshsize+1)*teller_xi+teller_zeta+1,...
            (meshsize+1)*teller_xi+teller_zeta...
            ];
    end
end

disp('Mesh defined')

%% Search for neighbouring facets to be able to connect them

% For every facet 
for ii = 1:length(faces.vertices) 
    % Compare to every other facet
    for ii2 = 1:length(faces.vertices)
        % But not itself
        if ii ~= ii2
            % Look each edge of the first facet
            for jj = 1:4
                % Check if it is next to the next facet
                for kk = 1:4
                    
                    % difference between the nodes meeting up
                    test2(:,:) = sort(faces.mesh(ii,edge_nodes(jj,:),:)) - sort(faces.mesh(ii2,edge_nodes(kk,:),:));
                    
                    % If this is numerically 0 then they are co-indicent
                    if norm(abs(test2)) < mean(mean(mean(faces.mesh)))*1e-8  
                        % Locations define which faces (ii) and edge on that
                        % face (jj) are connected to the first element
                        % (ii2) (location) and second element (edge) (kk)
                        join_locs(ii,jj,:) = [ii2 kk]; 
                        
                        % Plot for clarity and checking
%                         for jj2 = 1:6
%                             x_plot(:) = faces.mesh(ii2,edge_nodes(kk,jj2),1);
%                             y_plot(:) = faces.mesh(ii2,edge_nodes(kk,jj2),2);
%                             z_plot(:) = faces.mesh(ii2,edge_nodes(kk,jj2),3);
%                             hold on
%                             scatter3(x_plot,y_plot,z_plot,100,...
%                             'MarkerEdgeColor','k',...
%                             'MarkerFaceColor',[0.75 0 .75])
%                         end
                
                    end
                end
            end 
        end
    end
end

% Prettify
xlim([min(min(faces.mesh(:,:,1))) max(max(faces.mesh(:,:,1)))])
ylim([min(min(faces.mesh(:,:,2))) max(max(faces.mesh(:,:,2)))])
zlim([min(min(faces.mesh(:,:,3))) max(max(faces.mesh(:,:,3)))])

disp('Folds located and location saved')

%% Write header to file
disp('### Writing to file ###')
foldfile = fullfile(rundocfolder,strcat(filename,'.inp'));
fid = fopen(foldfile,'wt');
fprintf(fid,'*HEADING\n');
fprintf(fid,['Source file: ' filename '\n']);
fprintf(fid,'S.I. Units (N, kg, m, s)\n');
fprintf(fid,'**\n');
text_temp = sprintf('** PARTS: %i FACETS\n',length(faces.vertices));
fprintf(fid,text_temp);
clear text_temp
fprintf(fid,'**\n');

disp('Header written')

%% Write the parts to the file

for ii = 1:length(faces.vertices)
    % Name the facet
    fprintf(fid,['*Part, Name=FacetPart' num2str(ii) '\n']);

    % Write the nodes
    fprintf(fid,'*NODE\n');
    for jj = 1:total_nodes
        fprintf(fid,'%u, %.6f, %.6f, %.6f\n',jj,faces.mesh(ii,jj,1),faces.mesh(ii,jj,2),faces.mesh(ii,jj,3));
    end
    
    % Write the elements
    fprintf(fid,'%s%s%s \n','*ELEMENT, TYPE=',eltype,', ELSET=ElPtFt');
    for jj = 1:size(elFts,1)
        fprintf(fid,'%u, %u, %u, %u, %u\n',jj,elFts(jj,1),elFts(jj,2),elFts(jj,3),elFts(jj,4));
    end
    
    % Define as a shell and define the material
    fprintf(fid,'*SHELL SECTION, ELSET=ElPtFt, MATERIAL=CARD\n');
    % Define the thickness
    fprintf(fid,'%.3G\n',t);
    % Finish
    fprintf(fid,'*END PART\n');
    
end

fprintf(fid,'\n');
disp('Parts written')

%% Write Assembly
fprintf(fid,'*ASSEMBLY, Name=Assembly\n');
fprintf(fid,'\n');

% Write all the faces together into the assembly
for ii = 1:length(faces.vertices)
   
    % Write the facet name
    fprintf(fid,['*Instance, Name=Facet' num2str(ii) ', PART=FacetPart' num2str(ii) ' \n']);
    
    % Need to include the offset which is 0 because the position in space
    % is defined in the part
    fprintf(fid,'0,0,0\n');
    
    % Finish
    fprintf(fid,'*End Instance\n');
    
    % For readability of the input file break every 4 facets (works best
    % for unit cells which have 4 facets so feel free to change)
    if rem(ii,4) == 0
        fprintf(fid,'\n');
    end
end

disp('Assembly written')

%% Write connectors

join_counter = 0;
spring_count = 0;
% For every facet
for ii = 1:length(faces.vertices)
   
    % All the edges of that facet
    for kk = 1:4
        % only look at connecting to facets greater than the current facet then
        % there are no repeats
        if (join_locs(ii,kk,1) > ii)
            
            % Title the join
            fprintf(fid,['*Element, type=CONN3D2, elset=join-' num2str(ii) '-' num2str(join_locs(ii,kk,1)) '\n']);
            
            for jj = 2:size(edge_nodes,2)-1
               
                join_counter = join_counter + 1;
                fprintf(fid,[num2str(join_counter) ', Facet' num2str(ii) '.' num2str(edge_nodes(kk,jj)) ', Facet' num2str(join_locs(ii,kk,1)) '.' num2str(edge_nodes(join_locs(ii,kk,2),jj)) '\n']);
                
%                 % Plot the hinge to check everything has been written
%                 % properly
%                 x_plot = faces.mesh(ii,edge_nodes(kk,jj),1);
%                 y_plot = faces.mesh(ii,edge_nodes(kk,jj),2);
%                 z_plot = faces.mesh(ii,edge_nodes(kk,jj),3);
%                 hold on
%                 scatter3(x_plot,y_plot,z_plot,...
%                     'MarkerEdgeColor','k',...
%                     'MarkerFaceColor',[0 .75 .75])
            end
            
            spring_count = spring_count + 1;
            % Conclude the join
            fprintf(fid,['*Connector Section, elset=join-' num2str(ii) '-' num2str(join_locs(ii,kk,1)) ', behavior=TorsionSpring' num2str(spring_count) '\n']);
            fprintf(fid,'hinge\n');
            fprintf(fid,[num2str(ii) '-' num2str(join_locs(ii,kk,1)) '-CoOrd\n\n']);            
        end
    end
end

disp('Connectors written')

%% Write Node sets

%Write edge node sets as required for current load case
fprintf(fid,'\n**\n');
fprintf(fid,'** Node Sets\n');
fprintf(fid,'**\n\n');

%Left side facets
LFacets = [1:2:(4*m-1)];

%Front side facets
for ii = 0:n-1
    newadd = [1+((n+m)*2*ii),2+((n+m)*2*ii)];
    if ii==0
        FFacets = newadd;
    else
        FFacets = [FFacets,newadd];
    end
end

%Right side facets
for ii = 0:m-1
    newadd = [2+(4*m*(n-1))+(4*ii),4+(4*m*(n-1))+(4*ii)];
    if ii==0
        RFacets = newadd;
    else
        RFacets = [RFacets,newadd];
    end
end

for ii = 1:length(faces.vertices)
    
    instancename = strcat('Facet',num2str(ii));
    
    %Write Left Face Nodes where facet number is correct
    if ismember(ii,LFacets)
        kk = 0;
        fprintf(fid,strcat('*NSET, NSET=LeftFace\n'));
        for jj = 1:meshsize+1:total_nodes-2*meshsize-1
            fprintf(fid, '%s.%u, ',instancename, jj);
            kk = kk+1;
            if kk == 16
                fprintf(fid,'\n');
                kk = 0;
            end
        end
        jj = (total_nodes)-meshsize;
        fprintf(fid, '%s.%u',instancename, jj);
        fprintf(fid,'\n');
    end
    
    %Write relevant Floor Face Nodes for all facets
    faceindex = mod(ii,4);
        fprintf(fid,strcat('*NSET, NSET=FloorFace\n'));
    if faceindex == 1 || faceindex == 2 || faceindex == 0
        kk = 0;
        for jj = 1:meshsize
           fprintf(fid,'%s.%u, ',instancename,jj);
           kk = kk+1;
            if kk == 16
                fprintf(fid,'\n');
                kk = 0;
            end
       end
       jj = meshsize+1;
       fprintf(fid,'%s.%u',instancename,jj);
    end
    if faceindex == 3
        kk = 0;
        for jj = total_nodes-meshsize:total_nodes-1
            fprintf(fid,'%s.%u, ',instancename,jj);
            kk = kk+1;
            if kk == 16
                fprintf(fid,'\n');
                kk = 0;
            end
        end
        jj = total_nodes;
        fprintf(fid,'%s.%u',instancename,jj);
    end
    fprintf(fid,'\n');
    
    %Write Front edge nodes where facet number is correct
    if ismember(ii,FFacets)
        fprintf(fid,strcat('*NSET, NSET=FrontEdge\n'));
        faceindex = mod(ii,2);
        kk = 0;
        if faceindex == 1
            fprintf(fid,'%s.%u\n',instancename,1);
            kk = kk+1;
            if kk == 16
                fprintf(fid,'\n');
                kk = 0;
            end
        end
        kk = 0;
        if faceindex == 0
            fprintf(fid,'%s.%u\n',instancename,meshsize+1);
            kk = kk+1;
            if kk == 16
                fprintf(fid,'\n');
                kk = 0;
            end
        end        
    end
    
    %Write right edge nodes where facet number is correct
    if ismember(ii,RFacets)
        fprintf(fid,strcat('*NSET, NSET=LoadEdge\n'));
        kk = 0;
        for jj = meshsize+1:meshsize+1:total_nodes-meshsize-1
            fprintf(fid,'%s.%u, ',instancename,jj);
            kk = kk+1;
            if kk == 16
                fprintf(fid,'\n');
                kk = 0;
            end
        end
        jj = total_nodes;
        fprintf(fid,'%s.%u',instancename,jj);
        fprintf(fid,'\n');
    end    
end

disp('Node sets written')

%% Write Coordinate systems

% Write header
fprintf(fid,'\n**\n');
fprintf(fid,'** Coordinate Systems\n');
fprintf(fid,'**\n\n');

% Run through every face
for ii = 1:length(faces.vertices)
    % And all 4 edges
    for kk = 1:4
        if (join_locs(ii,kk,1) > ii) && (join_locs(ii,kk,1) ~= 0)
            
            % Define the corners of a panel as a plane from which to create
            % a co-ordinate frame where 1->3 is along the fold
            point1(:) = faces.mesh(ii,edge_nodes(kk,end),:);
            if kk < 3
                point2(:) = faces.mesh(ii,edge_nodes(kk+2,1),:);
            else
                point2(:) = faces.mesh(ii,edge_nodes(kk-2,1),:);
            end
            point3(:) = faces.mesh(ii,edge_nodes(kk,1),:);
            
            % Write co-ordinate system
            fprintf(fid,['*Orientation, name=' num2str(ii) '-' num2str(join_locs(ii,kk,1)) '-CoOrd\n']);
            fprintf(fid,[num2str(point1(1)) ', ' num2str(point1(2)) ', ' num2str(point1(3)) ', ' num2str(point2(1)) ', ' num2str(point2(2)) ', ' num2str(point2(3)) ', ' num2str(point3(1)) ', ' num2str(point3(2)) ', ' num2str(point3(3)) '\n']);
            fprintf(fid,[num2str(1) ', ' num2str(0) '\n\n']); 
        end
    end
end
fprintf(fid,'\n*END ASSEMBLY\n');

disp('Coordinate systems written')

%% Define material properties

% Write header
fprintf(fid,'\n**\n');
fprintf(fid,'** Hinge Spring Stiffness\n');
fprintf(fid,'**\n');

spring_count = 0;
% For every facet
for ii = 1:length(faces.vertices)
   
    % Round all the edges of that facet
    for kk = 1:4
        FoldStiffnessPerLength = NominalFoldStiffnessPerLength + (stiffnessvary/100)*NominalFoldStiffnessPerLength*randn;
        if (join_locs(ii,kk,1) > ii)
            spring_count = spring_count + 1;
            length_vec(:) = faces.mesh(ii,edge_nodes(kk,1),:) - faces.mesh(ii,edge_nodes(kk,end),:);
            FoldStiffness = FoldStiffnessPerLength*(norm(length_vec)/(meshsize-2));
            
            % Write the spring in
            fprintf(fid,'*Connector Behavior, name=TorsionSpring%i\n',spring_count);
            fprintf(fid,'*Connector Elasticity, component=4\n');
            fprintf(fid,'%.6f\n',FoldStiffness);
            
        end
        
    end
end

% Write header
fprintf(fid,'**\n');
fprintf(fid,'** Material Properties\n');
fprintf(fid,'**\n');
fprintf(fid,'*MATERIAL, NAME=CARD\n');

% Write properties
fprintf(fid,'*Density\n');
fprintf(fid,'%.6g, %i\n',rho,30);
fprintf(fid,'*Elastic\n');
fprintf(fid,'%.6g, %.6f\n',E,nu);

disp('Material properties written')

%% Boundary Conditions
fprintf(fid,'*BOUNDARY\n');
fprintf(fid,'FLOORFACE, 3\n');
fprintf(fid,'LEFTFACE, 1\n');
fprintf(fid,'FRONTEDGE, 2\n');

%% Load Step
%
fprintf(fid,'** STEP: LOAD_STEP\n');
fprintf(fid,'** \n');
fprintf(fid,'*Step, name=LOAD_STEP, nlgeom=YES, inc=1000000\n');
fprintf(fid,'*Static\n');
fprintf(fid,'10., 100., 0.0000000000001, 10.\n');
fprintf(fid,'** \n');
fprintf(fid,'** BOUNDARY CONDITIONS\n');
fprintf(fid,'** \n');
fprintf(fid,'** Name: LoadBC Load Type: Displacement/Rotation\n');
fprintf(fid,'*Boundary\n');
%Third number on next line is current right hand edge displacement
fprintf(fid,'LoadEdge, 1, 1, %i\n',displacement);
fprintf(fid,'** \n');
fprintf(fid,'** CONTROLS\n');
fprintf(fid,'** \n');
fprintf(fid,'*Controls, reset\n');
fprintf(fid,'*Controls, parameters=time incrementation\n');
fprintf(fid,', , , , , , , 60, , , \n');
fprintf(fid,'** \n');
fprintf(fid,'** Database Output Information\n');
fprintf(fid,'**\n');
fprintf(fid,'*ENERGY PRINT\n');
fprintf(fid,'*End Step');

%% Close file
[~] = fclose(fid);

disp('### File written ###')
delete (fullfile(rundocfolder,filename));


end