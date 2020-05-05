% A single run generator for distorted Miura Ori Meshes, all necessary
% arguments are defined here, and used to call the various generation
% functions. Creates abaqus.inp files, batch files to run these and matlab
% scripts to extract the desired data. Tear runs require the
% communications toolbox to be installed.

close all
addpath('Functions');

%% Step One: Define facet properties

Side_Length_A = 30; % (mm)
Side_Length_B = 30; % (mm)
Angle_Alpha = 60; % (degrees)
Angle_Theta = 60; % (degrees)

%

%% Step Two: Define unit cell counts

Count_M = 2; % Unit cells in y direction
Count_N = 2; % Unit cells in x direction

%

%% Step Three: Define deviation method

Deviation_Method_Switch = 0;

% 0 = No Deviation
% 1 = Centre Point Deviation (LEGACY)
% 2 = All Points Deviation
% 3 = No Deviation with tears
% 4 = Deviated with tears
% 5 = No Deviation, fold stiffness varying
% 6 = Deviation with varying fold stiffness
% 7 = Tears with varying fold stiffness
% 8 = Deviation, tears and varying fold stiffness

%

%% Step Four: Define deviations

%Deviation = inputdev; % standard deviation of nodes function execution
Deviation = 5; %Standard deviation of nodes manual define
Tears = 0.04; %Proportion of folds to be torn
StiffnessVary = 20; %Standard deviation percentage of stiffness of edges

if Deviation_Method_Switch == 0 || Deviation_Method_Switch == 3 || Deviation_Method_Switch == 5
    Deviation = 0;
end

if Deviation_Method_Switch == 0 || Deviation_Method_Switch == 1 || Deviation_Method_Switch == 2 || Deviation_Method_Switch == 5 || Deviation_Method_Switch == 6
    Tears = 0;
end

if Deviation_Method_Switch == 0 || Deviation_Method_Switch == 1 || Deviation_Method_Switch == 2 || Deviation_Method_Switch == 3 || Deviation_Method_Switch == 4
    StiffnessVary = 0;
end

%

%% Step Five: Define physical, mesh and loading properties

Meshsize = 5; % elements per edge
Displacement = 0.4; % proportion of starting width to compress by
%Add a switch between absolute value displacement and proportional
%displacement?

%

%% Step Six: Repetition options

Repetition = 0; %Logical Switch
Repeat_Count = 4; %Number of Repeats

%% Generation and meshing

if Repetition == 0
    Repeat_Count = 1;
end

%% Naming

if(StiffnessVary ==0)
    if(Deviation ==0)
        if(Tears ==0)
            if(Repeat_Count ==1)
                %No stiffness vary, no deviation, no tears, no repeats
                rundocfolder = sprintf('runfolder_method%i',Deviation_Method_Switch);
                mkdir (rundocfolder);
                rundoc = fullfile(rundocfolder,sprintf('rundoc_method%i.bat',Deviation_Method_Switch));
                extractname = sprintf('extract_method%i.m',Deviation_Method_Switch);
                extractdoc = fullfile(rundocfolder,extractname);
                fid = fopen(rundoc,'w');
                fid1 = fopen(extractdoc,'w');
            else
                %No stiffness vary, no deviation, no tears, repeats
                rundocfolder = sprintf('runfolder_method%i_rep%i',Deviation_Method_Switch,Repeat_Count);
                mkdir (rundocfolder);
                rundoc = fullfile(rundocfolder,sprintf('rundoc_method%i_rep%i.bat',Deviation_Method_Switch,Repeat_Count));
                extractname = sprintf('extract_method%i_rep%i.m',Deviation_Method_Switch,Repeat_Count);
                extractdoc = fullfile(rundocfolder,extractname);
                fid = fopen(rundoc,'w');
                fid1 = fopen(extractdoc,'w');
            end
        else
            if(Repeat_Count==1)
                %No stiffness vary, no deviation, tears, no repeats
                rundocfolder = sprintf('runfolder_method%i_tears%g',Deviation_Method_Switch,Tears);
                mkdir (rundocfolder);
                rundoc = fullfile(rundocfolder,sprintf('rundoc_method%i_tears%g.bat',Deviation_Method_Switch,Tears));
                extractname = sprintf('extract_method%i_tears%g.m',Deviation_Method_Switch,Tears);
                extractdoc = fullfile(rundocfolder,extractname);
                fid = fopen(rundoc,'w');
                fid1 = fopen(extractdoc,'w');
            else
                %No stiffness vary, no deviation, tears, repeats
                rundocfolder = sprintf('runfolder_method%i_tears%g_rep%i',Deviation_Method_Switch,Tears,Repeat_Count);
                mkdir (rundocfolder);
                rundoc = fullfile(rundocfolder,sprintf('rundoc_method%i_tears%g_rep%i.bat',Deviation_Method_Switch,Tears,Repeat_Count));
                extractname = sprintf('extract_method%i_tears%g_rep%i.m',Deviation_Method_Switch,Tears,Repeat_Count);
                extractdoc = fullfile(rundocfolder,extractname);
                fid = fopen(rundoc,'w');
                fid1 = fopen(extractdoc,'w');
            end
        end
    else
        if(Tears==0)
            if(Repeat_Count==1)
                %No stiffness vary, deviation, no tears, no repeats
                rundocfolder = sprintf('runfolder_method%i_dev%g',Deviation_Method_Switch,Deviation);
                mkdir (rundocfolder);
                rundoc = fullfile(rundocfolder,sprintf('rundoc_method%i_dev%g.bat',Deviation_Method_Switch,Deviation));
                extractname = sprintf('extract_method%i_dev%g.m',Deviation_Method_Switch,Deviation);
                extractdoc = fullfile(rundocfolder,extractname);
                fid = fopen(rundoc,'w');
                fid1 = fopen(extractdoc,'w');
            else
                %No stiffness vary, deviation, no tears, repeats
                rundocfolder = sprintf('runfolder_method%i_dev%g_rep%i',Deviation_Method_Switch,Deviation,Repeat_Count);
                mkdir (rundocfolder);
                rundoc = fullfile(rundocfolder,sprintf('rundoc_method%i_dev%g_rep%i.bat',Deviation_Method_Switch,Deviation,Repeat_Count));
                extractname = sprintf('extract_method%i_dev%g_rep%i.m',Deviation_Method_Switch,Deviation,Repeat_Count);
                extractdoc = fullfile(rundocfolder,extractname);
                fid = fopen(rundoc,'w');
                fid1 = fopen(extractdoc,'w');
            end
        else
            if(Repeat_Count==1)
                %No stiffness vary, deviation, tears, no repeats
                rundocfolder = sprintf('runfolder_method%i_dev%g_tears%g',Deviation_Method_Switch,Deviation,Tears);
                mkdir (rundocfolder);
                rundoc = fullfile(rundocfolder,sprintf('rundoc_method%i_dev%g_tears%g.bat',Deviation_Method_Switch,Deviation,Tears));
                extractname = sprintf('extract_method%i_dev%g_tears%g.m',Deviation_Method_Switch,Deviation,Tears);
                extractdoc = fullfile(rundocfolder,extractname);
                fid = fopen(rundoc,'w');
                fid1 = fopen(extractdoc,'w');
            else
                %No stiffness vary, deviation, tears, repeats
                rundocfolder = sprintf('runfolder_method%i_dev%g_tears%g_rep%i',Deviation_Method_Switch,Deviation,Tears,Repeat_Count);
                mkdir (rundocfolder);
                rundoc = fullfile(rundocfolder,sprintf('rundoc_method%i_dev%g_tears%g_rep%i.bat',Deviation_Method_Switch,Deviation,Tears,Repeat_Count));
                extractname = sprintf('extract_method%i_dev%g_tears%g_rep%i.m',Deviation_Method_Switch,Deviation,Tears,Repeat_Count);
                extractdoc = fullfile(rundocfolder,extractname);
                fid = fopen(rundoc,'w');
                fid1 = fopen(extractdoc,'w');
            end
        end
    end
else %stiffness vary != 0
    if(Deviation==0)
        if(Tears==0)
            if(Repeat_Count==1)
                %Stiffness vary !=0, deviation = 0, tears = 0, repeats = 1
                rundocfolder = sprintf('runfolder_method%i_stiffvary%g',Deviation_Method_Switch,StiffnessVary);
                mkdir (rundocfolder);
                rundoc = fullfile(rundocfolder,sprintf('rundoc_method%i_stiffvary%g.bat',Deviation_Method_Switch,StiffnessVary));
                extractname = sprintf('extract_method%i_stiffvary%g.m',Deviation_Method_Switch,StiffnessVary);
                extractdoc = fullfile(rundocfolder,extractname);
                fid = fopen(rundoc,'w');
                fid1 = fopen(extractdoc,'w');
            else %repeats != 1
                %Stiffness vary !=0, deviation = 0, tears = 0, repeats !=1
                rundocfolder = sprintf('runfolder_method%i_stiffvary%g_rep%i',Deviation_Method_Switch,StiffnessVary,Repeat_Count);
                mkdir (rundocfolder);
                rundoc = fullfile(rundocfolder,sprintf('rundoc_method%i_stiffvary%g_rep%i.bat',Deviation_Method_Switch,StiffnessVary,Repeat_Count));
                extractname = sprintf('extract_method%i_stiffvary%g_rep%i.m',Deviation_Method_Switch,StiffnessVary,Repeat_Count);
                extractdoc = fullfile(rundocfolder,extractname);
                fid = fopen(rundoc,'w');
                fid1 = fopen(extractdoc,'w');
            end
        else % tears != 0
            if(Repeat_Count==1)
                %Stiffness vary !=1, deviation = 0, tears !=0, repeats =1
                rundocfolder = sprintf('runfolder_method%i_tears%g_stiffvary%g',Deviation_Method_Switch,Tears,StiffnessVary);
                mkdir (rundocfolder);
                rundoc = fullfile(rundocfolder,sprintf('rundoc_method%i_tears%g_stiffvary%g.bat',Deviation_Method_Switch,Tears,StiffnessVary));
                extractname = sprintf('extract_method%i_tears%g_stiffvary%g.m',Deviation_Method_Switch,Tears,StiffnessVary);
                extractdoc = fullfile(rundocfolder,extractname);
                fid = fopen(rundoc,'w');
                fid1 = fopen(extractdoc,'w');
            else
                %Stiffness vary !=1, deviation = 0, tears !=0, repeats !=1
                rundocfolder = sprintf('runfolder_method%i_tears%g_stiffvary%g_rep%i',Deviation_Method_Switch,Tears,StiffnessVary,Repeat_Count);
                mkdir (rundocfolder);
                rundoc = fullfile(rundocfolder,sprintf('rundoc_method%i_tears%g_stiffvary%g_rep%i.bat',Deviation_Method_Switch,Tears,StiffnessVary,Repeat_Count));
                extractname = sprintf('extract_method%i_tears%g_stiffvary%g_rep%i.m',Deviation_Method_Switch,Tears,StiffnessVary,Repeat_Count);
                extractdoc = fullfile(rundocfolder,extractname);
                fid = fopen(rundoc,'w');
                fid1 = fopen(extractdoc,'w');
            end
        end
    else % deviation != 0
        if(Tears==0)
            if(Repeat_Count==1)
                %stiffness vary !=1, devation !=0, tears = 0, repeats =1
                rundocfolder = sprintf('runfolder_method%i_dev%g_stiffvary%g',Deviation_Method_Switch,Deviation,StiffnessVary);
                mkdir (rundocfolder);
                rundoc = fullfile(rundocfolder,sprintf('rundoc_method%i_dev%g_stiffvary%g.bat',Deviation_Method_Switch,Deviation,StiffnessVary));
                extractname = sprintf('extract_method%i_dev%g_stiffvary%g.m',Deviation_Method_Switch,Deviation,StiffnessVary);
                extractdoc = fullfile(rundocfolder,extractname);
                fid = fopen(rundoc,'w');
                fid1 = fopen(extractdoc,'w');
            else % repeats != 1, deviation !=0, tears = 0, repeats !=1
                rundocfolder = sprintf('runfolder_method%i_dev%g_stiffvary%g_rep%i',Deviation_Method_Switch,Deviation,StiffnessVary,Repeat_Count);
                mkdir (rundocfolder);
                rundoc = fullfile(rundocfolder,sprintf('rundoc_method%i_dev%g_stiffvary%g_rep%i.bat',Deviation_Method_Switch,Deviation,StiffnessVary,Repeat_Count));
                extractname = sprintf('extract_method%i_dev%g_stiffvary%g_rep%i.m',Deviation_Method_Switch,Deviation,StiffnessVary,Repeat_Count);
                extractdoc = fullfile(rundocfolder,extractname);
                fid = fopen(rundoc,'w');
                fid1 = fopen(extractdoc,'w');
            end
        else % tears != 0
            if(Repeat_Count==1)
                %stiffness vary !=0, deviation !=0, tears !=0, repeats =1
                rundocfolder = sprintf('runfolder_method%i_dev%g_tears%g_stiffvary%g',Deviation_Method_Switch,Deviation,Tears,StiffnessVary);
                mkdir (rundocfolder);
                rundoc = fullfile(rundocfolder,sprintf('rundoc_method%i_dev%g_tears%g_stiffvary%g.bat',Deviation_Method_Switch,Deviation,Tears,StiffnessVary));
                extractname = sprintf('extract_method%i_dev%g_tears%g_stiffvary%g.m',Deviation_Method_Switch,Deviation,Tears,StiffnessVary);
                extractdoc = fullfile(rundocfolder,extractname);
                fid = fopen(rundoc,'w');
                fid1 = fopen(extractdoc,'w');
            else % repeats != 1
                %stiffness vary !=0, deviation !=0, tears !=0, repeats !=1
                rundocfolder = sprintf('runfolder_method%i_dev%g_tears%g_stiffvary%g_rep%i',Deviation_Method_Switch,Deviation,Tears,StiffnessVary,Repeat_Count);
                mkdir (rundocfolder);
                rundoc = fullfile(rundocfolder,sprintf('rundoc_method%i_dev%g_tears%g_stiffvary%g_rep%i.bat',Deviation_Method_Switch,Deviation,Tears,StiffnessVary,Repeat_Count));
                extractname = sprintf('extract_method%i_dev%g_tears%g_stiffvary%g_rep%i.m',Deviation_Method_Switch,Deviation,Tears,StiffnessVary,Repeat_Count);
                extractdoc = fullfile(rundocfolder,extractname);
                fid = fopen(rundoc,'w');
                fid1 = fopen(extractdoc,'w');
            end
        end
    end
end

%%

addpath(rundocfolder);

for ii = 1:Repeat_Count
    
    % Generation of .FOLD File
    
    switch Deviation_Method_Switch
        case 0
            % 0 = No Deviation
            Deviation = 0;
            [filename,width] = NoDefects(Side_Length_A,Side_Length_B,Angle_Alpha,Angle_Theta,Count_M,Count_N,Deviation_Method_Switch,Deviation,ii,Repeat_Count,Tears,StiffnessVary,rundocfolder);
            displacement = -Displacement*width/1000; % Adjustment of sign convention to global coordinates and unit change
            FEAfromFOLD(filename,Meshsize,displacement,Count_M,Count_N,rundocfolder);
            
        case 1
            % 1 = Centre Point Deviation (LEGACY)
            [filename,width] = CentreDefects(Side_Length_A,Side_Length_B,Angle_Alpha,Angle_Theta,Count_M,Count_N,Deviation,ii,Repeat_Count,rundocfolder);
            displacement = -Displacement*width/1000; % Adjustment of sign convention to global coordinates and unit change
            FEAfromFOLD(filename,Meshsize,displacement,Count_M,Count_N,rundocfolder);
            
        case 2
            % 2 = All Points Deviation
            [filename,width] = NoDefects(Side_Length_A,Side_Length_B,Angle_Alpha,Angle_Theta,Count_M,Count_N,Deviation_Method_Switch,Deviation,ii,Repeat_Count,Tears,StiffnessVary,rundocfolder);
            FOLD_Distorter(filename,Side_Length_A,Side_Length_B,Count_M,Count_N,Deviation,rundocfolder);
            displacement = -Displacement*width/1000; % Adjustment of sign convention to global coordinates and unit change
            FEAfromFOLD(filename,Meshsize,displacement,Count_M,Count_N,rundocfolder);
            
        case 3
            % 3 = No Deviation with tears
            [filename,width] = NoDefects(Side_Length_A,Side_Length_B,Angle_Alpha,Angle_Theta,Count_M,Count_N,Deviation_Method_Switch,Deviation,ii,Repeat_Count,Tears,StiffnessVary,rundocfolder);
            displacement = -Displacement*width/1000; % Adjustment of sign convention to global coordinates and unit change
            TornFEAfromFOLD(filename,Meshsize,displacement,Count_M,Count_N,Tears,rundocfolder);
            
        case 4
            % 4 = Deviated with tears
            [filename,width] = NoDefects(Side_Length_A,Side_Length_B,Angle_Alpha,Angle_Theta,Count_M,Count_N,Deviation_Method_Switch,Deviation,ii,Repeat_Count,Tears,StiffnessVary,rundocfolder);
            FOLD_Distorter(filename,Side_Length_A,Side_Length_B,Count_M,Count_N,Deviation);
            displacement = -Displacement*width/1000; % Adjustment of sign convention to global coordinates and unit change
            TornFEAfromFOLD(filename,Meshsize,displacement,Count_M,Count_N,Tears,rundocfolder);
            
        case 5
            % 5 = No Deviation, fold stiffness varying
            [filename,width] = NoDefects(Side_Length_A,Side_Length_B,Angle_Alpha,Angle_Theta,Count_M,Count_N,Deviation_Method_Switch,Deviation,ii,Repeat_Count,Tears,StiffnessVary,rundocfolder);
            displacement = -Displacement*width/1000; % Adjustment of sign convention to global coordinates and unit change
            FEAfromFOLDVaryingStiffness(filename,Meshsize,displacement,Count_M,Count_N,StiffnessVary,rundocfolder);
            
        case 6
            % 6 = Deviation with varying fold stiffness
            [filename,width] = NoDefects(Side_Length_A,Side_Length_B,Angle_Alpha,Angle_Theta,Count_M,Count_N,Deviation_Method_Switch,Deviation,ii,Repeat_Count,Tears,StiffnessVary,rundocfolder);
            FOLD_Distorter(filename,Side_Length_A,Side_Length_B,Count_M,Count_N,Deviation,rundocfolder);
            displacement = -Displacement*width/1000; % Adjustment of sign convention to global coordinates and unit change
            FEAfromFOLDVaryingStiffness(filename,Meshsize,displacement,Count_M,Count_N,StiffnessVary,rundocfolder);
            
        case 7
            % 7 = Tears with varying fold stiffness no Deviation
            [filename,width] = NoDefects(Side_Length_A,Side_Length_B,Angle_Alpha,Angle_Theta,Count_M,Count_N,Deviation_Method_Switch,Deviation,ii,Repeat_Count,Tears,StiffnessVary,rundocfolder);
            displacement = -Displacement*width/1000; % Adjustment of sign convention to global coordinates and unit change
            TornFEAfromFOLDVaryingStiffness(filename,Meshsize,displacement,Count_M,Count_N,Tears,rundocfolder);
            
        case 8
            % 8 = Deviation, tears and varying fold stiffness
            [filename,width] = NoDefects(Side_Length_A,Side_Length_B,Angle_Alpha,Angle_Theta,Count_M,Count_N,Deviation_Method_Switch,Deviation,ii,Repeat_Count,Tears,StiffnessVary,rundocfolder);
            FOLD_Distorter(filename,Side_Length_A,Side_Length_B,Count_M,Count_N,Deviation);
            displacement = -Displacement*width/1000; % Adjustment of sign convention to global coordinates and unit change
            TornFEAfromFOLDVaryingStiffness(filename,Meshsize,displacement,Count_M,Count_N,Tears,rundocfolder);
    end
    
    %if Repeat_Count ~= 1
    fprintf(fid,'call abaqus job=%s interactive',filename);
    fprintf(fid,'\n');
    %end
    
    if ii ==1
        fprintf(fid1,'filename = cell(%i,1);',Repeat_Count);
        fprintf(fid1,'\n');
        fprintf(fid1,'rundocfolder = cell(%i,1);',Repeat_Count);
        fprintf(fid1,'\n');
    end
        
    fprintf(fid1,'filename{%i} = ''%s.dat'';',ii,filename);
    fprintf(fid1,'\n');
    fprintf(fid1,'rundocfolder{%i} = ''%s'';',ii,rundocfolder);
    fprintf(fid1,'\n');
    
    if ii == Repeat_Count
        
        fprintf(fid1,'started = 0;');
        fprintf(fid1,'\n');
        fprintf(fid1,'for ii=1 : length(filename)');
        fprintf(fid1,'\n');
        fprintf(fid1,'if isfile(filename{ii})');
        fprintf(fid1,'\n');
        fprintf(fid1,'values = data_extract(rundocfolder{ii},filename{ii});');
        fprintf(fid1,'\n');
        fprintf(fid1,'if ismatrix(values) && started == 0');
        fprintf(fid1,'\n');
        fprintf(fid1,'timeenergy = [values];');
        fprintf(fid1,'\n');
        fprintf(fid1,'started = 1;');
        fprintf(fid1,'\n');
        fprintf(fid1,'elseif ismatrix(values) && started == 1');
        fprintf(fid1,'\n');
        fprintf(fid1,'timeenergy = [timeenergy,values];');
        fprintf(fid1,'\n');
        fprintf(fid1,'end');
        fprintf(fid1,'\n');
        fprintf(fid1,'end');
        fprintf(fid1,'\n');
        fprintf(fid1,'end');
        fprintf(fid1,'\n');
        fprintf(fid1,'if started == 1');
        fprintf(fid1,'\n');
        fprintf(fid1,'writematrix(timeenergy,strcat(rundocfolder{1},''.xls''));');
        fprintf(fid1,'\n');
        fprintf(fid1,'end');
        fprintf(fid1,'\n');
        fprintf(fid1,'clear all');
        fprintf(fid1,'\n');
    
    end

end
fclose('all');
%