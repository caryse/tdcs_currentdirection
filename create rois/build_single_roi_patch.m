%% Create ROIs on GM surface
% (c) Carys Evans & Catharina Zich, UCL 
% carys.evans@ucl.ac.uk
% July 2022

%This script:
% Builds bank and crown ROI patch at the same time using ROI centre
% coordinates taken from roi_centre_coords.m
% Removes overlap between ROI patches

%Script also requires library for calculating geodesic distance:
% library for geodesic distance
% https://uk.mathworks.com/matlabcentral/fileexchange/18168-exact-geodesic-for-triangular-meshes
% download from here, add to path and run example1

%% Define directories and parameters

global geodesic_library;
geodesic_library = 'geodesic_debug';      %"release" is faster and "debug" does additional checks
addpath(genpath('\PATH\TO\TOOLBOX\geodesic_matlab-master\')); %geodesic toolbox
surfPATH = 'D:\PATH\TO\SURFACE\DATA\';         %surface path
cfmPATH = 'D:\PATH\TO\ROAST\MODELS\';          %ROAST models

% params
%Region = {'S1','M1'};
Region = {'M1'};
roiSize = 5;
cd(surfPATH)
%k = dir('1*'); subj={k.name}'; clear k
slashIdx = strfind(surfPATH, '/');
subj = extractBetween(surfPATH, (slashIdx(1,end-1)+1),(slashIdx(1,end)-1));
%subj = cell2mat(extractBetween(surfPATH, (slashIdx(1,end-1)+1),(slashIdx(1,end)-1)));
clear slashIdx

vis = 1; % 1 = visualise patches; 0 = don't visualise patches

%% ROI - CREATE ROI PATCH

% loop over subjects
for s = 1:length(subj)
    
  % get subj. path and load surface
  %  SUBJ_PATH = [surfPATH, sprintf('%s\surf\',subj{s})]; %subject folder within surface path   
    SUBJ_PATH = [surfPATH, 'surf/']; %subject folder within surface path   
    sl = gifti([SUBJ_PATH 'ro_white-pial.gii']);
    
    % loop over Regions (e.g. M1/S1)
    for i = 1:length(Region)
        % load coordinates
        cd(SUBJ_PATH)
        k=dir([str2mat(Region),'*coord.mat']); myroi = k.name;
        load([SUBJ_PATH, myroi]); %Coordinates for bank ROI (e.g. M1 bank)
        
        % find index of closest vertex for ROI centre coords
        ROI_centre_ind = knnsearch(sl.vertices,ROI_centre_coord_approx);

        % initialise ROI_vert_ind
        ROI_vert_ind = [ROI_centre_ind];

        % loop over roisize iterations
        for r = 1:roiSize % defines the 'size' of ROI patch (approximate at the minute. Finds neighbouring vertex 'r' times (i.e. find neighbour of centre, find neighbour of neighbour etc)
           % BUILD PATCH
            ROI_vert_ind_temp = [];
            
            % for each vertex in ROI_vert_ind
            for v = 1:length(ROI_vert_ind)
                % get the indices of all faces this vertex is part of
                [face_ind,~] = ind2sub(size(sl.faces),find(sl.faces == ROI_vert_ind(v)));
                
                % for each face
                for f = 1:length(face_ind)
                    % get other two verticies that are part of the face
                    ROI_vert_ind_temp = [ROI_vert_ind_temp,sl.faces(face_ind(f),:)];
                end % f
            end % v            
            
            ROI_vert_ind = [ROI_vert_ind,ROI_vert_ind_temp];
            ROI_vert_ind = unique(ROI_vert_ind);
            
        end % r
        
        %plot original patches if vis = 1
        if vis
            f1 = figure('units','normalized','outerposition',[0 0 1 1]);hold on;axis equal;
            t = patch('Faces',sl.faces,'Vertices',sl.vertices,'FaceVertexCData',ones(size(sl.vertices,1),1),'EdgeColor',[0.8 0.8 0.8],'FaceColor','flat');
            colormap([0.95 0.95 0.95]);
            scatter3(sl.vertices(ROI_vert_ind,1),sl.vertices(ROI_vert_ind,2),sl.vertices(ROI_vert_ind,3),60,'o','MarkerFacecolor','b','MarkerEdgeColor','b');
 %           xlim([ROI_centre_coord_approx(1)-10 ROI_centre_coord_approx(1)+10]);ylim([ROI_centre_coord_approx(2)-10 ROI_centre_coord_approx(2)+10]);zlim([ROI_centre_coord_approx(3)-10 ROI_centre_coord_approx(3)+10]);
            savefig([SUBJ_PATH Region{i} '_pre.fig']); 
            
        end
  
        ROI_vert_ind = ROI_vert_ind; save([SUBJ_PATH, Region{i}, '_patch.mat'],'ROI_vert_ind');
            
    %    close(f1);
        patch_size(s,i,1) = length(ROI_vert_ind);
        
    end % region (S1/M1)
end % subject
