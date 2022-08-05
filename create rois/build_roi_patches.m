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
Region = {'S1','M1'};
roiSize = 5;
cd(surfPATH)
k = dir('1*'); subj={k.name}'; clear k
vis = 1; % if visualisation on 

%% ROI - CREATE ROI PATCH

% loop over subjects
for s = 1:length(subj)
    
    % get subj. path and load surface
    SUBJ_PATH = [surfPATH, sprintf('%s\%s\surf\', subj{sub},subj{sub})]; %subject folder within surface path   
    sl = gifti([SUBJ_PATH 'ro_white-pial.gii']);
    
    % loop over Regions (S1,M1)
    for i = 1:length(Region)
        % load coordinates
        load([SUBJ_PATH, Region{i}, 'coord_bank.mat']); %Coordinates for bank ROI (e.g. M1 bank)
        ROI_centre_coord_approx_bank = ROI_centre_coord_approx;clear ROI_centre_coord_approx
        load([SUBJ_PATH, Region{i}, 'coord_crown.mat']); %Coordinates for crown ROI (e.g. M1 crown)
        ROI_centre_coord_approx_crown = ROI_centre_coord_approx;clear ROI_centre_coord_approx
        
        % find index of closest vertex for ROI centre coords
        ROI_centre_ind_bank = knnsearch(sl.vertices,ROI_centre_coord_approx_bank);
        ROI_centre_ind_crown = knnsearch(sl.vertices,ROI_centre_coord_approx_crown);
        
        % initialise ROI_vert_ind
        ROI_vert_ind_bank = [ROI_centre_ind_bank];
        ROI_vert_ind_crown = [ROI_centre_ind_crown];
        
        % loop over roisize iterations
        for r = 1:roiSize % defines the 'size' of ROI patch (approximate at the minute. Finds neighbouring vertex 'r' times (i.e. find neighbour of centre, find neighbour of neighbour etc)
            % BUILD CROWN
            ROI_vert_ind_crown_temp = [];
            
            % for each vertex in ROI_vert_ind
            for v = 1:length(ROI_vert_ind_crown)
                % get the indices of all faces this vertex is part of
                [face_ind,~] = ind2sub(size(sl.faces),find(sl.faces == ROI_vert_ind_crown(v)));
                
                % for each face
                for f = 1:length(face_ind)
                    % get other two verticies that are part of the face
                    ROI_vert_ind_crown_temp = [ROI_vert_ind_crown_temp,sl.faces(face_ind(f),:)];
                end % f
            end % v
            
            % BUILD BANK
            ROI_vert_ind_bank_temp = [];
            
            % for each vertex in ROI_vert_ind
            for v = 1:length(ROI_vert_ind_bank)
                % get the indices of all faces this vertex is part of
                [face_ind,~] = ind2sub(size(sl.faces),find(sl.faces == ROI_vert_ind_bank(v)));
                
                % for each face
                for f = 1:length(face_ind)
                    % get other two verticies that are part of the face
                    ROI_vert_ind_bank_temp = [ROI_vert_ind_bank_temp,sl.faces(face_ind(f),:)];
                end % f
            end % v
            
            ROI_vert_ind_crown = [ROI_vert_ind_crown,ROI_vert_ind_crown_temp];
            ROI_vert_ind_crown = unique(ROI_vert_ind_crown);
            
            ROI_vert_ind_bank = [ROI_vert_ind_bank,ROI_vert_ind_bank_temp];
            ROI_vert_ind_bank = unique(ROI_vert_ind_bank);
        end % r
        
        %plot original patches if vis = 1
        if vis
            f1 = figure('units','normalized','outerposition',[0 0 1 1]);hold on;axis equal;
            t = patch('Faces',sl.faces,'Vertices',sl.vertices,'FaceVertexCData',ones(size(sl.vertices,1),1),'EdgeColor',[0.8 0.8 0.8],'FaceColor','flat');
            colormap([0.95 0.95 0.95]);
            scatter3(sl.vertices(ROI_vert_ind_bank,1),sl.vertices(ROI_vert_ind_bank,2),sl.vertices(ROI_vert_ind_bank,3),60,'o','MarkerFacecolor','b','MarkerEdgeColor','none');
            scatter3(sl.vertices(ROI_vert_ind_crown,1),sl.vertices(ROI_vert_ind_crown,2),sl.vertices(ROI_vert_ind_crown,3),60,'o','MarkerFacecolor','r','MarkerEdgeColor','none');
            xlim([ROI_centre_coord_approx_bank(1)-10 ROI_centre_coord_approx_bank(1)+10]);ylim([ROI_centre_coord_approx_bank(2)-10 ROI_centre_coord_approx_bank(2)+10]);zlim([ROI_centre_coord_approx_bank(3)-10 ROI_centre_coord_approx_bank(3)+10]);
            savefig([SUBJ_PATH Region{i} '_pre.fig']); 
        end
        
        % check if crown and bank patch are overlapping
        overlap = intersect(ROI_vert_ind_crown,ROI_vert_ind_bank);
        
        % if there is overlap: solve 
        if length(overlap)>0
            
            if vis
                % add overlap (also add original center points again, in
                % case they were overdrawn)
                scatter3(sl.vertices(overlap,1),sl.vertices(overlap,2),sl.vertices(overlap,3),80,'o','MarkerFacecolor','g','MarkerEdgeColor','none');
                scatter3(sl.vertices(ROI_centre_ind_bank,1),sl.vertices(ROI_centre_ind_bank,2),sl.vertices(ROI_centre_ind_bank,3),60,'o','MarkerFacecolor','k','MarkerEdgeColor','none');
                scatter3(sl.vertices(ROI_centre_ind_crown,1),sl.vertices(ROI_centre_ind_crown,2),sl.vertices(ROI_centre_ind_crown,3),60,'o','MarkerFacecolor','k','MarkerEdgeColor','none');
                savefig([SUBJ_PATH Region{i} '_pre.fig']);
            end
            
            % for each overlapping vertex check if the distance (geodesic)
            % is closest to crown or bank 
            clear sl_small org2small
            sl_small.vertices = double(sl.vertices(unique([ROI_vert_ind_crown,ROI_vert_ind_bank]),:));
            parfor s = 1:size(sl.vertices,1)
                [~, org2small(s)] = ismember(sl.vertices(s,:),sl_small.vertices,'rows');
            end
            sl_small.faces(:,1) = org2small(sl.faces(:,1));sl_small.faces(:,2) = org2small(sl.faces(:,2));sl_small.faces(:,3) = org2small(sl.faces(:,3));
            sl_small.faces = sl_small.faces(all(sl_small.faces,2),:);
            
            % initialise mesh and algorithm for geodesic
            mesh = geodesic_new_mesh(sl_small.vertices,sl_small.faces);
            algorithm = geodesic_new_algorithm(mesh, 'dijkstra');
            
            % get geodesic distance from each centre (bank and crown) to
            % each vertex
            clear distances
            f2 = figure('units','normalized','outerposition',[0 0 1 1]);
            for y = 1:2
                % get index of original center coords in small mesh
                if y==1
                    [~, LOCB] = ismember(sl.vertices(ROI_centre_ind_crown,:),sl_small.vertices,'rows');
                else
                    [~, LOCB] = ismember(sl.vertices(ROI_centre_ind_bank,:),sl_small.vertices,'rows');
                end
                source_points = {geodesic_create_surface_point('vertex',LOCB,sl_small.vertices(LOCB,:))};
                geodesic_propagate(algorithm, source_points);   %propagation stage of the algorithm (the most time-consuming)
                distances{y} = zeros(size(sl_small.vertices,1),1);              %find distances to all vertices of the mesh (actual pathes are not computed)
                [source_id, distances{y}] = geodesic_distance_and_source(algorithm);     %find distances to all vertices of the mesh; in this example we have a single source, so source_id is always equal to 1
            
                if vis
                    % plot distance
                    subplot(1,2,y);axis equal;hold on;
                    t = patch('Faces',sl_small.faces,'Vertices',sl_small.vertices,'FaceVertexCData',distances{y},'EdgeColor','none','FaceColor','flat');
                    scatter3(sl_small.vertices(LOCB,1),sl_small.vertices(LOCB,2),sl_small.vertices(LOCB,3),'ko','filled');
                    colormap parula;colorbar;
                    savefig([SUBJ_PATH Region{i} '_dist.fig']);
                end
            end
            
            % solve overlap (if vertex is closer to center of bank, remove from
            % crown, and vice versa)
            for o = 1:length(overlap)
                % get index of original coords in small mesh
                [~, LOCB] = ismember(sl.vertices(overlap(o),:),sl_small.vertices,'rows');
                
                if distances{2}(LOCB) < distances{1}(LOCB)
                    ROI_vert_ind_crown(ROI_vert_ind_crown==overlap(o)) = [];
                    if vis
                        % re-colour vertex from overlap to newly assigned
                        % patch
                        figure(f1);
                        scatter3(sl.vertices(overlap(o),1),sl.vertices(overlap(o),2),sl.vertices(overlap(o),3),50,'o','MarkerFacecolor','b','MarkerEdgeColor','none');
                    end
                else
                    ROI_vert_ind_bank(ROI_vert_ind_bank==overlap(o)) = [];
                    if vis
                        figure(f1);
                        scatter3(sl.vertices(overlap(o),1),sl.vertices(overlap(o),2),sl.vertices(overlap(o),3),50,'o','MarkerFacecolor','r','MarkerEdgeColor','none');
                    end
                end
            end
            
            if vis
                %(add original center points again, in case they were overdrawn)
                scatter3(sl.vertices(ROI_centre_ind_bank,1),sl.vertices(ROI_centre_ind_bank,2),sl.vertices(ROI_centre_ind_bank,3),60,'o','MarkerFacecolor','k','MarkerEdgeColor','none');
                scatter3(sl.vertices(ROI_centre_ind_crown,1),sl.vertices(ROI_centre_ind_crown,2),sl.vertices(ROI_centre_ind_crown,3),60,'o','MarkerFacecolor','k','MarkerEdgeColor','none');
                savefig([SUBJ_PATH Region{i} '_post.fig']);
            end
            
            % double check
            overlap2 = intersect(ROI_vert_ind_crown,ROI_vert_ind_bank);
            if length(overlap2)>0
                keyboard;
            end
            
            %keyboard;
            close(f2);
        end % if there is overlap
        
        ROI_vert_ind = ROI_vert_ind_crown; save([SUBJ_PATH, Region{i}, '_crown_patch.mat'],'ROI_vert_ind');
        ROI_vert_ind = ROI_vert_ind_bank; save([SUBJ_PATH, Region{i}, '_bank_patch.mat'],'ROI_vert_ind');
            
        close(f1);
        patch_size(s,i,1) = length(ROI_vert_ind_bank);
        patch_size(s,i,2) = length(ROI_vert_ind_crown);
        
    end % region (S1/M1)
end % subject
