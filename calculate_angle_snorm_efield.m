%% Extract surface norms and E-field norms from CFM
% (c) Carys Evans, UCL
% carys.evans@ucl.ac.uk
% July 2022

%This script:
%Creates GM surface and extracts surface norm and E-field direction
%Loads ROI patch and extracts E-field
%Calculates DISTANCE between snorm and EF norm
%Calculates ANGLE between snorm and EF norm
%Saves relevant data to subjectID_surfacedata.mat file

%% Define file names and directories

%set directories
% addpath(genpath('C:\Matlab2018b\spm12\'));           %spm
% surfPATH = 'D:\SURFACE\DATA\';                       %surface pathway
% cfmPATH = 'D:\ROAST\MODELS\';                        %roast models
% ROIPATH = 'D:\ROI\DATA\';                            %ROI pathway

addpath(genpath('/Users/carysevans/Documents/MATLAB/spm12/'));           %spm
surfPATH = '/Users/carysevans/Desktop/directionWlkthru/Joyce/';                       %surface pathway
cfmPATH = '/Users/carysevans/Desktop/directionWlkthru/tDCS_simulation/';                        %roast models
ROIPATH = '/Users/carysevans/Desktop/directionWlkthru/Joyce/';                            %ROI pathway

%filenames
%Roast result after 'subjectID_'
cfmResultFile = {'CP3FCZ_2mA_roastResult.mat'};

%ROI
ROIfile = 'M1coord.mat'; %centre coordinates
ROIpatch = 'M1_patch.mat'; %ROI patch (e.g. M1 bank)
%ROI + MONTAGE label - also becomes result output filename
ROIlabel = {'M1_CP3FCZ'}; %M1_CP5FC1

roiSize = 5;    %ROI size
%surface filename
gmfilename = 'ro_white-pial'; %surface filename


%SUBJECT MRI.NII
%=================
%subjects
cd(surfPATH)
slashIdx = strfind(surfPATH, '/');
subj = extractBetween(surfPATH, (slashIdx(1,end-1)+1),(slashIdx(1,end)-1));
%subj = cell2mat(extractBetween(surfPATH, (slashIdx(1,end-1)+1),(slashIdx(1,end)-1)));
clear slashIdx

%subjects name in roast model
cd(cfmPATH) 
k = dir(['*',cell2mat(cfmResultFile)]);
subj_cfm = {extractBefore(k.name,['_',cell2mat(cfmResultFile)])};
clear k



%% ==============================================
% CREATE GM SURFACE & EXTRACT SURFACE NORM AND E-FIELD DIRECTION
%%===============================================

for mont = 1:length(cfmResultFile)
    for sub = 1:length(subj_cfm)
        
        savefileName = sprintf('%s_surfacedata_%s.mat', subj{sub},ROIlabel{mont});
        
        %direct = [ROIPATH, sprintf('%s/', subj{sub})];
        if ~exist(ROIPATH, 'dir')
            continue
        else
        cd(ROIPATH);
        end
        
        if isfile(savefileName)
            fprintf('%s skipped. Already done \n', subj_cfm{sub});
            clearvars -except surfPATH cfmPATH ROIPATH mont subj subj_cfm cfmResultFile ROIpatch ROIfile roiSize ROIlabel filename
            continue
        else
        end
        
        %E-FIELD RESULT.MAT from Roast
        %=================
        cd(cfmPATH)
        %cfmData = sprintf('%s_CP5FC1_2mA_result.mat', subj_cfm{sub});
        cfmData = [sprintf('%s_', subj_cfm{sub}), cfmResultFile{mont}];
        if isfile(cfmData)
            load([cfmPATH, sprintf('%s_', subj_cfm{sub}),cfmResultFile{mont}]);
            %load([cfmPATH,cfmData])
        else
            fprintf('%s skipped. No ROAST model \n', subj_cfm{sub});
            clearvars -except surfPATH cfmPATH ROIPATH mont subj subj_cfm cfmResultFile ROIpatch ROIfile roiSize ROIlabel gmfilename
            continue
        end
              
        %MRI .nii volume (ras transformed by roast)
        if  isfile([cfmPATH, sprintf('%s_ras.nii',subj_cfm{sub})])
            ef_vol = spm_vol([cfmPATH, sprintf('%s_ras.nii',subj_cfm{sub})]);
        else
            ef_vol = spm_vol([cfmPATH, sprintf('%s.nii',subj_cfm{sub})]);
        end

        
        %GREY MATTER SURFACE.GII
        %=================
        surface = [surfPATH, sprintf('surf/%s.gii', gmfilename)];
        
        %generate gifti from surface file
        sl=gifti(surface);
        
        %plot(sl)
        
        %% Calculate surface norms from grey matter surface
        
        %calculate norms
        snorm=spm_mesh_normals(sl);
        
        %% Transform Roast voxel indices into Surface space
        
        %Step 1) get i,j,k coords from ef_all voxel matrix
        %==============================================
        %only uses non-NaNs from 'ef_all' from here on...
        
        %finds coordinates in ef_all where EF is present (vectors ~= NaN)
        %outputs ijk which are equivalent to xyz but in ef_all coordinate system (voxel indices)
        [ef_i,ef_j,ef_k] = ind2sub([size(ef_all,1),size(ef_all,2),size(ef_all,3)],find(~isnan(ef_all(:,:,:,1))));
        %keep NaNs:
        %[ef_i,ef_j,ef_k] = ind2sub([size(ef_all,1),size(ef_all,2),size(ef_all,3)],find(ef_all(:,:,:,1)));
        
        %combines ijk into 1 matrix
        %size variesperson to person because removed NaN
        ef_ijk = [ef_i,ef_j,ef_k];
        
        %cleanup: clear separate variables i j k
        clear ef_i ef_j ef_k
        
        %Step 2) transform ijk coords to xyz coords
        %==============================================
        %ijk = voxel; xyz = MNI; needs subject MRI
        
        %transform coordinates 
        ef_xyz = ef_vol.mat(1:3,1:3) * ef_ijk' + ef_vol.mat(1:3,4);
        %transpose result
        ef_xyz = ef_xyz';
  
        %Step 3) for each sl.vertices get closest ef vector
        %==============================================
        
        %find nearest neighbour in X (ef_xyz) for each query point in Y (sl.vertices)
        efield_surfmap=knnsearch(ef_xyz,sl.vertices);
        
        %find coordinates for each ef vector using roast's voxel space
        efsurf_coord=ef_ijk(efield_surfmap,:);
        
        %get ef vector from ef_all
        for v_in = 1:length(efsurf_coord)
            % get ef vector information from ef_all (4th dim of ef_all)
            ef_vec(v_in,:) = ef_all(efsurf_coord(v_in,1),efsurf_coord(v_in,2),efsurf_coord(v_in,3),:);
        end
        
        clear v_in
        
        %% Calculate E-field vectors on surface norm coordinates
        
        for v_in = 1:length(efsurf_coord)
            % get ef mag for ef vectors
            ef_vec_mag(v_in,:) = ef_mag(efsurf_coord(v_in,1),efsurf_coord(v_in,2),efsurf_coord(v_in,3));
        end
        
        %calculate ef normal component (magnitude removed)
        ef_norm = ef_vec./sqrt(sum(ef_vec.^2,2));
        
        clear v_in
        
        %% ==============================================
        % LOAD ROI PATCH & EXTRACT E-FIELD IN ROI
        %%===============================================
         
        %LOAD ROI CENTRE COORDS
        %=================
        
        %cd([ROIPATH, sprintf('%s/', subj{sub})]);
        cd(ROIPATH);

        if isfile(ROIfile)
            ROI_centre_coord_approx = importdata(ROIfile);
        else %pause;
            fprintf('%s skipped. No ROI file \n', subj_cfm{sub});
            clearvars -except surfPATH cfmPATH ROIPATH subj subj_cfm cfmResultFile ROIpatch ROIfile roiSize ROIlabel gmfilename
            continue
        end
        
        % find index of closest vertex for ROI centre coords
        ROI_centre_ind = knnsearch(sl.vertices,ROI_centre_coord_approx);
        
        %ROI coords
        %=================       
        ROI_vert_ind = importdata(ROIpatch);        
        
        %% ROI - Find SNORM & EF in ROI        
        %%to reverse snorm (so arrows point outward from surface): snormRev = snorm*-1;
        
        %find coordinates of snorm in ROI in surface space
        snormROI_coord = sl.vertices(ROI_vert_ind,:);
        
        %find vectors of snorm in ROI
        snormROI = snorm(ROI_vert_ind,:);
        
        %find index of EF in ROI patch (result e.g. 126x1)
        efield_ROImap=knnsearch(ef_xyz,sl.vertices(ROI_vert_ind,:));
        
        %find coordinates for these EF indices in roast voxel space (result e.g. 126x3)
        efROI_coord=ef_ijk(efield_ROImap,:);
        
        %find coordinates for these EF indices in surface space (result e.g. 126x3)
        efROI_coordxyz=ef_xyz(efield_ROImap,:);        
        
        %% ROI - Calculate E-field vectors on surface norm coordinates
        
        %get ef vector from ef_all
        for v_in = 1:length(efROI_coord)
            
            % get ef vector information from ef_all (4th dim of ef_all), i.e. ef
            % vectors including magnitude
            efvec_inROI(v_in,:) = ef_all(efROI_coord(v_in,1),efROI_coord(v_in,2),efROI_coord(v_in,3),:);
        end
        
        %calculate ef unit vector (i.e. remove magnitude)
        ef_norm_inROI = efvec_inROI./sqrt(sum(efvec_inROI.^2,2));
        
        %extract ef_mag from ef vectors
        for v_in = 1:length(efROI_coord)
            % get ef mag for ef vectors
            ef_vec_mag_inROI(v_in,:) = ef_mag(efROI_coord(v_in,1),efROI_coord(v_in,2),efROI_coord(v_in,3));
        end
        
        %% ROI - Mean snorm / EF vectors in ROI
        
        %mean EF vector in ROI (includes magnitude info, efvec)
        meanEFVecinROI = mean(efvec_inROI);
        
        %mean EF unit vector in ROI (no magnitude, efnorm)
        meanEFNORMinROI = mean(ef_norm_inROI);
        
        %mean snorm in ROI
        meanSNORMinROI = mean(snorm(ROI_vert_ind,:)); %meanSNORMinROI = mean(snormROI)
        
        %mean EF magnitude in ROI
        meanEFMAGinROI = sqrt((meanEFVecinROI(:,1).^2)+(meanEFVecinROI(:,2).^2)+(meanEFVecinROI(:,3).^2));
        
        %% ==============================================
        % CALCULATE DISTANCE BETWEEN SNORM AND EF NORM
        %%===============================================
        
        %% GREY MATTER - Calculate DISTANCE between surface norm and NORM E-field vectors (Direction: snorms to e-field)
        %Distance = sqrt(x1-x2)2 + (y1-y2)2 + (z1-z2)2
        
        %calculate distance between snorm and ef
        for i=1:size(ef_norm,1)
            
            SE_Distance(i,1) = sqrt((snorm(i,1)- ef_norm(i,1)).^2 ...
                +(snorm(i,2)-ef_norm(i,2)).^2 ...
                +(snorm(i,3)-ef_norm(i,3)).^2);
            
        end
        
        clear i
        
        %INCLUDE MAGNITUDE        
        %maintain magnitude
        snorm_efmag = ef_vec;
        
        for i=1:size(ef_vec,1)            
            SE_DistanceMAG(i,1) = sqrt((snorm(i,1)- ef_vec(i,1)).^2 ...
                +(snorm(i,2)-ef_vec(i,2)).^2 ...
                +(snorm(i,3)-ef_vec(i,3)).^2);
        end        
        
        %% ROI - CALCULATE DISTANCE between surface norm and NORM E-field vectors (Direction: snorms to e-field)
        %Distance = sqrt(x1-x2)2 + (y1-y2)2 + (z1-z2)2
        
        % ========================
        % WHOLE ROI
        % ========================
        
        %calculate distance between snorm and ef vector - whole ROI
        for i=1:size(snormROI,1)
            
            SE_DistanceROI(i,1) = sqrt((snormROI(i,1)- efvec_inROI(i,1)).^2 ...
                +(snormROI(i,2)-efvec_inROI(i,2)).^2 ...
                +(snormROI(i,3)-efvec_inROI(i,3)).^2);
        end
        clear i
        
        %calculate distance between snorm and ef vector - whole ROI (NO MAGNITUDE)
        for i=1:size(snormROI,1)
            
            SE_DistanceROInorm(i,1) = sqrt((snormROI(i,1)- ef_norm_inROI(i,1)).^2 ...
                +(snormROI(i,2)-ef_norm_inROI(i,2)).^2 ...
                +(snormROI(i,3)-ef_norm_inROI(i,3)).^2);
        end
        clear i
        
        % ========================
        % MEANS W/IN ROI
        % ========================
        
        %MEANS - calculate distance between snorm and ef vector
        %for i=1:size(ef_norm,1)
        i = 1;
        MeanSE_Distance_inROI(i,1) = sqrt((meanSNORMinROI(i,1)- meanEFVecinROI(i,1)).^2 ...
            +(meanSNORMinROI(i,2)-meanEFVecinROI(i,2)).^2 ...
            +(meanSNORMinROI(i,3)-meanEFVecinROI(i,3)).^2);
        clear i
        
        %end
        
        %MEANS(NO MAGNITUDE) - calculate distance between snorm and efnorm
        
       % for i=1:size(ef_norm,1)
           i = 1;
            MeanSE_Distance_inROInorm(i,1) = sqrt((meanSNORMinROI(i,1)- meanEFNORMinROI(i,1)).^2 ...
                +(meanSNORMinROI(i,2)-meanEFNORMinROI(i,2)).^2 ...
                +(meanSNORMinROI(i,3)-meanEFNORMinROI(i,3)).^2);            
   %     end
        
        clear i
        
        %% ==============================================
        % CALCULATE ANGLE BETWEEN SNORM AND EF NORM
        %%===============================================
        
        %% GREY MATTER - Calculate ANGLE between snorm and ef vectors: cos(?) = a · b / (|a| x |b|)
        
        % a · b = snorm · ef_vec (surface norm vectors & ef vectors)
        normDotP = dot(snorm,ef_vec,2); %dot product: snorm · ef_vec
        
        %|a| = magnitude of snorm
        snorm_mag(:,1) = sqrt((snorm(:,1).^2)+(snorm(:,2).^2)+(snorm(:,3).^2));
        
        %|b| = magnitude of ef vectors
        ef_vec_mag; %efnorm_efmag;
        
        %cos(?). NB:snorm_mag is always 1, so (snorm_mag.*efnorm_efmag) is equivalent to efnorm_efmag
        cos_ang = (normDotP./(snorm_mag.*ef_vec_mag));
        
        %? = angle between snorm and ef vectors. acosd = inverse cosine in degrees
        angle = acosd(cos_ang);
        
        %angle in radians rather than degrees (useful for Polar Plots)
        anglerad = acos(cos_ang);
                
        %% ROI - Calculate angle between meanSNORM and meanEFNORM vectors:  cos(?) = a · b / (|a| x |b|)
        
        % ========================
        % WHOLE ROI
        % ========================
        
        % a · b = snorm · ef_vec (surface norm vectors & ef vectors)
        normDotP_ROI = dot(snormROI,efvec_inROI,2); %dot product: snorm · ef_norm (could also do ef_vec)
        
        %|a| = magnitude of snorm
        snormROI_mag= sqrt((snormROI(:,1).^2)+(snormROI(:,2).^2)+(snormROI(:,3).^2));        
        
        %|b| = magnitude of ef vectors
        efROI_mag= sqrt((efvec_inROI(:,1).^2)+(efvec_inROI(:,2).^2)+(efvec_inROI(:,3).^2));
        
        %cos(?). NB:snorm_mag is always 1, so (snorm_mag.*efnorm_efmag) is equivalent to efnorm_efmag
        cos_angROI = (normDotP_ROI./(snormROI_mag.*efROI_mag));
        
        %? = angle between snorm and ef vectors. acosd = inverse cosine in degrees
        ROIangle = acosd(cos_angROI);
        
        %angle in radians rather than degrees (useful for Polar Plots)
        ROIanglerad = acos(cos_angROI);
        
        % ========================
        % MEANS W/IN ROI
        % ========================
        
        % a · b = snorm · ef_vec (surface norm vectors & ef vectors)
        normDotP_meanROI = dot(meanSNORMinROI,meanEFNORMinROI,2); %dot product: snorm · ef_norm (could also do ef_vec)
        
        %|a| = magnitude of snorm
        snormROImean_mag(:,1)= sqrt((meanSNORMinROI(:,1).^2)+(meanSNORMinROI(:,2).^2)+(meanSNORMinROI(:,3).^2));
        
        %|b| = magnitude of ef vectors
        efnormROImean_mag(:,1)= sqrt((meanEFNORMinROI(:,1).^2)+(meanEFNORMinROI(:,2).^2)+(meanEFNORMinROI(:,3).^2));
        
        %cos(?). NB:snorm_mag is always 1, so (snorm_mag.*efnorm_efmag) is equivalent to efnorm_efmag
        cos_ang_meanROI = (normDotP_meanROI./(snormROImean_mag.*efnormROImean_mag));
        
        %? = angle between snorm and ef vectors. acosd = inverse cosine in degrees
        MeanROIangle = acosd(cos_ang_meanROI);
        
        %angle in radians rather than degrees (useful for Polar Plots)
        MeanROIanglerad = acos(cos_ang_meanROI);
        
        
        %% ==============================================
        % SAVE DATA
        %%===============================================
        
        %Data info
        Data_info = struct('NIH_ID', subj{sub},'CFM_ID',subj_cfm{sub},...
            'ROIsize',roiSize, 'ROIcoord', ROI_centre_coord_approx,...
            'cfmData', cfmData, 'MRIscan', ef_vol, 'surfaceData', surface);
        sl_info = 'sl = GM surface';
        
        % Grey matter
        GM_vectors = struct('EFvectors', ef_vec, 'EFnorm', ef_norm, 'snorm', snorm);
        GM_distance = struct('EFvectors', SE_Distance, 'EFnorm', SE_DistanceMAG);
        GM_angle = struct('Degrees', angle, 'Radians', anglerad);
        GM_efmag = ef_vec_mag;
        
        % ROI
        ROI_vertices = sl.vertices(ROI_vert_ind,:);
        ROI_Index = ROI_vert_ind;
        ROI_vectors = struct('EFvectors', efvec_inROI, 'EFnorm', ef_norm_inROI, 'snorm', snormROI);
        ROI_distance = struct('EFvectors', SE_DistanceROI, 'EFnorm', SE_DistanceROInorm);
        ROI_angle = struct('Degrees', ROIangle, 'Radians', ROIanglerad);
        ROI_efmag = ef_vec_mag_inROI;
        
        % Means w/in ROI
        meanROI_vertices = sl.vertices(ROI_centre_ind,:);
        meanROI_Index = ROI_centre_ind;
        meanROI_vectors = struct('EFvectors', meanEFVecinROI,'EFnorm', meanEFNORMinROI,...
            'snorm', meanSNORMinROI);
        meanROI_distance = struct('EFvectors',MeanSE_Distance_inROI, 'EFnorm', MeanSE_Distance_inROInorm);
        meanROI_angle = struct('Degrees', MeanROIangle, 'Radians', MeanROIanglerad);
        meanROI_efmag = meanEFMAGinROI;
        
        %Save relevant variables to file: surface data
        savefileName = sprintf('%s_surfacedata_%s.mat', subj{sub},ROIlabel{mont});
%        cd([ROIPATH, sprintf('%s/', subj{sub})]);
        cd(ROIPATH);
        save(savefileName, 'Data_info', 'sl', 'sl_info', 'GM_vectors', 'GM_distance', 'GM_angle', 'GM_efmag',...
            'ROI_vertices', 'ROI_Index', 'ROI_vectors', 'ROI_distance', 'ROI_angle', 'ROI_efmag',...
            'meanROI_vertices', 'meanROI_Index', 'meanROI_vectors', 'meanROI_distance', 'meanROI_angle', 'meanROI_efmag');
        %movefile(fileName,PATH);
        
        %clearvars -except subject sl sl_info GM_vertices GM_vectors GM_distance GM_angle ROI_vertices...
        %ROI_vectors ROI_distance ROI_angle meanROI_vertices meanROI_vectors meanROI_distance meanROI_angle
        
        fprintf('%s finished \n', subj_cfm{sub});
        clearvars -except surfPATH cfmPATH ROIPATH mont subj subj_cfm cfmResultFile ROIpatch ROIfile roiSize ROIlabel gmfilename        
        
    end
end
