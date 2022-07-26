%% Calculate ANGLE between electrode coordinates and motor strip / M1 bank ROI
% (c) Carys Evans, UCL
% carys.evans@ucl.ac.uk
% July 2022

%Calculates ANGLE between motor strip and electrode coordinates
%Calculates ANGLE between M1 bank ROI and electrode coodinates
%Saves relevant data to All_angle_MS_elecs.mat/.csv file


%% Define file names and directories

%set directories
MotorStripDataPATH = ('D:\PATH\TO\MOTORSTRIP\');            %motor strip data
datadir = 'D:\PATH\TO\ELECTRODE\DATA\';                     %electrode data

%set files
CoordFile = 'M1_motorstrip_coords.mat';          %motor strip coordinates
elecCoords = 'allPTs_electrodeCoords.mat';       %electrode coordinates 
montage = {'CP3FCZ','CPZFC3','C1FP2'};           %electrode montages
savefileName = 'All_angle_MS_elecs';             %output filename 
M1bankROIcoords = 'M1bankcoord.mat';             %M1bank coordinate file

%check if motor strip and M1 coordinate files exist, else skip subject
cd(MotorStripDataPATH)
k = dir('1*'); k={k.name}';         %subject folders

%define subjects
for i = 1:length(k)
    cd([MotorStripDataPATH, sprintf('%s/%s', k{i},k{i})]);
    
    if isfile(CoordFile) && isfile(M1bankROIcoords)
        subj(i,1) = k(i);
    else
        continue
    end
end

%load electrode coordinates
cd(datadir);
load(elecCoords);

%define subjects (necessary if filenames differ) 
subj_cfm = electrodeCoord_data.CP3FCZ.subjs; %(~cellfun('isempty',electrodeCoord_data.CP3FCZ.subjs));

%% CALCULATE ANGLE BETWEEN SNORM AND EF NORM
%%===============================================

for s = 1:length(subj)
    for m = 1:length(montage)
        
        if isempty(subj{s})
            continue
        else
        end
 
        cd([MotorStripDataPATH, sprintf('%s/%s', subj{s},subj{s})]);
        
        %load M1 ROI data
        surfdata = strcat(subj{s},'_surfacedata_M1_',montage{m},'.mat');
        surfdata = load(surfdata);
        ROIvect = surfdata.meanROI_vectors.snorm;
        
        %load motor strip coordinates
        load(CoordFile)

        %motor strip coordinates / vector
        stripCoords_xyz = [M1medial; M1lateral];
        stripVect = stripCoords_xyz(2,:)-stripCoords_xyz(1,:); %motor strip vector: lateral minus medial coords        
        stripVectNorm = flip([0,(-stripVect(3)),stripVect(2)]); %creates angle orthogonal to stripVect, pointing posterior lateral to anterior medial (in line with electrode direction anode-cathode)
        % stripVectNorm = [-stripVect(2),stripVect(1),0];%creates angle orthogonal to stripVect, pointing anterior medial to posterior lateral
        
        %get electrode coordinates / vector
        anCoord = electrodeCoord_data.(montage{m}).anodeCoordsxyz(s,:); %anode
        caCoord = electrodeCoord_data.(montage{m}).cathodeCoordsxyz(s,:); %cathode
        elecVect = caCoord-anCoord; %electrode vector: cathode minus anode coords (anterior minus posterior most electrode)
        
        
        %% Calculate ANGLE between Motor Strip and Electrode vectors: cos(?) = a · b / (|a| x |b|)
        
        % a · b = stripVect · elecVect (motor strip vectors & electrode vectors)
        normDotP = dot(stripVect,elecVect,2); %dot product
        
        %|a| = magnitude of motor strip vector
        stripV_mag = sqrt((stripVect(:,1).^2)+(stripVect(:,2).^2)+(stripVect(:,3).^2));
        
        %|b| = magnitude of electrode vector
        elecV_mag = sqrt((elecVect(:,1).^2)+(elecVect(:,2).^2)+(elecVect(:,3).^2));
        
        %cos(?).
        cos_ang = (normDotP./(stripV_mag.*elecV_mag));
        
        %? = angle between motor strip and electrode vectors. acosd = inverse cosine in degrees
        angle = acosd(cos_ang);
        
        %angle in radians rather than degrees (useful for Polar Plots)
        anglerad = acos(cos_ang);
        
        %% Calculate ANGLE between Motor Strip NORMAL and Electrode vectors: cos(?) = a · b / (|a| x |b|)
        
        % a · b = snorm · elecVect (motor strip vectors & electrode vectors)
        normDotPn = dot(stripVectNorm,elecVect,2); %dot product
        
        %|a| = magnitude of motor strip vector
        stripVn_mag = sqrt((stripVectNorm(:,1).^2)+(stripVectNorm(:,2).^2)+(stripVectNorm(:,3).^2));
        
        %|b| = magnitude of electrode vector
        %elecV_mag = sqrt((elecVect(:,1).^2)+(elecVect(:,2).^2)+(elecVect(:,3).^2));
        
        %cos(?).
        cos_angN = (normDotPn./(stripVn_mag.*elecV_mag));
        
        %? = angle between motor strip and electrode vectors. acosd = inverse cosine in degrees
        angleN = acosd(cos_angN);
        
        %angle in radians rather than degrees (useful for Polar Plots)
        angleradN = acos(cos_angN);
        
        
        %% Calculate ANGLE between M1 ROI and Electrode vectors: cos(?) = a · b / (|a| x |b|)
        
        % a · b = ROIvect · ef_vec (ROI norm vectors & ef vectors)
        normDotP = dot(ROIvect,elecVect,2); %dot product
        
        %|a| = magnitude of ROI vector
        ROIV_mag = sqrt((ROIvect(:,1).^2)+(ROIvect(:,2).^2)+(ROIvect(:,3).^2));
        
        %|b| = magnitude of electrode vector
        %elecV_mag = sqrt((elecVect(:,1).^2)+(elecVect(:,2).^2)+(elecVect(:,3).^2));
        
        %cos(?).
        cos_ang = (normDotP./(ROIV_mag.*elecV_mag));
        
        %? = angle between ROI and electrode vectors. acosd = inverse cosine in degrees
        angleROI = acosd(cos_ang);
        
        %angle in radians rather than degrees (useful for Polar Plots)
        angleradROI = acos(cos_ang);
        
        %% Save data
        
        strip_data = struct('strip_coords',stripCoords_xyz,...
            'strip_vectors',stripVect,'strip_vectorsNorm', stripVectNorm,...
            'stripV_magnitude', stripV_mag);
        elec_data = struct('anode_coords', anCoord, 'cathode_coords', caCoord,...
            'elec_vectors', elecVect, 'elecV_magnitude', elecV_mag);

        subjdata.strip_data = strip_data;
        subjdata.elec_data = elec_data;
        subjdata.M1_ROIvectors = ROIvect;
        
        subjdata.angle = angle;
        subjdata.angleN = angleN;
        subjdata.angleROI = angleROI;
       
        alldata.angle(s,m) = angle;
        alldata.angleN(s,m) = angleN;
        alldata.angleROI(s,m) = angleROI;
        alldata.angle_montage = montage;
        alldata.subj(s,1) = subj(s);
        alldata.subj(s,2) = subj_cfm(s);
               
        alldata.(subj_cfm{s}).(montage{m}) = subjdata;
        
        clearvars -except MotorStripDataPATH CoordFile datadir subj subj_cfm montage s m electrodeCoord_data alldata savefileName
        
    end
    
end

%% SAVE DATA
%save as .mat
cd(datadir)
matfile = [savefileName,'.mat'];
save(matfile,'alldata');

%save as .csv
% csvfile = [savefileName,'.csv'];
% mymont = [repmat({'PA'},(length(alldata.subj)),1);...
%     repmat({'ML'},(length(alldata.subj)),1);...
%     repmat({'C'},(length(alldata.subj)),1)];
% tabdat = table([alldata.angle(:,1);alldata.angle(:,2);alldata.angle(:,3)],...
% [alldata.angleROI(:,1);alldata.angleROI(:,2);alldata.angleROI(:,3)],...    
% mymont,(repmat(alldata.subj(:,1),3,1)),(repmat(alldata.subj(:,2),3,1)));
% tabdat.Properties.VariableNames = ...
%     {'Angle','AngleROI','Montage','SubjNum','SubjName'};
% writetable(tabdat,csvfile);

%%
% %% Optional: plot validation - Motor Strip
% %%===============================================
% %%
% % 
% cd('D:\DATA\LOCATION');
% sl = gifti('white-pial.gii');
%  mymontage = 'CPZFC3';
% s=10; %vector length
%  
% 
% myanglePA = nanmean(round(alldata.angleN(:,1)));
% myanglePA = strcat((num2str(myanglePA)),char(176));
% myangleML = nanmean(round(alldata.angleN(:,2)));
% myangleML = strcat((num2str(myangleML)),char(176));
% 
% h = plot(sl);set(h,'facecolor',[0.7 0.7 0.7],'facealpha',0.3); 
% %xlim([-100,100]),ylim([-100,100]),zlim([-100,100]);
% hold on;
% 
% for i = 1:50
%     if isempty(alldata.subj{i,2})
%         continue
%     else
%     end  
%     
%     s2 =     scatter3(alldata.(alldata.subj{i,2}).(mymontage).elec_data.cathode_coords(1),...
%         alldata.(alldata.subj{i,2}).(mymontage).elec_data.cathode_coords(2),...
%         alldata.(alldata.subj{i,2}).(mymontage).elec_data.cathode_coords(3),s,'b',...
%         'linewidth',1);
%     hold on;
% end

% %% Optional: Plot all coordinates on one head
% %%===============================================
%    
% for i = 1:50
%     if isempty(alldata.subj{i,2})
%         continue
%     else
%     end  
% 
% anode
% s1 =     scatter3(alldata.(alldata.subj{i,2}).(mymontage).elec_data.anode_coords(1),...
%         alldata.(alldata.subj{i,2}).(mymontage).elec_data.anode_coords(2),...
%         alldata.(alldata.subj{i,2}).(mymontage).elec_data.anode_coords(3),s,'r',...
%         'linewidth',1);
%     hold on;
%
% cathode
% s2 =     scatter3(alldata.(alldata.subj{i,2}).(mymontage).elec_data.cathode_coords(1),...
%         alldata.(alldata.subj{i,2}).(mymontage).elec_data.cathode_coords(2),...
%         alldata.(alldata.subj{i,2}).(mymontage).elec_data.cathode_coords(3),s,'b',...
%         'linewidth',1);
% 
% %motor strip endpoint 1
% s3 =    scatter3(alldata.(alldata.subj{i,2}).(mymontage).strip_data.strip_coords(1,1),...
%         alldata.(alldata.subj{i,2}).(mymontage).strip_data.strip_coords(1,2),...
%         alldata.(alldata.subj{i,2}).(mymontage).strip_data.strip_coords(1,3),s,...
%         [0.2 0.2 0.2],'linewidth',1);
%
% %motor strip endpoint 2
% s4 =    scatter3(alldata.(alldata.subj{i,2}).(mymontage).strip_data.strip_coords(2,1),...
%         alldata.(alldata.subj{i,2}).(mymontage).strip_data.strip_coords(2,2),...
%         alldata.(alldata.subj{i,2}).(mymontage).strip_data.strip_coords(2,3),s,...
%         [0.2 0.2 0.2],'linewidth',1);    
%  
% %motor strip
% q1 =     quiver3(alldata.(alldata.subj{i,2}).(mymontage).strip_data.strip_coords(1,1),...
%         alldata.(alldata.subj{i,2}).(mymontage).strip_data.strip_coords(1,2),...
%         alldata.(alldata.subj{i,2}).(mymontage).strip_data.strip_coords(1,3),...
%         alldata.(alldata.subj{i,2}).(mymontage).strip_data.strip_vectorsNorm(1,1),...
%         alldata.(alldata.subj{i,2}).(mymontage).strip_data.strip_vectorsNorm(1,2),...
%         alldata.(alldata.subj{i,2}).(mymontage).strip_data.strip_vectorsNorm(1,3),0,...
%         'linewidth',1, 'color', [0 0 0]);
% 
% %electrodes: CP3FCZ (black)
% q2 =     quiver3(alldata.(alldata.subj{i,2}).CP3FCZ.elec_data.anode_coords(1),...
%         alldata.(alldata.subj{i,2}).CP3FCZ.elec_data.anode_coords(2),...
%         alldata.(alldata.subj{i,2}).CP3FCZ.elec_data.anode_coords(3),...
%         alldata.(alldata.subj{i,2}).CP3FCZ.elec_data.elec_vectors(1),...
%         alldata.(alldata.subj{i,2}).CP3FCZ.elec_data.elec_vectors(2),...
%         alldata.(alldata.subj{i,2}).CP3FCZ.elec_data.elec_vectors(3),0,'k',...
%         'linewidth',1,'linestyle','--');
% 
% %electrodes: CP3FCZ
% q3 =     quiver3(alldata.(alldata.subj{i,2}).CP3FCZ.elec_data.anode_coords(1),...
%         alldata.(alldata.subj{i,2}).CP3FCZ.elec_data.anode_coords(2),...
%         alldata.(alldata.subj{i,2}).CP3FCZ.elec_data.anode_coords(3),...
%         alldata.(alldata.subj{i,2}).CP3FCZ.elec_data.elec_vectors(1),...
%         alldata.(alldata.subj{i,2}).CP3FCZ.elec_data.elec_vectors(2),...
%         alldata.(alldata.subj{i,2}).CP3FCZ.elec_data.elec_vectors(3),0,...
%         'linewidth',1,'linestyle','--','color',[0 0 0]);   
%  
% %electrodes: CPZFC3
% q4 =     quiver3(alldata.(alldata.subj{i,2}).CPZFC3.elec_data.anode_coords(1),...
%         alldata.(alldata.subj{i,2}).CPZFC3.elec_data.anode_coords(2),...
%         alldata.(alldata.subj{i,2}).CPZFC3.elec_data.anode_coords(3),...
%         alldata.(alldata.subj{i,2}).CPZFC3.elec_data.elec_vectors(1),...
%         alldata.(alldata.subj{i,2}).CPZFC3.elec_data.elec_vectors(2),...
%         alldata.(alldata.subj{i,2}).CPZFC3.elec_data.elec_vectors(3),0,...
%         'linewidth',1,'linestyle',':', 'color', [0 0 0]);               
% 
%     legend([s1 s2 s3 q1 q2 q3 q4],...
%         {'anode','cathode','motor strip endpoints','motor strip orientation',...
%         'electrode orientation', ['PA: ',myanglePA], ['ML: ',myangleML]},'location','southeastoutside'); legend boxoff
%     set(gcf,'color','white')
% end
