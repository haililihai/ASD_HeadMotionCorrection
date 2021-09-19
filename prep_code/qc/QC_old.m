% ------------------------------------------------------------------------------
% This script performs the primary analyses for
% Parkes, Fulcher, Yucel, Fornito. An evaluation of the efficacy, reliability, 
% and sensitivity of motion correction strategies for resting-state functional MRI
% 
% Copyright (C) 2017, Linden Parkes <lindenparkes@gmail.com>,
% ------------------------------------------------------------------------------
clear all; close all; clc

% ------------------------------------------------------------------------------
% Set string switches
% ------------------------------------------------------------------------------
% REQUIRED
Projects = {'TR2_200_01','TR2_200_02','TR2_240_01','TR2_240_02','TR3_120_01','TR3_120_02','TR3_240_01','TR3_240_02'};
WhichProject = Projects{8};

% Parcellation, 'Gordon' or 'Yeo'
% REQUIRED
WhichParc = 'Yeo';

% 'Motion' 'Diagnostic'
% Note, this only effect the NBS parts of the script.
WhichSplit = 'Motion';

% ------------------------------------------------------------------------------
% Set logical switches
% ------------------------------------------------------------------------------
% REQUIRED
runIndividual = false;
runPlot = true;


fprintf(1, 'Running rfMRI QC.\n\tDataset: %s\n\tParcelation: %s\n',WhichProject,WhichParc);

% ------------------------------------------------------------------------------
% Set project variables
% ------------------------------------------------------------------------------
% REQUIRED
parentdir = '/data/';
funcdir = [parentdir,'ASD/code/func/'];
addpath(funcdir);

switch WhichProject
    case 'TR2_200_01'
        projdir = parentdir;
        sublist = [projdir,'ASD/TR2_200/sub.csv'];
        datadir = [projdir,'ASD/TR2_200/rest/'];
        preprostr = '/01/prepro/';
        t1dir= [projdir,'ASD/TR2_200/anat/'];
        sessionid='01';
        plotdir = [projdir,'ASD/TR2_200/plotdir/01/'];
        TR = 2;
    case 'TR2_200_02'
        projdir = parentdir;
        sublist = [projdir,'ASD/TR2_200/sub.csv'];
        datadir = [projdir,'ASD/TR2_200/rest/'];
        preprostr = '/02/prepro/';
        t1dir= [projdir,'ASD/TR2_200/anat/'];
        sessionid='02';
        plotdir = [projdir,'ASD/TR2_200/plotdir/02/'];
        TR = 2;
    case 'TR2_240_01'
        projdir = parentdir;
        sublist = [projdir,'ASD/TR2_240/sub.csv'];
        datadir = [projdir,'ASD/TR2_240/rest/'];
        preprostr = '/01/prepro/';
        t1dir= [projdir,'ASD/TR2_240/anat/'];
        sessionid='01';
        plotdir = [projdir,'ASD/TR2_240/plotdir/01/'];
        TR = 2;
    case 'TR2_240_02'
        projdir = parentdir;
        sublist = [projdir,'ASD/TR2_240/sub.csv'];
        datadir = [projdir,'ASD/TR2_240/rest/'];
        preprostr = '/02/prepro/';
        t1dir= [projdir,'ASD/TR2_240/anat/'];
        sessionid='02';
        plotdir = [projdir,'ASD/TR2_240/plotdir/02/'];
        TR = 2;
    case 'TR3_120_01'
        projdir = parentdir;
        sublist = [projdir,'ASD/TR3_120/sub.csv'];
        datadir = [projdir,'ASD/TR3_120/rest/'];
        preprostr = '/01/prepro/';
        t1dir= [projdir,'ASD/TR3_120/anat/'];
        sessionid='01';
        plotdir = [projdir,'ASD/TR3_120/plotdir/01/'];
        TR = 3;
    case 'TR3_120_02'
        projdir = parentdir;
        sublist = [projdir,'ASD/TR3_120/sub.csv'];
        datadir = [projdir,'ASD/TR3_120/rest/'];
        preprostr = '/02/prepro/';
        t1dir= [projdir,'ASD/TR3_120/anat/'];
        sessionid='02';
        plotdir = [projdir,'ASD/TR3_120/plotdir/02/'];
        TR = 3;
    case 'TR3_240_01'
        projdir = parentdir;
        sublist = [projdir,'ASD/TR3_240/sub.csv'];
        datadir = [projdir,'ASD/TR3_240/rest/'];
        preprostr = '/01/prepro/';
        t1dir= [projdir,'ASD/TR3_240/anat/'];
        sessionid='01';
        plotdir = [projdir,'ASD/TR3_240/plotdir/01/'];
        TR = 3;
    case 'TR3_240_02'
        projdir = parentdir;
        sublist = [projdir,'ASD/TR3_240/sub.csv'];
        datadir = [projdir,'ASD/TR3_240/rest/'];
        preprostr = '/02/prepro/';
        t1dir= [projdir,'ASD/TR3_240/anat/'];
        sessionid='02';
        plotdir = [projdir,'ASD/TR3_240/plotdir/02/'];
        TR = 3;
end

% ------------------------------------------------------------------------------
% Set parcellation
% Note, code is not setup to process multiple parcellations concurrently.
% ------------------------------------------------------------------------------
% REQUIRED
ROIDir = [parentdir,'ASD/code/ROIs/'];
switch WhichParc
    case 'Gordon'
        Parc = 1;
        ROI_Coords = dlmread([ROIDir,'Gordon/Gordon_Centroids.txt']);
        fileName = [ROIDir,'Gordon/CommunityModified.txt'];
    case 'Yeo'
        Parc = 2;
        ROI_Coords = dlmread([ROIDir,'Yeo/Yeo_coordinates.txt']);
        fileName = [ROIDir,'Yeo/Yeo_community.txt'];
end

fileID = fopen(fileName);
ROIStruct = textscan(fileID,'%s'); ROIStruct = ROIStruct{1};
% rearrange by community
[ROIStruct_com,ROI_idx] = sort(ROIStruct);

% Convert text labels to unique integer values
ROILabels = unique(ROIStruct);
ROIStructID = zeros(size(ROIStruct));
numROIComms = length(ROILabels);
for i = 1:numROIComms
    x = find(strcmp(ROIStruct, ROILabels{i}));
    ROIStructID(x) = i;
end

% ------------------------------------------------------------------------------
% Load ROI coordinates
% ------------------------------------------------------------------------------
% Calculate pairwise euclidean distance
ROIDist = pdist2(ROI_Coords,ROI_Coords,'euclidean');

% Flatten distance matrix
ROIDistVec = LP_FlatMat(ROIDist);

% Calculate number of ROIs
numROIs = size(ROIDist,1);

% Calculate number of edges
numConnections = numROIs * (numROIs - 1) / 2;

% ------------------------------------------------------------------------------
% 							Preprocessing pipelines
% ------------------------------------------------------------------------------
switch WhichProject
    case 'TR2_200_01'
        noiseOptions = {'6P+2P',...
                        '6P+2P+GSR',...
                        'ICA-AROMA+2P',...
                        'ICA-AROMA+2P+GSR'};
        noiseOptionsNames = {'6HMP+2Phys',...
                             '6HMP+2Phys+GSR',...
                             'ICA-AROMA+2Phys',...
                             'ICA-AROMA+2Phys+GSR'};
    case 'TR2_200_02'
        noiseOptions = {'6P+2P',...
                        '6P+2P+GSR',...
                        'ICA-AROMA+2P',...
                        'ICA-AROMA+2P+GSR'};
        noiseOptionsNames = {'6HMP+2Phys',...
                             '6HMP+2Phys+GSR',...
                             'ICA-AROMA+2Phys',...
                             'ICA-AROMA+2Phys+GSR'};
    case 'TR2_240_01'
        noiseOptions = {'6P+2P',...
                        '6P+2P+GSR',...
                        'ICA-AROMA+2P',...
                        'ICA-AROMA+2P+GSR'};
        noiseOptionsNames = {'6HMP+2Phys',...
                             '6HMP+2Phys+GSR',...
                             'ICA-AROMA+2Phys',...
                             'ICA-AROMA+2Phys+GSR'};
    case 'TR2_240_02'
        noiseOptions = {'6P+2P',...
                        '6P+2P+GSR',...
                        'ICA-AROMA+2P',...
                        'ICA-AROMA+2P+GSR'};
        noiseOptionsNames = {'6HMP+2Phys',...
                             '6HMP+2Phys+GSR',...
                             'ICA-AROMA+2Phys',...
                             'ICA-AROMA+2Phys+GSR'};
    case 'TR3_120_01'
        noiseOptions = {'6P+2P',...
                        '6P+2P+GSR',...
                        'ICA-AROMA+2P',...
                        'ICA-AROMA+2P+GSR'};
        noiseOptionsNames = {'6HMP+2Phys',...
                             '6HMP+2Phys+GSR',...
                             'ICA-AROMA+2Phys',...
                             'ICA-AROMA+2Phys+GSR'};
    case 'TR3_120_02'
        noiseOptions = {'6P+2P',...
                        '6P+2P+GSR',...
                        'ICA-AROMA+2P',...
                        'ICA-AROMA+2P+GSR'};
        noiseOptionsNames = {'6HMP+2Phys',...
                             '6HMP+2Phys+GSR',...
                             'ICA-AROMA+2Phys',...
                             'ICA-AROMA+2Phys+GSR'};
    case 'TR3_240_01'
        noiseOptions = {'6P+2P',...
                        '6P+2P+GSR',...
                        'ICA-AROMA+2P',...
                        'ICA-AROMA+2P+GSR'};
        noiseOptionsNames = {'6HMP+2Phys',...
                             '6HMP+2Phys+GSR',...
                             'ICA-AROMA+2Phys',...
                             'ICA-AROMA+2Phys+GSR'};
    case 'TR3_240_02'
        noiseOptions = {'6P+2P',...
                        '6P+2P+GSR',...
                        'ICA-AROMA+2P',...
                        'ICA-AROMA+2P+GSR'};
        noiseOptionsNames = {'6HMP+2Phys',...
                             '6HMP+2Phys+GSR',...
                             'ICA-AROMA+2Phys',...
                             'ICA-AROMA+2Phys+GSR'};
end

numPrePro = length(noiseOptions);

% ------------------------------------------------------------------------------
% Subject list
% ------------------------------------------------------------------------------
fprintf(1, 'Loading metadata...\n');
metadata = readtable(sublist);
% convert participant IDs to strings
if ~iscellstr(metadata.ParticipantID)
    metadata.ParticipantID =  cellfun(@num2str, num2cell(metadata.ParticipantID), 'UniformOutput', false);
end

% Only do this if there is more than a single group
if numel(unique(metadata.Diagnosis)) > 1
    % Retain only group 1 (assumed to be HCs) and group 2 (assumed to be patients)
    % I do this because the OCDPG dataset has some PGs in it that need to be removed
    metadata(metadata.Diagnosis == 3,:) = [];
    
    % Retain only HCs
    metadata = metadata(metadata.Diagnosis == 1,:);
end

numGroups = numel(unique(metadata.Diagnosis));

% ------------------------------------------------------------------------------
% Exclusion and censoring stuff
% ------------------------------------------------------------------------------
[metadata.exclude,metadata.mov,metadata.fdJenk,metadata.fdJenk_m,metadata.fdPower,metadata.fdPower_m,metadata.dvars,metadata.JP12ScrubMask,metadata.JP14ScrubMask] = GetExcludeForSample(datadir,metadata.ParticipantID,TR,preprostr);
fprintf(1, 'done\n');

% compute number of volumes using the length of fdJenk
% note, this is assumed to be same for all subjects!
numVols = length(metadata.fdJenk{1});

% compute numsubs
numSubs = size(metadata,1);

% ------------------------------------------------------------------------------
% Plot Movement params
% ------------------------------------------------------------------------------
fd = metadata.fdJenk_m;
fprintf(1, 'Mean FD for sample w/ lenient: %f (%f) \n', round(mean(fd(~metadata.exclude(:,1))),2), round(std(fd(~metadata.exclude(:,1))),2))
fprintf(1, 'Mean FD for sample w/ stringent: %f (%f) \n', round(mean(fd(~metadata.exclude(:,2))),2), round(std(fd(~metadata.exclude(:,2))),2))
fprintf(1, 'Mean FD for sample w/ spike reg: %f (%f) \n', round(mean(fd(~metadata.exclude(:,3))),2), round(std(fd(~metadata.exclude(:,3))),2))
fprintf(1, 'Mean FD for sample w/ basic scrub: %f (%f) \n', round(mean(fd(~metadata.exclude(:,4))),2), round(std(fd(~metadata.exclude(:,4))),2))
fprintf(1, 'Mean FD for sample w/ optimized scrub: %f (%f) \n\n', round(mean(fd(~metadata.exclude(:,5))),2), round(std(fd(~metadata.exclude(:,5))),2))

if numGroups > 1
    fd_HC = fd(metadata.Diagnosis == 1);
    fprintf(1, 'Mean FD for HCs w/ lenient: %f (%f) \n', round(mean(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,1))),2), round(std(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,1))),2))
    fprintf(1, 'Mean FD for HCs w/ stringent: %f (%f) \n', round(mean(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,2))),2), round(std(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,2))),2))
    fprintf(1, 'Mean FD for HCs w/ spike reg: %f (%f) \n', round(mean(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,3))),2), round(std(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,3))),2))
    fprintf(1, 'Mean FD for HCs w/ basic scrub: %f (%f) \n', round(mean(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,4))),2), round(std(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,4))),2))
    fprintf(1, 'Mean FD for HCs w/ optimized scrub: %f (%f) \n\n', round(mean(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,5))),2), round(std(fd_HC(~metadata.exclude(metadata.Diagnosis == 1,5))),2))

    
    fd_Pat = fd(metadata.Diagnosis == 2);
    fprintf(1, 'Mean FD for patients w/ lenient: %f (%f) \n', round(mean(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,1))),2), round(std(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,1))),2))
    fprintf(1, 'Mean FD for patients w/ stringent: %f (%f) \n', round(mean(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,2))),2), round(std(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,2))),2))
    fprintf(1, 'Mean FD for patients w/ spike reg: %f (%f) \n', round(mean(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,3))),2), round(std(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,3))),2))
    fprintf(1, 'Mean FD for patients w/ basic scrub: %f (%f) \n', round(mean(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,4))),2), round(std(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,4))),2))
    fprintf(1, 'Mean FD for patients w/ optimized scrub: %f (%f) \n\n', round(mean(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,5))),2), round(std(fd_Pat(~metadata.exclude(metadata.Diagnosis == 2,5))),2))


    [h,p,~,stats] = ttest2(fd_HC,fd_Pat,'Vartype','unequal');
    if p < 0.05
        fprintf(1, 'There is a significant group difference in mean FD. t-value = %s. p-value = %s\n',num2str(stats.tstat),num2str(p));
    end

    BF_JitteredParallelScatter({fd_HC,fd_Pat});
    ax = gca;
    ax.XTick = [1,2];
    ax.XTickLabel = {'Controls','Patients'};
    xlabel('Group')
    ylabel('Mean FD')
end

% backup metadata
metadata_bak = metadata;

% plot dir
if exist(plotdir) == 0
    fprintf(1,'Initialising outdir\n')
    mkdir(plotdir)
elseif exist(plotdir) == 7
    fprintf(1,'Cleaning and re-initialising outdir\n')
    rmdir(plotdir,'s')
    mkdir(plotdir)
end

save(fullfile(plotdir,'metadata.mat'),'metadata','-v7.3');

% % ------------------------------------------------------------------------------
% % Variables
% % ------------------------------------------------------------------------------
% allData = struct('noiseOptions',noiseOptions,...
%                 'noiseOptionsNames',noiseOptionsNames,...
%                 'PercExcluded',[],...
%                 'cfg',struct,...
%                 'FC',[],...
%                 'FCVec',[],...
%                 'VarCovar',[],...
%                 'Var',[],...
%                 'GCOR',[],...
%                 'NaNFilter',[],...
%                 'QCFCVec',[],...
%                 'QCFC_PropSig_corr',[],...
%                 'QCFC_PropSig_unc',[],...
%                 'QCFC_AbsMed',[],...
%                 'QCFC_DistDep',[],...
%                 'QCFC_DistDep_Pval',[],...
%                 'MeanEdgeWeight',[],...
%                 'tDOF',[],...
%                 'tDOF_mean',[],...
%                 'tDOF_std',[],...
%                 'tDOF_gmean',[],...
%                 'tDOF_gstd',[]);

% % ------------------------------------------------------------------------------
% % 						Loop over preprocessing pipelines
% % ------------------------------------------------------------------------------
% for i = 1:numPrePro
%     allData(i).metadata=metadata;
%     WhichNoise = allData(i).noiseOptions;
%     WhichNoiseName = allData(i).noiseOptionsNames;
%     fprintf(1, '\nProcessing data: %s\n',WhichNoise);

%     WhichNoiseSplit = strsplit(WhichNoise,'+');

%     % ------------------------------------------------------------------------------
%     % Get time series and functional connectivity data
%     % ------------------------------------------------------------------------------
%     cfgFile = 'cfg.mat';
%     % WhichExclude = 5;
%     if any(strmatch('JP12Scrub',WhichNoiseSplit,'exact')) == 1
%         WhichNoise_temp = WhichNoiseSplit;
%         WhichNoise_temp(end) = [];
%         WhichNoise_temp = strjoin(WhichNoise_temp,'+');

%         WhichExclude = 4;
%         metadata = metadata_bak(~metadata_bak.exclude(:,WhichExclude),:);
%         numSubs = size(metadata,1);
%         fprintf(1, 'Excluded %u subjects \n', sum(metadata_bak.exclude(:,WhichExclude)));
%         [allData(i).cfg,allData(i).FC,allData(i).FCVec,allData(i).VarCovar,allData(i).Var,allData(i).GCOR] = GetFCForSample(datadir,metadata.ParticipantID,preprostr,WhichNoise_temp,cfgFile,Parc,numROIs,numConnections,metadata.JP12ScrubMask);
%     elseif any(strmatch('JP14Scrub',WhichNoiseSplit,'exact')) == 1
%         WhichExclude = 5;
%         metadata = metadata_bak(~metadata_bak.exclude(:,WhichExclude),:);
%         numSubs = size(metadata,1);
%         fprintf(1, 'Excluded %u subjects \n', sum(metadata_bak.exclude(:,WhichExclude)));
%         [allData(i).cfg,allData(i).FC,allData(i).FCVec,allData(i).VarCovar,allData(i).Var,allData(i).GCOR] = GetFCForSample(datadir,metadata.ParticipantID,preprostr,WhichNoise,cfgFile,Parc,numROIs,numConnections,metadata.JP14ScrubMask);
%     elseif any(strmatch('SpikeReg',WhichNoiseSplit,'exact')) == 1
%         WhichExclude = 3;
%         metadata = metadata_bak(~metadata_bak.exclude(:,WhichExclude),:);
%         numSubs = size(metadata,1);
%         fprintf(1, 'Excluded %u subjects \n', sum(metadata_bak.exclude(:,WhichExclude)));
%         [allData(i).cfg,allData(i).FC,allData(i).FCVec,allData(i).VarCovar,allData(i).Var,allData(i).GCOR] = GetFCForSample(datadir,metadata.ParticipantID,preprostr,WhichNoise,cfgFile,Parc,numROIs,numConnections);
%     else
%         WhichExclude = 1;
%         metadata = metadata_bak(~metadata_bak.exclude(:,WhichExclude),:);
%         numSubs = size(metadata,1);
%         fprintf(1, 'Excluded %u subjects \n', sum(metadata_bak.exclude(:,WhichExclude)));
%         [allData(i).cfg,allData(i).FC,allData(i).FCVec,allData(i).VarCovar,allData(i).Var,allData(i).GCOR] = GetFCForSample(datadir,metadata.ParticipantID,preprostr,WhichNoise,cfgFile,Parc,numROIs,numConnections);
%     end

%     % Store percentage of participants excluded
%     allData(i).PercExcluded = sum(metadata_bak.exclude(:,WhichExclude) / size(metadata_bak,1)) * 100;

%     % ------------------------------------------------------------------------------
%     % Compute QC-FC
%     % Note, make sure use z transformed r values
%     % ------------------------------------------------------------------------------
%     fprintf(1, 'Computing QC-FC: %s\n',WhichNoise);
%     [allData(i).QCFCVec,allData(i).NaNFilter,allData(i).QCFC_PropSig_corr,allData(i).QCFC_PropSig_unc,allData(i).QCFC_AbsMed,allData(i).QCFC_DistDep,allData(i).QCFC_DistDep_Pval] = RunQCFC(metadata.fdJenk_m,allData(i).FC,ROIDistVec);

%     % ------------------------------------------------------------------------------
%     % Compute mean edge weight
%     % ------------------------------------------------------------------------------
%     fprintf(1, 'Computing mean edge weights: %s\n',WhichNoise);
%     FCVec = [];
%     for j = 1:numSubs
%         % flatten FC for subject j
%         vec = LP_FlatMat(allData(i).FC(:,:,j));
%         % filter NaNs from QCFC analyses
%         vec = vec(allData(i).NaNFilter);
%         % store
%         FCVec(:,j) = vec;
%     end

%     % average across subject
%     allData(i).MeanEdgeWeight = mean(FCVec,2);

%     % ------------------------------------------------------------------------------
%     % Get tDOF
%     % ------------------------------------------------------------------------------
%     fprintf(1, 'Computing tDOF: %s\n',WhichNoise);
%     allData(i).tDOF = zeros(numSubs,1);
%     for j = 1:numSubs

%         % get tDOF
%         % First, find size of second dimension of noiseTS
%         allData(i).tDOF(j) = size(allData(i).cfg(j).noiseTS,2);

%         if any(strmatch('JP12Scrub',WhichNoiseSplit,'exact')) == 1
%             allData(i).tDOF(j) = allData(i).tDOF(j) + sum(metadata.JP12ScrubMask{j});
%         elseif any(strmatch('JP14Scrub',WhichNoiseSplit,'exact')) == 1
%             allData(i).tDOF(j) = allData(i).tDOF(j) + sum(metadata.JP14ScrubMask{j});
%         end

%         % Then, if ICA-AROMA pipeline, find number of ICs and add to tDOF
%         if ~isempty(strfind(WhichNoise,'ICA-AROMA'))
%             x = dlmread([datadir,metadata.ParticipantID{j},preprostr,'ICA-AROMA_output/classified_motion_ICs.txt']);
%             allData(i).tDOF(j) = allData(i).tDOF(j) + length(x);
%         end
%     end

%     % ------------------------------------------------------------------------------
%     % Calculate mean temporal degrees of freedom lost
%     % ------------------------------------------------------------------------------
%     % tDOF will be the same for most pipelines
%     % but some have variable regressor amounts
%     % so we take mean over subjects
%     allData(i).tDOF_mean = mean(allData(i).tDOF);
%     allData(i).tDOF_std = std(allData(i).tDOF);

%     % also get mean by diagnostic groups (metadata.Diagnosis)
%     if ismember('OCDPG',WhichProject,'rows') | ismember('UCLA',WhichProject,'rows') | ismember('COBRE',WhichProject,'rows')
%         allData(i).tDOF_gmean(1) = mean(allData(i).tDOF(metadata.Diagnosis == 1));
%         allData(i).tDOF_gmean(2) = mean(allData(i).tDOF(metadata.Diagnosis == 2));
        
%         allData(i).tDOF_gstd(1) = std(allData(i).tDOF(metadata.Diagnosis == 1));
%         allData(i).tDOF_gstd(2) = std(allData(i).tDOF(metadata.Diagnosis == 2));
%     end

%     % ------------------------------------------------------------------------------
%     % Perform t-test on tDOF-loss
%     % ------------------------------------------------------------------------------
%     if ismember('OCDPG',WhichProject,'rows') | ismember('UCLA',WhichProject,'rows') | ismember('COBRE',WhichProject,'rows')
%         x = allData(i).tDOF(metadata.Diagnosis == 1);
%         y = allData(i).tDOF(metadata.Diagnosis == 2);

%         [h,p,~,stats] = ttest2(x,y,'Vartype','unequal');
%         if p < 0.05
%             fprintf(1, 'Significant group difference in tDOF-loss. t-value(%s) = %s. p-value = %s\n',num2str(stats.df),num2str(stats.tstat),num2str(p));
%         else
%             fprintf(1, 'NO significant group difference in tDOF-loss. t-value(%s) = %s. p-value = %s\n',num2str(stats.df),num2str(stats.tstat),num2str(p));
%         end
%         fprintf(1, '\tMean tDOF-loss, group 1: %s\n', num2str(round(allData(i).tDOF_gmean(1),2)));
%         fprintf(1, '\tMean tDOF-loss, group 2: %s\n', num2str(round(allData(i).tDOF_gmean(2),2)));
%     end

%     % ------------------------------------------------------------------------------
%     % Get num components for aCompCor50
%     % ------------------------------------------------------------------------------
%     if any(strmatch('12P+aCC50',WhichNoise,'exact')) == 1
%         for j = 1:numSubs
%             aCC50_num_wm(j) = dlmread([datadir,metadata.ParticipantID{j},preprostr,'12P+aCC50/aCC_num_wm.txt']);
%             aCC50_num_csf(j) = dlmread([datadir,metadata.ParticipantID{j},preprostr,'12P+aCC50/aCC_num_csf.txt']);
%         end
%         aCC50_num_wm_mean = mean(aCC50_num_wm);
%         aCC50_num_wm_std = std(aCC50_num_wm);

%         aCC50_num_csf_mean = mean(aCC50_num_csf);
%         aCC50_num_csf_std = std(aCC50_num_csf);
%     end

%     if runIndividual
%         [srt,idx] = sort(metadata.fdJenk_m,'descend');
%         % ------------------------------------------------------------------------------
%         % Loop over subjects in order of descending movement issues
%         % ------------------------------------------------------------------------------
%         for i = 1:numSubs
%         % for i = [1,numSubs] % this will just do the highest and lowest motion subject
%             fprintf(1,'\tProcessing subject %u/%u: %s\n',i,numSubs,metadata.ParticipantID{idx(i)})

%             preprodir = [datadir,metadata.ParticipantID{idx(i)},preprostr,WhichNoise,'/'];
%             t1datadir = [t1dir,metadata.ParticipantID{idx(i)},'/',sessionid,'/'];

%             EPI = 'epi_prepro.nii.gz';
%             cfg = load([datadir,metadata.ParticipantID{idx(i)},preprostr,WhichNoise,'/',cfgFile]);
%             t1name = cfg.cfg.t1name;
%             gmMask = ['gm50_bin.nii.gz'];
%             wmMask = ['wbc2N4c',regexprep(t1name,'(\.nii)*(\.gz)?',''),'_e1.nii.gz'];

%             filesToCheck = {{EPI}, ...
%                             {gmMask,wmMask}};

%             dirsOfFiles = {preprodir, ...
%                             t1datadir};

%             % Check if the epi_prepro.nii exists at all, if it doesnt, then the subject was processed
%             % (this will happen if they were excluded according to optimized scrubbing procedures)
%             if exist([dirsOfFiles{1},filesToCheck{1}{1},'.gz']) == 2 | exist([dirsOfFiles{1},filesToCheck{1}{1}]) == 2

%                 fprintf(1, '\t\tDrawing figure...\n');
%                 % ------------------------------------------------------------------------------
%                 % GetTSCompartment
%                 % ------------------------------------------------------------------------------
%                 fsldir = '/usr/local/fsl/bin/'; % M3
%                 % fsldir = '/usr/share/fsl/5.0/bin/'; % local macbook
%                 [ts_compartment,key_compartment] = GetTSCompartment(fsldir,[preprodir,EPI],[t1datadir,gmMask],[t1datadir,wmMask]);
%                 % normalise
%                 ts_compartment = BF_NormalizeMatrix(ts_compartment,'maxmin');

%                 % ------------------------------------------------------------------------------
%                 % ThePlot
%                 % ------------------------------------------------------------------------------
%                 % Threshold for flagging problem volumes
%                 ThePlot([metadata.ParticipantID{idx(i)},' / ',WhichNoiseName],metadata.mov{idx(i)},metadata.fdPower{idx(i)},metadata.JP14ScrubMask{idx(i)},metadata.fdJenk{idx(i)},metadata.dvars{idx(i)},ts_compartment,key_compartment,TR)
                
%                 fig = gcf;
%                 set(fig,'PaperPositionMode','Auto')
%                 print(fig,[plotdir,num2str(i),'_',metadata.ParticipantID{idx(i)},'_',WhichNoise,'.bmp'],'-dbmp')
%                 close all

%             else
%                 fprintf(1, '\t\tData not present for %s / %s. Skipping...\n', metadata.ParticipantID{idx(i)}, WhichNoise);
%             end
%         end
%     end
% end

% % ------------------------------------------------------------------------------
% % Check whether NaNs filtered out of each pipeline are the same.
% % They should be since if there are NaNs in QCFC, the same subject(s) will cause them each time.
% % But, still worth checking explicitly
% % ------------------------------------------------------------------------------
% temp = [allData(:).NaNFilter]';
% allRowsEqual = size(unique(temp,'rows'),1) == 1;
% if allRowsEqual == 1
%     % If all nan filters are the same across pipelines
%     % Safe to use any of the nan filters to filter ROI dist vec permenantly
%     % ROIDistVec = ROIDistVec(allData(end).NaNFilter);
% elseif allRowsEqual ~= 1
%     warning('NaN filters are not the same across pipelines!');
%     % ROIDistVec = ROIDistVec(allData(end).NaNFilter);
% end

% % ------------------------------------------------------------------------------
% % Figures
% % ------------------------------------------------------------------------------
% FSize = 10;

% clear extraParams
% % ------------------------------------------------------------------------------
% % Chart colors and line styles
% % ------------------------------------------------------------------------------
% % tempColors = num2cell([255,105,97;97,168,255;178,223,138;117,112,179]./255,2);  
% tempColors = num2cell([255,105,97;97,168,255;178,223,138;117,112,179;255,179,71]./255,2);  
% for i = 1:numPrePro
%     strs = strsplit(allData(i).noiseOptions,'+');
%     % colors
%     if any(strmatch('6P',strs,'exact')) == 1; theColors{i} = tempColors{1}; end
%     if any(strmatch('24P',strs,'exact')) == 1; theColors{i} = tempColors{2}; end
%     if any(strmatch('aCC',strs,'exact')) == 1 || any(strmatch('aCC50',strs,'exact')) == 1; theColors{i} = tempColors{3}; end
%     if any(strmatch('ICA-AROMA',strs,'exact')) == 1; theColors{i} = tempColors{4}; end
    
%     if any(strmatch('SpikeReg',strs,'exact')) == 1 | ...
%         any(strmatch('JP12Scrub',strs,'exact')) == 1 | ...
%         any(strmatch('JP14Scrub',strs,'exact')) == 1
%         theColors{i} = tempColors{5};
%     end
    
%     % line style
%     if any(strmatch('GSR',strs,'exact')) == 1 | any(strmatch('4GSR',strs,'exact')) == 1 | any(strmatch('2GSR',strs,'exact')) == 1
%         theLines{i} = ':';
%         % theLines{i} = '--';
%     else
%         theLines{i} = '-';
%     end
% end

% % Plot labels
% x = {allData(:).noiseOptionsNames};
% y = num2cell(100 - [allData.PercExcluded]);
% xy = cell(x);
% for i = 1:length(x)
%     if y{i} ~= 100
%         xy{i} = strcat(x{i}, ' (',num2str(y{i},'%0.0f'),'%)');
%     elseif y{i} == 100
%         xy{i} = x{i};
%     end
% end

% % ------------------------------------------------------------------------------
% % Essential figures
% % ------------------------------------------------------------------------------
% if runPlot
%     % ------------------------------------------------------------------------------
%     % QC-FC significant proportion
%     % Corresponds to Figure 1 and Figure 4 in case of censoring
%     % ------------------------------------------------------------------------------
%     if ~exist(plotdir) mkdir(plotdir);end
    
%     Fig_QCFC_Dist = figure('color','w', 'units', 'centimeters', 'pos', [0 0 16 9], 'name',['Fig_QCFC_Dist']); box('on'); movegui(Fig_QCFC_Dist,'center');
%     sp = subplot(1,3,1);
%     pos = get(sp,'Position');
%     % set(gca,'Position',[pos(1)*2.75, pos(2)*1.2, pos(3)*1.4, pos(4)*1]); % [left bottom width height]
%     set(gca,'Position',[pos(1)*3.25, pos(2)*1.2, pos(3)*1.2, pos(4)*1]); % [left bottom width height]

%     % Create data
%     % data = {[allData(:).QCFC_PropSig_corr]'};
%     data = {[allData(:).QCFC_PropSig_unc]'};
%     data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));

%     % Create table
%     T = table(data{1},'RowNames',{allData(:).noiseOptionsNames}','VariableNames',{'QCFC_PropSig'})

%     % Create bar chart
%     clear extraParams
%     extraParams.xTickLabels = xy;
%     extraParams.xLabel = ''; % 'Pipeline'
%     % extraParams.yLabel = 'QC-FC (%)';
%     extraParams.yLabel = 'QC-FC uncorrected (%)';
%     extraParams.theColors = theColors;
%     extraParams.theLines = theLines;
%     extraParams.yLimits = [0 110];

%     TheBarChart(data,data_std,false,extraParams)

%     % ------------------------------------------------------------------------------
%     % QCFC distributions
%     % ------------------------------------------------------------------------------
%     sp = subplot(1,3,3);
%     pos = get(sp,'Position');
%     set(gca,'Position',[pos(1)*1, pos(2)*1.2, pos(3)*1.4, pos(4)*1]); % [left bottom width height]
%     data2 = {allData(:).QCFCVec};
%     clear extraParams
%     % extraParams.theLabels = {allData(:).noiseOptionsNames};
%     extraParams.customSpot = '';
%     extraParams.add0Line = true;
%     extraParams.theColors = theColors;
%     BF_JitteredParallelScatter(data2,1,1,0,extraParams);
%     ax = gca;

%     % Set axis stuff
%     ax.FontSize = FSize;
%     % ax.XTick = [1:size(data{1},1)];
%     ax.XTick = [];
%     ax.XTickLabel = [];
%     ax.XLim = ([0 numPrePro+1]);
%     if ismember('NYU_2',WhichProject,'rows') | ismember('OCDPG',WhichProject,'rows')
%         ax.YLim = ([-0.8 1.25]);
%     else
%         ax.YLim = ([-0.6 1]);
%     end
%     ylabel('QC-FC (Pearson''s r)')

%     % add text
%     TextRotation = 0;
%     strprec = '%0.2f';
%     data3 = {allData(:).QCFC_AbsMed};

%     text(1:size(data3,2),repmat(ax.YLim(2) - ax.YLim(2)*.05,1,size(data3,2)),num2str([data3{1,:}]',strprec),... 
%     'HorizontalAlignment','right',... 
%     'VerticalAlignment','middle',...
%     'Color','black',...
%     'FontSize', FSize,...
%     'Rotation',TextRotation)

%     view(90,90)

%     saveas(Fig_QCFC_Dist,fullfile(plotdir,'QCFC_Dist.png'));

%     % ------------------------------------------------------------------------------
%     % QC-FC corrected
%     % ------------------------------------------------------------------------------
%     % Create data
%     data = {[allData(:).QCFC_PropSig_corr]'};
%     data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));

%     % Create table
%     T = table(data{1},'RowNames',{allData(:).noiseOptionsNames}','VariableNames',{'QCFC_PropSig_corr'})

%     % Create bar chart
%     clear extraParams
%     % extraParams.xTickLabels = {allData(:).noiseOptionsNames};
%     extraParams.xTickLabels = xy;
%     extraParams.xLabel = ''; % 'Pipeline'
%     extraParams.yLabel = 'QC-FC FDR-corr. (%)';
%     extraParams.theColors = theColors;
%     extraParams.theLines = theLines;
%     extraParams.yLimits = [0 110];

%     Fig_QCFC_FDR = figure('color','w', 'units', 'centimeters', 'pos', [0 0 10.5 9], 'name',['Fig_QCFC_FDR']); box('on'); movegui(Fig_QCFC_FDR,'center');
%     sp = subplot(1,1,1);
%     pos = get(sp,'Position');
%     set(gca,'Position',[pos(1)*4.9, pos(2)*1.2, pos(3)*0.425, pos(4)*1]); % [left bottom width height]

%     TheBarChart(data,data_std,false,extraParams)

%     saveas(Fig_QCFC_FDR,fullfile(plotdir,'QCFC_FDR.bmp'));

%     % ------------------------------------------------------------------------------
%     % QC-FC distance dependence
%     % Corresponds to Figure 2 and Figure 4 in case of censoring
%     % ------------------------------------------------------------------------------
%     data = {[allData(:).QCFC_DistDep]'};

%     data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));
    
%     % Create table
%     T = table(data{1},'RowNames',{allData(:).noiseOptionsNames}','VariableNames',{'QCFC_DistDep'})

%     % Create bar chart
%     extraParams.xTickLabels = xy;
%     extraParams.xLabel = '';
%     extraParams.yLabel = 'QC-FC distance dependence (Spearman''s rho)';
%     extraParams.theColors = theColors;
%     extraParams.theLines = theLines;
%     extraParams.yLimits = [-0.5 0.5];

%     Fig_QCFC_DistDep = figure('color','w', 'units', 'centimeters', 'pos', [0 0 10.5 9], 'name',['Fig_QCFC_DistDep']); box('on'); movegui(Fig_QCFC_DistDep,'center');
%     sp = subplot(1,1,1);
%     pos = get(sp,'Position');
%     set(gca,'Position',[pos(1)*4.9, pos(2)*1.2, pos(3)*0.425, pos(4)*1]); % [left bottom width height]

%     TheBarChart(data,data_std,false,extraParams)

%     saveas(Fig_QCFC_DistDep,fullfile(plotdir,'QCFC_DistDep.bmp'));

%     % ------------------------------------------------------------------------------
%     % Figures: QCFC
%     % ------------------------------------------------------------------------------
%     % Initialise figures
%     Fig_QCFC_DistDepBig = figure('color','w', 'units', 'centimeters', 'pos', [0 0 25 27], 'name',['Fig_QCFC_DistDepBig']); box('on'); movegui(Fig_QCFC_DistDepBig,'center');

%     pipelines2Retain = logical([1 1 1 1]);
%     noiseOptions_temp = {allData(pipelines2Retain).noiseOptions};
%     noiseOptionsNames_temp = {allData(pipelines2Retain).noiseOptionsNames};
%     theColors_temp = theColors(pipelines2Retain);

%     numPrePro_temp = length(noiseOptions_temp);

%     for i = 1:numPrePro_temp

%         WhichNoise = noiseOptions_temp{i};
%         WhichNoiseName = noiseOptionsNames_temp{i};
%         idx = strmatch(WhichNoise,{allData(:).noiseOptions},'exact');

%         % ------------------------------------------------------------------------------
%         % Plot: distance dependence
%         % ------------------------------------------------------------------------------
%         figure(Fig_QCFC_DistDepBig)
%         subplot(4,ceil(numPrePro_temp/4),i)
%         set(gca,'FontSize',FSize)

%         % Bin QCFC data by distance and generate means and stds for each
%         numThresholds = numConnections; % bins
%         BF_PlotQuantiles(ROIDistVec(allData(idx).NaNFilter),allData(idx).QCFCVec,numThresholds,0,0,theColors_temp{i})
%         hold on
%         plot([0:200],zeros(1,201),'--','Color','k')
        
%         xlabel('Distance (mm)')
%         ylabel('QC-FC')

%         title(WhichNoiseName,'Interpreter', 'none','FontSize',FSize,'FontWeight','normal')
%     end

%     saveas(Fig_QCFC_DistDepBig,fullfile(plotdir,'QCFC_DistDepBig.bmp'));

%     % ------------------------------------------------------------------------------
%     % Figures: QCFC bins
%     % ------------------------------------------------------------------------------
%     % Initialise figures
%     Fig_QCFC_DistDepBig_bins = figure('color','w', 'units', 'centimeters', 'pos', [0 0 25 27], 'name',['Fig_QCFC_DistDepBig_bins']); box('on'); movegui(Fig_QCFC_DistDepBig_bins,'center');

%     pipelines2Retain = logical([1 1 1 1]);
%     noiseOptions_temp = {allData(pipelines2Retain).noiseOptions};
%     noiseOptionsNames_temp = {allData(pipelines2Retain).noiseOptionsNames};
%     theColors_temp = theColors(pipelines2Retain);

%     numPrePro_temp = length(noiseOptions_temp);

%     for i = 1:numPrePro_temp

%         WhichNoise = noiseOptions_temp{i};
%         WhichNoiseName = noiseOptionsNames_temp{i};
%         idx = strmatch(WhichNoise,{allData(:).noiseOptions},'exact');

%         % ------------------------------------------------------------------------------
%         % Plot: distance dependence
%         % ------------------------------------------------------------------------------
%         figure(Fig_QCFC_DistDepBig_bins)
%         subplot(4,ceil(numPrePro_temp/4),i)
%         set(gca,'FontSize',FSize)

%         % Bin QCFC data by distance and generate means and stds for each
%         numThresholds = 20; % bins
%         BF_PlotQuantiles(ROIDistVec(allData(idx).NaNFilter),allData(idx).QCFCVec,numThresholds,0,0,theColors_temp{i})
%         hold on
%         plot([0:200],zeros(1,201),'--','Color','k')
        
%         xlabel('Distance (mm)')
%         ylabel('QC-FC')

%         title(WhichNoiseName,'Interpreter', 'none','FontSize',FSize,'FontWeight','normal')
%     end

%     saveas(Fig_QCFC_DistDepBig_bins,fullfile(plotdir,'QCFC_DistDepBig_bins.bmp'));


%     % ------------------------------------------------------------------------------
%     % tDOF
%     % Corresponds to Figure 5
%     % ------------------------------------------------------------------------------
%     % Create data
%     data = {[allData(:).tDOF_mean]'};
%     data_std = {[allData(:).tDOF_std]'};

%     % Create table
%     T = table(data{1},'RowNames',{allData(:).noiseOptionsNames}','VariableNames',{'tDOF_mean'})

%     % Create bar chart
%     extraParams.xTickLabels = xy;
%     extraParams.xLabel = '';
%     extraParams.yLabel = 'tDOF-loss (# regressors)';
%     extraParams.theColors = theColors;
%     extraParams.theLines = theLines;
%     extraParams.yLimits = [0 140];

%     Fig_tDOF = figure('color','w', 'units', 'centimeters', 'pos', [0 0 10.5 9], 'name',['Fig_tDOF']); box('on'); movegui(Fig_tDOF,'center');
%     sp = subplot(1,1,1);
%     pos = get(sp,'Position');
%     set(gca,'Position',[pos(1)*4.9, pos(2)*1.2, pos(3)*0.425, pos(4)*1]); % [left bottom width height]

%     TheBarChart(data,data_std,false,extraParams)

%     saveas(Fig_tDOF,fullfile(plotdir,'tDOF.bmp'));

% end

% % ------------------------------------------------------------------------------
% % Overlap: sig edges across pipelines
% % ------------------------------------------------------------------------------
% runOverlapPlots = false;
% if runOverlapPlots
%     % ------------------------------------------------------------------------------
%     % 1) Pairwise overlap between significant networks
%     % ------------------------------------------------------------------------------
%     WhichOverlap = 'jaccard'; % 'jaccard' 'phi' 
%     f = figure('color','w', 'units', 'centimeters', 'pos', [0 0 23 12], 'name',['NBS Overlap']); box('on'); movegui(f,'center');
%     for i = 1:numContrasts
%         if i == 1
%             data = cellfun(@(x) x(2),{allData(:).NBS_sigVec}); data = cell2mat(data);
%         elseif i == 2
%             data = cellfun(@(x) x(1),{allData(:).NBS_sigVec}); data = cell2mat(data);
%         end

%         % Retain only pipelines showing an effect
%         pipeFilter = any(data);
%         data = data(:,pipeFilter);
%         data = full(data);

%         % Filter NaNs
%         NaNFilter = [allData(1).NaNFilter];
%         data = data(NaNFilter,:); 

%         switch WhichOverlap
%             case 'jaccard'
%                 mat = [];
%                 for j = 1:size(data,2)
%                     for k = 1:size(data,2)
%                         x = data(:,j);
%                         y = data(:,k);
%                         % Num of signifcant edges in intersection
%                         Ci = sum(x + y == 2);
%                         % Num of signifcant edges in union
%                         Cu = sum(x + y > 0);

%                         mat(j,k) = Ci/Cu;
%                     end
%                 end
%             case 'phi'
%                 mat = corr(data,'type','Pearson');
%         end

%         % ------------------------------------------------------------------------------
%         % Plot
%         % ------------------------------------------------------------------------------
%         subplot(1,2,i)
%         imagesc(mat)
%         axis square
%         axis tight
%         % colormap([flipud(BF_getcmap('blues',9,0));1,1,1;BF_getcmap('reds',9,0)])
%         colormap(BF_getcmap('reds',9,0))
%         caxis([0 1])
%         colorbar
%         ax = gca;
%         ax.FontSize = FSize;
%         x = [1:length(mat)];
%         ax.XTick = x; ax.YTick = x;
%         % ax.XTickLabel = '';
%         ax.XTickLabel = xy(:,pipeFilter);
%         ax.XTickLabelRotation = 45;
%         ax.YTickLabel = xy(:,pipeFilter);
%         ax.TickLength = [0,0];

%         % plot values in lower triangle of matrix
%         for i = 1:size(mat,1)
%             for j = i:size(mat,2)
%                 if i ~= j
%                     text(i,j,num2str(mat(i,j),'%0.1f'),'HorizontalAlignment','center',...
%                         'Color','k','FontSize',FSize,'FontWeight','normal');
%                 end
%             end
%         end

%         if i == 1
%             title('A) SCZ>HC','FontSize',FSize,'FontWeight','normal')
%         elseif i == 2
%             title('B) HC>SCZ','FontSize',FSize,'FontWeight','normal')
%         end
%     end

%     % ------------------------------------------------------------------------------
%     % Figure 10: neuromarvl
%     % Neuromarvl colours
%     % BLUE | GREEN | YELLOW | RED
%     % #0029ff | #00ff75 | #ffff00 | #ff0000
%     % ------------------------------------------------------------------------------
%     pipelines2Retain = logical([0 0 0 0 0 1 0 1 0 0 0 0 0 1 0 0 1 0 0]);
%     noiseOptions_temp = {allData(pipelines2Retain).noiseOptions};

%     % node coords
%     T = table(ROI_Coords(:,1),ROI_Coords(:,2),ROI_Coords(:,3),'VariableNames',{'x','y','z'});
%     writetable(T,'coordinates.txt')
%     clear T

%     % node attributes
%     T = table(ones(numROIs,1),ROIStructID,'VariableNames',{'dummy','ROIStructID'});
%     writetable(T,'attributes.txt')

%     % Neuromarvl files
%     for i = 1:length(noiseOptions_temp)
%         idx = strmatch(noiseOptions_temp{i},{allData(:).noiseOptions},'exact');
%         for j = 1:numContrasts
%             SigMatrix = full(allData(idx).NBS_sigMat{j});
%             ConMatrix = full(allData(idx).NBS_statMat{j});
%             ConMatrix(SigMatrix == 0) = 0;
%             dlmwrite([noiseOptions_temp{i},'_con',num2str(j),'.txt'],ConMatrix)

%             % Get t-value colourbar upper and lower (for illustrator)
%             x = ConMatrix(:);
%             x = sort(x,'descend');
%             fprintf(1, '%s %s: T-value upper: %s lower: %s \n', noiseOptions_temp{i}, num2str(j), num2str(round(x(1),2)),num2str(round(x(200),2)));
%         end
%     end

%     % ------------------------------------------------------------------------------
%     % Figure 10: matrix plots
%     % ------------------------------------------------------------------------------
%     for i = 1:length(noiseOptions_temp)
%         idx = strmatch(noiseOptions_temp{i},{allData(:).noiseOptions},'exact');
        
%         % Contrast 1
%         SigMatrix = full(allData(idx).NBS_sigMat{1});
%         [~,plotMat,~] = plotClassifiedEdges(SigMatrix,ROIStructID);
%         % add constant
%         plotMat = plotMat + 1;
%         % expand
%         x = triu(plotMat);
%         x = [zeros(numROIComms,1), x]; x = [x; zeros(1,numROIComms+1)];
%         x = [zeros(numROIComms+1,1), x]; x = [x; zeros(1,numROIComms+2)];

%         % Contrast 2
%         SigMatrix = full(allData(idx).NBS_sigMat{2});
%         [~,plotMat,~] = plotClassifiedEdges(SigMatrix,ROIStructID);
%         % add constant
%         plotMat = plotMat + 1;
%         % expand
%         y = tril(plotMat);
%         y = [y, zeros(numROIComms,1)]; y = [zeros(1,numROIComms+1); y];
%         y = [y, zeros(numROIComms+1,1)]; y = [zeros(1,numROIComms+2); y];

%         % Combine
%         plotMat = x + y;

%         % Create mask
%         imAlpha = ones(size(plotMat));
%         imAlpha(plotMat == 0) = 0;
        
%         % remove the constant
%         plotMat = plotMat - 1;
%         % convert to %
%         plotMat = plotMat * 100;

%         f = figure('color','w', 'units', 'centimeters', 'pos', [0 0 7 7], 'name',[noiseOptions_temp{i},'_con',num2str(j)]); box('off');
%         cmax = 12; 
%         % cmax = ceil(max(plotMat(:)));
%         imagesc(plotMat,'AlphaData',imAlpha);
%             axis square
%             axis tight
%             colormap(BF_getcmap('reds',6,0))
%             caxis([0 cmax])
%         set(gca,'TickLength', [0 0])
%         % set(gca,'XTick',''); % Set x tick vals
%         set(gca,'XTick',1:1:numROIComms); % Set x tick vals
%         set(gca,'XTickLabel',ROILabels,'FontSize',8,'FontWeight','normal','XTickLabelRotation',90); % set x tick labels. these will be off and need to be corrected in illustrator
%         set(gca,'YTick',1:1:numROIComms); % set y tick vals
%         set(gca,'YTickLabel',ROILabels,'FontSize',8,'FontWeight','normal'); % set y tick labels
%         colorbar; % show colour bar

%         set(f, 'PaperPositionMode', 'auto')
%         print(f,noiseOptions_temp{i},'-dsvg')
%     end
% end

% save(fullfile(plotdir, 'allData.mat'),'allData','-v7.3');

% clear pval* vec temp FC*
