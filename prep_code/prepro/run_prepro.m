%% run_prepro: 
% Copyright (C) 2017, Linden Parkes <lindenparkes@gmail.com>,

function [] = run_prepro(WhichProject,WhichSessScan,subject,smoothing,discard,slicetime,despike,detr,intnorm,bandpass)
    cfg.WhichProject = WhichProject;
    cfg.WhichSessScan = WhichSessScan;
    cfg.subject = subject;

    smoothing='after';

    % preprocessing options

    if nargin < 5
        cfg.discard = 1;
    else
        cfg.discard = discard;
    end

    if nargin < 6
        cfg.slicetime = 1;
    else
        cfg.slicetime = slicetime;
    end

    if nargin < 7
        cfg.despike = 0;
    else
        cfg.despike = despike;
    end

    if nargin < 8
        cfg.detr = 1;
    else
        cfg.detr = detr;
    end

    if nargin < 9
        cfg.intnorm = 1;
    else
        cfg.intnorm = intnorm;
    end

    if nargin < 10
        cfg.runBandpass = 1;
    else
        cfg.runBandpass = bandpass;
    end

    % ------------------------------------------------------------------------------
    % Store date and time
    % ------------------------------------------------------------------------------
    cfg.DateTime = datetime('now');

    % ------------------------------------------------------------------------------
    % Parent dir
    % ------------------------------------------------------------------------------
    % REQUIRED
    cfg.parentdir = '/data/';

    % ------------------------------------------------------------------------------
    % Add paths - edit this section
    % ------------------------------------------------------------------------------
    % where the prepro scripts are
    % REQUIRED
    cfg.scriptdir = [cfg.parentdir,'ASD/code/prepro/'];
    addpath(cfg.scriptdir)
    cfg.funcdir = [cfg.parentdir,'ASD/code/func/'];
    addpath(cfg.funcdir)

    % where spm is
    % REQUIRED
    cfg.spmdir = '/opt/spm8/';
    addpath(cfg.spmdir)

    % set FSL environments
    % REQUIRED
    cfg.fsldir = '/usr/share/fsl/5.0/bin/';
    setenv('FSLDIR',cfg.fsldir(1:end-4));
    setenv('FSLOUTPUTTYPE','NIFTI');
    setenv('LD_LIBRARY_PATH',[getenv('LD_LIBRARY_PATH'),':/usr/lib/fsl/5.0'])

    % where ICA-AROMA scripts are
    % REQUIRED
    cfg.scriptdir_ICA = '/opt/ICA-AROMA/';

    % ANTs
    % REQUIRED
    cfg.antsdir = '/opt/ANTs-1.9.v4-Linux/bin/';
    setenv('ANTSPATH',cfg.antsdir);

    % Directory to AFNI functions
    % REQUIRED
    cfg.afnidir = '/usr/lib/afni/bin/';

    % add by Hai Li
    % REQUIRED
    cfg.shell='/bin/bash';
    setenv('MATLAB_SHELL',cfg.shell);

    % add by Hai Li
    % REQUIRED
    addpath('/opt/REST_V1.8_130615/');

    % ------------------------------------------------------------------------------
    % Set project settings and parameters
    % Use WhichProject if you're juggling multiple datasets
            % cfg.order = [1:1:cfg.numSlices]; % ascending
            % cfg.order = [cfg.numSlices:-1:1]; % descending
            % cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices]; % interleaved
            % cfg.order = [2:2:cfg.numSlices, 1:2:cfg.numSlices-1]; %interleave alt
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define. 
            % cfg.refSlice = round(cfg.numSlices/2); % use for sequential order (e.g., ascending or descending)
            % cfg.refSlice = cfg.numSlices-1; % use for interleaved order
            % cfg.refSlice = cfg.numSlices; % use for interleaved alt order
    % ------------------------------------------------------------------------------
    % REQUIRED
    switch cfg.WhichProject
        case 'TR2_200'
            % Where the subjects' directories are
            cfg.datadir = [cfg.parentdir,'ASD/TR2_200/'];

            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,'rest/',cfg.subject,'/01/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,'anat/',cfg.subject,'/01/'];
                case 'Sess1_Scan2'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,'rest/',cfg.subject,'/02/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,'anat/',cfg.subject,'/02/'];
            end

            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/']; 
            
            % file name of EPI 4d file
            cfg.EPI = ['rest.nii.gz'];
            % name of t1 file.
            cfg.t1name = ['anat.nii.gz'];

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];    

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 200; 
            % Repetition time of acquistion in secs
            cfg.TR = 2;
            % Number of slices in EPI volumes.
            cfg.numSlices = 33;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices-1]; % interleaved
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define. 
            cfg.refSlice = cfg.numSlices; % use for interleaved order
            
            % Scalar value indicating spatial smoothing kernal size in mm. 
            cfg.kernel = 6;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08) 
            cfg.LowPass = 0.08;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
        case 'TR2_240'
            % Where the subjects' directories are
            cfg.datadir = [cfg.parentdir,'ASD/TR2_240/'];

            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,'rest/',cfg.subject,'/01/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,'anat/',cfg.subject,'/01/'];
                case 'Sess1_Scan2'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,'rest/',cfg.subject,'/02/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,'anat/',cfg.subject,'/02/'];
            end

            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/']; 
            
            % file name of EPI 4d file
            cfg.EPI = ['rest.nii.gz'];
            % name of t1 file.
            cfg.t1name = ['anat.nii.gz'];

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];    

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 240; 
            % Repetition time of acquistion in secs
            cfg.TR = 2;
            % Number of slices in EPI volumes.
            cfg.numSlices = 33;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            cfg.order = [1:2:cfg.numSlices,2:2:cfg.numSlices-1]; % interleaved
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define. 
            cfg.refSlice = cfg.numSlices; % use for interleaved order
            
            % Scalar value indicating spatial smoothing kernal size in mm. 
            cfg.kernel = 6;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08) 
            cfg.LowPass = 0.08;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
        case 'TR3_120'
            % Where the subjects' directories are
            cfg.datadir = [cfg.parentdir,'ASD/TR3_120/'];

            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,'rest/',cfg.subject,'/01/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,'anat/',cfg.subject,'/01/'];
                case 'Sess1_Scan2'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,'rest/',cfg.subject,'/02/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,'anat/',cfg.subject,'/02/'];
            end

            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/']; 
            
            % file name of EPI 4d file
            cfg.EPI = ['rest.nii.gz'];
            % name of t1 file.
            cfg.t1name = ['anat.nii.gz'];

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];    

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 120; 
            % Repetition time of acquistion in secs
            cfg.TR = 3;
            % Number of slices in EPI volumes.
            cfg.numSlices = 48;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            cfg.order = [1:2:cfg.numSlices-1,2:2:cfg.numSlices]; % interleaved
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define. 
            cfg.refSlice = cfg.numSlices-1; % use for interleaved order
            
            % Scalar value indicating spatial smoothing kernal size in mm. 
            cfg.kernel = 6;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08) 
            cfg.LowPass = 0.08;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;
        case 'TR3_240'
            % Where the subjects' directories are
            cfg.datadir = [cfg.parentdir,'ASD/TR3_240/'];

            switch cfg.WhichSessScan
                case 'Sess1_Scan1'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,'rest/',cfg.subject,'/01/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,'anat/',cfg.subject,'/01/'];
                case 'Sess1_Scan2'
                    % where the unprocessed EPI 4d file is
                    cfg.rawdir = [cfg.datadir,'rest/',cfg.subject,'/02/'];
                    % Directory where the t1 is
                    cfg.t1dir = [cfg.datadir,'anat/',cfg.subject,'/02/'];
            end

            % where the processed epi 4d files will be output to from prepro_base
            cfg.preprodir = [cfg.rawdir,'prepro/']; 
            
            % file name of EPI 4d file
            cfg.EPI = ['rest.nii.gz'];
            % name of t1 file.
            cfg.t1name = ['anat.nii.gz'];

            % the path and filename of the template in MNI space to which everything
            % will be normalized
            cfg.mni_template = [cfg.spmdir,'templates/T1.nii'];    

            % preprocessing settings
            % length of time series (no. vols)
            cfg.N = 240; 
            % Repetition time of acquistion in secs
            cfg.TR = 3;
            % Number of slices in EPI volumes.
            cfg.numSlices = 48;
            % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
            % See help of slicetime_epis.m for guidance on how to define)
            cfg.order = [1:2:cfg.numSlices-1,2:2:cfg.numSlices]; % interleaved
            % Reference slice for slice timing acquisition. See SlicetimeEPI.m for guidance on how to define. 
            cfg.refSlice = cfg.numSlices-1; % use for interleaved order
            
            % Scalar value indicating spatial smoothing kernal size in mm. 
            cfg.kernel = 6;
            % Low-pass cut-off for bandpass filter in Hz (e.g., .08) 
            cfg.LowPass = 0.08;
            % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
            cfg.HighPass = 0.008;

    end

    % ------------------------------------------------------------------------------
    % run prepro_base
    % REQUIRED
    runBase = 1;
    % ------------------------------------------------------------------------------
    if runBase == 1
        [cfg.tN,cfg.gm,cfg.wm,cfg.csf,cfg.epiBrainMask,cfg.t1BrainMask,cfg.BrainMask,cfg.gmmask,cfg.wmmask,cfg.csfmask,cfg.dvars,cfg.dvarsExtract,cfg.fdThr,cfg.dvarsThr,cfg.exclude,cfg.outEPI] = prepro_base(cfg);
    elseif runBase == 0
        % assumes 6P has been run
        temp = load(fullfile(cfg.preprodir,[cfg.WhichProject,'_',cfg.subject,'_',cfg.WhichSessScan,'.mat']));
        cfg.tN = temp.cfg.tN;
        cfg.gm = temp.cfg.gm;
        cfg.wm = temp.cfg.wm;
        cfg.csf = temp.cfg.csf;
        cfg.epiBrainMask = temp.cfg.epiBrainMask;
        cfg.t1BrainMask = temp.cfg.t1BrainMask;
        cfg.BrainMask = temp.cfg.BrainMask;
        cfg.gmmask = temp.cfg.gmmask;
        cfg.wmmask = temp.cfg.wmmask;
        cfg.csfmask = temp.cfg.csfmask;
        cfg.dvars = temp.cfg.dvars;
        cfg.dvarsExtract = temp.cfg.dvarsExtract;
        cfg.fdThr = temp.cfg.fdThr;
        cfg.dvarsThr = temp.cfg.dvarsThr;
        cfg.exclude = temp.cfg.exclude;
        cfg.outEPI = temp.cfg.outEPI;
    end

    % ------------------------------------------------------------------------------
    % noise correction and time series
    % ------------------------------------------------------------------------------
    % Enter the names of the noise corrections strategies you want to use into the cell 'noiseOptions'
    % This script will then loop over each strategy
    % If you only want to run one then just use something like: noiseOptions = {'24P+aCC'};
    
    % All pipelines
    % noiseOptions = {'6P','6P+2P','6P+2P+GSR','24P','24P+8P','24P+8P+4GSR','12P+aCC','24P+aCC','12P+aCC50','24P+aCC50','24P+aCC+4GSR','24P+aCC50+4GSR','ICA-AROMA+2P','ICA-AROMA+2P+GSR','ICA-AROMA+8P','ICA-AROMA+8P+4GSR','24P+8P+4GSR+SpikeReg'};
    % REQUIRED
    noiseOptions = {'6P+2P','6P+2P+GSR','ICA-AROMA+2P','ICA-AROMA+2P+GSR'};
    

    % If subject was not marked for exclusion for scrubbing, then append the JP14 pipelines
    if cfg.exclude == 0 & cfg.intnorm == 1 & cfg.runBandpass == 1
        fprintf(1, '\n\t\t Adding JP14 pipelines \n\n');
        noiseOptions2Append = {'24P+4P+2GSR+JP14Scrub'};
        noiseOptions = [noiseOptions,noiseOptions2Append];
    end

    % Loop over noise correction options
    for i = 1:length(noiseOptions)

        % Set noise correction
        cfg.removeNoise = noiseOptions{i};

        % Override smoothing order if using ICA-AROMA
        if any(~cellfun('isempty',strfind({cfg.removeNoise},'ICA-AROMA'))) == 1
            fprintf(1, '\n\t\t !!!! Forcing smoothing order to ''before'' for ICA-AROMA !!!! \n\n');
            cfg.smoothing = 'before';
        else
            cfg.smoothing = smoothing;
        end

        % define inputs to noise correction
        switch cfg.smoothing
            case {'after','none'}
                if any(~cellfun('isempty',strfind({cfg.removeNoise},'JP14Scrub'))) == 0
                    % If we run smoothing AFTER noise correction, then the input file is unsmoothed epi from prepro_base
                    cfg.CleanIn = cfg.outEPI{1};
                    % and, the nuisance inputs are the same image
                    cfg.NuisanceIn_wm = cfg.outEPI{1};
                    cfg.NuisanceIn_csf = cfg.outEPI{1};
                elseif any(~cellfun('isempty',strfind({cfg.removeNoise},'JP14Scrub'))) == 1
                    % If Power's scrubbing then its the power equivalent instead
                    cfg.CleanIn = cfg.outEPI{3};
                    cfg.NuisanceIn_wm = cfg.outEPI{3};
                    cfg.NuisanceIn_csf = cfg.outEPI{3};
                end
            case 'before'
                if any(~cellfun('isempty',strfind({cfg.removeNoise},'JP14Scrub'))) == 0
                    % If we run smoothing BEFORE noise correction, then the input file is smoothed epi from prepro_base
                    cfg.CleanIn = cfg.outEPI{2};
                    % and, the nuisance inputs are the same image
                    cfg.NuisanceIn_wm = cfg.outEPI{2};
                    cfg.NuisanceIn_csf = cfg.outEPI{2};
                elseif any(~cellfun('isempty',strfind({cfg.removeNoise},'JP14Scrub'))) == 1
                    % If Power's scrubbing then its the power equivalent instead
                    cfg.CleanIn = cfg.outEPI{4};
                    cfg.NuisanceIn_wm = cfg.outEPI{4};
                    cfg.NuisanceIn_csf = cfg.outEPI{4};
                end
        end

        % ------------------------------------------------------------------------------
        % run prepro_noise
        % ------------------------------------------------------------------------------

        [cfg.noiseTS,cfg.outdir,cfg.noiseTSz] = prepro_noise(cfg);

        % ------------------------------------------------------------------------------
        % extract time series
        % REQUIRED
        runTS = 1;
        % ------------------------------------------------------------------------------
        if runTS == 1
            cd(cfg.outdir)
            
            % Parcellation, default Gordon & Yeo, you can also add your parcellation here
            cfg.parcFiles = {[cfg.parentdir,'ASD/code/ROIs/Gordon/Parcels_MNI_222.nii'],...
                             [cfg.parentdir,'ASD/code/ROIs/Yeo/Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask_resliced.nii']};

            cfg.parcWeightGM = {'yes',...
                                'yes'};

            % Set input image for time series extraction
            cfg.ExtractIn = 'epi_prepro.nii';

            % Initialise roi time series variable
            cfg.roiTS = [];

            % Loop over parcellation files
            for i = 1:length(cfg.parcFiles)
                % Set parcellation file
                cfg.parcFile = cfg.parcFiles{i};
                % set GM weight
                cfg.weightGM = cfg.parcWeightGM{i};
                % extract time series 
                cfg.roiTS{i} = prepro_extractTS_FSL(cfg);
            end

            cfg = rmfield(cfg,{'parcFile','weightGM'});
        end

        % Save data
        save('cfg.mat','cfg')

        % ------------------------------------------------------------------------------
        % Compress final outputs
        % ------------------------------------------------------------------------------
        fprintf('\n\t\t ----- Compressing %s outputs ----- \n\n',cfg.removeNoise);
        cd(cfg.outdir)
        gzip('*.nii')
        pause(5)
        delete('*.nii')
        fprintf('\n\t\t ----- Finished. ----- \n\n');
    end

    % ------------------------------------------------------------------------------
    % Compress base & t1 outputs
    % ------------------------------------------------------------------------------
    fprintf('\n\t\t ----- Compressing base outputs ----- \n\n');
    cd(cfg.preprodir)
    gzip('*.nii')
    pause(5)
    delete('*.nii')

    cd(cfg.t1dir)
    gzip('*.nii')
    pause(5)
    delete('*.nii')

    if runBase == 1
        cd(cfg.rawdir)
        gzip('*.nii')
        pause(5)
        delete('*.nii')
    end
    fprintf('\n\t\t ----- Finished. ----- \n\n');
end
