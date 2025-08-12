function BK_statistical_analysis_pipeline(subject)

% I modified VPF's script to run rwls on the smoothed main data. The
% following main steps were applied for each subject:
%   1. copy the data from the external HD to locally and smooth with a
%   predefined kernel.
%   2. define motion and compcor parameters for the GLM
%   3. fit the first level GLM on the whole dataset and create contrast
%   images of the main effect of cogntition, pain, and interaction.
%   4. Move ß, contrasts, and SPM.mat to the external HD and delete locally
%  
% About the storage space of the analysis and the interim steps:
% A temporary path is set to copy the raw files (motion corrected and
% registered to the anatomical image), therefore the calculations are
% faster.
% One run is 10GB of zipped data, but 40gb of unzipped data. And I have three runs per participant--> ~120GB.
% The output is 1GB of zipped data, which is moved back to the external HD.
% 
%
% 
% 
% Previous runs:
%   - 10.2024: 34 subjects data were analyzed and saved in the external HD (...ELHID/ses-02\func\layers\rwls_stats_(compcor_)smoothed_hpf180)
%       used arguments:
%           - smooting kernel size (4mm)
%           - with and without compcor regressors
%           - HPfilter cutoff: 180sec
%   - 08.2024: 33 subjects data were analyzed and saved in the external HD (...ELHID/ses-02\func\layers\rwls_stats_compcor_smoothed)
%       used arguments:
%           - smooting kernel size (4mm)
%           - compcor regressors
%           - HPfilter cutoff: 400sec
%
% The followings need to be specified before starting this script:
%   1. The participan ELH ID. One can specifiy multiple participant, they
%   will analyzed sequentially.
%   2. smoothing kernel size
%   3. use of compcor regressors
%   4. temporary working directory (pipepath) - copy data, perform smoothing and 1st lvl GLM fit 
% 
% 
% 
% 
% Update 29.10.2024:
%   Whether use compcor as a nuissance regressor need be specified at the
%   beginning (see compcor var)
% 
% 
% 
% Balint Kincses
% 07.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%this script performes robust weighted linear least squares analysis on
%multiple subjects and estimates the correspoding SPM.mat.
%
%The output will be a folder named rwls_stats for each experiment (WM_localizer,
%pain_calib, post_calib, layers) and subject. In this folder, the SPM.mat and
%beta maps are saved.
%
%In case of the main experiment, the design matrix is different. Hence, it
%is treated seperately. To speed things up, data are copied to a working
%directory before running the analysis. The resulting SPM.mat file can then
%be used to run GLM on layers. See VPF_layer_analysis_stats.m.
%
% by VPF 03.2024



    
%% specifying arguments
% Subjects:
% We invited 35 subjects in total in the scanner, but 7483 has no second session(too big head for the coil).
% We work with 34 subjects (rdwl estimation did not work for subj 7376,but no issue with the classic spm estimation)
% List of all subjects:
% 7349,  7356(!),  7361,  7375,  7376(!!!),  7383,  7402,  7403,  7404,  7405
% 7408,  7415,  7426,  7434,  7443,  7445,  7449,  7453,  7455,  7457,  7469,  7485  
% 7414,  7425,  7433,  7435,  7444,  7448,  7452,  7454,  7456,  7468,  7482,  7484  
% 
% 
% 7356 - the anatomical image has different name, needed to adjust the batch call
% 7376 - the design matrix is not the same across runs (second run was stopped earlier, and the rwls estimation did not work) 
%
%
%The following subjects were run for compcor/nocompcor:
% 7375
%  7408,  7415,  7426,  7434,  7443,  7445,  7449, 7453, 7455,   7457,  7469, 7485 
%  7414,  7425,  7433,  7435,  7444,  7448,  7452,  7454, 7456,  7468,  7482,  7484
    smoothingkernel=[4];
    
    compcor=[0 ,1]; %0-no compcor is done, 1-compcor is done
    
    pipepath = 'C:\Users\lenov\Documents\DATA_tmp_egqaandorg_beforeupload\H_layer';
    
    hpfiltcutoff=180;
    
    % for subject = subjects
    tic
    subject
%IMPORTANT that the two external HD plugged in the "right" position as
%there is sligthly different path on those.
% define subjectpath on the external HD and the path of the anatomical image as
% well.
    if subject{:} <= 7405
        structpath_base = 'E:\pain_layers\main_project\derivatives\pipeline';
    else
        structpath_base = 'D:\main_project\derivatives\pipeline';
    end

    subjpath = [structpath_base '/' num2str(subject{:})];
    structpath = [structpath_base '/' num2str(subject{:}) '/ses-01/anat/presurf_MPRAGEise/presurf_UNI'];

    %% 1. copy data locally and smooth
    % do copying and smoothing ~40min/subject
    % hardcoded the number of runs here as it is always 3    
    for run = 1:3%runs
        %create temporary folder locally
        tmpfolder = [pipepath '/run' num2str(run)];
        if ~exist(tmpfolder, 'dir')
            mkdir(tmpfolder);
        end



        
        BK_prepare_working_directory(tmpfolder,subjpath,run,smoothingkernel);
        sprintf('The time for copy and smooth one run:')
        toc
    end %end of run for copying data
            
            %% 2. define motion and compcor params
            %
    for run=1:3
        for cmpcor=compcor
            if cmpcor
                typofcompcor='rwls_stats_compcor';
            else
                typofcompcor='rwls_stats';
            end

            % define the path to the images in each run:
            runpath=[subjpath '/ses-02/func/layers/run' num2str(run) '/func/'];
            % define the SPM.mat file, the rationale of this is that
            % the compcor file is not for all participants available,
            % therefore I use the SPM from the previous run by VPF.
            spmrunpath=[subjpath '/ses-02/func/layers/run' num2str(run) '/func/' typofcompcor '/SPM.mat'];

            % Create/copy the motion and compcor regressor file in the
            % temp folder. In some cases the compcor txt file is
            % missing, so I need to save that from the SPM.Sess.C.C
            %todo when we do not want to use copmcor 
            motion_file = dir([runpath '*_MoCorr.txt']);
            if isempty(motion_file)
                myspm= load(spmrunpath);
                covmatrix = myspm.SPM.Sess.C.C;
                motionmatrix=covmatrix(:,1:6);
                writematrix(motionmatrix,[tmpfolder '/MoCor_r' num2str(run) '.txt'])
                motion_file=[tmpfolder '/MoCor_r' num2str(run) '.txt'];
            else
                motion_file = [motion_file.folder '/' motion_file.name];
            end
            

            multi_regfile = [cellstr(motion_file)];
            %only load compcor when we need it.
            if cmpcor
                compcorfile = dir([runpath 'compcor_regressors*.txt']);
                if isempty(compcorfile)
                    myspm= load(spmrunpath);
                    covmatrix = myspm.SPM.Sess.C.C;
                    if width(covmatrix)==11
                        %hardcoded here the number of PC from the compcor
                        %(always 5)
                        copmcormatrix=covmatrix(:,end-4:end);

                        writematrix(copmcormatrix,[tmpfolder '/compcor_regressors_r' num2str(run) '.txt'])
                        compcorfile=[tmpfolder '/compcor_regressors_r' num2str(run) '.txt'];
                    end
                    
                else
                    compcorfile = [compcorfile.folder '/' compcorfile.name];
                    
                end
                multi_regfile = [multi_regfile; cellstr(compcorfile)];
            end
            all_multireg_file.(typofcompcor){run}=multi_regfile;
            maki{run}=multi_regfile;
        end %end of nuissance regression file creation (with and without compcor)
    end %end run for copmcorfile creation per run

    %% 3. fits first lvl GLM on the smoothed data
    %the current implementation runs with and without compcor nuissance
    %regressors
    for cmpcor=compcor
        if cmpcor
            typofcompcor='rwls_stats_compcor';
        else
            typofcompcor='rwls_stats';
        end
        %the rwls estimation (BK_execute_SPM_batch_layers_main_allrun) is not recommended on smoothed data, therefore
        %the group estimation is done with the classic default pipeline in
        %SPM (ReML estimation - BK_execute_SPM_batch_smoothedmaindata)
%         BK_execute_SPM_batch_layers_main_allrun(pipepath,subjpath,structpath,all_multireg_file.(typofcompcor),typofcompcor,hpfiltcutoff);
        BK_execute_SPM_batch_smoothedmaindata(pipepath,subjpath,structpath,all_multireg_file.(typofcompcor),typofcompcor,hpfiltcutoff);
        sprintf('The time for run the GLM:')
        toc
    end


    %% 4. save data to the external HD
    for cmpcor=compcor
        if cmpcor
            typofcompcor='rwls_stats_compcor';
        else
            typofcompcor='rwls_stats';
        end
        datapath = dir([pipepath '/' typofcompcor '_smoothed/*.nii']);
        outpath=[subjpath '\ses-02\func\layers\' typofcompcor '_smoothed_hpf180'];
        if ~exist(outpath, 'dir')
            mkdir(outpath);
        end
        %only 23 sec.ok,do not need to adjust more.
    %             tic
        parfor vols=1:length(datapath)                
             gzip(cellstr([datapath(vols).folder '\' datapath(vols).name]),outpath);
    %                  gzip(cellstr([datapath(vols).folder '\' datapath(vols).name]));
        end
    %             sprintf('The time for copy the data back to the extHD:')
    %             toc
       movefile([pipepath '/' typofcompcor '_smoothed/SPM.mat'],[outpath '/SPM.mat'])
       
       rmdir([pipepath '/' typofcompcor '_smoothed'],'s')
       
    end
    rmdir([pipepath '/run*'],'s')
    sprintf('The time to process one subject:')
    toc
% end %end subjecct
% toc
    clear;%clc;
%% the old version modified from VPF original similarly named file. Above I removed all the comments to make it more concise.
% %this script performes robust weighted linear least squares analysis on
% %multiple subjects and estimates the correspoding SPM.mat.
% 
% %The output will be a folder named rwls_stats for each experiment (WM_localizer,
% %pain_calib, post_calib, layers) and subject. In this folder, the SPM.mat and
% %beta maps are saved.
% 
% %In case of the main experiment, the design matrix is different. Hence, it
% %is treated seperately. To speed things up, data are copied to a working
% %directory before running the analysis. The resulting SPM.mat file can then
% %be used to run GLM on layers. See VPF_layer_analysis_stats.m.
% clear;clc;
% tic
% % subjects = {7402,7403,7404,7405};
% subjects = {7408};
% smoothingkernel=[4];
% % experiments = [cellstr('WM_localizer'),cellstr('pain_calib'),cellstr('post_calib')]; 
% experiments = [cellstr('layers')];
% if any(strcmp(experiments,'layers'))
%     %this is a temporary path to copy the files, bc then the whol
%     %calculation would be faster. However, I do not need this.(do not have
%     %much space on my laptop...). One run is 10GB of zipped data, but 40gb
%     %of unzipped data. And I have three runs per participant.
%     %Instead, I would unzip it temporarily for one subject and delet it
%     %after that.
% %     pipepath = 'D:\main_project\derivatives\pipeline\tmp';
%     pipepath = 'C:\Users\lenov\Documents\DATA_tmp_egqaandorg_beforeupload\H_layer';
% %     pipepath = '/home/pfaffenrot/work/postdoc/projects/ANT_workdir/';
% else
%     pipepath = '/media/pfaffenrot/My Passport1/pain_layers/main_project/derivatives/pipeline/';
% 
% end
% 
% for subject = subjects
%     for experiment = experiments
%         if strcmp(experiment{:},'layers')
%             if subject{:} <= 7405
%                 structpath_base = 'E:\pain_layers\main_project\derivatives\pipeline';
% 
% %                 structpath_base = '/media/pfaffenrot/My Passport1/pain_layers/main_project/derivatives/pipeline/';
%             else
%                 structpath_base = 'D:\main_project\derivatives\pipeline';
% 
% %                 structpath_base = '/media/pfaffenrot/My Passport2/main_project/derivatives/pipeline/';
%             end
%             subjpath = [structpath_base '/' num2str(subject{:})];
%             structpath = [structpath_base '/' num2str(subject{:}) '/ses-01/anat/presurf_MPRAGEise/presurf_UNI'];
% 
%             %create WM mask for compcor
%             %I do not need this as the copmcor regressors are already
%             %prepared(see below)
% %             VPF_compcor_mask_creation_pain(structpath);
% 
%             
%             
% %             rundir = dir([pipepath '/run' run]);
% %             runs = length(rundir);
%             %hardcoded the number fo runs here as it is always 3
%             for run = 1:3%runs
%                 
%                 tmpfolder = [pipepath '/run' num2str(run)];
%                 if ~exist(tmpfolder, 'dir')
%                     mkdir(tmpfolder);
%                 end
%                 %todo now the compcor spm is used, but one can use without
%                 %that (rwls_stats).
%                 runpath=[subjpath '/ses-02/func/layers/run' num2str(run) '/func/'];
%                 spmrunpath=[subjpath '/ses-02/func/layers/run' num2str(run) '/func/rwls_stats_compcor/SPM.mat'];
%                 %copy data to working directory
%                 %I do not copy the data, bc I do not have enough space.
% %                 VPF_prepare_working_directory(inppath,structpath);
%                 %insted, I would copy the data to a temporary file in the
%                 %external hd and run the analysis for that.
%                 % I also run now the smoothing here.
%                 %this takes 1599.283595 seconds. per subject per run.
%                  %1759.050467 seconds. per subject
% %                 BK_prepare_working_directory(tmpfolder,subjpath,run,smoothingkernel);
%                 
%                 %run compcor
%                 %do not need it as from V's run theses regressors are
%                 %available
% %                 VPF_compcor_pain(inppath,structpath)
%                 
%                 % Create/copy the motion and compcor regressor file in the
%                 % temp folder. In some cases the compcor txt file is
%                 % missing, so I need to save that from the SPM.Sess.C.C
%                 %todo when copmcor is false
%                 motion_file = dir([runpath '*_MoCorr.txt']);
%                 if isempty(motion_file)
%                     myspm= load(spmrunpath);
%                     covmatrix = myspm.SPM.Sess.C.C;
%                     motionmatrix=covmatrix(:,1:6);
%                     writematrix(motionmatrix,[tmpfolder '/MoCor_r' num2str(run) '.txt'])
%                     motion_file=[tmpfolder '/MoCor_r' num2str(run) '.txt'];
%                 else
%                     motion_file = [motion_file.folder '/' motion_file.name];
%                 end
%                 
% 
%                 multi_regfile = [cellstr(motion_file)];
%                 compcorfile = dir([runpath 'compcor_regressors*.txt']);
%                 if isempty(compcorfile)
%                     myspm= load(spmrunpath);
%                     covmatrix = myspm.SPM.Sess.C.C;
%                     if width(covmatrix)==11
%                         %hardcoded here the number of PC from the compcor
%                         %(always 5)
%                         copmcormatrix=covmatrix(:,end-4:end);
% 
%                         writematrix(copmcormatrix,[tmpfolder '/compcor_regressors_r' num2str(run) '.txt'])
%                         compcorfile=[tmpfolder '/compcor_regressors_r' num2str(run) '.txt'];
%                     end
%                     
%                 else
%                     compcorfile = [compcorfile.folder '/' compcorfile.name];
%                     
%                 end
%                 multi_regfile = [multi_regfile; cellstr(compcorfile)];
%                 all_multireg_file{run}=multi_regfile;
% 
% %                 spmrunpath
% 
% 
%                 %run 1st level analysis on voxel space
% %                 for DO_COMPCOR = [true,false]
%                 %only with compcor now
%                 
%                 for DO_COMPCOR = [true]
%                     %I use for each run the SPM.Sess.C.C for motion and
%                     %compcor regressors (in the rwls_stats_compcor subfolder)
%                     %it took Elapsed time is 591.722424 seconds.
% %                     BK_execute_SPM_batch_layers_main(tmpfolder,subjpath,structpath,DO_COMPCOR,run,multi_regfile);
%                 end
%                 
%                 %clean working directory
% %                 VPF_clean_working_directory(inppath,[structpath_base '/' num2str(subject{:})]);
% %                 BK_clean_working_directory(tmpfolder,runpath);
%             end %end run
%                             %clean working directory
%             %remove subjct wise temporary folders
%             toc
%             %1800.883564 seconds.
% %             BK_execute_SPM_batch_layers_main_allrun(pipepath,subjpath,structpath,DO_COMPCOR,all_multireg_file);
%             toc
%             datapath = dir([pipepath '/rwls_stats_compcor_smoothed/*.nii']);
% 
% %             data=strcat(vertcat(datapath(:).folder),'\',vertcat(datapath(:).name));
%             outpath=[subjpath '\ses-02\func\layers\rwls_stats_compcor_smoothed'];
%             if ~exist(outpath, 'dir')
%                 mkdir(outpath);
%             end
%             parfor vols=1:length(datapath)
%                 
%                  gzip(cellstr([datapath(vols).folder '\' datapath(vols).name]),outpath);
% %             all_multireg_file
% 
%             end
%            movefile([pipepath '/rwls_stats_compcor_smoothed/SPM.mat'],[outpath '/SPM.mat'])
%            rmdir([pipepath '/run*'],'s')
%            rmdir([pipepath '/rwls_stats_compcor_smoothed'],'s')
% %             rmdir([pipepath '/run1'])
% %             rmdir([pipepath '/run2'])
% %             rmdir([pipepath '/run3'])
% %             rmdir([pipepath '/rwls_stats_compcor_smoothed'])
%             %run 1st level analysis on layer level
% %             BK_layer_analysis_stats([structpath_base '/' num2str(subject{:})]);
%             all_multireg_file
%         else %all but layers. Done on harddrive, not working drive
% 
%             %run 1st level analysis on voxel space
%             VPF_execute_SPM_batch(pipepath,subject{:},experiment{:});
% 
%             %compress the used nifies to save space
%             experiment_path = [pipepath num2str(subject{:}) '/ses-02/func/' experiment{:}];
%             Nruns = length(dir([experiment_path '/run*']));
% 
%             for run = 1:Nruns
%                 dat_path = [experiment_path '/run' num2str(run) '/func/'];
%                 files = dir([dat_path '/*Warped-to-Anat.nii']);
% 
%                 for ii = 1:length(files)
%                     file = [files(ii).folder '/' files(ii).name];
%                     system(['pigz ' '"' file '"']);
%                 end
%             end
%         end
%     end %layer or not layer data
% end %end subjecct
% toc

