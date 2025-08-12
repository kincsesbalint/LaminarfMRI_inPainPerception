function BK_statistical_analysis_pipeline_unsmoothedmaindata(subject)
% I modified VPF's script to estimate ß/con maps on the UNsmoothed main data. This is about ~40 min/subject. The
% following main steps were applied for each subject:
%
% 
%   1. copy the data from the external HD to locally.
%   2. define motion and compcor parameters for the GLM.
%   3. fit the first level GLM on the whole dataset and create contrast
%   images of the main effect of cognition, pain, and interaction.
%   4. Move ß, contrasts, and SPM.mat to the external HD and delete
%   locally.
%  
% About the STORAGE SPACE of the analysis and the interim steps:
% A temporary path is set to copy the raw files (motion corrected and
% registered to the anatomical image), therefore the calculations are
% faster.
% One run is 10GB of zipped data, but 40gb of unzipped data. And I have three runs per participant--> ~120GB.
% The output is 1GB of zipped data, which is moved back to the external HD.
% 
%
% The followings need to be specified before starting the script:
%   1. The mode of parameter estimation: classical reml OR rwls. see below typeofestimation
%   2. smoothing kernel size (this should be set to 0, as the copying
%   function then skip the smoothing step)
%   3. use of compcor regressors. We use compcor as that was the result
%   from a previous interim result.
%   4. temporary working directory (pipepath) - copy data (unzip), estimate 1st lvl GLM fit 
% 
%   How to call:for subj={7349}; BK_statistical_analysis_pipeline_unsmoothedmaindata(subj); end
% 
% 
% 
% 
% Balint Kincses
% 07.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 




    
%% specifying arguments
% Subjects:
% We invited 35 subjects in total in the scanner, but 7483 has no second session(too big head for the coil).
% We work with 34 subjects (rswl estimation did not work for subj 7376,but no issue with the classic spm estimation)
% List of all subjects:
% 7349,  7356(!),  7361,  7375,  7376(!!!),  7383,  7402,  7403,  7404,  7405
% 7408,  7415,  7426,  7434,  7443,  7445,  7449,  7453,  7455,  7457,  7469,  7485  
% 7414,  7425,  7433,  7435,  
% 7444,  7448,  7452,  7454,  7456,  7468,  7482,  7484  
% 
% 
% 7356 - the anatomical image has different name, needed to adjust the
% batch call (it works now)
% 7376 - the design matrix is not the same across runs (second run was stopped earlier, and the rwls estimation on smoothed data did not work, BUT classical spm estimation work just fine) 
%
%7453
smoothingkernel=0; %it is needed for backward copmatibility with previous version, so I do not need to modify to much thinkgs. 
compcor=[1]; %0-no compcor is done, 1-compcor is done

pipepath = 'C:\Users\lenov\Documents\DATA_tmp_egqaandorg_beforeupload\H_layer';

hpfiltcutoff=180;

typeofestimation=['rwls'];%['reml', 'rwls'];
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
    % do copying ~9min/subject
    % hardcoded the number of runs here as it is always 3    
    for run = 1:3%runs
        %create temporary folder locally
        tmpfolder = [pipepath '/run' num2str(run)];
        if ~exist(tmpfolder, 'dir')
            mkdir(tmpfolder);
        end
        %smoothing kernel set to 0--> skip smoothing step, BUT keep
        %unzipping still.
        BK_prepare_working_directory(tmpfolder,subjpath,run,smoothingkernel);
        sprintf('The time for copy and smooth one run:')
        toc
%         subjpathforfunctionaldata=[subjpath '/ses-02/func/layers/'];
    end %end of run for copying data
            
    %% 2. define motion and compcor params
    %
    for run=1:3
        for cmpcor=compcor
            if cmpcor
%                 typofcompcor='rwls_stats_compcor';
                typofcompcor='stats_compcor';
            else
%                 typofcompcor='rwls_stats';
                typofcompcor='stats';
            end
            
            % define the path to the images in each run:
            runpath=[subjpath '/ses-02/func/layers/run' num2str(run) '/func/'];
            % define the SPM.mat file, the rationale of this is that
            % the compcor file is not for all participants available,
            % therefore I use the SPM from the previous run by VPF.
            % the reml estimation used the same motion and copmcor files so
            % it fine here to hardcode the rwsl_func folder for looking
            % these files
            spmrunpath=[subjpath '/ses-02/func/layers/run' num2str(run) '/func/rwls_' typofcompcor '/SPM.mat'];

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
            typofcompcor='stats_compcor';
        else
            typofcompcor='stats';
        end
        
        %creat output folder
        %figure out this part
        spmoutputpath = [pipepath '/' typeofestimation '_' typofcompcor '_UNsmoothed_hpf' num2str(hpfiltcutoff)];
        if ~exist(spmoutputpath, 'dir')
            mkdir(spmoutputpath)
        end
        
        %the rwls estimation (BK_execute_SPM_batch_layers_main_allrun) is not recommended on smoothed data, therefore
        %the group estimation is done with the classic default pipeline.
        %SPM (ReML estimation - BK_execute_SPM_batch_smoothedmaindata)
        %However, the unsmoothed data can be run any of it. I need to
        %decide which one to select.
        
        if strcmp(typeofestimation,'rwls')
            
            BK_execute_SPM_batch_layers_main_allrun(pipepath,subjpath,structpath,all_multireg_file.(typofcompcor),spmoutputpath,hpfiltcutoff); %rwls estimation
            %change back that the files first locally copied and (unzipped
            %is needed)
        elseif strcmp(typeofestimation,'reml')
            %TODO: change this function as well ti fit with the current
            %implementation in terms of hte arfuments.
            BK_execute_SPM_batch_UNsmoothedmaindata(pipepath,subjpath,structpath,all_multireg_file.(typofcompcor),typofcompcor,hpfiltcutoff);%reml estimation
        end
         
        sprintf('The time for run the GLM:')
        toc
    end


    %% 4. save data to the external HD
    for cmpcor=compcor
        if cmpcor
            typofcompcor='stats_compcor';
        else
            typofcompcor='stats';
        end
        spmoutputpath = [pipepath '/' typeofestimation '_' typofcompcor '_UNsmoothed_hpf' num2str(hpfiltcutoff)];
%         spmoutput = dir(spmoutputpath);
        spmoutputniftis = dir([spmoutputpath '/*.nii']);
        outpath=[subjpath '\ses-02\func\layers\' typeofestimation '_' typofcompcor '_UNsmoothed_hpf' num2str(hpfiltcutoff)];
        if ~exist(outpath, 'dir')
            mkdir(outpath);
        end
        %only 23 sec.ok,do not need to adjust more.
    %             tic
        parfor vols=1:length(spmoutputniftis)                
             gzip(cellstr([spmoutputniftis(vols).folder '\' spmoutputniftis(vols).name]),outpath);
    %                  gzip(cellstr([datapath(vols).folder '\' datapath(vols).name]));
        end
    %             sprintf('The time for copy the data back to the extHD:')
    %             toc
       movefile([spmoutputpath '/SPM.mat'],[outpath '/SPM.mat'])
       
       rmdir(spmoutputpath,'s')
       
    end
    rmdir([pipepath '/run*'],'s')
    sprintf('The time to process one subject:')
    toc

    clear;%clc;


