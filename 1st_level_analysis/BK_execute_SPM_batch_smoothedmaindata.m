function BK_execute_SPM_batch_smoothedmaindata(tmppath,subjpath,structpath,multi_regfile,typofcompcor,hpfiltcutoff)


%load the corresponding timing file. Used to fill in the onset and duration fields of the SPM batch
for run=1:3
    timingfilepath=[subjpath '\ses-02\func\layers\sub-*\eachcondition\' '*_run-0' num2str(run) '.csv'];
    timing_csv = dir([timingfilepath]);
    timing_info_file{run} = readtable([timing_csv.folder '/' timing_csv.name]);
end

%create the statistics output path
stats_path = [tmppath '/' typofcompcor '_smoothed'];

if ~(isfolder(stats_path))
    mkdir(stats_path)
end



%load the data as filenames in a struct
for run=1:3
    tmp=dir([tmppath '/run' num2str(run) '/*_Warped-to-Anat.nii']);
% tmp = dir([tmppath '/*_Warped-to-Anat.nii']);
    scanfiles = cell(size(tmp));
    for jj = 1:numel(scanfiles)
        scanfiles(jj) = cellstr([tmp(jj).folder '/' tmp(jj).name ',1']);
    end
    scans{run} = scanfiles;
end



%run the SPM batch
inputs = cell(0, 1);
spm('defaults', 'FMRI');
spm_jobman('run', batch_to_execute(stats_path,scans,timing_info_file,structpath,multi_regfile,hpfiltcutoff), inputs{:});
end

function matlabbatch = batch_to_execute(stats_path,scans,timing_info,structpath,multi_regfile,hpfiltcutoff)

%-----------------------------------------------------------------------
% Job saved on 31-Jan-2023 12:18:03 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

%% 1st lvl model setup
matlabbatch{1}.spm.stats.fmri_spec.dir = {stats_path};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 3;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
for run=1:3
    trial_types = sort(unique(timing_info{run}.trial_type));
    %the order of the trial tpyes are the following:
%     {'anticipation_high_cognition'}
%     {'anticipation_low_cognition' }
%     {'pain_high_cogn_high_pain'   }
%     {'pain_high_cogn_low_pain'    }
%     {'pain_low_cogn_high_pain'    }
%     {'pain_low_cogn_low_pain'     }
%     {'rating'                     }
    matlabbatch{1}.spm.stats.fmri_spec.sess(run).scans = scans{run};
    last_idx = 0;
    for ii = 1:length(trial_types)
        trial_type = trial_types{ii};
        onset = [timing_info{run}.onset(strcmp(timing_info{run}.trial_type,trial_type))];
        duration = [timing_info{run}.duration(strcmp(timing_info{run}.trial_type,trial_type))];
    
        matlabbatch{1}.spm.stats.fmri_spec.sess(run).cond(ii).name = trial_type;
        matlabbatch{1}.spm.stats.fmri_spec.sess(run).cond(ii).onset = onset;
        matlabbatch{1}.spm.stats.fmri_spec.sess(run).cond(ii).duration = duration;
        matlabbatch{1}.spm.stats.fmri_spec.sess(run).cond(ii).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(run).cond(ii).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(run).cond(ii).orth = 1;
 
    end

    matlabbatch{1}.spm.stats.fmri_spec.sess(run).multi_reg = multi_regfile{run};
    matlabbatch{1}.spm.stats.fmri_spec.sess(run).hpf = hpfiltcutoff;
end
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0;
subjid= extractBefore(structpath,'/ses-');
if str2double(subjid(end-3:end))==7356
    matlabbatch{1}.spm.stats.fmri_spec.mask = {[structpath '/UNI_MPRAGEised_brainmask.nii']};
else
    matlabbatch{1}.spm.stats.fmri_spec.mask = {[structpath '/UNI_MoCo_MPRAGEised_brainmask.nii']};
end
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; %or None? - see Tor's recommendation in which cases. or FAST? https://github.com/wiktorolszowy/fMRI_temporal_autocorrelation/blob/master/analysis_for_one_subject_SPM.m
%matlabbatch{1}.spm.tools.rwls.fmri_rwls_spec.cvi = 'wls'; %this is not recommended by the developers if the data is smoothed
%% Model estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0; % set to 1 to keep residuals
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
%% Contrast estimation
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('fMRI model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'maineffectofcognition';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [0 0 1 1 -1 -1 0];
% matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [0 0 0.5 0.5 -0.5 -0.5 0];%never used this part 
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'both';

matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'maineffectofpain';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 0 1 -1 1 -1 0];
% matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 0 0.5 -0.5 0.5 -0.5 0];%never used this part 
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'both';

matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'interactioneffect';
% matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 1 -1 -1 1 0]; %this was defined in the previous run with hpfilter=400sec
matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 -1 1 1 -1 0]; %this would mean a similar direction as in the behavior, that is the positive value would mean higher cognitive effect in the high pain condition
% matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 0.5 -0.5 -0.5 0.5 0]; %never used this part 
matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'both';
    %the order of the trial tpyes are the following:
%     {'anticipation_high_cognition'}
%     {'anticipation_low_cognition' }
%     {'pain_high_cogn_high_pain'   }
%     {'pain_high_cogn_low_pain'    }
%     {'pain_low_cogn_high_pain'    }
%     {'pain_low_cogn_low_pain'     }
%     {'rating'                     }
matlabbatch{3}.spm.stats.con.delete = 0;
%% diagnostic plot for - it can only be used if RWSL toolbox was used for estimation
% matlabbatch{4}.spm.tools.rwls.fmri_rwls_plotres.spmmat_plot(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
% matlabbatch{4}.spm.tools.rwls.fmri_rwls_plotres.plot_subset = [1 Inf];
% matlabbatch{4}.spm.tools.rwls.fmri_rwls_plotres.movparam = vertcat(multi_regfile{:});

end