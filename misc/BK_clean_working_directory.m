function BK_clean_working_directory(tmppath,outpath)
outpath=[outpath '/rwls_stats_compcor_smoothed/'];
if ~exist(outpath, 'dir')
    mkdir(outpath);
end
datapath=dir([tmppath '/rwls_stats_compcor_smoothed/*.nii']);
% datapath=dir([tmppath '/rwls_stats_compcor/*.nii']);
for file=1:length(datapath)
    data(file,1)=cellstr(strcat(datapath(file).folder,'\',datapath(file).name));
end
% stats_path = [tmppath '/rwls_stats_compcor_smoothed'];
gzip(cellstr(data),outpath);
movefile([tmppath '/rwls_stats_compcor_smoothed/SPM.mat'],[outpath '/SPM.mat']);
% movefile tmppath
% outpath = [outpath '/ses-02/func/layers'];
% run = extractBefore(tmppath,'_split');
% run = run(end);
% 
% stats_folders = dir([tmppath '/rwls_*']);
% 
% %copy every folder with 'rwls' in its name to hard drive
% for folder = 1:length(stats_folders)
%     source = [stats_folders(folder).folder '/' stats_folders(folder).name];
%     target = [outpath '/run' run '/func/' stats_folders(folder).name];
%     system(['mv ' source ' "' target '"']);
% end
% 
% %copy compcor file to hard drive
% compcor_file = dir([tmppath '/compcor_regressors*.txt']);
% system(['cp ' compcor_file.folder '/' compcor_file.name ' "' outpath '/run' run '/func"']);
% 
% %delete content of working directory
% system(['rm -r ' tmppath '/*']);
end