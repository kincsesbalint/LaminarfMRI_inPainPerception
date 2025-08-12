function BK_prepare_working_directory(outpath,subjpath,run,smoothingkernel)

%copy data and decompress
% first define list of files and the output file based on the inputs
datapath = dir([subjpath '/ses-02/func/layers/run' num2str(run) '/func/*Warped-to-Anat.*']);
rawfilenms=vertcat(datapath(:).name);
data=strcat(vertcat(datapath(:).folder),'\',(rawfilenms));

outputdata=strcat(outpath,'\',(rawfilenms(:,1:end-3)));

gunzip(cellstr(data),outpath);
%smooth the data locally.
if smoothingkernel~=0
    parfor vol=1:height(outputdata)
        spm_smooth(outputdata(vol,:),outputdata(vol,:),smoothingkernel);
    end
end

end

%% the old version modified from VPF original similarly named file. Above I removed all the comments to make it more concise.

% 
% function BK_prepare_working_directory(outpath,subjpath,run,smoothingkernel)
% 
% %subpath = extractBefore(structpath,'/ses-');
% % run = extractBefore(outpath,'_split');
% % run = run(end);
% 
% % inppath = [subjpath '/ses-02/func/layers'];
% 
% % %copy timing file .csv of each condition from sub-* to run folder in working dir
% % csv_file = dir([inppath '/sub*/eachcondition/*run-0' run '.csv']);
% % system(['cp -n "' csv_file.folder '/' csv_file.name '" ' outpath]);
% % 
% % %copy motion file
% % motion_file = dir([inppath '/run' run '/func/*_MoCorr.txt']);
% % system(['cp -n "' motion_file.folder '/' motion_file.name '" ' outpath]);
% 
% %copy data and decompress
% datapath = dir([subjpath '/ses-02/func/layers/run' num2str(run) '/func/*Warped-to-Anat.*']);
% rawfilenms=vertcat(datapath(:).name);
% % parfor vols=1:length(datapath)
% %     gunzip(cellstr([datapath(vols).folder,'\',rawfilenms(vols,:)]),outpath)
% %     
% % end
% data=strcat(vertcat(datapath(:).folder),'\',(rawfilenms));
% outputdata=strcat(outpath,'\',(rawfilenms(:,1:end-3)));
% % outputdata=outputdata(1:40,:)
% %check which is faster, when the tmp output file is locally(so copy
% %locally), or when it stays in the external HD.
% gunzip(cellstr(data),outpath);
% %think about smoothing
% parfor vol=1:height(outputdata)
%     spm_smooth(outputdata(vol,:),outputdata(vol,:),smoothingkernel);
% end
% % for file = 1:length(datapath)
% %     filename = [datapath(file).folder '/' datapath(file).name];
% %     system(['cp -n "' filename '" ' outpath]);
% %     [~,~,ext] = fileparts(filename);
% %     if strcmp(ext,'.gz')
% % %         maki="D:\main_project\tmp\mag_POCS_r1_1406_Warped-to-Anat.nii.gz"
% % %         outpath="D:\main_project\tmp"
% %         hami=sprintf('%s e "%s" -o"%s" -y', '"C:/Program Files/7-Zip/7z.exe"', maki, outpath);
% %         system(hami);
% %         system(['pigz -df ' outpath '/' datapath(file).name]);
% % 
% %         %OR using 
% %         gunzip(maki)
% %         %and smooth as well
% %         spm_smooth([WM_path '/spmT_0001.nii'],[WM_path '/sspmT_0001.nii'],[4.5 4.5 4.5]);
% %     end
% % end
% end
