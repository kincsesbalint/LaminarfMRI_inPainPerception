plotidx=0;
subplotsize=[8 4];
allsubjlayer={};
topcolumnnumber=200; %used this previously
layermeans = zeros(20, length(allsubj));

myroiname='S1';
hemisphere={'left' 'right'};
figuretype='raw'; %possible input: 'raw', 'tSNR',
focusofimgs='laminar'; %two possibliities: conditions; laminar
figurelocationpath='C:\Users\lenov\Documents\layerfMRI\E_stats\images\tomanuscipt\firstdrafts\supplementary\rawboldlayerdependentlTS\';% this is now the layer dependency

clustersize=200;%200; %possible values: 200 - top 200 vertices, 0 - all vertices

roipath='C:\Users\lenov\Documents\layerfMRI_DATA\groupavg_correctBET\';

for hemishperes=1:2
%     figure('Position',[100, 100, 2400, 1200]);
    plotidx=0;
    for subj=allsubj
    %     roipath='C:\Users\lenov\Documents\layerfMRI_DATA\groupavg_correctBET\';
        rwdtmtrx=loadrawdatamatrix(subj,myroiname);
    %     columnartsfile=[roipath num2str(subj{:}) '\functionalmasks\interimdata_rwls_' myroiname 'raw.mat'];
    %     if exist(columnartsfile,"file")
    %         disp([num2str(subj{:}) ' subject columnar info r'])
    %         plotidx=plotidx+1;
    plotidx=plotidx+1;
%     subplot(5,7,plotidx);
    if clustersize==200
%             [columnwisestat,~,~]=BK_layer_sampling_pain_study_pipeline(num2cell(subj),'glmestimate',myroiname,'raw','no','maineff','200'); 
            %(subj,'glmestimate',region,'raw',visualize,typeofcon,selectedcolumns);
            %[S2_200vertices_wholecolumn,procesubjid,columndistr]=layervalues('S2','maineff+int','visualization','no');
    %         allsubjid(plotidx)=subj{:};
    %         allsubjlayer{plotidx}=layerwisestat;
    %         allsubjcoldistr{plotidx}=columndistribution;
    %         allsubjstat{plotidx}=stat;
%             [top200Values, top200Indices] =maxk(columnwisestat(hemishperes).T(:,3),topcolumnnumber);
            
            % Filter the indices where the values are positive
%             positiveIndices = top200Indices(top200Values > 0);
            % Check if any negative values were removed
%             if length(positiveIndices) < topcolumnnumber
%                 warning('The top 200 values contained negative numbers. Only %d positive values are retained.', length(positiveIndices));
%             end
            outputpath=[roipath num2str(subj) '\smootheddata\'];
            vertexwiseparamfilenm=[outputpath 'vertexwiseparamestimation_' myroiname '.mat'];
            load(vertexwiseparamfilenm,'columnwisestat')


            tmpdata=rwdtmtrx.interimdata_columns{1,2}{hemishperes,1};
            tmpdata=tmpdata(positiveIndices,:,:);

   elseif clustersize==0
            tmpdata=rwdtmtrx.interimdata_columns{1,2}{hemishperes,1};
    end
            myspm=loadtimingfile(subj);
            trialonset=findtrialonset(myspm);
            conditionnms=fieldnames(trialonset);

            %find the values corresponding for trial onset and focus on the
            %-1 to end of trial (which can be between 19 to 22 volumes
            tmpdata_avgroi=squeeze(mean(tmpdata,1));
            idx = {1:7, 8:13, 14:20};
            layersmean = zeros(3, size(tmpdata_avgroi,2), 'like', tmpdata);
            for k = 1:3
                layersmean(k,:) = mean(tmpdata_avgroi(idx{k},:), 1);
            end
            
            %create a layer by timeseries matrix (this is averaged across
            %the 200 vertex pairs)
            roimean_bylayer=layersmean;
            sequencelength=size(roimean_bylayer,2);
            for condition=1:length(conditionnms)
                mycond=conditionnms{condition};
                myvols=trialonset.(mycond);
                condts=cell(3,1);
                for trial=1:length(myvols)
                    myvol=myvols(trial);
                    for mylayers=1:3
                        if (myvol+17)< sequencelength %this has to be put here as one subject has shorter data because of the failure of the measurement.
                            if myvol==1
                                mimi=roimean_bylayer(mylayers,myvol:(myvol+19));
                                mimi=[NaN,mimi];
                                restingperiod=roimean_bylayer(mylayers,(myvol+13):(myvol+19));
                            else
                                %original, no baseline, no percent signal
                                %changecalculation
                                if myvol+19>sequencelength
                                    mimi=roimean_bylayer(mylayers,(myvol-1):end);
                                    while length(mimi)~=21
                                        mimi=[mimi,NaN];
                                    end
                                else
                                    mimi=roimean_bylayer(mylayers,(myvol-1):(myvol+19));
                                end
                                restingperiod=roimean_bylayer(mylayers,(myvol-6):(myvol-1));
                            end
                            %this is the use of a baseline before adding
                            %together the conditions
    %                         restingperiod=roimean((myvol-6):(myvol-1));
                            baslinedmimi=mimi-mean(restingperiod);
                            %this is calculation percent signal change
                            baslinedmimi_percchange=100*(mimi-mean(restingperiod))/std(restingperiod);
                            condts{mylayers}=[condts{mylayers}; baslinedmimi];
                        end
                    end
                end
                
                trialts.(mycond)=condts;

            end
%             noise_std = std(vertexmean, 0, 2);
%             layermean=mean(vertexmean,2);
%             tSNR=layermean./noise_std;
%             if strcmp(figuretype,'raw')
%                 imagedata=layermean;
%             elseif strcmp(figuretype,'tSNR')
%                 imagedata=tSNR;
%             end
%             colourcode.pain_high_cogn_high_pain='-r';
%             colourcode.pain_high_cogn_low_pain='-b';
%             colourcode.pain_low_cogn_high_pain='-.r';
%             colourcode.pain_low_cogn_low_pain='-.b';
%             legend_entries = cell(length(conditionnms), 1);
            for condition=1:length(conditionnms)
                mycond=conditionnms{condition};
                for layeridx=1:3
                    mydata(layeridx,:)=mean(trialts.(mycond){layeridx,1});
                end
%                 line_handles(condition) = plot(mydata, colourcode.(mycond), 'LineWidth', 2,'DisplayName', mycond);
                
%                 hold on
%                 legend_entries{condition} = strrep(mycond(6:end), '_', '');
                grouplvldata(plotidx).(mycond)=mydata;
            end 
%             vline(2,'r'); %this is the start of the pain stimulation, as 0sec is the anticipation (3sec is the forst start of the pain stimulation). The stimulation reaches the plateu with some delay(2sec?), also the canonical HRFis with 5.6 secdelay in the SPM.mat file, so it would mean, the max should be seen 3volume away. 
%             vline(5,'b');
%             vline(9,'r');
%             vline(14,'r');
%             legend(line_handles,legend_entries, 'Location', 'bestoutside', 'Position', [0.92, 0.15, 0.05, 0.7]);
%             
%             title([ 'Subject ' num2str(subj)]);  % add a title to each subplot
%             hold off
%             layermeans(:, plotidx) = imagedata';
    %     else
    %         disp([num2str(subj{:}) ' subject has  no sampled cluster'])
    %     end
    
    end
%     saveas(gcf, [figurelocationpath myroiname hemisphere{hemishperes} 'trialts_' num2str(clustersize) 'vertices_mean_eachsubj_baselined.jpg'], 'jpeg'); % _baselined OR _percsignchange
    
    
    colourcode.deep='-g';
    colourcode.middle='-b';
    colourcode.superficial='-y';
    
    figure;  % create a new figure to visualize the layermean matrix
    layernames={'deep','middle','superficial'};
    legend_entries = cell(length(layernames), 1);
    for condition=1:length(conditionnms)
        mycond=conditionnms{condition};
        for mylayer=1:3
            for mysubj=1:size(grouplvldata,2)
                allsubjdata(mysubj,:)=grouplvldata(mysubj).(mycond)(mylayer,:);
            end
        


            mydata=mean(allsubjdata,1,"omitnan");
            mydatastd=std(allsubjdata,1,"omitnan")/sqrt(size(allsubjdata,1));
            subplot(2,2,condition)
            line_handles(mylayer) = plot(mydata, colourcode.(layernames{mylayer}), 'LineWidth', 2,'DisplayName', mycond);
            hold on 
            errorbar(1:21, mydata, mydatastd, colourcode.(layernames{mylayer}));
    %         errorbar(mydata, std(vertcat(grouplvldata.(mycond)),[],1),  colourcode.(mycond));
            hold on
            legend_entries{mylayer} = layernames{mylayer};
        end
        title(strrep(mycond(6:end), '_', ''))
        vline(2,'r');
        vline(5,'b');
        vline(9,'r');
        vline(14,'r');
        
    end 
    
    legend(line_handles,legend_entries, 'Location', 'best');

%     plot(1:20, mean(layermeans, 2), '-r', 'LineWidth', 1);
%     hold on
%     errorbar(1:20, mean(layermeans, 2), std(layermeans, [], 2), 'k');
%     hold off
%     xlabel('laminar depth (WM to CSF)');
    %     ylabel(' raw BOLD signal');
%     title('Grand average across subjects');
    saveas(gcf, [figurelocationpath myroiname hemisphere{hemishperes} 'trialts_' num2str(clustersize) 'vertices_mean_subjavg_baselined.jpg'], 'jpeg');
end