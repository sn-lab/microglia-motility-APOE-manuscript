clear; close all;

%% Get directories, settings
currentfolder = pwd;

% Get directories using dialog box
newfolder = questdlg('Add folder of images?', 'New Folder', 'Yes', 'No', 'Yes');

if strcmp(newfolder, 'Yes')
    n = 1;
    inpfiledir{n} = uigetdir;
  
    % Loop to add more folders
    morefiles = questdlg('Add another folder?', 'Continue', 'Yes', 'No', 'No');

    while strcmp(morefiles, 'Yes')
        n = n + 1;
        inpfiledir{n} = uigetdir;
        morefiles = questdlg('Add another folder?', 'Continue', 'Yes', 'No', 'No');
    end
end

% Check if user selected any folders
if isempty(inpfiledir)
    error('No folders selected. Exiting.');
end
%% User input for analysis parameters

prompt = {'Enter number of channels in these images:', ...
    'Enter the channel number that contains microglia:', ...
    'Enter numeric values of timepoints separated by a space:'};
dlgtitle = 'Input';
fieldsize = [1 45; 1 45; 1 45];
definput = {'4', '3', '0 5 10 15 20 25 30'};
answer = inputdlg(prompt, dlgtitle, fieldsize, definput);
channels = str2double(answer{1, 1});
m_channel = str2double(answer{2, 1});
tp = str2num(answer{3, 1});

%% User input for image filtering parameters


processing_params = 'No - change values';
gaussblur = 100;
p = 1;
q = 10;
k = 0.25;
neighborhood = 10;
medfiltrad = 1;
medfiltrounds = 1;
minfilt =2;

while strcmp(processing_params, 'No - change values')
     prompt = {'p:', 'q:','k:', ...
        'Neighborhood size:', 'Median filter (radius):','#Times to perform sequential median filters:','Gaussblur:'};
    dlgtitle = 'Image filtering';
    fieldsize = [1 45; 1 45; 1 45; 1 45;1 45; 1 45; 1 45];
    definput = {num2str(p), num2str(q), num2str(k), num2str(neighborhood), num2str(medfiltrad), num2str(medfiltrounds),num2str(gaussblur)};
    answer = inputdlg(prompt, dlgtitle, fieldsize, definput);
    p = str2double(answer{1, 1});
    q = str2double(answer{2, 1});
    k = str2double(answer{3, 1});
    neighborhood = str2double(answer{4, 1});
    medfiltrad = str2double(answer{5, 1});
    medfiltrounds = str2double(answer{6, 1});
    gaussblur = str2double(answer{7, 1});
    medfilt = 2*medfiltrad+1;
    
        matfiles = dir(fullfile(inpfiledir{1}, '*.tif'));
        nfiles = numel(matfiles);
        FileName = matfiles(1).name;
        
        imgtest = (imread(fullfile(inpfiledir{1}, FileName),(channels) - (channels - m_channel)));
        imgtest = im2uint16(rescale((double(imgtest))));
   
        imgtest_blur = imgaussfilt(double(imgtest),gaussblur,'Padding','symmetric');
        imgtest_sub_back = double((imgtest)) - double((imgtest_blur));

        imgtest_sub_back(imgtest == 0) = NaN;
        imgtest_sub_back(imgtest_sub_back < 0) = 0;

        [imgtest_thresholded, p, q] = phansalkar((imgtest_sub_back),p,q,[neighborhood neighborhood], k, 'circular');
        imgtest_thresholded(imgtest_sub_back == 0) = 0;
        
        imgtest_filt = imgtest_thresholded;

        m = 0;
        while m < medfiltrounds
        imgtest_filt = ordfilt2(imgtest_filt,round((medfilt^2)/2),ones(medfilt, medfilt));
        m = m+1;
        end
 
        fig = figure;
        t = tiledlayout(4,8);
        nexttile(1,[4 4])
        imagesc(imgtest,[prctile(imgtest,50,'all') prctile(imgtest,97,'all')])
        colormap gray
        title('Original Image')
        nexttile(5, [4 4])
        imshow(imgtest_filt)
        title('Background Subtract - Median Filtered - Thresholded')
        fig.Position = [10 100 1250 600];

        processing_params = questdlg('Do filtering settings look good?', 'Processing Parameters', 'Yes', 'No - change values', 'Yes');
        clear local_thresh
end
   
%close all
% Create a table to store the input parameters
T = table(channels, m_channel, tp, p,q,k, neighborhood, medfiltrad, medfiltrounds,gaussblur);
%% Process each folder of images

for f = 1:numel(inpfiledir)
    close all;
    date_now = datestr(now, 'dd-mm-yy-HH-MM-SS');
    cd(inpfiledir{f});
    foldername = ['Motility_', date_now];
    mkdir(foldername);
    dirresults = fullfile(inpfiledir{f}, foldername);

    %% Read image files and perform pre-processing

    matfiles = dir(fullfile(inpfiledir{f}, '*.tif'));
    nfiles = numel(matfiles);
    FileName = cell(nfiles, 1);

    for i = 1:nfiles
        i
        FileName{i} = matfiles(i).name;
        img = cell(1, numel(imfinfo(append(inpfiledir{f},'/',string(FileName(i,1)))))/channels);
        for n = 1:numel(imfinfo(append(inpfiledir{f},'/',string(FileName(i,1)))))/channels
            img{1, n} = (imread(fullfile(inpfiledir{f}, FileName{i}), (channels * n) - (channels - m_channel)));
        end
  
  %% Gaussian Blur

    img_blur  = img;
    img_fill_zeros = img;
 
        for n = 1:numel(imfinfo(append(inpfiledir{f},'/',string(FileName(i,1)))))/channels
           img_fill_zeros{1,n}(img{1,n} == 0) = 1.1*mean(img{1,n},'all');
           img_blur{1,n} = imgaussfilt(double(img_fill_zeros{1,n}),gaussblur,'Padding','symmetric');
        end

    %% Subtract background

        for n = 1:numel(imfinfo(append(inpfiledir{f},'/',string(FileName(i,1)))))/channels
           img_sub_back{1,n} = double(img{1,n})-double(img_blur{1,n});
           img_sub_back{1,n}(img{1,n} == 0) = NaN;
           img_sub_back{1,n}(img_sub_back{1,n} < 0) = 0;
           img_nan{1,n} = img_sub_back{1,n};
        end

    %% Threshold

    cd(currentfolder)

        for n = 1:numel(imfinfo(append(inpfiledir{f},'/',string(FileName(i,1)))))/channels  
         img_thresholded{1,n} = phansalkar((img_nan{1,n}),p,q,[neighborhood neighborhood], k,'symmetric');
         img_thresholded{1,n}(img_sub_back{1,n} == 0) = 0;

        end
      
     %% Img filt for prctile of only non-zero pixels

        for n = 1:numel(imfinfo(append(inpfiledir{f},'/',string(FileName(i,1)))))/channels
           img_thresh_nan{1,n} = double(img_thresholded{1,n});
           img_thresh_nan{1,n}(img{1,n} == 0) = NaN;
        end

     %% Median filter 
    
    img_thresh_filt = img_thresh_nan;
    medfilt = medfiltrad*2+1;

        for n = 1:numel(imfinfo(append(inpfiledir{f},'/',string(FileName(i,1)))))/channels
            m = 0;
            while m < medfiltrounds
           img_thresh_filt{1,n} = ordfilt2(img_thresh_filt{1,n},round((medfilt^2)/2),ones(medfilt, medfilt));
           m = m + 1;
           end
        end
     
      %% Exclude missing timepoints

        for n = 1:numel(imfinfo(append(inpfiledir{f}, '/', string(FileName(i, 1))))) / channels
            if mean(mean([img{1, n}])) == 0
                img_thresh_filt{1, n} = NaN(size(img{1, n}));

            elseif mean(mean([img{1, n}])) < 1
                img_thresh_filt{1, n} = NaN(size(img{1, n}));

            elseif mean(mean([img{1, n}])) == 1
                img_thresh_filt{1, n} = NaN(size(img{1, n}));
            end
        end
   
    %% Difference between Frames    
    
    img_diff = cell(1,numel(imfinfo(append(inpfiledir{f},'/',string(FileName(i,1)))))/channels-1);

        for n = 1:numel(imfinfo(append(inpfiledir{f},'/',string(FileName(i,1)))))/channels-1
            img_diff{1,n} = double(img_thresh_filt{1,n}) + 2*double(img_thresh_filt{1,n+1});
        end
    
    %% Distinguish Extended, Retracted, Stable, Empty pixels

    img_px_retracted  = img_diff;
    img_px_extended  = img_diff;
    img_px_stable  = img_diff;
    img_px_zero  = img_diff;

        for n = 1:width(img_diff)
        img_px_extended{1,n}(img_diff{1,n} == 2) = 255;
        img_px_extended{1,n}(img_diff{1,n} ~= 2) = 0;

        img_px_retracted{1,n}(img_diff{1,n} == 1) = 255;
        img_px_retracted{1,n}(img_diff{1,n} ~= 1) = 0;

        img_px_stable{1,n}(img_diff{1,n} == 3) = 255;
        img_px_stable{1,n}(img_diff{1,n} ~= 3) = 0;

        img_px_zero{1,n}(img_diff{1,n} == 0) = 255;
        img_px_zero{1,n}(img_diff{1,n} ~= 0) = 0;
       end
     
    %% Stability/Instability index

    img_stability = img_diff;

        for n = 1:numel(imfinfo(append(inpfiledir{f},'/',string(FileName(i,1)))))/channels-2
           img_stability{1,n} = double(img_diff{1,n}) + 5*double(img_diff{1,n+1});
        end

        img_px_stability =  img_stability;
        img_px_instability =  img_stability;

        for n = 1:width(img_stability)
        img_px_stability{1,n}(img_stability{1,n} == 17) = 255;
        img_px_stability{1,n}(img_stability{1,n} ~= 17) = 0;

        img_px_instability{1,n}(img_stability{1,n} == 8) = 255;
        img_px_instability{1,n}(img_stability{1,n} ~= 8) = 0;

       end

        for n = 1:width(img_stability)
           stability_index(i,n) = (sum(img_px_stability{1,n},'all'))/(sum(img_px_extended{1,n},'all'));
           instability_index(i,n) = (sum(img_px_instability{1,n},'all'))/(sum(img_px_stable{1,n},'all')); 
        end
           stability_index(stability_index == 0) = NaN;
           instability_index(instability_index == 0) = NaN;

           mean_stability_index(i,:) = mean(stability_index(i,:),"omitnan");
           mean_instability_index(i,:) = mean(instability_index(i,:),"omitnan");
  
    %% Calculate Motility Index

        for n = 1:numel(imfinfo(append(inpfiledir{f},'/',string(FileName(i,1)))))/channels-1
           motility_index(i,n) = (sum(img_px_extended{1,n},'all')+sum(img_px_retracted{1,n},'all'))/...
               (sum(img_thresh_filt{1,n},'all')*255); %motility index calculated by (extendedPx + retractedPx)/(total thresholdedPx in tp)
           s_motility_index(i,n) = (sum(img_px_extended{1,n},'all')+sum(img_px_retracted{1,n},'all'))/...
               (sum(img_px_stable{1,n},'all')); %motility index calculated by (extendedPx + retractedPx)/(total stablePx in tp)
           extension_index(i,n) = (sum(img_px_extended{1,n},'all'))/(sum(img_px_stable{1,n},'all')); 
           retraction_index(i,n) = (sum(img_px_retracted{1,n},'all'))/(sum(img_px_stable{1,n},'all')); 
        end

           motility_index(motility_index == 0) = NaN;
           s_motility_index(s_motility_index == 0) = NaN;
           extension_index(extension_index == 0) = NaN;
           retraction_index(retraction_index == 0) = NaN;

           mean_extension_index(i,:) = mean(extension_index(i,:),"omitnan");
           mean_retraction_index(i,:) = mean(retraction_index(i,:),"omitnan");
           mean_motility_index(i,:) = mean(motility_index(i,:),"omitnan");
           mean_s_motility_index(i,:) = mean(s_motility_index(i,:),"omitnan");

    %% Save Threshold, Extended, Retracted, Stable Images

    dirresults_threshold = fullfile(dirresults, 'img_thresholded');
    mkdir(dirresults_threshold);

    dirresults_pos = fullfile(dirresults, 'extended_px');
    mkdir(dirresults_pos);

    dirresults_neg = fullfile(dirresults, 'retracted_px');
    mkdir(dirresults_neg);

    dirresults_stable = fullfile(dirresults, 'stable_px');
    mkdir(dirresults_stable);

        filename_thresh = ['thresh_', FileName{i}];
        imwrite(uint8(img_thresh_filt{1, 1} * 255), fullfile(dirresults_threshold, filename_thresh));
        for n = 2:numel(imfinfo(append(inpfiledir{f},'/',string(FileName(i,1)))))/channels
            imwrite(uint8(img_thresh_filt{1, n} * 255), fullfile(dirresults_threshold, filename_thresh), "WriteMode", "append");
        end

        filename_pos = ['extended_', FileName{i}];
        imwrite(uint8(zeros(size(img_px_extended{1, 1}))), fullfile(dirresults_pos, filename_pos), "tif");
        for n = 1:numel(imfinfo(append(inpfiledir{f},'/',string(FileName(i,1)))))/channels-1
            imwrite(uint8(img_px_extended{1, n}), fullfile(dirresults_pos, filename_pos), "WriteMode", "append");
        end
        
        filename_neg = ['retracted_', FileName{i}];
        imwrite(uint8(zeros(size(img_px_retracted{1, 1}))), fullfile(dirresults_neg, filename_neg));
        for n = 1:numel(imfinfo(append(inpfiledir{f},'/',string(FileName(i,1)))))/channels-1
            imwrite(uint8(img_px_retracted{1, n}), fullfile(dirresults_neg, filename_neg), "WriteMode", "append");
        end

        filename_stable = ['stable_', FileName{i}];
        imwrite(uint8(zeros(size(img_px_stable{1, 1}))), fullfile(dirresults_stable, filename_stable));
        for n = 1:numel(imfinfo(append(inpfiledir{f},'/',string(FileName(i,1)))))/channels-1
            imwrite(uint8(img_px_stable{1, n}), fullfile(dirresults_stable, filename_stable), "WriteMode", "append");
        end
    end
     %% Save Excel sheets

    cd(dirresults);

    filename = ['Processing_parameters_',date_now,'.xlsx'];
    writetable(T,fullfile(dirresults,filename))

    filename = ['Mean_motility_index_',date_now,'.xlsx'];
    TT = cell2table(num2cell(mean_motility_index)', 'VariableNames', FileName');
    writetable(TT,fullfile(dirresults,filename))
    clear TT
    
    filename = ['Mean_s_motility_index_',date_now,'.xlsx'];
    TT = cell2table(num2cell(mean_s_motility_index)', 'VariableNames', FileName');
    writetable(TT,fullfile(dirresults,filename))
    clear TT

    filename = ['Mean_retraction_index_',date_now,'.xlsx'];
    TT = cell2table(num2cell(mean_retraction_index)', 'VariableNames', FileName');
    writetable(TT,fullfile(dirresults,filename))
    clear TT

    filename = ['Mean_extension_index_',date_now,'.xlsx'];
    TT = cell2table(num2cell(mean_extension_index)', 'VariableNames', FileName');
    writetable(TT,fullfile(dirresults,filename))
    clear TT

    filename = ['Mean_stability_index_',date_now,'.xlsx'];
    TT = cell2table(num2cell(mean_stability_index)', 'VariableNames', FileName');
    writetable(TT,fullfile(dirresults,filename))
    clear TT

    filename = ['Mean_instability_index_',date_now,'.xlsx'];
    TT = cell2table(num2cell(mean_instability_index)', 'VariableNames', FileName');
    writetable(TT,fullfile(dirresults,filename))
    clear TT

    close all 
    clear extension_index retraction_index motility_index stability_index instability_index
    clear mean_extension_index mean_retraction_index mean_motility_index mean_s_motility_index mean_stability_index mean_instability_index
end

disp('Finished motility analysis!')