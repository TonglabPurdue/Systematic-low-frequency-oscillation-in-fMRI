% (C) Copyright 2018                TonglabPurdue
% All rights reserved               Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
% James Hao Wang & Jinxa (Fiona) Yao, 2018

clc
clear
close all

%Setting FSL environment (may not be necessary)
setenv('FSLDIR','/usr/local/fsl');
setenv('FSLOUTPUTTYPE','NIFTI_GZ');

list = {'201111'}; %array of subject ID's
grandun_folder = '/Users/aurous/Desktop'; %grandunwarp folder location
subject_loc = '/Volumes/4T/20_good_subject'; %location of subjects
start_run = 1; %First image per subject
end_run = 2; %Last image per subject
LRRL = {'LR','RL'}; %Image scan types

for i = 1:length(list)
    for j = start_run:end_run
        run = num2str(j);
        ID = char(list(i));      %file number
        
        % Constants/Parameters
        pixel_expansion = 0    %number of pixels to add on surface (dilating masks)
        Value = 15000   %threshold for ICA
        low_value = 8000    %threshold for SSS
        cutting = 1     %removing largest mask (ideally removing fat/skull), [0 --> do not remove, 1 --> remove]
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Blood Vessel Extraction and Masking
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% - Calling grandunwarp program, reoreinting T1w and T2w images, and BET'ing the T1w image
        command = ['python ',grandun_folder,'/gradunwarp-master-2/gradunwarp/core/gradient_unwarp.py  ',subject_loc,'/',ID,'_strc/unprocessed/3T/T1w_MPR',run,'/',ID,'_3T_T1w_MPR',run,'.nii.gz ',subject_loc,'/',ID,'_strc/unprocessed/3T/T1w_MPR',run,'/T1w_warped.nii.gz siemens -g ',grandun_folder,'/gradunwarp-master-2/gradunwarp/core/grad.grad -n'];
        status = system(command)
        
        command = ['python ',grandun_folder,'/gradunwarp-master-2/gradunwarp/core/gradient_unwarp.py ',subject_loc,'/',ID,'_strc/unprocessed/3T/T2w_SPC',run,'/',ID,'_3T_T2w_SPC',run,'.nii.gz ',subject_loc,'/',ID,'_strc/unprocessed/3T/T2w_SPC',run,'/T2w_warped.nii.gz siemens -g ',grandun_folder,'/gradunwarp-master-2/gradunwarp/core/grad.grad -n'];
        status = system(command)
        
        command = ['${FSLDIR}/bin/fslreorient2std ',subject_loc,'/',ID,'_strc/unprocessed/3T/T1w_MPR',run,'/T1w_warped.nii ',subject_loc,'/',ID,'_strc/unprocessed/3T/T1w_MPR',run,'/sub-T1w_reorient.nii.gz'];
        status = system(command)
        
        command = ['${FSLDIR}/bin/fslreorient2std ',subject_loc,'/',ID,'_strc/unprocessed/3T/T2w_SPC',run,'/T2w_warped.nii ',subject_loc,'/',ID,'_strc/unprocessed/3T/T2w_SPC',run,'/sub-T2w_reorient.nii.gz'];
        status = system(command)
        
        command = ['${FSLDIR}/bin/bet ',subject_loc,'/',ID,'_strc/unprocessed/3T/T1w_MPR',run,'/sub-T1w_reorient.nii.gz ',subject_loc,'/',ID,'_strc/unprocessed/3T/T1w_MPR',run,'/sub-T1w_brain.nii.gz'];
        status = system(command)
        
        % loads files and creates copies of the images
        cd('/Users/aurous/Desktop/Nifti')
        T1 = load_untouch_nii(strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/T1w_MPR',run,'/sub-T1w_reorient.nii.gz'));
        T2 = load_untouch_nii(strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/T2w_SPC',run,'/sub-T2w_reorient.nii.gz'));
        
        %% - Initializations
        dimension = size(T1.img)
        bg = zeros(dimension);
        working = T1;
        crop_left = T1;
        crop_right = T1;
        crop_vein = T1;
        working.img = bg;
        crop_left.img = bg;
        crop_right.img = bg;
        crop_vein.img = bg;
        cropping_left = working;
        cropping_right = working;
        cropping_vein = working;
        cut = working;
        cut_temp = cut;
        ICA_left = cropping_left;
        ICA_right = cropping_right;
        ICA_vein = cropping_vein;
        squad_div = working;
        div_squared = working;
        calculated_matrix = working;
        
        %% - Use T1w and T2w to create image with increased contrast
        squad_div.img = (double(T1.img) .^ 2) ./ double(T2.img);
        div_squared.img = double(T1.img) ./ (double(T2.img).^2);
        calculated_matrix.img = squad_div.img - div_squared.img;

        working.img = calculated_matrix.img;
        save_untouch_nii(working,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/T1T2calc1_',run,'_',ID,'.nii.gz'));
        
        %% - Cut out largest mask
        if cutting == 1
            idx1_cut = find(working.img>Value);
            cut_temp.img(idx1_cut) = working.img(idx1_cut);
            cc = bwconncomp(cut_temp.img);
            numPixels = cellfun(@numel,cc.PixelIdxList);
            maximum = max(numPixels);
            location_cut = find(numPixels == maximum);
            cut.img(cc.PixelIdxList{location_cut}) = working.img(cc.PixelIdxList{location_cut});
            working.img = working.img - cut.img;
            save_untouch_nii(working,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/T1T2calc2_',run,'_',ID,'.nii.gz'));
        end
        
        %% - Cropping out evaluation zone for ICA
        threshold = Value;
        top = dimension(3);
        bot = 1;
        left = 1;
        right = dimension(1);
        front = 1;
        back = dimension(2);
        
        % set cropped regions
        LR_Val = round((right - left)/10);
        FB_Val = round((back - front)/10);
        LR_range_left = (left+3*LR_Val):(left+5*LR_Val);
        LR_range_right = (left+5*LR_Val):(left+7*LR_Val);
        LR_range_vein = (left+ceil(4.7*LR_Val)):(left+ceil(5.3*LR_Val));
        FB_range = (front+ceil(4.5*FB_Val)):(front+ceil(5.5*FB_Val));
        TB_range = (bot:ceil(0.4*top));
        
        crop_left.img(LR_range_left,FB_range,TB_range) = working.img(LR_range_left,FB_range,TB_range);
        crop_right.img(LR_range_right,FB_range,TB_range) = working.img(LR_range_right,FB_range,TB_range);
        crop_vein.img(LR_range_vein,[front:back],[ceil(0.6*top):top]) = working.img(LR_range_vein,[front:back],[ceil(0.6*top):top]);
        
        save_untouch_nii(crop_left,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/left_area_',run,'_',ID,'.nii.gz'));
        save_untouch_nii(crop_right,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/right_area_',run,'_',ID,'.nii.gz'));
        save_untouch_nii(crop_vein,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/vein_area_',run,'_',ID,'.nii.gz'));
        
        % value search regions
        idx1_left = find(crop_left.img>Value);
        idx1_right = find(crop_right.img>Value);
        
        % create binary image
        cropping_left.img(idx1_left) = 1;
        cropping_right.img(idx1_right) = 1;
        
        %% - Find ICA's for left and right side
        cc = bwconncomp(cropping_left.img,18);
        numPixels = cellfun(@numel,cc.PixelIdxList);
        maximum = max(numPixels);
        location_left = find(numPixels == maximum);
        left_count = numPixels(location_left)
        ICA_left.img(cc.PixelIdxList{location_left}) = cropping_left.img(cc.PixelIdxList{location_left});
        
        cc = bwconncomp(cropping_right.img,18);
        numPixels = cellfun(@numel,cc.PixelIdxList);
        maximum = max(numPixels);
        location_right = find(numPixels == maximum);
        right_count = numPixels(location_right)
        ICA_right.img(cc.PixelIdxList{location_right}) = cropping_right.img(cc.PixelIdxList{location_right});
        
        save_untouch_nii(ICA_left,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/cleft_',run,'_',ID,'.nii.gz'));
        save_untouch_nii(ICA_right,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/cright_',run,'_',ID,'.nii.gz'));
        
        % mask dilation
        if pixel_expansion ~= 0
            save_untouch_nii(ICA_left,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/dileft0_',run,'_',ID,'.nii.gz'));
            save_untouch_nii(ICA_right,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/diright0_',run,'_',ID,'.nii.gz'));
            ICA_leftDi = ICA_left;
            ICA_rightDi = ICA_right;
            for i=1:pixel_expansion
                nu = num2str(i);
                SE = strel('sphere', i);
                ICA_leftDi.img = imdilate(ICA_left.img,SE);
                ICA_rightDi.img = imdilate(ICA_right.img,SE);
                save_untouch_nii(ICA_leftDi,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/dileft',nu,'_',run,'_',ID,'.nii.gz'));
                save_untouch_nii(ICA_rightDi,strcat('/media/aurous/4T/ABCD_M_110/fmriresults01/sub-NDARINV',ID,'/ses-baselineYear1Arm1/diright',nu,'_',run,'_',ID,'.nii.gz'));
            end
        end
        
        %% - Find SSS
        idx1_vein = find(crop_vein.img>low_value);
        cropping_vein.img(idx1_vein) = 1;
        cc = bwconncomp(cropping_vein.img,18);
        numPixels = cellfun(@numel,cc.PixelIdxList);
        maximum = max(numPixels);
        location_vein = find(numPixels == maximum);
        ICA_vein.img(cc.PixelIdxList{location_vein}) = cropping_vein.img(cc.PixelIdxList{location_vein});
        save_untouch_nii(ICA_vein,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/cvein_',run,'_',ID,'.nii.gz'));
        
        bet = load_untouch_nii(strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/T1w_MPR',run,'/sub-T1w_brain.nii.gz'));
        bet.img(find(bet.img ~= 0)) = 1;
        
        save_untouch_nii(bet,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/cGM_',run,'_',ID,'.nii.gz'));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Motion Correction, Registration, and Timeseries Extraction
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1:length(LRRL)
            position = char(LRRL(k));
            
            command = ['${FSLDIR}/bin/fslmaths ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/fmri_corrected.nii.gz ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/prefiltered_func_data -odt float'];
            status = system(command)
            
            command = ['${FSLDIR}/bin/fslroi ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/prefiltered_func_data ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/example_func 210 1'];
            status = system(command)
            
            command = ['${FSLDIR}/bin/mcflirt -in ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/prefiltered_func_data -out ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/prefiltered_func_data_mcf -mats -plots -reffile ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/example_func -rmsrel -rmsabs -spline_final'];
            status = system(command)
            
            command = ['${FSLDIR}/bin/fsl_tsplot -i ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/prefiltered_func_data_mcf.par -t ''MCFLIRT estimated rotations (radians)'' -u 1 --start=1 --finish=3 -a x,y,z -w 640 -h 144 -o ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/rot.png'];
            status = system(command)
            
            command = ['${FSLDIR}/bin/fsl_tsplot -i ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/prefiltered_func_data_mcf.par -t ''MCFLIRT estimated translations (mm)'' -u 1 --start=4 --finish=6 -a x,y,z -w 640 -h 144 -o ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/trans.png'];
            status = system(command)
            
            command = ['${FSLDIR}/bin/fsl_tsplot -i ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/prefiltered_func_data_mcf_abs.rms -i ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/prefiltered_func_data_mcf_rel.rms -t ''MCFLIRT estimated mean displacement (mm)'' -u 1 -w 640 -h 144 -a absolute,relative -o ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/disp.png'];
            status = system(command)
            
            command = ['${FSLDIR}/bin/flirt -in ',subject_loc,'/',ID,'_strc/unprocessed/3T/T1w_MPR',run,'/sub-T1w_brain.nii.gz -ref ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/example_func.nii -out ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/highres2lowres -omat ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/highres2lowres.mat -bins 256 -cost corratio -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear'];
            status = system(command)
            
            command = ['${FSLDIR}/bin/flirt -in ',subject_loc,'/',ID,'_strc/unprocessed/3T/cleft_',run,'_',ID,'.nii.gz -ref ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/example_func.nii.gz -out ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/low_res_left.nii -applyxfm -init ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/highres2lowres.mat -interp trilinear'];
            status = system(command)
            
            command = ['${FSLDIR}/bin/flirt -in ',subject_loc,'/',ID,'_strc/unprocessed/3T/cright_',run,'_',ID,'.nii.gz -ref ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/example_func.nii.gz -out ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/low_res_right.nii -applyxfm -init ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/highres2lowres.mat -interp trilinear'];
            status = system(command)
            
            command = ['${FSLDIR}/bin/flirt -in ',subject_loc,'/',ID,'_strc/unprocessed/3T/cvein_',run,'_',ID,'.nii.gz -ref ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/example_func.nii.gz -out ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/low_res_SSS.nii -applyxfm -init ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/highres2lowres.mat -interp trilinear'];
            status = system(command)
            
            command = ['${FSLDIR}/bin/flirt -in ',subject_loc,'/',ID,'_strc/unprocessed/3T/cGM_',run,'_',ID,'.nii.gz -ref ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/example_func.nii.gz -out ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/low_res_GM.nii -applyxfm -init ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/highres2lowres.mat -interp trilinear'];
            status = system(command)
            
            command = ['${FSLDIR}/bin/fslmeants -i ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/prefiltered_func_data_mcf.nii.gz -m ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/low_res_left.nii -w -o ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/left_TS.txt'];
            status = system(command)
            
            command = ['${FSLDIR}/bin/fslmeants -i ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/prefiltered_func_data_mcf.nii.gz -m ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/low_res_right.nii -w -o ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/right_TS.txt'];
            status = system(command)
            
            command = ['${FSLDIR}/bin/fslmeants -i ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/prefiltered_func_data_mcf.nii.gz -m ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/low_res_SSS.nii -w -o ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/SSS_TS.txt'];
            status = system(command)
            
            command = ['${FSLDIR}/bin/fslmeants -i ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/prefiltered_func_data_mcf.nii.gz -m ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/low_res_GM.nii -w -o ',subject_loc,'/',ID,'_',run,'/unprocessed/3T/rfMRI_REST',run,'_',position,'/GM_TS.txt'];
            status = system(command)
        end
    end
end
