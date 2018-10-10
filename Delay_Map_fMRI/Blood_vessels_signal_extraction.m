% (C) Copyright 2018                TonglabPurdue
% All rights reserved               Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
% James Hao Wang & Jinxa (Fiona) Yao, 2018

clc
clear
close all

setenv('FSLDIR','/usr/local/fsl');
setenv('FSLOUTPUTTYPE','NIFTI_GZ');

tic
% list = {'09ZE6UUK','1R7RYN1P','2C7F1286','2NLERKN7','3C0FHDGM','4LUBR4AE','6MH7AY14','7N9P7K46','8UF58XPK','9AEE483X','9X7P9U2G'};
% list = {'test','test1','test2','test3','test4','test5','test6','test7'};
list = {'146432','147737','148840','149539','151223','153025','154734','159340','160123','163129','176542','190031','192540','214423','221319','414229','672756','751348','856766','899885'};
% 147737 153025_2, 159340_2, 160123 
list = {'153025','159340'};
grandun_folder = '/Users/aurous/Desktop';
subject_loc = '/Volumes/4T/20_good_subject';
start_run = 1;
end_run = 2;
LRRL = {'LR','RL'};
for i = 1:length(list)
    for j = start_run:end_run
        run = num2str(j);
        
        ID = char(list(i));      %file number
        % number = '211417'
        pixel_expansion = 0    %number of pixels to add on surface
        Value = 15000
        low_value = 8000
        cutting = 1
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Artery Calc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        command = ['python ',grandun_folder,'/gradunwarp-master-2/gradunwarp/core/gradient_unwarp.py  ',subject_loc,'/',ID,'_strc/unprocessed/3T/T1w_MPR',run,'/',ID,'_3T_T1w_MPR',run,'.nii.gz ',subject_loc,'/',ID,'_strc/unprocessed/3T/T1w_MPR',run,'/T1w_warped.nii.gz siemens -g ',grandun_folder,'/gradunwarp-master-2/gradunwarp/core/grad.grad -n'];
        %     clipboard('copy',command)
        status = system(command)
        %
        command = ['python ',grandun_folder,'/gradunwarp-master-2/gradunwarp/core/gradient_unwarp.py ',subject_loc,'/',ID,'_strc/unprocessed/3T/T2w_SPC',run,'/',ID,'_3T_T2w_SPC',run,'.nii.gz ',subject_loc,'/',ID,'_strc/unprocessed/3T/T2w_SPC',run,'/T2w_warped.nii.gz siemens -g ',grandun_folder,'/gradunwarp-master-2/gradunwarp/core/grad.grad -n'];
        %     clipboard('copy',command)
        status = system(command)
        
        command = ['${FSLDIR}/bin/fslreorient2std ',subject_loc,'/',ID,'_strc/unprocessed/3T/T1w_MPR',run,'/T1w_warped.nii ',subject_loc,'/',ID,'_strc/unprocessed/3T/T1w_MPR',run,'/sub-T1w_reorient.nii.gz'];
        status = system(command)
        
        command = ['${FSLDIR}/bin/fslreorient2std ',subject_loc,'/',ID,'_strc/unprocessed/3T/T2w_SPC',run,'/T2w_warped.nii ',subject_loc,'/',ID,'_strc/unprocessed/3T/T2w_SPC',run,'/sub-T2w_reorient.nii.gz'];
        status = system(command)
        
        command = ['${FSLDIR}/bin/bet ',subject_loc,'/',ID,'_strc/unprocessed/3T/T1w_MPR',run,'/sub-T1w_reorient.nii.gz ',subject_loc,'/',ID,'_strc/unprocessed/3T/T1w_MPR',run,'/sub-T1w_brain.nii.gz'];
        status = system(command)
        % command = strcat('fsl5.0-fslchfiletype NIFTI /media/aurous/4T/ABCD_M_110/fmriresults01/sub-NDARINV',ID,'/ses-baselineYear1Arm1/anat/sub-T1w_brain');
        % status = system(command)
        
        %loads files and creates copies of the images for editing
        cd('/Users/aurous/Desktop/Nifti')
        T1 = load_untouch_nii(strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/T1w_MPR',run,'/sub-T1w_reorient.nii.gz'));
        T2 = load_untouch_nii(strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/T2w_SPC',run,'/sub-T2w_reorient.nii.gz'));
        
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
        
        
        squad_div.img = (double(T1.img) .^ 2) ./ double(T2.img);
        div_squared.img = double(T1.img) ./ (double(T2.img).^2);
        
        % calculated_matrix.img = (double(T1.img) - double(T2.img)) / (double(T1.img) + double(T2.img));
        calculated_matrix.img = squad_div.img - div_squared.img;
        % save_untouch_nii(calculated_matrix,strcat('/media/aurous/4T/LLHCP2/Calc/T1T2calc_',number,'.nii.gz'));
        
        
        %readjusting values
        % matrix = reshape(calculated_matrix.img,[],1);
        % maximum = max(matrix);
        % minimum = min(matrix);
        % desired_max = 32767;
        %
        % % square_img = (original.img - minimum) .^ 2;
        % working.img = (desired_max / (maximum - minimum)) * (calculated_matrix.img - minimum);
        % % working.img = ceil((double(desired_max) / double(max(reshape(square_img,[],1)))) * (square_img));
        % matrix = (desired_max / (maximum - minimum)) * (matrix - minimum);
        % avg = mean(matrix);
        
        working.img = calculated_matrix.img;
        save_untouch_nii(working,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/T1T2calc1_',run,'_',ID,'.nii.gz'));
        
        if cutting == 1
            idx1_cut = find(working.img>Value);
            cut_temp.img(idx1_cut) = working.img(idx1_cut);
            cc = bwconncomp(cut_temp.img);
            numPixels = cellfun(@numel,cc.PixelIdxList);
            maximum = max(numPixels);
            location_cut = find(numPixels == maximum);
            cut.img(cc.PixelIdxList{location_cut}) = working.img(cc.PixelIdxList{location_cut});
            working.img = working.img - cut.img;
            % save_untouch_nii(cut,strcat('/media/aurous/4T/LLHCP2/Calc/V2T1T2calc_',number,'.nii.gz'));
            save_untouch_nii(working,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/T1T2calc2_',run,'_',ID,'.nii.gz'));
        end
        
        
        % figure(1)
        % histogram(working.img)
        
        % find min of top 10% of values
        % Ms = sort(matrix,'descend');
        
        threshold = Value;
        top = dimension(3);
        bot = 1;
        left = 1;
        right = dimension(1);
        front = 1;
        back = dimension(2);
        
        % crop
        LR_Val = round((right - left)/10);
        FB_Val = round((back - front)/10);
        LR_range_left = (left+3*LR_Val):(left+5*LR_Val);
        LR_range_right = (left+5*LR_Val):(left+7*LR_Val);
        LR_range_vein = (left+ceil(4.7*LR_Val)):(left+ceil(5.3*LR_Val));
        %LR_range = left:right;
        FB_range = (front+ceil(4.5*FB_Val)):(front+ceil(5.5*FB_Val));
        TB_range = (bot:ceil(0.4*top));
        
        % for n=1:1:round(top/2)
        %     bg(LR_range,FB_range,n)=imcrop(working.img(:,:,n),[(back+FB_Val),(left+LR_Val),FB_Val,LR_Val]);
        %     crop.img(:,:,n) = bg(:,:,n);
        % end
        % crop_matrix_left = reshape(working.img(LR_range_left,FB_range,TB_range),[],1);
        % crop_matrix_right = reshape(working.img(LR_range_right,FB_range,TB_range),[],1);
        crop_left.img(LR_range_left,FB_range,TB_range) = working.img(LR_range_left,FB_range,TB_range);
        crop_right.img(LR_range_right,FB_range,TB_range) = working.img(LR_range_right,FB_range,TB_range);
        crop_vein.img(LR_range_vein,[front:back],[ceil(0.6*top):top]) = working.img(LR_range_vein,[front:back],[ceil(0.6*top):top]);
        
        save_untouch_nii(crop_left,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/left_area_',run,'_',ID,'.nii.gz'));
        save_untouch_nii(crop_right,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/right_area_',run,'_',ID,'.nii.gz'));
        save_untouch_nii(crop_vein,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/vein_area_',run,'_',ID,'.nii.gz'));
        
        
        % value search regions
        idx1_left = find(crop_left.img>Value);
        idx1_right = find(crop_right.img>Value);
        
        % binary image
        cropping_left.img(idx1_left) = 1;
        cropping_right.img(idx1_right) = 1;
        
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
        
        % figure(2)
        % histogram(crop_matrix)
        % Ms = sort(crop_matrix,'descend');
        % crop_Value = min(Ms(1:ceil(length(Ms(:))*0.05)))
        
        % ICA_left.img = smooth3(ICA_left.img);
        % ICA_right.img = smooth3(ICA_right.img);
        
        
        save_untouch_nii(ICA_left,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/cleft_',run,'_',ID,'.nii.gz'));
        save_untouch_nii(ICA_right,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/cright_',run,'_',ID,'.nii.gz'));
        
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
        
        idx1_vein = find(crop_vein.img>low_value);
        cropping_vein.img(idx1_vein) = 1;
        cc = bwconncomp(cropping_vein.img,18);
        numPixels = cellfun(@numel,cc.PixelIdxList);
        maximum = max(numPixels);
        location_vein = find(numPixels == maximum);
        ICA_vein.img(cc.PixelIdxList{location_vein}) = cropping_vein.img(cc.PixelIdxList{location_vein});
        %ICA_vein.img = smooth3(ICA_vein.img);
        save_untouch_nii(ICA_vein,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/cvein_',run,'_',ID,'.nii.gz'));
        
        % save_untouch_nii(crop,strcat('/media/aurous/4T/cropsquare_',number,'.nii.gz'));
        % save_untouch_nii(working,strcat('/home/james/Desktop/Patient Files/LS',number,'/adjustedT1_',number,'.nii.gz'));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SSS calc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        bet = load_untouch_nii(strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/T1w_MPR',run,'/sub-T1w_brain.nii.gz'));
        bet.img(find(bet.img ~= 0)) = 1;
        
        save_untouch_nii(bet,strcat(subject_loc,'/',ID,'_strc/unprocessed/3T/cGM_',run,'_',ID,'.nii.gz'));
        
        
        %        mkdir(strcat('/media/aurous/4T/ABCD_M_110/fmriresults01/sub-NDARINV',ID,'/ses-baselineYear1Arm1/run_0',run));
        
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
            
            %         fid = fopen(strcat('/media/aurous/4T/ABCD_M_110/fmriresults01/sub-NDARINV',ID,'/ses-baselineYear1Arm1/func/sub-NDARINV',ID,'_ses-baselineYear1Arm1_task-rest_run-0',run,'_motion.tsv'));
            %         C = textscan(fid, '%f %f %f %f %f %f %f', 'HeaderLines', 1);
            %         fclose(fid);
            %         mat = [C{1,2}(:),C{1,3}(:),C{1,4}(:),C{1,5}(:),C{1,6}(:),C{1,7}(:)];
            %         dlmwrite(strcat('/media/aurous/4T/ABCD_M_110/fmriresults01/sub-NDARINV',ID,'/ses-baselineYear1Arm1/run_0',run,'/motion.mat'),mat,'delimiter','\t','precision',2);
            %
            % %            command = strcat('fsl5.0-Text2Vest /media/aurous/4T/ABCD_M_110/fmriresults01/sub-NDARINV',ID,'/ses-baselineYear1Arm1/run_0',run,'/motion.txt /media/aurous/4T/ABCD_M_110/fmriresults01/sub-NDARINV',ID,'/ses-baselineYear1Arm1/run_0',run,'/motion.mat');
            % %            status = system(command)
            %
            %         command = strcat('fsl5.0-fsl_glm -i /media/aurous/4T/ABCD_M_110/fmriresults01/sub-NDARINV',ID,'/ses-baselineYear1Arm1/run_0',run,'/filtered_func_data.nii.gz -d /media/aurous/4T/ABCD_M_110/fmriresults01/sub-NDARINV',ID,'/ses-baselineYear1Arm1/run_0',run,'/motion.mat --out_res=/media/aurous/4T/ABCD_M_110/fmriresults01/sub-NDARINV',ID,'/ses-baselineYear1Arm1/run_0',run,'/filtered_func_data.nii.gz');
            %         status = system(command)
            
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

% open_FSLeyes = strcat('fsleyes /home/james/Desktop/Patient\ Files/LS',number,'/adjustedT1_',number,'.nii.gz');
% system(open_FSLeyes)
% system(strcat('fsleyes /home/james/Desktop/Patient\ Files/LS',number,'/test_',number,'.nii.gz'));
