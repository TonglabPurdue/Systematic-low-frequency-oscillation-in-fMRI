% (C) Copyright 2018                TonglabPurdue
% All rights reserved               Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
% Yunjie Tong, 2012

function [output]= Delay_map(imag_nii,GM,o_filename,TR, nTR, thresh)

%you will need the function "oversample_ts.m" and the function in NIfTI_30250233.zip to run this function.

%imag_nii should be an input 4D fMRI scan.

%GM is global mean of the imag_nii. You can applied FSL function: fslmeants
%to calculate the global mean. And then put into this function.

%o_filename is the grouping name for output files.

%TR is the temporal resolution of fMRI scan. TR=2 means that it took MRI 2
%seconds to complete one scan of brain.

%nTR is the step, n is the ratio to the TR; if the step is twice the TR,
%then nTRis 2

%thresh is the threshold for maximum max correlation coefficient (MCCC) after calculation of cross-correlation. 


ratioNtoOld=floor(TR/0.15); %upsample ratio; TR/ratioNtoOld is the new temporal resolution

xcorr_range=floor(6/(TR/ratioNtoOld));% for most subject in CBV (cerebral blood velocity) the time range is 6s. FOr Sub10 is 10s.


disp(['The values of xcorr_range: ', num2str(xcorr_range)]);
disp(['The values of ratioNtoOld: ', num2str(ratioNtoOld)]);
disp(['The values of step: ', num2str(nTR*TR)]);

%xcorr_range in TRs
ref_ts=detrend(oversample_ts(GM,ratioNtoOld),'linear');

gz=strfind(imag_nii,'.gz');
if (gz~=0)
    F1=gunzip(imag_nii);
    image1= str2mat(F1);
    imag=load_untouch_nii(image1);
    delete(image1);
else
    imag=load_untouch_nii(imag_nii);
end

AA1=size(imag.img);

%thresh=1.96/sqrt(length(NIRS_ts));
%thresh=0;
groupfactor=1;

%%11 
%New_img=zeros(AA1(1),AA1(2),AA1(3),numofreg);

 
New_img=zeros(AA1(1),AA1(2),AA1(3));
New_img_r=zeros(AA1(1),AA1(2),AA1(3));


New_img_movie=zeros(AA1(1),AA1(2),AA1(3), floor((xcorr_range*2+1)/ratioNtoOld/nTR)+1);% (xcorr_range*2+1)/ratioNtoOld is num of TRs in the whole range;



%% trial
for i=1:(AA1(1))
    for j=1:(AA1(2))
        for k=1:(AA1(3))
            
            if (imag.img(i,j,k,1)~=0)
                
                ts=reshape(imag.img(i,j,k,:),AA1(4),1);
                ts_detrend=detrend(oversample_ts(ts,ratioNtoOld),'linear');

%% calculate xcorr coef with ONE NIRS ts   
                
                [r,p]=xcorr(ts_detrend,ref_ts,xcorr_range,'coeff');
                
                if (max(r)>thresh) %max r = MCCC
                T=p(find(r==max(r))); %T = corresponding time delay of MCCC
                New_img(i,j,k)= T(1)*(TR/ratioNtoOld); %Delay map
                
                New_img_r(i,j,k)= max(r); %MCCC map
                Bin=floor((xcorr_range+T(1)+groupfactor)/ratioNtoOld/nTR)+1;
                New_img_movie(i,j,k,Bin)=1;
                end
                
            end
        end
    end
end

%% end of trial
output=New_img;
fname=strcat(o_filename,'_delayMap.nii')
imag.hdr.dime.dim(5)=1;
imag.img=New_img;
save_untouch_nii(imag,fname);
gzip(fname);
delete(fname);

fname=strcat(o_filename,'_maxccMap.nii')
imag.img=New_img_r;
save_untouch_nii(imag,fname);
gzip(fname);
delete(fname);

fname=strcat(o_filename,'_movie.nii')
imag.hdr.dime.dim(5)=floor((xcorr_range*2+1)/ratioNtoOld/nTR)+1;
imag.hdr.dime.pixdim(5)=nTR*TR;
imag.img=New_img_movie;
save_untouch_nii(imag,fname);
gzip(fname);
delete(fname);

