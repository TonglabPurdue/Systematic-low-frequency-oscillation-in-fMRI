% (C) Copyright 2018                TonglabPurdue
% All rights reserved               Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
% Yunjie Tong, 2012

function TS=oversample_ts(ts,ratioNtoOld)


%% oversample
%ratioNtoOld=5 mean oversample 5 times
[m,n]=size(ts);
kk=round(m*ratioNtoOld);
TS=zeros(kk,n);
for i=1:n
    TS(:,i)=interpft (ts(:,i), kk);
end
