% (C) Copyright 2018                TonglabPurdue
% All rights reserved               Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
% Yunjie Tong, 2012


function [demean_ts]=demean(ts)
demean_ts=detrend(ts)/std(detrend(ts));
