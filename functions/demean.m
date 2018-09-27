function [demean_ts]=demean(ts)
demean_ts=detrend(ts)/std(detrend(ts));
