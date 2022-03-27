function [ts_q,ts_dq,ts_x] = solpath2ts(file,h)
% READ FILE
A = importdata(file,' ',2);

q = unwrap(A.data); %unwrap to resample
dq = zeros(size(q));
x = [q dq];

tf = length(q)*h-h;
t = 0:h:tf;

%CREATE TIMESERIES 
ts_q = timeseries(q,t);
ts_dq = timeseries(dq,t);
ts_x = timeseries(x,t);

end