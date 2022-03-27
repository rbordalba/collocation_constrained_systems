function [ts_q, ts_dq, ts_u, ts_x, ts_ddq] = sol2ts(file,fs)
% READ FILE
A = importdata(file,' ',3);
info = str2num(A.textdata{1});
n_samples = info(1);
n_x = info(2);
n_u = info(3);
n_t = info(4);

% if no dummy variables
if rem(n_x,2)==0    
    x = A.data(1:3:end-1,1:n_x); %read states
else
    x = A.data(1:3:end-1,[1:(n_x-1)/2 (n_x-1)/2+2:n_x]); %read states
end
% end state
x = [x; A.data(end,1:n_x)]; %read states  
x(:,1:size(x,2)/2) = unwrap(x(:,1:size(x,2)/2)); % unwrap angles 

u = A.data(2:3:end,1:n_u); %read inputs
u = [u; u(end,:)]; % repeat input to match state size
q = x(:,1:size(x,2)/2);
dq = x(:,size(x,2)/2+1:end);
dt = [0; A.data(3:3:end,n_t)];
t=cumsum(dt(1:end));

%CREATE TIMESERIES 
ts_q = timeseries(q,t);
ts_dq = timeseries(dq,t);
ts_u = timeseries(u,t);
ts_x = timeseries(x,t);


% numerically differentiate for acceleration?
out = diff(ts_dq.Data)./diff(ts_dq.Time);
out(end+1,:) = out(end,:);
ts_ddq = timeseries(out,ts_dq.Time);
% remove NaN
[rows,cols] = find(isnan(ts_ddq.Data));
remove = unique(rows);
ts_ddq = delsample(ts_ddq,'Index',remove);

if nargin==2
    % RESAMPLE DATA AT FIXED TIME STEP 
    n_samples = floor(t(end)*fs);

    ts_x = resample(ts_x,linspace(t(1),t(end),n_samples));
    ts_q = resample(ts_q,linspace(t(1),t(end),n_samples));
    ts_dq = resample(ts_dq,linspace(t(1),t(end),n_samples));
    ts_ddq = resample(ts_ddq,linspace(t(1),t(end),n_samples));
    ts_u = resample(ts_u,linspace(t(1),t(end),n_samples));
end

end