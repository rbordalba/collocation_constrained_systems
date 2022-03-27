function ts2sol(ts_x,ts_u,name)
N = size(ts_x.Data, 1);
n_x = size(ts_x.Data, 2);
n_q = n_x/2;
n_u = size(ts_u.Data, 2);
t = [0; diff(ts_x.Time)];
fid = fopen(name, 'wt');

% write trajectory info
out =[num2str(N) ' ' num2str(n_x) ' ' num2str(n_u) ' 1' '\n']; 
fprintf(fid, out);
    
% write trajectory (u,t,x)
for i=1:N
    % input
    u = getdatasamples(ts_u, i);
    out =[num2str(u,12) '\n']; 
    fprintf(fid, out);
    % time increment from previous
    dt = t(i);
    out =[num2str(dt,6) '\n']; 
    fprintf(fid, out);
    % state
    x = getdatasamples(ts_x, i);
    out =[num2str(x,12) '\n']; 
    fprintf(fid, out);
end

% triple empty line
fprintf(fid, '\n\n\n');

% write again for boxes
for i=1:N
    out = ['{' num2str(n_q*2)];
    for j=1:n_q*2
        out = [out ' [' num2str(ts_x.Data(i,j),12) ',' num2str(ts_x.Data(i,j),12) ']'];
    end
    out =[out '}\n']; 
    fprintf( fid, out);
end
fclose(fid);
