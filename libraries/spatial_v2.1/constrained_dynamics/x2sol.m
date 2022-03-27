function x2sol(x,name)
n_q = size(x,1)/2;

fid = fopen( name, 'wt' );

for i=1:size(x,2)
    out = ['{' num2str(n_q*2)];
    for j=1:n_q*2
        out = [out ' [' num2str(x(j,i),12) ',' num2str(x(j,i),12) ']'];
    end
    out =[out '}\n']; 
    fprintf( fid, out);
end
fclose(fid);
