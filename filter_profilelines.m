function valid = filter_profilelines(angles, spans)

da = [0,diff(angles)];
tol = 10*pi/180;
dl = [0,diff(spans)];
mean_ang = mean(angles);
dm_ang = angles - mean_ang;

mean_sp = mean(spans);
valid = (spans>0) & (da>-tol) & (da<tol) & (abs(dl)<20) & (spans<mean_sp+20);

% first could be outlier - diff will not catch 
if abs(dm_ang(1))>tol
    valid(1)=false;
end

