function stat = calc_stat_only_vedba(x, y, z)
    fs = 12;
    tm = [1:length(x)]';
    samp = [6];
    
    dba_x = nanmean(calc_vedba(x,samp(1)));
    dba_y = nanmean(calc_vedba(y,samp(1)));
    dba_z = nanmean(calc_vedba(z,samp(1)));
    
    val = nanmean(sqrt(dba_x^2+...
                dba_y^2+...
                dba_z^2));
            
    clear samp x y z dba_x dba_y dba_z;
    vars = who;
    stat = [];

    for yy = 1:length(vars)
         stat = [stat  eval(char(vars(yy)))];
    end
   