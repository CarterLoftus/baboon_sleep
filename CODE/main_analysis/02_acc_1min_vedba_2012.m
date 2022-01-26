function acc_1min_veba_2012()
% The code reads acc data downloaded from movebank and calculates VeDBA, using the codes calc_stat_only_vedba and calc_vedba
%% read GPS data

ds = datastore('Papio Anubis Mpala 2012.csv'); % acc data downloaded from movebank
ds.SelectedFormats(3) = {'%{yyyy-MM-dd HH:mm:ss.S}D'};  % changing formats
ds.SelectedVariableNames = {'timestamp',... 3
    'eobs_accelerations_raw',... 6
    'tag_local_identifier',... 11
    'eobs_acceleration_sampling_frequency_per_axis',...
    'individual_local_identifier'}; % 12
tt_acc = tall(ds);
tt_acc.tag_local_identifier = categorical(tt_acc.tag_local_identifier);

tt_filtered_acc = table2timetable(tt_acc);
raw_acc = tt_filtered_acc.eobs_accelerations_raw;
tag_acc = tt_filtered_acc.tag_local_identifier;
time_acc = tt_filtered_acc.timestamp;
freq_acc = tt_filtered_acc.eobs_acceleration_sampling_frequency_per_axis;
[acc_vals, tag, freq, dt_time] = gather(raw_acc, tag_acc, freq_acc, time_acc);

acc_all = table(acc_vals, tag, freq, dt_time);
%% open acc bursts to single measurements
tag_acc = categorical(acc_all.tag);
raw_acc = acc_all.acc_vals;
time_acc = acc_all.dt_time;

tag_samp={};
time_samp={};
time_samp_st = {};
x = {};
y = {};
z = {};
% seperate bursts to different rows
parfor pp = 1:length(raw_acc)
    vals = cell2mat(cellfun(@str2num, strsplit(raw_acc{pp})', 'UniformOutput', false));
    samp_num = length(vals)/3;
    %             acc_mat(pp,1:length(vals)) =  vals;
    tmp_time = time_acc(pp):...
        seconds(1/freq(pp)):...
        time_acc(pp)+seconds((samp_num)/(freq(pp)));
    
    tag_samp{pp} = repmat(tag_acc(pp),1,samp_num);
    time_samp{pp} = tmp_time(1:end-1);
    time_samp_st{pp} = repmat(time_acc(pp),1,samp_num);
    x{pp} = vals(1:3:end)';
    y{pp} = vals(2:3:end)';
    z{pp} = vals(3:3:end)';
end

tag_samp1 = [tag_samp{:}]';
time_samp1 = [time_samp{:}]';
time_samp_st1 = [time_samp_st{:}]';
x1 = [x{:}]';
y1 = [y{:}]';
z1 = [z{:}]';

% sorting data by tag \ time
[a, v] = sortrows([double(tag_samp1) datenum(time_samp1)]);
tag_samp1 = tag_samp1(v);
time_samp1 = time_samp1(v);
time_samp_st1 = time_samp_st1(v);
x1 = x1(v);
y1 = y1(v);
z1 = z1(v);

% default acc vals from calibrations
acc_const.X0 = 2036;
acc_const.Y0 = 1925;
acc_const.Z0 = 1893;
acc_const.Cx = 0.002257;
acc_const.Cy = 0.004673;
acc_const.Cz = 0.002281;

acc = table(tag_samp1, time_samp_st1, time_samp1, x1, y1, z1);
acc.Properties.VariableNames = ...
    {'tag_samp' , 'time_samp_st' , 'time_samp', 'x', 'y', 'z'};

acc.time_samp_st_num = datenum(acc.time_samp_st);
% save(['acc_2012_stage1'],'acc','-v7.3');
%% caculate VeDBA
vedba = table();
tag_names = unique(acc.tag_samp);
for yy = 1:length(tag_names)
    disp(num2str(yy))
    tag_indexes = tag_names(yy) == acc.tag_samp;
    time_acc = acc.time_samp(tag_indexes);
    tag_acc = acc.tag_samp(tag_indexes); 
    x = acc.x(tag_indexes); 
    y = acc.y(tag_indexes); 
    z = acc.z(tag_indexes); 
    
    x_cal=[(x'-acc_const.X0)*acc_const.Cx*9.81]'; % Calculate X-acceleration AccData m/sec2
    y_cal=[-1*(y'-acc_const.Y0)*acc_const.Cy*-9.81]'; % Calculate Y-acceleration AccData m/sec2
    z_cal=[(z'-acc_const.Z0)*acc_const.Cz*9.81]'; % Calculate Z-acceleration AccData m/sec2
    
    acc_time_min = (dateshift(time_acc,'start','minute'));
    st_min_times = unique(dateshift(time_acc,'start','minute'));
    tag = tag_acc(1);
    timestamp = acc_time_min(1);
    val = [];
    parfor mm = 1:length(st_min_times)
        tag(mm,1) = tag_acc(mm);
        timestamp(mm,1) = st_min_times(mm);
        min_indexes = find(st_min_times(mm) == acc_time_min);
        if length(indexes) > 2        
            stat = calc_stat_only_vedba(x_cal(min_indexes), y_cal(min_indexes), z_cal(min_indexes));
            val(mm,1) = stat(1);

        else
            val(mm,1) = NaN;
        end
    end
    vedba = [vedba ; table(tag,timestamp,val)];
end
%% save VeDBA table
writetable(vedba, 'vedba_mean_2012.csv')