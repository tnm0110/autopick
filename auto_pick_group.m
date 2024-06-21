function [tp_picked, ts_picked, tp_picked_2,ts_picked_2, t_Pn, t_Sn,coment] = auto_pick_group(sfile)

% % ----------------------------------------------------------------- % %
%   This script automatically pick the P and S wave arrival time. It 
%   applies sequential triggnering to search for P and S waves in the
%   group velocity windows 
%   Optimize the parameters accordingly to improve the pick accuracy  
% 
% % -------------------------------------------------------------------% %
%                           Good Luck! 
%                      Md Mohimanul Islam
%                University of Missouri Columbia
%                       mibhk@umsystem.edu
%                      updated : 03/12/2024
%
% % -------------------------------------------------------------------% %

%clear;
close all
%ddir = 'C:\Users\tnm\OneDrive - University of Missouri\MUSE\autopick\data_new\2023_07_01_01_00_06.0/';
%ddir = '/data/mica/hjbf5b/Eventlist/afad/GOOD/2023_07_08_14_02_34.0';
%sfile = '2023.206.07.57.54.0000.YB.EF016.00.DPZ.D.SAC';

%cd (ddir);

%


% reading the sac file 

try
    sac=rdsac(sfile);
    dist=sac.HEADER.GCARC;
    fs=1./sac.HEADER.DELTA;
    edp=sac.HEADER.EVDP;
catch
    disp(['Error reading SAC file: ', sfile]);
    tp_picked = NaN; tp_picked_2 = NaN;
    ts_picked = NaN; ts_picked_2 = NaN;
    t_Pn = NaN; t_Sn = NaN;
    coment = 'ERROR reading SAC file';
    return;
end
    
if isfield(sac, 'T1') && isfield(sac,'T2')
    T1=sac.HEADER.T1; T2 = sac.HEADER.T2;   % picked arrvial times
end
 
%t_Pn = get_ttime(dist,edp,'Pn');
%t_Sn = get_ttime(dist,edp,'Sn');

t_Pn = get_ttime(dist,edp,'P,p,pn');
t_Sn = get_ttime(dist,edp,'S,s,sn');

% if delta < 1 is too colse
if isnan(t_Pn)
    t_Pn = get_ttime(dist,edp,'p');
end

if isnan(t_Sn)
    t_Sn = get_ttime(dist,edp,'s');
end

n_tpn = t_Pn*fs; n_tsn = t_Sn*fs;
t=sac.HEADER.B:sac.HEADER.DELTA:sac.HEADER.E;

if length(t) ~= length(sac.d)
    t(end+1) = t(end) + sac.HEADER.DELTA;
end

% remove mean
sac2 = sac.d - mean(sac.d);
% remove trend
sac2=detrend(sac2);

% filter the seismogram
c1 = 0.05;   % Lower corner frequency in Hz
c2 = 3.0;    % Upper corner frequency in Hz
n = 3;       % Filter order (number of poles)

% Design a bandpass Butterworth filter

try
    [b, a] = butter(n, [c1, c2]/(fs/2), 'bandpass');
    sac_filt = filtfilt(b, a, sac2); % apply zero-phase bp filter
catch 
    disp(['License Error: ', sfile]);
    error('stop execution')
end

% apply STA/LTA for P wave picking

%  p_win1 = 10; % P wave window start time (relative to theoritical arrival)
%  p_win2 = 5;  % P wave window end time

% widow using group velocity 
V_pn = deg2km(dist)/t_Pn;
V_pn1 = V_pn - 1.5 ; V_pn2 = V_pn + 6.5;
t_pn1 = deg2km(dist)/V_pn1;  
t_pn2 = deg2km(dist)/V_pn2;  

% for colse events (*win2 is in the slower end )
if dist < 0.5
    p_win1 = (t_pn1 - t_Pn) + ((3/dist) - (2/dist));   % right window limit 
    p_win2 = (t_Pn - t_pn2) + ((6/dist) - (2.5/dist)); % left window limit
else
    p_win1 = (t_pn1 - t_Pn); 
    p_win2 = (t_Pn - t_pn2);
end

% if begin window is smaller than sac start widow
if t_Pn < p_win2
    p_win2 = t_Pn;
    sac_p = sac_filt(ceil(n_tpn-round(p_win2*fs))+1:round(n_tpn+round(p_win1*fs)));
else
    sac_p = sac_filt(ceil(n_tpn-round(p_win2*fs)):round(n_tpn+round(p_win1*fs)));
end

% implement sequential triggering
thd = 2.5 : -0.05 : 1.8;
trig = zeros(1,length(thd));
for i = 1 : length (thd)
    trig(i) = pickps(sac_p,fs,1,10,thd(i));
end

trig = rmmissing(trig);      % remove the nan values
if isempty(trig)
    trig = NaN;
end

%tp = median(trig); % take the median value

% control the outlier
tpdiff = max(trig) - min(trig);
if tpdiff < 4.5 + dist*0.5
    tp = min(trig);
else
    tp = median(trig);
end

% convert pick to the absolute time scale
if isnan(tp)
    tp_picked = NaN;
else
    tp_picked = (t_Pn - p_win2) + tp;
end

% get end window
if ~isnan(tp_picked)
    Vp_picked = (deg2km(dist)/tp_picked);
    Vp_w2 = Vp_picked - 1.5;     % taking 1km/sec slower 
    tp_picked_2=deg2km(dist)/Vp_w2;
else
    tp_picked_2 = NaN;
end
coment=' ';

% check the S wave arrival time from P wave
% assuming vs ~ 3.5km/s
if ~isnan(tp_picked)
    %ts_pre = deg2km(dist) / 3.5;
    %t_Sn_pre = tp_picked + (t_Sn - t_Pn);
    if isnan(tp_picked)
        t_Sn_pre = t_Sn;
    else
        t_Sn_pre = tp_picked + (t_Sn - t_Pn);
    end
    
    if abs(t_Sn - t_Sn_pre) > (8 + dist*2) || Vp_picked < 4.5 || Vp_picked > 10 
        disp(['Event misassocation of wrong P wave arrival picked'])
        coment='P/S';
    else
        coment = 'ok';
    end
end


%% For S wave  
V_sn = deg2km(dist)/t_Sn;                 % prediceted S wave group velocity
V_sn1 = V_sn - 0.8 ; V_sn2 = V_sn + 0.6;  % velocity window
t_sn1 = deg2km(dist)/V_sn1; t_sn2 = deg2km(dist)/V_sn2;  % time window 

%s_win1 = (t_sn1 - t_Sn) + 2; s_win2 = t_Sn - t_sn2;
%disp(s_win1 + s_win2)

%disp(s_win2)

% for colse events
if dist < 0.7
    s_win1 = (t_sn1 - t_Sn) + ((5/dist) - (2/dist));     % right widow lim
    s_win2 = (t_Sn - t_sn2) + ((10/dist) - (2.5/dist));  % left window lim
else
    s_win1 = (t_sn1 - t_Sn); 
    s_win2 = (t_Sn - t_sn2);
end

%disp(s_win1 + s_win2)
%disp(s_win2)

% when S window exceeds the P arrivals
if t_Sn - t_Pn < s_win2
    s_win2 = ceil((t_Sn - t_Pn)*0.9);
    s_win1 = (15- s_win2);
end

if round(n_tsn+round(s_win1*fs)) > length(sac_filt)
    sac_s = sac_filt(ceil(n_tsn-(round(s_win2*fs))):end);
else
    sac_s=sac_filt(ceil(n_tsn-(round(s_win2*fs))):round(n_tsn+round(s_win1*fs)));
end



% implement sequential triggering
thds = 2.5 : -0.025 : 1.5;
trigs = zeros(1,length(thds));
for j = 1 : length (thds)
    trigs(j) = pickps(sac_s,fs,1,5,thds(j));
end


trigs = rmmissing(trigs);      % remove the nan values
if isempty(trigs)
    trigs = NaN;
end

% control the outlier
tsdiff = max(trigs) - min(trigs);
if tsdiff < 4 + dist*0.3
    ts = min(trigs);
else
    ts = median(trigs);
end

% use a wide filter bandpass Butterworth filter if not properly triggered
ts_diff = abs((t_Sn - s_win2) + ts - t_Sn);
ts_pick1 = (t_Sn - s_win2) + ts;

% check if ts is very close to tp
tstp_diff = 100;
if ~isnan(tp)
    tp_picked = (t_Pn - p_win2) + tp;
    tstp_diff = ts_pick1 - tp_picked;
end

if isnan(ts) || ts_diff > 3.0 + dist*0.3 || tstp_diff < 5 + dist*0.3
    c2 = 5;
    [b, a] = butter(n, [c1, c2]/(fs/2), 'bandpass');
    sac_filt2 = filtfilt(b, a, sac2); % apply zero-phase bp filter
    
    if round(n_tsn+round(s_win1*fs)) > length(sac_filt)
        sac_s = sac_filt2(ceil(n_tsn-(round(s_win2*fs))):end);
        else
        sac_s=sac_filt2(ceil(n_tsn-(round(s_win2*fs))):round(n_tsn+round(s_win1*fs)));
    end

    %sac_s=sac_filt2(ceil(n_tsn-(round(s_win2*fs))):round(n_tsn+round(s_win1*fs)));
    
    % implement sequential triggering
    thds = 2.5 : -0.025 : 1.5;
    trigs2 = zeros(1,length(thds));
    for j = 1 : length (thds)
        trigs2(j) = pickps(sac_s,fs,1,5,thds(j));
    end

    trigs2 = rmmissing(trigs2);      % remove the nan values
    if isempty(trigs2)
        trigs2 = NaN;
    end
    
    if ~isnan(trigs2)

        tsdiff = max(trigs2) - min(trigs2);
        if tsdiff < 3 + dist *0.5 && tstp_diff > 5 + dist*0.5
            ts2 = min(trigs2);
        elseif tstp_diff < 5 + dist*0.5
            ts2 = max(trigs);
        else
            ts2 = median(trigs2);
        end

        ts_pick2 = (t_Sn - s_win2) + ts2;
        ts_pick1 = (t_Sn - s_win2) + ts;

        if dist < 0.2
            ts = ts2;
        else
            diff1 = abs((ts_pick1 - t_Sn));
            diff2 = abs((ts_pick2 - t_Sn));
                if diff2 <= diff1 || isnan(diff1) || tstp_diff < 100
                    ts = ts2;
                end
        end
    end
end

if isnan(ts)
    ts_picked = NaN;
    
else
    ts_picked = (t_Sn - s_win2) + ts;
end

%% condition when S wave is not trigged but p traveltime is accurate 

res_p = abs(tp_picked - t_Pn);
if isnan(ts_picked) && res_p < 2 + dist*0.2
    ts_picked = t_Sn + (tp_picked - t_Pn);

end

if ~isnan(ts_picked)
    Vs_picked = (deg2km(dist)/ts_picked);
    Vs_w2 = Vs_picked - 0.9;        % group velocity of the end window
    ts_picked_2=deg2km(dist)/Vs_w2;
else
    ts_picked_2 = NaN;
    coment ='Nopick';
end

% condition when no data point  
if length(sac_filt) < ceil(fs*ts_picked)
    ts_picked = NaN;
    coment ='Nodata';
end
if length(sac_filt) < ceil(fs*ts_picked_2)
    ts_picked_2 = NaN;
    coment ='Nodata';
end

if ~isnan(ts_picked)
    %ts_pre = deg2km(dist) / 3.5;
    %t_Sn_pre = tp_picked + (t_Sn - t_Pn);
    if isnan(tp_picked)
        t_Sn_pre = t_Sn;
    else
        t_Sn_pre = tp_picked + (t_Sn - t_Pn);
    end
   
    if abs(t_Sn - t_Sn_pre) > (8 + dist*2) || Vs_picked < 2.4 || Vs_picked > 4.0 
        disp(['Event misassocation of wrong S wave arrival picked'])
        coment='P/S';
    else
        coment ='ok';
    end
end

%%  plot the pick window
fig1 = figure('position',[200, 650, 500, 280]);
plot(0:1/fs:(length(sac_p)-1)/fs,sac_p)
hold on
if ~isnan(tp)
    xline(tp,'LineWidth',2,'Color','g')
end
title('P - wave Pick', 'FontSize',14, 'FontWeight','bold')

fig2 = figure('position',[1000, 650, 500, 280]);
plot(0:1/fs:(length(sac_s)-1)/fs,sac_s)
hold on
if ~isnan(ts)
    xline(ts,'LineWidth',2,'Color','g')
    hold on
end
title('S - wave Pick', 'FontSize',14, 'FontWeight','bold')


%% plot the sac file

figure(3)
% filtered seismogram
subplot(2,1,1) 
plot(t,sac.d)
if ceil(t_Pn -20) < 0 
    xlim1 = 0;
    xlim2 = t_Sn + 60;
else
    xlim1 = ceil(t_Pn -20);
    xlim2 = t_Sn + 80;
end

xlim([xlim1 xlim2])

hold on
% plot predicted arrvial time
if ~isnan(t_Pn)
    xline(t_Pn,'LineWidth',1.5,'Color','r')
    hold on
end
if ~isnan(t_Sn)
    xline(t_Sn, 'LineWidth',1.5,'Color','r')
end
hold on

% plot picked arrvial time
if exist('T1', 'var')
    if ~isnan(T1) && ~isnan(T2)
        xline(T1,'LineWidth',1.5,'Color',[0 0 0])
        hold on
        xline(T2, 'LineWidth',1.5,'Color',[0 0 0])
        hold on
    end
end

tstring ='Original Waveform';
title(tstring, 'FontSize',12,'FontWeight','bold');
%legend('','Predicted Arrival', '', 'Picked Arrival')

subplot(2,1,2)
plot(t,sac_filt)
%xlim([sac.HEADER.B sac.HEADER.E])

if ceil(t_Pn -20) < 0 
    xlim1 = 0;
    xlim2 = t_Sn + 60;
else
    xlim1 = ceil(t_Pn -20);
    xlim2 = t_Sn + 100;
end

xlim([xlim1 xlim2])

hold on

% plot the predicted arrival
if ~isnan(t_Pn)
    xline(t_Pn,'--','LineWidth',1.5,'Color','r')
    hold on
end
if ~isnan(t_Sn)
    xline(t_Sn,'--','LineWidth',1.5,'Color','r')
    hold on
end

% plot the autopick
if ~isnan(tp_picked)
    xline(tp_picked,'LineWidth',1.5,'Color','g')
    hold on
    xline(tp_picked_2,'--','LineWidth',1.5,'Color','g')
    hold on
end

if ~isnan(ts_picked) 
    xline(ts_picked,'LineWidth',1.5,'Color','g')
    hold on
end
if ~isnan(ts_picked_2) 
    xline(ts_picked_2,'--','LineWidth',1.5,'Color','g')
    hold on
end

% plot picked arrvial time
if exist('T1', 'var')
    xline(T1,'--','LineWidth',1.5,'Color',[0 0 0])
    hold on
    xline(T2,'--', 'LineWidth',1.5,'Color',[0 0 0])
end

% plot end time window
tstring ='Filtered Waveform';
xlabel('time (s)','FontSize',12,'FontWeight','bold')
%xpos=0.9; ypos = 0.9;
% Add text to the top right corner
%annotation('textbox', [xpos, ypos, 0.1, 0.1], 'string', ...
%    tstring, 'EdgeColor', 'none', 'HorizontalAlignment', 'right', ...
%    'VerticalAlignment', 'top');
title(tstring, 'FontSize',12,'FontWeight','bold');
%legend('','Predicted Arrival', '', '','', 'Auto-picked Arrival', '', 'Picked Arrival')
tit=['Event : ', sfile(1:20), '    STN:   ', sfile(27:31), '  C: ', coment];
sgtitle(tit, 'FontSize',12,'FontWeight','bold')

set(gcf, 'Position',  [500, 100, 800, 480]);
filename=[sfile(1:31) '.png'];
%exportgraphics(gcf,filename,'Resolution',300) % expg does not work 
saveas(gcf, filename);
end
