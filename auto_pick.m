%function [tp1, ts1] = auto_pick(sfile)
%clear;
close all
sfile='2023.132.08.33.21.0800.YB.EF040.00.DPZ.D.SAC';
sac=rdsac(sfile);
dist=sac.HEADER.GCARC;
fs=1./sac.HEADER.DELTA;

edp=sac.HEADER.EVDP;

T1=sac.HEADER.T1; T2 = sac.HEADER.T2;   % picked arrvial times
 
%t_Pn = get_ttime(dist,edp,'Pn');
%t_Sn = get_ttime(dist,edp,'Sn');

t_Pn = get_ttime(dist,edp,'P');
t_Sn = get_ttime(dist,edp,'S');

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
c1 = 0.05; % Lower corner frequency in Hz
c2 = 2.0; % Upper corner frequency in Hz
n = 3; % Filter order (number of poles)


% Design a bandpass Butterworth filter
[b, a] = butter(n, [c1, c2]/(fs/2), 'bandpass');
sac_filt = filtfilt(b, a, sac2); % apply zero-phase bp filter


% apply STA/LTA for P wave picking

p_win1 = 10; % P wave window start time (relative to theoritical arrival)
p_win2 = 5;  % P wave window end time

% widow using group velocity 



% if begin window is smaller than start widow
if t_Pn < p_win1
    p_win1 = t_Pn;
    sac_p = sac_filt(ceil(n_tpn-round(p_win1*fs))+1: round(n_tpn+round(p_win2*fs)));
else
    sac_p = sac_filt(ceil(n_tpn-round(p_win1*fs)): round(n_tpn+round(p_win2*fs)));
end

tp = pickps(sac_p,fs,1,10,2.5);


% lower the threshold value if not triggered
if isnan(tp)
    tp = pickps(sac_p,fs,1,10,2.0);
end


% for S wave

s_win1 = 5;  % S wave window start time
s_win2 = 10; % S wave window end time

% when S window exceeds the P arrivals
if t_Sn - t_Pn < s_win1
    s_win1 = ceil((t_Sn - t_Pn)*0.6);
end

sac_s = sac_filt(ceil(n_tsn-round(s_win1*fs)): round(n_tsn+round(s_win2*fs)));
ts = pickps(sac_s,fs,0.5,10,2.0);

% lower the threshold value if not triggered
if isnan(ts)
    ts = pickps(sac_s,fs,0.5,10,1.5);
end



% convert pick to the absolute time scale
tp1 = (t_Pn - p_win1) + tp;
ts1 = (t_Sn - s_win1) + ts;

%%
fig = figure('position',[200, 650, 500, 280]);
plot(0:1/fs:(length(sac_p)-1)/fs,sac_p)
hold on
xline(tp,'LineWidth',2,'Color','g')
title('P - wave Pick', 'FontSize',14, 'FontWeight','bold')
fig = figure('position',[1000, 650, 500, 280]);
plot(0:1/fs:(length(sac_s)-1)/fs,sac_s)
hold on
xline(ts,'LineWidth',2,'Color','g')
title('S - wave Pick', 'FontSize',14, 'FontWeight','bold')


%% plot the sac file

figure(3)
subplot(2,1,1)
plot(t,sac.d)
xlim([ceil(t_Pn -10) ceil(t_Sn + 60)])
hold on
% plot predicted arrvial time
xline(t_Pn,'LineWidth',1.5,'Color','r')
hold on
xline(t_Sn, 'LineWidth',1.5,'Color','r')
hold on
% plot picked arrvial time

xline(T1,'LineWidth',1.5,'Color',[0 0 0])
hold on
xline(T2, 'LineWidth',1.5,'Color',[0 0 0])

tstring ='Original Waveform';
title(tstring, 'FontSize',12,'FontWeight','bold');
legend('','Predicted Arrival', '', 'Picked Arrival')

subplot(2,1,2)
plot(t,sac_filt)
%xlim([sac.HEADER.B sac.HEADER.E])
xlim([ceil(t_Pn -10) ceil(t_Sn + 60)])
hold on
xline(t_Pn,'LineWidth',1.5,'Color','r')
hold on
xline(t_Sn,'LineWidth',1.5,'Color','r')
hold on
% plot the autopick
xline(tp1,'LineWidth',1.5,'Color','g')
hold on
xline(ts1,'LineWidth',1.5,'Color','g')
tstring ='Filtered Waveform';
xlabel('time (s)','FontSize',12,'FontWeight','bold')
%xpos=0.9; ypos = 0.9;
% Add text to the top right corner
%annotation('textbox', [xpos, ypos, 0.1, 0.1], 'string', ...
%    tstring, 'EdgeColor', 'none', 'HorizontalAlignment', 'right', ...
%    'VerticalAlignment', 'top');
title(tstring, 'FontSize',12,'FontWeight','bold');
legend('','Predicted Arrival', '', 'Auto-picked Arrival')
tit=['Event : ', sfile(1:20), '    Stn : ', sfile(27:31)];
sgtitle(tit, 'FontSize',12,'FontWeight','bold')

set(gcf, 'Position',  [500, 150, 800, 400])
filename=[sfile(1:31) '.png'];
exportgraphics(gcf,filename,'Resolution',300)

%end