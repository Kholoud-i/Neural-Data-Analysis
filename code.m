lfp = LoadBinary([prefix base_name '\' base_name '.eeg'],1); % getslfpchannel #1


sr = 1.25e3;


bpFilt = designfilt('bandpassfir','FilterOrder',20, …
'CutoffFrequency1',100,'CutoffFrequency2',250, …
'SampleRate',sr);
lfp_filt = filtfilt(bpFilt,lfp);


%% ripple detection
clc


params.zscore_cutoff = 2;
params.minimum_duration = 15; % ms
params.minimum_interval = 15; % ms


ripples = ripple_detection(lfp_filt,sr,params);


%% plot ripples
% grey: original lfp
% blue: filtered lfp
% thickened blue: ripple
figure,
for i = 1:20


subplot(4,5,i)
stime = ripples(i).start_time;
etime = ripples(i).end_time;
% 50ms before and after the ripple
tt = 50/1000;
y1 = lfp_filt(floor(stime*sr-tt*sr+1):floor(etime*sr+tt*sr));
t1 = linspace(0,length(y1)/sr,length(y1));
plot(t1,y1,'color',[0, 0.4470, 0.7410,.5])
hold on
% plot original lfp
y2 = lfp(floor(stime*sr-tt*sr+1):floor(etime*sr+tt*sr));
plot(t1,y2,'color',[0 0 0 .2])
hold on
% plot the ripple itself
y3 = ripples(i).data;
t3 = linspace(tt,tt+etime-stime,length(y3));
plot(t3,y3,'color',[0, 0.4470, 0.7410],'linewidth',1)
xlim([0 max(t1)])
title(num2str(i))
end




set(gcf,'position',[100,100,900,600])
han=axes(gcf,'visible','off');
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'Time [Hz]');
% ylabel(han,'Power (dB)');


%% hilbert transform
figure,
for i = 1:20
subplot(4,5,i)
y = hilbert(ripples(i).data);
[pxx,f] = pwelch(y,[],0,1:300,sr,'power');
p = pow2db(pxx);
% [p,f] = pspectrum(y,sr,'Power','FrequencyLimits',[0 300]);
plot(f,p)
title(num2str(i))
end
han=axes(gcf,'visible','off');
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'Frequency [Hz]');
ylabel(han,'Power [dB]');


function ripples = ripple_detection(lfp_filt,sr,params)


% envelop
[up,~] = envelope(lfp_filt);


% z score
zscore = (up.^2 - mean(up.^2))/std(up.^2);


% parameters for ripple detection
min_duration = params.minimum_duration/1000; % second
cutoff = params.zscore_cutoff;
min_interval = params.minimum_interval/1000; % second




% detect ripples
sig_index = find(zscore>=cutoff);
i = 1;n=1;
ripples = struct();
while i<=length(sig_index)
rip_start = sig_index(i);
zscore_rest = zscore(rip_start:end);
rip_end = rip_start + find(zscore_rest<cutoff,1)-1-1;
if rip_end-rip_start+1 > min_duration*sr
ripples(n).start_time = rip_start/sr;
ripples(n).end_time = rip_end/sr;
n = n+1;
end
i = find(sig_index>rip_end,1);
end


% merge ripples
i = 1;
while i < length(ripples)
interval = ripples(i+1).start_time - ripples(i).end_time;
if interval > min_interval
i = i+1;
else
ripples(i).end_time = ripples(i+1).end_time;
ripples(i+1) = [];
end
end


%
for i = 1:length(ripples)
ripples(i).duration = ripples(i).end_time - ripples(i).start_time;ripples(i).data =
lfp_filt(floor(ripples(i).start_time*sr):floor(ripples(i).end_time*sr));
end


