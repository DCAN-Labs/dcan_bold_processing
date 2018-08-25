function PMUstructmat = TransformPMUOutput(PMUdir,varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
respfreqmax = 25;
pulsfreqmax = 200;
if isempty(varargin) == 0
    for i = 1:size(varargin,2)
        if isstruct(varargin{i}) == 0
            switch(varargin{i})
                case('MaxRespFreq')
                    respfreqmax = varargin{i+1};
                case('MaxPulsFreq')
                    pulsfreqmax = varargin{i+1};
            end
        end
    end
end
%identify files to load
volfile = dir(strcat(PMUdir,'/*VolumeSliceAcq*.txt'));
respfile = dir(strcat(PMUdir,'/*RESP.txt'));
pulsfile = dir(strcat(PMUdir,'/*PULS.txt'));
volfile = [volfile.folder filesep volfile.name];
respfile = [respfile.folder filesep respfile.name];
pulsfile = [pulsfile.folder filesep pulsfile.name];
%get data from files
fid = fopen(volfile);
header = textscan(fid,'%s%s%s',1);
raw_data = textscan(fid,'%d%d%d\n');
volume_number = raw_data{1};
slice_number = double(raw_data{2});
data_acq_time = double(raw_data{3});
fclose('all');
clear header
clear raw_data
clear fid
%dump the first volume because the timing is screwed up and
%immeasurable -- determine the slice sampling rate
data_acq_time = data_acq_time(volume_number > 0);
slice_number = slice_number(volume_number > 0);
volume_number = volume_number(volume_number > 0);
slice_rate = 1/(mode(diff(data_acq_time(1:2:end)))/0.4)*1000; %interleaved acquisition means we have to take every other to calculate the slice rate
disp(strcat('slice acquisition rate identified --',num2str(slice_rate),' Hz'))
if data_acq_time(2) > data_acq_time(1)
    data_acq_order = 1;
else
    data_acq_order = 2;
end
fid = fopen(respfile);
header = textscan(fid,'%s%s%s',1);
raw_data = textscan(fid,'%d%d%d\n');
time_resp = double(raw_data{1});
signal_resp = double(raw_data{2});
fclose('all');
clear header
clear raw_data
clear fid
%match the start of the bold data to the respiration acquisition --
%chuck the noise and determine the respiration sampling rate
if isempty(find(time_resp == data_acq_time(data_acq_order)))
    phase_shift = '';
    if (max(data_acq_time) > max(time_resp))
        disp('acquisition of EPI and respiration are out of phase!! Shifting respiration FORWARD in time')
        phase_shift = 'forward';
    else
        disp('acquisition of EPI and respiration are out of phase!! Shifting respiration BACKWARD in time')
        phase_shift = 'backward';
    end
    time_resp_point = max(find(time_resp < data_acq_time(data_acq_order)));
    switch(phase_shift)
        case('forward')
            time_resp = time_resp(time_resp_point:end);
            signal_resp = signal_resp(time_resp_point:end);
            time_resp_diff = data_acq_time(data_acq_order) - time_resp(1);
            time_resp = time_resp + time_resp_diff;
        case('backward')
            time_resp = time_resp(time_resp_point+1:end);
            signal_resp = signal_resp(time_resp_point+1:end);
            time_resp_diff = time_resp(1) - data_acq_time(data_acq_order);
            time_resp = time_resp + time_resp_diff;
    end
else
    signal_resp = signal_resp(find(time_resp == data_acq_time(data_acq_order)):end);
    time_resp = time_resp(find(time_resp == data_acq_time(data_acq_order)):end);
end
resp_rate = 1/(mode(diff(time_resp))/0.4)*1000;
disp(strcat('respiration rate identified --',num2str(resp_rate),' Hz'))
fid = fopen(pulsfile);
header = textscan(fid,'%s%s%s',1);
raw_data = textscan(fid,'%d%d%d\n');
time_puls = double(raw_data{1});
signal_puls = double(raw_data{2});
fclose('all');
clear header
clear raw_data
clear fid
%match the start of the bold data to the respiration acquisition --
%chuck the noise and determine the respiration sampling rate
if isempty(find(time_puls == data_acq_time(data_acq_order)))
    phase_shift = '';
    if (max(data_acq_time) > max(time_puls))
        disp('acquisition of EPI and heart rate are out of phase!! Shifting heart rate FORWARD in time')
        phase_shift = 'forward';
    else
        disp('acquisition of EPI and heart rate are out of phase!! Shifting heart rate BACKWARD in time')
        phase_shift = 'backward';
    end
    time_puls_point = max(find(time_puls < data_acq_time(data_acq_order)));
    switch(phase_shift)
        case('forward')
            time_puls = time_puls(time_puls_point:end);
            signal_puls = signal_puls(time_puls_point:end);
            time_puls_diff = data_acq_time(data_acq_order) - time_puls(1);
            time_puls = time_puls + time_puls_diff;
        case('backward')
            time_puls = time_puls(time_puls_point+1:end);
            signal_puls = signal_puls(time_puls_point+1:end);
            time_puls_diff = time_puls(1) - data_acq_time(data_acq_order);
            time_puls = time_puls + time_puls_diff;
    end
else
    signal_puls = signal_puls(find(time_puls == data_acq_time(data_acq_order)):end);
    time_puls = time_puls(find(time_puls == data_acq_time(data_acq_order)):end);
end
puls_rate = 1/(mode(diff(time_puls))/0.4)*1000;
disp(strcat('heart rate identified --',num2str(puls_rate),' Hz'))
% create struct matrix for data and rates
PMUstructmat.slice = [data_acq_time./0.4 slice_number];
PMUstructmat.respiration = [time_resp./0.4 signal_resp];
PMUstructmat.pulse = [time_puls./0.4 signal_puls];
PMUstructmat.resprate = resp_rate;
PMUstructmat.pulsrate = puls_rate;
PMUstructmat.respduration = (max(time_resp) - min(time_resp))/0.4;
PMUstructmat.pulsduration = (max(time_puls) - min(time_puls))/0.4;
% use FFTs to produce frequency-power plots
pulsFFT = fft(PMUstructmat.pulse(:,2));
respFFT= fft(PMUstructmat.respiration(:,2));
P2_puls = abs(pulsFFT/length(pulsFFT));
P2_resp = abs(respFFT/length(respFFT));
P1_puls = P2_puls(1:length(pulsFFT)/2 + 1);
P1_resp = P2_resp(1:length(respFFT)/2 + 1);
P1_puls(2:end-1) = 2*P1_puls(2:end-1);
P1_resp(2:end-1) = 2*P1_resp(2:end-1);
f_resp = PMUstructmat.resprate*(0:(length(respFFT)/2))/length(respFFT);
f_puls = PMUstructmat.pulsrate*(0:(length(pulsFFT)/2))/length(pulsFFT);
h = figure(1);
subplot(2,1,1);
plot(f_resp(2:end),P1_resp(2:end));
xlim([0 respfreqmax])
xlabel('frequency (Hz)')
ylabel('amplitude')
title('Frequency Power plot for respiration');
set(gca,'FontSize',16,'FontWeight','Bold','FontName','Arial');
subplot(2,1,2);
plot(f_puls(2:end),P1_puls(2:end));
xlim([0 pulsfreqmax])
xlabel('frequency (Hz)')
ylabel('amplitude')
title('Frequency Power plot for heart rate');
set(gca,'FontSize',16,'FontWeight','Bold','FontName','Arial');
PMUstructmat.respFFT = respFFT;
PMUstructmat.pulsFFT = pulsFFT;
PMUstructmat.f_resp = f_resp;
PMUstructmat.f_puls = f_puls;
saveas(h,strcat(PMUdir,'/freq_power.tif'));
save(strcat(PMUdir,'/PMUextracted.mat'),'PMUstructmat');
close all
%quit

