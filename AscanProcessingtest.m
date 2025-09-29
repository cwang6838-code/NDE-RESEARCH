clc; clear; close all;
%finds LS_PANEL_11_3_PE anywhere on your computer and sets a path to it
scriptFullPath = mfilename('fullpath'); 
repoRoot = fileparts(scriptFullPath);
dataFolder = fullfile(repoRoot);
%uses found path to extract test data
filepath = (dataFolder);
radialPos = ['t', 'r', 'u', 's', "c1", 'n', "d1", 'o', 'p', 'q', 'f', 'g', 'h', 'i', 'j', 'a', 'b', 'c', 'd', 'e', 'w', 'a1', 'e1', 'e2', 'e3', 'e4', 'e5', 'l', 'v', 'k', 'x', 'y', 'm', 'z', 'b1', 'ref1', 'ref2', 'ref3', 'ref4', 'ref5', 'ref6', 'ref7', 'ref0'];
tests = ["ref1", 'ref2', 'ref3', 'ref4', 'ref5', 'ref6', 'ref7', 'a', 'a1', 'b', 'b1', 'b2', 'c', 'c1', 'd', 'd1', 'e', 'e1', 'e2', 'e3', 'e4', 'e5', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'];
% 1=time, 2=Ch1 (trigger), 3=Ch2 (data), 4=Ref (=ref0)
%% Loading raw csv
filenames = strings([1, 43]);
numTests = 43;
testMatrix = zeros(4999, 43);
for ii = 1:length(radialPos)
    filename = filepath + "\" + radialPos(ii) + ".csv";
    filenames(1, ii) = filename;
end
for kk = 1:numTests-1
    testSignal = readmatrix(filenames(kk));
    testMatrix(:, kk) = testSignal(3:end-1, 3);
end
ref0signal = readmatrix(filenames(36));
testMatrix(:, 43) = ref0signal(3:end-1, 4);

Ns = size(testMatrix,1);
tvec = ref0signal(3:end-1,1);
dt = tvec(10)-tvec(9); fs = 1/dt;
time = 0:dt:(Ns-1)*dt;

trigger = ref0signal(3:end-1,2);
%% Loading saved .mat file
% References: ref0-7 (ref0=saveds as Ref1 on Oscilloscope), b1 (new!)
figure(1);
plot(time*1e6,trigger);grid on;hold on;
for g = 35:42
plot(time*1e6,testMatrix(:,g));grid on;hold on;
end
ylim([-0.5 0.5])
title('Reference');xlabel('time [\mus]')

% % References: b2 (new!, old b1) on no paint edge
% figure(7);
% plot(time*1e6,trigger);grid on;hold on; 
% plot(time*1e6,b2(:,1));grid on;hold on; 
% ylim([-0.5 0.5])
% title('no paint edge');xlabel('time [\mus]')

% Pristine close to paint chip: w, a1, l, e1-5,v,k,x,y,m
% in ascanprocessing3 a1 has time not multiplied by 1e6, not sure if its a
% mistake or not but left it out here
figure(2);
plot(time*1e6,trigger);grid on;hold on;
for h = 21:33
    plot(time*1e6,testMatrix(:,h));grid on;hold on;
end 
ylim([-0.5 0.5])
title('Pristine (close to edge)');xlabel('time [\mus]')

% No paint close to edge: a, b, c, d, e 
figure(3);
plot(time*1e6,trigger);grid on;hold on; 
plot(time*1e6,testMatrix(:,43),'k');grid on;hold on; 
for v = 16:20
    plot(time*1e6,testMatrix(:,v));grid on;hold on;
end
ylim([-0.5 0.5])
title('No paint (further from LS center)');xlabel('time [\mus]')

% No paint further from edge: f, g, h, i, j
figure(4);
plot(time*1e6,trigger);grid on;hold on; 
plot(time*1e6,testMatrix(:,43),'k');grid on;hold on; 
for j = 11:15
    plot(time*1e6,testMatrix(:,j));grid on;hold on;
end
ylim([-0.5 0.5])
title('No paint (closer to LS center)');xlabel('time [\mus]')

% Damage (black area): n, d1,o,p,q
figure(5)
plot(time*1e6,trigger);grid on;hold on; 
plot(time*1e6,testMatrix(:,43),'k');grid on;hold on; 
for o = 6:10
    plot(time*1e6,testMatrix(:,o));grid on;hold on;
end
ylim([-0.5 0.5])
title('Severe Damage');xlabel('time [\mus]')

% Damage (delam and severe): r,s,t,u,c1
figure(6)
plot(time*1e6,trigger);grid on;hold on; 
plot(time*1e6,testMatrix(:,43),'k');grid on;hold on; 
for i = 1:5
    plot(time*1e6,testMatrix(:,i));grid on;hold on;
end
ylim([-0.5 0.5])
title('Severe Damage/Delam');xlabel('time [\mus]')

%% Calculate TOF and AMP of max
jj = 37; % ref0=1, ref1=2, ref7=8; 
% i=9,w; 10-l, 11-e1, 12-e2. 13-e3, 14-e4, 15-e5, 16-k, 17-x
% i=18-a, 19-b, 20-c, 21-d, 22-e
% i=23-f, 24-g, 25-h, 26-ii, 27-j
% i=28-d1, 29-n, 30-o, 31-p, 32-q
% i=33-c1, 34-r, 35-s, 36-t, 37-u
signal = testMatrix(:,3);
[ampl(jj) ind] = max(abs(signal(1201:1701)),[],1); % max bw 2.4mus-3.4mus
tof(jj) = time(1201+ind)*1e6; % mus


% Pristine close to paint chip: w, a1, l, e1-5,v,k,x,y,m
%pmatrix = [t, r, u, s, c1, n, d1, o, p, q, f, g, h, i, j, a, b, c, d, e, w, a1, e1, e2, e3, e4, e5, l, v, k, x, y, m, z, b1, ref0, ref1, ref2, ref3, ref4, ref5, ref6, ref7];
% pmatrix is the 2d matrix with x being position and y being voltage in
% order of the rms graph before timegating
tgatemat = testMatrix;
tgate1 = testMatrix;
tgate1(1: 1201,:) = 0;tgate1(1701:end,:) = 0;
tgate2 = testMatrix;
tgate2(1:1771, :) = 0;tgate2(2391:end,:) = 0;
hann1 = tgatemat(1202:1700, :) .* hann(499);
hann2 = tgatemat(1772:2390, :) .* hann(619);
testZeros = zeros(4999, 43);
testZeros2 = zeros(4999, 43);
testZeros(1202:1700, :) = hann1;
testZeros2(1772:2390, :) = hann2;
bWall = testZeros;
ref2 = testZeros2;

% tgatemat is the 2d pmatrix with the timegating for the rms graph
RMSvec1 = rms(tgate1);
RMSvec2 = rms(tgate2);
labels = {'t','r','u','s','c1','n','d1','o','p','q','f','g','h','i','j','a','b','c','d','e','w','a1','e1','e2','e3','e4','e5','l','v','k','x','y','m','z','b1','ref1','ref2','ref3','ref4','ref5','ref6','ref7', 'ref0'};

figure(10)
subplot(2, 1, 1);
plot(1:43, RMSvec1, "LineWidth", 2); hold on; grid on;
text(1:43, RMSvec1, labels,'VerticalAlignment','bottom');
xl = xline(20.5,'-.','Paint Edge','DisplayName','Paint Edge'); hold on;
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
title('BW RMS');
subplot(2, 1, 2);
plot(1:43, RMSvec2, "LineWidth", 2); hold on; grid on;
text(1:43, RMSvec2, labels,'VerticalAlignment','bottom');
xl = xline(20.5,'-.','Paint Edge','DisplayName','Paint Edge'); hold on;
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
title('REF2 RMS');

% RMS POSITIONS PROPAGATE RADIALLY: Left side = plate center, right side = far from center
% Pos 1-5 = severe / delam, 6-10 = black area, 11-15 = no paint further
% from edge, 16-20 = no paint close to edge, 21-34 = pristine close to paint chip, 35-42 = reference locations 

% Frequency Domain Analysis
%freq = (fs.*((0:Ns/2)*(fs/(Ns)))/1000)';
freq = linspace(-fs/2, fs/2, Ns);
bWallF = fftshift(abs(fft(bWall, Ns)./(Ns)));
bWallF2 = fftshift(abs(fft(ref2, Ns)./(Ns)));

figure(12)
plot(freq/1e6, bWallF, 'LineWidth', 1.2); hold on; grid on;
xlabel('Frequency (MHz)');
ylabel('Amplitude');
title('BW Two-Sided Amplitude Spectrum (Normalized)');
xlim([0, 12]);



figure(13)
plot(freq/1e6, bWallF2, 'LineWidth', 1.2); hold on; grid on;
xlabel('Frequency (MHz)');
ylabel('Amplitude');
title('2nd Ref Two-Sided Amplitude Spectrum (Normalized)');
xlim([0, 12]);


maxAmp1st = zeros(1, length(radialPos));  % Vector of frequencies corresponding to max Ampl.
maxAmp2nd = maxAmp1st;
for zz = 1:length(radialPos)
    [maxFw, indexFw] = max(bWallF);
    [maxBw, indexBw] = max(bWallF2);
end
    maxAmp1st(:) = freq(indexFw(:));
    maxAmp2nd(:) = freq(indexBw(:));

intFw = zeros(1, length(radialPos));
intBw = intFw;
for zz = 1:length(radialPos)
    intFw(zz) = trapz(freq, bWallF(:, zz).^2);
    intBw(zz) = trapz(freq,bWallF2(:, zz).^2);
end

figure(15);
subplot(2, 1, 1);
plot(1:43, intFw, 'LineWidth', 1.2); grid on;
xlabel('Radial Position');
ylabel('Power (V^2)');
title('BW Total Power');
xl = xline(20.5,'-.','Paint Edge','DisplayName','Paint Edge'); hold on;
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
subplot(2, 1, 2);
plot(1:43, intBw, 'LineWidth', 1.2); grid on;
xlabel('Radial Position');
ylabel('Power (V^2)');
title('2nd Ref Total Power');
xl = xline(20.5,'-.','Paint Edge','DisplayName','Paint Edge'); hold on;
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';

figure(16);
subplot(2, 1, 1);
plot(1:43, maxAmp1st, 'LineWidth', 1.2); grid on;
xlabel('Radial Position');
ylabel('Frequency (MHz)');
title('BW Freq. of Max Ampl.');
xl = xline(20.5,'-.','Paint Edge','DisplayName','Paint Edge'); hold on;
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
subplot(2, 1, 2);
plot(1:43, maxAmp2nd, 'LineWidth', 1.2); grid on;
xlabel('Radial Position');
ylabel('Frequency (MHz)');
title('2nd Ref Freq. of Max Ampl.');
xl = xline(20.5,'-.','Paint Edge','DisplayName','Paint Edge'); hold on;
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';


