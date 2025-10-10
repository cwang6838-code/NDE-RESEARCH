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

testMatrix = testMatrix - mean(testMatrix);

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

figure(7)
plot(time*1e6, trigger); grid on; hold on;
for ll = 1:43
plot(time*1e6, testMatrix(:, ll));
end
ylim([-0.5 0.5])

%% Calculate TOF and AMP of max
jj = 37; % ref0=1, ref1=2, ref7=8; 
% i=9,w; 10-l, 11-e1, 12-e2. 13-e3, 14-e4, 15-e5, 16-k, 17-x
% i=18-a, 19-b, 20-c, 21-d, 22-e
% i=23-f, 24-g, 25-h, 26-ii, 27-j
% i=28-d1, 29-n, 30-o, 31-p, 32-q
% i=33-c1, 34-r, 35-s, 36-t, 37-u

% Pristine close to paint chip: w, a1, l, e1-5,v,k,x,y,m
%pmatrix = [t, r, u, s, c1, n, d1, o, p, q, f, g, h, i, j, a, b, c, d, e, w, a1, e1, e2, e3, e4, e5, l, v, k, x, y, m, z, b1, ref0, ref1, ref2, ref3, ref4, ref5, ref6, ref7];
% pmatrix is the 2d matrix with x being position and y being voltage in
% order of the rms graph before timegating

% tgatemat = testMatrix;
% tgate1 = testMatrix;
% tgate1(1: 1201,:) = 0;tgate1(1701:end,:) = 0;
% tgate2 = testMatrix;
% tgate2(1:1771, :) = 0;tgate2(2391:end,:) = 0;
% hann1 = tgatemat(1202:1700, :) .* hann(499);
% hann2 = tgatemat(1772:2390, :) .* hann(619);
% testZeros = zeros(4999, 43);
% testZeros2 = zeros(4999, 43);
% testZeros(1202:1700, :) = hann1;
% testZeros2(1772:2390, :) = hann2;

adaptiveTgate1 = testMatrix;
adaptiveTgate2 = testMatrix;

tg.boundleft1 = 2.3e-6 / dt + 1; % Left bound for MAX signal search (s)
tg.boundright1 = 3.5e-6 / dt + 1;

tg.left = 0.35e-6  / dt;
tg.right = 0.75e-6 / dt;
tgateamplitude = zeros(2,43);
maxIndex = zeros(1,43);

for zz = 1:numTests
    [Amax, imax] = max(abs(testMatrix(tg.boundleft1:tg.boundright1, zz)));
    imax = imax + tg.boundleft1;
    maxIndex(1,zz) = -1*imax;
    windowStart = round(imax - tg.left);
    windowEnd = round(imax + tg.right);
    hannAdapt = hann(round(windowEnd - windowStart + 1));
    adaptiveTgate1(1:windowStart-1, zz) = 0;
    adaptiveTgate1(windowEnd+1:end, zz) = 0;
    adaptiveTgate1(windowStart:windowEnd, zz) = adaptiveTgate1(windowStart:windowEnd, zz) .* hannAdapt;
    tgateamplitude(1,zz) = max(abs(adaptiveTgate1(:,zz)));
end

tg.boundleft2 = 3.5e-6 / dt + 1; % Left bound for MAX signal search (s)
tg.boundright2 = 5e-6 / dt + 1;

for zz = 1:numTests
    [Amax, imax] = max(abs(testMatrix(tg.boundleft2:tg.boundright2, zz)));
    imax = imax + tg.boundleft2;
    maxIndex(1,zz) = maxIndex(1,zz) + imax;
    windowStart = round(imax - tg.left);
    windowEnd = round(imax + tg.right);
    hannAdapt = hann(round(windowEnd - windowStart + 1));
    adaptiveTgate2(1:windowStart-1, zz) = 0;
    adaptiveTgate2(windowEnd+1:end, zz) = 0;
    adaptiveTgate2(windowStart:windowEnd, zz) = adaptiveTgate2(windowStart:windowEnd, zz) .* hannAdapt;
    tgateamplitude(2,zz) = max(abs(adaptiveTgate2(:,zz)));
end

adaptiveTgate1 = adaptiveTgate1 - mean(adaptiveTgate1);
adaptiveTgate2 = adaptiveTgate2 - mean(adaptiveTgate2);

figure(8);
plot(time*1e6, adaptiveTgate1(:, :));
grid on;
title('Time Gated Signals');

figure(9);
plot(time*1e6, adaptiveTgate2(:, :));
grid on;
title('Time Gated Signals');

% adaptiveTgate1 = adaptiveTgate1 - mean(adaptiveTgate1);
% adaptiveTgate2 = adaptiveTgate2 - mean(adaptiveTgate2);

% tgatemat is the 2d pmatrix with the timegating for the rms graph
RMSvec1 = rms(adaptiveTgate1);
RMSvec2 = rms(adaptiveTgate2);
labels = {'t','r','u','s','c1','n','d1','o','p','q','f','g','h','i','j','a','b','c','d','e','w','a1','e1','e2','e3','e4','e5','l','v','k','x','y','m','z','b1','ref1','ref2','ref3','ref4','ref5','ref6','ref7', 'ref0'};

figure(10)
subplot(2, 1, 1);
plot(1:43, RMSvec1, "LineWidth", 2); hold on; grid on;
text(1:43, RMSvec1, labels,'VerticalAlignment','bottom');
xl = xline(20.5,'-.','Paint Edge','DisplayName','Paint Edge'); hold on;
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
title('RMS (reflection 1)');
subplot(2, 1, 2);
plot(1:43, RMSvec2, "LineWidth", 2); hold on; grid on;
text(1:43, RMSvec2, labels,'VerticalAlignment','bottom');
xl = xline(20.5,'-.','Paint Edge','DisplayName','Paint Edge'); hold on;
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
title('RMS (reflection 2)');

% RMS POSITIONS PROPAGATE RADIALLY: Left side = plate center, right side = far from center
% Pos 1-5 = severe / delam, 6-10 = black area, 11-15 = no paint further
% from edge, 16-20 = no paint close to edge, 21-34 = pristine close to paint chip, 35-42 = reference locations 

% Frequency Domain Analysis
%freq = (fs.*((0:Ns/2)*(fs/(Ns)))/1000)';
freq = linspace(-fs/2, fs/2, Ns);
reflection1F = fft(adaptiveTgate1, Ns)./(Ns);
reflection2F = fft(adaptiveTgate2, Ns)./(Ns);

figure(12)
plot(freq/1e6, fftshift(abs(reflection1F)), 'LineWidth', 1.2); hold on; grid on;
xlabel('Frequency (MHz)');
ylabel('Amplitude');
title('Two-Sided Amplitude Spectrum (reflection 1)');
xlim([0, 12]);

figure(13)
plot(freq/1e6, fftshift(abs(reflection2F)), 'LineWidth', 1.2); hold on; grid on;
xlabel('Frequency (MHz)');
ylabel('Amplitude');
title('Two-Sided Amplitude Spectrum (reflection 2)');
xlim([0, 12]);

% Freq. of Max Amplitude
maxAmp1 = zeros(1, length(radialPos));  % Vector of frequencies corresponding to max Ampl.
maxAmp2 = maxAmp1;
for zz = 1:numTests
    [max1, index1] = max(abs(fftshift(reflection1F)));
    [max2, index2] = max(abs(fftshift(reflection2F)));
end
    maxAmp1(:) = freq(index1(:));
    maxAmp2(:) = freq(index2(:));

    % POWER
int1 = zeros(1, length(radialPos));
int2 = int1;
for zz = 1:length(radialPos)
    int1(zz) = trapz(freq, abs(fftshift(reflection1F(:, zz))).^2 / fs);
    int2(zz) = trapz(freq,abs(fftshift(reflection2F(:, zz))).^2 / fs);
end

figure(15);
subplot(2, 1, 1);
plot(1:43, int1, 'LineWidth', 1.2); grid on;
xlabel('Radial Position');
ylabel('Power (V^2)');
title('Total Power (reflection 1)');
xl = xline(20.5,'-.','Paint Edge','DisplayName','Paint Edge'); hold on;
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
subplot(2, 1, 2);
plot(1:43, int2, 'LineWidth', 1.2); grid on;
xlabel('Radial Position');
ylabel('Power (V^2)');
title('Total Power (reflection 2)');
xl = xline(20.5,'-.','Paint Edge','DisplayName','Paint Edge'); hold on;
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';

figure(16);
subplot(2, 1, 1);
plot(1:43, maxAmp1, 'LineWidth', 1.2); grid on;
xlabel('Radial Position');
ylabel('Frequency (MHz)');
title('Freq. of Max Ampl. (reflection 1)');
xl = xline(20.5,'-.','Paint Edge','DisplayName','Paint Edge'); hold on;
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
subplot(2, 1, 2);
plot(1:43, maxAmp2, 'LineWidth', 1.2); grid on;
xlabel('Radial Position');
ylabel('Frequency (MHz)');
title('Freq. of Max Ampl. (reflection 2)');
xl = xline(20.5,'-.','Paint Edge','DisplayName','Paint Edge'); hold on;
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';

figure(17);
subplot(2, 1, 1);
plot(time*1e6, abs(testMatrix(:, 1))); grid on; hold on;
plot(time*1e6, abs(testMatrix(:, 43))); grid on; hold on;
legend('(t)Severe/Delam', '(Ref0) Reference')
title('Non-timegated Severe and Reference');xlabel('time [\mus]')
subplot(2, 1, 2);
plot(time*1e6, abs(adaptiveTgate1(:, 1))); grid on; hold on;
plot(time*1e6, abs(adaptiveTgate2(:, 1))); grid on; hold on;
plot(time*1e6, abs(adaptiveTgate1(:, 43))); grid on; hold on;
plot(time*1e6, abs(adaptiveTgate2(:, 43))); grid on; hold on;
legend('(t)Severe/Delam First Reflection','(t)Severe/Delam Second Reflection', '(Ref0) Reference First Reflection', '(Ref0) Reference Second Reflection')
title('Timegated Severe and Reference');

figure(18);
%subplot(2, 1, 1);
plot(1:43,tgateamplitude(1,:),"Marker","o","LineStyle","-"); grid on; hold on;
plot(1:43,tgateamplitude(2,:),"Marker","o","LineStyle",'-'); grid on; hold on;
xl = xline(20.5,'-.','Paint Edge','DisplayName','Paint Edge'); hold on;
legend('First Reflection','Second Reflection', 'Location','northwest')
title('Max Amplitude of first and second reflection vs Radial Position');
%subplot(2, 1, 2);
%plot(1:43,tgateamplitude(2,:),"Marker","o","LineStyle",'-');
%xl = xline(20.5,'-.','Paint Edge','DisplayName','Paint Edge'); hold on;

%time of flight  \/ \/
tof = maxIndex*dt*1e6;
figure(19);
plot(1:43,tof,'LineStyle','-','Marker','o'); grid on; hold on;
xl = xline(20.5,'-.','Paint Edge','DisplayName','Paint Edge'); hold on;
xlabel('Radial Samples')
ylabel('TOF(mus)')
title('Time of Flight Vs Radial Position')
