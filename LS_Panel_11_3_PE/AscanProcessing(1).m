clc; clear; close all;
filepath = ("C:\Users\cb232\Research\NDE Research\Maio Code\LS_Panel_11_3_PE");
tests = ["ref1", 'ref2', 'ref3', 'ref4', 'ref5', 'ref6', 'ref7', 'a', 'a1', 'b', 'b1', 'b2', 'c', 'c1', 'd', 'd1', 'e', 'e1', 'e2', 'e3', 'e4', 'e5', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'];
% 1=time, 2=Ch1 (trigger), 3=Ch2 (data), 4=Ref (=ref0)
%% Loading raw csv
filenames = strings([1, 43]);
for ii = 1:length(tests)
    filename = filepath + "\" + tests(ii) + ".csv";
    filenames(1, ii) = filename;
end
ref1 = readmatrix(filenames(1));
ref2 = readmatrix(filenames(2));
ref3 = readmatrix(filenames(3));
ref4 = readmatrix(filenames(4));
ref5 = readmatrix(filenames(5));
ref6 = readmatrix(filenames(6));
ref7 = readmatrix(filenames(7));
a = readmatrix(filenames(8));
a1 = readmatrix(filenames(9));
b = readmatrix(filenames(10));
b1 = readmatrix(filenames(11));
b2 = readmatrix(filenames(12));
c = readmatrix(filenames(13));
c1 = readmatrix(filenames(14));
d = readmatrix(filenames(15));
d1 = readmatrix(filenames(16));
e = readmatrix(filenames(17));
e1 = readmatrix(filenames(18));
e2 = readmatrix(filenames(19));
e3 = readmatrix(filenames(20));
e4 = readmatrix(filenames(21));
e5 = readmatrix(filenames(22));
f = readmatrix(filenames(23));
g = readmatrix(filenames(24));
h = readmatrix(filenames(25));
i = readmatrix(filenames(26));
j = readmatrix(filenames(27));
k = readmatrix(filenames(28));
l = readmatrix(filenames(29));
m = readmatrix(filenames(30));
n = readmatrix(filenames(31));
o = readmatrix(filenames(32));
p = readmatrix(filenames(33));
q = readmatrix(filenames(34));
r = readmatrix(filenames(35));
s = readmatrix(filenames(36));
t = readmatrix(filenames(37));
u = readmatrix(filenames(38));
v = readmatrix(filenames(39));
w = readmatrix(filenames(40));
x = readmatrix(filenames(41));
y = readmatrix(filenames(42));
z = readmatrix(filenames(43));

Ns = size(ref1,1)-2;
tvec = ref1(3:end-1,1);
dt = tvec(10)-tvec(9); fs = 1/dt;
time = 0:dt:(Ns-1)*dt;

trigger = ref1(3:end,2);
ref0 = ref1(3:end,4);
ref1 = ref1(3:end,3);
ref2 = ref2(3:end,3);
ref3 = ref3(3:end,3);
ref4 = ref4(3:end,3);
ref5 = ref5(3:end,3);
ref6 = ref6(3:end,3);
ref7 = ref7(3:end,3);
a = a(3:end,3); a1 = a1(3:end,3);
b = b(3:end,3); b1 = b1(3:end,3); b2 = b2(3:end,3);
c = c(3:end,3); c1 = c1(3:end,3);
d = d(3:end,3); d1 = d1(3:end,3);
e = e(3:end,3);
e1 = e1(3:end,3);e2 = e2(3:end,3);e3 = e3(3:end,3);e4 = e4(3:end,3);e5 = e5(3:end,3);
f = f(3:end,3); g = g(3:end,3);
h = h(3:end,3); i = i(3:end,3);
j = j(3:end,3); k = k(3:end,3);
l = l(3:end,3); m = m(3:end,3);
n = n(3:end,3); o = o(3:end,3);
p = p(3:end,3); q = q(3:end,3);
r = r(3:end,3); s = s(3:end,3);
t = t(3:end,3); u = u(3:end,3);
v = v(3:end,3); w = w(3:end,3);
x = x(3:end,3); y = y(3:end,3);
z = z(3:end,3); 
%% Loading saved .mat file
% References: ref0-7 (ref0=saveds as Ref1 on Oscilloscope), b1 (new!)
figure(1);
plot(time*1e6,trigger);grid on;hold on; 

%plot(time*1e6,ref0(:,1),'k');grid on;hold on; 
plot(time*1e6,ref1(:,1));grid on;hold on; 
plot(time*1e6,ref2(:,1));grid on;hold on; 
plot(time*1e6,ref3(:,1));grid on;hold on; 
plot(time*1e6,ref4(:,1));grid on;hold on; 
plot(time*1e6,ref5(:,1));grid on;hold on; 
plot(time*1e6,ref6(:,1));grid on;hold on; 
plot(time*1e6,ref7(:,1));grid on;hold on; 
ylim([-0.5 0.5])
title('Reference');xlabel('time [\mus]')

% % References: b2 (new!, old b1) on no paint edge
% figure(7);
% plot(time*1e6,trigger);grid on;hold on; 
% plot(time*1e6,b2(:,1));grid on;hold on; 
% ylim([-0.5 0.5])
% title('no paint edge');xlabel('time [\mus]')

% Pristine close to paint chip: w, a1, l, e1-5,v,k,x,y,m
figure(2);
plot(time*1e6,trigger);grid on;hold on; 
plot(time*1e6,w(:,1));grid on;hold on;
plot(time,a1(:,1));grid on;hold on; 
plot(time*1e6,l(:,1));grid on;hold on; 
plot(time*1e6,m(:,1));grid on;hold on; 
plot(time*1e6,e1(:,1));grid on;hold on; 
plot(time*1e6,e2(:,1));grid on;hold on; 
plot(time*1e6,e3(:,1));grid on;hold on; 
plot(time*1e6,e4(:,1));grid on;hold on;
plot(time*1e6,e5(:,1));grid on;hold on; 
plot(time*1e6,v(:,1));grid on;hold on; 
plot(time*1e6,k(:,1));grid on;hold on; 
plot(time*1e6,x(:,1));grid on;hold on; 
plot(time*1e6,y(:,1));grid on;hold on; 
ylim([-0.5 0.5])
title('Pristine (close to edge)');xlabel('time [\mus]')

% No paint close to edge: a, b, c, d, e 
figure(3);
plot(time*1e6,trigger);grid on;hold on; 
plot(time*1e6,ref0(:,1),'k');grid on;hold on; 

plot(time*1e6,a(:,1));grid on;hold on;
plot(time*1e6,b(:,1));grid on;hold on;
plot(time*1e6,c(:,1));grid on;hold on;
plot(time*1e6,d(:,1));grid on;hold on;
plot(time*1e6,e(:,1));grid on;hold on;
ylim([-0.5 0.5])
title('No paint (further from LS center)');xlabel('time [\mus]')

% No paint further from edge: f, g, h, i, j
figure(4);
plot(time*1e6,trigger);grid on;hold on; 
plot(time*1e6,ref0(:,1),'k');grid on;hold on; 

plot(time*1e6,f(:,1));grid on;hold on;
plot(time*1e6,g(:,1));grid on;hold on;
plot(time*1e6,h(:,1));grid on;hold on;
plot(time*1e6,i(:,1));grid on;hold on;
plot(time*1e6,j(:,1));grid on;hold on;
ylim([-0.5 0.5])
title('No paint (closer to LS center)');xlabel('time [\mus]')

% Damage (black area): n, d1,o,p,q
figure(5)
plot(time*1e6,trigger);grid on;hold on; 
plot(time*1e6,ref0(:,1),'k');grid on;hold on; 

plot(time*1e6,d1(:,1));grid on;hold on; 
plot(time*1e6,n(:,1));grid on;hold on; 
plot(time*1e6,o(:,1));grid on;hold on; 
plot(time*1e6,p(:,1));grid on;hold on; 
plot(time*1e6,q(:,1));grid on;hold on; 
ylim([-0.5 0.5])
title('Severe Damage');xlabel('time [\mus]')

% Damage (delam and severe): r,s,t,u,c1
figure(6)
plot(time*1e6,trigger);grid on;hold on; 
plot(time*1e6,ref0(:,1),'k');grid on;hold on; 

plot(time*1e6,c1(:,1));grid on;hold on; 
plot(time*1e6,r(:,1));grid on;hold on; 
plot(time*1e6,s(:,1));grid on;hold on; 
plot(time*1e6,t(:,1));grid on;hold on; 
plot(time*1e6,u(:,1));grid on;hold on; 
ylim([-0.5 0.5])
title('Severe Damage/Delam');xlabel('time [\mus]')

%% Calculate TOF and AMP of max
jj = 37; % ref0=1, ref1=2, ref7=8; 
% i=9,w; 10-l, 11-e1, 12-e2. 13-e3, 14-e4, 15-e5, 16-k, 17-x
% i=18-a, 19-b, 20-c, 21-d, 22-e
% i=23-f, 24-g, 25-h, 26-ii, 27-j
% i=28-d1, 29-n, 30-o, 31-p, 32-q
% i=33-c1, 34-r, 35-s, 36-t, 37-u
signal = u;
[ampl(jj) ind] = max(abs(signal(1201:1701)),[],1); % max bw 2.4mus-3.4mus
tof(jj) = time(1201+ind)*1e6; % mus

ref0(1:1201) = 0; ref0(1701:end) = 0;
ref1(1:1201) = 0; ref1(1701:end) = 0;
ref2(1:1201) = 0; ref2(1701:end) = 0;
ref3(1:1201) = 0; ref3(1701:end) = 0;
ref4(1:1201) = 0; ref4(1701:end) = 0;
ref5(1:1201) = 0; ref5(1701:end) = 0;
ref6(1:1201) = 0; ref6(1701:end) = 0;
ref7(1:1201) = 0; ref7(1701:end) = 0;
a(1:1201) = 0; a(1701:end) = 0;
a1(1:1201) = 0; a1(1701:end) = 0;
b(1:1201) = 0; b(1701:end) = 0;
b1(1:1201) = 0; b1(1701:end) = 0;
b2(1:1201) = 0; b2(1701:end) = 0;
c(1:1201) = 0; c(1701:end) = 0;
c1(1:1201) = 0; c1(1701:end) = 0;
d(1:1201) = 0; d(1701:end) = 0;
d1(1:1201) = 0; d1(1701:end) = 0;
e(1:1201) = 0; e(1701:end) = 0;
e1(1:1201) = 0; e1(1701:end) = 0;
e2(1:1201) = 0; e2(1701:end) = 0;
e3(1:1201) = 0; e3(1701:end) = 0;
e4(1:1201) = 0; e4(1701:end) = 0;
e5(1:1201) = 0; e5(1701:end) = 0;
f(1:1201) = 0; f(1701:end) = 0;
g(1:1201) = 0; g(1701:end) = 0;
h(1:1201) = 0; h(1701:end) = 0;
i(1:1201) = 0; i(1701:end) = 0;
j(1:1201) = 0; j(1701:end) = 0;
k(1:1201) = 0; k(1701:end) = 0;
l(1:1201) = 0; l(1701:end) = 0;
m(1:1201) = 0; m(1701:end) = 0;
n(1:1201) = 0; n(1701:end) = 0;
o(1:1201) = 0; o(1701:end) = 0;
p(1:1201) = 0; p(1701:end) = 0;
q(1:1201) = 0; q(1701:end) = 0;
r(1:1201) = 0; r(1701:end) = 0;
s(1:1201) = 0; s(1701:end) = 0;
t(1:1201) = 0; t(1701:end) = 0;
u(1:1201) = 0; u(1701:end) = 0;
v(1:1201) = 0; v(1701:end) = 0;
w(1:1201) = 0; w(1701:end) = 0;
x(1:1201) = 0; x(1701:end) = 0;
y(1:1201) = 0; y(1701:end) = 0;
z(1:1201) = 0; z(1701:end) = 0;

% Pristine close to paint chip: w, a1, l, e1-5,v,k,x,y,m
RMSvec = [rms(t), rms(r),rms(u),rms(s),rms(c1),rms(n),rms(d1),rms(o),rms(p),rms(q),rms(f),rms(g),rms(h),rms(i),rms(j),rms(a),rms(b),rms(c),rms(d),rms(e),rms(w),rms(a1),rms(e1),rms(e2),rms(e3),rms(e4),rms(e5),rms(l),rms(v),rms(k),rms(x),rms(y),rms(m),rms(z),rms(b1),rms(ref0),rms(ref1),rms(ref2),rms(ref3), rms(ref4), rms(ref5), rms(ref6), rms(ref7)];
    
figure(10)
plot(1:43, RMSvec, "LineWidth", 2); hold on; grid on;
% RMS POSITIONS PROPAGATE RADIALLY: Left side = plate center, right side = far from center
% Pos 1-5 = severe / delam, 6-10 = black area, 11-15 = no paint further
% from edge, 16-20 = no paint close to edge, 21-34 = pristine close to paint chip, 35-42 = reference locations 
title('RMS (All positions');xlabel('Position from center'); ylabel('RMS');

% Frequency Domain Analysis
%freq = (0:Ns-1)*(fs/Ns);
freqB=(0:fs/Ns:fs - fs/Ns)';
freq = freqB(1:Ns/2)/1000;  
zzz = fft(ref1);
figure(11)
plot(freq(1:Ns/2), abs(zzz(1:Ns/2))); hold on; grid on;












