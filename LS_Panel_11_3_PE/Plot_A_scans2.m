
% % filepath = ('C:\Users\mcapriotti\Documents\SDSU\Research\Projects\Lightning Strike\LS_Panel_11_3_PE');
% % filename = ([filepath,'\a']);
% % 1=time, 2=Ch1 (trigger), 3=Ch2 (data), 4=Ref (=ref0)
%% Loading raw csv
% Ns = size(ref1,1);
% t = ref1(3:end,1);
% dt = t(10)-t(9); fs = 1/dt;
% time = 0:dt:(Ns-1)*dt;

% trigger = ref1(3:end,2);
% ref0 = ref1(3:end,4);
% ref1 = ref1(3:end,3);
% ref2 = ref2(3:end,3);
% ref3 = ref3(3:end,3);
% ref4 = ref4(3:end,3);
% ref5 = ref5(3:end,3);
% ref6 = ref6(3:end,3);
% ref7 = ref7(3:end,3);
% a = a(3:end,3); a1 = a1(3:end,3);
% b = b(3:end,3); b1 = b1(3:end,3); b2 = b2(3:end,3);
% c = c(3:end,3); c1 = c1(3:end,3);
% d = d(3:end,3); d1 = d1(3:end,3);
% e = e(3:end,3);
% e1 = e1(3:end,3);e2 = e2(3:end,3);e3 = e3(3:end,3);e4 = e4(3:end,3);e5 = e5(3:end,3);
% f = f(3:end,3); g = g(3:end,3);
% h = h(3:end,3); ii = i1(3:end,3);
% j = j(3:end,3); k = k(3:end,3);
% l = l(3:end,3); m = m(3:end,3);
% n = n(3:end,3); o = o(3:end,3);
% p = p(3:end,3); q = q(3:end,3);
% r = r(3:end,3); s = s(3:end,3);
% t = t(3:end,3); u = u(3:end,3);
% v = v(3:end,3); w = w(3:end,3);
% x = x(3:end,3); y = y(3:end,3);
% z = z(3:end,3); 

%% Loading saved .mat file
% References: ref0-7 (ref0=saveds as Ref1 on Oscilloscope), b1 (new!)
figure(1);
plot(time*1e6,trigger);grid on;hold on; 

plot(time*1e6,ref0(:,1),'k');grid on;hold on; 
plot(time*1e6,ref1(:,1));grid on;hold on; 
plot(time*1e6,ref2(:,1));grid on;hold on; 
plot(time*1e6,ref3(:,1));grid on;hold on; 
plot(time*1e6,ref4(:,1));grid on;hold on; 
plot(time*1e6,ref5(:,1));grid on;hold on; 
plot(time*1e6,ref6(:,1));grid on;hold on; 
% plot(time*1e6,ref7(:,1));grid on;hold on; 
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
% plot(time,a1(:,1));grid on;hold on; 
plot(time*1e6,l(:,1));grid on;hold on; 
% plot(time*1e6,m(:,1));grid on;hold on; 
plot(time*1e6,e1(:,1));grid on;hold on; 
plot(time*1e6,e2(:,1));grid on;hold on; 
plot(time*1e6,e3(:,1));grid on;hold on; 
% plot(time*1e6,e4(:,1));grid on;hold on;
plot(time*1e6,e5(:,1));grid on;hold on; 
% plot(time*1e6,v(:,1));grid on;hold on; 
plot(time*1e6,k(:,1));grid on;hold on; 
plot(time*1e6,x(:,1));grid on;hold on; 
% plot(time*1e6,y(:,1));grid on;hold on; 
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
plot(time*1e6,ii(:,1));grid on;hold on;
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
i = 37; % ref0=1, ref1=2, ref7=8; 
% i=9,w; 10-l, 11-e1, 12-e2. 13-e3, 14-e4, 15-e5, 16-k, 17-x
% i=18-a, 19-b, 20-c, 21-d, 22-e
% i=23-f, 24-g, 25-h, 26-ii, 27-j
% i=28-d1, 29-n, 30-o, 31-p, 32-q
% i=33-c1, 34-r, 35-s, 36-t, 37-u
signal = u;
[ampl(i) ind] = max(abs(signal(1201:1701)),[],1); % max bw 2.4mus-3.4mus
tof(i) = time(1201+ind)*1e6; % mus

