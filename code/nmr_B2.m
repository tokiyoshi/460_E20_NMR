close all
clear
sampleNames = {'Doped Water','Ethanol','Rubber'};
counter = 1;
for fileNum = 3:5
filePath = char("C:\Users\Benjamin's PC\Dropbox\School Work\PHYS460A\Experiment 20\460_E20_NMR\rawdata\");
% fileNum = 3;
fileName = ['T',sprintf('%04.0f',fileNum),'CH1.CSV'];
data = csvread([filePath,fileName],16,0);

t = data(:,1)*1e3;
vRaw = data(:,2);

%% Detect background and subtract
bgStart = 1;
bgEnd = 300;
bg = mean(vRaw(bgStart:bgEnd));
vBG = vRaw - bg;

%% Filter the voltage signal
% winSize = 5;
% vFilt = medfilt1(vBG,winSize,'truncate');
vFilt = vBG;

%% Exponential fit
% x = t;
% y = vBG;

% figure(1); hold on;
% plot(x,y)
% [peakX,peakY] = ginput(1);
% close

% fitType = fittype('heaviside(x-x0).*M0.*exp(-(x-x0)/T2)',...
%     'dependent',{'y'},'independent',{'x'},...
%     'coefficients',{'M0','T2','x0'});
% vFitObj = fit(  x,y,fitType,...
%                 'Exclude',abs(y)>13,...
%                 'Lower',[-2*peakY,0.5,peakX-0.01],...
%                 'Upper',[ 2*peakY+1,2,peakX+0.01])
%             
% vFit = feval(vFitObj,t);

%% Semi-log fit
figure; hold on;
plot(t,log(vFilt))
[peakX,peakY] = ginput(2); hold off;
close

[~,xLeft] = min(abs(t-peakX(1))); [~,xRight] = min(abs(t-peakX(2)));
x = t(xLeft:xRight);
y = log(vBG(xLeft:xRight));

fitType = fittype('a1 - x/T2s',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a1','T2s'});
[vFitLogObj,gof] = fit(x,y,fitType)
          
vLog = feval(vFitLogObj,t);

figure(30); hold on;
subplot(1,3,counter); hold on;
plot(t,log(vFilt))
plot(t,vLog)
xlabel('Time (ms)')
ylabel('Log(M)')
title(char(sampleNames(counter)));
axis([0,10,-4,5])
hold off;


vFitLog = exp(vFitLogObj.a1).*exp(-t./vFitLogObj.T2s);

vFitLog(t<0)=NaN;

vFit(:,fileNum-2) = vFitLog;

v(:,fileNum-2) = vBG;
counter = counter+1;
end


figure; hold on;
for ii = 1:3
%% Plot final data 
subplot(3,1,ii); hold on;
plot(t,v(:,ii));
% plot(t,vFit);
plot(t,vFit(:,ii));
xlabel('Time, ms')
ylabel('Voltage, V')
legend('Experimental Curve','Fit')
end