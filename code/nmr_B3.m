close all
clear

filePath = char("C:\Users\Benjamin's PC\Dropbox\School Work\PHYS460A\Experiment 20\460_E20_NMR\rawdata\");
fileNum = 8;
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
winSize = 5;
vFilt = medfilt1(vBG,winSize,'truncate');

%% Fit to function
x = t;
y = vBG;

figure(1); hold on;
plot(x,y)
[peakX,peakY] = ginput(2);
close

fitType = fittype('heaviside(x-x0_1).*M0_1.*exp(-(x-x0_1)/T1_1) + heaviside(x-x0_2).*M0_2.*exp(-(x-x0_2)/T1_2)',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'M0_1','T1_1','x0_1','M0_2','T1_2','x0_2'});
vFitObj = fit(  x,y,fitType,...
                'Exclude',abs(y)>13,...
                'Lower',[-2*peakY(1),0.5,peakX(1)-0.05,-2*peakY(2),0.5,peakX(2)-0.05],...
                'Upper',[ 2*peakY(1)+1,2,peakX(1)+0.05, 2*peakY(2)+1,2,peakX(2)+0.05])
            
vFit = feval(vFitObj,t);

%% Find time difference
% [pks,loc] = findpeaks(abs(vBG),'SortStr','descend');
% numPeaks = 2;
% pkLocs = loc(1:numPeaks);
% dT = abs(t(pkLocs(2)) - t(pkLocs(1)))

%% Plot final data 
figure(2); hold on;
plot(t,vBG);
plot(t,vFit);
% plot(t,heaviside(t-0.05).*9.*exp(-t./1))
xlabel('Time, ms')
ylabel('Voltage, V')
