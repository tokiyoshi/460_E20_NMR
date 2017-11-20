close all
clear
startNum = 17;
endNum = 19;
C = cbrewer('qual', 'Paired', 6);

for fileNum = startNum:endNum
filePath = char("C:\Users\Benjamin's PC\Dropbox\School Work\PHYS460A\Experiment 20\460_E20_NMR\rawdata\");
% fileNum = 7;
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
winSize = 50;
vFilt = medfilt1(vBG,winSize,'truncate');

%% Find time difference
[pks,loc] = findpeaks(abs(vBG),'SortStr','descend');
numPeaks = 2;
pkLocs = loc(1:numPeaks);
dT(fileNum-startNum+1) = abs(t(pkLocs(2)) - t(pkLocs(1)));

%% Plot final data 
figure(1); hold on;
plot(t,vBG,'Color',C(2*(fileNum-startNum+1),:));
% plot(t,vFit);
scatter(t(pkLocs),vBG(pkLocs),[],C(2*(fileNum-startNum+1)-1,:));
% plot(t,heaviside(t-0.05).*9.*exp(-t./1))
xlabel('Time, ms')
ylabel('Voltage, V')
end

tau = mean(dT);
tauE = max(mean(dT)-dT);

T1 = tau/log(2)
T1E = tauE/log(2)
