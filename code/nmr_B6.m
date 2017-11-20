close all
clear

waterPath = char("C:\Users\Benjamin's PC\Dropbox\School Work\PHYS460A\Experiment 20\460_E20_NMR\rawdata\");
fileNum = 3;
fileName = ['T',sprintf('%04.0f',fileNum),'CH1.CSV'];
waterData = csvread([waterPath,fileName],16,0);
waterT = waterData(:,1)*1e3;
waterV = waterData(:,2);
bgStart = 1;
bgEnd = 300;
bg = mean(waterV(bgStart:bgEnd));
waterV = waterV - bg;
waterV = medfilt1(waterV,10);


icePath = char("C:\Users\Benjamin's PC\Dropbox\School Work\PHYS460A\Experiment 20\460_E20_NMR\rawdata\");
fileNum = 89;
fileName = ['T',sprintf('%04.0f',fileNum),'CH1.CSV'];
iceData = csvread([icePath,fileName],16,0);
iceT = iceData(:,1)*1e3;
iceV = iceData(:,2);
bgStart = 1;
bgEnd = 300;
bg = mean(iceV(bgStart:bgEnd));
iceV = iceV - bg;
iceV = medfilt1(iceV,10);

C = cbrewer('qual', 'Set1', 3);

figure(1); hold on;
plot(waterT,waterV,'Color',C(:,1));
plot(iceT,iceV,'Color',C(:,2));
legend('Doped Water','Frozen Doped Water')
xlabel('Time, ms')
ylabel('Magnetization, V')
axis([-1,10,-1,10])

