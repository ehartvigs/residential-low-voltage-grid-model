clear all
clc
close all

steplength =1000;



%% First part is to repeat the driving patterns to get a full year of data.
%% The original data from Stens project was measured for 30-60 days and thus
%% does not cover a whole year.


% load BRDChargingPowerAndSoCsBatterySize100.mat
% 
% 
% %load BRDChargingData.mat
% 
% 
% 
% BRDChargingPower69_Home_ED17 = chargingPowerBRD.chargingPowerHome_ED017_CR69;
% NewBRDCharging = zeros(52704,429);
% 
% % First make it into full year data
% DayIndices = 0:6*24:52704;
% 
% for k = 1:429
%     start = find(BRDChargingPower69_Home_ED17(:,k)==1,1);
%     finish = find(BRDChargingPower69_Home_ED17(:,k)==1,1,'last');
%     
%     if finish-start > 18000
%         start = find(BRDChargingPower69_Home_ED17(1:26000,k)==1,1,'last');
%         finish = find(BRDChargingPower69_Home_ED17(26000:end,k)==1,1,'first')+26000;
%         [tmp startpos] = min(abs(DayIndices-start));
%         [tmp endpos] = min(abs(DayIndices-finish));
%         
%         FirstHalf = BRDChargingPower69_Home_ED17(1:startpos*6*24,k);
%         ThirdHalf = BRDChargingPower69_Home_ED17(endpos*6*24:end,k);
%         
%         NoCopies = floor((endpos-startpos)/((366-endpos)+startpos));
%         SecHalf = repmat([BRDChargingPower69_Home_ED17((endpos)*6*24:end,k) ;BRDChargingPower69_Home_ED17(1:startpos*24*6,k)],NoCopies,1);
%         DaysToTakeFromData = int16((((endpos-startpos)/((366-endpos)+startpos))-floor((endpos-startpos)/((366-endpos)+startpos)))*((366-endpos)+startpos));
%         Set = [BRDChargingPower69_Home_ED17((endpos)*6*24:end,k) ;BRDChargingPower69_Home_ED17(1:startpos*24*6,k)];
%         SecHalf = [SecHalf;Set(1:DaysToTakeFromData*6*24)];
%         
%         tmp = [FirstHalf;SecHalf;ThirdHalf];
%         NewBRDCharging(:,k) = tmp(1:52704);
%     else
%         [tmp startpos] = min(abs(DayIndices-start));
%         [tmp endpos] = min(abs(DayIndices-finish));
%         
%         NoCopies = floor((366-endpos)/(endpos-startpos))+1;
%         SecHalf = repmat(BRDChargingPower69_Home_ED17((startpos)*6*24:(endpos)*6*24,k),NoCopies,1);
%         DaysToTakeFromData = (366-endpos)-floor((366-endpos)/(endpos-startpos))*(endpos-startpos);
%         SecHalf = [SecHalf;BRDChargingPower69_Home_ED17((startpos)*6*24:(startpos+DaysToTakeFromData)*6*24,k)];
%         
%         NoCopies = floor((startpos)/(endpos-startpos));
%         FirstHalf = repmat(BRDChargingPower69_Home_ED17((startpos)*6*24:(endpos)*6*24,k),NoCopies,1);
%         DaysToTakeFromData = startpos-floor((startpos)/(endpos-startpos))*(endpos-startpos);
%         FirstHalf = [FirstHalf;BRDChargingPower69_Home_ED17((startpos)*6*24:(startpos+DaysToTakeFromData)*6*24,k)];
%         
%         tmp = [FirstHalf;SecHalf];
%         NewBRDCharging(:,k) = tmp(1:52704);
%         
%     end
%     
% 
% end
% 
% BRDChargingAnnualHomeED17CR69 = NewBRDCharging(1:52704,:);
% 
% save('BRDChargingAnnualHomeED17CR69100kWh.mat','BRDChargingAnnualHomeED17CR69')



%% This part creates the coincidence file using the annual profile of charging
% demand.


load BRDChargingAnnualHomeED17CR69100kWh.mat

BRDChargingPower69_Home_ED17 = BRDChargingAnnualHomeED17CR69;

coincidence_ED017_CR69 = zeros(31,100,52704);


% I need to change something to get the coincidence (the way it
% is now is because data has not been extrapolated to cover all days


tic
for k = 9:31
   clear tmpSum
    parfor j = 1:1000
        ix = randperm(429,k);
        ix = ix';

        tmpSum(:,j) = sum(BRDChargingPower69_Home_ED17(1:52704,ix)')./(k);
        
    end
    [ii,kk] = sort(max(tmpSum));
    bb = tmpSum(:,kk);
    bb = bb(:,1:10:1000);
    coincidence_ED017_CR69(k,:,:) = reshape(bb',[1 100 52704]);
end
toc



%save('CoincidenceLarge100kWh.mat','coincidence_ED017_CR69','-v7.3')









