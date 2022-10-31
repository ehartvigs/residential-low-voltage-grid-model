clear all
clc
tic

Limit = [];
Likelihood = [];
PopDensity = [];
AptNumber = [];
HHNumber = [];
FuseLimit = [];
xCord = [];
yCord = [];
LP = reshape(1:35,5,7);

% for mn = 1:7
%     parfor p=1:5
%         [FuseL Lim Like Pop HHNum AptNum XCord YCord] = MainModelEV_SWE_2022();
% 
%         Limit = [Limit;Lim];
%         Likelihood = [Likelihood;Like];
%         PopDensity = [PopDensity;Pop];
%         HHNumber = [HHNumber HHNum];
%         AptNumber = [AptNumber AptNum];
%         FuseLimit = [FuseLimit FuseL'];
%         xCord = [xCord; XCord];
%         yCord = [yCord; YCord];
%     end
% end

for d = 1:1
    [tmp] = Main_model_EV();
    Likelihood = [Likelihood tmp];
end
toc

