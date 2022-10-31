clear all
clc
tic

Limit = [];
CapPerCust = [];
Energy = [];
AptNumber = [];
HHNumber = [];
FuseLimit = [];
Cap = [];

for p=1:100
    [CapPC En HHNum AptNum cap] = Main_model_SWE(p);
    
    %Limit = [Limit;Lim];
    CapPerCust = [CapPerCust CapPC];
    Cap = [Cap cap];
    Energy = [Energy En];
    HHNumber = [HHNumber HHNum];
    AptNumber = [AptNumber AptNum];
    
    %FuseLimit = [FuseLimit;FuseL'];
end

save('Sweden_LimitsPVResults.mat')

toc

