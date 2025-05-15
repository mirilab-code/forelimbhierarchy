%oneway CCM
function [sugiSource, source_cut] = ccm_qoneway(source, target, tau, E, tp)
%%
%tp is time delay added, which is set to 0 if not specified.
%tp is the extended TE 
%matlab is dumb!    
if tp >0
    error('Bad input')
end
%%
%shift time series based on timelags. 
source=reshape(source,[],1);
target=reshape(target,[],1);

source=source(1-tp:end);
target=target(1:end+tp);

%%
L=length(source); %total length
T=1+(E-1)*tau; %samples to remove/not look at
target_sm=zeros((L-T+1),E);
SugiN=E+1; % number of nearest neighbors to look for
N = L-T+1; %number of timesteps
%% RECONTRUCTIONS OF ORIGINAL SYSTEMS

for t=1:(L-T+1)
    target_sm(t,:)=target((T+t-1):-tau:(T+t-1-(E-1)*tau));
end
%%
sugiSource=zeros(N,1);

%% CCM

%a little confused here, but this weird indexing seems to only look
%for nearest neighbors in HALF the data. (not all of it). 


dat=floor((L-T+1)/2); %i think we need to figure out this later

%i think this just iterates through all timesteps, but for some reason
%very mildly confusing
%pb = progressbar('searching for nearest neighbors n stuff... ');
for ii=(dat+1):(L-T+1)
    %pb.print(ii, L-T+1);
    [n1s,d1s]=knnsearch(target_sm((ii-dat):(ii-1),:),target_sm(ii,:),'k',SugiN);
    %same for y

    u1s=exp(-d1s/d1s(1)); %some weird normalizing
    w1s=u1s/sum(u1s); %also some weird normalizing
    sugiSource(ii)= w1s*source(n1s+T-1+ii-(dat+1)); %lets get the og point activity, weighted
end

source_cut=source(T:end);
source_cut=source_cut((dat+1):(L-T+1));
sugiSource=sugiSource((dat+1):(L-T+1));

%%
%SugiCorr=zeros(2,1);

%SugiCorr1=corrcoef(origY,SugiY,'Rows','complete');
%SugiCorr(2,1)=SugiCorr1(1,2);

%SugiCorr2=corrcoef(origX,SugiX,'Rows','complete');
%SugiCorr(1,1)=SugiCorr2(1,2);

% SugiR(2,1)=sqrt((sum((origY-SugiY).^2)/numel(origY)))/std(origY);
% SugiR(1,1)=sqrt((sum((origX-SugiX).^2)/numel(origX)))/std(origX);

end