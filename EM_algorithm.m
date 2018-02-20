close all;
clear all;

data=load('data1.txt'); %Loading data
%Initializing distributions
gaussNum=3;
dataSize=length(data);
priors=ones(1,gaussNum)*1/gaussNum;
means=[0.1 0.4 0.8];%0.3 0.8 0.9]*1;%rand(1,3);
stDev=ones(1,gaussNum)*2;
pCX=zeros(gaussNum,dataSize);
pXC=zeros(dataSize,gaussNum);
pX=zeros(dataSize,1);
N=zeros(1,gaussNum);
pCandX=pCX;lLL=-1000;oldlLL=-10000;
cnt=0;
while (lLL(end)-oldlLL)>0.001
        for i=1:gaussNum
            pXC(:,i)=normpdf(data,means(i),stDev(i));
        end
    cnt=cnt+1;
    pCandX=(pXC.*priors)';
    pX=(sum(pCandX,1))';
    oldlLL=lLL(end);
    lLL=[lLL;sum(log10(pX),1)];
    pCX=pCandX./pX';%pX';
    N=(sum(pCX,2))';
    means=(sum(pCX.*(data'),2))'./N;
    stDev=sqrt((sum(pCX.*((data').^2),2))'./N-means.^2);
    priors=N/sum(N);
end
hold on;
plot(2:cnt,lLL(3:end),'DisplayName','2 Gaussians');
xlabel('iritations');
ylabel('log likelihood');
title('log likelihood vs iritations');
lg=legend('show');
lg.FontSize=10;
