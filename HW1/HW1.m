clear all; 
close all; clc;

% initialization
dF=0.001;   % Frequency resolution
F=0:dF:0.5; % Frequency axis from 0 to 0.5
Hd=zeros(1,length(F)); 
fs=6000;
ftrans=1350;
trans_band_low=1200;
trans_band_high=1500;
Hd(1:(ftrans/fs)/dF+1)=1;                              % Ideal desired filter Hd(F)
trans_band=[(trans_band_low/fs) (trans_band_high/fs)]; % Transition band boundary

% Weighting function 
W=zeros(1,length(F)); 
W(1:trans_band(1)/dF+1)=1;    % W(F)=1 for 0<F<trans_band_lFt
W(trans_band(2)/dF+1:end)=0.6;    % W(F)=0.6 for trans_band_right0<F<0.5
N=17;       % Filter length N
k=(N-1)/2;
Fi=[0 0.04 0.1 0.12 0.15 0.19 0.3 0.38 0.44 0.5]; %initial guess
delta=0.0001;                                  
%%%% MiniMaxFIR %%%%
dF=F(2)-F(1); % Frequency axis
k=(N-1)/2; 
E=[]; 
for i=1:17
    %%%%%%   Find Filter Parameter s[n]  %%%%%%
    n=0:k;  
    m=0:k+1;
    A=[cos(2*pi*Fi.'*n) ((-1).^m./W(round(Fi/dF+1))).'];
    H=Hd(round(Fi/dF+1)).';
    Sn=inv(A)*H;       % inverse matrix
    s=Sn(1:k+1).';
    R = s * cos(2*pi*n.'*F);
    err = (R - Hd) .* W;    % compute err(F)
    % new local maximal & check the number of the extreme points
    prev=err(1:end-2); 
    curr=err(2:end-1); 
    next=err(3:end);
    Fi=F(find(prev<curr & next<curr | prev>curr & next>curr)+1);
    if sign(err(1))~=sign(err(2)-err(1)) 
        Fi=[F(1) Fi];
    end
    if sign(err(end))~=sign(err(end-1)-err(end))
        Fi=[Fi F(end)];
    end      
    if length(Fi) == k+3
        BoundinFi=find(Fi==F(1) | Fi==trans_band(1) | Fi==trans_band(2) | Fi==F(end));  
        BoundErr=err(Fi(BoundinFi)/dF+1);                                               
        [minerr,index]=min(abs(BoundErr));                                              
        Fi(BoundinFi(index))=[];                                                        
    elseif length(Fi) > k+3
        disp('Number of extreme frequency exceed k+3');
    end

    maxerr=max(abs(err(round(Fi/dF+1))));
    E=[E maxerr];
    if length(E) >= 2
        if abs(E(end-1)-E(end)) < delta   % |E1-E0| < delta
            break;
        end
    end
end

figure; 
plot(F,Hd,F,R); 
hold on;
for m=1:k+2
    plot([Fi(m) Fi(m)],[-0.5 1.5],'m:');
end
legend('Hd(F)','H(F)','Fm');
title('Frequency Response'); 
xlabel('Normalized Frequency');  

%Maximal error for each iteration
disp(sprintf('Maximal error :'));disp(E);
% impulse response h
h=zeros(1,N);
h(k+1)=s(1);
h(1:k)=s(k+1:-1:2)/2;
h(k+2:N)=s(2:k+1)/2;
n=0:N-1;
figure; 
stem(n,h);
title('Impulse Response'); 
xlabel('n');
disp(sprintf('Impulse Response :'));disp(h);