function []=TimeResolvedNets(TR,K,thr,SNC_size,band,varargin)

% Marika Strindberg, Karolinska Institutet 2018-2020. 
% marika.strindberg@ki.se, marikastrindberg@gmail.com

% The ouput of this code is described in "Spatio-temporally flexible subnetworks reveal the quasi-cyclic nature of
% integration and segregation in the human brain", Strindberg et al, 2021

% (Row 87-113 and 282-298 are adapted from Cabral et al 2017) 

% TR = Timeresolution of aquisition (repetition time)
% K = Expansion in each iteration to find SNCs, default K=10
% thr = is by default 0, if set higher it only searches trough SNCs with qint
% above this thr
% SNC_size = size of final SNC i.e. number of areas
% band = if 1 then bandpass filter will be used, if 0 IMFs will be calculated
% varargin = If predefined seed areas are prefered, add areas as a vector 

% Subrutines:
    % SNCs.m  calculates the subnetworks
    % PhaseIntegrationSNC.m calculates the phase coherence within a subset
    % of empirically derived SNCs and randomly composed SNCs
     
% Example use: TimeResolvedNets(0.72,10,0,8,0) 

% Not that the directories on lines
% 61,67,1331,1341,1590,1599,1644,1653,2538,2547,2570,2579 needs to be utdated to fit the
% local environement.

% BrainNet viewer needs to be installed for the visualizations of networks,
% see ReadMe for further details.
    
% Main calculations: 

% 1a. In the case of IMFs, Decompose signal into IMFs and decide which IMF to use based on powerspectra 
% 1b. In the case of band-pass filter, the user will be prompted to set the
% upper and lower filter thresholds
% 2. Louvain communities are calculuated for all areas for all timepoints and subjects.
% 3. Louvain run is chosen
% 4. Calculation of SNCs
% 5. Grouping of SNCs into SNs
% 6. Calculation of SNC timeseries and SNC phase 
% 7. Calculation of SN timeseries and SN phase
% 8. Calculate integration and segregation of SNs relative to a) time of overlap b) all
% time as well as relative activation/deactivation
% 9. Selection of hierachical cluster solution that maximizes integration within MNs
% 10. Meta-network timeseries and maps
% 11. Naming of MNs
% 12. SN maps
% 13. Video of MN and SN maps
% 17. State vector recurrence plots
% 14. Phase coherence within and between SNs in the same community and
% between SNs in different commmunities
% 15. Phase coherence for a sample of SNC an d random components prior to and during  integration
% 16. Duration of SNs
% 17. Calculation of number of SNs per timepoint,average number of communities SNs are integrated each timepoint, number of areas per SN
% 18. Calculation relative flexibility/modularity of areas

% % 
% % Set the directory where the BOLD timeseries of all subjects are located
%  BOLD=dir('/Users/Documents/MATLAB/BOLDtimeseries/*.txt');
% 
%  S=length(BOLD); 
%  
% %Load BOLD signals into matrix BOLDdata
% for s=1:S
%     BOLDdata(s,:,:)=load(['/Users/Documents/MATLAB/BOLDtimeseries/' BOLD(s).name]);
% end
% 
% save BOLD BOLDdata

load BOLD

[S N T]= size(BOLDdata)


if band == 1
 %%%% 1.Bandpass instead of EMD%%   
    
j=1;
j_common=1;
b1=input("Choose lower threshold for bandpass filter ")
b2=input("Choose upper threshold for bandpass filter ")

[S, N, T]=size(BOLDdata);

% Bandpass filter settings
    fnq=1/(2*TR);                 % Nyquist frequency
    flp = b1;                    % lowpass frequency of filter (Hz)
    fhi = b2;                    % highpass
    Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
    k=2;                          % 2nd order butterworth filter
    [bfilt,afilt]=butter(k,Wn);   % construct the filter
    clear fnq flp fhi Wn k
       
BOLD_band=zeros(S,N,T);
    
for s=1:S  
    BOLD = squeeze(BOLDdata(s,:,:)); 
    % Filter signal and then get phase using the Hilbert transform
    for n=1:N
        BOLD(n,:)=BOLD(n,:)-mean(BOLD(n,:));
        signal_filt =filtfilt(bfilt,afilt,BOLD(n,:));
        BOLD_band(s,n,:) = angle(hilbert(signal_filt));
    end
end


for s=1:S
    for p=1:N
        Phase_imf(s,p,1,:)=angle(hilbert(BOLD_band(s,p,:)));
    end
end

save PhaseAndImf  Phase_imf BOLD_band

    
else 
 %%%% 1.Decompose all timeseries into imperical mode functions%%.    
    
fs=1/TR;
imf_freq=zeros(S,N,T,10);
BOLD_imf=zeros(S,N,T,10);
for s=1:S
    %Decompose signals fo all parcels into imfs
    
    x=squeeze(BOLDdata(s,:,:));
    
    for p=1:N
        
        y=squeeze(x(p,:));
        y=y(:)-mean(y(:));
        [imf,~]=emd(y);
        [hs,~,~,imfinsf,~]=hht(imf,fs); %to get the freq content for power spectrum
        [a b]=size(imf);
        IMFnr(s,p)=b;
        
        clear hs
        
        
        imf_freq(s,p,:,1:b)=imfinsf;
        BOLD_imf(s,p,:,1:b)=imf;
        
        imf=[];
        imfinsf=[];
    end
end
%
j_common=min(IMFnr(:))

disp(sprintf('Minimum number of IMFs common to all areas and all participants: %d',j_common ))

for s=1:S
    for p=1:N
        for i=1:j_common
            Phase_imf(s,p,i,:)=angle(hilbert(BOLD_imf(s,p,:,i)));
        end
    end
end
save PhaseAndImf Phase_imf BOLD_imf imf_freq IMFnr



%Powerspectra for shared IMFs

% Plot Power Spectrum for all nodes&subjects in each IMF

w = T;
Ts = w*TR; % define total duration of sample in seconds
freq = (0:w/2-1)/Ts; % frequency range in Hz

PS=zeros(j_common,numel(freq));
PW=zeros(j_common,numel(freq));
t=0
for imf=1:j_common
    for s=1:S
        for n=1:N
            X=squeeze(BOLD_imf(s,n,:,imf));
            pw = abs(fft(X-mean(X)));
            ps_X=pw(1:floor(w/2)).^2/(w/2);
            PS(imf,:)=PS(imf,:)+ps_X';
            PW=PW(imf,:)+pw;
            t=t+1;
        end
    end
end
  
figure
a=[1:1:j_common]
for i=1:j_common
    hold on
    plot(freq,PS(i,:),'LineWidth',2)
    Psum=sum(PS(i,:));
    [r c]=find(PS(i,:)==max(PS(i,:)));
    Peak(i)=freq(c)
    
end
xlabel('Frequency (Hz)')
ylabel('Power')
%%legend('IMF 1','IMF 2','IMF 3','IMF 4','IMF 5','IMF 6')


%Power spectrum for BOLD signal
PSBold=zeros(1,numel(freq));
for s=1:S
        for n=1:N
            X=squeeze(BOLDdata(s,n,:));
            pw = abs(fft(X-mean(X)));
            ps_X=pw(1:floor(w/2)).^2/(w/2);
            PSBold=PSBold+ps_X';  
        end     
end
hold on
plot(freq,PSBold,'LineWidth',2)
legend


%Range of power spectra for the IMFs: column1=peak, column 2=99% of power
%below this frequency, column 3=99% of power above this frequency
for i=1:j_common
    [r c]=find(PS(i,:)==max(PS(i,:)))
    PsN(i,1)=freq(c)
    a=sum(PS(i,:))
    for j=1:T/2
        a2=sum(PS(i,1:j));
        if a2>=a*0.99 %approx max freq
            PsN(i,2)=freq(j)
            break
        end
    end
end

for i=1:j_common
    a=sum(PS(i,:));
    for j=1:Ts
        a2=sum(PS(i,1:j));
        if a2>=a*0.01 %approx min freq
            PsN(i,3)=freq(j)
            break
        end
    end
end


%How much of the power <0.1 Hz is contained in each IMF
% Find threshold 0.1Hz
F=freq
f=F(freq<=0.1);
for i=1:j_common
    powerBOLDclassic(i)=sum(sum(PS(i,1:numel(f))))/sum(sum(PS(:,1:numel(f))))
end


format short g
fileID = fopen('IMFsummary.txt','w')
fprintf(fileID,'Min number of common IMFs %d \n',j_common)
for i=1:j_common
    a1=powerBOLDclassic(i)
    a2=PsN(i,1)
    a3=PsN(i,3)
    a4=PsN(i,2)
    fprintf(fileID,'IMF %d  \n', i)
    fprintf(fileID,'Power <=0.1Hz %d \n', a1)
    fprintf(fileID,'Peak freq (hz) %d \n', a2 )
    fprintf(fileID,'Range (99 percentil of distribution) %d - %d \n',a3, a4)
   
end
fclose(fileID)

j=input("Choose IMF ")

end

%%%

%for each timepoint compute dPC for the band pass filtered signal or the IMFs(j_common), Louvain communites
%division as well as Leida weights
    L=100 %Number of iterations of the Louvain alogrithm
    Louvain_s=zeros(L,S,T-10,N);
   
    %Calculate the dPCs
    for s=1:S
        
        iFCalls=zeros(T-10,N,N);
     
        parfor t=1:T-10 % Five first and last timepoints are not used due to artefacts resulting from the hilbert transformation
            % Calculate phase difference matrixes dFC
            dPC=zeros(N,N);
            for m=1:N
                for n=1:N
                    if n>m
                        dPC(m,n)=cos(Phase_imf(s,m,j,t+5)-Phase_imf(s,n,j,t+5));
                    else
                    end
                end
            end
            dPC=dPC+tril(dPC',1);
            iFCalls(t,:,:)=dPC;
           
           
            % Calculate the Louvain communities
           
            for k=1:L
                CLPs=community_louvain(dPC,1,[],'negative_sym');
                Louvain_s(k,s,t,:)=CLPs;
            end
            
        end
        
        save(sprintf('iFCall_s%d',s), 'iFCalls', '-v7.3')
    end
    save(sprintf('Louvain%d_%d',j,K), 'Louvain_s',  '-v7.3')
    
  
    L=100
    Tmax=T-10;
    PairW=zeros(L,S,N,N);
    
    for k=1:L
        
        for s=1:S
            
            A=squeeze(Louvain_s(k,s,:,:));
            
            for n=1:N
                for m=1:N
                    if m>=n
                        PairW(k,s,n,m)=numel(nonzeros(A(:,n)==A(:,m)))/Tmax;
                        PairW(k,s,m,n)=numel(nonzeros(A(:,n)==A(:,m)))/Tmax;
                        
                    else
                    end
                end
            end
        end
    end
    
    for k=1:L
        for n=1:N
            for m=1:N
                Mean_matrix(k,n,m)=mean(PairW(k,:,n,m));
                Min_matrix(k,n,m)=min(PairW(k,:,n,m));
                Max_matrix(k,n,m)=max(PairW(k,:,n,m));
                Var_matrix(k,n,m)=std(PairW(k,:,n,m));
            end
        end
    end
    

    save(sprintf('PairWise_K%d_%d',K,j), 'Mean_matrix', 'Min_matrix', 'Max_matrix' ,'Var_matrix', '-v7.3')


    
% Find the most representative Louvain solution ie the one that contains most pairs with highest frequency accross the iterations    
    
QM=zeros(N,N);
OrdAll=[];
for k=1:L
    
    V=squeeze(Mean_matrix(k,:,:));
    V=V-diag(diag(V));
    V=triu(V);
    [row col]=find(V>0);
    R=[row col];
    
    R=sort(R,2);
    R=unique(R,'rows');
    
    for i=1:length(R)
        Qr(i)=V(R(i,1),R(i,2));
    end
    
    QR=[R Qr'];
    QRs=sortrows(QR,3,'descend');
    
    Ord=[];
    for i=1:N
        temp=[];r=[];tempmax=[];
        [r c]=find(QRs==i); 
        temp=QRs(r,:);
        tempmax=max((QRs(r,3))); % BUGfix added 6 july
        [rx ~]=find(temp==tempmax) %
        Ord(i)=r(rx(1));   %
    end
    Ord=unique(Ord); % Remove dublicates
    rr=find(Ord>length(R));
    
    if rr>0
        for i=1:numel(rr)
            r2(i)=Ord(rr(i))-length(R);
        end
        Ord(rr)=r2;
    end
    QRtop=QRs(Ord,:);
    
    %test that all areas are included
   
    QRtopS(k,1:length(Ord),:)=sortrows(QRtop,3,'descend');
    Q=squeeze(QRtop(:,1:2));
    
    %Count number of top pair frequency
    
    for m=1:length(Ord)
        
        QM(Q(m,1),Q(m,2))=QM(Q(m,1),Q(m,2))+1;
    end
    OrdAll(k)=numel(Ord); %Bugfix 18 July
end
save QRsTOP_IMFiter QRtopS QM OrdAll
figure; imagesc(QM)

QT=QRtopS(:,:,1:2);

MAX=max(QM(:));
[row col]=find(QM==MAX)
f1=[row col];

Hh=unique(nonzeros(QM)); %Different frequencies
HH=sort(Hh,'descend');

CG=zeros(L,N);

if numel(HH)>1 % Added to accomodate the possibility that all Louvain solutions create exactly the same TOP pairs

kk=[];
TOP=[];
TOPs=[];
for n=1:N
    
    [row col]=find(QM==HH(n));
    f=[row col];
    kk=[kk;f];
    if (length(unique(kk,'rows')))<N
        
        for i=1:L
            d=squeeze(QT(i,:,:));
            CG(i,n)= length(unique([kk;d],'rows'))-N;
        end
    else
        
        final=n
        ll=length(kk);
        break
    end
end
SumCG=max(CG');
[r d]=find(SumCG'==min(SumCG));
rof=r(1);% Bugfix 18 July

else  % if all Louvain solutions give exactly the same TOP pairs
   rof=1;% Bugfix 18 July
end

TOP=squeeze(QT(rof,1:OrdAll(rof),:)); % Bugfix 18 July
TOPS=squeeze(QRtopS(rof,1:OrdAll(rof),:)); %Bugfix 18 July
if S>1
    Louv=squeeze(Louvain_s(rof,:,:,:)) % Bugfix 18 July
else
    Louv(1,:,:)=squeeze(Louvain_s(rof,:,:,:)); % Bugfix 18 July
end

V=squeeze(Mean_matrix(rof,:,:)); % Bugfix 18 July
C=mean(mean(V))*diag(diag(V));
V=V-diag(diag(V));
V1=V+C;
figure; imagesc(V1)
colorbar
colormap('jet')
title(['Qint Matrix ',num2str(k),' mean ', num2str(round(mean(mean(V1)),3)),', max ',num2str(round(max(max(V1)),3)),', min ',num2str(round(min(min(V1)),3)) ])

Image = getframe(gcf);
imwrite(Image.cdata, sprintf('PairWiseQint.jpg'));
if size(varargin)>0
    TOP=cell2mat(varargin);
    tempo=size(TOP); 
    if tempo(1,1)==1% Make inte vertical vector if horizontal
        TOP=TOP';
    end
    save(sprintf('Top'), 'TOP','Louv','CG','V','r','-v7.3')
else
    save(sprintf('Top'), 'TOP', 'TOPS','Louv','CG','V','r','-v7.3')
end

%%%%%% Get subnetworks (SNCs and SNs) %%%%%
%Decide if seed is even or odd
xc=size(TOP)
if min(xc)==1
    evenodd=1
else
    evenodd=0
end

%%%%%%%%
if evenodd==0
    %Choose If PDFs should be calculated
    aa=input('Do you want to background probability distributions? if yes type 1 if NO type 0')
    
    if aa==1
        X=input('How many random samples?')
        %Create probability distributions:
        
        
        load Top
        
        pp=[4 6 8 10 12]
        
        for nn=1:numel(pp)
            nn
            a=0;
            
            n=pp(nn)
            randComp=zeros(X,n);
            Qint=zeros(X,S);
            MeanInt=[];
            p=0;
            
            for r=1:X
                r
                ran = round(1 + (N-1).*rand(n,1),0);
                if numel(unique(ran))==n % ~(randomA(1,1)==randomA(2,1))
                    p=p+1
                    randComp(p,:)=ran;
                    if S>1
                        parfor s=1:S
                            tempB=squeeze(Louv(s,:,:));% Time and community
                            
                            InCom=0;
                            y=tempB(:,ran);
                            z=prod(y,2); %Row with same community nr have product 1^6=1,2^6=64,3^6=729
                            if n==4
                                zz=(z==[1 16 81 ]);
                            elseif n==6
                                zz=(z==[1 64 729]);
                            elseif n==8
                                zz=(z==[1 256 6561]);
                            elseif n==10
                                zz=(z==[1 1024 59049]);
                            else
                                zz=(z==[1 4096 531441]);
                            end
                            zzz=(zz(:,1)+zz(:,2)+zz(:,3));
                            TimeOfInt(s,p,:)=zzz;
                            InCom=numel(nonzeros(zzz));
                            Qint(p,s)=InCom/Tmax;
                            
                        end
                    else
                        for s=1:S
                            tempB=squeeze(Louv(s,:,:));% Time and community
                            
                            InCom=0;
                            y=tempB(:,ran);
                            z=prod(y,2); 
                            if n==4
                                zz=(z==[1 16 81 ]);
                            elseif n==6
                                zz=(z==[1 64 729]);
                            elseif n==8
                                zz=(z==[1 256 6561]);
                            elseif n==10
                                zz=(z==[1 1024 59049]);
                            else
                                zz=(z==[1 4096 531441]);
                            end
                            zzz=(zz(:,1)+zz(:,2)+zz(:,3));
                            TimeOfInt(s,p,:)=zzz;
                            InCom=numel(nonzeros(zzz));
                            Qint(p,s)=InCom/Tmax;
                            
                        end
                    end
                    
                    MeanInt(p)=mean(Qint(p,:));
                    
                    randComp(p,:)=ran;
                    
                end
                
            end
            
            MeanPDF(nn,1)=mean(MeanInt);
            MeanPDF(nn,2)=std(MeanInt);
            
            save(sprintf('PDF_IMF_%d',n),'randComp','MeanInt','Qint','MeanPDF','TimeOfInt''-v7.3')
            
        end
        
        figure
        p=[4 6 8 10 12];
        for i=1:5
            
            load(sprintf('PDF_IMF_%d',p(i)))
            ME(i,1)=numel(nonzeros(MeanInt))/numel(MeanInt) %percent non-zero entries
            ME(i,2)=mean(MeanInt)
            ME(i,3)=std(MeanInt)
            ME(i,4)=max(MeanInt)
            ME(i,5)=min(MeanInt)
            hold on;histogram(MeanInt,50,'Normalization','probability')
            
        end
        
        title('PDFs Random components')
        legend('size 4','size 6','size 8','size 10', 'size 12')
        xlabel('qint')
        ylabel('probability')
    else
    end
end

%%%%Get SNCs
Tmax=T-10;
%SNC_size=input('Choose size of SNCs (default=8, other options 10 or 12 if initital seed was a pair, otherwise odd numbers up to 11)')

SNCs(K,Louv,V,TOP,Tmax,S,N,SNC_size,thr,evenodd) % calculatest the SNCs 


%%%% Get SNs

OV=input('Choose minimum area overlap between SNCs in SN (default=5 if SNC_size=8)')
if SNC_size==8||SNC_size==7
    load(sprintf('P%d_I8',K))
    H=R8_IDXsUfreqSort;
elseif SNC_size==6||SNC_size==5
    load(sprintf('P%d_I6',K))
    H=R6_IDXsUfreqSort;
elseif SNC_size==4||SNC_size==3
    load(sprintf('P%d_I4',K))
    H=R4_IDXsUfreqSort;
elseif SNC_size==10||SNC_size==9
    load(sprintf('P%d_I10',K))
    H=R10_IDXsUfreqSort;
elseif SNC_size==12||SNC_size==11
    load(sprintf('P%d_I12',K))
    H=R12_IDXsUfreqSort;
end

if evenodd==1
    for i=1:length(H)
        Temp(i,:)=unique(H(i,:));
    end
    H=Temp;
else
    Temp=H;
end

H=H(:,2:SNC_size+1); %Only the areas

k=0
for p=1:length(Temp(:,1))
    A=nonzeros(unique(Temp(1,2:SNC_size+1))); %Areas of first row.

    TT=zeros(length(Temp(:,1)),SNC_size);
     k=k+1
    for i=1:SNC_size
        Ttemp=zeros(length(Temp(:,1)),SNC_size);
        Ttemp(Temp(:,2:SNC_size+1)==A(i))=1;
        Hv=sum(Ttemp,2); %make into vector
        TT(:,i)=Hv; %rows that overlap with A(i)
    end
    z=sum(TT,2); %Rows with overlap
    
    r=find(z>(OV-1)); % Rows beloing to SN
    D(k,1:numel(r),:)=Temp(r,:);
    Temp(r,:)=[]; %Delet rows that are assigned to SN
    length(Temp(:,1)) %Length of reminder
    if length(Temp(:,1))<=1
        break
    else
    end
end

for i=1:length(D(:,1))
    
    bb=squeeze(D(i,:,2:SNC_size));
    tem=unique(nonzeros(bb(:,:)));
    SNs(i,1:numel(tem))=tem;
    Areas_SNs(i)=numel(tem);
    clear tem
    
end
max(Areas_SNs)
min(Areas_SNs)

%Spatial Overlap between SNs
NK=zeros(length(D(:,1)),length(D(:,1)));
for i=1:length(D(:,1))
    for n=1:length(D(:,1))
        if n>i
            M(i,n)=numel(intersect(nonzeros(SNs(i,:)),nonzeros(SNs(n,:))));
            NK(i,n)=2*numel(intersect(nonzeros(SNs(i,:)),nonzeros(SNs(n,:))))/(Areas_SNs(i)+Areas_SNs(n));
        end
    end
end
figure; imagesc(M)
figure; imagesc(NK)

save(sprintf('Subnetworks_overlap%d',OV), 'SNs', 'Areas_SNs', 'NK',  'D') % 

%%%% Get convergence of SNCs between iterations %%%%

load(sprintf('P%d_I4',K))
P(1,1)=mean(FreqIDXsU4)
P(1,2)=iqr(FreqIDXsU4)
P(1,3)=max(FreqIDXsU4)
P(1,4)=min(FreqIDXsU4)
Konvergence(1)=1-length(R4_IDXsU(:,1))/length(R4_IDXsALL(:,1))

if SNC_size >4
    load(sprintf('P%d_I6',K))
    P(2,1)=mean(FreqIDXsU6)
    P(2,2)=iqr(FreqIDXsU6)
    P(2,3)=max(FreqIDXsU6)
    P(2,4)=min(FreqIDXsU6)
    Konvergence(2)=1-length(R6_IDXsU(:,1))/length(R6_IDXsALL(:,1))
    
    if SNC_size>6
        load(sprintf('P%d_I8',K))
        
        P(3,1)=mean(FreqIDXsU8)
        P(3,2)=iqr(FreqIDXsU8)
        P(3,3)=max(FreqIDXsU8)
        P(3,4)=min(FreqIDXsU8)
        Konvergence(3)=1-length(R8_IDXsU(:,1))/length(R8_IDXsALL(:,1))
    end
end



%%%%%%

%1) Calculate SNC-timeseries and Leida Weights

load Top Louv
load(sprintf('Louvain%d_%d',j,K)) % Third IMF

load(sprintf('Subnetworks_overlap%d',OV))
load PhaseAndImf

Phase=squeeze(Phase_imf(:,:,j,6:T-5));
Tmax=T-10;
W=length(D(1,1,:))

for w=1:length(D(:,1,1))
    Qint_SNC=[];
    w
    SNCbinTS=[]; SNCcomTS =[]; SNCmeanPhase=[];
    L=length(nonzeros(D(w,:,1)))
    SNCbinTS=zeros(S,L,Tmax);
    SNCcomTS=zeros(S,L,Tmax);
    SNCmeanPhase=zeros(S,L,Tmax);
    
    
    NN=squeeze(D(w,1:L,2:W));
    
    parfor k=1:L
        
        if L>1
            snc=NN(k,:); %the SNC areas
        else
            snc=NN';
        end
        for s=1:S
            
            tempB=squeeze(Louv(s,:,:));% Time and community
            InCom=0;
            y=tempB(:,snc);
            z=prod(y,2);
	    if SNC_size==8 
                zz=(z==[1 256 6561]); %Row with same community nr have product 1^8=1,2^8=256,3^8=6561
            elseif SNC_size==7
                zz=(z==[1 128 2187]);
            elseif SNC_size==6
                zz=(z==[1 64 729]);
            elseif SNC_size==5
                zz=(z==[1 32 243]);
            elseif SNC_size==9
                zz=(z==[1 512 19683]);
            elseif SNC_size==10
                zz=(z==[1 1024 59049]);
            elseif SNC_size==11
                zz=(z==[1 2048 177147]);
            elseif SNC_size==12
                zz=(z==[1 4096 531441]);
            elseif SNC_size==4
                zz=(z==[1 16 81]);
            elseif SNC_size==3
                zz=(z==[1 8 27]);  
            end
            zzz=(zz(:,1)+zz(:,2)+zz(:,3));
            SNCbinTS(s,k,:)=zzz; %Binary timeseries
            temp3=tempB(:,snc(1)); % get the community nr of the SNC
            temp3(~logical(zzz))=0;
            SNCcomTS(s,k,:)=temp3;% timeseries with community beloning
           
            temp5=squeeze(Phase(s,:,:));
           
            SNCmeanPhase(s,k,:)=mean(temp5(snc,:),1);% mean Phase value of the SNC at each timepoint
            InCom=numel(nonzeros(zzz));
            Qint_SNC(k,s)=InCom/Tmax;
        end
        
    end
    save(sprintf('SNC8_TS_%d',w),'Qint_SNC', 'SNCbinTS','SNCcomTS','SNCmeanPhase','-v7.3')    
end

% test=squeeze(SNCbinTS(2,1,:,:));
% figure; plot(test')
% test=squeeze(SNCcomTS(2,1,:,:));
% figure; plot(test')

%%%% Calculate SN-timeseries and SN mean phase %%%%
% Check the proportion of times SNCs in the same SN are split into
% different communities

Y=length(D(:,1,1))
SNcomTS=zeros(Y,S,Tmax);
ww=zeros(Y,S,Tmax);
T2=zeros(Y,S,Tmax);
Bb=zeros(Y,S,3);

SNmeanPhase=zeros(Y,S,Tmax);

for y=1:Y
   y 
    clear Qint_SNC SNCbinTS SNCcomTS  SNCmeanPhase
    
   load(sprintf('SNC8_TS_%d',y))
    p=0;
    pp=0;
    
    for s=1:S
        for t=1:Tmax
            
            a=nonzeros(unique(SNCcomTS(s,:,t)));
            if a>0 %if any SNC belonging to the SN is integrated
                if numel(a)==1 % if all in the same community
                    p=p+1;
                    SNcomTS(y,s,t)=a; %community timeseries
                   
                    SNmeanPhase(y,s,t)=std(nonzeros(SNCmeanPhase(s,:,t)));
                    
                elseif numel(a)>1 %if split into more than one community
                    pp=pp+1;
                    ww(y,s,t)=1;
                    Bb(y,s,1:numel(a))=a;
                    if numel(a)==2
                        a1=numel((find(SNCcomTS(s,:,t)==a(1))));%nr of instance of com 1
                        a2=numel((find(SNCcomTS(s,:,t)==a(2))));%nr of instance of com 2
                        if a1>a2
                            SNcomTS(y,s,t)=a(1);
                            r=find(SNCcomTS(s,:,t)==a(1));
                            T2(y,s,t)=a2/(a1+a2);
                            SNmeanPhase(y,s,t)=std(nonzeros(SNCmeanPhase(s,r,t)));
                            
                        else
                            SNcomTS(y,s,t)=a(2);
                            r=find(SNCcomTS(s,:,t)==a(2));
                            T2(y,s,t)=a1/(a1+a2);
                           SNmeanPhase(y,s,t)=std(nonzeros(SNCmeanPhase(s,r,t)));
                            
                        end
                        
                    else
                        a1=numel((find(SNCcomTS(s,:,t)==a(1))));
                        a2=numel((find(SNCcomTS(s,:,t)==a(2))));
                        a3=numel((find(SNCcomTS(s,:,t)==a(3))));
                        if a1>a2 && a1>a3
                            SNcomTS(y,s,t)=a(1);
                            r=find(SNCcomTS(s,:,t)==a(1));
                            T2(y,s,t)=(a2/(a1+a2)+a3/(a1+a3))/2;
                            SNmeanPhase(y,s,t)=std(nonzeros(SNCmeanPhase(s,r,t)));
                        elseif a2>a1 && a2>a3
                            SNcomTS(y,s,t)=a(2);
                            r=find(SNCcomTS(s,:,t)==a(2));
                            T2(y,s,t)=(a1/(a1+a2)+a3/(a3+a2))/2;
                            SNmeanPhase(y,s,t)=std(nonzeros(SNCmeanPhase(s,r,t)));
                        else
                            SNcomTS(y,s,t)=a(3);
                            r=find(SNCcomTS(s,:,t)==a(3));
                            T2(y,s,t)=(a2/(a2+a3)+a1/(a3+a1))/2;%Proportion of clusters in the non dominant community
                            SNmeanPhase(y,s,t)=std(nonzeros(SNCmeanPhase(s,r,t)));
                        end
                    end
                end
            end
        end
        QintSN(y,s)=numel(find(SNcomTS(y,s,:)>0))/Tmax;
        
    end
    MeanQintSN(y,1)=mean(QintSN(y,:));
    MeanQintSN(y,2)=std(QintSN(y,:));
   
    CommunitySplit(y)=pp/(p+pp);
end
save SN_Time T2 CommunitySplit QintSN MeanQintSN SNmeanPhase SNcomTS ww Bb

% Plot example of all SN_timeseries for subject 1
temp=squeeze(SNcomTS(:,1,:));
figure;
for i=1:Y
    plot(temp(i,:)+i)
    hold on;
end

%%%Get mean IMF activation

load PhaseAndImf.mat
if band==0
[S N T J]=size(BOLD_imf)

for s=1:S
    
        imf=squeeze(BOLD_imf(s,:,6:T-5,j));
        for a=1:N
            IMFdm(s,a,:)=imf(a,:)-mean(imf(a,:)); %normalize imf-timeseries
            amp(s,a,:)=IMFdm(s,a,:)/max(abs(IMFdm(s,a,:))); %normalize against max amplitude within imf
            %amp(s,a,:)=IMFdm(s,a,:); %only demean
        end
   
end

else 
    
 [S N T J]=size(BOLD_band)

for s=1:S
  
        imf=squeeze(BOLD_band(s,:,6:T-5));
        for a=1:N
            IMFdm(s,a,:)=imf(a,:)-mean(imf(a,:)); %normalize imf-timeseries
            amp(s,a,:)=IMFdm(s,a,:)/max(abs(IMFdm(s,a,:))); %normalize against max amplitude within imf
            %amp(s,a,:)=IMFdm(s,a,:); %only demean
        end
    
end   
    
end

save Activation amp



%%%%Activation/deactivation of each area in SNs as well as full SN %%%%

load(sprintf('SN_Time'))
load(sprintf('Subnetworks_overlap%d',OV))
load Activation

E=Areas_SNs; 
Y=numel(Areas_SNs)
TimeArea=zeros(Y,S,Tmax,max(E));
TimeAreaProp=zeros(Y,S,Tmax);
mean_amp=zeros(Y,S,Tmax);

W=length(D(1,1,:));
for y=1:Y
    
    clear  Qint_SN SNCbinTS SNCcomTS DD F e DomCom
    load([sprintf('SNC8_TS_%d',y)])
    
    DD=squeeze(D(y,:,2:W));
    l=length(nonzeros(DD(:,1)));
    DD=DD(1:l,:);
    
    
    for s=1:S
        
        F=squeeze(SNCcomTS(s,:,:));
        if l==1
           
            F=F';
        end
        ImfA=squeeze(amp(s,:,:)); %get relativea ctivation/deactivation of areas
        
        for t=1:Tmax
            
            if SNcomTS(y,s,t)>0 % if SN is assembled
                [r c]=find(F(:,t)==SNcomTS(y,s,t));
                u=unique(DD(r,:)); %Get the unique areas
                
                TimeArea(y,s,t,1:numel(u))=u; %unique areas at the timepoint
                TimeAreaC(y,s,t,1:numel(u))=SNcomTS(y,s,t); %the community
                TimeAreaProp(y,s,t)=numel(u)/E(y);
                mean_amp(y,s,t)=mean(ImfA(u,t));
                
            end
        end
    end
end

mean_amp(isnan(mean_amp))=0;

TimeAreaProp(isnan(TimeAreaProp))=0;
TimeArea(isnan(TimeArea))=0;
TimeAreaC(isnan(TimeAreaC))=0;

save Time_normalized_IMF_C  mean_amp  TimeArea TimeAreaC  TimeAreaProp 




%%%%Integration and segregation of SNs relative to a) time of overlap b) all
%time
[Y S Tmax]=size(SNcomTS)
sync=zeros(S,Y,Y);
overlap=zeros(S,Y,Y);
for s=1:S
    temp=squeeze(SNcomTS(:,s,:));
    ex1=(temp==1);%in same communities
    ex2=(temp==2);
    ex3=(temp==3);
    temp(temp>0)=1;
    for i=1:Tmax
        sync(s,ex1(:,i),ex1(:,i))=sync(s,ex1(:,i),ex1(:,i))+1;
        sync(s,ex2(:,i),ex2(:,i))=sync(s,ex2(:,i),ex2(:,i))+1;
        sync(s,ex3(:,i),ex3(:,i))=sync(s,ex3(:,i),ex3(:,i))+1;
        overlap(s,logical(temp(:,i)),logical(temp(:,i)))=overlap(s,logical(temp(:,i)),logical(temp(:,i)))+1;
    end
    
    syncALLtime(s,:,:)=sync(s,:,:)/Tmax;
    syncOver(s,:,:)=sync(s,:,:)./overlap(s,:,:);
    clear temp
end
sync(isnan(sync))=0;
syncALLtime(isnan(syncALLtime))=0;
overlap(isnan(overlap))=0;
if S>1
    SynC=mean(syncALLtime);
    SynC=squeeze(SynC(1,:,:));
    syncOver(isnan(syncOver))=0;
    SynCover=mean(syncOver,1);
    SynCover=squeeze(SynCover(1,:,:));
else
    SynC=squeeze(syncALLtime(1,:,:));
    syncOver(isnan(syncOver))=0;
    SynCover=squeeze(syncOver(1,:,:));
end
save SYNC SynC SynCover syncALLtime sync overlap

figure; imagesc(SynC)
title('Sync relative to all time')
colorbar
colormap(copper)
figure; imagesc(SynCover)
title('Sync relative to time of simultaneous integration')
colorbar
colormap(copper)


%%%% Find hierachical cluster solution that maximizes integration within cluster and segregation between MNs

load SYNC

F=ones(Y,Y)-SynCover;%Distance measure
Z=linkage(F,'single')
NumC=round(Y-1); %Number of cluster solutions to try
for i=1:NumC
    
    c = cluster(Z,'Maxclust',i);
    
    M=zeros(i,Y);
    HierOrd=[];HierOrdNum=[];
    q=0; p=0;
    for k=1:i %for each of clusters in the solution
        
        r=find(c==k); %Find the MNs belonging to the cluster
        M(k,1:numel(r))=r;
        
        E=[];
        if numel(r)>1 % if mor than 1 MN in the cluster
            for a=1:numel(r) %calculate mean integration within the MN ie between the SN that consitues it
                for b=1:numel(r)
                    if b>a
                        E(a,b)=SynCover(r(a),r(b));
                    else
                    end
                end
            end
            MeanINT(i,k)=mean(nonzeros(E)); %Mean integration within MN k for total number of clusters i
            MinINT(i,k)=min(nonzeros(E)); %Min integration within MN k for total number of clusters i 
            Diff(i,k)=mean(nonzeros(E))-min(nonzeros(E)); 
        else
            q=q+1;
            MeanINT(i,k)=0; %if only one SN in the MC set to zeros to try to maximize to solutions with few single clusters
            MinINT(i,k)=0;
            Diff(i,k)=0;
        end
        % get degree o sync with other SNs outside this MN
        r2=setdiff(1:Y,r);
        F=[];
        if numel(r2)>1
            for a=1:numel(r2) %calculate mean segregation with non-MN SNs
                for b=1:numel(r)
                    F(a,b)=1-SynCover(r2(a),r(b));
                end
            end
            MNseg(i,k)=mean(nonzeros(F));
            MinSeg(i,k)=min(nonzeros(F));
        else
            p=p+1;
            MNseg(i,k)=0; %if only one SN in the MC set to zeros to try to maximize to solutions with few single clusters
            MinSeg(i,k)=0;
        end  
        
    end
    %M=sortrows(M,1,'ascend'); %Sort rows such that row with SN1 is first
    HierOrd(:,1)=nonzeros(M');
    NK=[];
    for e=1:length(M(:,1))
        NK(e,1:numel(nonzeros(M(e,:))))= ones(numel(nonzeros(M(e,:))),1)*e;
    end
    HierOrd(:,2)=nonzeros(NK'); % Column 1 contains SNs, column 2 contains cluster nr
    CLUSTERsummary(i,10)=mean(nonzeros(Diff(i,:)));
    CLUSTERsummary(i,5)=mean(nonzeros(MinINT(i,:)));%Mean of minimum integration with MNs
    CLUSTERsummary(i,6)=mean(nonzeros(MinSeg(i,:)));%Mean of min segretation With other MNs
    CLUSTERsummary(i,1)=mean(nonzeros(MeanINT(i,:)));%Mean of Mean integration within MN, Consider only the clusters with more than one SN
    CLUSTERsummary(i,2)=std(nonzeros(MeanINT(i,:)));
    CLUSTERsummary(i,3)=max(nonzeros(MeanINT(i,:)));% max mean integration among clusters
    CLUSTERsummary(i,4)=min(nonzeros(MeanINT(i,:)));% minimum mean integration within composed MNs per cluster solution
    CLUSTERsummary(i,7)= q;%number of single clusters
    if i>1
        CLUSTERsummary(i,8)=mean(nonzeros(MNseg(i,:)))% mean segregation with other clusters
        CLUSTERsummary(i,9)=min(nonzeros(MNseg(i,:)))% Min of min segregation with other SNs
    else
        CLUSTERsummary(i,8)=0;
        CLUSTERsummary(i,9)=0;
    end
  
    CLUSTERsummary(i,7)= q;%number of single clusters
    
%     figure;imagesc(SynCover(HierOrd(:,1),HierOrd(:,1)))
%     set(gca,'FontSize',10)
%     set(gca,'FontWeight','bold')
%     colorbar
%     colormap('spring')
%     xticks([1:1:Y])
%     xticklabels(HierOrd(:,2))
%     yticks([1:1:Y])
%     yticklabels(HierOrd(:,2))
    
%     title({['HierC Ver ',num2str(i),'. Mean of mean int within MNs', num2str(round(CLUSTERsummary(i,1),3)),' (sd ',num2str(round(CLUSTERsummary(i,2),3)),')']
%         ['Mean of Min int within MNs ',num2str(round(CLUSTERsummary(i,4),3)) ] },'Fontsize',15)
%     Image = getframe(gcf);
%     imwrite(Image.cdata, sprintf('ClustVer_%d.jpg',i));
    save(sprintf('HierClustV_%d',i),'M','HierOrd','HierOrdNum')
end

figure; plot(CLUSTERsummary(:,1),'*','Linewidth',2); hold on; plot(CLUSTERsummary(:,5),'*','Linewidth',2)
hold on;  plot(CLUSTERsummary(:,8),'*','Linewidth',2);
title('Mean and minimum integration within MetaNetwork')
xlabel('Number of Clusters')

legend('Mean of mean integration', 'Min of Mean integration','Mean of segregation other SNs')

% Prompt to choose cluster version
 x = input('Standard cluster thresholds (Mean integration within MNs >=0.95 and Min integrattion within MNs >=0.80) (0) or set different thresholds (1)?')

% Atomatic selection 
if x==0
% Find first cluster solution with Mean of Mean integration>=0.95 and Min
% of Mean integrattion >=0.80
a1=0.95;
b1=0.80
ra=find(CLUSTERsummary(:,1)>=a1); 
rb=find(CLUSTERsummary(:,4)>=b1);

if numel(ra) && numel(rb)>0
    if min(ra)<min(rb)
        C=min(rb)
    else
        C=min(ra)
    end
else
    C=NumC;
    sprintf('Cluster thresholds were not met, cluster with higest number choosen')
end

else 
    a1=input('Mean of mean integration within MNs')
    b1=input('Min of mean integration within MNs')
    ra=find(CLUSTERsummary(:,1)>=a1);
    rb=find(CLUSTERsummary(:,4)>=b1);

if numel(ra) && numel(rb)>0
    if min(ra)<min(rb)
        C=min(rb)
    else
        C=min(ra)
    end
else
    C=NumC;
    sprintf('Cluster thresholds were not met, cluster with higest number choosen')
end
    
end
save ClusterSUM CLUSTERsummary MeanINT MinINT C a1 b1


 %%%%Meta network timeseries%%%%

load(sprintf('SN_Time'))
load(sprintf('Subnetworks_overlap%d',OV))
load(sprintf('HierClustV_%d',C))

SNcomTS(SNcomTS>0)=1; %binarize ts

MN=zeros(C,max(Areas_SNs),S,Tmax);
for i=1:C
    
    m=nonzeros(M(i,:));
    nMN(i)=numel(m);
    MN(i,1:numel(m),:,:)=SNcomTS(m,:,:);
    MN_TS(i,:,:)=squeeze(sum(MN(i,:,:,:),2));
    
end
MN_TS(MN_TS>0)=1;
if S>1
    QintMNsubj=mean(MN_TS,3);
    MeanQintMN(:,1)=mean(mean(MN_TS,3),2) ;
    mm=mean(MN_TS,3);
    MeanQintMN(:,2)=std(mm') ;
else
    QintMNsubj=mean(MN_TS,2);
    MeanQintMN(:,1)=mean(MN_TS,2) ;
    MeanQintMN(:,2)=0 ;  %only one subj
end
save(sprintf('MNQint_%d',C), 'MN_TS', 'QintMNsubj', 'MeanQintMN') %MNs in unsorted order


%Sort MNs from most frequent to least

[MNsort MN_I]=sortrows(MeanQintMN',1,'descend')


% Sort HierOrd according to Qint of the MNs
[MeanQintMNSorted, Isort] =sortrows(MeanQintMN,1,'descend');

Hsort=[];
for i=1:C
    r=find(HierOrd(:,2)==Isort(i))
    Hsort=[Hsort;(HierOrd(r,:))]
end

save SortedHierOrd  Hsort MeanQintMNSorted Isort



%%%%%%%%%%%%%%
%%%%Get Meta cluster maps weighted by probability of occurence per area%%%%
load SN_Time 
load(sprintf('Subnetworks_overlap%d',OV))
load(sprintf('HierClustV_%d',C))  % Sorted
load Time_normalized_IMF_C
[SN S Tmax Ar]=size(TimeArea)
TimeAreaMetaW=zeros(S,C,N);
SNcomTS(SNcomTS>0)=1;

for s=1:S
    Tem=squeeze(TimeArea(:,s,:,:));
    for i=1:C
        AA=zeros(N,1);
        Temp=squeeze(Tem(HierOrd(HierOrd(:,2)==i,1),:,:));
        if size(Temp,3)>1 %more than one SN
            for t=1:Tmax
                [aa bb]=hist(nonzeros(Temp(:,t,:)),unique(nonzeros(Temp(:,t,:))));
                AA(bb)=AA(bb)+1;
            end
        else
            for t=1:Tmax
                [aa bb]=hist(nonzeros(Temp(t,:)),unique(nonzeros(Temp(t,:))));
                AA(bb)=AA(bb)+1;
            end
        end
        TimeAreaMetaW(s,i,1:N)=AA/Tmax;%reltive freq (all time) per area,metanetwork and subj
        B=squeeze(SNcomTS(HierOrd(HierOrd(:,2)==i,1),s,:)); % Timecourses of SNs belonging to i:th MN 
        if size(Temp,3)>1
            B=sum(B(:,:)); %Add timeseries of the SNs belongning to the MNs
            B(B>0)=1; %Common TS for MN
        end
        MeanQintMNsubj(s,i)=sum(B)/Tmax; %relative freq of MN in relation to all time
        TimeAreaMetaInt(s,i,1:N)=AA/sum(B);% Relative freq of Areas in relation to number of timepoints of MN integration
    end
end

TimeAreaMetaW(isnan(TimeAreaMetaW))=0;
TimeAreaMetaInt(isnan(TimeAreaMetaInt))=0;

if S>1
    for i=1:C
        MetaMeanArea(i,:)=mean(TimeAreaMetaW(:,i,:));%relative to all time
        MetaMeanInt(i,:)=mean(TimeAreaMetaInt(:,i,:));
    end
else
    MetaMeanArea=squeeze(TimeAreaMetaW(1,:,:));%relative to all time
    MetaMeanInt=squeeze(TimeAreaMetaInt(1,:,:));
end

save MetaProbMaps MetaMeanArea MetaMeanInt TimeAreaMetaInt TimeAreaMetaW MeanQintMNsubj %MeanMetaQint%Int is relative to networks own time of assembly

% SNs belonging to each MN
for i=1:C
    temp=HierOrd(HierOrd(:,2)==i,1);
    SNnrMN(i,1:numel(temp))=temp;
end

% Areas beloning to each MN
MNmap=zeros(length(SNnrMN(:,1)),N);
for i=1:length(SNnrMN(:,1))
    MNmap(i,1:numel(nonzeros(unique(SNs(nonzeros(SNnrMN(i,:)),:)))))=nonzeros(unique(SNs(nonzeros(SNnrMN(i,:)),:)));
end

%Sort MNmap accordint to Qint

MNmapSort=MNmap(Isort,:);
SNnrMNsort=(SNnrMN(Isort,:));
  
save MN_MAP MNmap SNnrMN SNnrMNsort MNmapSort %MeanMetaQintSorted



% Make MN maps

setenv('FSLOUTPUTTYPE','NIFTI_GZ')
getenv('FSLOUTPUTTYPE')

Weight=MetaMeanInt(Isort,:);

%Set map scale
maxWeight=round(max(100*Weight(:)));
minWeight=round(min(100*Weight(:)));
load MaxVox
EC.vol.px=maxWeight;
EC.vol.pn=minWeight;
save MaxVox EC


MNC=C
MN=MNmapSort;
for y=1:length(MN(:,1))
    AP=nonzeros(Weight(y,:))*100;
    A=nonzeros(MN(y,:))
    a1=A(1);
    cs = sprintf('cp /Users/marist1/Documents/MATLAB/Parcellation236/%dSH_7mask.nii.gz MN%d_OV%d_%d.nii.gz',a1,MNC,OV,y);
    system(cs)
    %weigthed:
    cmd=sprintf('/usr/local/fsl/bin/fslmaths  MN%d_OV%d_%d.nii.gz -mul %d  MN%d_OV%d_%d.nii.gz ',MNC,OV,y,AP(1),MNC,OV,y);
    system(cmd)
    
    for a=2:numel(A)
        
        n=A(a);
        %with weight
        cmd = sprintf('/usr/local/fsl/bin/fslmaths  /Users/marist1/Documents/MATLAB/Parcellation236/%dSH_7mask.nii.gz -mul %d temp_%d.nii.gz',n,AP(a),a);
        system(cmd)
        
        cmd = sprintf('/usr/local/fsl/bin/fslmaths MN%d_OV%d_%d.nii.gz -add temp_%d.nii.gz MN%d_OV%d_%d.nii.gz',MNC,OV,y,a,MNC,OV,y);
        
        system(cmd);
        delete(sprintf('temp_%d.nii.gz',a))
    end
    
    BrainNet_MapCfg('BrainMesh_ICBM152.nv',sprintf('MN%d_OV%d_%d.nii.gz',MNC,OV,y),'MAxVox.mat',sprintf('MN%d_OV%d_%d.jpg',MNC,OV,y));
  
    clear r A
    y
end



close all

v = VideoWriter(sprintf('MNs_K%d_OV%d_%d_ALL.avi',K,OV,C));
v.FrameRate=1.5

for i=1:length(MNmap(:,1))
    imshow(sprintf('MN%d_OV%d_%d.jpg',C,OV,i),'InitialMagnification','fit')
    title(['MN ', num2str(i), ', Qint ',num2str(round(MeanQintMNSorted(i,1),2))], 'FontSize',30, 'FontWeight', 'bold')
  
    frame=getframe(gcf);
    im=frame2im(frame);
    imwrite(im,sprintf('MN_%d.jpg',i))
    F = getframe(gcf) ;
    open(v)
    writeVideo(v,F)
    close(gcf)
    
end

 


%%%%Find most dominant Network per MN to name them accordingly

%Get how Shaefer parcell Nets contribute to different MNs
load(sprintf('HierClustV_%d',C))
load(sprintf('Subnetworks_overlap%d',OV))
load SortedHierord Isort
MNarea=zeros(C,N);
for i=1:C
    
    r=find(HierOrd(:,2)==i)
    aa=numel(unique(nonzeros(SNs(HierOrd(r,1),:))))
    MNarea(i,1:aa)=unique(nonzeros(SNs(HierOrd(r,1),:)));
    MNaCount(i)=aa; %areas per MN
    
end

MNarea=MNarea(:,1:max(MNaCount));
MNareaSort=MNarea(Isort,:);
MNaCountSort=MNaCount(Isort);

save MN_AREAS MNarea MNaCount MNareaSort MNaCountSort

%Number of areas per MN
%Weights contributed from each MN
%MN number, total nr of areas of the MN, percent of total netw areas,
%weight for relative activation
load MN_AREAS
load NetName
load MetaProbMaps
load AreaParcelA236_N9
load SortedHierOrd
Nettot=table2array(NetName(:,3));
for i = 1:C
    clear k k2 w w2 w3 re re2 re3
    nets=AreaNetw(nonzeros(MNarea(i,:)),1); %Network beleoning of the areas
    netsWe=MetaMeanInt(i,nonzeros(MNarea(i,:))); %the weighting of the areas in the MN
    
    
    netu=unique(nets);% networks
    c=histc(nets,netu); % count area contribution per networks
    k=[netu c c./Nettot(netu)]; % uniq nets, count, rel freq of total of that network
    
    % how much each network contribute to the MN in terms of weight
    for p=1:numel(netu)
        row=find(nets==netu(p));
        w(p)=sum(netsWe(row));
    end
    
    
    MNclass{i}=[k w']
    re=find(w==max(w))
    if numel(re)==1 %if one unique max value
        MN_Net_Nr(i,1)=netu(re);% Number of largest contributor
        
        w2=setdiff(w,max(w)); %identify second largest weight
        if  max(w2)/max(w)>1/3
            re2=find(w==max(w2));
            if numel(re2)==1
                MN_Net_Nr(i,2)=netu(re2);
            else % in small networks weights can be whole numer and equal for more than one netw
                MN_Net_Nr(i,2:3)=netu(re2);
            end
            
            if numel(re2)==1
                w3=setdiff(w2,max(w2)); %identify second largest weight
                if max(w3)/max(w)>1/3
                    re3=find(w==max(w3))
                    MN_Net_Nr(i,3)=netu(re3);
                else
                    MN_Net_Nr(i,3)=0;
                end
            else
            end
        else
            MN_Net_Nr(i,2)=0;
        end
    elseif numel(re)==2
        MN_Net_Nr(i,1:2)=netu(re);% Number of largest contributor
        w2=setdiff(w,max(w)); %identify second largest weight
        if  max(w2)/max(w)>1/3
            re2=find(w==max(w2));
            if numel(re2)==1
                MN_Net_Nr(i,3)=netu(re2);
            else % in small networks weights can be whole numer and equal for more than one netw
                MN_Net_Nr(i,3:4)=netu(re2);
            end
            
            if numel(re2)==1 % Add more only if there is only 2 names before
                w3=setdiff(w2,max(w2)); %identifythird largest weight
                if max(w3)/max(w)>1/3
                    re3=find(w==max(w3))
                    if numel(re3)==1
                        MN_Net_Nr(i,4)=netu(re3);
                    else
                        MN_Net_Nr(i,4:5)=netu(re3);
                    end
                end
            else
            end
        else
            MN_Net_Nr(i,3)=0;
        end
    else
        MN_Net_Nr(i,1:3)=netu(re);
        
    end
end

MNclassSorted=MNclass(Isort);
MN_Net_Nr_Sorted=MN_Net_Nr(Isort,:);

for i=1:C
    temp=MN_Net_Nr_Sorted(i,:)
   a= numel(nonzeros(temp))
   for k=1:a
       MNnamesSorted(i,k)=NetName(temp(k),2)
   end
end


save MNnaming MNclass MNnamesSorted MN_Net_Nr MN_Net_Nr_Sorted MNclassSorted

% Calculate all RSNs that contribute to each MN
 load MN_AREAS
 load AreaParcelA236_N9
 [C y ]=size(MNareaSort)
 R=zeros(C,y);
 for p=1:9
   a=find(AreaNetw(:,1)==p); %areas per network
   temp=[];
   for i=1:length(a)
      
           [r c]= find(MNareaSort==a(i));
           for j=1:numel(r)
               MN_RSN(r(j),c(j))=p; % replace by netwnr
           end
       end
  
 end

 for i=1:C
     temp=unique(nonzeros(MN_RSN(i,:)));
     MN_RSNunique(i,1:numel(temp))=temp;
 end
 
 save RSNsPerMN MN_RSN MN_RSNunique
 

%%%% Weigted SNs maps %%%%
load SN_Time
load Time_normalized_IMF_C

%Binarize SNcomTS
SNcomTS(SNcomTS>0)=1;
[SN S Tmax Ar]=size(TimeArea);
for s=1:S
    Tem=squeeze(TimeArea(:,s,:,:));
    for i=1:SN
        AA=zeros(N,1);
        Temp=squeeze(Tem(i,:,:));
        
        for t=1:Tmax
            [aa bb]=hist(nonzeros(Temp(t,:)),unique(nonzeros(Temp(t,:))));
            AA(bb)=AA(bb)+1;
        end
        
        TimeAreaSNW(s,i,1:N)=AA/Tmax;%reltive freq(all time)
        B=squeeze(SNcomTS(i,s,:)); %
        Mean_SNTS(s,i)=sum(B)/Tmax;
        TimeAreaSNInt(s,i,1:N)=AA/sum(B);%
    end
end
TimeAreaSNW(isnan(TimeAreaSNW))=0;
TimeAreaSNInt(isnan(TimeAreaSNInt))=0;

if S>1
    for i=1:SN
        SNMean(i,:)=mean(TimeAreaSNW(:,i,:));
        SNMeanInt(i,:)=mean(TimeAreaSNInt(:,i,:));
    end
else
    SNMean=squeeze(TimeAreaSNW(1,:,:));
    SNMeanInt=squeeze(TimeAreaSNInt(1,:,:));
end
save SnWeigtMaps SNMean SNMeanInt TimeAreaSNInt TimeAreaSNW %Int is relative to networks own time of assembly


% Plot SN-maps per MN(sorted by highest Qint) and such that SN with highest
% qint within MN comes first
setenv('FSLOUTPUTTYPE','NIFTI_GZ')
getenv('FSLOUTPUTTYPE')
load SN_Time
load MN_MAP 
load SnWeigtMaps
load Subnetworks_overlap5 SNs %areas per SN

 
for f=1:C
    
    sn=numel(nonzeros(SNnrMNsort(f,:))); %number of SNs making up the MN
    temp=nonzeros(SNnrMNsort(f,:));
    temp2=MeanQintSN(temp);
    [t3 sortSN]=sort(temp2,'descend')
    SnsortedQint=temp(sortSN); %order of SN according to highest qint
    
    for y=1:sn%length(SNs(:,1))
        AP=nonzeros(SNMeanInt(SnsortedQint(y),:))*100;
        A=nonzeros(SNs(SnsortedQint(y),:)); %
        a1=A(1);
        
        cs = sprintf('cp /Users/marist1/Documents/MATLAB/Parcellation236/%dSH_7mask.nii.gz SN%d_OV%d_%d.nii.gz',a1,f,OV,y);
        system(cs)
        %%weigthed:
        cmd=sprintf('/usr/local/fsl/bin/fslmaths SN%d_OV%d_%d.nii.gz -mul %d SN%d_OV%d_%d.nii.gz ',f,OV,y,AP(1),f,OV,y);
        system(cmd)
        
        for a=2:numel(A)
            n=A(a);
            %with weight
            cmd = sprintf('/usr/local/fsl/bin/fslmaths  /Users/marist1/Documents/MATLAB/Parcellation236/%dSH_7mask.nii.gz -mul %d temp_%d.nii.gz',n,AP(a),a);
            system(cmd)
            cmd = sprintf('/usr/local/fsl/bin/fslmaths SN%d_OV%d_%d.nii.gz -add temp_%d.nii.gz SN%d_OV%d_%d.nii.gz',f,OV,y,a,f,OV,y);
            
            system(cmd);
            delete(sprintf('temp_%d.nii.gz',a))
        end
        BrainNet_MapCfg('BrainMesh_ICBM152.nv',sprintf('SN%d_OV%d_%d.nii.gz',f,OV,y),'MaxVox.mat',sprintf('SN%d_OV%d_%d.jpg',f,OV,y));
        clear r A AP
        
    end
    
end


close all

v = VideoWriter(sprintf('SNsPerMN_K%d_OV%d_%d_ALL.avi',K,OV,C));
v.FrameRate=1.5

for i=1:length(MNmap(:,1))
    imshow(sprintf('MN%d_OV%d_%d.jpg',C,OV,i),'InitialMagnification','fit')
    title(['MN ', num2str(i), ', Qint ',num2str(round(MeanQintMNSorted(i,1),2))], 'FontSize',30, 'FontWeight', 'bold')
   
    frame=getframe(gcf);
    im=frame2im(frame);
    imwrite(im,sprintf('MNS_%d.jpg',i))
    F = getframe(gcf) ;
    open(v)
    writeVideo(v,F)
    close(gcf)
    
    %Each of the SNs making up the i:th MN
    sn=numel(nonzeros(SNnrMNsort(i,:))); %number of SNs making up the MN
    tempA=nonzeros(SNnrMNsort(i,:)); % The SNs numbered in original order that make up i:th MN
    tempB=MeanQintSN(tempA); % The qint of SNs in TempA
    [t3 sortSN]=sort(tempB,'descend'); %
    tempC=tempA(sortSN); % SNs in i:th MN in sorted order according to qint
    SnsortedQint=tempB(sortSN); %Qint value of SNs in sorted order
    
    for y=1:sn %length(SNs(:,1))
        AP=nonzeros(SNMeanInt(tempC(y),:))*100;
        A=nonzeros(SNs(tempC(y),:)); %
        a1=A(1);
        
        cs = sprintf('cp /Users/marist1/Documents/MATLAB/Parcellation236/%dSH_7mask.nii.gz MN%d_SN%d_OV%d.nii.gz',a1,i,y,OV);
        system(cs)
        %%weigthed:
        cmd=sprintf('/usr/local/fsl/bin/fslmaths MN%d_SN%d_OV%d.nii.gz -mul %d MN%d_SN%d_OV%d.nii.gz ',i,y,OV,AP(1),i,y,OV);
        system(cmd)
        
        for a=2:numel(A)
            n=A(a);
            %with weight
            cmd = sprintf('/usr/local/fsl/bin/fslmaths  /Users/marist1/Documents/MATLAB/Parcellation236/%dSH_7mask.nii.gz -mul %d temp_%d.nii.gz',n,AP(a),a);
            system(cmd)
            cmd = sprintf('/usr/local/fsl/bin/fslmaths MN%d_SN%d_OV%d.nii.gz -add temp_%d.nii.gz MN%d_SN%d_OV%d.nii.gz',i,y,OV,a,i,y,OV);
            
            system(cmd);
            delete(sprintf('temp_%d.nii.gz',a))
        end
        BrainNet_MapCfg('BrainMesh_ICBM152.nv',sprintf('MN%d_SN%d_OV%d.nii.gz',i,y,OV),'MaxVox.mat',sprintf('MN%d_SN%d_OV%d.jpg',i,y,OV));
        clear r A AP
        
        close all
        imshow(sprintf('MN%d_SN%d_OV%d.jpg',i,y,OV),'InitialMagnification','fit')
        title(['SN ', num2str(y), ' (MN ',num2str(i), ') , Qint ',num2str(round(SnsortedQint(y),2))], 'FontSize',30, 'FontWeight', 'bold')
        
        frame=getframe(gcf);
        im=frame2im(frame);
        imwrite(im,sprintf('MN_%d_SN_%d.jpg',i,y))
        F = getframe(gcf) ;
        open(v)
        writeVideo(v,F)
        close(gcf)
        
    end
end



% Sync and SyncOv sorted
load SortedHierOrd
load SYNC

figure; imagesc(SynCover(Hsort(:,1),Hsort(:,1)))
title('Sync relative to time of simultaneous integration, sorted according to Qint of MN')
colorbar
colormap(copper)
figure; imagesc(SynC(Hsort(:,1),Hsort(:,1)))
title('Sync relative to all time, sorted according to Qint of MN')
colorbar
colormap(copper)


%%%% State vector recurrence plots %%%%
%Relative activation/deactivation of MNs
load(sprintf('HierClustV_%d',C),'HierOrd')
load Time_normalized_IMF_C mean_amp %SN mean IMF activation
load SortedHierOrd 
u=unique(HierOrd(:,2));
[nc l]=histc(HierOrd(:,2), u); %SNs per MN
meanIMFC=mean_amp(HierOrd(:,1),:,:); %sort in MN order

for s=1:S
    
    meanIMF=squeeze(meanIMFC(:,s,:));
    for t=1:Tmax
        r=1;
        for n=1:C
            StateNonZero(n,s,t)=mean(nonzeros(meanIMF(r:r+nc(n)-1,t))); %Mean of activation of SNs beloning to the MN
            r=r+nc(n);
        end
    end
end

StateNonZero(isnan(StateNonZero))=0;
StateNonZeroSorted=StateNonZero(Isort,:,:);
save StateIMFnonz StateNonZero StateNonZeroSorted

%example one subject
MNstate=squeeze(StateNonZeroSorted(:,3,:));
figure;imagesc(MNstate)

title('Relative activation/deactivation of MNs accross time for one subject')



%Plot an example first MNs

example=squeeze(StateNonZeroSorted(:,3,:));
figure; hold on
t=10

for i=1:t
    
    plot(example(i,:)+t,'LineWidth',1)
    t=t-1;
end



%%Create recurrence plots of statevectors for all subjects and save as
%%video

load StateIMFnonz 
load Activation
v = VideoWriter('StatesNew.avi');
v.FrameRate=1.5

t=0;

for s=1:S
    t=t+1;
    %If blank space between subjects:
%     m=ones(Tmax,Tmax);
%     
%     figure; imagesc(m)
%     colorbar
%     colormap cool
%     F =getframe(gcf) ;
%     open(v)
%     writeVideo(v,F)
%     close(gcf)
    t=t+1;
  
    MNstate=squeeze(StateNonZeroSorted(:,s,:)); %MN state vector
    SNstate=squeeze(mean_amp(:,s,:)); %SN state vector
    rawIMF=squeeze(amp(s,:,:)); %raw IMF normalized values per area
    
    
    MN=[];SN=[]; Raw=[];
    for m=1:Tmax
        for n=1:Tmax
            MN(m,n)=corr2(MNstate(:,m),MNstate(:,n));
            SN(m,n)=corr2(SNstate(:,m),SNstate(:,n));
            Raw(m,n)=corr2(rawIMF(:,m),rawIMF(:,n));
            
        end
    end
    
    
    % RawIMF state vector
    Raw=Raw-diag(diag(Raw));
    Raw(isnan(Raw))=0;
    MeanStateRaw(s)=mean(nonzeros(Raw));
    Mtrn=triu(Raw,1);
    Mt2rn=triu(Raw,2);
    MDrn=Mtrn-Mt2rn;
    MeanStateChangeRaw(s,1)=mean(nonzeros(MDrn));
    MeanStateChangeRaw(s,2)=std(nonzeros(MDrn));
    StatesAllRaw(s,:,:)=Raw;
    MaxminRaw(s,1)=min(nonzeros(Raw));
    MaxminRaw(s,2)=max(nonzeros(Raw));
    
    figure;
    imagesc(Raw)
    caxis([-1 1])
    colorbar
    colormap cool
    xticks(250:250:1000)
    xticklabels({'3','6','9','12'})
    xlabel('Time (min)')
    yticks(250:250:1000)
    yticklabels({'3','6','9','12'})
    ylabel('Time (min)')
    
    set(gca,'FontSize',20,'FontWeight','bold')
    title(['Subject Raw ',num2str(s) ])
   
    Image = getframe(gcf);
    imwrite(Image.cdata, sprintf('RawIMFstate_%d.jpg',s));
    saveas(gcf,sprintf('RawIMFstate_%d.fig',s))
    %load(sprintf('RawIMFstate_%d.jpg',s))
    
    F = getframe(gcf) ;
    open(v)
    writeVideo(v,F)
    close(gcf)
    
    % SN-state vectors
    SN=SN-diag(diag(SN));
    SN(isnan(SN))=0;
    MeanStateSN(s)=mean(nonzeros(SN));
    Mtsn=triu(SN,1);
    Mt2sn=triu(SN,2);
    MDsn=Mtsn-Mt2sn;
    MeanStateChangeSN(s,1)=mean(nonzeros(MDsn));
    MeanStateChangeSN(s,2)=std(nonzeros(MDsn));
    StatesAllSN(s,:,:)=SN;
    MaxminSN(s,1)=min(nonzeros(SN));
    MaxminSN(s,2)=max(nonzeros(SN));
    
    figure;
    imagesc(SN)
    caxis([-1 1])
    colorbar
    colormap cool
    xticks(250:250:1000)
    xticklabels({'3','6','9','12'})
    xlabel('Time (min)')
    yticks(250:250:1000)
    yticklabels({'3','6','9','12'})
    ylabel('Time (min)')
    set(gca,'FontSize',20,'FontWeight','bold')
    title(['Subject SN ',num2str(s) ])
    Image = getframe(gcf);
    imwrite(Image.cdata, sprintf('SNstate_%d.jpg',s));
    saveas(gcf,sprintf('SNstate_%d.fig',s))
    %load(sprintf('SNstate_%d.jpg',s))
    F = getframe(gcf) ;
    open(v)
    writeVideo(v,F)
    close(gcf)
    
    %MN state vectors
    MN=MN-diag(diag(MN));
    MN(isnan(MN))=0;
    MeanStateMN(s)=mean(nonzeros(MN));
    Mtmn=triu(MN,1);
    Mt2mn=triu(MN,2);
    MDmn=Mtmn-Mt2mn;
    MeanStateChangeMN(s,1)=mean(nonzeros(MDmn));
    MeanStateChangeMN(s,2)=std(nonzeros(MDmn));
    StatesAllMN(s,:,:)=MN;
   
    MaxminMN(s,1)=min(nonzeros(MN));
    MaxminMN(s,2)=max(nonzeros(MN));
    
    figure;
    imagesc(MN)
    caxis([-1 1])
    colorbar
    colormap cool
    xticks(250:250:1000)
    xticklabels({'3','6','9','12'})
    xlabel('Time (min)')
    yticks(250:250:1000)
    yticklabels({'3','6','9','12'})
    ylabel('Time (min)')
    
    set(gca,'FontSize',20,'FontWeight','bold')
    title(['Subject MN ',num2str(s) ])
   
    Image = getframe(gcf);
    imwrite(Image.cdata, sprintf('MNstate_%d.jpg',s));
    saveas(gcf,sprintf('MNstate_%d.fig',s))
    %load(sprintf('MNstate_%d.jpg',s))
    
    F = getframe(gcf) ;
    open(v)
    writeVideo(v,F)
    close(gcf)
    % Mantel test
    [r,p]=bramila_mantel(MN,SN,100,'spearman');
    CorrMNSN(s,:)= [r,p];%MN state vector versus SN vector (spearman's r and p-value)
    [r,p]= bramila_mantel(MN,Raw,100,'spearman');
    CorrMNRaw(s,:)=[r,p]; %MN state vector versus Raw vector
    [r,p] = bramila_mantel(SN,Raw,100,'spearman');
    CorrSNraw(s,:)=[r,p]; %Raw state vector versus SN vector
    
end
save States StatesAllMN MeanStateMN MeanStateChangeMN MaxminMN StatesAllSN MeanStateSN MeanStateChangeSN MaxminSN StatesAllRaw MeanStateRaw MeanStateChangeRaw MaxminRaw 



%Calculate average constrast witin recurrence plots

for s=1:S
    e1=abs(squeeze(StatesAllRaw(s,:,:))); 
    e2=abs(squeeze(StatesAllSN(s,:,:)));
    e3=abs(squeeze(StatesAllMN(s,:,:)));
    
    Contrast(s,1)=mean(e1(:));
    Contrast(s,2)=mean(e2(:));
    Contrast(s,3)=mean(e3(:));
end

save ContrastStates Contrast


%%%%%%%%


%%Duration of SNs
load(sprintf('SN_Time')) %SN-timeseries
load(sprintf('Subnetworks_overlap%d',OV))
[SN S Tmax]=size(SNcomTS)
 
TempOr=zeros(S,SN,Tmax+1);

  %%Create timeing matrix of all SNs for all subject that track before during and after integration of SNs   
    
        CT=SNcomTS;
        CT(CT>1)=1; %Binarize timeseris
        for w=1:SN
            for s=1:S
                
                d=squeeze(CT(w,s,:));
                
                temp=d;  
                temp=[0 temp'];% add artifical 0 at the end to make catch tail
                tt=length(temp);
                %
                for t=1:length(temp)
                    if temp(t)==1 && temp(t-1)==0
                        TempOr(s,w,t)=1;
                    elseif  temp(t)==1 && temp(t-1)==-1 % added in this version to accomodate negative
                        TempOr(s,w,t)=1;
                    elseif temp(t)==1 && temp(t-1)==1
                        TempOr(s,w,t)=TempOr(s,w,t-1)+1;
                        
                    elseif t< tt-2 && temp(t)==0 && temp(t+1)==1
                        TempOr(s,w,t)=-1;
                    elseif t< tt-2 && temp(t)==0 && temp(t+1)==0 &&temp(t+2)==1
                        TempOr(s,w,t)=-2;
                    elseif t< tt-3 && temp(t)==0 && temp(t+1)==0 &&temp(t+2)==0 &&temp(t+3)==1
                        TempOr(s,w,t)=-3;
                    elseif t< tt-4 && temp(t)==0 && temp(t+1)==0 &&temp(t+2)==0 &&temp(t+3)==0 &&temp(t+4)==1
                        TempOr(s,w,t)=-4;
                    elseif t< tt-5 && temp(t)==0 && temp(t+1)==0 &&temp(t+2)==0 &&temp(t+3)==0 &&temp(t+4)==0 &&temp(t+5)==1
                        TempOr(s,w,t)=-5;
                    elseif t< tt-6 && temp(t)==0 && temp(t+1)==0 &&temp(t+2)==0 &&temp(t+3)==0 &&temp(t+4)==0 &&temp(t+5)==0 &&temp(t+6)==1
                        TempOr(s,w,t)=-6;
                    elseif t< tt-7 && temp(t)==0 && temp(t+1)==0 &&temp(t+2)==0 &&temp(t+3)==0 &&temp(t+4)==0 &&temp(t+5)==0 &&temp(t+6)==0  &&temp(t+7)==1
                        TempOr(s,w,t)=-7;
                    elseif t< tt-8 && temp(t)==0 && temp(t+1)==0 &&temp(t+2)==0 &&temp(t+3)==0 &&temp(t+4)==0 &&temp(t+5)==0 &&temp(t+6)==0  &&temp(t+7)==0 &&temp(t+8)==1
                        TempOr(s,w,t)=-8;
                    elseif t< tt-9 && temp(t)==0 && temp(t+1)==0 &&temp(t+2)==0 &&temp(t+3)==0 &&temp(t+4)==0 &&temp(t+5)==0 &&temp(t+6)==0  &&temp(t+7)==0 &&temp(t+8)==0&&temp(t+9)==1
                        TempOr(s,w,t)=-9;
                    elseif t< tt-10 && temp(t)==0 && temp(t+1)==0 &&temp(t+2)==0 &&temp(t+3)==0 &&temp(t+4)==0 &&temp(t+5)==0 &&temp(t+6)==0  &&temp(t+7)==0 &&temp(t+8)==0&&temp(t+9)==0 &&temp(t+10)==1
                        TempOr(s,w,t)=-10;
                    elseif t< tt-11 && temp(t)==0 && temp(t+1)==0 &&temp(t+2)==0 &&temp(t+3)==0 &&temp(t+4)==0 &&temp(t+5)==0 &&temp(t+6)==0  &&temp(t+7)==0 &&temp(t+8)==0&&temp(t+9)==0 &&temp(t+10)==0 && temp(t+11)==1
                        TempOr(s,w,t)=-11;
                    elseif t< tt-12 && temp(t)==0 && temp(t+1)==0 &&temp(t+2)==0 &&temp(t+3)==0 &&temp(t+4)==0 &&temp(t+5)==0 &&temp(t+6)==0  &&temp(t+7)==0 &&temp(t+8)==0&&temp(t+9)==0 &&temp(t+10)==0 && temp(t+11)==0 && temp(t+12)==1
                        TempOr(s,w,t)=-12;
                    elseif t< tt-13 && temp(t)==0 && temp(t+1)==0 &&temp(t+2)==0 &&temp(t+3)==0 &&temp(t+4)==0 &&temp(t+5)==0 &&temp(t+6)==0  &&temp(t+7)==0 &&temp(t+8)==0&&temp(t+9)==0 &&temp(t+10)==0 && temp(t+11)==0 && temp(t+12)==0 && temp(t+13)==1
                        TempOr(s,w,t)=-13;
                    elseif t< tt-14 && temp(t)==0 && temp(t+1)==0 &&temp(t+2)==0 &&temp(t+3)==0 &&temp(t+4)==0 &&temp(t+5)==0 &&temp(t+6)==0  &&temp(t+7)==0 &&temp(t+8)==0&&temp(t+9)==0 &&temp(t+10)==0 && temp(t+11)==0 && temp(t+12)==0 && temp(t+13)==0  && temp(t+14)==1
                        TempOr(s,w,t)=-14;
                    elseif t< tt-15 && temp(t)==0 && temp(t+1)==0 &&temp(t+2)==0 &&temp(t+3)==0 &&temp(t+4)==0 &&temp(t+5)==0 &&temp(t+6)==0  &&temp(t+7)==0 &&temp(t+8)==0&&temp(t+9)==0 &&temp(t+10)==0 && temp(t+11)==0 && temp(t+12)==0 && temp(t+13)==0  && temp(t+14)==0 && temp(t+15)==1
                        TempOr(s,w,t)=-15;
                    elseif t< tt-16 && temp(t)==0 && temp(t+1)==0 &&temp(t+2)==0 &&temp(t+3)==0 &&temp(t+4)==0 &&temp(t+5)==0 &&temp(t+6)==0  &&temp(t+7)==0 &&temp(t+8)==0&&temp(t+9)==0 &&temp(t+10)==0 && temp(t+11)==0 && temp(t+12)==0 && temp(t+13)==0  && temp(t+14)==0 && temp(t+15)==0 && temp(t+16)==1
                        TempOr(s,w,t)=-16;
                    elseif t< tt-17 && temp(t)==0 && temp(t+1)==0 &&temp(t+2)==0 &&temp(t+3)==0 &&temp(t+4)==0 &&temp(t+5)==0 &&temp(t+6)==0  &&temp(t+7)==0 &&temp(t+8)==0&&temp(t+9)==0 &&temp(t+10)==0 && temp(t+11)==0 && temp(t+12)==0 && temp(t+13)==0  && temp(t+14)==0 && temp(t+15)==0 && temp(t+16)==0 && temp(t+17)==1
                        TempOr(s,w,t)=-17;
                    elseif t< tt-18 && temp(t)==0 && temp(t+1)==0 &&temp(t+2)==0 &&temp(t+3)==0 &&temp(t+4)==0 &&temp(t+5)==0 &&temp(t+6)==0  &&temp(t+7)==0 &&temp(t+8)==0&&temp(t+9)==0 &&temp(t+10)==0 && temp(t+11)==0 && temp(t+12)==0 && temp(t+13)==0  && temp(t+14)==0 && temp(t+15)==0 && temp(t+16)==0 && temp(t+17)==0 && temp(t+18)==1
                        TempOr(s,w,t)=-18;
                    elseif t>1 && t<(Tmax-16) &&temp(t)==0 &&temp(t+6)==0 && TempOr(s,w,t-1)>0 && TempOr(s,w,t-1)<100%THis line adds after cluster assembly has ended (280 before when -12 was max)
                        TempOr(s,w,t)=TempOr(s,w,t-1)+101;
                    elseif t>1 && t<(Tmax-16) &&temp(t)==0 &&temp(t+6)==0 && TempOr(s,w,t-1)>0 && TempOr(s,w,t-1)>100%THis line adds after cluster assembly has ended
                        TempOr(s,w,t)=TempOr(s,w,t-1)+1;
                    elseif t>1 && t<(Tmax-17) &&temp(t)==0 &&temp(t+6)~=1 && temp(t+6)~=1 && TempOr(s,w,t-1)>0 && TempOr(s,w,t-1)>100%THis line adds after cluster assembly has ended
                        TempOr(s,w,t)=TempOr(s,w,t-1)+1;
                    elseif t>1 && t<(Tmax-18) &&temp(t)==0 &&temp(t+6)~=1 && temp(t+6)~=1 && TempOr(s,w,t-1)>0 && TempOr(s,w,t-1)>100%THis line adds after cluster assembly has ended
                        TempOr(s,w,t)=TempOr(s,w,t-1)+1;
                    else
                        TempOr(s,w,t)=0;
                        
                    end
                    
                end
            end
            
        end
 
 TempOr=TempOr(:,:,2:Tmax+1); % Remove artifical zero
   
 save TempOR_Timing TempOr
 % find avarage duration of integration, 
 % zero assemblies are not calculated
 load SN_Time QintSN
 
 TempOr(TempOr<0)=0;
 TempOr(TempOr>99)=0;
 OnlyOneDur=0;
 
     for w=1:SN
         for s=1:S
             
             D=squeeze(TempOr(s,w,:));
             b=unique(nonzeros(D)); %nr of different uniqe assembly durations
             c=histc(nonzeros(D),b); %instances per duration
             if numel(b)>1
                %for every unique length calculate the instance minus all instances with a duration 
                 % one shorter 
                 for i=1:length(b)-1
                     ry(s,w,i)=c(i)-c(i+1);
                 end
                 ry(s,w,i+1)=c(i+1); %Vector with nr of assemblies with a particular length
                 r=find(ry(s,w,:)>0);
                 aa=squeeze(ry(s,w,1:max(r)));
                
                 rr=[1:1:max(r)]'; %this will alwaysbe same as b
                 Wid(w,s)=sum(aa.*rr)/sum(aa); %mean duration (duration multiplied by nr of instances)
                 MaxWid(w,s)=max(r);
                 
                 %Calculate median
                 A=[];
                 for k=1:length(b)
                     if length(aa)>0
                         A=[A; ones(aa(k),1)*k];
                     end
                 end
                 WidMedian(w,s)=median(A);
                 
             elseif numel(b)==1 %ie if only one duration
                 Wid(w,s)=b;
                 MaxWid(w,s)=b;
                 WidMedian(w,s)=b;
                 OnlyOneDur=OnlyOneDur+1
             else
                  end
             aa=[];
         end
         
         WidInt(w,1)=mean(Wid(w,:));
         WidInt(w,2)=iqr(Wid(w,:));
         WidInt(w,3)=max(Wid(w,:));
         WidInt(w,4)=min(Wid(w,:));
         WidInt(w,5)=median(WidMedian(w,:));
         WidInt(w,6)=iqr(WidMedian(w,:));
         WidInt(w,7)=mean(WidMedian(w,:));
         WidInt(w,8)=std(WidMedian(w,:));
         WidInt(w,9)=median(Wid(w,:));
         WidInt(w,10)=iqr(Wid(w,:));
     end
save DurationSN WidInt TempOr Wid WidMedian 

load(sprintf('Subnetworks_overlap%d',OV))


figure;
a=WidInt(:,1)*TR
hold on; scatter(a,QintSN(:,1),100,'filled' )
title('Qint as a function of duration')
xlabel('Median duration (s)')
ylabel('Qint')
legend('SNs')
saveas(gcf,'SN_QintVsDuration.jpg')
saveas(gcf,'SN_QintVsDuration.fig')

%Duration as a function of total nr of areas?
figure; 

a=WidInt(:,1)*TR;
hold on; scatter(Areas_SNs,a,100,'filled' )
title('Duration as a function of total nr of areas?')
saveas(gcf,'SN_NrAreasVsDuration.jpg')
saveas(gcf,'SN_NrAreasDuration.fig')

%Duration as a function of nr of SNCs?
for i=1:SN
    temp=squeeze(D(i,:,1));
    SNCnr(i)=length(nonzeros(temp))
end
figure; 

a=WidInt(:,1)*TR;
hold on; scatter(SNCnr,a,100,'filled' )
title('Duration as a function of total nr of SNCs?')
[roh pval]=corr(SNCnr',a,'Type','Spearman')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Phase Diff Within and Between SNs

load(sprintf('Subnetworks_overlap%d',OV))
for y=1:SN
    
    load(sprintf('SNC8_TS_%d',y))
    
    parfor s=1:S
        
        if length(Qint_SNC(:,1))>1 % Only calculate difference if nr of SNCs >1
            
            e=squeeze(SNCmeanPhase(s,:,:));
            for t=1:Tmax
                k=nonzeros(e(:,t)); % Mean Leida values of SNCs belongin to SN y at time t in subj s
                M=zeros(numel(k),numel(k));
                for f=1:numel(k)
                    for g=1:numel(k)
                        if g>f
                            M(f,g)=cos(k(f)-k(g));
                        end
                    end
                end
              
                MeanDiffPhaseSN(y,s,t)=mean(nonzeros(M));
            end
        end
    end

end
MeanDiffPhaseSN(isnan(MeanDiffPhaseSN))=0;
save('DiffPhaseSNmean','MeanDiffPhaseSN', '-v7.3')
 
for s=1:S
    for t=1:Tmax
        E(s,t)=mean(nonzeros(MeanDiffPhaseSN(:,s,t)));
    end
E(isnan(E))=0;
E1(s)=mean(E(s,:));
end
mean(E1)
std(E1)


% Difference in mean phase between SNs in different communities

load(sprintf('SN_Time')) %SN-timeseries
DiffiCom=zeros(S,Tmax);
for s=1:S
    Tem=squeeze(SNmeanPhase(:,s,:)); %Mean SN phase weigt
    Tem2=squeeze(SNcomTS(:,s,:)); % SN community
    for t=1:Tmax
        u=unique(nonzeros(Tem2(:,t))); %number of commmunities
        for i=1:numel(u)
            c1=(Tem2==u(i));
            c2=logical(c1);
            
            Ac(s,t,i)= mean(Tem(c1(:,t),t)); %mean of each com
            
        end
        Nd=zeros(numel(u),numel(u));
        for i=1:numel(u)
            for k=1:numel(u)
                if k>i
                    Nd(i,k)=cos(Ac(s,t,i)-Ac(s,t,k));
                end
            end
        end
        DiffiCom(s,t)=mean(nonzeros(Nd));
    end
end

DiffiCom(isnan(DiffiCom))=0;
figure; histogram(nonzeros(DiffiCom),'Normalization','probability')
for s=1:S
    F(s)=mean(nonzeros(DiffiCom(s,:)));
end
mean(F)
std(F)
save DiffComSNsmean DiffiCom


% Difference in phase in same communitity for
%1) SNs of same MN, 2) SNs of different MN
load SN_Time
load(sprintf('HierClustV_%d',C))
H=HierOrd; % SNs per MNs sorted in MN order
H2=sortrows(HierOrd,1,'ascend'); % Sort H in SN order
SameMN=zeros(3,C,S,Tmax); %nr of communities, MN, subj,Timepoints
DiffMN=zeros(3,C,S,Tmax);
for s=1:S
    Tem=squeeze(SNmeanPhase(:,s,:)); %Median SN Phase weigt
    Tem2=squeeze(SNcomTS(:,s,:)); % SN community
    
    for t=1:Tmax
        
        u=unique(nonzeros(Tem2(:,t))); %Different communities with SNs at time t
        for q=1:numel(u) %for each community
            d1=(Tem2(:,t)==u(q)); % rows of SNs in the community
            
            A= Tem(d1,t); %Phase angle
            
            B= H2(d1,2); % MNs
            
            k=unique(B); %MNs in the same community
            
            for i=1:numel(k)% for each of the MNs
                
                r=find(B==k(i));
                rn=numel(r);
                v=ones(numel(B),1);
                v(r)=0; %rows of SNs of sam MN
                vn=B(logical(v)); % rows of SNs of other MNs
                CK=A(logical(v));
                %calculate diff within MN ie between SNs of same MN
                if rn>1
                    M=zeros(rn,rn);
                    
                    for l=1:rn
                        for  x=1:rn
                            if x>l
                                M(x,l)=cos(A(r(l))-A(r(x)));
                            end
                        end
                    end
                    
                    SameMN(q,k(i),s,t)=mean(nonzeros(M));
                else
                    SameMN(q,k(i),s,t)=0; %if MN consists of only one SN
                end
                %diff between MNs
                if numel(vn)>0
                    NL=zeros(rn,sum(v));
                    for l=1:rn
                        for x=1:sum(v)
                            if x>l
                                NL(x,l)=cos(A(r(l))-CK(x) );
                            end
                        end
                    end
                    
                    DiffMN(q,k(i),s,t)=mean(nonzeros(NL));
                else
                    DiffMN(q,k(i),s,t)=0;
                end
            end
        end
    end
end

SameMN(isnan(SameMN))=0;
DiffMN(isnan(DiffMN))=0;

for s=1:S
    for t=1:Tmax
        for q=1:3
            A1(s,t,q)=mean(nonzeros(SameMN(q,:,s,t)));
            B1(s,t,q)=mean(nonzeros(DiffMN(q,:,s,t)));
        end
        A1(isnan(A1))=0;
        A2(s,t)=mean(nonzeros(A1(s,t,:)));
        B1(isnan(B1))=0;
        B2(s,t)=mean(nonzeros(B1(s,t,:)));
    end
    A2(isnan(A2))=0;
    A3(s)=mean(nonzeros(A2(s,:)));
    B2(isnan(B2))=0;
    B3(s)=mean(nonzeros(B2(s,:)));
    
end
A3(isnan(A3))=0;
B3(isnan(B3))=0;
SumMN(1,1)=mean(A3)
SumMN(1,2)=std(A3)
SumMN(1,3)=mean(B3)
SumMN(1,4)=std(B3)
SumMN=round(SumMN,4)


save(sprintf('MNdiff_K%d',K), 'DiffMN', 'SameMN', 'SumMN')

figure; histogram(nonzeros(A3),'Normalization','probability')
hold on; histogram(nonzeros(B3),'Normalization','probability')

[p,h,stats] = ranksum(nonzeros(SameMN),nonzeros(DiffMN))


 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNC integration and random component integration
% 1)first get L random SNCs and random components with timeseries
% 2)calculate phase integration prior to and during integration into the
% same community

load Top

L=input('How many SNC and random components for the phase integration calculation')
b=[0 1]
for i=1:2
    a=b(i)
    PhaseIntegrationSNC(a,Tmax,Louv,TR,S,N,K,L,SNC_size)
end

%Permutation t-test to test if there is any difference in max phase
%synchrony between SNCs and random components
load SNCLines
c1=MaxMean;
load RandLines
c2=MaxMean;
%onesided
[p2, observeddifference, effectsize2] = permutationTest(c2,c1, 1000,'sidedness','smaller','plotresult',1 )  
%twosided
[p2, observeddifference, effectsize2] = permutationTest(c2,c1, 1000,'plotresult',1 )  
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Plot an example five first SNs (idenpendent of Qint order)
example=squeeze(mean_amp(:,1,:));
figure; hold on
t=0
tt=[5 4 3 2 1]
for i=1:5
   plot(example(tt(i),:)+t,'LineWidth',1) 
   t=t+1;  
end
%if Schaefer 200 + BrainNetome 236 subcortical
% xticks(250:250:1000)
% xticklabels({'3','6','9','12'})
% xlabel('Time (min)')
% yticklabels({'','SN5','SN4','SN3','SN2','SN1'})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Calculate number of SNs per timepoint  
%calculate nr of Communities SNs are integrated in
%calculate nr of areas per SN
%calculate nr of areas per tp
load SN_Time
for s=1:S
    d=squeeze(SNcomTS(:,s,:));
    d1=d; d1(d1>0)=1;
    Kn(s,:)=sum(d1);
    Kmin(s)=numel(nonzeros(sum(d1)==0));
    % find nr of unique communities:
    c1 =(d==1);
    c11=sum(c1,1);
    c11(c11>0)=1;
    c2 =(d==2);
    c22=sum(c2,1);
    c22(c22>0)=1;
    c3=(d==3);
    c33=sum(c3,1);
    c33(c33>0)=1;
    KD(s,:)=c11+c22+c33;
end

% sum(Kmin)/(S*Tmax)
% mean(Kmin) 
% std(Kmin)
MeanCom=mean(KD,2);
MeanMeanCom(1)=mean(MeanCom)
MeanMeanCom(2)=std(MeanCom)

MeanSN=mean(Kn,2);
MeanMeanSN(1)=mean(MeanSN)
MeanMeanSN(2)=std(MeanSN)

%timepoint with no SN5/(Tmax*S)
tt=squeeze(sum(SNcomTS));
sum(sum(tt'==0))


%%%%Calculate nr of unique areas at each timepoint, percentage of total areas
%at each timepoint
load('Time_normalized_IMF_C.mat')
load(sprintf('Subnetworks_overlap%d',OV))
[a1 a2]=size(SNs)

for s=1:S
  
   a=squeeze(TimeArea(:,s,:,:)); 
 
   for t=1:Tmax
      Uni(s,t)=numel(unique(nonzeros(a(:,t,:))))/N; 
   end
   meanArea(s,1)=mean(Uni(s,:));
   meanArea(s,2)=std(Uni(s,:));
   
end
mean(meanArea(:,1))
std(meanArea(:,1))
min((meanArea(:,1)))
max(meanArea(:,1))

%%%%Number of Louvain communities independent of if SNs are integrated in them or
%not

L=Louv;
for s=1:S
    d=squeeze(L(s,:,:))';
    c1 =(d==1);
    c11=sum(c1,1);
    c11(c11>0)=1;
    c2 =(d==2);
    c22=sum(c2,1);
    c22(c22>0)=1;
    c3=(d==3);
    c33=sum(c3,1);
    c33(c33>0)=1;
    ComD(s,:)=c11+c22+c33;
    meanComD(s)=mean(ComD(s,:));
end
mean(meanComD)
std(meanComD)
numel(nonzeros(ComD==1))/(S*Tmax)

%%
%%%%%Flexibility/modularity%%%%%%%%%%%

% Calculate number of SNCs  per area 
% Percentage of total areas that each area share an SNC with - Pdiv
load(sprintf('P%d_I%d',K,SNC_size))
if SNC_size==8||SNC_size==7
    load(sprintf('P%d_I8',K))
    H=R8_IDXsU;
    H2=R8_IDXsUfreqSort;
elseif SNC_size==6||SNC_size==5
    load(sprintf('P%d_I6',K))
    H=R6_IDXsU;
    H2=R6_IDXsUfreqSort;
elseif SNC_size==4||SNC_size==3
    load(sprintf('P%d_I4',K))
    H=R4_IDXsU;
    H2=R4_IDXsUfreqSort;
elseif SNC_size==10||SNC_size==9
    load(sprintf('P%d_I10',K))
    H=R10_IDXsU;
    H2=R10_IDXsUfreqSort;
elseif SNC_size==12||SNC_size==11
    load(sprintf('P%d_I12',K))
    H=R12_IDXsU;
    H2=R12_IDXsUfreqSort;
end
tot=length(H);
for i=1:N
    NumSNCarea(i,1)=sum(sum(H==i)); %nr of SNCs with that area
    NumSNCarea(i,2)=sum(sum(H==i))/tot; %percent of total SNCs with that area
    %weight by qint of the SNC
    SumQintSNCarea(i)=sum(H(logical(sum(H==i,2)),1));
    %total nr of other areas in SNCs
    Pdiv(i)=numel(unique(H(logical(sum(H==i,2)),:)))-1;
end

Pdiv=Pdiv/N;

%%%%Average freq of areas as part of any SNC: Pint
load Time_normalized_IMF_C
[SN S Tmax Ar]=size(TimeArea)

AreaTime=zeros(S,Tmax,Ar);
AreaTimeCount=zeros(S,N);
for s=1:S
    Time=squeeze(TimeArea(:,s,:,:));
    
    TA=permute(Time,[2 1 3 ]);
   
    AA=zeros(1,N);
    
    for t=1:Tmax
        unA=nonzeros(unique(TA(t,:,:)));
        AreaTime(s,t,1:numel(unA))=unA; %Unique areas in s at time t
        AA(unA)=AA(unA)+1;
    end
    AreaTimeCount(s,:)=AA; % Number of Timepoints areas appear
    
end
if S>1
    Pint=mean(AreaTimeCount/Tmax); %Percentage of area is part of any SN
else
    Pint=AreaTimeCount/Tmax;
end
save(sprintf('AreaTIMEconK%d',K), 'AreaTimeCount','AreaTime', 'Pint', 'Pdiv')



load(sprintf('AreaTIMEconK%d',K))
load AreaParcelA236_N9 AreaNetw

figure; hold on
P=9 % Nr of static networks in parcellation
for i=1:P
    r=find(AreaNetw(:,1)==i);
    CO(i)=corr2(Pint(r),Pdiv(r));
    scatter(Pint(r),Pdiv(r),120,'filled')
end
ax = gca;
ax.XAxis.FontSize = 24;
ax.YAxis.FontSize = 24;
ax.FontWeight = 'bold';
ax.Title.FontSize = 24;
thr1=mean(Pint);
thr2=mean(Pdiv);
xlabel('Pint')
ylabel('Pdiv')
hold on;
line([thr1 thr1],[0 1],'LineWidth',2,'Color','k');
line([0 1],[thr2 thr2], 'LineWidth',2);
legend ({'VIS','SOM','DAN','VAN','Limbic','FPN','DMN','BG','Thal','M Pint','M Pdiv'},'Location','northwest','NumColumns',1,'Fontsize',20)



% Brain maps of  Pdiv and Pint

% You might need to change the range
setenv('FSLOUTPUTTYPE','NIFTI_GZ')
getenv('FSLOUTPUTTYPE')

%Set map scale
maxPdiv=round(max(100*Pdiv));
minPdiv=round(min(100*Pdiv));
load Pdiv
EC.vol.px=maxPdiv;
EC.vol.pn=minPdiv;
save Pdiv EC

AP=(Pdiv)*100; %the weights
A=[1:N];
a1=A(1);
cs = sprintf('cp /Users/marist1/Documents/MATLAB/Parcellation236/%dSH_7mask.nii.gz Pdiv.nii.gz',a1);
system(cs)
%weigthed:
cmd=sprintf('/usr/local/fsl/bin/fslmaths  Pdiv.nii.gz -mul %d  Pdiv.nii.gz ',AP(1));
system(cmd)

for a=2:numel(A)
    n=A(a);
    %with weight
    cmd = sprintf('/usr/local/fsl/bin/fslmaths  /Users/marist1/Documents/MATLAB/Parcellation236/%dSH_7mask.nii.gz -mul %d temp_%d.nii.gz',n,AP(a),a);
    system(cmd)
    k(a)=AP(a);
    cmd = sprintf('/usr/local/fsl/bin/fslmaths Pdiv.nii.gz -add temp_%d.nii.gz Pdiv.nii.gz',a);
    system(cmd);
    delete(sprintf('temp_%d.nii.gz',a))
end

BrainNet_MapCfg('BrainMesh_ICBM152.nv','Pdiv.nii.gz','Pdiv.mat',sprintf('Pdiv%d.jpg',K));


% Pint map
%Set map scale
maxPint=round(max(100*Pint));
minPint=round(min(100*Pint));
load Pint
EC.vol.px=maxPint;
EC.vol.pn=minPint;
save Pint EC

AP=Pint*100; %the weights
A=[1:N];
a1=A(1);
cs = sprintf('cp /Users/marist1/Documents/MATLAB/Parcellation236/%dSH_7mask.nii.gz Pint.nii.gz',a1);
system(cs)
%weigthed:
cmd=sprintf('/usr/local/fsl/bin/fslmaths  Pint.nii.gz -mul %d  Pint.nii.gz ',AP(1));
system(cmd)

for a=2:numel(A)
    n=A(a);
    %with weight
    cmd = sprintf('/usr/local/fsl/bin/fslmaths  /Users/marist1/Documents/MATLAB/Parcellation236/%dSH_7mask.nii.gz -mul %d temp_%d.nii.gz',n,AP(a),a);
    system(cmd)
    k(a)=AP(a);
    cmd = sprintf('/usr/local/fsl/bin/fslmaths Pint.nii.gz -add temp_%d.nii.gz Pint.nii.gz',a);
    system(cmd);
    delete(sprintf('temp_%d.nii.gz',a))
end

BrainNet_MapCfg('BrainMesh_ICBM152.nv','Pint.nii.gz','Pint.mat',sprintf('Pint%d.jpg',K));



end
