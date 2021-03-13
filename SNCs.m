function []=SNCs(K,Louv,V,Q,Tmax,S,N,SNC_size,thr,evenodd)
% Marika Strindberg, Karolinska Institutet 2018-2020

% marika.strindberg@ki.se, marikastrindberg@gmail.com
% Subrutine to TimeResolvedNets that calculates SNCs of size n = 4, n = 6, n = 8, n = 10,
% n = 12 areas as specified by SNC_size (alternatively size n = 5, n = 7, n = 9, n = 11 if single seeds are specified)

% At each iteration the mean accross participants are calculated and the
% top K constellations for each seed is used as seed in the next
% iteration. 
% The runtime is sensitive to number of subjects, number of areas in the parcellation, number of seed pairs/areas but much less to total number of time-points.

% Q = contains the initial seed pairs (or seed if only on area is chosen as seed)
% K = the expansion factor ie how many top pairs in terms of qint is kept
% for each seed(SNC) in each iteration. This number is allowed to be as
% great as K+add to accommodate multiple pairs with same qint
% V = the qint-matrix for all pairs
% Louv = contains the community assignment of all areas at all timepoints. The
% numbers range from 1-3 and does not have any meaning in them selves beyond
%markers of different groups.
% Tmax =  the number of timepoints
% thr = by default 0, if set higher it only searches trough SNCs with qint
% above this thr
% S = number of subjects
% N = number of areas in parcellation
% SNC_size = the size of the final SNC ie number of areas they contain,
% default = 8. If one area is used instead of a pair it will be an odd
% number.


% Max number of components that have the same qint as the K:th highest
add=40;
%Divide computations into smaller units
Div=100;
%a=0;
if length(Q)>Div
    
    X=floor(length(Q)/Div)
    B=[Div:Div:X*Div];
    B=[0 B];
    
    A=ones(1,X)*Div; % If A is even multiplication of div
    if (length(Q)-X*Div)>0 % If A is not an even multiplication of div
        A=[A length(Q)-X*Div];
    else
    end
    Y=numel(A)
else
    
    X=floor(length(Q)/Div)
    B=[Div:Div:X*Div];
    B=[0 B];
    
    Div=length(Q)
    A=length(Q)
    Y=numel(A)
end

%Identify pairs with a qint above threshold a
[row col]=find(V>thr);
R=[row col];

R=sort(R,2);
F=unique(R,'rows');


if evenodd==0
    for p=1:Y
        meanNet_2_4=zeros(A(p),S,N,N);
        parfor h=1:A(p)
            
            mN_2_4=zeros(S,N,N);
            
            for s=1:S
                
                temp2=squeeze(Louv(s,:,:));% Time and community
                for m=1:length(F) % only pairs above thr a is seached through
                    
                    if ~ismember(Q(B(p)+h,:),F(m,:))
                        
                        InCom=0;
                        assembly=[Q(B(p)+h,:) F(m,:)];
                        
                        y=temp2(:,assembly);
                        z=prod(y,2); %Row with same community nr have product 1^4=1,2^4=16,3^4=81
                        zz=(z==[1 16 81]);
                        zzz=(zz(:,1)+zz(:,2)+zz(:,3));
                        InCom=numel(nonzeros(zzz));
                        mN_2_4(s,F(m,1),F(m,2))=InCom/Tmax;
                        
                    end
                end
                
                meanNet_2_4(h,s,:,:)=mN_2_4(s,:,:);
                
            end
            
        end
        MeanALLNetw2_4=zeros(A(p),N,N);
        if S>1
            for h=1:A(p) %calculate mean for each matrix in each row
                G=squeeze(meanNet_2_4(h,:,:,:));
                MeanALLNetw2_4(h,:,:)=mean(G);
            end
        else
            
            G=squeeze(meanNet_2_4(:,1,:,:));
            MeanALLNetw2_4(:,:,:)=G;
        end
        
        
        save(sprintf('P%d_MeanAllSNC4_%d',K,p),'MeanALLNetw2_4','-v7.3')
        
    end
else
    
    for p=1:Y
        meanNet_2_4=zeros(A(p),S,N,N);
        parfor h=1:A(p)
            
            mN_2_4=zeros(S,N,N);
            
            for s=1:S
                
                temp2=squeeze(Louv(s,:,:));% Time and community
                for m=1:length(F) % only pairs above threshold a is seached through
                    
                    if ~ismember(Q(B(p)+h,:),F(m,:))
                        
                        InCom=0;
                        assembly=[Q(B(p)+h,:) F(m,:)];
                        
                        y=temp2(:,assembly);
                        z=prod(y,2); %Row with same community nr have product 1^3=1,2^^=8,3^3=27
                        zz=(z==[1 8 27]);
                        zzz=(zz(:,1)+zz(:,2)+zz(:,3));
                        InCom=numel(nonzeros(zzz));
                        mN_2_4(s,F(m,1),F(m,2))=InCom/Tmax;
                        
                    end
                end
                
                meanNet_2_4(h,s,:,:)=mN_2_4(s,:,:);
                
            end
            
        end
        MeanALLNetw2_4=zeros(A(p),N,N);
        if S>1
            for h=1:A(p) %calculate mean for each matrix in each row
                G=squeeze(meanNet_2_4(h,:,:,:));
                MeanALLNetw2_4(h,:,:)=mean(G);
            end
        else
            
            G=squeeze(meanNet_2_4(:,1,:,:));
            MeanALLNetw2_4(:,:,:)=G;
        end
        
        
        save(sprintf('P%d_MeanAllSNC4_%d',K,p),'MeanALLNetw2_4','-v7.3')
        
    end
    
end

% For each pair calculate the the top 10 most likely 4 node components

DnumelAll=[];
for p=1:Y
    clear MeanAllNetw2_4
    load(sprintf('P%d_MeanAllSNC4_%d',K,p)) % Load search space matrix
    R4_IDXt=zeros(A(p),K+add,4);
    
    R4_length_t=zeros(1,A(p));
    
    Freq_Y=zeros(A(p),K+add); %Save all qints for division Y
    for h=1:A(p)   %For each 2 area seed
        R4_t=zeros(K+add,4); % Stores the 4 area seeds
        
        D=squeeze(MeanALLNetw2_4(h,:,:));
        Dnumel(h)=numel(nonzeros(D));
        if Dnumel(h)>0 && Dnumel(h)>=K %if more than K pairs
            a=maxk(nonzeros(D),K);
            [row col]=find(D>=a(K));%pairs above or equal to 10th highest value
            R4_t(1:numel(row),:)=[Q(B(p)+h,:).*ones(numel(row),2) row col]; %Save the top 4 area seeds
            Freq_Y(h,1:numel(row))=D(D>=a(K)); % The qint of the pairs
            
        elseif Dnumel(h)>0
            [row col]=find(D>0);
            R4_t(1:numel(row),:)=[Q(B(p)+h,:).*ones(numel(row),2) row col];%Save the top 4 area seeds
            Freq_Y(h,1:numel(row))=D(D>0);
            
        else
            R4_t(1,:)=[Q(B(p)+h,:) 0 0];
            Freq_Y(h,1)=0;
            
        end
        R4_IDXt(h,:,:)=R4_t(:,:);
         
    end
    
    save(sprintf('P%d_R4prob%d',K,p),'Dnumel','Freq_Y', 'R4_IDXt','-v7.3') 
    DnumelAll=[DnumelAll Dnumel];
    
end
save SizeS2  DnumelAll

%dubblecheck that the SNCs have unique areas

for p=1:Y
    load(sprintf('P%d_R4prob%d',K,p))
    h=1;
    Unika=[];
    for i=1:length(nonzeros(R4_IDXt(:,1,1)))
        temp= squeeze(R4_IDXt(i,:,:));
        for j=1:length(nonzeros(temp(:,1)))
            Unika(h)=numel(unique(temp(j,:)));
            h=h+1;
        end
    end
    MinR4IDX(p)=min(Unika)
    
end
save(sprintf('UniqueR4IDXt_%d',p), 'MinR4IDX')

%Add to one matrix
R4_IDX=[];R4_length=[];Freq=[];
for p=1:Y
    
    load(sprintf('P%d_R4prob%d',K,p))
    
    R4_IDX=[R4_IDX; R4_IDXt];
    
    Freq=[Freq; Freq_Y]; %Probability of the SNC components
    save(sprintf('P%d_R4idx',K), 'R4_IDX',  'Freq')
    clear R R4_IDXt  Freq_Y R4_length_t
    delete(sprintf('P%d_R4prob%d.mat',K,p))
end


%Sort from lowest to highest area nr in each row in order to later remove
%dubplicates
for h=1:length(R4_IDX(:,1,1))
    for k=1:length(R4_IDX(h,:,1))
        R4_IDXs(h,k,:)=sort(R4_IDX(h,k,:));
    end
end



R4_IDXsALL=[];
R4_IDXsALL_ORD=[];
R4_Freq_ALL=[];

for h=1:length(nonzeros(R4_IDXs(:,1,1)))
    
    T=squeeze(R4_IDXs(h,:,:));
    [H_i I]=unique(T,'rows','stable'); %zero row at the end in most cases
    Nu=ones(numel(I),1).*h; % saves the origin, ie which row the seed is derived from in R4_IDXsU
    R4_Freq_ALL=[R4_Freq_ALL; Freq(h,I)']; % get the accompanying qint for each of the SNCs
    R4_IDXsALL=[R4_IDXsALL;H_i];
    R4_IDXsALL_ORD=[R4_IDXsALL_ORD; Nu ];
    
end

%delete rows with zeros
[ro col]=find(R4_IDXsALL==0);
zeRO=unique(ro); numel(zeRO)
R4_IDXsALL(zeRO,:)=[];
R4_IDXsALL_ORD(zeRO,:)=[];
R4_Freq_ALL(zeRO)=[]; 


%Get unique rows
[R4_IDXsU I4]=unique(R4_IDXsALL,'rows','stable'); % keep the order, dont sort
FreqIDXsU4=R4_Freq_ALL(I4);%%%%
R4_IDXsUfreq=[FreqIDXsU4(:,1) R4_IDXsU];
R4_IDXsUfreqSort=sortrows(R4_IDXsUfreq,1,'descend');

save(sprintf('P%d_I4',K),'R4_IDXsUfreq','R4_IDXsUfreqSort', 'R4_Freq_ALL','FreqIDXsU4', 'R4_IDXsU', 'I4', 'R4_IDXsALL', 'R4_IDXsALL_ORD', 'R4_IDX', 'R4_IDXs', '-v7.3') % 'ProbDR4_ALL',


   

% Dubbelcheck no rows with duplicate area
duplicate=0;
U=0;
for hh=1:length(R4_IDXsU)
    U(hh)=numel(unique(R4_IDXsU(hh,:)));
    if evenodd==0
        if U(hh)<4
            duplicate=1
            U(hh)
            save(sprintf('P%d_4D',K), 'U')
            break
        end
    else
        if U(hh)<3
            duplicate=1
            save(sprintf('P%d_4D',K), 'U')
            break
        end
        
    end
end


fileID = fopen('Step1.txt','w');
fprintf(fileID,'Number of rows in step 2 is %d\n',size(R4_IDXsU,1));
fprintf(fileID,'Duplicate %d \n', duplicate);
fprintf(fileID,'Number of divisions %d', Y);
fclose(fileID);
Dnumel=[];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Step2 ie from SNC=4 to SNC=6
if SNC_size>4
    
    Div=100
    
    X=floor(numel(I4)/Div); % Division not to make matrix too big
    B=[Div:Div:X*Div];
    B=[0 B];
    
    A=ones(1,X)*Div;
    if (length(I4)-X*Div)>0 % If A is not an even multiplication of div
        A=[A numel(I4)-X*Div]; % Nr of rows  in each division
    else
    end
    Y=numel(A) %Nr of divisions
    
    if evenodd==0
        for p=1:Y
            
            clear  TempMeanALL E MeanALLNetw2_4_6
            meanNet_2_4_6=zeros(A(p),S,N,N);
            if p< Y
                R4ORDtemp=R4_IDXsALL_ORD(I4(B(p)+1:B(p+1))); % length decided by A(p)
                
            else
                R4ORDtemp=R4_IDXsALL_ORD(I4(B(p)+1:B(p)+A(p)));
                
            end
            
            
            TempMeanALL=zeros(length(R4ORDtemp),N,N);
            for i=1:length(R4ORDtemp) % length decided by A(p) ie max=div
                MeanALLNetw2_4=[];
                F=R4ORDtemp(i);
                b1=floor(F/Div);
                b2=mod(F,Div);
                if b2==0
                    load(sprintf('P%d_MeanAllSNC4_%d',K,b1));
                    TempMeanALL(i,:,:)=MeanALLNetw2_4(Div,:,:);
                else
                    load(sprintf('P%d_MeanAllSNC4_%d',K,b1+1));
                    TempMeanALL(i,:,:)=MeanALLNetw2_4(b2,:,:);
                end
                
                %Identify pairs with a qint above threshold a to search through
                temp=squeeze(TempMeanALL(i,:,:));
                [row col]=find(temp>thr);
                Rx=[row col];
                Rx=sort(Rx,2);
                Rx=unique(Rx,'rows');
                E{i}=Rx; %E has max length=div
            end
            
            
            parfor h=1:A(p)
                mN_2_4_6=zeros(S,N,N);
                F=cell2mat(E(h));
                for s=1:S
                    temp2=squeeze(Louv(s,:,:));% Time and community
                    for m=1:length(F(:,1)) %number of pairs to search through
                        
                        if ~ismember(R4_IDXsU(h+B(p),:),F(m,:))
                            InCom=0;
                            assembly=[R4_IDXsU(h+B(p),:) F(m,:)];
                            y=temp2(:,assembly);
                            z=prod(y,2); %Row with same community nr have product 1^6=1,2^6=64,3^729
                            zz=(z==[1 64 729]);
                            zzz=(zz(:,1)+zz(:,2)+zz(:,3));
                            InCom=numel(nonzeros(zzz));
                            mN_2_4_6(s,F(m,1),F(m,2))=InCom/Tmax;
                            
                        end
                    end
                    
                    meanNet_2_4_6(h,s,:,:)=mN_2_4_6(s,:,:);
                    
                end
                
            end
            MeanALLNetw2_4_6=zeros(A(p),N,N);
            if S>1
                for h=1:A(p) %calculate mean for each matrix in each row
                    G=squeeze(meanNet_2_4_6(h,:,:,:));
                    MeanALLNetw2_4_6(h,:,:)=mean(G);
                end
            else
                G=squeeze(meanNet_2_4_6(:,1,:,:));
                MeanALLNetw2_4_6(:,:,:)=G;
            end
            
            save(sprintf('P%d_MeanAllSNC6_%d',K,p),'MeanALLNetw2_4_6','-v7.3')
        end
        
    else
        for p=1:Y
            
            clear  TempMeanALL E MeanALLNetw2_4_6
            meanNet_2_4_6=zeros(A(p),S,N,N);
            if p< Y
                R4ORDtemp=R4_IDXsALL_ORD(I4(B(p)+1:B(p+1))); % length decided by A(p)
                
            else
                R4ORDtemp=R4_IDXsALL_ORD(I4(B(p)+1:B(p)+A(p)));
                
            end
            
            
            TempMeanALL=zeros(length(R4ORDtemp),N,N);
            for i=1:length(R4ORDtemp) % length decided by A(p) ie max=div
                MeanALLNetw2_4=[];
                F=R4ORDtemp(i);
                b1=floor(F/Div);
                b2=mod(F,Div);
                if b2==0
                    load(sprintf('P%d_MeanAllSNC4_%d',K,b1));
                    TempMeanALL(i,:,:)=MeanALLNetw2_4(Div,:,:);
                else
                    load(sprintf('P%d_MeanAllSNC4_%d',K,b1+1));
                    TempMeanALL(i,:,:)=MeanALLNetw2_4(b2,:,:);
                end
                
                %Identify pairs with a qint above threshold a to search through
                temp=squeeze(TempMeanALL(i,:,:));
                [row col]=find(temp>thr);
                Rx=[row col];
                Rx=sort(Rx,2);
                Rx=unique(Rx,'rows');
                E{i}=Rx; %E has max length=div
            end
            
            
            parfor h=1:A(p)
                mN_2_4_6=zeros(S,N,N);
                F=cell2mat(E(h));
                for s=1:S
                    temp2=squeeze(Louv(s,:,:));% Time and community
                    for m=1:length(F(:,1)) %number of pairs to search through
                        
                        if ~ismember(R4_IDXsU(h+B(p),:),F(m,:))
                            InCom=0;
                            assembly=[R4_IDXsU(h+B(p),:) F(m,:)];
                            y=temp2(:,assembly);
                            z=prod(y,2); %Row with same community nr have product 1^5=1,2^5=32,3^5=243
                            zz=(z==[1 32 243]);
                            zzz=(zz(:,1)+zz(:,2)+zz(:,3));
                            InCom=numel(nonzeros(zzz));
                            mN_2_4_6(s,F(m,1),F(m,2))=InCom/Tmax;
                            
                        end
                    end
                    
                    meanNet_2_4_6(h,s,:,:)=mN_2_4_6(s,:,:);
                    
                end
                
            end
            MeanALLNetw2_4_6=zeros(A(p),N,N);
            if S>1
                for h=1:A(p) %calculate mean for each matrix in each row
                    G=squeeze(meanNet_2_4_6(h,:,:,:));
                    MeanALLNetw2_4_6(h,:,:)=mean(G);
                end
            else
                G=squeeze(meanNet_2_4_6(:,1,:,:));
                MeanALLNetw2_4_6(:,:,:)=G;
            end
            
            save(sprintf('P%d_MeanAllSNC6_%d',K,p),'MeanALLNetw2_4_6','-v7.3')
        end
     
        
    end
    
    
    DnumelAll=[];
    for p=1:Y
        load(sprintf('P%d_MeanAllSNC6_%d',K,p)) % Load search space matrix
        
        R6_IDXt=zeros(A(p),K+add,6);
        R6_freqt=zeros(A(p),K+add);
        R6_length_t=zeros(1,A(p));
        
        Freq_Y=zeros(A(p),K+add);
        for h=1:A(p)   %For each 4 area seed
            R6_t=zeros(K+add,6); %Stores the 6 area seed
            
            
            D=squeeze(MeanALLNetw2_4_6(h,:,:));
            Dnumel(h)=numel(nonzeros(D));
            if Dnumel(h)>0 && Dnumel(h)>=K %if more than 10 pairs
                a=maxk(nonzeros(D),K); %The pairs above or equal to 10th highest value
                [row col]=find(D>=a(K));
                R6_t(1:numel(row),:)=[R4_IDXsU(B(p)+h,:).*ones(numel(row),4) row col];
                Freq_Y(h,1:numel(row))=D(D>=a(K)); % The qint of the pairs
                
            elseif Dnumel(h)>0
                [row col]=find(D>0);
                R6_t(1:numel(row),:)=[R4_IDXsU(B(p)+h,:).*ones(numel(row),4) row col];
                Freq_Y(h,1:numel(row))=D(D>0);
                
            else
                R6_t(1,:)=[R4_IDXsU(B(p)+h,:) 0 0];
                Freq_Y(h,1)=0;
                
            end
            R6_IDXt(h,:,:)=R6_t(:,:);
           
        end
        
        save(sprintf('P%d_R6prob%d',K,p),'Dnumel','Freq_Y', 'R6_IDXt','-v7.3')
        DnumelAll=[DnumelAll Dnumel];
    end
    save SizeS2 DnumelAll
    
    %Add to one matrix
    R6_IDX=[];R6_Freq=[];
    for p=1:Y
        
        load(sprintf('P%d_R6prob%d',K,p))
        
        R6_IDX=[R6_IDX; R6_IDXt];
        
        R6_Freq=[R6_Freq; Freq_Y]; %Probability of the SNC component
        
        save(sprintf('P%d_R6idx',K), 'R6_IDX', 'R6_Freq')
        clear R R6_IDXt  ProbDR6t ProbIn6t
        delete(sprintf('P%d_R6prob%d.mat',K,p))
        
    end
    
    
       
        
       %Sort each SNC from lowest to highest area nr in order to later
        %delete dublicate rows that have converged in this  step
        R6_IDXs=[];
        R6_IDXsALL=[];
        R6_IDXsALL_ORD=[];
        R6_Freq_ALL=[];
        
        q=mod(length(R6_IDX),Div)
        qa=floor(length(R6_IDX)/Div) %Rows in each division
        
        
        if qa>1000
            Av=0;
            for i=1:Div
                Av(i+1,1)=qa*i;
            end
            
            Le=ones(Div,1)*qa;
            if q>0
                Av(Div+1)=q;
                Le=[Le; q];
                Di=Div+1; %Number of divisions
            else
                Di=Div;
            end
            

            for g=1:Di
                
                temp=[];
                R6_IDXsALLtemp=[];
                R6_IDXsALL_ORDtemp=[];
                R6_Freq_ALLtemp=[];
                for h=1:Le(g)
                    
                    for k=1:length(R6_IDX(h+Av(g),:,1))
                        temp(h,k,:)=sort(R6_IDX(h+Av(g),k,:));
                    end
                    H_i=[]; I=[];Nu=[];
                    T=squeeze(temp(h,:,:));
                    [H_i I]=unique(T,'rows','stable');
                    
                    Nu=ones(numel(I),1).*h;
                    
                    R6_Freq_ALLtemp=[R6_Freq_ALLtemp; R6_Freq(h+Av(g),I)'];
                    R6_IDXsALLtemp=[R6_IDXsALLtemp;H_i];
                    R6_IDXsALL_ORDtemp=[R6_IDXsALL_ORDtemp;Nu];
                    
                end
                
                R6_IDXs=[R6_IDXs; temp];
                R6_IDXsALL=[R6_IDXsALL;R6_IDXsALLtemp];
                R6_IDXsALL_ORD=[R6_IDXsALL_ORD;R6_IDXsALL_ORDtemp]
                R6_Freq_ALL=[R6_Freq_ALL;R6_Freq_ALLtemp];
                
            end 
            
        else
            for h=1:length(R6_IDX)
                for k=1:length((R6_IDX(h,:,1)))
                    R6_IDXs(h,k,:)=sort(R6_IDX(h,k,:));
                end
                H_i=[]; I=[];Nu=[];
                T=squeeze(R6_IDXs(h,:,:));
                [H_i I]=unique(T,'rows','stable');
                
                Nu=ones(numel(I),1).*h;
                
                R6_Freq_ALL=[R6_Freq_ALL; R6_Freq(h,I)'];
                R6_IDXsALL=[R6_IDXsALL;H_i];
                R6_IDXsALL_ORD=[R6_IDXsALL_ORD;Nu];
                
            end
            
        end
  
    % %delete rows with zeros
    [ro col]=find(R6_IDXsALL==0);
    zeRO=unique(ro); numel(zeRO)
    R6_IDXsALL(zeRO,:)=[];
    R6_IDXsALL_ORD(zeRO)=[];
    R6_Freq_ALL(zeRO)=[];
    
    [R6_IDXsU I6]=unique(R6_IDXsALL,'rows');
    FreqIDXsU6=R6_Freq_ALL(I6);
    R6_IDXsUfreq=[FreqIDXsU6 R6_IDXsU];
    R6_IDXsUfreqSort=sortrows(R6_IDXsUfreq,1,'descend');
    save(sprintf('P%d_I6',K), 'R6_IDXsUfreqSort', 'R6_IDXsU', 'FreqIDXsU6', 'R6_IDXs', 'I6',  'R6_IDXsALL','R6_IDXsALL_ORD', 'R6_IDXsUfreq','R6_Freq_ALL','-v7.3') 
    
    %Dubbelcheck no dublicates
    duplicate=0;
    for hh=1:length(R6_IDXsU)
        U(hh)=numel(unique(R6_IDXsU(hh,:)));
        if evenodd==0
            if U(hh)<6
                duplicate=1
                save(sprintf('P%d_6',K), 'U')
                break
            end
        else
            if U(hh)<5
                duplicate=1
                save(sprintf('P%d_6',K), 'U')
                break
            end
        end
    end
    
    
    fileID = fopen('Step2.txt','w');
    fprintf(fileID,'Number of rows in step 3 is %d\n',size(R6_IDXsU,1)); 
    fprintf(fileID,'area dubplicate in component %d \n', duplicate);
    fprintf(fileID,'min search space %d ',min(DnumelAll));
    fprintf(fileID,'Number of divisions %d', Y);
    fclose(fileID);
    
    if SNC_size>6
       
        % From SNC=6 to SNC=8 (or from 5 to 7)
        R6_IDXsU=[];
        load(sprintf('P%d_I6',K))
        
        X=floor(numel(I6)/Div)
        B=[Div:Div:X*Div];
        B=[0 B];
        
        A=ones(1,X)*Div;
        if (length(I6)-X*Div)>0 % If A is not an even multiplication of div
            A=[A numel(I6)-X*Div];
        else
        end
        Y=numel(A)
        
        if evenodd==0
            for p=1:Y
                clear  TempMeanALL E  MeanALLNetw2_4_6_8 
                meanNet2_4_6_8=zeros(A(p),S,N,N);
                if p< Y
                    R6ORDtemp=R6_IDXsALL_ORD(I6(B(p)+1:B(p+1)));
                else
                    R6ORDtemp=R6_IDXsALL_ORD(I6(B(p)+1:B(p)+A(p)));
                end
                
                TempMeanALL=zeros(length(R6ORDtemp),N,N);
                for i=1:length(R6ORDtemp)
                    MeanALLNetw2_4_6=[];
                    F=R6ORDtemp(i);
                    b1=floor(F/Div);
                    b2=mod(F,Div);
                    if b2==0
                        load(sprintf('P%d_MeanAllSNC6_%d',K,b1))
                        TempMeanALL(i,:,:)=MeanALLNetw2_4_6(Div,:,:);
                    else
                        load(sprintf('P%d_MeanAllSNC6_%d',K,b1+1))
                        TempMeanALL(i,:,:)=MeanALLNetw2_4_6(b2,:,:);
                    end
                    %Identify pairs with a qint above threshold a to search through
                    temp=squeeze(TempMeanALL(i,:,:));
                    [row col]=find(temp>thr);
                    Rx=[row col];
                    Rx=sort(Rx,2);
                    Rx=unique(Rx,'rows');
                    E{i}=Rx; %E has max length=div
                end
                
                parfor h=1:A(p)
                    mN_2_4_6_8=zeros(S,N,N);
                    F=cell2mat(E(h));
                    for s=1:S
                        temp2=squeeze(Louv(s,:,:));% Time and community
                        
                        for m=1:length(F(:,1))
                            if ~ismember(R6_IDXsU(h+B(p),:),F(m,:)) % if area is not already part of the seed
                                
                                
                                InCom=0;
                                assembly=[R6_IDXsU(h+B(p),:) F(m,:)];
                                y=temp2(:,assembly);
                                z=prod(y,2); %Row with same community nr have product 1^8=1,2^8=256,3^8=6561
                                zz=(z==[1 256 6561]);
                                zzz=(zz(:,1)+zz(:,2)+zz(:,3));
                                InCom=numel(nonzeros(zzz));
                                mN_2_4_6_8(s,F(m,1),F(m,2))=InCom/Tmax;
                            end
                        end
                        
                        meanNet2_4_6_8(h,s,:,:)=mN_2_4_6_8(s,:,:);
                    end
                end
                
                MeanALLNetw2_4_6_8=zeros(A(p),N,N);
                if S>1
                    for h=1:A(p)
                        
                        G=squeeze(meanNet2_4_6_8(h,:,:,:));
                        MeanALLNetw2_4_6_8(h,:,:)=mean(G);
                        
                    end
                else
                    G=squeeze(meanNet2_4_6_8(:,1,:,:));
                    MeanALLNetw2_4_6_8(:,:,:)=G;
                end
                save(sprintf('P%d_MeanAllSNC8_%d',K,p),'MeanALLNetw2_4_6_8','-v7.3')
                
                
            end
            
        else
            
            for p=1:Y
                clear  TempMeanALL E  MeanALLNetw2_4_6_8 
                meanNet2_4_6_8=zeros(A(p),S,N,N);
                if p< Y
                    R6ORDtemp=R6_IDXsALL_ORD(I6(B(p)+1:B(p+1)));
                else
                    R6ORDtemp=R6_IDXsALL_ORD(I6(B(p)+1:B(p)+A(p)));
                end
                
                TempMeanALL=zeros(length(R6ORDtemp),N,N);
                for i=1:length(R6ORDtemp)
                    MeanALLNetw2_4_6=[];
                    F=R6ORDtemp(i);
                    b1=floor(F/Div);
                    b2=mod(F,Div);
                    if b2==0
                        load(sprintf('P%d_MeanAllSNC6_%d',K,b1))
                        TempMeanALL(i,:,:)=MeanALLNetw2_4_6(Div,:,:);
                    else
                        load(sprintf('P%d_MeanAllSNC6_%d',K,b1+1))
                        TempMeanALL(i,:,:)=MeanALLNetw2_4_6(b2,:,:);
                    end
                    %Identify pairs with a qint above threshold a to search through
                    temp=squeeze(TempMeanALL(i,:,:));
                    [row col]=find(temp>thr);
                    Rx=[row col];
                    Rx=sort(Rx,2);
                    Rx=unique(Rx,'rows');
                    E{i}=Rx; %E has max length=div
                end
                
                parfor h=1:A(p)
                    mN_2_4_6_8=zeros(S,N,N);
                    F=cell2mat(E(h));
                    for s=1:S
                        temp2=squeeze(Louv(s,:,:));% Time and community
                        
                        for m=1:length(F(:,1))
                            if ~ismember(R6_IDXsU(h+B(p),:),F(m,:)) % if area is not already part of the seed
                                
                                
                                InCom=0;
                                assembly=[R6_IDXsU(h+B(p),:) F(m,:)];
                                y=temp2(:,assembly);
                                z=prod(y,2); %Row with same community nr have product 1^7=1,2^7=128,3^8=2187
                                zz=(z==[1 128 2187]);
                                zzz=(zz(:,1)+zz(:,2)+zz(:,3));
                                InCom=numel(nonzeros(zzz));
                                mN_2_4_6_8(s,F(m,1),F(m,2))=InCom/Tmax;
                            end
                        end
                        
                        meanNet2_4_6_8(h,s,:,:)=mN_2_4_6_8(s,:,:);
                    end
                end
                
                MeanALLNetw2_4_6_8=zeros(A(p),N,N);
                if S>1
                    for h=1:A(p)
                        
                        G=squeeze(meanNet2_4_6_8(h,:,:,:));
                        MeanALLNetw2_4_6_8(h,:,:)=mean(G);
                        
                    end
                else
                    G=squeeze(meanNet2_4_6_8(:,1,:,:));
                    MeanALLNetw2_4_6_8(:,:,:)=G;
                end
                save(sprintf('P%d_MeanAllSNC8_%d',K,p),'MeanALLNetw2_4_6_8','-v7.3')
                
                
            end
            
            
        end
        
        
        DnumelAll=[];
        for p=1:Y
            load(sprintf('P%d_MeanAllSNC8_%d',K,p))
            
            R8_IDXt=zeros(A(p),K+add,8);
            
            
            
            Freq_Y=zeros(A(p),K+add);
            
            for h=1:A(p)
                R8_t=zeros(K+add,8);
                
                D=squeeze(MeanALLNetw2_4_6_8(h,:,:));
                Dnumel(h)=numel(nonzeros(D));
                if Dnumel(h)>0 && Dnumel(h)>=K
                    a=maxk(nonzeros(D),K);
                    [row col]=find(D>=a(K));
                    R8_t(1:numel(row),:)=[R6_IDXsU(B(p)+h,:).*ones(numel(row),6) row col];
                    Freq_Y(h,1:numel(row))=D(D>=a(K));
                    
                elseif Dnumel(h)>0
                    [row col]=find(D>0);
                    R8_t(1:numel(row),:)=[R6_IDXsU(B(p)+h,:).*ones(numel(row),6) row col]; % The SNCs
                    Freq_Y(h,1:numel(row))=D(D>0); % Stores qint of the SNCs
                    
                else
                    R8_t(1,:)=[R6_IDXsU(B(p)+h,:) 0 0];
                    Freq_Y(h,1)=0;
                    
                end
                R8_IDXt(h,:,:)=R8_t(:,:);
               
            end
            save(sprintf('P%d_R8prob%d',K,p),'Freq_Y', 'R8_IDXt','-v7.3')
            
            DnumelAll=[DnumelAll Dnumel];
        end
        save SizeS3  DnumelAll
        
        %Add to one matrix
        R8_IDX=[];Freq=[];   
        for p=1:Y
            
            load(sprintf('P%d_R8prob%d',K,p))
            
            R8_IDX=[R8_IDX; R8_IDXt];

            Freq=[Freq; Freq_Y];
            save(sprintf('P%d_R8idx',K), 'R8_IDX',  'Freq')
            clear R R8_IDXt Freq_Y
            delete(sprintf('P%d_R8prob%d.mat',K,p))
            
        end
        
        %Sort each SNC from lowest to highest area nr in order to later
        %delete dublicate rows that have converged in this  step
        R8_IDXs=[];
        R8_IDXsALL=[];
        R8_IDXsALL_ORD=[];
        R8_Freq_ALL=[];
        
        q=mod(length(R8_IDX),Div)
        qa=floor(length(R8_IDX)/Div) %Rows in each division
        
        
        if qa>1000
            Av=0;
            for i=1:Div
                Av(i+1,1)=qa*i;
            end
            
            Le=ones(Div,1)*qa;
            if q>0
                Av(Div+1)=q;
                Le=[Le; q];
                Di=Div+1; %Number of divisions
            else
                Di=Div;
            end
            
           
            for g=1:Di
                g
                temp=[];
                R8_IDXsALLtemp=[];
                R8_IDXsALL_ORDtemp=[];
                Freq_ALLtemp=[];
                for h=1:Le(g)
                    
                    for k=1:length(R8_IDX(h+Av(g),:,1))
                        temp(h,k,:)=sort(R8_IDX(h+Av(g),k,:));
                    end
                    H_i=[]; I=[];Nu=[];
                    T=squeeze(temp(h,:,:));
                    [H_i I]=unique(T,'rows','stable');
                    
                    Nu=ones(numel(I),1).*h;
                    
                    Freq_ALLtemp=[Freq_ALLtemp; Freq(h+Av(g),I)'];
                    R8_IDXsALLtemp=[R8_IDXsALLtemp;H_i];
                    R8_IDXsALL_ORDtemp=[R8_IDXsALL_ORDtemp;Nu];
                    
                end
                
                R8_IDXs=[R8_IDXs; temp];
                R8_IDXsALL=[R8_IDXsALL;R8_IDXsALLtemp];
                R8_IDXsALL_ORD=[R8_IDXsALL_ORD;R8_IDXsALL_ORDtemp];
                R8_Freq_ALL=[R8_Freq_ALL;Freq_ALLtemp];
                
            end 
            
        else
            for h=1:length(R8_IDX)
                for k=1:length((R8_IDX(h,:,1)))
                    R8_IDXs(h,k,:)=sort(R8_IDX(h,k,:));
                end
                H_i=[]; I=[];Nu=[];
                T=squeeze(R8_IDXs(h,:,:));
                [H_i I]=unique(T,'rows','stable');
                
                Nu=ones(numel(I),1).*h;
                
                R8_Freq_ALL=[R8_Freq_ALL; Freq(h,I)'];
                R8_IDXsALL=[R8_IDXsALL;H_i];
                R8_IDXsALL_ORD=[R8_IDXsALL_ORD;Nu];
                
            end
            
        end
        
      
        %delete rows with zeros
        [ro col]=find(R8_IDXsALL==0);
        zeRO=unique(ro); numel(zeRO)
        R8_IDXsALL(zeRO,:)=[];
        R8_IDXsALL_ORD(zeRO)=[];
        R8_Freq_ALL(zeRO)=[];
        
        
        [R8_IDXsU I8]=unique(R8_IDXsALL,'rows');
        FreqIDXsU8=R8_Freq_ALL(I8);
        R8_IDXsUfreq=[FreqIDXsU8 R8_IDXsU];
        R8_IDXsUfreqSort=sortrows(R8_IDXsUfreq,1,'descend');
        save(sprintf('P%d_I8',K),'R8_IDXsALL_ORD', 'R8_IDXsUfreqSort', 'R8_IDXsU', 'FreqIDXsU8', 'R8_IDXs', 'I8',  'R8_IDXsALL', 'R8_IDXsUfreq','-v7.3') 
        
        %Dubbelcheck there are no area duplicates
        duplicate=0;
        for hh=1:length(R8_IDXsU)
            U(hh)=numel(unique(R8_IDXsU(hh,:)));
            if evenodd==0
                if U(hh)<8
                    hh
                    duplicate=1
                    save(sprintf('P%d_8U',K),'U')
                    break
                end
            else
                if U(hh)<7
                    hh
                    duplicate=1
                    save(sprintf('P%d_8U',K),'U')
                    break
                end
            end
        end
        
        fileID = fopen('Step3.txt','w');
        fprintf(fileID,'Number of rows in step 4 is %d\n',size(R8_IDXsU,1));
        fprintf(fileID,'area dubplicate in component %d \n', duplicate);
        fprintf(fileID, 'Min search space %d ', min(DnumelAll));
        fprintf(fileID,'Number of divisions %d', Y);
        fclose(fileID);
        
        
        if SNC_size>8
            
            
            Div=100
            X=floor(numel(I8)/Div)
            B=[Div:Div:X*Div];
            B=[0 B];
            
            A=ones(1,X)*Div;
            if (length(I8)-X*Div)>0 % If A is not an even multiplication of div
                A=[A numel(I8)-X*Div];
            else
            end
            Y=numel(A)
            %
            if evenodd==0
                for  p=1:Y
                    clear  TempMeanALL E MeanALLNetw2_4_6_8_10 
                    meanNet2_4_6_8_10=zeros(A(p),S,N,N);
                    
                    if p< Y
                        R8ORDtemp=R8_IDXsALL_ORD(I8(B(p)+1:B(p+1)));
                        
                    else
                        R8ORDtemp=R8_IDXsALL_ORD(I8(B(p)+1:B(p)+A(p)));
                    end
                    
                    
                    
                    TempMeanALL=zeros(length(R8ORDtemp),N,N);
                    for i=1:length(R8ORDtemp)
                        MeanALLNetw2_4_6_8=[];
                        F=R8ORDtemp(i);
                        b1=floor(F/Div);
                        b2=mod(F,Div);
                        if b2==0
                            load(sprintf('P%d_MeanAllSNC8_%d',K,b1))
                            TempMeanALL(i,:,:)=MeanALLNetw2_4_6_8(Div,:,:);
                        else
                            load(sprintf('P%d_MeanAllSNC8_%d',K,b1+1))
                            TempMeanALL(i,:,:)=MeanALLNetw2_4_6_8(b2,:,:);
                        end
                        %Identify pairs with a qint above threshold a to search through
                        temp=squeeze(TempMeanALL(i,:,:));
                        [row col]=find(temp>thr);
                        Rx=[row col];
                        Rx=sort(Rx,2);
                        Rx=unique(Rx,'rows');
                        E{i}=Rx; %E has max length=div
                    end
                    
                    
                    parfor h=1:A(p) 
                        mN_2_4_6_8_10=zeros(S,N,N);
                        F=cell2mat(E(h));
                        for s=1:S
                            temp2=squeeze(Louv(s,:,:));% Time and community
                            
                            for m=1:length(F(:,1))
                                if ~ismember(R8_IDXsU(h+B(p),:),F(m,:))
                                    InCom=0;
                                    assembly=[R8_IDXsU(h+B(p),:) F(m,:)];
                                    y=temp2(:,assembly);
                                    z=prod(y,2);
                                    zz=(z==[1 1024 59049]);
                                    zzz=(zz(:,1)+zz(:,2)+zz(:,3));
                                    InCom=numel(nonzeros(zzz));
                                    mN_2_4_6_8_10(s,F(m,1),F(m,2))=InCom/Tmax;
                                end
                            end
                            
                            meanNet2_4_6_8_10(h,s,:,:)=mN_2_4_6_8_10(s,:,:);
                            
                        end
                    end
                    
                    MeanALLNetw2_4_6_8_10=zeros(A(p),N,N);
                    if S>1
                        for h=1:A(p)
                            G=squeeze(meanNet2_4_6_8_10(h,:,:,:));
                            MeanALLNetw2_4_6_8_10(h,:,:)=mean(G);
                            
                        end
                        
                    else
                        G=squeeze(meanNet2_4_6_8_10(:,1,:,:));
                        MeanALLNetw2_4_6_8_10(:,:,:)=G;
                    end
                    save(sprintf('P%d_MeanAllSNC10_%d',K,p),'MeanALLNetw2_4_6_8_10','-v7.3')
                    
                end
                
            else
                for  p=1:Y
                    clear  TempMeanALL E MeanALLNetw2_4_6_8_10  
                    meanNet2_4_6_8_10=zeros(A(p),S,N,N);
                    
                    if p< Y
                        R8ORDtemp=R8_IDXsALL_ORD(I8(B(p)+1:B(p+1)));
                        
                    else
                        R8ORDtemp=R8_IDXsALL_ORD(I8(B(p)+1:B(p)+A(p)));
                    end
                    
                    
                    
                    TempMeanALL=zeros(length(R8ORDtemp),N,N);
                    for i=1:length(R8ORDtemp)
                        MeanALLNetw2_4_6_8=[];
                        F=R8ORDtemp(i);
                        b1=floor(F/Div);
                        b2=mod(F,Div);
                        if b2==0
                            load(sprintf('P%d_MeanAllSNC8_%d',K,b1))
                            TempMeanALL(i,:,:)=MeanALLNetw2_4_6_8(Div,:,:);
                        else
                            load(sprintf('P%d_MeanAllSNC8_%d',K,b1+1))
                            TempMeanALL(i,:,:)=MeanALLNetw2_4_6_8(b2,:,:);
                        end
                        %Identify pairs with a qint above threshold a to search through
                        temp=squeeze(TempMeanALL(i,:,:));
                        [row col]=find(temp>thr);
                        Rx=[row col];
                        Rx=sort(Rx,2);
                        Rx=unique(Rx,'rows');
                        E{i}=Rx; %E has max length=div
                    end
                    
                    
                    parfor h=1:A(p) 
                        mN_2_4_6_8_10=zeros(S,N,N);
                        F=cell2mat(E(h));
                        for s=1:S
                            temp2=squeeze(Louv(s,:,:));% Time and community
                            
                            for m=1:length(F(:,1))
                                if ~ismember(R8_IDXsU(h+B(p),:),F(m,:))
                                    InCom=0;
                                    assembly=[R8_IDXsU(h+B(p),:) F(m,:)];
                                    y=temp2(:,assembly);
                                    z=prod(y,2); %Row with same community nr have product 1^9=1,2^9=512,3^9=19683
                                    zz=(z==[1 512 19683]);
                                    zzz=(zz(:,1)+zz(:,2)+zz(:,3));
                                    InCom=numel(nonzeros(zzz));
                                    mN_2_4_6_8_10(s,F(m,1),F(m,2))=InCom/Tmax;
                                end
                            end
                            
                            meanNet2_4_6_8_10(h,s,:,:)=mN_2_4_6_8_10(s,:,:);
                            
                        end
                    end
                    
                    MeanALLNetw2_4_6_8_10=zeros(A(p),N,N);
                    if S>1
                        for h=1:A(p)
                            G=squeeze(meanNet2_4_6_8_10(h,:,:,:));
                            MeanALLNetw2_4_6_8_10(h,:,:)=mean(G);
                            
                        end
                        
                    else
                        G=squeeze(meanNet2_4_6_8_10(:,1,:,:));
                        MeanALLNetw2_4_6_8_10(:,:,:)=G;
                    end
                    save(sprintf('P%d_MeanAllSNC10_%d',K,p),'MeanALLNetw2_4_6_8_10','-v7.3')
                    
                end

            end
            
            SizeSearch4=[];
            DnumelAll=[];
            for p=1:Y
                load(sprintf('P%d_MeanAllSNC10_%d',K,p))
                
                R10_IDXt=zeros(A(p),K+add,10);
                
            
                
                Freq_Y=zeros(A(p),K+add);
               
                for h=1:A(p)
                    R10_t=zeros(K+add,10);
                    
                    
                    D=squeeze(MeanALLNetw2_4_6_8_10(h,:,:));
                    Dnumel(h)=numel(nonzeros(D));
                    if Dnumel(h)>0 && Dnumel(h)>=K
                        a=maxk(nonzeros(D),K);
                        [row col]=find(D>=a(K));
                        R10_t(1:numel(row),:)=[R8_IDXsU(B(p)+h,:).*ones(numel(row),8) row col];
                        Freq_Y(h,1:numel(row))=D(D>=a(K));
                        
                    elseif Dnumel(h)>0
                        [row col]=find(D>0);
                        R10_t(1:numel(row),:)=[R8_IDXsU(B(p)+h,:).*ones(numel(row),8) row col];
                        Freq_Y(h,1:numel(row))=D(D>0);
                        
                    else
                        R10_t(1,:)=[R8_IDXsU(B(p)+h,:) 0 0];
                        Freq_Y(h,1)=0;
                        
                    end
                    R10_IDXt(h,:,:)=R10_t(:,:);
                    
                    
                    
                end
                save(sprintf('P%d_R10prob%d',K,p),'SizeSearch4','Freq_Y', 'R10_IDXt', '-v7.3') 
                
                DnumelAll=[DnumelAll Dnumel];
            end
            save SizeS4  DnumelAll
            %Add to one matrix
            R10_IDX=[];R10_length=[]; Freq=[];
            for p=1:Y
                
                load(sprintf('P%d_R10prob%d',K,p))
                
                R10_IDX=[R10_IDX; R10_IDXt];
                
                
                Freq=[Freq; Freq_Y];
                save(sprintf('P%d_R10idx',K), 'R10_IDX', 'Freq') 
                clear  R10_IDXt Freq
                delete(sprintf('P%d_R10prob%d.mat',K,p))
            end
            
            
            for h=1:length(R10_IDX)
                for k=1:length((R10_IDX(h,:,1)))
                    R10_IDXs(h,k,:)=sort(R10_IDX(h,k,:));
                end
            end
            save(sprintf('P%d_R10idxS',K), 'R10_IDXs')
            
            R10_IDXsALL=[];
            R10_IDXsALL_ORD=[];
            R10_Freq_ALL=[];
            
            for h=1:length(nonzeros(R10_IDXs(:,1,1)))
                T=squeeze(R10_IDXs(h,:,:));
                [H_i I]=unique(T,'rows','stable');
                
                Nu=ones(numel(I),1).*h;
                
                R10_Freq_ALL=[R10_Freq_ALL; Freq(h,I)'];
                R10_IDXsALL=[R10_IDXsALL;H_i];
                R10_IDXsALL_ORD=[R10_IDXsALL_ORD;Nu];
                
            end
            
            
            
            %delete rows with zeros
            [ro col]=find(R10_IDXsALL==0);
            zeRO=unique(ro); numel(zeRO)
            R10_IDXsALL(zeRO,:)=[];
            R10_IDXsALL_ORD(zeRO)=[];
            R10_Freq_ALL(zeRO)=[];
            
            [R10_IDXsU I10]=unique(R10_IDXsALL,'rows');
            FreqIDXsU10=R10_Freq_ALL(I10);
            R10_IDXsUfreq=[FreqIDXsU10 R10_IDXsU];

            R10_IDXsUfreqSort=sortrows(R10_IDXsUfreq,1,'descend');
            
            save(sprintf('P%d_I10',K), 'R10_IDXsUfreqSort', 'R10_IDXsU', 'FreqIDXsU10', 'R10_IDXs', 'I10',  'R10_IDXsALL', 'R10_IDXsUfreq','-v7.3')
            
            %Dubbelcheck there are no area duplicates  
            duplicate=0;
            for hh=1:length(R10_IDXsU)
                U(hh)=numel(unique(R10_IDXsU(hh,:)));
                if evenodd==0
                    if U(hh)<10
                        hh
                        duplicate=1
                        save P10_10U U
                        break
                    end
                else
                    if U(hh)<9
                        hh
                        duplicate=1
                        save P10_10U U
                        break
                    end
                end
            end
            
            fileID = fopen('Step8_10.txt','w');
            fprintf(fileID,'Number of rows in step 4 is %d\n',size(R10_IDXsU,1));
            fprintf(fileID,'area dubplicate in component %d \n', duplicate);
            fprintf(fileID, 'Min search space %d ', min(DnumelAll));
            fprintf(fileID,'Number of divisions %d', Y);
            fclose(fileID);
            
            if SNC_size==12
                
                X=floor(numel(I10)/Div)
                B=[Div:Div:X*Div];
                B=[0 B];
                
                A=ones(1,X)*Div;
                if (length(I10)-X*Div)>0 % If A is not an even multiplication of div
                    A=[A numel(I10)-X*Div];
                else
                end
                Y=numel(A)
                %
                if evenodd==0
                    for  p=1:Y
                        clear  TempMeanALL E meanNet2_4_6_8_10_12
                        meanNet2_4_6_8_10_12=zeros(A(p),S,N,N);
                        if p< Y
                            R10ORDtemp=R10_IDXsALL_ORD(I10(B(p)+1:B(p+1)));
                            
                        else
                            R10ORDtemp=R10_IDXsALL_ORD(I10(B(p)+1:B(p)+A(p)));
                        end
                        
                        
                        TempMeanALL=zeros(length(R10ORDtemp),N,N);
                        for i=1:length(R10ORDtemp)
                            MeanALLNetw2_4_6_8_10=[];
                            F=R10ORDtemp(i);
                            b1=floor(F/Div);
                            b2=mod(F,Div);
                            if b2==0
                                load(sprintf('P%d_MeanAllSNC10_%d',K,b1))
                                TempMeanALL(i,:,:)=MeanALLNetw2_4_6_8_10(Div,:,:);
                            else
                                load(sprintf('P%d_MeanAllSNC10_%d',K,b1+1))
                                TempMeanALL(i,:,:)=MeanALLNetw2_4_6_8_10(b2,:,:);
                            end
                            %Identify pairs with a qint above threshold a to search through
                            temp=squeeze(TempMeanALL(i,:,:));
                            [row col]=find(temp>thr);
                            Rx=[row col];
                            Rx=sort(Rx,2);
                            Rx=unique(Rx,'rows');
                            E{i}=Rx; %E has max length=div
                        end
                        
                      
                        
                        parfor h=1:A(p)
                            mN_2_4_6_8_10_12=zeros(S,N,N);
                            F=cell2mat(E(h));
                            for s=1:S
                                temp2=squeeze(Louv(s,:,:));% Time and community
                                
                                for m=1:length(F(:,1))
                                    if ~ismember(R10_IDXsU(h+B(p),:),F(m,:))
                                        
                                        InCom=0;
                                        assembly=[R10_IDXsU(h+B(p),:) F(m,:)];
                                        y=temp2(:,assembly);
                                        z=prod(y,2); 
                                        zz=(z==[1 4096 531441]);
                                        zzz=(zz(:,1)+zz(:,2)+zz(:,3));
                                        InCom=numel(nonzeros(zzz));
                                        mN_2_4_6_8_10_12(s,F(m,1),F(m,2))=InCom/Tmax;
                                        
                                    end
                                end
                                
                                meanNet2_4_6_8_10_12(h,s,:,:)=mN_2_4_6_8_10_12(s,:,:);
                                
                            end
                        end
                        
                        MeanALLNetw2_4_6_8_10_12=zeros(A(p),N,N);
                        if S>1
                            for h=1:A(p)
                                G=squeeze(meanNet2_4_6_8_10_12(h,:,:,:));
                                MeanALLNetw2_4_6_8_10_12(h,:,:)=mean(G);
                                
                            end
                            
                        else
                            G=squeeze(meanNet2_4_6_8_10_12(:,1,:,:));
                            MeanALLNetw2_4_6_8_10_12(:,:,:)=G;
                        end
                        save(sprintf('P%d_MeanAllSNC12_%d',K,p),'MeanALLNetw2_4_6_8_10_12','-v7.3')
                    end
                    
                else
                    
                    for  p=1:Y
                        clear  TempMeanALL E meanNet2_4_6_8_10_12
                        meanNet2_4_6_8_10_12=zeros(A(p),S,N,N);
                        if p< Y
                            R10ORDtemp=R10_IDXsALL_ORD(I10(B(p)+1:B(p+1)));
                            
                        else
                            R10ORDtemp=R10_IDXsALL_ORD(I10(B(p)+1:B(p)+A(p)));
                        end
                        
                        
                        TempMeanALL=zeros(length(R10ORDtemp),N,N);
                        for i=1:length(R10ORDtemp)
                            MeanALLNetw2_4_6_8_10=[];
                            F=R10ORDtemp(i);
                            b1=floor(F/Div);
                            b2=mod(F,Div);
                            if b2==0
                                load(sprintf('P%d_MeanAllSNC10_%d',K,b1))
                                TempMeanALL(i,:,:)=MeanALLNetw2_4_6_8_10(Div,:,:);
                            else
                                load(sprintf('P%d_MeanAllSNC10_%d',K,b1+1))
                                TempMeanALL(i,:,:)=MeanALLNetw2_4_6_8_10(b2,:,:);
                            end
                            %Identify pairs with a qint above threshold a to search through
                            temp=squeeze(TempMeanALL(i,:,:));
                            [row col]=find(temp>thr);
                            Rx=[row col];
                            Rx=sort(Rx,2);
                            Rx=unique(Rx,'rows');
                            E{i}=Rx; %E has max length=div
                        end
                        
                        % end
                        
                        
                        parfor h=1:A(p)
                            mN_2_4_6_8_10_12=zeros(S,N,N);
                            F=cell2mat(E(h));
                            for s=1:S
                                temp2=squeeze(Louv(s,:,:));% Time and community
                                
                                for m=1:length(F(:,1))
                                    if ~ismember(R10_IDXsU(h+B(p),:),F(m,:))
                                        
                                        InCom=0;
                                        assembly=[R10_IDXsU(h+B(p),:) F(m,:)];
                                        y=temp2(:,assembly);
                                        z=prod(y,2); %Row with same community nr have product 1^11=1,2^11=2048,3^11=177147
                                        zz=(z==[1 2048 177147]);
                                        zzz=(zz(:,1)+zz(:,2)+zz(:,3));
                                        InCom=numel(nonzeros(zzz));
                                        mN_2_4_6_8_10_12(s,F(m,1),F(m,2))=InCom/Tmax;
                                        
                                    end
                                end
                                
                                meanNet2_4_6_8_10_12(h,s,:,:)=mN_2_4_6_8_10_12(s,:,:);
                                
                            end
                        end
                        
                        MeanALLNetw2_4_6_8_10_12=zeros(A(p),N,N);
                        if S>1
                            for h=1:A(p)
                                G=squeeze(meanNet2_4_6_8_10_12(h,:,:,:));
                                MeanALLNetw2_4_6_8_10_12(h,:,:)=mean(G);
                            end 
                        else
                            G=squeeze(meanNet2_4_6_8_10_12(:,1,:,:));
                            MeanALLNetw2_4_6_8_10_12(:,:,:)=G;
                        end
                        save(sprintf('P%d_MeanAllSNC12_%d',K,p),'MeanALLNetw2_4_6_8_10_12','-v7.3')
                    end
                   
                end
               
                DnumelAll=[];
                for p=1:Y
                    load(sprintf('P%d_MeanAllSNC12_%d',K,p))
                    
                    R12_IDXt=zeros(A(p),K+add,12);
                    
                   
                    
                    Freq_Y=zeros(A(p),K+add);
                    SizeSearch5=[];
                    for h=1:A(p)
                        R12_t=zeros(K+add,12);
                        
                        D=squeeze(MeanALLNetw2_4_6_8_10_12(h,:,:));
                        Dnumel(h)=numel(nonzeros(D));
                        if Dnumel(h)>0 && Dnumel(h)>=K
                            a=maxk(nonzeros(D),K);
                            [row col]=find(D>=a(K));
                            R12_t(1:numel(row),:)=[R10_IDXsU(B(p)+h,:).*ones(numel(row),10) row col];
                            Freq_Y(h,1:numel(row))=D(D>=a(K));
                            
                        elseif Dnumel(h)>0
                            [row col]=find(D>0);
                            R12_t(1:numel(row),:)=[R10_IDXsU(B(p)+h,:).*ones(numel(row),10) row col];
                            Freq_Y(h,1:numel(row))=D(D>0);
                            
                        else
                            R12_t(1,:)=[R10_IDXsU(B(p)+h,:) 0 0];
                            Freq_Y(h,1)=0;
                            
                        end
                        R12_IDXt(h,:,:)=R12_t(:,:);
                        
                        
                        
                    end
                    save(sprintf('P%d_R12prob%d',K,p),'SizeSearch5','Freq_Y', 'R12_IDXt', '-v7.3')
                    
                    DnumelAll=[DnumelAll Dnumel];
                end
                save SizeS5  DnumelAll
                %Add to one matrix
                R12_IDX=[];Freq=[];
                for p=1:Y
                    
                    load(sprintf('P%d_R12prob%d',K,p))
                    
                    R12_IDX=[R12_IDX; R12_IDXt];
                    
                    
                    Freq=[Freq; Freq_Y];
                    save(sprintf('P%d_R12idx',K), 'R12_IDX',   'Freq')
                    clear R R12_IDXt  Freq_Y
                end
                
                %Sort each assembly from lowest to highest
                for h=1:length(R12_IDX)
                    for k=1:length((R12_IDX(h,:,1)))
                        R12_IDXs(h,k,:)=sort(R12_IDX(h,k,:));
                    end
                end
                save(sprintf('P%d_R12idxS',K), 'R12_IDXs')
                
                R12_IDXsALL=[];
                R12_IDXsALL_ORD=[];
                R12_Freq_ALL=[];
                
                for h=1:length(nonzeros(R12_IDXs(:,1,1)))
                    T=squeeze(R12_IDXs(h,:,:));
                    [H_i I]=unique(T,'rows','stable');
                    
                    Nu=ones(numel(I),1).*h;
                    
                    R12_Freq_ALL=[R12_Freq_ALL; Freq(h,I)'];
                    R12_IDXsALL=[R12_IDXsALL;H_i];
                    R12_IDXsALL_ORD=[R12_IDXsALL_ORD;Nu];
                    
                end
                
                
                %delete rows with zeros
                [ro col]=find(R12_IDXsALL==0);
                zeRO=unique(ro); numel(zeRO)
                R12_IDXsALL(zeRO,:)=[];
                R12_IDXsALL_ORD(zeRO)=[];
                R12_Freq_ALL(zeRO)=[];
                
                
                [R12_IDXsU I12]=unique(R12_IDXsALL,'rows');
                FreqIDXsU12=R12_Freq_ALL(I12);
                R12_IDXsUfreq=[FreqIDXsU12 R12_IDXsU];
                
                
                R12_IDXsUfreqSort=sortrows(R12_IDXsUfreq,1,'descend');
                save(sprintf('P%d_I12',K), 'R12_IDXsUfreqSort', 'R12_IDXsU', 'FreqIDXsU12', 'R12_IDXs', 'I12', 'R12_IDXsALL_ORD', 'R12_IDXsALL', 'R12_IDXsUfreq','-v7.3')
                
                %Dubbelcheck there are no area duplicates  
                duplicate=0;
                for hh=1:length(R12_IDXsU)
                    U(hh)=numel(unique(R12_IDXsU(hh,:)));
                    if evenodd==0
                        if U(hh)<12
                            hh
                            duplicate=1
                            save P10_12U U
                            break
                        end
                    else
                        if U(hh)<11
                            hh
                            duplicate=1
                            save P10_12U U
                            break
                        end
                        
                    end
                end
                
                fileID = fopen('Step5.txt','w');
                fprintf(fileID,'Number of rows in step 4 is %d\n',size(R12_IDXsU,1));
                fprintf(fileID,'area dubplicate in component %d \n', duplicate);
                fprintf(fileID,'Min search space %d ', min(DnumelAll));
                fprintf(fileID,'Number of divisions %d', Y);
                fclose(fileID);
            else
                
            end
        else
            
            
        end
    else
      
    end
else
   
end

end
