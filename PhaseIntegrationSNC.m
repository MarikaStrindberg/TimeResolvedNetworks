function []=PhaseIntegrationSNC(AX,Tmax,Louv,TR,S,N,K,L,SNC_size)

% Marika Strindberg, Karolinska Institutet 2018-2020
% marika.strindberg@ki.se, marikastrindberg@gmail.com

% Subrutine to TimeResolvedNets that calculates dPC-values for a sample of
% SNCs and random komponents prior to and during integration (assignement to the same community)


NA=N;

if AX==1
    
    %Draw 500 random SNCs
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
    
    i=length(H);
    sncSample=round(1+(i-1).*rand(L,1),0);
    SNCSample=H(sncSample,2:SNC_size +1);
    QintSample=H(sncSample,1);
    save SNCSamples SNCSample sncSample QintSample
    
    %Get the timeseries
    NN=SNCSample;
   
        Qint_SNC=[];
        SNCSamplebinTS=[]; SNCSamplecomTS =[];
        SNCSamplebinTS=zeros(S,L,Tmax);
        SNCSamplecomTS=zeros(S,L,Tmax);
        for k=1:L
            
            snc=NN(k,:); %the SNC areas
            
            for s=1:S
                temp2=squeeze(Louv(s,:,:));% Time and community
                InCom=0;
                y=temp2(:,snc);
                z=prod(y,2); %Row with same community nr have product
                if SNC_size==8
                    zz=(z==[1 256 6561]);
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
                SNCSamplebinTS(s,k,:)=zzz; %Binary timeseries
                temp3=temp2(:,snc(1)); % get the community nr of the SNC
                temp3(~logical(zzz))=0;
                SNCSamplecomTS(s,k,:)=temp3;% timeseries with community beloning
                InCom=numel(nonzeros(zzz));
                Qint_SNC(k,s)=InCom/Tmax;
            end
        end
   
    save SNCSampleTS SNCSamplecomTS SNCSamplebinTS
    
   %Create the timing matrix
   TempOr=zeros(S,L,Tmax);
   CT=SNCSamplebinTS;
    
    for w=1:L
        for s=1:S
            
            d=squeeze(CT(s,w,:));
            
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
    
    
    save  TempOR_SampleSNC TempOr
    
    
else
    
    % Random components
    n=SNC_size;
    randSNC=zeros(L,n);
    Qint_rand=zeros(L,S);
    MeanInt=[];
    p=0;
    for r=1:L
        r
        ran = round(1 + (NA-1).*rand(n,1),0);
        if numel(unique(ran))==n
            p=p+1
            randSNC(p,:)=ran;
            
            for s=1:S
                temp2=squeeze(Louv(s,:,:));% Time and community
                InCom=0;
                y=temp2(:,ran);
                z=prod(y,2);
                if SNC_size==8
                    zz=(z==[1 256 6561]);
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
                RandSNC_TS(s,p,:)=zzz; %Timeseries
                InCom=numel(nonzeros(zzz));
                Qint_rand(p,s)=InCom/Tmax;
                
            end
            MeanInt(p)=mean(Qint_rand(p,:));
            randSNC(p,:)=ran;
        end
        
    end
    
    MeanRand(1,1)=mean(MeanInt);
    MeanRand(1,2)=std(MeanInt);
    L=p;
    
    % Check that the random SNCs do not belong to the SNCs
    load(sprintf('P%d_I%d',K,SNC_size))
    if SNC_size==8||SNC_size==7
        load(sprintf('P%d_I8',K))
        H=R8_IDXsU;
    elseif SNC_size==6||SNC_size==5
        load(sprintf('P%d_I6',K))
        H=R6_IDXsU;
    elseif SNC_size==4||SNC_size==3
        load(sprintf('P%d_I4',K))
        H=R4_IDXsU;
    elseif SNC_size==10||SNC_size==9
        load(sprintf('P%d_I10',K))
        H=R10_IDXsU;
    elseif SNC_size==12||SNC_size==11
        load(sprintf('P%d_I12',K))
        H=R12_IDXsU;
    end
    
    
    TEMP=[H;randSNC];
    l=length(H);
    [RR I]=unique(TEMP,'rows');
    dubli=L-(length(RR)-length(H));
    rando=[l+1:1:l+L];
    RaSNC=intersect(rando,I);
    RaSNC=RaSNC-l;
    RandomSNC=randSNC(RaSNC,:);
    RandomSNC_TS=RandSNC_TS(:,RaSNC,:);
    NN=RandomSNC;
    
    save('RandSample','randSNC','MeanInt','Qint_rand','MeanRand','RandSNC_TS','InCom', 'RandomSNC','RandomSNC_TS','NN','-v7.3')
    
    %%Create timeing matrix of all random SNCs for all subject that track before during and after integration
    
    TempOr=zeros(S,L,Tmax);
    CT=RandomSNC_TS;
    
    for w=1:L
        for s=1:S
            
            d=squeeze(CT(s,w,:));
            
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
    
    save  TempOR_SampleRand TempOr
end



j=1;
y=18;
%x=20;
x=18
Z=[101:1:150];

AB=zeros(L,S,y,28,30); %from integration including disintegration
% nr of snc, subjects, pairdiff
ABasPeakdis=zeros(L,S,y,28,30);%only integration
Ba12=zeros(L,S,y,28,30);%preintegration
ABmean=zeros(L,S,y,30);
ABmedian=zeros(L,S,y,30);
Ba12mean=zeros(L,S,y,30);
Ba12median=zeros(L,S,y,30);
ABPeakmean=zeros(L,S,y,30);
ABPeakmedian=zeros(L,S,y,30);

for s=1:S
    s
    load(sprintf('iFCall_s%d',s))
    d=iFCalls; % Phase diffmatrix
    d=permute(d, [2 3 1]);
    
    for w=1:length(NN(:,1))
        
        
        N=NN(w,:);
        clear A
        M=zeros(NA,NA);
        b=1;
        
        for i=1:SNC_size
            for k=1:SNC_size
                if i~=k && M(N(k),N(i))==0
                    M(N(i),N(k))=1;
                    A(b,:)=[N(i) N(k)];
                    b=b+1;
                end
            end
        end
        
        q=squeeze(TempOr(:,w,:));   
        q(q>100)=0; %remove the indicators that counts timepoints after assembly disintegration
        m(w)=max(max(max(q))); %Max length of w among the kids
        if S>1
        r=max(q(s,:)); %Max length for individual child
        else
            r=m(w);
        end
        a=0;
        H=zeros(100,1);
        k=zeros(y,1);
        E=zeros(length(Z)+x,1); 
        z=x-m(w);
        
        for t=1:length(TempOr(1,1,:))
            
            if TempOr(s,w,t)>0 %
                
                for h=1:r
                    
                    if TempOr(s,w,t)==h
                        H(h)=H(h)+1;
                        for i=1:length(A)%for each cluster-node pair combo
                            AB(w,s,h,i,H(h))=d(A(i,1),A(i,2),t);
                            
                            ABasPeakdis(w,s,h,i,H(h))=d(A(i,1),A(i,2),t);%Use to calculte exact assembly peak and disembly phase diff
                        end
                        
                        ABmean(w,s,h,H(h))=mean(AB(w,s,h,:,H(h)));
                        ABmedian(w,s,h,H(h))=median(AB(w,s,h,:,H(h)));
                        
                        ABPeakmean(w,s,h,H(h))=mean(ABasPeakdis(w,s,h,:,H(h)));
                        ABPeakmedian(w,s,h,H(h))=median(ABasPeakdis(w,s,h,:,H(h)));
                    end
                end
                
                for h=r+1:x
                    if TempOr(s,w,t)==Z(h-r+1)
                        E(h-r)=E(h-r)+1;
                        for i=1:length(A)%for each SNC area pair
                            AB(w,s,h,i,E(h-r))=d(A(i,1),A(i,2),t);
                        end
                        ABmean(w,s,h,E(h-r))=mean(AB(w,s,h,:,E(h-r)));
                        ABmedian(w,s,h,E(h-r))=median(AB(w,s,h,:,E(h-r)));
                    end
                end
                %
            elseif TempOr(s,w,t)<0
                if  TempOr(s,w,t)==-18
                    u=1;
                elseif  TempOr(s,w,t)==-17
                    u=2;
                elseif  TempOr(s,w,t)==-16
                    u=3;
                elseif  TempOr(s,w,t)==-15
                    u=4;
                elseif  TempOr(s,w,t)==-14
                    u=5;
                elseif  TempOr(s,w,t)==-13
                    u=6;
                elseif  TempOr(s,w,t)==-12
                    u=7;
                elseif  TempOr(s,w,t)==-11
                    u=8;
                elseif  TempOr(s,w,t)==-10
                    u=9;
                elseif  TempOr(s,w,t)==-9
                    u=10;
                elseif  TempOr(s,w,t)==-8
                    u=11;
                elseif  TempOr(s,w,t)==-7
                    u=12;
                elseif  TempOr(s,w,t)==-6
                    u=13;
                elseif  TempOr(s,w,t)==-5
                    u=14;
                elseif  TempOr(s,w,t)==-4
                    u=15;
                    
                elseif  TempOr(s,w,t)==-3
                    u=16;
                elseif  TempOr(s,w,t)==-2
                    u=17;
                elseif  TempOr(s,w,t)==-1
                    u=18;
                end
                k(u)=k(u)+1;
                for i=1:length(A)%for each area pair comb
                    Ba12(w,s,u,i,k(u))=d(A(i,1),A(i,2),t);
                end
                Ba12mean(w,s,u,k(u))=mean(Ba12(w,s,u,:,k(u)));
                Ba12median(w,s,u,k(u))=median(Ba12(w,s,u,:,k(u)));
            else
            end
        end 
    end    
end 

if AX==1
    save PhaseDiffTS Ba12mean Ba12median Ba12 AB ABmean ABmedian  ABasPeakdis  ABPeakmean ABPeakmedian
    
else
    save PhaseDiffTS_rand Ba12mean Ba12median Ba12 AB ABmean ABmedian  ABasPeakdis  ABPeakmean ABPeakmedian
end

%Take mean accross SNCs
size(AB)
size(Ba12)


y=18

for w=1:L
    w
    for s=1:S
        
        gg=squeeze(ABmean(w,s,:,:));
        cc=squeeze(Ba12mean(w,s,:,:));
        dd=squeeze(ABPeakmean(w,s,:,:));
        gg2=squeeze(ABmedian(w,s,:,:));
        cc2=squeeze(Ba12median(w,s,:,:));
        dd2=squeeze(ABPeakmedian(w,s,:,:));
        
        for b=1:y %before integration
            
            MeanAB(w,s,b)=mean(nonzeros((cc(b,:)))); %insert abs if vanted
            MedianAB(w,s,b)=median(nonzeros((cc2(b,:)))); %
            MeanAB(isnan(MeanAB))=0;
            MedianAB(isnan(MedianAB))=0;
            
            
        end
        for h=1:x %from timepoint of integration
            %Calculate the mean phasediff i at that position
            
            MeanAB(w,s,h+y)=mean(nonzeros((gg(h,:))));%accomodate two preassembly
            MedianAB(w,s,h+y)=median(nonzeros((gg2(h,:)))); %
            MeanABPeak(w,s,h+y)=median(nonzeros((dd(h,:)))); %
            MedianABPeak(w,s,h+y)=median(nonzeros((dd2(h,:)))); %
            
        end
        
    end
end
MeanAB(isnan(MeanAB))=0;
MedianAB(isnan(MedianAB))=0;
MeanABPeak(isnan(MeanABPeak))=0;
MedianABPeak(isnan(MedianABPeak))=0;

if AX==1
    save PhaseDiffMean MeanAB MedianAB MedianABPeak MeanABPeak
else
    save PhaseDiffMean_Rand MeanAB MedianAB MedianABPeak MeanABPeak
end


for w=1:L
    for h=1:(x+y)%
        meanPh(w,h)=mean(nonzeros(MeanAB(w,:,h)));
        meamPhSd(w,h)=std(nonzeros(MeanAB(w,:,h)));
        medianPh(w,h)=mean(nonzeros(MedianAB(w,:,h)));
        medianPhSd(w,h)=std(nonzeros(MedianAB(w,:,h)));
    end   
end


figure
for w=1:L
    %plot(nonzeros(medianPh(w,:)),'linewidth',1.5)
    plot(nonzeros(meanPh(w,:)),'linewidth',1.5)
    hold on
end

yL = get(gca,'YLim');
line([19 19],yL,'Color','k','Linewidth',1);
xticks([1:6:(x+y)])
xticklabels([ round(-18*TR,1) round(-12*TR,1) round(-6*TR,1)  0 round(6*TR,1 ) round(12*TR,1)  round(18*TR,1)])


%Create line of mean +-1sd

for  h=1:(x+y)
    %
    meanLine(h,1)=mean(nonzeros(meanPh(:,h)));
    meanLine(h,2)=std(nonzeros(meanPh(:,h)));
    meanLine(h,3)=max(nonzeros(meanPh(:,h)));
    meanLine(h,4)=min(nonzeros(meanPh(:,h)));
    medianLine(h,1)=mean(nonzeros(medianPh(:,h)));
    medianLine(h,2)=std(nonzeros(medianPh(:,h)));
    medianLine(h,3)=max(nonzeros(medianPh(:,h)));
    medianLine(h,4)=min(nonzeros(medianPh(:,h)));
    
end


x=18;y=18;
figure; %mean Plot
plot(meanLine(:,1),'linewidth',2)
hold on;
plot(nonzeros(meanLine(:,1)+meanLine(:,2)),'linewidth', 2) %+1sd
hold on
plot(nonzeros(meanLine(:,1)-meanLine(:,2)),'linewidth', 2) %-1sd
hold on
plot(nonzeros(meanLine(:,3)),'linewidth', 2) %max 
hold on
plot(nonzeros(meanLine(:,4)),'linewidth', 2) %min

yL = get(gca,'YLim');
line([19 19],yL,'Color','k','Linewidth',1);

xticks([1:6:(x+y)])

xticklabels([ round(-18*TR,1) round(-12*TR,1) round(-6*TR,1)  0 round(6*TR,1 ) round(12*TR,1)  round(18*TR,1)])
xlabel('Time (sec)')
ylabel('Phase cohorence (cos(\Delta\phi)')



%median plot
figure;

plot(medianLine(:,1),'linewidth',2)
hold on;
plot(nonzeros(medianLine(:,1)+medianLine(:,2)),'linewidth', 2)
hold on
plot(nonzeros(medianLine(:,1)-medianLine(:,2)),'linewidth', 2)
hold on
plot(nonzeros(medianLine(:,3)),'linewidth', 2)
hold on
plot(nonzeros(medianLine(:,4)),'linewidth', 2)

yL = get(gca,'YLim');
line([19 19],yL,'Color','k','Linewidth',1);
xticks([1:6:(x+y)])

xticklabels([ round(-18*TR,1) round(-12*TR,1) round(-6*TR,1)  0 round(6*TR,1 ) round(12*TR,1)  round(18*TR,1)])


figure;  plot(medianLine(:,1),'linewidth',2)
hold on;
plot(meanLine(:,1),'linewidth',2)
hold on;
yL = get(gca,'YLim');
line([19 19],yL,'Color','k','Linewidth',1);
xticks([1:6:(x+y)])

xticklabels([ round(-18*TR,1) round(-12*TR,1) round(-6*TR,1)  0 round(6*TR,1 ) round(12*TR,1)  round(18*TR,1)])
legend('median',  'mean')
title('Phase integration')

MaxMedian=max(medianLine');
MaxMean=max(meanLine');

if AX==0
    save RandLines meanPh medianPh MaxMean MaxMedian meanLine medianLine
else
    save SNCLines meanPh medianPh MaxMean MaxMedian meanLine medianLine
end





%DUration of integration
if AX==1
    load TempOR_SampleSNC
else
    load TempOR_SampleRand
end

TempOr(TempOr<0)=0;
TempOr(TempOr>99)=0;
OnlyOneDur=0;
SNC=L

    for w=1:L
        for s=1:S
            
            D=squeeze(TempOr(s,w,:));
            b=unique(nonzeros(D)); %nr of different uniqe assembly durations
            c=histc(nonzeros(D),b); %instances per duration
            if numel(b)>1
                %for every unique length calculate the instance minus all instances with a duration
                % one shorter (could also be accomplised by dist function)
                for i=1:length(b)-1
                    ry(s,w,i)=c(i)-c(i+1);
                end
                ry(s,w,i+1)=c(i+1); %Vector with nr of assemblies with a particular length
                r=find(ry(s,w,:)>0);
                aa=squeeze(ry(s,w,1:max(r)));
                
                rr=[1:1:max(r)]'; %this will alwaysbe same as b
                Wid(w,s)=sum(aa.*rr)/sum(aa); %mean duration (duration multiplied by nr of instances)
                MaxWid(w,s)=max(r);
                
                %Calculate medium
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
        WidInt(w,2)=std(Wid(w,:));
        WidInt(w,3)=max(Wid(w,:));
        WidInt(w,4)=min(Wid(w,:));
        WidInt(w,5)=median(WidMedian(w,:));
        WidInt(w,6)=iqr(WidMedian(w,:))
        WidInt(w,7)=median(Wid(w,:));
        WidInt(w,8)=iqr(Wid(w,:));
    end

if AX==1
    save Duration500SNC WidInt TempOr Wid WidMedian 
else
    save Duration500rand WidInt TempOr Wid WidMedian 
end

mean(WidInt(:,1))*TR
std(WidInt(:,1))*TR


% relationship between mean duration of SNC and max phase coherence

MaxMeanPh=max(meanPh');


figure; scatter(WidInt(:,1)*TR,MaxMeanPh')

if AX==1
title('Mean Duration vs Max Phase coherence SNC')

else
  title('Mean Duration vs Max Phase coherence RandComp')  
end
ylabel('Max phase cohorence (cos(\Delta\phi) ')
xlabel('Mean duration (s)')
end