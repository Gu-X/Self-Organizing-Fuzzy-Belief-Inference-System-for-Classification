function [output]=SOFBISplus(input,mode,task)
disttype='minkowski';
P0=1;P1=2;
if strcmp(mode,'learning')==1
    gamma=input.gamma;
    data0=input.data;
    y0=input.y;
    ck0=input.chunksize;
    [L,W]=size(data0);
    seqck=[1:ck0:L,L+1];
    Lck=length(seqck);
    averdist=0;
    Npd=0;
    CEN=[];
    OUT=[];
    MEM=[];
    if strcmp(task,'c')==1
        lc=length(unique(y0));
    end
    for tt=1:1:(Lck-1)
        data=data0(seqck(tt):1:seqck(tt+1)-1,:);
        [data,temp]=unique(data,'rows');
        y=y0(seqck(tt):1:seqck(tt+1)-1,:);
        if strcmp(task,'c')==1
            y=full(ind2vec(y(temp,:)',lc))';
        elseif strcmp(task,'r')==1
            y=y(temp,:);
        end
        LC=length(y(:,1));
        dist00=pdist(data,disttype,P0).^P1;
        dist01=squareform(dist00);
        g0=mean(dist00);
        pp=1;
        F=[];
        tempseq4=[];
        while true
            dist00(dist00>mean(dist00))=[];
            tempB=mean(dist00);
            averdist1=tempB;
            tempdist=zeros(LC);
            tempdist(dist01<=averdist1)=1;
            tempdist0=tempdist;
            tempseq=sum(exp(-1*tempdist.*dist01./averdist1),2);
            tempdist=tempdist.*repmat(tempseq',LC,1);
            tempseq=max(tempdist,[],2);
            tempseq1=diag(tempdist);
            tempseq2=find(tempseq1-tempseq==0);
            purity=mean(min(dist01(tempseq2,:),[],1))/g0;
            LCC=length(tempseq2);
            SCratio=LCC/LC;
            F(pp)=purity+SCratio*gamma;
            if pp>=3
                [kpt]=kneept(F);
                if kpt~=pp
                    tempseq2=tempseq4;
                    tempdist0=tempdist1;
                    LCC=length(tempseq2);
                    centre0=zeros(LCC,W);
                    outp0=zeros(LCC,size(y,2));
                    member0=zeros(LCC,1);
                    for ii=1:1:LCC
                        tempseq3=find(tempdist0(tempseq2(ii),:)==1);
                        centre0(ii,:)=mean(data(tempseq3,:),1);
                        outp0(ii,:)=mean(y(tempseq3,:),1);
                        member0(ii)=length(tempseq3);
                    end
                    break
                end
            end
            tempdist1=tempdist0;
            tempseq4=tempseq2;
            pp=pp+1;
        end
           centre=centre0;
            outp=outp0;
            member=member0;
        averdist=(Npd*averdist+LC*averdist1)/(Npd+LC);
        Npd=Npd+LC;
        if tt==1
            CEN=centre;
            OUT=outp;
            MEM=member;
        else
            dist4=pdist2(centre,CEN,disttype,P0).^P1;
            [distseq,idx]=min(dist4,[],2);
            distseq=distseq.^P1;
            seq2=find(distseq>averdist);
            seq3=1:1:length(member(:,1));
            seq3(seq2)=[];
            CEN=[CEN;centre(seq2,:)];
            OUT=[OUT;outp(seq2,:)];
            MEM=[MEM;member(seq2)];
            for jj=seq3
                temp=MEM(idx(jj))+member(jj);
                CEN(idx(jj),:)=((CEN(idx(jj),:)*MEM(idx(jj))+centre(jj,:)*member(jj))/temp);
                OUT(idx(jj),:)=((OUT(idx(jj),:)*MEM(idx(jj))+outp(jj,:)*member(jj))/temp);
                MEM(idx(jj),:)=temp;
            end
        end
    end
    output.CEN=CEN;
    output.OUT=OUT;
    output.MEM=MEM;
    output.AVD=averdist;
end
if strcmp(mode,'testing')==1
    data=input.data;
    ck0=input.chunksize;
    averdist=input.AVD;
    CEN=input.CEN;
    OUT=input.OUT;
    MEM=(input.MEM)';
    [L,W]=size(data);
    seqck=[1:ck0:L,L+1];
    Lck=length(seqck);
    CL=length(CEN(:,1));
    for tt=1:1:(Lck-1)
        seq1=seqck(tt):1:seqck(tt+1)-1;
        dist30=pdist2(data(seq1,:),CEN,disttype,P0).^P1;
        dist3=exp(-1*(dist30./averdist));
        seq=find(sum(dist3,2)==0);
        dist3(seq,:)=repmat(ones(1,CL),length(seq),1);
%          dist3=dist3.*repmat(MEM,length(seq1),1);
        dist3=dist3./repmat(sum(dist3,2),1,CL);
        PredY1(seq1,:)=dist3*OUT;
    end
    if strcmp(task,'c')==1
        [~,output.pred] = max(PredY1,[],2);
    elseif strcmp(task,'r')==1
        output.pred=PredY1;
    end
end
end
function [kpt]=kneept(y)
[~,kpt]=min(y);
end