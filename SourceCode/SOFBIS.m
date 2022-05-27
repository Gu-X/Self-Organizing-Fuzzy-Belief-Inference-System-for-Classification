function [output]=SOFBIS(input,mode,task)
disttype='minkowski';
P0=1;P1=2;
if strcmp(mode,'learning')==1
    data0=input.data;
    y0=input.y;
    ck0=input.chunksize;
    granlevel=input.granlevel;
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
        [data,temp] = unique(data,'rows');
        y=y0(seqck(tt):1:seqck(tt+1)-1,:);
        if strcmp(task,'c')==1
            y=full(ind2vec(y(temp,:)',lc))';
        elseif strcmp(task,'r')==1
            y=y(temp,:);
        end
        LC=length(y(:,1));
        dist00=pdist(data,disttype,P0).^P1;
        dist01=squareform(dist00);
        tempB=mean(dist00);
        for qq=1:granlevel
            tempA=tempB;
            dist00(dist00>mean(dist00))=[];
            tempB=mean(dist00);
            if tempB==0
                granlevel=qq-1;
                tempB=tempA;
                break
            end
        end
        tempdist=zeros(LC);
        tempdist(dist01<=tempB)=1;
        tempdist0=tempdist;
        tempseq=sum(exp(-1*tempdist.*dist01./tempB),2);
        tempdist=tempdist.*repmat(tempseq',LC,1);
        tempseq=max(tempdist,[],2);
        tempseq1=diag(tempdist);
        tempseq2=find(tempseq1-tempseq==0);
        LCC=length(tempseq2);
        centre=zeros(LCC,W);
        outp=zeros(LCC,size(y,2));
        member=zeros(LCC,1);
        for ii=1:1:LCC
            tempseq3=find(tempdist0(tempseq2(ii),:)==1);
            centre(ii,:)=mean(data(tempseq3,:),1);
            outp(ii,:)=mean(y(tempseq3,:),1);
            member(ii)=length(tempseq3);
        end
        averdist=(Npd*averdist+LC*tempB)/(Npd+LC);
        Npd=Npd+LC;
        if tt==1
            CEN=centre;
            OUT=outp;
            MEM=member;
        else
            dist4=pdist2(centre,CEN,disttype,P0);
            [distseq,idx]=min(dist4,[],2);
            distseq=distseq.^P1;
            seq2=find(distseq>averdist);
            seq3=1:1:length(outp(:,1));
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
        dist3=dist3./repmat(sum(dist3,2),1,CL);
        PredY1(seq1,:)=dist3*OUT;
    end
    if strcmp(task,'c')==1
        [~,output.pred] = max(PredY1,[],2);
    elseif strcmp(task,'r')==1
        output.pred=PredY1;
    end
end
