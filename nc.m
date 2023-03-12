%x=input('give the data  ');
%sf=input('give the sampling frequency   ');
%function [afr,aco,SP1,SP2,con_lim]=nc(x,sf);
function V=nc(x,sf,winlen,mi,ma); 
L=winlen;
lx=length(x);sp1=zeros(L,1);sp2=zeros(L,1);cs=zeros(L,1);co=zeros(L,1);fi=zeros(L,1);
m=0;
for i=1:L:lx-L+1
    m=m+1;
    % if(i+L-1)>lx;break;end
    x1=x(i:i+L-1,1);
    x2=x(i:i+L-1,2);
%     h=hanning(L);
%     Q=1/L*sum(h);
%     x1=(x1.*h)/sqrt(Q);
%     x2=(x2.*h)/sqrt(Q);
    x1=(x1-mean(x1))/std(x1);
    x2=(x2-mean(x2))/std(x2);
    f1=fft(x1)/length(x1);
    f2=fft(x2)/length(x2);
    cs1=f1.*conj(f2);
    FI=unwrap(angle(cs1));
    cs=cs+cs1;
    p1=f1.*conj(f1);
    p2=f2.*conj(f2);
    sp1=sp1+p1;
    sp2=sp2+p2;
    fi=fi+FI;
  %  co=co+(abs(cs1)./sqrt(p1)./sqrt(p2));
end
cs=cs/m;
sp1=sp1/m;
sp2=sp2/m;
fi=angle(cs);
%co=abs(cs)./sqrt(sp1)./sqrt(sp2);
%co=co/m;
co=(abs(cs).^2)./sp1./sp2;
fr=(0:L-1)/L*sf;
q1=find(fr>=0);
q2=find(fr>=250);
aco=co(q1(1):q2(1));
afr=fr(q1(1):q2(1));
cf=0.01^(1/(m-1));% sinignificance of zero coherence 
con_lim=(1-cf)*ones(length(afr),1);
SP1=sp1(q1(1):q2(1));
SP2=sp2(q1(1):q2(1));
fi1=fi(q1(1):q2(1));
%figure
%plot(afr,aco)% to calculate the 8 parameters it has 8 columns
%('freq,powersp1,powersp2,coherence,conf-limit,phase,errobar for phase estimate,
%crossspectrum')
V(:,1)=afr';
V(:,2)=SP1;
V(:,3)=SP2;
V(:,4)=aco;
V(:,5)=con_lim;
V(:,6)=fi1;
V(:,7)=1.96*(1/(2*m)*(1./aco-1)).^0.5;
V(:,8)=cs(q1(1):q2(1));
return