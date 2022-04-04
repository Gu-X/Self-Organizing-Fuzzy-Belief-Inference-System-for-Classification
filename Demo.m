clear all
clc
close all

load("data.mat")
color1=[[0, 0, 1];[1, 0, 0];[0, 0.4, 0]];
color=[[0.25, 0.6970, 0.9910];[1.00, 0.5750, 0.3480];[0.4660, 0.6740, 0.1880]];
color2=color*1.8;
color2(color2>1)=1;
color=color/1.2;
y=-1.95:0.01:3.45;
L=length(y);
data1=[];
for x=-2.45:0.01:2.45
    data1=[data1;[ones(L,1)*x,y']];
end


figure
plot(data0(Y0==1,1),data0(Y0==1,2),'.','markersize',16,'color',color(1,:));
hold on
plot(data0(Y0==2,1),data0(Y0==2,2),'.','markersize',16,'color',color(2,:));
plot(data0(Y0==3,1),data0(Y0==3,2),'.','markersize',16,'color',color(3,:));
axis([-2.5,2.5,-2,3.5])
xlabel('x_1')
ylabel('x_2')
% grid on
set(gca,'fontsize',16)
hold off

GranLevel=2;
mode='learning';

task='c';
L=length(Y0);
input.y=Y0;
input.data=data0;
input.chunksize=L;
input.granlevel=GranLevel;
[output0]=SOFBIS(input,mode,task);
center1=output0.CEN;

mode='testing';
input=output0;
input.chunksize=100;
input.data=data0;
[output]=SOFBIS(input,mode,task);
label =output.pred;

mode='testing';
input=output0;
input.chunksize=100;
input.data=data1;
[output]=SOFBIS(input,mode,task);
label2 =output.pred;

figure
hold on
scatter(data1(label2==1,1),data1(label2==1,2),16,'.','MarkerEdgeColor',color2(1,:),'MarkerEdgeAlpha',0.8);
scatter(data1(label2==2,1),data1(label2==2,2),16,'.','MarkerEdgeColor',color2(2,:),'MarkerEdgeAlpha',0.8);
scatter(data1(label2==3,1),data1(label2==3,2),16,'.','MarkerEdgeColor',color2(3,:),'MarkerEdgeAlpha',0.8);
plot(data0(Y0==1,1),data0(Y0==1,2),'.','markersize',16,'color',color(1,:));
plot(data0(Y0==2,1),data0(Y0==2,2),'.','markersize',16,'color',color(2,:));
plot(data0(Y0==3,1),data0(Y0==3,2),'.','markersize',16,'color',color(3,:));
xlabel('x_1')
ylabel('x_2')
grid on
voronoi(center1(:,1),center1(:,2),'k--');
set(gca,'fontsize',16)
plot(center1(:,1),center1(:,2),'.','markersize',24,'linewidth',2,'color',[0,0.0,0]);
axis([-2.5,2.5,-2,3.5])
hold off
set(gca, 'box', 'on')



%%
GranLevel=2;
mode='learning';
task='c';
L=length(Y0);
input.y=Y0;
input.data=data0;
input.chunksize=L;
input.gamma=1;
[output0]=SOFBISplus(input,mode,task);
center1=output0.CEN;

mode='testing';
input=output0;
input.chunksize=100;
input.data=data0;
[output]=SOFBISplus(input,mode,task);
label =output.pred;

mode='testing';
input=output0;
input.chunksize=100;
input.data=data1;
[output]=SOFBISplus(input,mode,task);
label2 =output.pred;

figure
hold on
scatter(data1(label2==1,1),data1(label2==1,2),16,'.','MarkerEdgeColor',color2(1,:),'MarkerEdgeAlpha',0.8);
scatter(data1(label2==2,1),data1(label2==2,2),16,'.','MarkerEdgeColor',color2(2,:),'MarkerEdgeAlpha',0.8);
scatter(data1(label2==3,1),data1(label2==3,2),16,'.','MarkerEdgeColor',color2(3,:),'MarkerEdgeAlpha',0.8);
plot(data0(Y0==1,1),data0(Y0==1,2),'.','markersize',16,'color',color(1,:));
plot(data0(Y0==2,1),data0(Y0==2,2),'.','markersize',16,'color',color(2,:));
plot(data0(Y0==3,1),data0(Y0==3,2),'.','markersize',16,'color',color(3,:));
xlabel('x_1')
ylabel('x_2')
grid on
voronoi(center1(:,1),center1(:,2),'k--');
set(gca,'fontsize',16)
plot(center1(:,1),center1(:,2),'.','markersize',24,'linewidth',2,'color',[0,0.0,0]);
axis([-2.5,2.5,-2,3.5])
hold off
set(gca, 'box', 'on')