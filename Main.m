clc,clear
close all
%��������
load('data1.mat')
% load('data2.mat')

Dist=zeros(size(location,1),size(location,1));
for i=1:size(location,1)
    for j=i+1:size(location,1)
        Dist(i,j)=((location(i,1)-location(j,1))^2+(location(i,2)-location(j,2))^2)^0.5;
        Dist(j,i)=Dist(i,j);
    end
end
%�㷨ͨ�ò���
SizePop=100;             %��Ⱥ��ģ
MaxGen=500;            	%������������
StopGen=50;             %���Ž����ٱ��ִ���
%%
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('��������ò�ͬ���˹������㷨���TSP����')
disp('�˹������㷨ѡ��')
fprintf(' 1.�Ŵ��㷨    2.�������Ⱥ�㷨    3.ģ���˻��㷨    4.��Ⱥ�㷨 \n')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
number_ANN=input('�������˹������㷨���=');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if number_ANN==1            %�Ŵ��㷨

    %�㷨����
    GGAP=0.9;               %ѡ����ʣ���ѡ������Ĭ��Ϊ1��
    Pc=0.9;                 %�������
    Pm=0.05;                %�������
    %�㷨����
    for i=1:SizePop
        Chrom(i,:)=randperm(size(Dist,1)); %��ʼ����Ⱥ
        ObjV(i,1)=ObjF(Chrom(i,:),Dist); 	%�����ʼ��Ⱥ�����Ŀ�꺯��ֵ
    end
    [bestObjV,bestindex]=min(ObjV);
    bestsol=Chrom(bestindex,:);
    %������ʼ
    gen=1;                              %��ʼ�Ŵ�����
    gen0=0;                             %��¼���ֲ���Ĵ���   
    while gen0<StopGen && gen<=MaxGen
        disp(gen)
        FitnV=ranking(ObjV);          	%������Ӧ��ֵ(Assign fitness values)
        %ranking������ԽС�����Խ��
        SelCh=TSP_Select(Chrom,FitnV,GGAP);%ѡ��
        SelCh=TSP_Recombin(SelCh,Pc);   %����
        SelCh=TSP_Mutate(SelCh,Pm);     %����
        SelCh=TSP_Reverse(SelCh,Dist);     %��ת
        Chrom=TSP_Reins(Chrom,SelCh,ObjV);%�ز����Ӵ�������Ⱥ
         for i=1:SizePop
            ObjV(i,1)=ObjF(Chrom(i,:),Dist);%�����ʼ��Ⱥ�����Ŀ�꺯��ֵ
        end
        [newbestObjV,newbestindex]=min(ObjV);
        [worestObjV,worestindex]=max(ObjV);
        if newbestObjV<bestObjV
            bestObjV=newbestObjV;
            bestsol=Chrom(newbestindex,:);
            gen0=0;
        else
            gen0=gen0+1;                %����ֵ���ִ�����1
        end
        Chrom(worestindex,:)=bestsol;
        ObjV(worestindex,:)=bestObjV;
        trace(gen,1)=bestObjV;         	%��¼ÿ��������ֵ
        trace(gen,2)=sum(ObjV)/SizePop; %��¼ÿһ��������ƽ����Ӧ��
        gen=gen+1;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
elseif number_ANN==2      %�������Ⱥ�㷨
    
    num=size(location,1);                 %����λ������
    for i=1:SizePop
        individual(i,:)=randperm(num);    %���ɳ�ʼ��
    end
    individual_ObjV=ones(SizePop,1)*Inf;%��¼��������ֵ
    bestObjV=Inf;                       %��¼Ⱥ������ֵ
    gen=1;                              %��ʼ�Ŵ�����
    gen0=0;                             %��¼���ֲ���Ĵ���   
    while gen0<StopGen && gen<=MaxGen
        disp(gen)
        newindividual=individual;
        %������Ӧ��ֵ
        for i=1:SizePop
            ObjV(i,:)=ObjF(individual(i,:),Dist);
        end
        %���µ�ǰ���ź���ʷ����
        for i=1:SizePop
            if ObjV(i)<individual_ObjV(i)
                individual_ObjV(i)=ObjV(i);
                individual_sol(i,:)=individual(i,:);
            end
        end
        [newbestObjV,index]=min(individual_ObjV);
        if newbestObjV<bestObjV
            bestObjV=newbestObjV;
            gen0=0;
            bestsol=individual_sol(index,:);
        else 
            gen0=gen0+1;
        end
        trace(gen,1)=bestObjV;
        trace(gen,2)=mean(individual_ObjV);
        %����ͱ������
        for i=1:SizePop
            %��������Ž��н���
            temp1=randperm(num);
            c1=temp1(1);
            c2=temp1(2);
            cros=individual_sol(i,min(c1,c2):max(c1,c2));
            %ɾ���뽻��������ͬԪ��
            temp2=newindividual(i,:);
            for j=1:size(cros,2)
                temp2(temp2==cros(j))=[];
            end
            newindividual(i,:)=[temp2 cros];
            %��·�����ȱ�������
            newObjV=ObjF(newindividual(i,:),Dist);
            if newObjV<ObjV(i)
                individual(i,:)=newindividual(i,:);
                ObjV(i)=newObjV;
            end
            % ��ȫ�����Ž��н���
            temp3=randperm(num);
            c1=temp3(1);
            c2=temp3(2);
            cros=bestsol(min(c1,c2):max(c1,c2));
            %ɾ���뽻��������ͬԪ��
            temp4=newindividual(i,:);
            for j=1:size(cros,2)
                temp4(temp4==cros(j))=[];
            end
            newindividual(i,:)=[temp4 cros];
            %��·�����ȱ�������
            newObjV=ObjF(newindividual(i,:),Dist);
            if newObjV<ObjV(i)
                individual(i,:)=newindividual(i,:);
                ObjV(i)=newObjV;
            end
            %�������
            temp5=randperm(num);
            c1=temp5(1);
            c2=temp5(2); 
            temp6=newindividual(i,c1);
            newindividual(i,c1)=newindividual(i,c2);
            newindividual(i,c2)=temp6;
            %��·�����ȱ�������
            newObjV=ObjF(newindividual(i,:),Dist);
            if newObjV<ObjV(i)
                individual(i,:)=newindividual(i,:);
                ObjV(i)=newObjV;
            end
        end
        gen=gen+1;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
elseif number_ANN==3      %ģ���˻��㷨
    
    %ģ���˻��㷨����
    T0=1000;   %��ʼ�¶�
    Tend=1e-6;  %��ֹ�¶�
    q=0.97;    %��������
    num=size(location,1);    %���еĸ���
    %���ɳ�ʼ��
    S1=randperm(num);
    bestObjV=ObjF(S1,Dist);                    	%��¼����ֵ
    bestsol=S1;                                 %��¼���Ž�
    %��ʼ����
    gen=1;                                      %��¼��������
    gen2=0;                                     %��¼�¶��½�������Ŀ�꺯��ֵ����Ĵ���  
    while T0>Tend && gen2<StopGen*2
        disp(gen)                               %�����ǰ��������
        temp=[];
        gen3=1;
        S2=TSP_NewAnswer(S1);
        [S1,R]=TSP_Metropolis(S1,S2,Dist,T0);  %Metropolis �����㷨
        temp(gen3,:)=[S1 R];
        gen4=0;                                 %��¼�̶��¶���Ŀ�꺯��ֵ����Ĵ���
        while  gen3<=MaxGen && gen4<StopGen*5
            gen3=gen3+1;
            S2=TSP_NewAnswer(S1);%�����½�
            [S1,R]=TSP_Metropolis(S1,S2,Dist,T0);%Metropolis�����ж��Ƿ�����½�
            temp(gen3,:)=[S1 R];          %��¼��һ·�ߵļ���·��
            if temp(gen3,end)<temp(gen3-1,end)
                gen4=0;
            else
                gen4=gen4+1;
            end
        end
        %��¼ÿ�ε������̵�����·��
        [newbestObjV,index]=min(temp(:,end)); 	%�ҳ���ǰ�¶�������ֵ
        if newbestObjV<bestObjV
            bestObjV=newbestObjV;             	%�����ǰ�¶�������ֵС����һ�¶�������ֵ�����¼��ǰ����ֵ
            bestsol=temp(index,1:end-1);
            gen2=0;
        else
            gen2=gen2+1;
        end
        trace(gen,1)=bestObjV;                	%��¼��ǰ�¶ȵ�����ֵ
        T0=q*T0;                               	%����
        gen=gen+1;                              %���µ�������
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
elseif number_ANN==4      %��Ⱥ�㷨
    
    %��ʼ������
    alpha = 1;                                  %��Ϣ����Ҫ�̶�����
    beta = 5;                                  	%����������Ҫ�̶�����
    rho = 0.1;                                 	%��Ϣ�ػӷ�����
    Q = 1;                                     	%��ϵ��
    Eta = 1./Dist;                             	%��������
    Tau = ones(size(location,1));             	%��Ϣ�ؾ���
    %�㷨����
	gen=1;                            % ����������ֵ
    disp(gen)
    Table=zeros(SizePop,size(location,1));       	% ·����¼��
	for i=1:SizePop
        Table(i,1)=randperm(size(location,1),1);
	end
	% ������ռ�
	location_index=1:size(location,1);
	% �������·��ѡ��
	for i=1:SizePop
        % �������·��ѡ��
        for j=2:size(location,1)
            tabu=Table(i,1:(j - 1));           % �ѷ��ʵĳ��м���(���ɱ�)
            allow=location_index(~ismember(location_index,tabu));  % �����ʵĳ��м���
            P=ones(1,length(allow));
            % ������м�ת�Ƹ���
            for k=1:length(allow)
                P(k)=Tau(tabu(end),allow(k))^alpha * Eta(tabu(end),allow(k))^beta;%ת�Ƹ��ʵķ���
            end
            P=P/sum(P);%ת�Ƹ���
            % ���̶ķ�ѡ����һ�����ʳ���
            Pc=cumsum(P);     
            Table(i,j)=allow(find(Pc>=rand,1));
        end
	end
	% ����������ϵ�·������
	ObjV=zeros(SizePop,1);
	for i=1:SizePop
        ObjV(i)=ObjF(Table(i,:),Dist);
	end
	% ������Ϣ��
	Delta_Tau=zeros(size(location,1));
	% ������ϼ���
	for i=1:SizePop
        % ������м���
        for j=1:size(location,1)-1
            Delta_Tau(Table(i,j),Table(i,j+1))=Delta_Tau(Table(i,j),Table(i,j+1))+Q/ObjV(i);
        end
        Delta_Tau(Table(i,size(location,1)),Table(i,1))=Delta_Tau(Table(i,size(location,1)),Table(i,1))+Q/ObjV(i);
	end
	Tau=(1-rho)*Tau+Delta_Tau;
    [bestObjV,bestindex]=min(ObjV);
    trace(gen,1)=bestObjV;  
    trace(gen,2)=mean(ObjV);
    bestsol=Table(bestindex,:);
    % ������ʼ
    gen0=0;                             %��¼���ֲ���Ĵ���   
    while gen0<StopGen && gen<=MaxGen
        gen=gen+1;
        disp(gen)
        Table=zeros(SizePop,size(location,1));
        % ��������������ϵ�������
        for i=1:SizePop
            Table(i,1)=randperm(size(location,1),1);
        end
        % �������·��ѡ��
        for i=1:SizePop
            % �������·��ѡ��
            for j=2:size(location,1)
                tabu=Table(i,1:(j - 1));           % �ѷ��ʵĳ��м���(���ɱ�)
                allow=location_index(~ismember(location_index,tabu));  % �����ʵĳ��м���
                P=ones(1,length(allow));
                % ������м�ת�Ƹ���
                for k=1:length(allow)
                    P(k)=Tau(tabu(end),allow(k))^alpha * Eta(tabu(end),allow(k))^beta;%ת�Ƹ��ʵķ���
                end
                P=P/sum(P);%ת�Ƹ���
                % ���̶ķ�ѡ����һ�����ʳ���
                Pc=cumsum(P);     
                Table(i,j)=allow(find(Pc>=rand,1));
            end
        end
        % ����������ϵ�·������
        ObjV=zeros(SizePop,1);
        for i=1:SizePop
            ObjV(i)=ObjF(Table(i,:),Dist);
        end
        % ������Ϣ��
        Delta_Tau=zeros(size(location,1));
        % ������ϼ���
        for i=1:SizePop
            % ������м���
            for j=1:size(location,1)-1
                Delta_Tau(Table(i,j),Table(i,j+1))=Delta_Tau(Table(i,j),Table(i,j+1))+Q/ObjV(i);
            end
            Delta_Tau(Table(i,size(location,1)),Table(i,1))=Delta_Tau(Table(i,size(location,1)),Table(i,1))+Q/ObjV(i);
        end
        Tau=(1-rho)*Tau+Delta_Tau;
        [newbestObjV,newbestindex]=min(ObjV);
        if newbestObjV<bestObjV
            bestObjV=newbestObjV;
            bestsol=Table(newbestindex,:);
            gen0=0;
        else
            gen0=gen0+1;                %����ֵ���ִ�����1
        end
        trace(gen,1)=bestObjV;         	%��¼ÿ��������ֵ
        trace(gen,2)=sum(ObjV)/SizePop; %��¼ÿһ��������ƽ����Ӧ��
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    disp('�������˹������㷨��ų�����Χ������')
end
%%
%�����ͼ
%ͼ1���Ż�����ͼ
if size(trace,2)==1
    plot(1:size(trace,1),trace(:,1));
    legend('��������ֵ');
elseif size(trace,2)==2
    plot(1:size(trace,1),trace(:,1),'r-',1:size(trace,1),trace(:,2),'b--');
    legend('��������ֵ','����ƽ��ֵ');
end
title(['Ŀ�꺯��ֵ�Ż�����  ' '��ֹ������' num2str(size(trace,1))],'fontsize',12);
xlabel('��������','fontsize',12);
ylabel('Ŀ�꺯��ֵ','fontsize',12);
xlim([1,size(trace,1)])
%ͼ2��·�߹滮ͼ
figure;
hold on
plot(location(:,1),location(:,2),'r.', 'MarkerSize', 20)
for i=1:size(location,1)
    text(location(i,1)+0.1,location(i,2)+0.1,num2str(i),'color',[1,0,0]);
end
plot(location(bestsol(1),1),location(bestsol(1),2),'rp','MarkerSize',20)
R=[bestsol bestsol(1)]; %һ�������(����)
A=location(R,:);
for i=2:size(A,1)
    [arrowx,arrowy] = dsxy2figxy(gca,A(i-1:i,1),A(i-1:i,2));%����ת��
    annotation('textarrow',arrowx,arrowy,'HeadWidth',8,'color',[0,0,1]);
end
hold off
xlabel('������')
ylabel('������')
title(['�Ż�·��(����ܾ���:' num2str(bestObjV) ')'])
box on
grid on
%% 
%������
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('TSP����·�ߣ�')
R=[bestsol,bestsol(1)];
p=num2str(R(1));
for i=2:length(R)
    p=[p,'->',num2str(R(i))];
end
disp(p)
disp(['����ܾ��룺',num2str(bestObjV)]);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')