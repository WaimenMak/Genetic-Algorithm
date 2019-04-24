function GeneticA_tsp(popsize,pc,pm) %实数编码
if nargin<1
    popsize = 30; % 种群大小
end
if nargin<2
    pc = 0.1;%杂交概率
end
if nargin<3
    pm = 0.7;%变异概率
end

% pos = [0,0;1,1;0,1;1,0];
% pos=[0.4 0.4439;0.2439 0.1463;0.1707 0.2293;0.2293 0.761;0.5171 0.9414;
%         0.8732 0.6536;0.6878 0.5219;0.8488 0.3609;0.6683 0.2536;0.6195 0.2634];%10 cities d'=2.691
%    pos = [41 94;37 84;54 67;25 62;7 64;2 99;68 58;71 44;54 62;83 69;64 60;18 54;22 60;
%         83 46;91 38;25 38;24 42;58 69;71 71;74 78;87 76;18 40;13 40;82 7;62 32;58 35;45 21;41 26;44 35;4 50];%30 cities d'=423.741 by D B Fogel
    pos=[48 21;52 26;55 50;50 50;41 46;51 42;55 45;38 33;33 34;45 35;40 37;50 30;
        55 34;54 38;26 13;15 5;21 48;29 39;33 44;15 19;16 19;12 17;50 40;22 53;21 36;
        20 30;26 29;40 20;36 26;62 48;67 41;62 35;65 27;62 24;55 20;35 51;30 50;
        45 42;21 45;36 6;6 25;11 28;26 59;30 60;22 22;27 24;30 20;35 16;54 10;50 15;
        44 13;35 60;40 60;40 66;31 76;47 66;50 70;57 72;55 65;2 38;7 43;9 56;15 56;
        10 70;17 64;55 57;62 57;70 64;64 4;59 5;50 4;60 15;66 14;66 8;43 26];%75 cities d'=549.18 by D B Fogel
%     
    
% city_num = 10;
city_num = 75;
% city_num = 30;


cost = 1e5;   %计算成本

X_pop = zeros(popsize+1,city_num);     

Distance = zeros(city_num,city_num);

fit_val = zeros(popsize+1,1);
for i = 1:city_num
    for j = 1:city_num
        Distance(i,j) = ((pos(i,1) - pos(j,1))^2 + (pos(i,2) - pos(j,2))^2)^0.5;
    end
end

T = ceil(cost/popsize);

average_len = zeros(1,T);

best_len = zeros(1,T);



for generation = 1:T
    if generation == 1
%-----------------------------------
% %随机打乱
%         %initialize pop
%         for i = 1:popsize+1
%             X_pop(i,:) = randperm(city_num) ;
%         end
%-----------------------------------
%%K均值
    K = 5;
%     K = 4;
%     K = 3;
    class_num = zeros(1,K+1);
    interval = zeros(1,K);
    class_num(1,1) = 0;
        Idx=kmeans(pos,K);
        k = 1;
        for j = 1:K
            for i2 = 1:city_num
                if Idx(i2) == j
                    X_pop(1,k) = i2;
                    k = k+1;
                end
            end
            class_num(j+1) = k;
        end
        class_num(1,end) = city_num;
        
        for i3 = 1:K
            interval(i3) = class_num(i3+1) - class_num(i3);
        end
%         X_pop(2:10,:) = InitRandPop(X_pop(1,:),interval,9);
            
%         for q = 2:popsize+1                  %把类内顺序打乱
%              for p = 1:size(class_num,2)-1
%                 temp = X_pop(1,class_num(p)+1:class_num(p+1));
%                 temp_ = randperm(size(temp,2));
%                 X_pop(q,class_num(p)+1:class_num(p+1)) = temp(temp_);
%              end
%         end

        for q = 2:popsize+1            %类内顺序不乱，打乱类间
            X_pop(q,:) = X_pop(1,:);
            if rand > 0.5
                temp_2 = X_pop(q,:);
                X_pop(q,:) = InitRandPop_copy(temp_2,interval,1);
            end
        end
%    
        
        
%         for q = 2:popsize+1
% %             X_pop(q,:) = randperm(city_num);    %后面全部打乱顺序
%             X_pop(q,:) = X_pop(1,:);              %把后面全部赋值
%    
%                 mpoint_1 = ceil(rand*(city_num));
%                 mpoint_2 = ceil(rand*(city_num));
%                 temp = X_pop(q,mpoint_1);
%                 X_pop(q,mpoint_1) = X_pop(q,mpoint_2);
%                 X_pop(q,mpoint_2) = temp;
%         end        

    else
        father_pop = fathers(fit_val,popsize+1,X_pop);
        new_pop = cross_over(father_pop,pc,pm,popsize,city_num,generation,T,rol(1),class_num);
        X_pop(2:popsize+1,:) = new_pop;    
    end
      
    for i = 1:popsize+1
        [fit_val(i)] = object_val(X_pop(i,:),Distance);
    end
    optimal = min(fit_val);
    average = mean(fit_val);
    [rol,~] = find(fit_val == optimal);
    if fit_val(1) > optimal
        X_pop(1,:) = X_pop(rol(1),:);
    end
    
    average_len(1,generation) = average;
    best_len(1,generation) = optimal;  %绘图用
    if generation == 1
        variable = [optimal, X_pop(rol(1),:)];
    elseif optimal < variable(1)
        variable = [optimal, X_pop(rol(1),:)];
    end
    if mod(generation,10)== 1
    visualization(variable,pos,city_num);
    title(['迭代次数为n=' num2str(generation)]);
    pause(0.05);
    hold off;
    end
end

    variable
    fprintf('The optimal is --->>%5.4f\n',variable(1));
%     fprintf('The last optimal is --->>%5.4f\n',optimal);
    
    i = 1:T;
    figure(2);
    plot(i,best_len);
    hold on
    plot(i,average_len,'g');
    xlabel('Generation');
    ylabel('Optimal');
    title('Convergency');
    
end

function [f_value] = object_val(entity,Dis)

% X_ = entity;
f_value = 0;
for i = 1:size(entity,2) - 1
    f_value = f_value + Dis(entity(i),entity(i+1));
end
f_value = f_value + Dis(entity(1),entity(size(entity,2))); %与首部相连

end

%轮盘
function fatherpop = fathers(fit_val,popsize,pop)
% fatherpop = zeros(size(pop));
fatherpop = pop;
% min_f = min(fit_val);
max_f = max(fit_val);
adaption_value = cumsum((max_f*ones(popsize,1) - fit_val)/(popsize*max_f - sum(fit_val)));

point = sort(rand(popsize,1));      %取值0，1之间
enti_num = 1;
poin_num = 1;
while poin_num <= popsize && enti_num <= popsize
    if adaption_value(enti_num) >= point(poin_num)
        fatherpop(poin_num,:) = pop(enti_num,:);
        poin_num = poin_num + 1;
    else
        enti_num = enti_num + 1;
    end
end

end

function newpop = cross_over(fatherpop,pc,pm,popsize,dim,t,T,rol,class_num)

%随机打乱
newpop = fatherpop(2:popsize+1,:);

% %不打乱
% newpop = fatherpop;
% pos = round(rand*dim);

%杂交
for k = 1:popsize
    r = rand;  
    if r < pc
        prob_1 = ceil(rand*dim);
        prob_2 = dim+1 - prob_1;
        if prob_1 == prob_2   %城市奇数才有的情况
           newpop(k,:) = randperm(dim);
        elseif prob_1 < prob_2
            temp = [newpop(k,1:prob_1) newpop(k,prob_2:dim)];
            c = randperm(size(temp,2));
            temp_2 = temp(c);
            newpop(k,prob_2:dim) = temp_2(1:size(temp_2,2)/2);
            newpop(k,1:prob_1) = temp_2(size(temp_2,2)/2+1:end);
        else
            temp = [newpop(k,1:prob_2) newpop(k,prob_1:dim)];
            c = randperm(size(temp,2));
            temp_2 = temp(c);
            newpop(k,prob_1:dim) = temp_2(1:size(temp_2,2)/2);
            newpop(k,1:prob_2) = temp_2(size(temp_2,2)/2+1:end);
%             temp = newpop(k,prob_1:dim);
%             newpop(k,prob_1:dim) = newpop(k,1:prob_2);
%             newpop(k,1:prob_2) = temp;
        end
    end
           
end

%%变异

for q = 1:popsize
% ------------------------------
%全部打乱
%     if (1 - rand*(1-t/T)) <=pm
%         newpop(q,:) = randperm(dim);
%     end
%-------------------------------
%部分打乱，类内打乱，效果不好
%         if (1-rand*(t/T))<=pm && (t/T) > (1/3)
%             for p = 1:size(class_num,2)-1
%                 if (rand*(t/T))<=pm
%                 temp = fatherpop(rol,class_num(p)+1:class_num(p+1));
%                 temp_ = randperm(size(temp,2));
%                 newpop(q,class_num(p)+1:class_num(p+1)) = temp(temp_);
%                 end
%             end
%         end
%-------------------------------
%交换类的组合顺序，效果没怎么改进

% if (rand*(t/T))<=pm && (t/T) > (2/3)
%     part = 2;
%     temp = randperm(part);      
%     temp_ = newpop(q,:);       %效果好一点
% %     temp_ = fatherpop(rol,:);  %效果不好
%     for o = 1:part          %-1:dim - dim/part+1
%         newpop(q,(o-1)*dim/part+1:o*dim/part) = temp_(1,(temp(1,o)-1)*dim/part+1:temp(1,o)*dim/part);
%     end
% end


%-------------------------------
    if (1 - rand*(1-t/T))<=pm
%     if rand <= pm
        mpoint_1 = ceil(rand*(dim));
        mpoint_2 = ceil(rand*(dim));
        temp = newpop(q,mpoint_1);
        newpop(q,mpoint_1) = newpop(q,mpoint_2);
        newpop(q,mpoint_2) = temp;
    end
end
end

function visualization(variable,position,city_num)
figure(1);
plot(position(:,1),position(:,2),'b*');
hold on
temp = zeros(city_num,2);
for i = 1:city_num
    temp(i,:) = position(variable(i+1),:);
end
plot(temp(:,1),temp(:,2),'rs-','LineWidth',1);
plot([temp(1,1),temp(city_num,1)],[temp(1,2),temp(city_num,2)],'rs-','LineWidth',1);

end


%InitRandPop(数据向量,各类长度向量,种群数量)%%%打乱每一类的顺序
function newpop = InitRandPop_copy(data,len,popsize)
[~,px] = size(len);
a = 1;
for i = 2:px
    len(i) = len(i-1) + len(i);
end;

for i = px+1:-1:2
    len(i) = len(i-1);
end;
len(1) = 0;
newpop = zeros(popsize,len(px+1));
for i = 1:popsize
    R = randperm(px);
    theta = 1;
    for j = 1:px
        b = len(R(j)+1);
        a = len(R(j));
        newpop(i,theta:(theta + b - a - 1)) = data(a+1:b);
        theta = theta + b - a;
    end;
end;
end
