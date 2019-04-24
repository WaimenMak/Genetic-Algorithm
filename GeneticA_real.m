function GeneticA_real(popsize,pc,pm) %实数编码
if nargin<1
    popsize = 10; % 种群大小
end
if nargin<2
    pc = 0.6;%杂交概率
end
if nargin<3
    pm = 0.01;%变异概率
end

% bounds = [-5*ones(5,1) 10*ones(5,1)];
% bounds = [-3.0 12.1;4.1 5.8];
% bounds = [-3.0 12.1;0 10];
% bounds = [-1,1;-1,1];
% bounds = [-15*ones(10,1) 30*ones(10,1)];%ackley
% bounds = [0*ones(2,1) pi*ones(2,1)]; %mich
% bounds = [-5*ones(10,1) 10*ones(10,1)];%rosen
% bounds = [-4.1*ones(10,1) 6.4*ones(10,1)];%rast
% bounds = [-480*ones(5,1) 750*ones(5,1)];%griewank
% bounds = [-8*ones(4,1) 12.5*ones(4,1)];
% bounds = [-100*ones(2,1) 100*ones(2,1)];%easom
bounds = [0 1;0 1];
% bounds = [0 1200;0 1200;-0.55 0.55;-0.55 0.55];

[dim,~] = size(bounds);
cost = 1e4;   %计算成本
decpop = zeros(popsize,1);   % 化为10进制（因变量）
X_pop = zeros(popsize,dim);        %十进制种群（自变量）
variable = zeros(1,dim+1);     %保存的变量

T = ceil(cost/popsize);

best_y = zeros(1,T);


for generation = 1:T
    if generation == 1
        %initialize pop
        for i = 1:popsize
            prob = rand;
            X_pop(i,:) =  ((bounds(:,2) - bounds(:,1))*prob + bounds(:,1))';
        end
%         for i = 1:popsize
%             for j = 1:dim
%                 prob = rand;
%                 X_pop(i,j) =  ((bounds(j,2) - bounds(j,1))*prob + bounds(j,1));
%             end
%         end
    else
        father_pop = fathers(decpop,popsize,X_pop);
        new_pop = cross_over(father_pop,pc,pm,popsize,dim,bounds,generation,T);
        X_pop = new_pop;    
    end
      
    for i = 1:popsize
        [decpop(i)] = object_function(X_pop(i,:),generation);
    end
    optimal = max(decpop);
    [rol,~] = find(decpop == optimal);
    best_y(1,generation) = optimal;
    if generation == 1
        variable = [optimal, X_pop(rol(1),:)];
    elseif optimal > variable(1)
        variable = [optimal, X_pop(rol(1),:)];
    end
%     if mod(generation,10)== 1 && dim <=2 
%     visualization(X_pop(:,1),X_pop(:,2),decpop,bounds);
%     title(['迭代次数为n=' num2str(generation)]);
%     plot3(X_pop(rol(1),1),X_pop(rol(1),2),optimal,'o','LineWidth',2);
%     pause;
%     hold off;
%     end
end

    variable
    fprintf('The optimal is --->>%5.4f\n',variable(1));
%     fprintf('The last optimal is --->>%5.4f\n',optimal);
    
    i = 1:T;
    plot(i,best_y);
    xlabel('Generation');
    ylabel('Optimal');
    title('Convergency');
    
end

function [f_value] = object_function(entity,iter)

X_ = entity;
% f_value = -(X_(1)^2 + X_(2)^2);
% f_value = 21.5+X_(1)*sin(4*pi*X_(1))+X_(2)*sin(20*pi*X_(2));
% f_value = -1*griewank(X_);
% f_value = -1*easom(X_);
% f_value = -1*mich(X_);
% f_value = -1*ackley(X_);
% f_value = -1*rosen(X_);
% f_value = -1*rast(X_);
% f_value = colville(X_);
% f_value = -1*Fun1(X_,1);
f_value = Fun2(X_,5);

end

%轮盘
function fatherpop = fathers(decpop,popsize,pop)
% fatherpop = zeros(size(pop));
fatherpop = pop;
min_f = min(decpop);
if min_f < 0              %检查是否有负数
    decpop = decpop - min_f*ones(popsize,1)+ones(popsize,1); %%不能直接加最小值，因为当所有值相同时，计算概率的分母和为0
end

adaption_value = cumsum(decpop/sum(decpop));
% adaption_value = cumsum(decpop + ones(size(decpop))/(sum(decpop) + popsize));
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

function newpop = cross_over(fatherpop,pc,pm,popsize,dim,bounds,t,T)

%随机打乱
randrow = randperm(popsize);
newpop = fatherpop(randrow,:);

% %不打乱
% newpop = fatherpop;
% pos = round(rand*dim);

%杂交
for k = 1:2:popsize - 1
    r = rand;  
    if r < pc
        fa_1 = newpop(k,:);
        for p = 1:dim
            r_1 = rand;
            newpop(k,p) = r_1*newpop(k,p) + (1-r_1)*newpop(k+1,p);
            newpop(k+1,p) = r_1*newpop(k+1,p) + (1-r_1)*fa_1(1,p);
        end
    end
end
%%变异

C = rand;%变异函数参数
b = 2;%变异函数参数
for q = 1:popsize
    mpoint = round(rand*(dim));
    r = rand(1);
    if r < pm
        if mpoint == 0
            mpoint = 1;
        end
        newpop(q,mpoint) = bounds(mpoint,1) + (bounds(mpoint,2) - bounds(mpoint,1))*C*(1-t/T)^b;
    end
end
end

function visualization(dim_1,dim_2,decpop,bounds)

[x,y] = meshgrid(bounds(1,1):0.1:bounds(1,2),bounds(2,1):0.1:bounds(2,2));
z = -(x.^2 + y.^2);
% [row,col] = size(x);
% z = zeros(row,col);

% for i = 1:row
%     for j = 1:col
%         z(i,j) = branin([x(i,j),y(i,j)]);
%     end
% end
% z = 21.5+x.*sin(4*pi.*x)+y.*sin(20*pi.*y);

surfl(x,y,z);
% view(0,0);
hold on;
plot3(dim_1,dim_2,decpop,'r*');

end

%% 惩罚以后的目标函数（适应值函数）
function y = Fun1(x,iteration)
% fun1 = 2*x(1)*x(2);
fun1 = 3*x(1) + 0.000001*x(1)^3 + 2*x(2) + 0.000002/(3*x(2)^3);
%惩罚函数很难选取好
% flag = x(1)^2 + x(2)^2 - 1; %对约束的违反程度
flag1 = 1000*sin(-x(3)-0.25) + 1000*sin(-x(4) - 0.25) + 894.8 - x(1);
flag2 = 1000*sin(x(3)-0.25) + 1000*sin(x(3)-x(4) - 0.25) + 894.8 - x(2);
flag3 = 1000*sin(x(4)-0.25) + 1000*sin(x(4)-x(3) - 0.25) + 1294.8;
flag4 = x(4) - x(3) + 0.55;flag5 = x(3) - x(4) + 0.55;
c_1 = 0;c_2 = 0;c_3 = 0;c_4 = 0;c_5 = 0;
% pun = zeros(size(x,1),1); % 惩罚函数值
% for ii=1:size(x,1)
    if abs(flag1)>0 && abs(flag1)<1e-4
        c_1 = 20;
    elseif abs(flag1)>1e-4 && abs(flag1)<0.1
        c_1 = 50;
    else
        c_1 = 100;
    end
    if abs(flag2)>0 && abs(flag2)<1e-4
        c_2 = 20;
    elseif abs(flag2)>1e-4 && abs(flag2)<0.1
        c_2 = 50;
    else
        c_2 = 100;
    end
    if abs(flag3)>0 && abs(flag3)<1e-4
        c_3 = 20;
    elseif abs(flag3)>1e-4 && abs(flag3)<0.1
        c_3 = 50;
    else
        c_3 = 100;
    end
    if flag4<0 && flag4>-1e-4
        c_4 = 20;
    elseif flag4<-1e-4
        c_4 = 100;
    end
    if flag5<0 && flag5>-1e-4
        c_5 = 20;
    elseif flag5<-1e-4
        c_5 = 100;
    end
    pun = c_1*abs(flag1)^2+c_2*abs(flag2)^2+c_3*abs(flag3)^2+c_4*abs(flag4)^2+c_5*abs(flag5)^2;
% end
y = fun1 + pun;
end

function y = Fun2(x,iteration)
fun1 = 2*x(1)*x(2);

%惩罚函数很难选取好
flag = x(1)^2 + x(2)^2 - 1; %对约束的违反程度
pun = 0;
if abs(flag)>1e-4
    pun = iteration^2*flag^2;
end

% end
y = fun1 - pun;
end
