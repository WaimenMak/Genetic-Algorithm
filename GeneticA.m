function GeneticA(popsize,pc,pm)
if nargin<1
    popsize = 200; % 种群大小
end
if nargin<2
    pc = 0.6;%杂交概率
end
if nargin<3
    pm = 0.01;%变异概率
end

% bound1 = 1;    %边界   %1是大的，2是小的
% bound2 = -1;
% bound3 = 1;
% bound4 = -1;
bound1 = 12.1;    %边界   %1是大的，2是小的
bound2 = -3.0;
bound3 = 5.8;
bound4 = 4.1;
bound = [bound1,bound2,bound3,bound4];

q = 4 ; %%保留小数两位;
chrom_length1 = ceil(log2((bound(1) - bound(2))*10^q));  %染色体长度
chrom_length2 = ceil(log2((bound(3) - bound(4))*10^q));
chrom_length = chrom_length1 + chrom_length2;

cost = 1e4;   %计算成本

% popsize = 20;  %种群个体数
pop = zeros(popsize,chrom_length);
decpop = zeros(popsize,1);   % 化为10进制
x1 = zeros(popsize,1);x2 = zeros(popsize,1);
variable = zeros(1,3);
for generation = 1:ceil(cost/popsize)
    if generation == 1
        for i = 1:popsize
            pop(i,:) = round(rand(1,chrom_length));
        end
    else
        father_pop = fathers(decpop,popsize,pop);
        new_pop = cross_over(father_pop,pc,pm,chrom_length,popsize);
        pop = new_pop;    
    end
      
    for i = 1:popsize
        [decpop(i),x1(i),x2(i)] = object_function(pop(i,:),chrom_length1,chrom_length2,bound);
    end
    optimal = max(decpop);
    [rol,~] = find(decpop == optimal);
    if generation == 1
        variable = [optimal, x1(rol(1)),x2(rol(1))];
    elseif optimal > variable(1)
        variable = [optimal, x1(rol(1)),x2(rol(1))];
    end
    if mod(generation,10) == 1
        visualization(x1,x2,decpop,bound);
        title(['迭代次数为n=' num2str(generation)]);
        plot3(x1(rol(1)),x2(rol(1)),optimal,'o','LineWidth',2);
        pause;
        hold off;
    end
end
fprintf('var1 --->>%5.2f\n',variable(2));
fprintf('var2 --->>%5.2f\n',variable(3));
fprintf('The optimal is --->>%5.2f\n',variable(1));
fprintf('The last optimal is --->>%5.4f\n',optimal);

end

function visualization(dim_1,dim_2,decpop,bound)

[x,y] = meshgrid(bound(2):0.1:bound(1),bound(4):0.1:bound(3));
% z = -(x.^2 + y.^2);
z = 21.5+x.*sin(4*pi.*x)+y.*sin(20*pi.*y);

surfl(x,y,z);
% view(0,0);
hold on;
plot3(dim_1,dim_2,decpop,'r*');

end

function [f_value,x_1,x_2] = object_function(entity,chrom_length1,chrom_length2,bound)
chrom_length = chrom_length1 + chrom_length2;
temp = entity(1:chrom_length1);
x_1 = bound(2) + mybin2dec(temp)*(bound(1) - bound(2))/(2^chrom_length1 - 1);
temp2 = entity(chrom_length1 + 1:chrom_length);
x_2 = bound(4) + mybin2dec(temp2)*(bound(3) - bound(4))/(2^chrom_length2 - 1);

% f_value = -(x_1^2 + x_2^2);
f_value = 21.5+x_1*sin(4*pi*x_1)+x_2*sin(20*pi*x_2);
end



function dec = mybin2dec(bin)
dec = 0;
for o = 1:size(bin,2)
    dec = dec + bin(o)*2^(size(bin,2) - o);
end
end



%轮盘
function fatherpop = fathers(decpop,popsize,pop)
fatherpop = zeros(size(pop));
min_f = min(decpop);
if min_f < 0              %检查是否有负数
    decpop = decpop - min_f*ones(popsize,1)+ones(popsize,1); %%不能直接加最小值，因为当所有值相同时，计算概率的分母和为0
end

adaption_value = cumsum(decpop/sum(decpop));
point = sort(rand(popsize,1));      %有问题
enti_num = 1;
poin_num = 1;
while poin_num <= popsize
    if adaption_value(enti_num) >= point(poin_num)
        fatherpop(poin_num,:) = pop(enti_num,:);
        poin_num = poin_num + 1;
    else
        enti_num = enti_num + 1;
    end
end
end

function newpop = cross_over(fatherpop,pc,pm,chromlength,popsize)

%随机打乱
randrow = randperm(popsize);
newpop = fatherpop(randrow,:);

% %不打乱
% newpop = fatherpop;
pos = round(rand*chromlength);
%杂交
for k = 1:2:popsize - 1
    r = rand(1);
    if r < pc
        temp = newpop(k,1:pos);
        newpop(k,1:pos) = newpop(k+1,1:pos);
        newpop(k+1,1:pos) = temp;
    end
end
%%变异
for i = 1:popsize
    mpoint = round(rand*chromlength);
    r = rand(1);
    if r < pm
        if mpoint == 0
            mpoint = 1;
        else
            if newpop(i,mpoint) == 1
                newpop(i,mpoint) = 0;
            else
                newpop(i,mpoint) = 1;
            end
        end
    end
end
end
        
        

        
        
        
        
    





            
    