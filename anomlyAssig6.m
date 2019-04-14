% 1-calc mean of data (col)
% 2-calc std of data 
% 3-calc pdf by using Normpdf(x,mean,std)- each one of the inputs should be
% 18(do it in a loop) and put the o/p in a vector
% for i=1:18
  % a(i) = Normpdf(x,mean,std)
  %end
  
% 4- deh betdrb prod(a) > threshold (0.999), multiply all the pdfs to know
% wether the data is anomly or not
% if there is correlation bet the features then the anomly will not be accurate 

clc
clear all
close all

[num,text] = xlsread('house_prices_data_training_data.csv');
x = num(:,4:end);
n = 18;
m = mean(x);
s = std(x);
a = [];
counter = 0;
[r, c]= size(x);
for j=1:r
for i = 1:c
a(i) = normcdf(x(j,i),m(i),s(i));
end
bol = prod(a);
if bol < 10^-7 || bol > 0.999
   
   counter = counter +1;

  
end
end
counter