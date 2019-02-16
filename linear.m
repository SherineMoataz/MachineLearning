clc
clear all 
close all 
Filename='house_prices_data_training_data.csv';
Variablest=xlsread(Filename);
Price_all=Variablest(:,3);
m_all = length(Price_all);
alpha=0.003;
MSE=[];
lamda=0.02;
%removing the extra variables from the list:
n = length(Variablest(1,:));
Variabless =Variablest(:,[4:n]);

% for i=1:n
% figure(i)
% stem(Price_all,Variables(:,i))
% end
%generation of diff hyp:
%x1= [Variabless];%best hypothesis
x1= [Variabless Variabless(:,[3,9,10,16]).^3];
%x1= [Variabless Variabless(:,[9,16]).^2 Variabless(:,[3,10]).^3];
%x1= [Variabless Variabless(:,[3,9]).^2]; 
x0= ones(m_all,1);
%generation of x0;
x1 = [x0 x1];
%x1= [x0 x1];
%x1= [x0 x1];
%x1=[x0 x1];

% %generation of thetas
% theta_old= ones(1,n);
%Cross_validation training set:
xtraining = x1([1:(0.6*m_all)],:);
Price= Price_all([1:(0.6*m_all)],:);
%generation of thetas
n= length(xtraining(1,:));
theta_old= ones(1,n);
m=length(xtraining(:,1));
%Scaling:
avg_pr=mean(Price);
st_pr=std(Price);
Price=(Price-avg_pr)/st_pr;
avg_training=zeros(1,n);
st_training=ones(1,n);
for i=2:n
avg_training(i)= mean(xtraining(:,i));
st_training(i) = std(xtraining(:,i));
i
temp=(xtraining(:,i)-avg_training(i))/st_training(i);
xtraining(:,i)= temp;

end
%generation of h and cost function
h=[];
for i=1:m
h(i) = theta_old*xtraining(i,:)';
end
J_old = ((1/(2*m))*sum((h' -Price).^2));
%generation of new thetas:
for i=1:n
theta_new= theta_old - (alpha*(1/m)*sum( (h' - Price).*xtraining(:,i)));
end
for i=1:m
h(i) = theta_new*xtraining(i,:)';
end
J_training = ((1/(2*m))*sum((h' -Price).^2));
iterations =1;
 mse1 = J_training ;
 hu=[];
 mse=[];
while J_old - J_training > 10^-5

theta_old = theta_new;
%generation of h and cost function

for i=1:m
h(i) = theta_old*xtraining(i,:)';
end
J_old = ((1/(2*m))*sum((h' -Price).^2));
for i=1:n
theta_new= theta_old - (alpha*(1/m)*sum( (h' - Price).*xtraining(:,i)));
end
for i=1:m
hu(i) = theta_new*xtraining(i,:)';
end
J_training = ((1/(2*m))*sum((hu' -Price).^2));
mse(iterations) = J_training;
iterations= iterations+1
end
mse = [mse1 mse];
iterations =length(mse);
figure(1)
plot([1:iterations],mse)
xlabel('No. of iterations')
ylabel('MSE')
% Cross Validation set:
p=(0.6*m_all);
d=(0.8*m_all);
xcross = x1(p:d,:);
Price_cross= Price_all(p:d,:);
n_cross= length(xcross(1,:));
m_cross=length(xcross(:,1));
h_cross=[];
%scaling:
Price_cross=(Price_cross - avg_pr)/st_pr;
for i=2:n_cross
temp=(xcross(:,i)-avg_training(i))/st_training(i);
xcross(:,i)= temp;
end
for i=1:m_cross
h_cross(i) = theta_new*xcross(i,:)';
end
J_cross = (1/(2*m_cross))*sum((h_cross' -Price_cross).^2);
%Test sets:
r=m_all;
d=(0.8*m_all);
xtest = x1(d:r,:);
Price_test= Price_all(d:r,:);
n_test= length(xtest(1,:));
m_test=length(xtest(:,1));
% Scaling:
Price_test=(Price_test-avg_pr)/st_pr;
for i=2:n_test
temp=(xtest(:,i)-avg_training(i))/st_training(i);
xtest(:,i)= temp;
end
%hypothesis
h_test=[];
for i=1:m_test
h_test(i) = theta_new*xtest(i,:)';
end
J_test = (1/(2*m_test))*sum((h_test' -Price_test).^2);

%Normal method
Price_normal=(Variablest(:,3)-mean(Variablest(:,3)))/std(Variablest(:,3));
 m_normal = length(Price_normal);
n_normal = length(Variablest(1,:));
for i=1:n_normal
temp=(Variablest(:,i)-mean(Variablest(:,i)))/std(Variablest(:,i));
Variablest(:,i)= temp;
end
%removing the unwanted variables 
Variablesw =Variablest(:,[4:n_normal]);
%x0
Variables = [ones(m_normal,1) Variablesw];
%thetas
theta_normal=inv(Variables'*Variables)*(Variables')*Price_normal;
for i=1:m_normal
h_normal(i) = theta_normal'*Variables(i,:)';
end
 J_normal = (1/(2*m_normal))*sum((h_normal' -Price_normal).^2);
