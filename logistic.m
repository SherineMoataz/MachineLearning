clc
clear all 
close all 
Filename='heart_DD.csv';
Variablest=xlsread(Filename);
target_all=Variablest(:,14);
m_all = length(target_all);
%removing target from the list:
Variabless =Variablest(:,[1:13]);
n = length(Variabless(1,:));
alpha=0.003;
MSE=[];
lamda=0.02;
n = length(Variabless(1,:));
% for i=1:n
% figure(i)
% stem(Variables(:,i),target_all)
% end
%generation of diff hyp:
%x1= [Variabless(:,1:2:n)];
%x1= [Variabless Variabless(:,[3,9,10,13]).^2];
x1= [Variabless]; %best hypothesis
%x1= [Variabless Variabless(:,[10,13]).^2 Variabless(:,[4,6]).^3 ];
x0= ones(m_all,1);
%generation of x0;
x1 = [x0 x1];
%x1= [x0 x1];
%x1= [x0 x1];
%x1=[x0 x1];
n= length(x1(1,:));
% %generation of thetas
% theta_old= ones(1,n);
%Cross_validation training set:
xtraining = x1([1:(0.6*m_all)],:);
target= target_all([1:(0.6*m_all)],:);
%generation of thetas
n= length(xtraining(1,:));
theta_old= ones(1,n);
m=length(xtraining(:,1));

avg_tr=mean(target);
 st_tr=1;

 target=(target - avg_tr);
avg_training=zeros(1,n);
st_training=ones(1,n);
for i=2:n
avg_training(i)= mean(xtraining(:,i));
st_training(i) = std(xtraining(:,i));
temp=(xtraining(:,i)-avg_training(i))/st_training(i);
xtraining(:,i)= temp;

end

%generation of h and cost function
h=[];
phase=[];
for i=1:m
    i
phase(i)= theta_old*xtraining(i,:)';
h(i)= 1/(1+exp(-phase(i)));
end
J_old = ((1/(2*m))*sum((h' -target).^2));
%generation of new thetas:
for i=1:n
theta_new= theta_old - (alpha*(1/m)*sum( (h' - target).*xtraining(:,i)));
end
for i=1:m
phase = theta_new*xtraining(i,:)';
h(i)=1/(1+exp(-phase));
end
J_training = ((1/(2*m))*sum((h' -target).^2));
iterations =1;
 mse1 = J_training ;
 hu=[];
 mse=[];
while J_old - J_training > 10^-6

theta_old = theta_new;
%generation of h and cost function

for i=1:m
phase = theta_old*xtraining(i,:)';
h(i)=1/(1+exp(-phase));
end
J_old = ((1/(2*m))*sum((h' -target).^2));
for i=1:n
theta_new= theta_old - (alpha*(1/m)*sum( (h' - target).*xtraining(:,i)));
end
for i=1:m
phase = theta_new*xtraining(i,:)';
hu(i)=1/(1+exp(-phase));
end
J_training = ((1/(2*m))*sum((hu' -target).^2));
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
target_cross= target_all(p:d,:);
%generation of thetas
n_cross= length(xcross(1,:));
m_cross=length(xcross(:,1));
%scaling:
target_cross=(target_cross - avg_tr)/st_tr;
for i=2:n_cross
temp=(xcross(:,i)-avg_training(i))/st_training(i);
xcross(:,i)= temp;
end
h_cross=[];
for i=1:m_cross
phase = theta_new*xtraining(i,:)';
h_cross(i)=1/(1+exp(-phase));
end
J_cross = ((1/(2*m_cross))*sum((h_cross' -target_cross).^2));
 %Test sets:
r=m_all;
d=(0.8*m_all);
xtest = x1(d:r,:);
target_test= target_all(d:r,:);
n_test= length(xtest(1,:));
m_test=length(xtest(:,1));
%generation of thetas
% Scaling:
target_test=(target_test-avg_tr)/st_tr;
for i=2:n_test
temp=(xtest(:,i)-avg_training(i))/st_training(i);
xtest(:,i)= temp;
end
%hypothesis
h_test=[];
for i=1:m_test
h_test(i) = theta_new*xtest(i,:)';
end
J_test = (1/(2*m_test))*sum((h_test' -target_test).^2);

%Normal method
target_normal=(Variablest(:,3)-mean(Variablest(:,3)))/std(Variablest(:,3));
 m_normal = length(target_normal);
n_normal = length(Variablest(1,:));
for i=1:n_normal
temp=(Variablest(:,i)-mean(Variablest(:,i)))/std(Variablest(:,i));
Variablest(:,i)= temp;
end
%removing the unwanted variables 
Variablesw =Variablest(:,[4:n_normal]);
%x0
Variables = [ones(m_normal,1) Variablesw];
%thetas:
theta_normal=inv(Variables'*Variables)*(Variables')*target_normal;
for i=1:m_all
h_normal(i) = theta_normal'*Variables(i,:)';
end
 J_normal = (1/(2*m_normal))*sum((h_normal' -target_normal).^2);
