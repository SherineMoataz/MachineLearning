clc
clear all
close all

DataTable = readtable('house_prices_data_training_data.csv');
X = table2array(DataTable(1:17999,4:21));

for i=1:length(X(1,:))
    X(:,i) = (X(:,i)-mean(X(:,i)))/std(X(:,i)); 
end

k = 10;

[ rows columns ] = size(X);


costVec = zeros(1,10);


for k = 1:10
    
    centroids = zeros( rows , columns);
    initial_index = randperm(rows);
    centroids = X(initial_index(1:k),:);
    
    prevCentroids = zeros(size(centroids));
    
    indices = zeros(size(X,1), 1);
    distance = zeros(rows,k);
    
    flag = true;
    iterations = 0;
    
    while(flag)

               
        
        
        for i = 1:rows
            for j = 1:k
                distance(i, j) = sum((X(i,:) - centroids(j, :)).^2);
            end
        end
        
        for i = 1:rows
            indices(i) = find(distance(i,:)==min(distance(i,:)));
        end
        
        for i = 1 : k
            
            clustering = X(find(indices == i), :);
            centroids(i, :) = mean(clustering);
            cost = 0;
            
            for z = 1 : size(clustering,1)
                cost = cost + sum((clustering(z,:) - centroids(i,:)).^2)/rows;
            end
            
            costVec(1,k) = cost;
            
        end
        
        if prevCentroids == centroids
            flag = false;
        end
        
        prevCentroids = centroids;
        
        iterations = iterations + 1;
        
    end
end
[ o bestKvalue ] = min(costVec);

noClusters = 1:10;

plot(noClusters, costVec)