
L_Q4 = data1(:, end) 
L_Q3 = data1(:, end-1) 
L_Q2 = data1(:, end-2) % 2022 year loss degree
L_Q1 = data1(:, end-3) 

L_Q5=L_Q3+L_Q4 % 1-10Y increase degree
T=data1(:, end-5) % total output
R=(data1(:, end-4)./T)' % emission intensity
X=data1(:, end-6) % final demand
A=data1(1:n,1:n)./(T') 


%%
S=zeros(n)
% forward 
for i=1:n
    b=['F',num2str(i)];
    eval([b,'=S;'])
end
% backward
for j=1:n
    a=['B',num2str(j)];
    eval([a,'=S;'])
end
    
for m=1:n
    for k=m
    eval(['F',num2str(m),'(k,:)','=','1',';']);
        end
end

for f=1:n
   eval(['F',num2str(f),'(f,f)','=','0',';']);
        end
    
for m1=1:n
    for k1=m1
    eval(['B',num2str(m1),'(:,k1)','=','1',';']);
        end
end
for f=1:n
    
    eval(['B',num2str(f),'(f,f)','=','0',';']);
end


num_matrices = n;
XFp = cell(1, num_matrices);


for o = 1:num_matrices
    
    matrix_name = ['F', num2str(o)];
    
    
    XFp{o} = eval(matrix_name);
end

XBp = cell(1, num_matrices);

for o = 1:num_matrices

    matrix_name = ['B', num2str(o)];
    
    XBp{o} = eval(matrix_name);
end


AFNew=cellfun(@(x)(ones(n)-x.*repmat(L_Q5,1,n)).*A,XFp,'UniformOutput',false)
ABNew=cellfun(@(x)(ones(n)-x.*repmat(L_Q5,1,n)).*A,XBp,'UniformOutput',false)
%%
m=12 


ResultList1 = zeros(m,3);
i=12 % The first (forward) or last (backward) sector number at the beginning of the path



visited = false(n);


for k = 1:m
    
    max_val = -inf;
    max_i1 = 0;
    max_j = 0;
    
 
    for i1=1:n
        for j = 1:n
           
            if visited(i1, j)
                continue;
            end
            
       
            val = -(R(i1) * A(i1,j) * X(j)-R(i1) * AFNew{i}(i1,j) * X(j));
            
           
            if val > max_val
                max_val = val;
                max_i1 = i1;
                max_j = j;
            end
        end
    end
        
        
       
        ResultList1(k, 1) = max_val;
        ResultList1(k, 2) = max_i1;
        ResultList1(k, 3) = max_j;
        
        
        visited(max_i1, max_j) = true;
end


ResultList1_1 = zeros(m, 3);

visited = false(n,n);


for k = 1:m
    
    max_val = -inf;
    max_i1 = 0;
    max_j = 0;
    
    
    for i1=1:n
        for j = 1:n
          
            if visited(j, i1)
                continue;
            end
            
           
            val = -(R(j) * A(j,i1) * X(i1)-R(j) * ABNew{i}(j,i1) * X(i1));
            
            
            if val > max_val
                max_val = val;
                max_i1 = i1;
                max_j = j;
            end
        end
    end
        
        
       
        ResultList1_1(k, 1) = max_val;
        ResultList1_1(k, 2) = max_j;
        ResultList1_1(k, 3) = max_i1;
        
       
        visited(max_j, max_i1) = true;
end

%%
% layer 2
m=20
ResultList2 = zeros(m, 4);

visited = false(n, n, n);


for k = 1:m
    
    max_val = -inf;
    max_i1 = 0;
    max_j1 = 0;
    max_j2 = 0;
    
    
    for i1=1:n
        for j1 = 1:n
            for j2 = 1:n
               
                if visited(i1, j1, j2)
                    continue;
                end
                
               
                val = -(R(i1) * A(i1,j1) * A(j1,j2) * X(j2)-R(i1) * AFNew{i}(i1,j1) * AFNew{i}(j1,j2) * X(j2));
                
               
                if val > max_val
                    max_val = val;
                    max_i1 = i1;
                    max_j1 = j1;
                    max_j2 = j2;
                end
            end
        end
    end
        
        

        ResultList2(k, 1) = max_val;
        ResultList2(k, 2) = max_i1;
        ResultList2(k, 3) = max_j1;
        ResultList2(k, 4) = max_j2;
        
        
        visited(max_i1, max_j1, max_j2) = true;
end

% backward : ResultList'n'_1, n means the layer number
ResultList2_1 = zeros(m, 4);

visited = false(n, n, n);


for k = 1:m
    
    max_val = -inf;
    max_i1 = 0;
    max_j1 = 0;
    max_j2 = 0;
    
    
    for i1=1:n
        for j1 = 1:n
            for j2 = 1:n
                
                if visited(j1, j2,i1)
                    continue;
                end
                
              
                val = -(R(j1) * A(j1,j2) * A(j2,i1) * X(i1)-R(j1) * ABNew{i}(j1,j2) * ABNew{i}(j2,i1) * X(i1));
                
                
                if val > max_val
                    max_val = val;
                    max_i1 = i1;
                    max_j1 = j1;
                    max_j2 = j2;
                end
            end
        end
    end
        
        
      
        ResultList2_1(k, 1) = max_val;
        ResultList2_1(k, 2) = max_j1;
        ResultList2_1(k, 3) = max_j2;
        ResultList2_1(k, 4) = max_i1;
        
        
        visited(max_j1, max_j2, max_i1) = true;
end
%%
% layer 3


ResultList3 = zeros(m, 5);


visited = false(n, n, n, n);

for k = 1:m
    
    max_val = -inf;
    max_indices = zeros(1, 4);
    
   
    for i1=1:n
        for j1 = 1:n
            for j2 = 1:n
                for j3 = 1:n
                
                    if visited(i1, j1, j2, j3)
                        continue;
                    end
                    
                  
                    val = -(R(i1) * A(i1,j1) * A(j1,j2) * A(j2,j3) * X(j3)-R(i1) * AFNew{i}(i1,j1) * AFNew{i}(j1,j2) * AFNew{i}(j2,j3) * X(j3));
                    
                    
                    if val > max_val
                        max_val = val;
                        max_indices = [i1, j1, j2, j3];
                    end
                end
            end
        end
    end
    
    
    ResultList3(k, :) = [max_val, max_indices];
    
   
    visited(max_indices(1), max_indices(2), max_indices(3), max_indices(4)) = true;
end


ResultList3_1 = zeros(m, 5);


visited = false(n, n, n, n);


for k = 1:m
    
    max_val = -inf;
    max_indices = zeros(1, 4);
    
   
    for i1=1:n
        for j1 = 1:n
            for j2 = 1:n
                for j3 = 1:n
                   
                    if visited(j1, j2, j3,i1)
                        continue;
                    end
                    
                
                    val = -(R(j1) * A(j1,j2) * A(j2,j3) * A(j3,i1) * X(i1)-R(j1) * ABNew{i}(j1,j2) * ABNew{i}(j2,j3) * ABNew{i}(j3,i1) * X(i1));
                    
                  
                    if val > max_val
                        max_val = val;
                        max_indices = [j1, j2, j3,i1];
                    end
                end
            end
        end
    end
  
    ResultList3_1(k, :) = [max_val, max_indices];
    
   
    visited(max_indices(1), max_indices(2), max_indices(3), max_indices(4)) = true;
end

%%
% layer 4


ResultList4 = zeros(m, 6);

visited = false(n, n, n, n, n);


for k = 1:m
    
    max_val = -inf;
    max_indices = zeros(1, 5);
    
    
    for i1=1:n
        for j1 = 1:n
            for j2 = 1:n
                for j3 = 1:n
                    for j4 = 1:n
                       
                        if visited(i1, j1, j2, j3, j4)
                            continue;
                        end
                        
                        
                        val = -(R(i1) * A(i1,j1) * A(j1,j2) * A(j2,j3) * A(j3,j4) * X(j4)-R(i1) * AFNew{i}(i1,j1) * AFNew{i}(j1,j2) * AFNew{i}(j2,j3) * AFNew{i}(j3,j4) * X(j4));
                        
                        
                        if val > max_val
                            max_val = val;
                            max_indices = [i1, j1, j2, j3, j4];
                        end
                    end
                end
            end
        end
    end

    
    
    ResultList4(k, :) = [max_val, max_indices];
    
    
    visited(max_indices(1), max_indices(2), max_indices(3), max_indices(4), max_indices(5)) = true;
end


ResultList4_1 = zeros(m, 6);

visited = false(n, n, n, n, n);


for k = 1:m
    
    max_val = -inf;
    max_indices = zeros(1, 5);
    
 
    for i1=1:n
        for j1 = 1:n
            for j2 = 1:n
                for j3 = 1:n
                    for j4 = 1:n
                        
                        if visited(j1, j2, j3, j4,i1)
                            continue;
                        end
                        
                      
                        val = -(R(j1) * A(j1,j2) * A(j2,j3) * A(j3,j4) * A(j4,i1) * X(i1)-R(j1) * ABNew{i}(j1,j2) * ABNew{i}(j2,j3) * ABNew{i}(j3,j4) * ABNew{i}(j4,i1) * X(i1));
                        
                        
                        if val > max_val
                            max_val = val;
                            max_indices = [j1, j2, j3, j4, i1];
                        end
                    end
                end
            end
        end
    end

    
    
    ResultList4_1(k, :) = [max_val, max_indices];
    
    
    visited(max_indices(1), max_indices(2), max_indices(3), max_indices(4), max_indices(5)) = true;
end

%%
% RF1: 1st layer forward path
RF1=[ResultList1,zeros(12,3)]
RF2=[ResultList2,zeros(20,2)]
RF3=[ResultList3,zeros(20,1)]
RF=[RF1;RF2;RF3;ResultList4]

% RB1: 1st layer backward path
RB1=[ResultList1_1,zeros(12,3)]
RB2=[ResultList2_1,zeros(20,2)]
RB3=[ResultList3_1,zeros(20,1)]
RB=[RB1;RB2;RB3;ResultList4_1]
