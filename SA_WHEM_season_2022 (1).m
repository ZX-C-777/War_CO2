S_Q4_0= data1(:,end)
S_Q3_0= data1(:,end-1)
S_Q2_0= data1(:,end-2)
S_Q1_0= data1(:,end-3)

L_Q4 = data1(:, end-4).* S_Q4_0
L_Q3 = data1(:, end-5).* S_Q3_0 
L_Q2 = data1(:, end-6).* S_Q2_0 
L_Q1 = data1(:, end-7).* S_Q1_0 
X=data1(:, end-9) % total output
E_I=diag(data1(:, end-8)./X) % emission intensity
Y=data1(:, end-10) % final demand
A=data1(1:n,1:n)./(X') % direct requirement coefficient 
Linv=inv(eye(n)-A) % Leontief inverse matrix
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
    
%%

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

%%
WA=cellfun(@(x)Linv-inv(eye(n)-((ones(n)-repmat(L_Q1,1,n).*x).*A)),XFp,'UniformOutput',false)
WB=cellfun(@(x)Linv-inv(eye(n)-((ones(n)-repmat(L_Q2,1,n).*x).*A)),XFp,'UniformOutput',false)
WC=cellfun(@(x)Linv-inv(eye(n)-((ones(n)-repmat(L_Q3,1,n).*x).*A)),XFp,'UniformOutput',false)
WD=cellfun(@(x)Linv-inv(eye(n)-((ones(n)-repmat(L_Q4,1,n).*x).*A)),XFp,'UniformOutput',false)

WA2=cellfun(@(x)Linv-inv(eye(n)-((ones(n)-repmat(L_Q1,1,n).*x).*A)),XBp,'UniformOutput',false)
WB2=cellfun(@(x)Linv-inv(eye(n)-((ones(n)-repmat(L_Q2,1,n).*x).*A)),XBp,'UniformOutput',false)
WC2=cellfun(@(x)Linv-inv(eye(n)-((ones(n)-repmat(L_Q3,1,n).*x).*A)),XBp,'UniformOutput',false)
WD2=cellfun(@(x)Linv-inv(eye(n)-((ones(n)-repmat(L_Q4,1,n).*x).*A)),XBp,'UniformOutput',false)

% CO2 reduction effects��forward��
CO2F1=cellfun(@(x)sum(E_I*x*Y),WA,'UniformOutput',false)
CO2F2=cellfun(@(x)sum(E_I*x*Y),WB,'UniformOutput',false)
CO2F3=cellfun(@(x)sum(E_I*x*Y),WC,'UniformOutput',false)
CO2F4=cellfun(@(x)sum(E_I*x*Y),WD,'UniformOutput',false)


% CO2 reduction effects��backward��
CO2B1=cellfun(@(x)sum(E_I*x*Y),WA2,'UniformOutput',false)
CO2B2=cellfun(@(x)sum(E_I*x*Y),WB2,'UniformOutput',false)
CO2B3=cellfun(@(x)sum(E_I*x*Y),WC2,'UniformOutput',false)
CO2B4=cellfun(@(x)sum(E_I*x*Y),WD2,'UniformOutput',false)

% Results: 
% 1: 2022Q1, 2:2022Q2, 3:Q3, 4:Q4.
CO2FB=[CO2F1',CO2F2',CO2F3',CO2F4',CO2B1',CO2B2',CO2B3',CO2B4']
