function imgOut = imgRecover(imgIn, blkSize, numSample)
%Reading image
C_Original=imgRead(imgIn);
%Initial the predicted final output
C_Final=zeros(length(C_Original(:,1)),length(C_Original(1,:)));
block_size=blkSize;
k=numSample;
%count=1;
if length(C_Original(:,1))==512
    count2=16;
else
    count2=5;
end
    

for iteration=1:count2
C=C_Original((iteration-1)*(length(C_Original(:,1))/count2)+1:length(C_Original(:,1))/count2*iteration,:);
%Initialize height and width
H=length(C(:,1));
W=length(C(1,:));
%Generate blocks
%Initiate predicted B
B_Final=zeros(H,W);

T=zeros(H*W);


         


for x=1:H
	for y=1:W
		for u=0:H-1
			for v=1:W
				T((x-1)*W+y,u*W+v)=sqrt((u==0)/H)*sqrt((v==1)/W)*cos(pi*(u)*(2*x-1)/(2*H))*cos(pi*(v-1)*(2*y-1)/(2*W))+sqrt(2/H)*sqrt((v==1)/W)*cos(pi*(u)*(2*x-1)/(2*H))*cos(pi*(v-1)*(2*y-1)/(2*W))+sqrt(2/H)*sqrt(2/W)*cos(pi*(u)*(2*x-1)/(2*H))*cos(pi*(v-1)*(2*y-1)/(2*W))+sqrt((u==0)/H)*sqrt(2/W)*cos(pi*(u)*(2*x-1)/(2*H))*cos(pi*(v-1)*(2*y-1)/(2*W));

			end
		end

	end
end
for o=1:H/block_size
for p=1:W/block_size
%Sampling from blocks
%o=14;p=10;
for i=1:block_size^2
	block_y(i,1)=rem(i,block_size);
end
block_y(block_y==0)=block_size;
block_y=block_y+(p-1)*block_size;
temp=0;
for i=1:block_size^2
	if rem(i,block_size)==1
		temp=temp+1;
	end
	block_x(i,1)=temp;
end
block_x(end,1)=block_size;
block_x=block_x+(o-1)*block_size;
rand_ind=randsample(block_size*block_size,k);
sample_x=ceil(rand_ind/block_size);
sample_y=rem(rand_ind,block_size);
sample_y(sample_y==0)=block_size;
sample_x=sample_x+(o-1)*block_size;
sample_y=sample_y+(p-1)*block_size;

T_block=T((block_x-1)*W+block_y,:);
T_sample=T((sample_x-1)*W+sample_y,:);

C_sample=zeros(k,1);
C_block=zeros(block_size^2,1);
for i=1:k
	C_sample(i)=C(sample_x(i),sample_y(i));
end

for i=1:block_size^2
	C_block(i)=C(block_x(i),block_y(i));
end
%Construct loops for sampling in each block
%for j=1:(H/block_size*W/block_size)
%Testing lambda
lambda=11;
error_lambda=[];

%Cross validation
%for n=1:30
N_s=length(C_sample(:,1));
m=floor(N_s/6);
%Count the rows of the sample
rows=(1:N_s)';
%Index for testing set
rand_cv=randsample(N_s,m);
%Index for training set
rows(rand_cv,:)=[];

%OMP initialize
%Used in CV
%B=C_sample(rows,:);
B=C_sample;
F=B;
p=1;
%Used in CV
%A=T_sample(rows,:);
A=T_sample;
theta=zeros(1,H*W);
omega=[];
%Used in CV
%B_Testing=C_sample(rand_cv,:);
%A_Testing=T_sample(rand_cv,:);

%Loop for convergence
while p<lambda
	for i=1:(H*W)
		theta(i)=dot(A(:,i),F);
	end
s=find(theta==max(theta));
omega=cat(2,omega,s);
alpha=A(:,omega(1,:))\B;
F=B-A(:,omega(1,:))*alpha;

p=p+1;
end

%Reconstruct alpha
alpha_final=zeros(H*W,1);
for i=1:length(omega(1,:))
	alpha_final(omega(1,i),1)=alpha(i,1);
end

%Test the alpha on testing set
%B_pred=A_Testing*alpha_final;
%MSE=mean((B_pred-B_Testing).^2);
%Record error of lambda
%error_lambda(n)=MSE;
%end
%average_error=mean(error_lambda);
%final_error(count)=average_error;
%count=count+1;


%Calculate the predicted B with sparse alpha solutions
B_Recovered=T_block*alpha_final;

%Potential problem here
for i=1:block_size^2
	B_Final(block_x(i),block_y(i))=B_Recovered(i);
end




end
end

C_Final((iteration-1)*(length(C_Original(:,1))/count2)+1:length(C_Original(:,1))/count2*iteration,:)=B_Final;

end
%Filter
C_Final=medfilt2(C_Final);

%Show image
imgOut=C_Final;
imgShow(imgOut);





end