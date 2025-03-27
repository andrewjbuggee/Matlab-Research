%% Random Varbiable with possibilities 1 and -1

% discrete uniform random variable with two possible outputs: 1 and -1


function X = rand_plusORminus_one(size_row,size_col)

X = rand(size_row,size_col);

index_LessThan = X>=0 & X<(1/2);
index_GreaterThan = X>=(1/2) & X<=1;


X(index_LessThan)=-1;           % half of the values are set to -1
X(index_GreaterThan)=1;        % the other half are set to 1


end
