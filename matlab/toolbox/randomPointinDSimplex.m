function p=randomPointinDSimplex(d,N)

if nargin==1,
    N=1;
end
y=rand(d-1,N);
u=sort([zeros(1,N);y;ones(1,N)],'ascend');
p=diff(u);

end