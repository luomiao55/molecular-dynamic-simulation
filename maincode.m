N=12;Lx=8;Ly=8;dt=0.02;vmax=1;t=200; ncal=t/dt;a=0;j=1;k=0;j1=3;
X=zeros(N,ncal);
Y=zeros(N,ncal);
Vx=zeros(N,ncal);
Vy=zeros(N,ncal);
Ax=zeros(N,N);
Ay=zeros(N,N);
ax=zeros(N,ncal);
ay=zeros(N,ncal);
R=zeros(N,N);
Dx=zeros(N,N);
Dy=zeros(N,N);
X(:,1)=[0.5*Lx/8
    1.5*Lx/8
    2.5*Lx/8
    3.5*Lx/8
    0.5*Lx/8
    1.5*Lx/8
    2.5*Lx/8
    3.5*Lx/8
    0.5*Lx/8
    1.5*Lx/8
    2.5*Lx/8
    3.5*Lx/8];
Y(:,1)=[0.5*Ly/3
    0.5*Ly/3
    0.5*Ly/3
    0.5*Ly/3
    1.5*Ly/3
    1.5*Ly/3
    1.5*Ly/3
    1.5*Ly/3
    2.5*Ly/3
    2.5*Ly/3
    2.5*Ly/3
    2.5*Ly/3];
for n=1:12
    Vx(n,1)=vmax*(2*rand(1,1)-1);
    Vy(n,1)=vmax*(2*rand(1,1)-1);
end
z=mean(Vx(:,1));
x=mean(Vy(:,1));
Vx(:,1)=Vx(:,1)-z;
Vy(:,1)=Vy(:,1)-x;
while (k<=N-1)
    k=k+1;l=0;
    while (l<=N-1)
        l=l+1;
         Dx(k,l)=X(k,1)-X(l,1);
         Dy(k,l)=Y(k,1)-Y(l,1);
           if (abs(Dx(k,l))>0.5*Lx)
              Dx(k,l)=Dx(k,l)-sign(Dx(k,l))*Lx;
           end
           if (abs(Dy(k,l))>0.5*Ly)
              Dy(k,l)=Dy(k,l)-sign(Dy(k,l))*Ly;
           end
           R(k,l)=(Dx(k,l)^2+Dy(k,l)^2)^0.5;
           if (k==l)
               Ax(k,l)=0;Ay(k,l)=0;
           else
               Ax(k,l)=24*Dx(k,l)*(2/(R(k,l)^6)-1)/(R(k,l)^8);
               Ay(k,l)=24*Dy(k,l)*(2/(R(k,l)^6)-1)/(R(k,l)^8);
               ax(k,1)=sum(Ax(k,:));
               ay(k,1)=sum(Ay(k,:));
           end 
    end     
end

while (j<=2)
    j=j+1;a=0;
        X(:,j)=X(:,j-1)+Vx(:,j-1)*dt+0.5*ax(:,j-1)*dt^2;
        Y(:,j)=Y(:,j-1)+Vy(:,j-1)*dt+0.5*ay(:,j-1)*dt^2;
        for q=1:N
           while (X(q,j)<0)
               X(q,j)=X(q,j)+Lx;
           end   
           while (X(q,j)>Lx)
               X(q,j)=X(q,j)-Lx;
           end    
           while (Y(q,j)<0)
              Y(q,j)=Y(q,j)+Ly;
           end
           while (Y(q,j)>Ly)
             Y(q,j)=Y(q,j)-Ly;
           end
        end
        while (a<=N-1)
           a=a+1;b=0;
           while (b<=N-1)
           b=b+1;
           Dx(a,b)=X(a,j)-X(b,j);
           Dy(a,b)=Y(a,j)-Y(b,j);
           if (abs(Dx(a,b))>0.5*Lx)
              Dx(a,b)=Dx(a,b)-sign(Dx(a,b))*Lx;
           end
           if (abs(Dy(a,b))>0.5*Ly)
              Dy(a,b)=Dy(a,b)-sign(Dy(a,b))*Ly;
           end
           R(a,b)=(Dx(a,b)^2+Dy(a,b)^2)^0.5;
           if (a==b)
               Ax(a,b)=0;Ay(a,b)=0;
           else
               Ax(a,b)=24*Dx(a,b)*(2/(R(a,b)^6)-1)/(R(a,b)^8);
               Ay(a,b)=24*Dy(a,b)*(2/(R(a,b)^6)-1)/(R(a,b)^8);
               ax(a,j)=sum(Ax(a,:));
               ay(a,j)=sum(Ay(a,:));
           end
           end
        end       
        Vx(:,j)=Vx(:,j-1)+0.5*dt*(ax(:,j)+ax(:,j-1));
        Vy(:,j)=Vy(:,j-1)+0.5*dt*(ay(:,j)+ay(:,j-1));
end
while (j1<=ncal-1)
    j1=j1+1;a1=0;
        X(:,j1)=2*X(:,j1-1)-X(:,j1-2)+ax(:,j1-1)*dt^2;
        Y(:,j1)=2*Y(:,j1-1)-Y(:,j1-2)+ay(:,j1-1)*dt^2;
        for q1=1:N
           while (X(q1,j1)<0)
               X(q1,j1)=X(q1,j1)+Lx;
           end   
           while (X(q1,j1)>Lx)
               X(q1,j1)=X(q1,j1)-Lx;
           end    
           while (Y(q1,j1)<0)
              Y(q1,j1)=Y(q1,j1)+Ly;
           end
           while (Y(q1,j1)>Ly)
             Y(q1,j1)=Y(q1,j1)-Ly;
           end
        end
        while (a1<=N-1)
           a1=a1+1;b1=0;
           while (b1<=N-1)
           b1=b1+1;
           Dx(a1,b1)=X(a1,j1)-X(b1,j1);
           Dy(a1,b1)=Y(a1,j1)-Y(b1,j1);
           if (abs(Dx(a1,b1))>0.5*Lx)
              Dx(a1,b1)=Dx(a1,b1)-sign(Dx(a1,b1))*Lx;
           end
           if (abs(Dy(a1,b1))>0.5*Ly)
              Dy(a1,b1)=Dy(a1,b1)-sign(Dy(a1,b1))*Ly;
           end
           R(a1,b1)=(Dx(a1,b1)^2+Dy(a1,b1)^2)^0.5;
           if (a1==b1)
               Ax(a1,b1)=0;Ay(a1,b1)=0;
           else
               Ax(a1,b1)=24*Dx(a1,b1)*(2/(R(a1,b1)^6)-1)/(R(a1,b1)^8);
               Ay(a1,b1)=24*Dy(a1,b1)*(2/(R(a1,b1)^6)-1)/(R(a1,b1)^8);
               ax(a1,j1)=sum(Ax(a1,:));
               ay(a1,j1)=sum(Ay(a1,:));
           end
           end
        end       
     
end
m=moviein(ncal);
figure(11);
plot(X(:,1),Y(:,1),'b.','Markersize',30);
axis([0 8 0 8]);
axis()
axis square
figure(12);
for n=1:ncal
h=plot(X(:,n),Y(:,n),'b.','Markersize',30);
axis([0 8 0 8]);
axis();
axis square;
grid off;
drawnow;
m(n)=getframe;
end
 for i=1:10
figure(i);
plot(X(:,1000*i),Y(:,1000*i),'b.','Markersize',30);
axis([0 8 0 8]);
axis()
axis square
grid off
 end

