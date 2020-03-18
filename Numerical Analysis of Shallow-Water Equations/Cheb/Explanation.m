
%polynomial interpolation in equispaced and Chebychev points

%number of grid points
%N=12; 
N=24;

%test both methods
for i = 3:4
    %equispaced points
    if i==3, s='Equispaced points'; x=-1+2*(0:N)/N;end
    %Chebychev points
    if i==4, s='Chebychev points'; x=cos(pi*(0:N)/N); end
   
    %define function
    u=1./(1+16*x.^2);
    
    %Polynomial curve fitting/interpolaton
    %use polyfit to fit a Nth-degree polynomial
    p=polyfit(x,u,N);  
    
    %define the space domain
    xx = -1.01:0.005:1.01;

    %evaluation of interpolant
    %evaluates the polynomial p at each point in xx
    pp=polyval(p,xx); 
     
    %plotting
    subplot(2,2,i)
    %plot the grid points of the fucntion u
    plot(x,u,'r.','markersize',13)
    hold on
    %plot the polynomial approximation
    plot(xx,pp,'b')
    axis([-1.1 1.1 -1 1.5]), title(s)
    
    %calculate error
    %define the function on each point of the space domain
    uu=1./(1+16*xx.^2);
    error = norm(uu-pp,inf);
    text(-0.5, -0.5, ['max error = ' num2str(error)])
    
end



