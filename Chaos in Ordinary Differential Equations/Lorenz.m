

function rhs = Lorenz(t,x,dummy,sigma,b,r)

rhs = [
        sigma*(-x(1)+x(2))
        -x(1)*x(3)+r*x(1)-x(2)
        x(1)*x(2)-b*x(3)
        ];

