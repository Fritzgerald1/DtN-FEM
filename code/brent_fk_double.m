function [ xmin, brent ] = brent_fk_double( ax,bx,cx,tol,kth1,kth2,mu1,mu2 )
CGOLD = 0.381966;
ZEPS = 1*10^(-10);
ITMAX = 100;

a = min(ax,cx);
b = max(ax,cx);
v = bx;
w = v;
x = v;
e = 0;
fx = double_fk(x,kth1,kth2,mu1,mu2);
fv = fx;
fw = fx;
for iter = 1:ITMAX
    xm = 0.5*(a+b);
    tol1 = tol*abs(x)+ZEPS;
    tol2 = 2*tol1;
    if abs(x-xm) <= (tol2-0.5*(b-a)) 
        break
    end
    if abs(e) > tol1
        r = (x-w)*(fx-fv);
        q = (x-v)*(fx-fw);
        p = (x-v)*q-(x-w)*r;
        q = 2*(q-r);
        if q > 0
            p = -p;
        end
            q = abs(q);
            etemp = e;
            e = d;
            if abs(p) >= abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)
                if x >= xm
                    e = a-x;
                else
                    e = b-x;
                end
                d = CGOLD*e;
            else
                d=p/q;
                u=x+d;
                if u-a < tol2 || b-u < tol2
                    d = sign(xm-x)*abs(tol1);
                end
            end
    else
        if x >= xm
            e = a-x;
        else
            e = b-x;
        end
        d = CGOLD*e;
    end
    if abs(d) >= tol1 
        u =x+d;
    else
        u = x + sign(d)*abs(tol1);
    end
    
    fu = double_fk(u,kth1,kth2,mu1,mu2);
    if fu <= fx
        if  u >= x
            a = x;
        else
            b = x;
        end
        v = w;
        fv =fw;
        w = x;
        fw = fx;
        x = u;
        fx = fu;
    else
        if u < x 
            a = u;
        else
            b = u;
        end
        if fu <= fw || w ==x 
            v =w;
            fv =fw;
            w = u;
            fw = fu;
        else
            if fu <= fv || v == x || v == w
                v = u;
                fv =fu;
            end
        end
    end                         
end
xmin = x;
brent = fx;

