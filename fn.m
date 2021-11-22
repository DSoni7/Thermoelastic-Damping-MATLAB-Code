function [y]=fn(A,B,C,Ak1,Ak2,Ak3,Ak4)
    y=-Ak1*A-Ak2*B-Ak3*(A^3)+Ak4*A*C;
end