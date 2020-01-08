function [par]=wylicz_param(C,par_set,a,b,c,d,k,l)

% choose input parameter set: 1 - standard nondimnesional numbers, 2 -
% dimensional, 3- non-standard nondimensional

% TO DO: add here some test to check if the parameters have sens

switch par_set
    case 1
        St=a;
        Sv=b;
        teta=c;
        A=d;
        
        tau_p=(C.nu*St*Sv^2/(C.g^2*sin(teta)^2*A))^(1/3);
        R=(9*C.nu*tau_p*C.ro_a/(2*C.ro_p))^(1/2);
        delta=(C.nu^2*Sv/(A^2*C.g*sin(teta)*St))^(1/3);
        
    case 2
        R=a;
        delta=b;
        teta=c;
        A=d;
        
        tau_p=(2*C.ro_p*R^2)/(9*C.nu*C.ro_a);
        St=tau_p*C.nu/(delta^2*A);
        Sv=C.g*delta*A*tau_p*sin(teta)/(C.nu);

end

gamma=2*C.nu/delta^2;
tau_f=delta^2*A/C.nu;
tau_g=tau_f/Sv;

Fr=St/Sv;
B=(C.rs*C.nu^2/(C.g*2^8*pi^3))*((sin(teta)*delta.^3)^(-1)); % 10.01.19 to zostalo zmienione z przeliczenia uzywajacego St,Sv i A! 
z_b=C.g*tau_p*cos(teta)./gamma;

D=k*delta; %radius
Z=l*D; % half-length [m]


par.St=St;
par.Sv=Sv;
par.A=A;
par.teta=teta;
par.gamma=gamma;
par.tau_f=tau_f;
par.tau_g=tau_g;
par.tau_p=tau_p;
par.R=R;
par.delta=delta;
par.v_f=delta/tau_f;
par.Fr=Fr;
par.B=B;
par.z_b=z_b;
par.D=D;
par.Z=Z;
end