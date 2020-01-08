function [points, varargout] = eq_points(Const,A,Sv,varargin)
%EQ_POINTS Calculates all eq. points of Burgers vortex and its stability
%   'points' is nx2 matrix where n is the number of eq. points and 2 gives X
%   and Y positions.
%   'stability' can have values: 0-unstable focus, 1-stable focus,2-saddle,
%   3-stable node



syms r
assume(r>=0)

if A>=0.02176
    
    r_eq0=double(vpasolve(r*sqrt(1+((1-exp(-r^2/2))/(2*pi*A*r^2))^2)-Sv/A==0,r));
    point_r=unique(r_eq0);
    points=zeros(numel(point_r),2);
    chi=(1-exp(-point_r.^2/2))./(2*pi*A*point_r.^2);
    points(:,1)=Sv*chi./(A*(1+chi.^2));
    points(:,2)=-Sv./(A*(1+chi.^2));
    
    if point_r>Const.rs
        stability=3;
    else
        if nargin==4
            St=varargin{1};
            fi_war=fi(point_r);
            war=A/(St*abs(fi_war));
            if war>=1
                stability=1;
            else
                stability=0;
            end
        else
            stability=NaN;
        end
    end
    varargout{1}=stability;
else
    rmax=200;
    r_bezw=0:0.0001:0.5*rmax;
    [~,kolejnosc]=sort(abs(krzywa(r_bezw,A,Sv)));
    r_pos=r_bezw(kolejnosc);
    r0=[0,r_pos(1:15),1.5,2,10,rmax,10*rmax];
    for n=1:numel(r0)
        x0=r0(n);
        r_equ=double(vpasolve(r.*sqrt(1+((1-exp(-r.^2/2))./(2*pi*A*r.^2)).^2)-Sv/A==0,r,x0));
        if numel(r_equ)==1
            r_eq0(n)=r_equ;
        elseif numel(r_equ)==0
            r_eq0(n)=NaN;
        else
            error('za duzo rozwiazan')
        end
    end
    r_eq0(isnan(r_eq0))=[];
    point_r=unique(r_eq0);
    if numel(point_r)>3
        error(['Sth wrong with nr of eqillibrium points:' num2str(numel(points))])
    end
    
    points=zeros(numel(point_r),2);
    chi=(1-exp(-point_r.^2/2))./(2*pi*A*point_r.^2);
    points(:,1)=Sv*chi./(A*(1+chi.^2));
    points(:,2)=-Sv./(A*(1+chi.^2));
    stability=zeros(size(point_r));
    
    for j=1:numel(point_r)
        if point_r(j)>Const.rs
            fi_war=fi(point_r(j));
            war=A/sqrt(fi_war);
            if war>=1
                stability(j)=3;
            else
                stability(j)=2;
            end
        else
            if nargin==4
                St=varargin{1};
                fi_war=fi(point_r(j));
                war=A/(St*abs(fi_war));
                if war>=1
                    stability(j)=1;
                else
                    stability(j)=0;
                end
            else
                stability(j)=NaN;
            end
        end
    end
    varargout{1}=stability;
end
end
% 25.08.2019 tu byl blad w definicji funkcji!!!
% wart=(1-exp(-r.^2/2)).*((1-exp(-r.^2/2))/r.^2-exp(-r.^2/2))./(r*(2*pi)).^2;

function wartosc=krzywa(r,A,Sv)
wartosc=r.*sqrt(1+((1-exp(-r.^2/2))./(2*pi*A*r.^2)).^2)-Sv/A;
end

function wart=fi(r)
part1=(1-exp(-r.^2/2))./r.^2;
part3=exp(-r.^2/2);
wart=part1.*(part1-part3)/(2*pi)^2;
end

