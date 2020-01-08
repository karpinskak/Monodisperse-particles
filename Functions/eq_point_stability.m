function  [stability] = eq_point_stability(r,A,varargin)
%EQ_POINT_STABILITY
% r is a vector of radial positions of equillibrium points
%'stability' can have values 1,2 which are stable and unstable, NaN means
%not enough data to regognize

fi_war=fi(r);
stability=zeros(size(r));

for j=1:numel(r)
    
    if fi_war(j)<0
        if nargin>2
            war=A/(St*abs(fi_war(j)));
        else
            war=NaN;
        end
    elseif fi_war(j)>0
        war=A/sqrt(fi_war(j));
    end
    
    if war>=1 && ~isnan(war)==1
        stability(j)=1;
    elseif war<1 && ~isnan(war)==1
        stability(j)=2;
    else
        stability(j)=NaN;
    end
end
end

function wart=fi(r)
part1=(1-exp(-r.^2/2))./r.^2;
part3=exp(-r.^2/2);
wart=part1.*(part1-part3)/(2*pi)^2;
end

