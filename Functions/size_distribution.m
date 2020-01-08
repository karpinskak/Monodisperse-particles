function [R] = size_distribution(N,Rmin,dispersion,varargin)
%SIZE_DISTRIBUTION Summary of this function goes here
%   varargin = Rmean or Rmean and Rdev or ProbabilityDistribution
%   
%   N - number of particles, dispersion is in {0,1,2,3} as parameter "disp".
%   If dispersion is 3 you need to declare probability distribution PD in a
%   form of Mx2 matrix, where 1st column gives and 2nd gives probabilities.
%   If dispersion is 2 you need to give Rmean and Rdev parameters.
%   For disp 1 you give Rmax.

switch dispersion
    case 0 % monodispersed
        if size(varargin)~=1
            error('Wrong nr of arguments')
        end
        R=zeros(N,1)+varargin{1};
    case 1 % uniform
        if size(varargin)~=1
            error('Wrong nr of arguments')
        end
        Rmax=varargin{1};
        R=rand(N,1)*(Rmax-Rmin)+Rmin;
        
    case 2 % gaussian
        if size(varargin)~=2
            error('Wrong nr of arguments')
        end
        Rmean=varargin{1};
        Rdev=varargin{2};
        R=randn(N,1)*Rdev+Rmean;
        %Usuniecie malych promieni
        while numel(nonzeros(R<=Rmin))~=0
            l_kropel_ujemnych=numel(nonzeros(R<=Rmin));
            R(R<=Rmin)=randn(l_kropel_ujemnych,1)*Rdev+Rmean;
        end
        
    case 3 % arbitrary
        if size(varargin)~=1
            error('Wrong nr of arguments')
        end
        PD=varargin{1};
     
        if abs(sum(PD(:,2))-1)>10^(-4)
            error('Wrong data for probability distribution')
        end
     
        D = cumsum(PD(:,2));
        D(end)=1; % So there is no error due to certain accuracy of the double number
        losowe = rand(N,1);
        prob = @(r) find(r<D,1,'first'); % find the 1st index s.t. r<D(i);
        % Now this are your results of the random trials
        R = arrayfun(prob,losowe)*PD(end,1)/size(PD,1);
         while (numel(nonzeros(R<=Rmin))+numel(nonzeros(R>max(PD(:,1)))))~=0
            l_kropel_ujemnych=numel(nonzeros(R<=Rmin));
            los = rand(l_kropel_ujemnych,1);
            if numel(nonzeros(R<=Rmin))~=0
            R(R<=Rmin)=arrayfun(prob,los)*PD(end,1)/size(PD,1);
            end
            l_kropel_duzych=numel(nonzeros(R>max(PD(:,1))));
            los2 = rand(l_kropel_duzych,1);
            if numel(nonzeros(R>max(PD(:,1))))~=0
            R(R>max(PD(:,1)))=arrayfun(prob,los2)*PD(end,1)/size(PD,1);
            end
        end
end

