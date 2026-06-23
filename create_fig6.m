% Set the experimental parameters
th0=[1e24 0.8 0.0025 2 2 2e4];

% nperturb=0.25;
% thini=(1+nperturb*randn(size(th0))).*th0;

fields={'DEL','g','z2','dydx','NyNx','blurs','kiso', 'taper'};
defstruct('params',fields,...
	    {[2670 630],9.81,35000,[20 20]*1e3,[128 128],3,NaN, 0});

% Do one SIMULROS simulation and MLEROS recovery
[Hx,Gx,th0,params,k,Hk,Gk,Sb,Lb]=simulros(th0,params);
[thhat,~,logli,th0,scl,p,e,o,gr,hs] = mleros(Hx,Gx,th0,params,[],[],[]);



%Plotting
clf
[~, ah, ha] = krijetem(subnum(2,5));
axes(ah(1))

axes(ah(2))

axes(ah(3))

axes(ah(4))

axes(ah(5))

axes(ah(6))

axes(ah(7))

axes(ah(8))

axes(ah(9))

axes(ah(10))



