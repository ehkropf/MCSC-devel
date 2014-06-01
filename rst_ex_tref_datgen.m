%% barbell resistor accuracy run
clear


NN = 4;
npts = 100;



%%
P = polygons(...
  polygon([ 0.25i -0.5+0.25i -0.5+0.5i -1.5+0.5i -1.5-0.5i ...
           -0.5-0.5i -0.5-0.25i 0.5-0.25i 0.5-0.5i 1.5-0.5i ...
           1.5+0.5i 0.5+0.5i 0.5+0.25i ]),...
  polygon([ -0.75 -0.75+0.25i -1.25+0.25i -1.25-0.25i -0.75-0.25i ]),...
  polygon([ 0.75 0.75-0.25i 1.25-0.25i 1.25+0.25i 0.75+0.25i ])...
);

[m vc vl beta] = polydat(P);
fixed = [0; 0];

X0 = [
   0.015086285340753; 0.015086267236461; -0.000166410519941;
   0.981624373360132; -0.000165795426656; -0.981624394245061;
   1.505892839263334; 1.539921029816123; 1.566343547709987;
   1.575590629565584; 1.602012635240045; 1.636039522502695;
   4.647145975676726; 4.681172812113288; 4.707594858922135;
   4.716841910053967; 4.743264213227991; 4.777292251033749;
   4.712716123001660; 6.400725960045659; 7.580050650064988;
   8.128381052774674; 9.307793339843695; 1.570507598767029;
   3.258592824151452; 4.437980801084158; 4.986310618112451;
   6.165658607077876;
];




%%
Rmean = zeros(length(NN),1);
Rmin = Rmean; Rmax = Rmean;
rsumSC = Rmean; rsumSRo = Rmean; rsumSRi = Rmean;
for n = 1:length(NN);
  N = NN(n);
  
  Xc = solve_parambd(P,N,fixed,X0);
  
  [f C z w rsumf] = make_mapbd(Xc,P,N,fixed);
  [c r t] = circDomain(Xc,vc);
  
  a = 0.2i;
  [g gp rsumgo rsumgi] = map_radslitring(c,r,a,N);
  
  zz = repmat(c(2:3),npts,1) + ...
    repmat(exp(2i*pi*(0:npts-1)'/(npts-1)),1,2)*diag(r(2:3));
  eta = g(zz);

  rout = abs(eta(:,1)); rin = abs(eta(:,2));
  Rmean(n) = 0.5*log(sum(rout)/sum(rin))/pi;
  Rmin(n) = 0.5*log(min(rout)/max(rin))/pi;
  Rmax(n) = 0.5*log(max(rout)/min(rin))/pi;
  
  rsumSC(n) = rsumf(end);
  rsumSRo(n) = rsumgo(end);
  rsumSRi(n) = rsumgi(end);
end

