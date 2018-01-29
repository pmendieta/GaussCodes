/* In this file I solve for the partial equilibrium 
   Hugget Economy with Endogenous Grid Method Here I explore incoporating a
   restriction to maintain always positive consumption
   in the arguments of marginal utility 
   --Sergio Salas Landeau, solved by Pablo Mendieta--*/

new;
clear all; 
cls;

Np   = 300;                                       @--number of evaluation points---@
Nk   = 1000;                                      @--maximum number of iterations--@

beta = 0.99322;                                   @--discount factor-------@
gama = 1.5;                                       @--param utility function@ 
q    = 1.0135023;                                 @--price of claims-------@
QQ     = (0.5~0.5)|(0.075~0.925);                 @--transition matrix-----@
yl   = 0.1;                                       @--low endowment---------@
yh   = 1;                                         @--high endowment--------@

alow = -2;                                     @--borrowing constraint--@
aup  = 4;                                        @--upper bound assets----@
agrid= grid(alow,aup,Np);                        @--grid for assets-------@

a0_l = yl.*ones(Np,1);//agrid;                             @--initial guess for asset policy, low shock---@
a0_h = yh.*ones(Np,1);//agrid;                             @--initial guess for asset policy, high shock--@



@--storing policy functions--@

al   = zeros(Np,Nk); al[.,1] = a0_l;
ah   = zeros(Np,Nk); ah[.,1] = a0_h;


als  = zeros(Np,Nk);                             @--storing a* low shock---@
ahs  = zeros(Np,Nk);                             @--storing a* high shock--@

Bl   = zeros(Np,Nk);                             @--storing RHS of Euler low endowment--@
Bh   = zeros(Np,Nk);                             @--storing RHS of Euler high endowment-@

for j(1,Nk-1,1);
    cl     = agrid+yl-q*al[.,j];
    ch     = agrid+yh-q*ah[.,j];
    dl     = cl .< 0;                               @--dummy, 1 if value negative-------@
    dh     = ch .< 0;                               @--dummy, 1 if value negative-------@

    Bl[.,j] = beta*(QQ[1,1]*substute(cl,dl,0.000001)^(-gama)+QQ[1,2]*substute(ch,dh,0.000001)^(-gama));
    Bh[.,j] = beta*(QQ[2,1]*substute(cl,dl,0.000001)^(-gama)+QQ[2,2]*substute(ch,dh,0.000001)^(-gama));
    als[.,j]= (Bl[.,j]/q)^(-1/gama)+q*agrid-yl;          @--solving for a^*(a'_i,yl)--@
    ahs[.,j]= (Bh[.,j]/q)^(-1/gama)+q*agrid-yh;          @--solving for a^*(a'_i,yh)--@

    for i(1,Np,1);

        if agrid[i] < als[1,j];
            al[i,j+1] = agrid[1];
        elseif agrid[i] >= als[Np,j];
            al[i,j+1] = agrid[Np];
        else;
            al[i,j+1] = lininter(als[.,j],agrid,agrid[i]);      
        endif;

        if agrid[i] < ahs[1,j];
            ah[i,j+1] = agrid[1];
        elseif agrid[i] >= als[Np,j];
            ah[i,j+1] = agrid[Np];
        else;
            ah[i,j+1] = lininter(ahs[.,j],agrid,agrid[i]);      
        endif;
    endfor;
    /* 
Here you should write the main code of the EGM

Guide from the code with EGM for a growth model with heterogeneous agents
    
     if agrid[i]<als[1];
                khat1_1[i] = agrid[1];
                lhat1_1[i] = 0;
                elseif agrid[i]>=als[Np];
                khat1_1[i] = agrid[Np];
                lhat1_1[i] = 0;
                else; 
                khat1_1[i] = lininter(als,agrid,agrid[i]);
                lhat1_1[i] = lininter(als,lam1,agrid[i]);
            endif;

            if agrid[i]<as0[1];
                khat0_1[i] = agrid[1];
                lhat0_1[i] = 0;
                elseif agrid[i]>as0[Np] ;
                khat0_1[i] = agrid[Np];
                lhat0_1[i] = 0;
                else;
                khat0_1[i] = lininter(as0,agrid,agrid[i]);
                lhat0_1[i] = lininter(as0,lam0,agrid[i]);
            endif;
        endfor;

*/

ac=abs((al[.,j+1]|ah[.,j+1])-(al[.,j]|ah[.,j]));

    if maxc(ac) <= 0.000001;
        print "convergence achieved";
        gl = al[.,j];
        gh = ah[.,j];
        print "number of iterations";
        lindex  = maxindc((Bl+beta*al)');            @--index for policy in low endowment enviroment--@
        hindex  = maxindc((Bh+beta*ah)');            @--index for policy in high endowment enviroment--@        
        print j;
        iter=j;
        break;
    endif;

endfor;

/*Here you need to apply EGM to find policies 

gh= a'(a,e_h), gl = a'(a,e_l)

And then use the code below */


for i(1,rows(gl),1);
    gl[i]=al[i,iter];
endfor;

for i(1,rows(gh),1);
    gh[i]=ah[i,iter];
endfor;

library pgraph;
graphset;
/*_pdate=2;*/
_pmcolor = {0,0,0,0,0,0,0,0,15};
_pcolor = {1,2,0,12};
_pgrid  = {1,0};
/*_plwidth= {7}; 
_plegctl= { 4, 5, 1.7, 4.5 };*/
_plegstr= " Low endowment\0"\
          " High endowment\0";
/*scale(0,0|12.5);*/
xlabel("Assets (a)");
ylabel("Next period assets (a')");

title("Policy Function (own replication based on Salas)");


xy(agrid,gl~gh);
//xy(kgrid,kpol~kpolana~(kgrid^alf+(1-delta)*kgrid));
waitc;
/*stop;*/


/* In this second part I compute the distributions */

p00=QQ[1,1];
p01=1-p00;
p10=QQ[2,1];
p11=1-p10;

Prob=(p00-1~1)|(p10~1);
b=(0|1);
epr = linsolve(b,Prob');

pr0=epr[1];                         @--ergodic probability low shock----------@                       
pr1=epr[2];                         @--ergodic probability high shock---------@


Md=3*Np;                            @--number of grid point for distribution--@
mgrid=grid(alow,aup,Md);            @--grid for distribution------------------@
Fini=(mgrid-alow)/(aup-alow);       @--initial guess for distribution---------@


Nk=300;                             @--maximum number iterations--------------@

F0=zeros(Md,Nk);                    @--matrix to store distributions----------@
F0[.,1]=pr0.*Fini;                  @--setting initial guess ergodic prob-----@

F1=zeros(Md,Nk);                    @--matrix to store distributions----------@
F1[.,1]=pr1.*Fini;                  @--setting initial guess------------------@

for i(1,Nk-1,1);
    for j(1,Md,1);
        if mgrid[j] < gl[1];
        m0 = lininter(mgrid,F0[.,i],lininter(gl,agrid,mgrid[j]));
        m1 = 0;
        F0[j,i+1] = p00*m0+p10*m1;
        F1[j,i+1] = p01*m0+p11*m1; 

        elseif mgrid[j] >= gh[1] and mgrid[j]<gl[Np];
        m0 = lininter(mgrid,F0[.,i],lininter(gl,agrid,mgrid[j]));
        m1 = lininter(mgrid,F1[.,i],lininter(gh,agrid,mgrid[j]));
        F0[j,i+1] = p00*m0+p10*m1;
        F1[j,i+1] = p01*m0+p11*m1; 

        elseif mgrid[j]>=gl[Np] and mgrid[j]<gh[Np];
        m0 = pr0;
        m1 = lininter(mgrid,F1[.,i],lininter(gh,agrid,mgrid[j]));
        F0[j,i+1] = p00*m0+p10*m1;
        F1[j,i+1] = p01*m0+p11*m1;
        
        elseif mgrid[j]>=gh[Np];
        m0 = pr0;
        m1 = pr1;
        F0[j,i+1] = p00*m0+p10*m1;
        F1[j,i+1] = p01*m0+p11*m1;

        endif;
    endfor;
    h = (F1[.,i+1]-F1[.,i])*(F1[.,i+1]-F1[.,i])';   
    if h<=0.00001;
    F1 = F1[.,i+1];
    F0 = F0[.,i+1];
    print "convergence achieved";
    print i;
    break;
    endif;
endfor;

F = F1+F0;

/* Excess demand of credit */ 

SA = (F1[2:Md]-F1[1:Md-1])'((mgrid[2:Md]+mgrid[1:Md-1])./2)+F1[1]*mgrid[1];
DA = (F0[2:Md]-F0[1:Md-1])'((mgrid[2:Md]+mgrid[1:Md-1])./2)+F0[1]*mgrid[1];

EXDA  = DA+SA;

print SA|DA|EXDA;
