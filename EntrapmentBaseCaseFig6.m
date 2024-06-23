%% EntrapmentBaseCaseFig6.m

%% Front Matter
%Simulated advancing simplified clinoform
%Related to the paper:-


%Onshore entrapment of seawater in coastal aquifers by rapid coastline progradation

%Vaughan R. Voller
%Nafis Sazeed
%Loc Luong
%Adrien Camille
%MichaelSteckler
%Kate Leary
%Mark Person

%THE CURRENT CODE IS SET TO RECOVER FIGURE 6B
%Computation elapse time on a Macbook is ~5mins

clear all


%% Basic Settings

%User Defined Data --see Table 2
kxval=5; % [m/day] Hydraulic conductivity in x-direction
aL=50; %  [m] Longitudinal dispersion
simtime=4911;; % [year] simulation time

%Physical properties and program variables held fixed in study see table 2
Breath =300000; % [m]  domain breath
cols=101; %number of equal spaced colums in domain
topslope=0.0005;% the water table slope
foreslope = 0.05;% the foreset slope 0.05

Delx=Breath/(cols-1); % [m] Fixed width between columns (300 m)
xcols=linspace(0,Breath,cols); %x-locations of columns
AR=80; %Aspect ratio of elements  AR=xdim/ydim
Dely=Delx/AR;% [m] Typical height of an element 37.5 m

delt=100; % [days] Time step

consea=1;  %Saturated Salt concentration in sea
rhorel=1.025; %Relative density of saturated saline rhosat/rhowater

eps=0.35; %Aquifer porosity
stov=0.001; %Aquifer storage

%Hydraulic conductivity m/day, orthotropic
% kxval--hydraulic conductivity in x-direction
kyval=kxval/100; %in y direction

%diffusion and dispersion coefficients to be used in anisotropic dispersion
Dmol= 0.00001;%  molecular diffusion  [m^2/day]
%aL=Longitudinal dispersion [m]
aT=aL/10; %Transverse

%% Geometrically Defined synthetic stratigraphy

%Models clinoform by assuming constant water-table and foreset slopes

sealevelin=300; %Initial Sea-level setting [m]
basement=250; % Depth of basement -- height of sea floor [m]
%foreslope--foreset slope, set above
%topslope--the Water table slope, set above

%initial position (m) of the location of delta toe, where the foreset intersects with the %basement
xtoein=(sealevelin-basement)/foreslope;
xtoe=xtoein;

toevel=.04;  %Toe velocity over basement [m/day]
seavel=25/16370/365; %Rate of sea-level change [m/day] 25m/16370years/365days
% = 4.1841e-6

%Create and store evolving margin and clinoform
%This done by, at each time step, calculating the bottom and top
%location of each numerical grid column (1001 cols), determined by the
%advance of the geometric delta

%simtime in [years] --set above
totstep=round(simtime*365/delt); %total number of time steps

for tstep=1:totstep+1
    sealevel=sealevelin+seavel*(tstep-1)*delt; %Change sea-level
    xtoe=xtoein+toevel*(tstep-1)*delt; %position of toe
    xshore=xtoe-(sealevel-basement)/foreslope; %position of shoreline
    %elevation of column tops at time tstep
    eta(tstep,:)=min ( max((xshore-xcols)*topslope+sealevel,sealevel), ...
        max((xtoe-xcols)*foreslope+basement,basement) );
    %elevation of column bottoms
    etabot(tstep,:)=zeros(1,cols);
    seastore(tstep)=sealevel;
    
end


%% Grid Initial domain

ytop=eta(1,:);  %initial elevation of top domain nodes
ybot=etabot(1,:); %initial elevation of bottom domain nodes
ncol=round((ytop-ybot)/Dely)+1; %number of nodes in each column
N=sum(round((ytop-ybot)/Dely)+1); %number of nodes in domain

%storage and initial values for head values
phi=max(ytop)*ones(N,1);
phi_insert=zeros(1,cols);

%storage and initial values for concentration
con_insert=zeros(1,cols);

%set initial salt concentrations
con=[];
for jj=1:cols
    if ytop(jj)>sealevelin;
        con=[con;zeros(ncol(jj),1)]; %Above sea-level C=0
    else
        con=[con;ones(ncol(jj),1)]; %Below sea-level C=1
    end
end



%% Time Stepping

tolh=1e-5; %convergence tolerance for head
tols=1e-7; %for concentration

printval=50;   %parameter to control plotting through time
printspace=50;


for tstep=1:totstep %Main time loop
    
    %% Adjustment of strat. and set sea-level
    
    %read in current top and bottom elevations
    if tstep<totstep+1
        ytop=eta(tstep,:);
        ybot=etabot(tstep,:);
    end
    
    %Interpolation for Inserted values
    %In the following we will check to see if a new node will be added to a
    %column. If this happens, we need to find a variable value of this node,
    %located Dely above the last but one node in the column, using
    %interpolation. Here we pre-calculate these values so they are available
    %when needed-- we do not need do this for the first time step
    
    if tstep>1
        for jj=1:cols
            mnode=Btop(1,jj); % node number of top node
            %interpolation ratio
            yrat=(y(mnode)-y(mnode-1)-Dely)/(y(mnode)-y(mnode-1));
            phi_insert(1,jj)=phi(mnode)-yrat*(phi(mnode)-phi(mnode-2));
            con_insert(1,jj)=con(mnode)-yrat*(con(mnode)-con(mnode-2));
        end
    end
    
    %number of nodes in each column
    ncolpre=ncol; %store previous value
    ncol=round((ytop-ybot)/Dely)+1; %calculate new value. Will only change if
    %top  node is at a distance > 1.5 Dely
    % above its neighbor below
    Btop=[]; %list of nodes on top surface
    Btopsea=[];  %list of top domain nodes under sea
    x=[]; %list of x locations of domain nodes
    y=[]; %list of y locations of domain nodes
    for jj=1:cols
        x=[x;Delx*(jj-1)*ones(ncol(jj),1)];
        y=[y;linspace(ybot(jj),ybot(jj)+(ncol(jj)-2)*Dely,ncol(jj)-1)';ytop(jj)];
        
        %Find Btop, nodes on top surface, Btopsea top surface nodes under sea
        Btop=[Btop,size(y,1)];
        if ytop(jj)<seastore(1,tstep)
            Btopsea=[Btopsea,size(y,1)];
        end
    end
    Bbot=[1,Btop(1:cols-1)+1]; %List of bottom nodes
    
    N=size(x,1); %updated number of node points
    
    %Resize and rebuild storage for head and concentration to account
    %for new nodes added to the grid.
    
    %store values from previous grid
    phipre=phi;
    conpre=con;
    
    
    %resize phi and con as grid expands
    con=[];
    phi=[];
    last=0;
    for jj=1:cols
        first=last+1;
        last=first+ncolpre(jj)-1;
        
        %No change in number of nodes in column
        if ncolpre(jj)==ncol(jj)
            phi=[phi;phipre(first:last)];
            con=[con;conpre(first:last)];
        end
        
        
        if ncolpre(jj)<ncol(jj)  %This condition indicates that
            % a node point has been added
            
            %We first assign the nodal values to new nodes
            %We copy the previous values to all nodes
            %And set the value at the added node to take the top node value
            phi=[phi;phipre(first:last);phinew(last)];
            con=[con;conpre(first:last);connew(last)];
            %We then overwrite the added node value---
            %using the stored interpolated values, calculated from
            %previous time step solution
            mnode=size(con,1); % gives top node position in current column
            phi(mnode-1)=phi_insert(1,jj);
            con(mnode-1)=con_insert(1,jj);
        end
        
    end
    
    
    %Set initial guess for new time step
    phinew=phi;
    connew=con;
    
    
    %% Make Mesh
    
    %We bring in two adjacent columns at a time and then use the
    %inbuilt MATLAB routine Delaunay to make the elements in the mesh
    t=[];
    first=1;
    for jj=1:cols-1
        last=first+ncol(jj)+ncol(jj+1)-1;
        xcc=x(first:last);
        ycc=y(first:last);
        tcc=delaunay(xcc,ycc)+first-1;
        t=[t;tcc];
        first=first+ncol(jj);
    end
    %The resulting array 't' has one line for each element in mesh.
    %this line contains a list of the three node indices
    %(in counter clockwise order) that are at the vertices on the element
    
    %NOTE the MATALB command
    %triplot(t,x,y, '-b');
    %will print out the grid at any time asked for
    
    Ntri=size(t,1);  %the size of t is the number of triangle elements
    
    %A COMMENT ON DATA STRUCTURE
    %The key data structure in CVFEM is a REGION OF SUPPORT
    %The region of support for a node p is  a list of all the nodes
    %that share an element side with node p--typically this will
    %be a list of the connected nodes in counter-clockwise order around p
    %Here we assemble the storage of the
    %region of support as follows. We have an array 'sup' with one row for
    %each node point. We sweep through all of our elements and each time
    %we find an element with node p as a vertex we sequentially add,
    %in a counter clockwise sense, the node numbers of the other two vertices
    %to row p of sup.
    %On completion of the sweep through all the elements row p will
    %contain all the information needed to identify the elements connected to
    %node p--each pair of entries in row p accounting for one element.
    
    %If there are np elements connected to node p then there will be
    %2*np entries on row p of array sup.
    %thus we size of sup as N X maxsup
    %N number of nodes in domain
    %maxsup=2Xmaxconnect (the maximum number of elements connected to a node)
    maxconnect=0;
    for ii=1:N
        maxconnect=max(maxconnect, size(find (t==ii),1));
        %NOTE size(find (t==ii),1)) finds number elements connected to node ii
    end
    maxsup=2*maxconnect; %double the nodes in support
    
    
    %END MESH XXXXXXXXXXXX
    
    %% Geometric properties of mesh
    %Here we calculate some geometric properties of the mesh
    %And assemble the support array sup
    
    Volp=zeros(N,1);    %CV volume
    xmid=zeros(Ntri,1); %elemnt mid point
    ymid=zeros(Ntri,1);
    Nx=zeros(Ntri,3);   %Derivatives of shape functions
    Ny=zeros(Ntri,3);
    %tsup=zeros(N,1); %triangles in support
    isup=ones(N,1);  %support index location
    sup=ones(N,maxsup); %support
    
    
    for itri=1:Ntri
        k1=t(itri,1); %global number of 1st node in triangle itri
        k2=t(itri,2); %2nd node
        k3=t(itri,3); %3rd node
        %element volume
        v=(x(k2)*y(k3)-x(k3)*y(k2)-x(k1)*y(k3)+x(k1)*y(k2) ...
            +y(k1)*x(k3)-y(k1)*x(k2))/2;
        
        
        Volp(k1)=Volp(k1)+v/3; %contribution to control volume
        Volp(k2)=Volp(k2)+v/3;
        Volp(k3)=Volp(k3)+v/3;
        
        xmid(itri)=(x(k1)+x(k2)+x(k3))/3; %mid point of element
        ymid(itri)=(y(k1)+y(k2)+y(k3))/3;
        
        %derivatives of shape functions
        Nx(itri,1)= (y(k2)-y(k3))/(2*v);   %Nx=shape function derivative
        Nx(itri,2)= (y(k3)-y(k1))/(2*v);   %the index 1, 2 or 3
        Nx(itri,3)= (y(k1)-y(k2))/(2*v);   %refers to local tri element node
        Ny(itri,1)=-(x(k2)-x(k3))/(2*v);   %Ny=shape function derivative
        Ny(itri,2)=-(x(k3)-x(k1))/(2*v);
        Ny(itri,3)=-(x(k1)-x(k2))/(2*v);
        
        
        %build support
        sup(k1,isup(k1))=k2;
        sup(k1,isup(k1)+1)=k3;
        isup(k1)=isup(k1)+2; %This sets position for next entry on row k1 of sup
        %tsup(k1)=tsup(k1)+1;
        
        sup(k2,isup(k2))=k3;
        sup(k2,isup(k2)+1)=k1;
        isup(k2)=isup(k2)+2;
        %tsup(k2)=tsup(k2)+1;
        
        sup(k3,isup(k3))=k1;
        sup(k3,isup(k3)+1)=k2;
        isup(k3)=isup(k3)+2;
        %tsup(k3)=tsup(k3)+1;
        
    end
    
    %% SET BOUNDARY CONDITIONS
    
    %BOUNDARY NODE POINTS
    
    %We are assuming that only no flow or prescribed boundary conditions
    %are applied.
    
    %By the manner in which our CVFEM discrete system is constructed no-flow
    %boundaries are the default setting.
    
    % The form of our discrete equation is give  in eq (10) V--is the Control volume Omega
    
    % (mu+bc_i)*V*u_i =mu*V*u^old_i + delt*[sum(a_j*u_j)-a_i*u_i+bb_i]
    
    %we rearrange to get
    %[(mu+ bc_i)*V+delt*a_i]*u_i =mu*Va^old_i + delt*[sum(a_j*u_j+bb_i]
    
    
    % To account for fixed value boundaries we add additional coefficients and terms
    %BCu_i ad BBu_i to get
    %[(mu+ bc_i)*V+delt*a_i +BCu_i]*u_i =mu*V*u^old_i + delt*[sum(a_j*u_j+bb_i]+BBu_i
    
    %The one-d arrays BCu and BBu are initialized with zeros
    %at nodes where a fixed value 'val' needs to be applied we set
    %BCu_i=1e18
    %BBu_i=1e18*val
    % ensuring, on solution, of the discrete equation that the correct value
    % is imposed at node i
    
    
    %Btop--stored node numbers on top boundary
    %Btopsea--stored top nodes at or below sea level.
    
    %Arrays for Boundary conditions
    %Recall
    %Btop--stored node numbers on top boundary
    %Btopsea--stored top nodes at or below sea level.
    
    Big=1e18;
    BCh=zeros(N,1); %boundary coeff values for head and solute
    BCs=zeros(N,1);
    BBh=zeros(N,1); %fixed head values
    BBs=zeros(N,1); %fixed solute concentration values
    
    %Fixed head value on top boundary
    BCh(Btop,1)=Big;
    BBh(Btop,1)=Big*(y(Btop));
    %correctiion for sea water
    BBh(Btopsea,1)=Big*(y(Btopsea)-(y(Btopsea)-seastore(1,tstep))*(rhorel));
    
    %concentration Fixed values 0 above sealevel 1 below
    BCs(Btop,1)=Big;
    BBs(Btopsea)=consea*Big;
    %Note for the submerged nodes stored in Btopsea we will
    %overwrite this condition when we detect outflow--see below
    
    %% Head Solution
    %Discrete Equation (see eq 10 in paper) has form
    %[(mu+ bc_i)*V+delt*a_i +BCu_i]*u_i =mu*V*u^old_i + delt*[sum(a_j*u_j+bb_i]+BB_i
    %Written as
    %[sto*V+delt*a_i +BCh_i]*phi_i =sto*V*phi^old_i + delt*[sum(a_j*phi_j+bb_i]+BBh_i
    %bc_i=0
    %bb_i -will be identified below--see variable BBvar (variable density source)
    %V=Omega=Volp
    
    %coefficents
    ap=zeros(N,1);        %coefficient at node
    asup=zeros(N,maxsup); %coefficients on nodes in support
    isup=ones(N,1);  %current point in support
    
    BBvar=zeros(N,1); %Variable density source THIS is bb_i is discrete equation
    
    kx=zeros(Ntri,1); %x direction conductivity value
    ky=zeros(Ntri,1); %y direction conductivity value
    
    for itri=1:Ntri
        
        %set elemnet conductivity
        kx(itri)=kxval;
        ky(itri)=kyval;
        
        %The element has three vertices nodes
        %Given local indices 1, 2 and 3
        %We will call at each node in the element in turn
        %and build the coefficients associated with that node
        %when we work with node 1 the neighbors are 2 and 3 (counter-clockwise)
        %when we work with node 2 the neighbors are 3 and 1 (counter-clockwise)
        %when we work with node 3 the neighbors are 1 and 2 (counter-clockwise)
        %We store this ordering in
        cyc=[1,2,3;2,3,1;3,1,2];
        
        for node=1:3
            ii=cyc(node,1);  %local node
            jj=cyc(node,2);
            kk=cyc(node,3);
            
            k1=t(itri,ii); %global node number
            k2=t(itri,jj);
            k3=t(itri,kk);
            
            Nx1=Nx(itri,ii);   %Nx=shape function derivative for element
            Nx2=Nx(itri,jj);
            Nx3=Nx(itri,kk);
            Ny1=Ny(itri,ii);   %Ny=shape function derivative
            Ny2=Ny(itri,jj);
            Ny3=Ny(itri,kk);
            
            
            %Each node in an element has two CV faces we determine the
            %coefficients by approximating flux across each face
            
            %Face1
            delx= (x(k1)+x(k2)+x(k3))/3-(x(k1)+x(k2))/2;
            dely= (y(k1)+y(k2)+y(k3))/3-(y(k1)+y(k2))/2;
            
            face1_k1=kx(itri)*Nx1*dely-ky(itri)*Ny1*delx;
            face1_k2=kx(itri)*Nx2*dely-ky(itri)*Ny2*delx;
            face1_k3=kx(itri)*Nx3*dely-ky(itri)*Ny3*delx;
            
            %variable density source(face value)
            %asssuming rel. density contribution is constant in element
            BBvar(k1)=BBvar(k1)-ky(itri)*((rhorel-1)/3)*(con(k1)+con(k2)+con(k3))*delx;
            
            
            %Face2
            delx= -(x(k1)+x(k2)+x(k3))/3+(x(k1)+x(k3))/2;
            dely= -(y(k1)+y(k2)+y(k3))/3+(y(k1)+y(k3))/2;
            
            
            face2_k1=kx(itri)*Nx1*dely-ky(itri)*Ny1*delx;
            face2_k2=kx(itri)*Nx2*dely-ky(itri)*Ny2*delx;
            face2_k3=kx(itri)*Nx3*dely-ky(itri)*Ny3*delx;
            
            
            ap(k1)=ap(k1)      -face1_k1-face2_k1;
            asup(k1,isup(k1))  =face1_k2+face2_k2;
            asup(k1,isup(k1)+1)=face1_k3+face2_k3;
            isup(k1)=isup(k1)+2;
            
            %varibale density source(face value)
            %asssuming rel. density contribution is constant in element
            BBvar(k1)=BBvar(k1)-ky(itri)*((rhorel-1)/3)*(con(k1)+con(k2)+con(k3))*delx;
        end
    end
    
    %Update Head
    %solution is for freshwater head
    %Uses Jacobi solver (works well with vectorization)
    %provides solution of equation 10 in text--written as
    %[sto*V+delt*a_i+BCh_i]*phinew_i =sto*V*phi_i + delt*[sum(a_j*phi_j+BBvar_i]+BBh_i
    %NOTE RHS=sum(asup.*phinew(sup),2)=sum(a_j*u_j)
    %bb_i=BBvar variable density source
    %V=Omega=Volp
    
    
    sto=stov; % sets storage value
    if tstep==1
        sto=0; % at tstep==1, to get initial conserved flow field storage=0
    end
    
    conver=1;
    phipre=phinew;
    
    while conver>tolh %convergence on iteration
        RHS=sum(asup.*phinew(sup),2);
        phinew=(sto*Volp.*phi+delt*(RHS+BBvar)+BBh)./(sto*Volp+BCh+ap*delt);
        conver=max(abs(phinew-phipre));
        phipre=phinew;
    end
    
    %calculate volume of flow stored per time step
    %need this term to account for storage in solute transport equation
    Vstore=sto*(Volp.*(phinew-phi));
    
    
    phi=phinew; %old for new
    
   
    
    %%  Solute solution
    %Discrete Equation (see eq 10 in paper) has form
    %[mu*V + bc_i +delt*a_i +BCu_i]*u_i =mu*V*u^old_i + delt*[sum(a_j*u_j+bb_i]+BB*u_i
    %Written as
    %[eps*V+ Vstore +delt*a_i +BCs_i]*connew_i =eps*V*con_i + delt*[sum(a_j*connew_j]+BBs_i
    %Vstore-- see above
    %bb_i=0
    %V=Omega=Volp
    
    %Coefficients
    
    ap=zeros(N,1);
    asup=zeros(N,maxsup);
    isup=ones(N,1);
    
    
    for itri=1:Ntri
        
        %Diffusion
        cyc=[1,2,3;2,3,1;3,1,2];
        for node=1:3
            ii=cyc(node,1);
            jj=cyc(node,2);
            kk=cyc(node,3);
            
            k1=t(itri,ii); %global node number of element vertices
            k2=t(itri,jj);
            k3=t(itri,kk);
            
            Nx1=Nx(itri,ii);   %Nx=shape function derivative
            Nx2=Nx(itri,jj);
            Nx3=Nx(itri,kk);
            Ny1=Ny(itri,ii);   %Ny=shape function derivative
            Ny2=Ny(itri,jj);
            Ny3=Ny(itri,kk);
            
            
            if node==1 %Dispersion for element
                
                % --contribution to discharge from fresh-water head
                qxval=-kx(itri)*(Nx1*phi(k1)+Nx2*phi(k2)+Nx3*phi(k3));
                qyval=-ky(itri)*(Ny1*phi(k1)+Ny2*phi(k2)+Ny3*phi(k3));
                
                %actual discharges at element midpoints
                %used to calculate dispersion tensor
                %also assumes rel. density contribution is constant in element
                qxmid=qxval;
                qymid=qyval-ky(itri)*((rhorel-1)/3)*(con(k1)+con(k2)+con(k3));
                
                qx(itri)=qxmid;
                qy(itri)=qymid;
                
                qx2=qxmid^2;
                qy2=qymid^2;
                qabs=sqrt(qx2+qy2);
                
                
                %Dispersion from Bear 1972 
                Dxx=aL*qx2/qabs+aT*qy2+Dmol*eps;
                Dyy=aT*qx2/qabs+aL*qy2+Dmol*eps;
                Dxy=(aL-aT)*(qxmid)*(qymid)/qabs;
                
                
            end
            
            %Face1
            %dischage across face 1--Note qyface included rel. density
            qxface=qxval;
            qyface=qyval-ky(itri)*((rhorel-1)/3)*(con(k1)+con(k2)+con(k3));
            delx= (x(k1)+x(k2)+x(k3))/3-(x(k1)+x(k2))/2;
            dely= (y(k1)+y(k2)+y(k3))/3-(y(k1)+y(k2))/2;
            
            %contribution to flux across face due to dispersion
            face1_k1=(Dxx*Nx1+Dxy*Ny1)*dely-(Dyy*Ny1+Dxy*Nx1)*delx;
            face1_k2=(Dxx*Nx2+Dxy*Ny2)*dely-(Dyy*Ny2+Dxy*Nx2)*delx;
            face1_k3=(Dxx*Nx3+Dxy*Ny3)*dely-(Dyy*Ny3+Dxy*Nx3)*delx;
            
            
            %contribution to flux across face due to flow
            %uses upwind, flow caries variable value in upstream direction
            qout=qxface*dely-qyface*delx; %flow out of vol k
            if qout>=0
                face1_k1=face1_k1-qout;
            else
                face1_k2=face1_k2-qout;
            end
            
            
            %Face2
            qxface=qxval;
            qyface=qyval-ky(itri)*((rhorel-1)/3)*(con(k1)+con(k2)+con(k3));
            
            delx= -(x(k1)+x(k2)+x(k3))/3+(x(k1)+x(k3))/2;
            dely= -(y(k1)+y(k2)+y(k3))/3+(y(k1)+y(k3))/2;
            
            face2_k1=(Dxx*Nx1+Dxy*Ny1)*dely-(Dyy*Ny1+Dxy*Nx1)*delx;
            face2_k2=(Dxx*Nx2+Dxy*Ny2)*dely-(Dyy*Ny2+Dxy*Nx2)*delx;
            face2_k3=(Dxx*Nx3+Dxy*Ny3)*dely-(Dyy*Ny3+Dxy*Nx3)*delx;
            
            %upwind
            qout=qxface*dely-qyface*delx; %flow out of vol k
            if qout>=0
                face2_k1=face2_k1-qout;
            else
                face2_k3=face2_k3-qout;
            end
            
            %update diagonal element and support
            
            ap(k1)=ap(k1)      -face1_k1-face2_k1;
            asup(k1,isup(k1))  =face1_k2+face2_k2;
            asup(k1,isup(k1)+1)=face1_k3+face2_k3;
            isup(k1)=isup(k1)+2;
            
        end
        
    end
    
    
    
    %Boundary Condition for outflows at submerged nodes
    %As a default we set the concentration at nodes on the
    %seafloor and foreset to  a prescribed constant value consea
    %this is achieved with the setting
    %BBs(Btopsea)=Big*consea;
    
    %But when there is outflow at a submerged surface node,
    %we reset the condition as follows:
    
    nsub=size(Btopsea,2); %submerged nodes
    %We can sweep over these nodes with the loop
    for ii=1:nsub
        %Obtaining for each ii the node index of the top surface node
        bnode=Btopsea(ii);
        %Recognizing the column storage nature of our grid
        %we known that the index of the node immediately below is bnode-1
        %Thus we can calculate the sign of the vertical discharge as
        qvertsig=-(phi(bnode,1)-phi(bnode-1,1))/(y(bnode,1)-y(bnode-1,1))...
                  -con(bnode-1)*(0.025);  
         if qvertsig>0 % we will have outflow and thus we set 
            BBs(bnode)=Big*con(bnode-1); 
            %i.e., with at outflow node we set the boundary concentration
            %equal to the value of the node immediately below; 
            %a Neumann condition 
            %Note this setting is explicit, lagging behind one time step
        end
    end
    
    
    
    %solution is for solute
    %Uses Jacobi solver (works well with vectorization)
    %provides solution of equation 10 in text--written as
    %[eps*V+ Vstore +delt*a_i +BCs_i]*connew_i =eps*V*con_i + delt*[sum(a_j*connew_j]+BBs_i
    %NOTE RHS=sum(asup.*connew(sup),2)=sum(a_j*connew_j)
    %V=Omega=Volp
    
    
    conver=1;
    conpre=connew;
    
    
    
    while conver>tols
        RHS=sum(asup.*connew(sup),2);
        connew=(eps*Volp.*con+delt*RHS+BBs)./(Vstore+eps*Volp+BCs+ap*delt);
        conver=max(abs(connew-conpre));
        conpre=connew;
    end
    
    
    con=connew;
    
    
    
    %Print out a flooded contour plot of concentrations at set and final times
    if tstep>printval || tstep >totstep-1 %will plot at last time step
        printval=printval+printspace;
        figure (1)
        %Plot outline of domain
        plot(xcols,ytop,'k-', 'LineWidth',2)
        ylim([0 450])
        hold on
        
        %Use trisurf MATLAB routine to print out
        %relative concentration field
        trisurf(t,x,y,con)
        map= [0 0 1
            .3 .3 1
            .5 .5 1
            0 1 0
            0 .5 0
            .7 .7 0
            1 .4 0
            1 0 0];
        colormap(map)
        shading interp
        view(0,90)
        colorbar
        daspect([157 1 1])
        xlabel('m')
        ylabel('m')
        
        
        hold off
        
        
        pause(.1)
        
    end
    
end