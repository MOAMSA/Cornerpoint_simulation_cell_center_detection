%%% This code findes the center points of a corner points structures in the
%%% 3D space. The algorithm of this code is as below:
%%% Step 1: Extracting the top and bottom coordinate of pillars and finding
%%% the pillar lines equation.
%%% Step 2: Extracting the Z-coordinates of grid corners and finding the X 
%%% and Y-coordinate of these points using the equations and Z-coordinate.
%%% Step 3: Finding the center points of top and bottom face of each cell
%%% by intersecting the digonal lines.
%%% Step 4: Finding the cell center using arithmetic average of two last
%%% points.
%
% Dimension of i,j and k-direction.
%
Nx = 138;
Ny = 106 ;
Nz = 19 ;
%
% Loading COORD file (contains pillar  coordinates and Z-coordinates of
% grid corners).
%
load COORD ;
%
% To remove negative Z-coordinates.
%
COORD_F = abs(COORD_F) ;
%
% Initialization output matrix.
%
Pcl = zeros(Nx*Ny*Nz,3) ;
Pgx = zeros(Nx,Ny,Nz,8) ;
Pgy = zeros(Nx,Ny,Nz,8) ;
Pgz = zeros(Nx,Ny,Nz,8) ;
%
% Constructing the center points of each cells.
%
for kk = 1:Nz
    for jj = 1:Ny
        for ii = 1:Nx
            %
            % Finding the bottom and points of near left pillars.
            %
            Pt = COORD_F((jj-1)*(Nx+1)+ii,1:3) ;
            Pb = COORD_F((jj-1)*(Nx+1)+ii,4:6) ;
            %
%             Znlt = ZCORN_F(2*(jj-1)*(Nx+1)+2*ii-1) ;
            % Finding the Z-coordinate of near left top point of each cell.
            %
            Znlt = ZCORN_F(4*(jj-1)*Nx+2*ii-1+8*(kk-1)*Nx*Ny) ;
            %
            % Finding the X and Y-coordinate of near left top point of each
            % cell using the Z-coordinate and line equation of near left
            % pillar.
            %
            t = (Znlt-Pb(3))/(Pt(3)-Pb(3)) ;
            Xnlt = (Pt(1)-Pb(1))*t+Pb(1) ;
            Ynlt = (Pt(2)-Pb(2))*t+Pb(2) ;
            %
            % [X Y Z] of near left top point.
            %
            Pnlt = [Xnlt Ynlt Znlt] ;
%             Pgx (ii,jj,kk,1) = Xnlt ;
%             Pgy (ii,jj,kk,1) = Ynlt ;
%             Pgz (ii,jj,kk,1) = Znlt ;
            %
            % Finding the bottom and points of near right pillars.
            %
            Pt = COORD_F((jj-1)*(Nx+1)+ii+1,1:3) ;
            Pb = COORD_F((jj-1)*(Nx+1)+ii+1,4:6) ;
%             Znrt = ZCORN_F(2*(jj-1)*(Nx+1)+2*ii) ;
            %
            % Finding the Z-coordinate of near right top point of each cell.
            %
            Znrt = ZCORN_F(4*(jj-1)*Nx+2*ii+8*(kk-1)*Nx*Ny) ;
            %
            % Finding the X and Y-coordinate of near right top point of each
            % cell using the Z-coordinate and line equation of near left
            % pillar.
            %
            t = (Znrt-Pb(3))/(Pt(3)-Pb(3)) ;
            Xnrt = (Pt(1)-Pb(1))*t+Pb(1) ;
            Ynrt = (Pt(2)-Pb(2))*t+Pb(2) ;
            %
            % [X Y Z] of near rigth top point.
            %
            Pnrt = [Xnrt Ynrt Znrt] ;
%             Pgx (ii,jj,kk,2) = Xnlt ;
%             Pgy (ii,jj,kk,2) = Ynlt ;
%             Pgz (ii,jj,kk,2) = Znlt ;
            %
            % Finding the bottom and points of far left pillars.
            %
            Pt = COORD_F(jj*(Nx+1)+ii,1:3) ;
            Pb = COORD_F(jj*(Nx+1)+ii,4:6) ;
            %
%             Zflt = ZCORN_F(2*jj*(Nx+1)+2*ii-1) ;
            % Finding the Z-coordinate of far left top point of each cell.
            %
            Zflt = ZCORN_F(6*(jj-1)*Nx+2*ii-1+8*(kk-1)*Nx*Ny) ;
            %
            % Finding the X and Y-coordinate of far left top point of each
            % cell using the Z-coordinate and line equation of near left
            % pillar.
            %
            t = (Zflt-Pb(3))/(Pt(3)-Pb(3)) ;
            Xflt = (Pt(1)-Pb(1))*t+Pb(1) ;
            Yflt = (Pt(2)-Pb(2))*t+Pb(2) ;
            %
            % [X Y Z] of far left top point.
            %
            Pflt = [Xflt Yflt Zflt] ;    
%             Pgx (ii,jj,kk,3) = Xnlt ;
%             Pgy (ii,jj,kk,3) = Ynlt ;
%             Pgz (ii,jj,kk,3) = Znlt ;
            %
            % Finding the bottom and points of far right pillars.
            %
            Pt = COORD_F(jj*(Nx+1)+ii+1,1:3) ;
            Pb = COORD_F(jj*(Nx+1)+ii+1,4:6) ;
%             Zfrt = ZCORN_F(2*jj*(Nx+1)+2*ii) ;
            %
            % Finding the Z-coordinate of far right top point of each cell.
            %
            Zfrt = ZCORN_F(6*(jj-1)*Nx+2*ii+8*(kk-1)*Nx*Ny) ;
            %
            % Finding the X and Y-coordinate of far right top point of each
            % cell using the Z-coordinate and line equation of near left
            % pillar.
            %
            t = (Zfrt-Pb(3))/(Pt(3)-Pb(3)) ;
            Xfrt = (Pt(1)-Pb(1))*t+Pb(1) ;
            Yfrt = (Pt(2)-Pb(2))*t+Pb(2) ;
            %
            % [X Y Z] of far right top point.
            %
            Pfrt = [Xfrt Yfrt Zfrt] ;
%             Pgx (ii,jj,kk,4) = Xnlt ;
%             Pgy (ii,jj,kk,4) = Ynlt ;
%             Pgz (ii,jj,kk,4) = Znlt ;
            % Top points of each cells.
            P1 = Pnlt ;
            P2 = Pnrt ;
            P3 = Pfrt ;
            P4 = Pflt ;
            %
            % Intersecting digonal lines of top face of cells to find the
            % mid point of top face.
            %
            t1 = (P1(1)*(P2(2)-P3(2))+P2(1)*(P3(2)-P1(2))+...
                P3(1)*(P1(2)-P2(2)))/((P1(2)-P3(2))*(P4(1)-...
                P2(1))+(P4(2)-P2(2))*(P3(1)-P1(1))) ;
             t2 = (P1(1)*(P2(3)-P3(3))+P2(1)*(P3(3)-P1(3))+P3(1)*(P1(3)-P2(3)))/((P1(3)-P3(3))*(P4(1)-P2(1))+(P4(3)-P2(3))*(P3(1)-P1(1))) ;
             t3 = (P1(3)*(P2(2)-P3(2))+P2(3)*(P3(2)-P1(2))+P3(3)*(P1(2)-P2(2)))/((P1(2)-P3(2))*(P4(3)-P2(3))+(P4(2)-P2(2))*(P3(3)-P1(3))) ;
             t1 = min(abs([t1 t2 t3])-0.5)+0.5 ;
             
             
            %
            %  Mid point of top face of each cell.
            %
           Pct(1,1) = (P4(1)-P2(1))*t1+P2(1) ;
           Pct(1,2) = (P4(2)-P2(2))*t1+P2(2) ;
           Pct(1,3) = (P4(3)-P2(3))*t1+P2(3) ;
           %
           % Finding the bottom and points of near left pillars.
           %
           Pt = COORD_F((jj-1)*(Nx+1)+ii,1:3) ;
           Pb = COORD_F((jj-1)*(Nx+1)+ii,4:6) ;
           %Znlb = ZCORN_F(2*(jj-1)*(Nx+1)+2*ii-1+4*Nx*Ny*Nz) ;
           Znlb = ZCORN_F(4*(jj-1)*Nx+2*ii-1+4*Nx*Ny+8*(kk-1)*Nx*Ny) ;
           t = (Znlb-Pb(3))/(Pb(3)-Pt(3)) ;
           Xnlb = (Pb(1)-Pt(1))*t+Pb(1) ;
           Ynlb = (Pb(2)-Pt(2))*t+Pb(2) ;
           Pnlb = [Xnlb Ynlb Znlb] ;

           Pgx (ii,jj,kk,5) = Xnlt ;
           Pgy (ii,jj,kk,5) = Ynlt ;
           Pgz (ii,jj,kk,5) = Znlt ;
           %
           % Finding the bottom and points of near rigth pillars.
           %
            Pt = COORD_F((jj-1)*(Nx+1)+ii+1,1:3) ;
            Pb = COORD_F((jj-1)*(Nx+1)+ii+1,4:6) ;
            %Znrb = ZCORN_F(2*(jj-1)*(Nx+1)+2*ii+4*Nx*Ny*Nz) ;
            Znrb = ZCORN_F(4*(jj-1)*Nx+2*ii+4*Nx*Ny+8*(kk-1)*Nx*Ny) ;
            t = (Znrb-Pb(3))/(Pb(3)-Pt(3)) ;
            Xnrb = (Pb(1)-Pt(1))*t+Pb(1) ;
            Ynrb = (Pb(2)-Pt(2))*t+Pb(2) ;
            Pnrb = [Xnrb Ynrb Znrb] ;

           Pgx (ii,jj,kk,6) = Xnlt ;
           Pgy (ii,jj,kk,6) = Ynlt ;
           Pgz (ii,jj,kk,6) = Znlt ;
            %
            % Finding the bottom and points of far left pillars.
            %
            Pt = COORD_F(jj*(Nx+1)+ii,1:3) ;
            Pb = COORD_F(jj*(Nx+1)+ii,4:6) ;
            %Zflb = ZCORN_F(2*jj*(Nx+1)+2*ii-1+4*Nx*Ny*Nz) ;
            Zflb = ZCORN_F(6*(jj-1)*Nx+2*ii-1+4*Nx*Ny+8*(kk-1)*Nx*Ny) ;
            t = (Zflb-Pb(3))/(Pb(3)-Pt(3)) ;
            Xflb = (Pb(1)-Pt(1))*t+Pb(1) ;
            Yflb = (Pb(2)-Pt(2))*t+Pb(2) ;
            Pflb = [Xflb Yflb Zflb] ;    

           Pgx (ii,jj,kk,7) = Xnlt ;
           Pgy (ii,jj,kk,7) = Ynlt ;
           Pgz (ii,jj,kk,7) = Znlt ;
            %
            % Finding the bottom and points of far rigth pillars.
            %
            Pt = COORD_F(jj*(Nx+1)+ii+1,1:3) ;
            Pb = COORD_F(jj*(Nx+1)+ii+1,4:6) ;
            Zfrb = ZCORN_F(6*(jj-1)*Nx+2*ii+4*Nx*Ny+8*(kk-1)*Nx*Ny) ;
            %Zfrb = ZCORN_F(2*jj*(Nx+1)+2*ii+4*Nx*Ny*Nz) ;
            t = (Zfrb-Pb(3))/(Pb(3)-Pt(3)) ;
            Xfrb = (Pb(1)-Pt(1))*t+Pb(1) ;
            Yfrb = (Pb(2)-Pt(2))*t+Pb(2) ;
            Pfrb = [Xfrb Yfrb Zfrb] ;

           Pgx (ii,jj,kk,8) = Xnlt ;
           Pgy (ii,jj,kk,8) = Ynlt ;
           Pgz (ii,jj,kk,8) = Znlt ;

            P1 = Pnlb ;
            P2 = Pnrb ;
            P3 = Pfrb ;
            P4 = Pflb ;

            t1 = (P1(1)*(P2(2)-P3(2))+P2(1)*(P3(2)-P1(2))+P3(1)*(P1(2)-P2(2)))/((P1(2)-P3(2))*(P4(1)-P2(1))+(P4(2)-P2(2))*(P3(1)-P1(1))) ;
            t2 = (P1(1)*(P2(3)-P3(3))+P2(1)*(P3(3)-P1(3))+P3(1)*(P1(3)-P2(3)))/((P1(3)-P3(3))*(P4(1)-P2(1))+(P4(3)-P2(3))*(P3(1)-P1(1))) ;
            t3 = (P1(3)*(P2(2)-P3(2))+P2(3)*(P3(2)-P1(2))+P3(3)*(P1(2)-P2(2)))/((P1(2)-P3(2))*(P4(3)-P2(3))+(P4(2)-P2(2))*(P3(3)-P1(3))) ;
              t1 = min(abs([t1 t2 t3])-0.5)+0.5 ;

           Pcb(1,1) = (P4(1)-P2(1))*t1+P2(1) ;
           Pcb(1,2) = (P4(2)-P2(2))*1+P2(2) ;
           Pcb(1,3) = (P4(3)-P2(3))*t1+P2(3) ;

           Pcl((kk-1)*Ny+(jj-1)*Nx+ii,1) = (Pct(1,1)+Pcb(1,1))/2 ;
           Pcl((kk-1)*Ny+(jj-1)*Nx+ii,2) = (Pct(1,2)+Pcb(1,2))/2 ;
           Pcl((kk-1)*Ny+(jj-1)*Nx+ii,3) = (Pct(1,3)+Pcb(1,3))/2 ;

%            figure(1) ;
                hold on ;
%            fill3([Pnlt(1) Pnrt(1) Pfrt(1) Pflt(1) Pnlt(1)],[Pnlt(2) Pnrt(2) Pfrt(2) Pflt(2) Pnlt(2)] ,-[Pnlt(3) Pnrt(3) Pfrt(3) Pflt(3) Pnlt(3)],'r')  ;
%            fill3([Pnlb(1) Pnrb(1) Pfrb(1) Pflb(1) Pnlb(1)],[Pnlb(2) Pnrb(2) Pfrb(2) Pflb(2) Pnlb(2)] ,-[Pnlb(3) Pnrb(3) Pfrb(3) Pflb(3) Pnlb(3)],'r')  ;
%            fill3([Pnlt(1) Pnlb(1) Pnrb(1) Pnrt(1)],[Pnlt(2) Pnlb(2) Pnrb(2) Pnrt(2)],-[Pnlt(3) Pnlb(3) Pnrb(3) Pnrt(3)],'r') ;
%            fill3([Pflt(1) Pflb(1) Pfrb(1) Pfrt(1)],[Pflt(2) Pflb(2) Pfrb(2) Pfrt(2)],-[Pflt(3) Pflb(3) Pfrb(3) Pfrt(3)],'r') ;
%            fill3([Pflt(1) Pflb(1) Pnlb(1) Pnlt(1)],[Pflt(2) Pflb(2) Pnlb(2) Pnlt(2)],-[Pflt(3) Pflb(3) Pnlb(3) Pnlt(3)],'r') ;
%            fill3([Pfrt(1) Pfrb(1) Pnrb(1) Pnrt(1)],[Pfrt(2) Pfrb(2) Pnrb(2) Pnrt(2)],-[Pfrt(3) Pfrb(3) Pnrb(3) Pnrt(3)],'r') ;
        end
    end
end
plot3(Pcl(1:(kk-1)*Ny+(jj-1)*Nx+ii,1),Pcl(1:(kk-1)*Ny+(jj-1)*Nx+ii,2),-Pcl(1:(kk-1)*Ny+(jj-1)*Nx+ii,3),'.') ;




