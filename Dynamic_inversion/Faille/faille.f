************************************************************************
* creation des points sources pour axitra
* entree au clavier ou par un fichier avec <
* sortie: source et axi.hist
************************************************************************

      program faille

      implicit none

      real t(100000),xs(100000),ys(100000),zs(100000)
      real pi
      real astk,adip,estk,edip
      real phi,delta,rake
      real slip
      real dstk,ddip
      real dxs,dxd,dys,dyd,dzd
      integer nstk,ndip,idip,istk,is,ns
      real zhypo
      real vr
      pi=3.14159265359

************************************************************************
* lecture de parametres d''entree
************************************************************************
      write(*,*)
      write(*,*)
      write(*,*) 'MAKe Sure The SOURCE file has been deleted first'
      write(*,*)
      write(*,*)
      open(55,form='formatted',file='faille.in')
       
      read(55,*) astk,adip
      write(6,*) 'longueur(stk) largeur(dip) en m        ',astk,adip
      write(*,*) astk,adip
      read(55,*) estk, edip
      write(6,*) 'position de l''epic. sur la faille en m ',estk,edip
      write(*,*)estk, edip
      read(55,*) zhypo
      write(6,*)"profondeur de l'hypocentre en m:         ", zhypo
      write(*,*) zhypo
      read(55,*) nstk,ndip
      write(6,*)"nb de sousfailles selon strk et dip:     ",nstk,ndip
      write(*,*) nstk,ndip
      read(55,*) phi
      write(6,*)"azim faille en degre par rapport au nord ",phi
      write(*,*)phi
      read(55,*) delta
      write(6,*)"pendage (dip) en degre:                  ",delta
      write(*,*) delta
      read(55,*) rake
      write(6,*)"rake en degre:                           ",rake
      write(*,*) rake
      read(55,*) slip
      write(6,*)"amplitude du glissement (m) ?            ", slip
      write(*,*) slip
      read(55,*) vr
      write(6,*)"vitesse de rupture en m/s:               ", vr
      write(*,*)vr
      write(*,*)

      close(55)
************************************************************************
* ouverture des fichiers de sortie
************************************************************************

      open(10,form='formatted',file='source')
      open(11,form='formatted',file='axi.hist')

************************************************************************
* nombre de points source (indexes de 1 a ns)
************************************************************************

      ns=nstk*ndip

************************************************************************
* angles en radian
************************************************************************

      phi=phi*pi/180.
      delta=delta*pi/180.

************************************************************************
* taille de chaque sous faille selon le strike et le dip
************************************************************************

      dstk=astk/nstk
      ddip=adip/ndip
      write(6,*)"dstk, dsty = ",dstk,ddip

************************************************************************
* decalage selon x,y ou z provoque entre sous-failles contigues
************************************************************************

      dxs=dstk*cos(phi)
      dxd=-ddip*cos(delta)*sin(phi)
      dys=dstk*sin(phi)
      dyd=ddip*cos(delta)*cos(phi)
      dzd=ddip*sin(delta)
c      write(6,*)"5"

************************************************************************
* coordonnees de la premiere sous-faille dans le repere x,y et z
************************************************************************

      xs(1)=(dxs+dxd)/2-estk*cos(phi)+edip*cos(delta)*sin(phi)
      ys(1)=(dys+dyd)/2-estk*sin(phi)-edip*cos(delta)*cos(phi)
      zs(1)=dzd/2-edip*sin(delta)+zhypo
      write(6,*)"first = (",xs(1), ys(1), zs(1)," )"

************************************************************************
* coordonnees des autres sous-failles et ecriture dans les fichiers
************************************************************************

      is=1
      do idip=1,ndip
        do istk=1,nstk
c        write(6,*)"1"
          xs(is)=xs(1)+(istk-1)*dxs+(idip-1)*dxd
          ys(is)=ys(1)+(istk-1)*dys+(idip-1)*dyd
          zs(is)=zs(1)+(idip-1)*dzd
          t(is)=(((xs(is)**2)+(ys(is)**2)+((zs(is)-zhypo)**2))**0.5)/vr
          write(10,*) is,xs(is),ys(is),zs(is)
          t(is)=t(is)
          write(11,'(I5,7G12.3)')is,slip,phi*180/pi,delta*180/pi,rake,
     $        ddip,dstk,t(is)
          is=is+1
       enddo
      enddo
      write(6,*)"last = (",xs(is-1), ys(is-1), zs(is-1)," )"

      close(10)
      close(11)

************************************************************************

      end

