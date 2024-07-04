cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                MODEL OF THE SCATTERING OF TWO PARTICLES (Deuteron, two glueballs...)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                THIS CODE IS INTENTIONALLY INCOMPLETE. 
C                In order to publish it, I decided not to provide the full code,
c                since it can be still used to do research.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                Trotti Enrico
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      
      
      implicit real*8(a-h,o-z)
      integer nmax,l
      parameter(nbasis=8)
c        
      parameter(nmax=nbasis)
      parameter(idim=nbasis)
      real*8 nnn(nmax,nmax),tnn(nmax,nmax),vnn(nmax,nmax),nun
     &,hnn(nmax,nmax),ener(nmax),vec(nmax,nmax),const,rlim,rlim2,cmg,vlim,vlim2
      common /ar1/a,r1
      common /param/cmr

 111  format(10(e9.3,2x))

c input parameters
      cmn=939d0 ! neutron mass (unit: MeV)
      cmp=938d0 ! proton mass (unit: MeV)
      cmr=cmn*cmp/(cmn+cmp) ! reduced mass (unit: MeV)

      pi=3.14159265358979d0
      cmg=1653d0
      do iconst=1,450
      const=50d0 + 1d0*dble(iconst)
      do irlim=1,100
      rlim=1d-2 + 1d-3*dble(irlim)
      
      vlim=(-const*dexp(-cmg*rlim/197d0))*197d0/(rlim*4d0*pi)      
      rlim2=1.5d-1
      vlim2=(-const*dexp(-cmg*rlim2/197d0))*197d0/(rlim2*4d0*pi)
      r1=0.005d0    ! unit : fm
      rmax=1.d0  ! unit : fm
      a  = (rmax/r1)**(1.d0/dble(nbasis-1))

      call mat(nmax,nnn,hnn,tnn,cmg,rlim,rlim2,const,vlim,vlim2)

      call doublediag(nmax,nnn,hnn,ener,vec)

c      write(*,*)'hamiltonian eigenvalues'
      write(*,*)'const , rlim , ener1 , ener2'
      write(*,111) const,rlim,ener(1),ener(2)
      write(*,*)
c I want the binding energy 1 to be mG-2mG
C and the energy 2 mG(0++*)-2mG      
      if((ener(1)+1730d0)**2.lt.4d0.and.(ener(2)+790d0)**2.lt.4d0)then
      write(*,*) 'this is ok'
      write(*,111) const,rlim,ener(1),ener(2)
      stop     
      end if

      end do
      end do

      if((ener(1)+1730)**2.gt.1.and.(ener(2)+790)**2.gt.1)then
      write(*,*) 'sorry, try again' 
      end if

c      if(ener(1).lt.emem)then
c         emem=ener(1)
c      endif


c      enddo


      open(919,file='a.txt',status='replace')


      write(919,*)'radius','     potential'     
      dr = 1.d-3  ! unit : fm
      do i=1,120000
         r=dble(i)*dr
      write(919,*)r
     & ,yukawa(r,cmg,rlim,rlim2,const,vlim,vlim2)
     & ,r*r*wavfct0(r/197.4d0,vec,nmax)**2/197.4d0**3*dr  
     & ,ener(2) 

      enddo

      close(919)


      t=0.d0
      do i=1,10000
         t=t
     &   +r*r*wavfct0(r,vec,nmax)**2*dr
      end do
      write(*,*)'All probability distribution',t


      stop
      end

ccccccccccccccccccccccccccc
c Can be used only in this program!
      real*8 function wavfctall(r,vec,nmax)
      implicit real*8(a-h,o-z)
      real*8 vec(nmax,nmax)
      pi=3.14159265358979d0

      wavfctall=0.d0
      do i=1,nmax
         wavfctall=wavfctall
     &     +vec(i,1)*basis(r,i,0)
      end do

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function retta(r,cmg,rlim,rlim2,const,vlim,vlim2)
      implicit real*8(a-h,o-z)
      pi=3.14159265358979d0
      retta=(-const*dexp(-cmg*rlim/197d0))*197d0/(rlim*4d0*pi)

      return
      end  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function yukawa(r,cmg,rlim,rlim2,const,vlim,vlim2)
      implicit real*8(a-h,o-z)
      pi=3.14159265358979d0
      yukawa=(-const*dexp(-cmg*r/197d0))*197d0/(r*4d0*pi) 

      return
      end  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function potential(r,cmg,rlim,rlim2,const,vlim,vlim2)
      implicit real*8(a-h,o-z)
      pi=3.14159265358979d0
      if(r.le.rlim)then
      potential=retta(r,cmg,rlim,rlim2,const,vlim,vlim2)
c      elseif(r.le.rlim2)then
c      potential=((-const*dexp(-cmg*rlim/197d0))*197d0/(rlim*4d0*pi))+
c     .    ((r-rlim)*(vlim2-vlim)/(rlim2-rlim))     
      else
      potential=yukawa(r,cmg,rlim,rlim2,const,vlim,vlim2)
      endif
      return
      end   
ccccccccccccccccccccccccccc
ccccccccccccccccccccccccccc
      real*8 function basis(r,n,l)
      implicit real*8 (a-h,o-z)
      real*8 nun

      pi=3.14159265358979d0

      lfac=1
      if(l.eq.0)then
         lfac=1
      else
  

      return
      end

ccccccccccccccccccccccccccccccccccccc
c
c Hamiltonian matrix elements are given
c
      subroutine mat(nmax,nnn,hnn,tnn,cmg,rlim,rlim2,const,vlim,vlim2)
      implicit real*8 (a-h,o-z)
      real*8 nnn(nmax,nmax),tnn(nmax,nmax)

     & ,ayuk(10)
      common /ar1/a,r1
      common /param/cmr

      pi=3.14159265358979d0


      n3=nmax


      redmn= cmr   ! reduced mass (unit: MeV)


c Norm matrix

      do i=1,nmax
         do j=1,nmax
            nnn(i,j)=0.d0
         end do
      end do

      l=0
      do 
         do 
            nnn(i,j)=

         end do
      end do


c Kinetic matrix

      do 
         do 
         end do
      end do

      l=0
      do 
         do 

         end do
      end do



c Coulomb force (for test)

      do i=1,nmax
         do j=1,nmax
            vdeut(i,j)=0.d0
         end do
      end do

      do i=1,n3
         do j=1,n3

!---- Numerical integral start
           mesh=100    ! number of integral mesh (how much integr discretized)
           temp=0.d0
c  dr is the spacing of integration (in fm)... here 1/100=0.01
           dr=1.d0/dble(mesh)
            do int=1,mesh *10
              r= dr*dble(int)

            enddo

!---- Numerical integral end

        vdeut(i,j)= dr *temp
     &       *4.d0*(4.d0*nun(i)*nun(j))**(3.d0/4.d0)/dsqrt(pi)

         end do
      end do




c Hamiltonian

      do i=1,nmax
         do j=1,nmax
          hnn(i,j)=tnn(i,j)+vdeut(i,j)
         end do
      end do
      

      return
      end
ccccccccccccccccccccccccccc
      real*8 function raku(i,j,l,mu)
      implicit real*8 (a-h,o-z)
      real*8 nun,mu
      common /ar1/a,r1

      raku=
     &     (2.d0*dsqrt(nun(i)*nun(j))
     &     /(nun(i)+nun(j)+mu) )**(dble(l)+1.5d0)

      return
      end
ccccccccccccccccccccccccccc
      real*8 function nun(n)
      implicit real*8 (a-h,o-z)


      return
      end

ccccccccccccccccccccccccccccccc
      real*8 function fac(l)

      fac=1.d0
      if(l.eq.0)then
         fac=1.d0
      else
         do i=1,l
            fac=dble(i)*fac
         end do
      endif

      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine doublediag(nmax,nnn,hnn,ener,vec)
      real*8 nnn(nmax,nmax),hnn(nmax,nmax),ener(nmax)
     &  ,vec(nmax,nmax),x(nmax,nmax),d(nmax)
     &  ,htbd(nmax,nmax)
     &  ,test(nmax,nmax)


      do i=1,nmax
         do j=1,nmax
            x(i,j)=nnn(i,j)
         end do
      end do

c      do i=1,nmax
c         write(*,111) (nnn(i,j),j=1,nmax)
c      end do
c      write(*,*)


      call realsymmatrixdiag(x,nmax,d)




c
c if norm matrix minus eigenvalues exist, then stop
c

      do i=1,nmax
         if(d(i).lt.0.d0)then
            write(*,*)'minus eigenvalue in norm matrix!'
            write(*,*)'Operation stopped'
            stop
         endif
      end do





      do i=1,nmax
         do j=1,nmax
            temp=0.d0
            do k1=1,nmax
               do k2=1,nmax
      temp=temp
     & +x(k1,i)*hnn(k1,k2)*x(k2,j)/(dsqrt(d(i)*d(j)))
               end do
            end do
            htbd(i,j)=temp
         end do
      end do


      call realsymmatrixdiag(htbd,nmax,ener)



      do i=1,nmax
         do j=1,nmax
            temp=0.d0
            do k=1,nmax
               temp=
            end do
            vec(i,j)=temp
         end do
      end do



      return
      end


c==================================================================
c
c     hermitian matrix diagonalization subroutine:
c
c     a  : input real symmetric matrix (becomes rotation matrix in output)
c     d : output eigenvalues 
c
      subroutine realsymmatrixdiag(a,nc,d)
      implicit real*8(a-h,o-z)
c      implicit complex*16(c)
c      complex*16 c(nc,nc),cz(nc,nc),czmin(nc)
      real*8 a(nc,nc),d(nc),e(nc),dc(nc),amin(nc)
      integer njogai(nc)


c     write the real symmetric matrix
c
c      do i=1,nc
c         write(*,*) (a(i,j),j=1,nc)
c      end do



c     first step: make the tridiagonal matrix

      do i=nc,2,-1
         l=i-1
         h=0.d0
         scale=0.d0
         if(l.gt.1) then
            do k=1,l
               scale= scale+ abs(a(i,k))
            end do
            if(scale.eq.0.d0) then
               e(i)=a(i,1)
            else
               do 
               f=0.d0
               do j=1,l
c     omit following line if finding only eigenvalues
                  a(j,i)=a(i,j)/h
                  g=0.d0
                  do k=1,j
                     g=g+a(j,k)*a(i,k)
                  end do
                  do k=j+1,l
                     g=g+a(k,j)*a(i,k)
                  end do
                  e(j)=g/h
                  f=f+e(j)*a(i,j)
               end do
               hh=f/(h+h)
               do j=1,l
                  
                  end do
               end do
            endif
            
         else
            e(i)=a(i,l)
         endif
         d(i)=h
      end do
      
c     Omit following line if finding only eigenvalues
      d(1)=0.d0
      e(1)=0.d0
      
      do i=1,nc
c     Delete lines from here ...
         l=i-1
         if(d(i).ne.0.d0)then
            do j=1,l
               g=0.d0
               do k=1,l
                  g=g+a(i,k)*a(k,j)
               end do
               do k=1,l
                  a(k,j)=a(k,j)-g*a(k,i)
               end do
            end do
         endif
c     .... to here when finding only eigenvaalues.
         d(i)=a(i,i)
c     also delete lines from here ....
         a(i,i)=1.d0
         do j=1,l
            a(i,j)=0.d0
            a(j,i)=0.d0
         end do
c     .... to here when finding only eigenvalues.
      end do

c     second step: diagonalize iteratively the tridiagonal matrix, take its rotation matrix

c     diagonalization

      do i=2,nc
         e(i-1)=e(i)
      end do
      e(nc)=0.d0
      do l=1,nc
         iter=0
 1              
         m=nc
 2              if(m.ne.l)then
            if(iter.eq.300)then
               write(*,*) 'too many iterations in tqli'
               stop
            endif
            iter=iter+1

c            write(*,*)'iter=',iter
c            write(*,*)'e(l)=',(e(k),k=1,nc)

            g=(d(l+1)-d(l))/(2.*e(l))
            r=pythag(g,1.d0)
            g=d(m)-d(l)+e(l)/(g+sign(r,g))
            s=1.d0
            bc=1.d0
            p=0.d0
            do i=m-1,l,-1
               f=s*e(i)
                  e(m)=0.d0
                  goto 1
               endif
               s=f/r
               bc=g/r
               g=d(i+1)-p
               r=(d(i)-g)*s+2.d0*bc*b
               p=s*r
               d(i+1)=g+p
               g=bc*r-b
c     Omit lines from here ....
               do k=1,nc
                  f=a(k,i+1)
                  a(k,i)=bc*a(k,i)-s*f
               end do
c     ....to here when finding only eigenvalues.
            end do
            d(l)=d(l)-p
            e(l)=g
            e(m)=0.d0
            goto 1
         endif
      end do

c     last step: sort eigenvalues in ascending order

      do i=1,nc-1
         dmin=d(i)
         mini=i
         do j=i+1,nc
            if(d(j).lt.dmin)then
               dmin=d(j)
               do k=1,nc
                  amin(k)=a(k,j)
               end do
               mini=j
            endif
               a(k,i)=amin(k)
            end do
         endif
      end do


      return
      end

c================================================================

      real*8 function pythag(a1,b1)
c      implicit none
      real*8 a1,b1
      real*8 absa,absb

      absa=abs(a1)
         pythag=absa*sqrt(1.d0+(absb/absa)**2)
      else
         if(absb.eq.0.d0)then
            pythag=0.d0
         else
            pythag=absb*sqrt(1.d0+(absa/absb)**2)
         endif
      endif
      return
      end
