Module mod_fit
implicit none

contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine optimize_param(x,y,sigma,n,parameters,m,converged,chisq,write_files)
  !! Input: x(n),y(n),sigma(n),parameters(m) (Initial Guess),write_files(1: generates files output,data.out)
  !! Input: subroutine funcs
  !! Output: parameters(m),converged(1:converged,0:unconverged),chisq
  !! Look for parameters tolerance and iter_max if having convergence issues
  implicit none
  integer n,m,converged,setup,lwork,iter_max,i,write_files
  real*8 x(n),y(n),sigma(n),sigma_inv(n),error(n)
  real*8 parameters(m),deriv(m)
  real*8 jacobian(n,m)
  real*8 lambda,chisq,tolerance,tol2,f
  real*8 saved_parameters(m)
  real*8,allocatable::work(:)
  character*3 ch1   !! Assumes less than 1000 parameters
  character*8 str
  integer j,k
  real*8 a(m), chisq_tmp, a_tmp(m)
  !! IMPORTANT PARAMETERS
  !---------------------------------------------------------- 
  tolerance=1.d-45  !! tolerance of chisq for convergence
  !tolerance=0.5d-52 !! tolerance of chisq for convergence
  !tol2=0.5d-42 ! tolerance of chisq for convergence
  tol2=1.d-45 !! tolerance of chisq for convergence
  !tol2=0.5d-52 !! tolerance of chisq for convergence
  iter_max=6000     !! Maximum iterations for convergence
  !iter_max=20000     !! Maximum iterations for convergence
  !--------------------------------------------------------- 
  !---------------------------------------------------------------------------------------

  if(dot_product(parameters,parameters)==0.d0) then
    write(*,*) "WARNING: All parameters initialized to zero"
    write(*,*) "Is that correct? Press Enter to continue"
    read(*,*)
  endif

  if(product(sigma)==0.d0) then
    write(*,*) "some of the sigma value is 0"
    stop
  endif

!  write(ch1,'(i3)')m
!  str=ch1//'f15.7)'

  sigma_inv=1.d0/sigma
  lambda=0.001d0;setup=1;converged=0

  saved_parameters=parameters
  chisq_tmp=0.d0
  do i=1,n
     call funcs(x(i),f,parameters,deriv,m)
    error(i)=y(i)-f
  enddo
  chisq_tmp = dot_product(error,error)/real(n)
 ! call levenberg_marquardt(x,y,sigma_inv,n,parameters,m,lambda,chisq,jacobian,error,work,lwork,converged,tolerance,setup)
  write(120,*) x,y !"entry",i,parameters, chisq!,lambda
  write(120,*) "entry",i,parameters, chisq_tmp!,lambda
!  if(dabs(chisq_tmp-chisq)<tol2)then
    !write(*,*) "chisq_temp2", chisq_tmp
!    parameters=saved_parameters
!    return! write(*,*) "chisq_temp2", chisq_tmp
!  endif

 
  !! First call sets up some variables, should be with setup=1
  call levenberg_marquardt(x,y,sigma_inv,n,parameters,m,lambda,chisq,jacobian,error,work,lwork,converged,tolerance,setup)

  if(write_files==1) then
    write(*,*) "Guess parameters=",parameters
    write(*,*) "Initial chisq=",chisq
    write(*,*) "##########################################"
  endif

  if(write_files==1)then
    open(10,file="output_mod_fit")
    write(10,*) "   #   parameters         chisq       lambda"
  endif
!------------------------------------------------------------------------------------
 !! Converging for optimal parameters
                  !     write(120,'(i5,21es20.11)') i,parameters, chisq!,lambda
  do i=1,iter_max
   ! if(write_files==1) !write(120,'(i5,21es14.17)') i,x,y,parameters!,chisq
                       !write(120,'(i5,21es20.11)') i,parameters, chisq!,lambda

                       !write(120,'(i5,2es30.18)') i,chisq!,lambda
   !saved_parameters=parameters
    call levenberg_marquardt(x,y,sigma_inv,n,parameters,m,lambda,chisq,jacobian,error,work,lwork,converged,tolerance,setup)
    !if(write_files==1) write(10,'(i5,3f15.7)') i,parameters,chisq,lambda
    !if(write_files==1) write(120,'(i5,20es15.7)') i,parameters,chisq,lambda
  if(dabs(chisq_tmp-chisq)<tol2)then
    !write(*,*) "chisq_temp2", chisq_tmp
     parameters=saved_parameters
    return! write(*,*) "chisq_temp2", chisq_tmp
  endif
    if(converged==1)exit
    !print *, "x after entering blahh", x
  enddo
  write(120,*)"exit", i,parameters, chisq!,lambda
  if(write_files==1)close(10)

  if(write_files==1) then
    if(converged==1) then
      write(*,*) "Converged parameters=",parameters
      write(*,*) "chisq=",chisq
    else
      write(*,*) "Unconverged parameters=",parameters
      write(*,*) "chisq=",chisq
      write(*,*) "Look for initial guess,tolerance,iter_max for convergence issues"
    endif

    open(10,file="data.out")
    write(10,*) "       x              fit             y"
    write(10,*) "###############################################"
    do i=1,n
      !call funcs(v1,c1,d1,x(i),f,parameters,deriv,m)
      call funcs(x(i),f,parameters,deriv,m)
      write(10,'(3es15.7)') x(i),f,y(i)
      write(110,'(3es15.7)') x(i),f,y(i)
    enddo
    close(10)

  endif

end subroutine optimize_param
!---------------------------------------------------------- 

subroutine levenberg_marquardt(x,y,sigma_inv,n,a,m,lambda,chisq,jacobian,error,work,lwork,converged,tolerance,setup)
  !! LEAST SQUARE SOLUTION OF CHISQ=SUM(F(X,A)-Y)**2
  !! Implementing algorithms from:
    !! Numerical Recipes in Fortan, Chapter 15, page 678-679
    !! http://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm
  !! Do not understand the technical details :(
  !! First call to it should be with setup=1
  implicit none
  integer i,n,m,converged,setup
  real*8 x(n),y(n),error(n),sigma_inv(n),v1,c1,d1
  real*8 a(m),deriv(m),a_new(m)
  real*8 lambda,chisq,chisq_new,f,lamb_plus_1,tolerance
  real*8 jacobian(n,m),transp_jacob(m,n)
  real*8 aa(m,m),bb(m)
  real*8,allocatable:: work(:)
  integer ipiv(m),lwork,info

  !! Default values of parameters if not provided
  if(lambda<=0.d0)lambda=0.001d0
  if(tolerance<=0.d0)tolerance=1.d-10
  lamb_plus_1=lambda+1.d0

  if(setup==0) then
    transp_jacob=transpose(jacobian)
  
    !! SOLVING LINEAR EQN AA*(DEL_A)=BB
    !! AA = J^T*J + LAMBDA*DIAG(J^T*J)
    !! BB = J^T*(Y-F)
    AA = matmul(transp_jacob,jacobian)
    do i=1,m
     aa(i,i) = aa(i,i)*lamb_plus_1
    enddo
    bb=matmul(transp_jacob,error)
  
    !! LAPACK subroutine to solve linear eqn.
    call dsysv('U', m, 1, AA, m, IPIV, BB, m, WORK, LWORK, INFO)
    if(info.ne.0) then
      write(*,*) "problem in dsysv",info
      stop
    endif

    a_new=a+bb

    !! Testing if chisq increased or decreased
    chisq_new=0.d0
    do i=1,n
      call funcs(x(i),f,a_new,deriv,m)
      error(i)=y(i)-f
    enddo
    chisq_new = dot_product(error,error)/real(n)

    converged=0
    if(dabs(chisq-chisq_new)<tolerance) then
   ! if(dabs(chisq_new)<tolerance) then
      converged=1
    else
      if(chisq_new>chisq) lambda=lambda*10
      if(chisq_new<chisq) then
        lambda=lambda/10.d0
        a=a_new
        chisq=chisq_new
      endif
    endif
  endif

  if(setup==1) then
    !! Done only on the first call
    chisq=0.d0
    do i=1,n
      call funcs(x(i),f,a,deriv,m)
      error(i)=(y(i)-f)*sigma_inv(i)
      jacobian(i,:)=deriv*sigma_inv(i)
    enddo
    chisq = dot_product(error,error)/real(n)
    lwork=-1;info=0
    allocate(work(m))
    call dsysv('U',m,1,AA,m,IPIV,BB,m,WORK,LWORK,INFO)
    if(info.ne.0) then
      write(*,*) "problem in dsysv"
      stop
    endif
    lwork=nint(work(1))
    deallocate(work)
    allocate(work(lwork))
    setup=0
  endif

end subroutine levenberg_marquardt
!---------------------------------------------------------- 

subroutine funcs(x,f,a,deriv,m)
!subroutine funcs(v1,c1,d1,x,f,a,deriv,m)
  !! provides the guess function(x,a), and its derivatives
  !! a(m) are the parameters
  implicit none
  integer m
  real*8 x,a(m),deriv(m),v1,c1,d1
  real*8,intent(out)::f

  f=a(1)+a(2)*(x-a(3))+a(4)*(x-a(3))**2
  deriv(1)=1.d0
  deriv(2)=x-a(3)
  deriv(3)=-a(2)-2.d0*a(4)*(x-a(3))
  deriv(4)=(x-a(3))**2
  !f=v1+c1*(x-a(1))+d1*(x-a(1))**2
  !deriv(1)=-c1-2.d0*d1*(x-a(1))
end subroutine funcs
!---------------------------------------------------------- 

end module mod_fit
