submodule (reconstruction) weno7
implicit none
! seventh order weno coefficients
double precision,parameter :: a400 =-1.d0/4.d0
double precision,parameter :: a401 = 13.d0/12.d0
double precision,parameter :: a402 =-23.d0/12.d0
double precision,parameter :: a403 = 25.d0/12.d0
double precision,parameter :: a410 = 1.d0/12.d0
double precision,parameter :: a411 =-5.d0/12.d0
double precision,parameter :: a412 = 13.d0/12.d0
double precision,parameter :: a413 = 3.d0/12.d0
double precision,parameter :: a420 =-1.d0/12.d0
double precision,parameter :: a421 = 7.d0/12.d0
double precision,parameter :: a422 = 7.d0/12.d0
double precision,parameter :: a423 =-1.d0/12.d0
double precision,parameter :: a430 = 1.d0/4.d0
double precision,parameter :: a431 = 13.d0/12.d0
double precision,parameter :: a432 =-5.d0/12.d0
double precision,parameter :: a433 = 1.d0/12.d0
double precision,parameter :: c40 = 1.d0/35.d0
double precision,parameter :: c41 = 12.d0/35.d0
double precision,parameter :: c42 = 18.d0/35.d0
double precision,parameter :: c43 = 4.d0/35.d0

double precision,parameter :: ss   = 1.0e-20

contains
end submodule