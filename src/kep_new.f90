module kepler_new
use iso_c_binding
implicit none

real*8, parameter :: pi = 4.d0 * Atan(1.d0), G = 8.887692445125634d-10
real*8, parameter :: o3 = 1.d0 / 3.d0

contains

subroutine kepler_solve(M, ecc, cosf, sinf, f_e, f_M, j) bind(C, name="kepler_solve")

    integer :: j, i
    real*8 :: tol, err, x, ecc
    real (c_double), bind(C), dimension(j) :: M
    real (c_double), bind(C), dimension(j), intent(out) :: cosf, sinf, f_e, f_M
    real*8, dimension(j) :: tanfhalf, tanfhalf2, denom, E, sE, cE
    
    if (ecc .lt. 1.d-10) then
        cE = Cos(M)
        sE = Sin(M)
        
        tanfhalf = sE / (1.d0 + cE)
        tanfhalf2 = tanfhalf * tanfhalf
        denom = 1.d0 / (tanfhalf2 + 1.d0)
        cosf = (1.d0 - tanfhalf2) * denom
        sinf = 2 * tanfhalf * denom
        
        f_M = 1.d0
        f_e = 2.d0 * sinf
    else
        tol = 1.d-10
        
        E = M
        x = Sqrt((1.d0 + ecc) / (1.d0 - ecc))
        
        do i=1,j,1
            err = 10.d0
            do while (abs(err) .gt. tol)
                err = - (E(i) - ecc * Sin(E(i)) - M(i)) / (1.d0 - ecc * Cos(E(i)))
                E(i) = E(i) + err
            end do
        end do    
        
        cE = Cos(E)
        tanfhalf = x * Sin(E) / (1.d0 + cE)
        tanfhalf2 = tanfhalf * tanfhalf
        denom = 1.d0 / (tanfhalf2 + 1.d0)
        cosf = (1.d0 - tanfhalf2) * denom
        sinf = 2 * tanfhalf * denom
        
        x = 1.d0 - ecc**2.d0
        f_M = (1.d0 + ecc * cosf)**2.d0 / x**1.5d0
        f_e = (2.d0 + ecc * cosf) * sinf / x
    end if
end

end