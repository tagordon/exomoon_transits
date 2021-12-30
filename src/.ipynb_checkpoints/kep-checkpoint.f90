module kepler
use iso_c_binding
implicit none

real*8, parameter :: pi = 4.d0 * Atan(1.d0), G = 8.8876413d-10
real*8, parameter :: o3 = 1.d0 / 3.d0
real*8, parameter :: days_in_year = 365.256, earths_in_sun = 332946.08
real*8, parameter :: au_rsun = 215.03215567054764

contains

subroutine solve_kepler(t, n, t0, ecc, a, r, cosf, sinf, j)

    integer :: j, i
    real*8 :: n, t0, ecc, a, tol, x, err
    real*8, dimension(j) :: t, M, E, sE, cE
    real*8, dimension(j), intent(out) :: r, cosf, sinf
    real*8, dimension(j) :: y, y2, denom
    
    M = n * (t - t0)
    if (ecc .lt. 1.d-5) then
        cosf = Cos(M)
        sinf = Sin(M)
        r = a
    else
        tol = 1.d-15
        
        E = M + ecc
        x = Sqrt((1.d0 + ecc) / (1.d0 - ecc))
        
        do i=1,j,1
            err = 1.d0
            do while (err .gt. tol)
                err = - (E(i) - ecc * Sin(E(i)) - M(i)) / (1.d0 - ecc * Cos(E(i)))
                E(i) = E(i) + err
            end do
        end do    
        
        cE = Cos(E)
        y = x * Sin(E) / (1.d0 + cE)
        y2 = y * y
        denom = 1.d0 / (y2 + 1.d0)
        cosf = (1.d0 - y2) * denom
        sinf = 2 * y * denom
        r = a * (1.d0 - ecc * cE)
    end if
end

subroutine coords(t, ms, t0p, ep, Pp, Op, wp, ip, mp, &
    &t0m, em, Pm, wm, Om, im, mm, j, &
    &xp, yp, zp, xm, ym, zm) bind(C, name="coords")

    real*8 :: np, nm, ap, am
    real*8 :: mrp, mrm, comegap, somegap, comegam, somegam
    real*8 :: cwp, swp, cip, cwm, swm, cim
    real*8, dimension(j) :: x, y, z, xbc, ybc, zbc
    real*8, dimension(j) :: r, cosf, sinf, cosfw, sinfw
    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: ms, t0p, ep, Pp, Op, wp, ip, mp 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im, mm
    real (c_double), bind(C), dimension(j) :: t
    real (c_double), bind(C), dimension(j), intent(out) :: xp, yp, zp, xm, ym, zm
    
    np = 2 * pi / Pp
    nm = 2 * pi / Pm
    ap = (G * (ms + mp + mm) / (np ** 2.d0)) ** o3
    am = (G * (mp + mm) / (nm ** 2.d0)) ** o3
    
    comegap = Cos(Op)
    somegap = Sin(Op)
    cwp = Cos(wp)
    swp = Sin(wp)
    cip = Cos(ip)
    
    comegam = Cos(Om)
    somegam = Sin(Om)
    cwm = Cos(wm)
    swm = Sin(wm)
    cim = Cos(im)
    
    mrp = -mm / (mp + mm)
    mrm = mp / (mp + mm)
    
    call solve_kepler(t, np, t0p, ep, ap, r, cosf, sinf, j)
    
    cosfw = cwp * cosf - swp * sinf
    sinfw = swp * cosf + sinf * cwp
    xbc = -r * (comegap * cosfw - somegap * sinfw * cip)
    ybc = -r * (somegap * cosfw + comegap * sinfw * cip)
    zbc = r * sinfw * Sin(ip)
    
    call solve_kepler(t, nm, t0m, em, am, r, cosf, sinf, j)
    
    cosfw = cwm * cosf - swm * sinf
    sinfw = swm * cosf + sinf * cwm
    x = -r * (comegam * cosfw - somegam * sinfw * cim)
    y = -r * (somegam * cosfw + comegam * sinfw * cim)
    z = r * sinfw * Sin(im)
    
    xp = xbc + x * mrp
    yp = ybc + y * mrp
    zp = zbc + z * mrp
    
    xm = xbc + x * mrm
    ym = ybc + y * mrm
    zm = zbc + z * mrm
    
end

subroutine input_coords(t, ms, t0p, ep, Pp, Op, wp, ip, mp, &
    &t0m, em, Pm, wm, Om, im, mm, j, &
    &xp, yp, zp, xm, ym, zm, bp, bpm, theta) bind(C, name="input_coords")

    integer :: i
    real*8 :: a, b, c, tmp, mu
    real*8, dimension(j) :: bm
    real*8, dimension(j), intent(out) :: xp, yp, zp, xm, ym, zm
    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: ms, t0p, ep, Pp, Op, wp, ip, mp 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im, mm
    real (c_double), bind(C), dimension(j) :: t
    real (c_double), bind(C), intent(out), dimension(j) :: bp, bpm, theta
    
    call coords(t, ms, t0p, ep, Pp, Op, wp, ip, mp, t0m, em, &
                Pm, wm, Om, im, mm, j, xp, yp, zp, xm, ym, zm)
    
    bm = Sqrt(xm**2.d0 + ym**2.d0) * au_rsun
    bp = Sqrt(xp**2.d0 + yp**2.d0) * au_rsun
    bpm = Sqrt((xm - xp)**2.d0 + (ym - yp)**2.d0) * au_rsun
    
    do i=1,j,1
    
        a = bp(i) 
        b = bpm(i)
        c = bm(i)
        
        if (a .GT. b) then
            tmp = b
            b = a
            a = tmp
        end if
        if (b .GT. c) then
            mu = c - (a - b)
        else
            mu = b - (a - c)
        end if
            
        theta(i) = 2 * Atan(Sqrt(((a - b) + c) * mu / ((a + (b + c)) * ((a - c) + b))))
    end do
    
end

end