module kepler_new
use iso_c_binding
implicit none

real*8, parameter :: pi = 4.d0 * Atan(1.d0), G = 8.887692445125634d-10
real*8, parameter :: o3 = 1.d0 / 3.d0
real*8, parameter :: au_rsun = 215.03215567054764

contains

subroutine kepler_solve(M, ecc, cosf, sinf, j) bind(C, name="kepler_solve")

    integer :: j, i
    real*8 :: tol, err, x, ecc
    real (c_double), bind(C), dimension(j) :: M
    real (c_double), bind(C), dimension(j), intent(out) :: cosf, sinf
    real*8, dimension(j) :: tanfhalf, tanfhalf2, denom, E, sE, cE
    
    if (ecc .lt. 1.d-10) then
        cE = Cos(M)
        sE = Sin(M)
        
        tanfhalf = sE / (1.d0 + cE)
        tanfhalf2 = tanfhalf * tanfhalf
        denom = 1.d0 / (tanfhalf2 + 1.d0)
        cosf = (1.d0 - tanfhalf2) * denom
        sinf = 2 * tanfhalf * denom
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
    end if
end

subroutine grad_kepler_solve(M, ecc, cosf, sinf, f_e, f_M, j) bind(C, name="grad_kepler_solve")

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

subroutine coords(t, ap, t0p, ep, Pp, wp, ip, am, &
    &t0m, em, Pm, Om, wm, im, mm, j, &
    &xp, yp, zp, xm, ym, zm) bind(C, name="coords")
    
    integer (c_int), bind(C) :: j
    real*8 :: np, nm
    real*8 :: mrp, mrm, mtot, comegam, somegam
    real*8 :: cwp, swp, cip, cwm, swm, cim
    real*8, dimension(j) :: x, y, z, xbc, ybc, zbc
    real*8, dimension(j) :: r, cosf, sinf, cosfw, sinfw
    real (c_double), bind(C) :: ap, t0p, ep, Pp, wp, ip, am 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im, mm
    real (c_double), bind(C), dimension(j) :: t
    real (c_double), bind(C), dimension(j), intent(out) :: xp, yp, zp, xm, ym, zm
            
    np = 2 * pi / Pp    
    nm = 2 * pi / Pm

    cwp = Cos(wp)
    swp = Sin(wp)
    cip = Cos(ip)
    
    comegam = Cos(Om)
    somegam = Sin(Om)
    cwm = Cos(wm)
    swm = Sin(wm)
    cim = Cos(im)

    mrp = - 1.d0 / (1.d0 + mm)
    mrm = - mm * mrp
    
    call kepler_solve(np * (t - t0p), ep, cosf, sinf, j)
    
    r = ap * (1.d0 - ep * ep) / (1.d0 + ep * cosf)
    cosfw = cwp * cosf - swp * sinf
    sinfw = swp * cosf + sinf * cwp
    xbc = -r * cosfw
    ybc = -r * sinfw * cip
    zbc = r * sinfw * Sin(ip)
    
    call kepler_solve(nm * (t - t0m), em, cosf, sinf, j)
    
    r = am * (1.d0 - em * em) / (1.d0 + em * cosf)
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

subroutine grad_coords(t, ap, t0p, ep, Pp, wp, ip, am, &
    &t0m, em, Pm, Om, wm, im, mm, j, &
    &xp, yp, xm, ym, dxp, dyp, dxm, dym) bind(C, name="grad_coords")

    integer (c_int), bind(C) :: j
    real*8 :: np, nm
    real*8 :: mrp, mrm, comegam, somegam
    real*8 :: cwp, swp, cip, cwm, swm, cim, sip, sim
    real*8, dimension(j) :: x, y, xbc, ybc
    real*8, dimension(j) :: r, cosf, sinf, cosfw, sinfw, denom
    real (c_double), bind(C) :: ap, t0p, ep, Pp, wp, ip, am 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im, mm
    real (c_double), bind(C), dimension(j) :: t
    real (c_double), bind(C), dimension(j), intent(out) :: xp, yp, xm, ym
    
    real (c_double), bind(C), dimension(14, j), intent(out) :: dxp, dyp, dxm, dym
    real*8, dimension(j) :: xbc_a, ybc_a, x_a, y_a
    
    real*8 :: np_Pp, nm_Pm
    real*8 :: np_ms, np_mp, np_mm, nm_mp, nm_mm
    real*8 :: ap_m, am_m, am_Pm, ap_Pp, am_nm, ap_np
    real*8 :: mtot, omtot2, mrp_mm, mrp_mp, mrm_mp, mrm_mm, mr
    real*8, dimension(j) :: r_t0, r_e, r_cosf
    real*8, dimension(j) :: f_M, f_e
    real*8, dimension(j) :: ccmss, cspsc, scpcs, ssmcc
    real*8, dimension(j) :: r_Pp, r_Pm
        
    np = 2 * pi / Pp
    np_Pp = -np / Pp
    
    nm = 2 * pi / Pm
    nm_Pm = -nm / Pm
    
    cwp = Cos(wp)
    swp = Sin(wp)
    cip = Cos(ip)
    sip = Sin(ip)
    
    comegam = Cos(Om)
    somegam = Sin(Om)
    cwm = Cos(wm)
    swm = Sin(wm)
    cim = Cos(im)
    sim = Sin(im)

    mrp = - 1.d0 / (1.d0 + mm)
    mrm = - mm * mrp
    
    mrp_mm = - mrp / (1.d0 + mm)
    mrm_mm = mrp_mm
    
    call grad_kepler_solve(np * (t - t0p), ep, cosf, sinf, f_e, f_M, j)
    
    denom = 1.d0 / (1.d0 + ep * cosf)
    r = ap * (1.d0 - ep * ep) * denom
    r_cosf = - r * ep * denom
    r_t0 = r_cosf * sinf * f_M * np
    r_e = -r_cosf * sinf * f_e - ap * (cosf + 2 * ep + cosf * ep * ep) * denom * denom
    r_Pp = -(r_cosf * sinf * f_M * (t - t0p) + 2 * r / (3.d0 * np)) * np_Pp
     
    cosfw = cwp * cosf - swp * sinf
    sinfw = swp * cosf + sinf * cwp
    
    ccmss = cosfw
    cspsc = sinfw
    scpcs = sinfw * cip
    ssmcc = - cosfw * cip
    xbc = -r * ccmss
    xbc_a = - r * ccmss / ap
    ybc = -r * scpcs
    ybc_a = - r * scpcs / ap
    
    dxp(1, :) = xbc_a
    dxm(1, :) = xbc_a
    dxp(2, :) = - r_t0 * ccmss - r * f_M * np * cspsc
    dxm(2, :) = dxp(2, :)
    dxp(3, :) = - r_e * ccmss + r * f_e * cspsc
    dxm(3, :) = dxp(3, :)
    dxp(4, :) = - r_Pp * ccmss + r * f_M * (t - t0p) * np_Pp * cspsc
    dxm(4, :) = dxp(4, :)
    dxp(5, :) = r * cspsc
    dxm(5, :) = dxp(5, :)
    dxp(6, :) = 0.d0
    dxm(6, :) = dxp(6, :)
        
    ! ms
    dyp(1, :) = ybc_a
    dym(1, :) = ybc_a
    ! t0p
    dyp(2, :) = - r_t0 * scpcs - r * f_M * np * ssmcc
    dym(2, :) = dyp(2, :)
    ! ep
    dyp(3, :) = - r_e * scpcs + r * f_e * ssmcc
    dym(3, :) = dyp(3, :)
    ! Pp
    dyp(4, :) = - r_Pp * scpcs + r * f_M * (t - t0p) * np_Pp * ssmcc
    dym(4, :) = dyp(4, :)
    ! wp
    dyp(5, :) = r * ssmcc
    dym(5, :) = dyp(5, :)
    ! ip
    dyp(6, :) = r * sinfw * sip
    dym(6, :) = dyp(6, :)
    
    call grad_kepler_solve(nm * (t - t0m), em, cosf, sinf, f_e, f_M, j)
    
    denom = 1.d0 / (1.d0 + em * cosf)
    r = am * (1.d0 - em * em) * denom
    r_cosf = - em * r * denom
    r_t0 = r_cosf * sinf * f_M * nm
    r_e = -r_cosf * sinf * f_e - am * (cosf + 2 * em + cosf * em * em) * denom * denom
    r_Pm = -(r_cosf * sinf * f_M * (t - t0m) + 2 * r / (3.d0 * nm)) * nm_Pm
     
    cosfw = cwm * cosf - swm * sinf
    sinfw = swm * cosf + sinf * cwm
    
    ccmss = comegam * cosfw - somegam * sinfw * cim
    cspsc = comegam * sinfw + somegam * cosfw * cim
    scpcs = somegam * cosfw + comegam * sinfw * cim
    ssmcc = somegam * sinfw - comegam * cosfw * cim
    x = -r * ccmss
    x_a = -r * ccmss / am
    y = -r * scpcs
    y_a = -r * scpcs / am
        
    ! ms, t0p, ep, Pp, Op, wp, ip, mp, t0m, em, Pm, wm, Om, im, mm
    xp = xbc + x * mrp
    dxp(7, :) = x_a * mrp
    dxp(8, :) = (- r_t0 * ccmss - r * f_M * nm * cspsc) * mrp
    dxp(9, :) = (- r_e * ccmss + r * f_e * cspsc) * mrp
    dxp(10, :) = (- r_Pm * ccmss + r * f_M * (t - t0m) * nm_Pm * cspsc) * mrp
    dxp(11, :) = r * cspsc * mrp
    dxp(12, :) = r * scpcs * mrp
    dxp(13, :) = -r * somegam * sinfw * sim * mrp
    dxp(14, :) = x * mrp_mm
    
    yp = ybc + y * mrp

    dyp(7, :) = y_a * mrp
    dyp(8, :) = (- r_t0 * scpcs - r * f_M * nm * ssmcc) * mrp
    dyp(9, :) = (- r_e * scpcs + r * f_e * ssmcc) * mrp
    dyp(10, :) = (- r_Pm * scpcs + r * f_M * (t - t0m) * nm_Pm * ssmcc) * mrp
    dyp(11, :) = r * ssmcc * mrp
    dyp(12, :) = - r * ccmss * mrp
    dyp(13, :) = r * comegam * sinfw * sim * mrp
    dyp(14, :) = y * mrp_mm

    xm = xbc + x * mrm

    dxm(7, :) = x_a * mrm
    dxm(8, :) = -dxp(8, :) * mm
    dxm(9, :) = -dxp(9, :) * mm
    dxm(10, :) = -dxp(10, :) * mm
    dxm(11, :) = -dxp(11, :) * mm
    dxm(12, :) = -dxp(12, :) * mm
    dxm(13, :) = -dxp(13, :) * mm
    dxm(14, :) = x * mrm_mm
    
    ym = ybc + y * mrm

    dym(7, :) = y_a * mrm
    dym(8, :) = -dyp(8, :) * mm
    dym(9, :) = -dyp(9, :) * mm
    dym(10, :) = -dyp(10, :) * mm
    dym(11, :) = -dyp(11, :) * mm
    dym(12, :) = -dyp(12, :) * mm
    dym(13, :) = -dyp(13, :) * mm
    dym(14, :) = y * mrm_mm
    
end

subroutine grad_impacts(t, ap, t0p, ep, Pp, wp, ip, am, &
    &t0m, em, Pm, Om, wm, im, mm, j, &
    &bp, bpm, theta, dbp, dbpm, dtheta) bind(C, name="grad_impacts")

    integer :: i
    real*8 :: a, b, c, tmp, mu
    real*8, dimension(j) :: bm, bm2, bp2, bpm2, sth
    real*8, dimension(j) :: xp, yp, xm, ym
    real*8, dimension(j) :: bm_xm, bm_ym, bp_xp, bp_yp, bpm_xm, bpm_ym, bpm_xp, bpm_yp
    real*8, dimension(j) :: theta_bp, theta_bpm, theta_bm
    real*8, dimension(14, j) :: dxp, dyp, dxm, dym, dbm
    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: ap, t0p, ep, Pp, wp, ip, am 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im, mm
    real (c_double), bind(C), dimension(j) :: t
    real (c_double), bind(C), intent(out), dimension(j) :: bp, bpm, theta
    real (c_double), bind(C), intent(out), dimension(14, j) :: dbp, dbpm, dtheta
    
    call grad_coords(t, ap, t0p, ep, Pp, wp, ip, am, t0m, em, &
                Pm, wm, Om, im, mm, j, xp, yp, xm, ym, &
                dxp, dyp, dxm, dym)
    
    ! ms, t0p, ep, Pp, Op, wp, ip, mp, t0m, em, Pm, wm, Om, im, mm
    bm2 = xm**2.d0 + ym**2.d0
    bm = Sqrt(bm2)
    bm_xm = xm / bm
    bm_ym = ym / bm
    
    bp2 = xp**2.d0 + yp**2.d0
    bp = Sqrt(bp2)
    bp_xp = xp / bp
    bp_yp = yp / bp
    
    bpm2 = (xm - xp)**2.d0 + (ym - yp)**2.d0
    bpm = Sqrt(bpm2)
    bpm_xm = (xm - xp) / bpm
    bpm_ym = (ym - yp) / bpm
    bpm_xp = - bpm_xm
    bpm_yp = - bpm_ym
    
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
    
    sth = Sin(theta)
    theta_bm = bm / (bp * bpm * sth)
    theta_bp = ((bpm - bm) * (bpm + bm) - bp2) / (2 * bp2 * bpm * sth)
    theta_bpm = ((bp - bm) * (bp + bm) - bpm2) / (2 * bpm2 * bp * sth)
    
    do i=1,14,1
        dbm(i, :) = bm_xm * dxm(i, :) + bm_ym * dym(i, :)
        dbp(i, :) = bp_xp * dxp(i, :) + bp_yp * dyp(i, :)
        dbpm(i, :) = bpm_xm * dxm(i, :) + bpm_ym * dym(i, :) &
                   + bpm_xp * dxp(i, :) + bpm_yp * dyp(i, :)
        dtheta(i, :) = theta_bm * dbm(i, :) &
                     + theta_bp * dbp(i, :) &
                     + theta_bpm * dbpm(i, :)
    end do

    bp = bp * au_rsun
    bpm = bpm * au_rsun
    dbp = dbp * au_rsun
    dbpm = dbpm * au_rsun
            
end

subroutine impacts(t, ap, t0p, ep, Pp, wp, ip, am, &
    &t0m, em, Pm, Om, wm, im, mm, j, bp, bpm, theta) bind(C, name="impacts")

    integer :: i
    real*8 :: a, b, c, tmp, mu
    real*8, dimension(j) :: bm, bm2, bp2, bpm2
    real*8, dimension(j) :: xp, yp, xm, ym, zp, zm
    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: ap, t0p, ep, Pp, wp, ip, am 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im, mm
    real (c_double), bind(C), dimension(j) :: t
    real (c_double), bind(C), intent(out), dimension(j) :: bp, bpm, theta
    
    call coords(t, ap, t0p, ep, Pp, wp, ip, am, t0m, em, &
                Pm, wm, Om, im, mm, j, xp, yp, zp, xm, ym, zm)
    
    bm2 = xm**2.d0 + ym**2.d0
    bm = Sqrt(bm2)
    
    bp2 = xp**2.d0 + yp**2.d0
    bp = Sqrt(bp2)
    
    bpm2 = (xm - xp)**2.d0 + (ym - yp)**2.d0
    bpm = Sqrt(bpm2)
    
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

    bp = bp * au_rsun
    bpm = bpm * au_rsun
            
end

end