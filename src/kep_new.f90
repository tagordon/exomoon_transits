module kepler_new
use iso_c_binding
implicit none

real*8, parameter :: pi = 4.d0 * Atan(1.d0), G = 8.887692445125634d-10
real*8, parameter :: o3 = 1.d0 / 3.d0
real*8, parameter :: au_rsun = 215.03215567054764

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

subroutine grad_coords(t, ms, t0p, ep, Pp, Op, wp, ip, mp, &
    &t0m, em, Pm, wm, Om, im, mm, j, &
    &xp, yp, xm, ym, dxp, dyp, dxm, dym) bind(C, name="grad_coords")

    integer (c_int), bind(C) :: j
    real*8 :: np, nm, ap, am
    real*8 :: mrp, mrm, comegap, somegap, comegam, somegam
    real*8 :: cwp, swp, cip, cwm, swm, cim, sip, sim
    real*8, dimension(j) :: x, y, xbc, ybc
    real*8, dimension(j) :: r, cosf, sinf, cosfw, sinfw
    real (c_double), bind(C) :: ms, t0p, ep, Pp, Op, wp, ip, mp 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im, mm
    real (c_double), bind(C), dimension(j) :: t
    real (c_double), bind(C), dimension(j), intent(out) :: xp, yp, xm, ym
    
    real (c_double), bind(C), dimension(15, j), intent(out) :: dxp, dyp, dxm, dym
    real*8, dimension(j) :: xbc_m, ybc_m, x_m, y_m
    
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
    
    ap = (G * (ms + mp + mm) / (np ** 2.d0)) ** o3
    ap_m = G / (3 * (ap * np) ** 2.d0)
    ap_np = - 2 * ap / (3 * np)
    ap_Pp = ap_np * np_Pp

    am = (G * (mp + mm) / (nm ** 2.d0)) ** o3
    am_m = G / (3 * (am * nm) ** 2.d0)
    am_nm = - 2 * am / (3 * nm)
    am_Pm = am_nm * nm_Pm
    
    comegap = Cos(Op)
    somegap = Sin(Op)
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

    mtot = mp + mm
    omtot2 = 1.d0 / (mtot * mtot)
    
    mrp = -mm / mtot
    mrm = mp / mtot
    
    mrp_mm = -mp * omtot2
    mrp_mp = mm * omtot2
    mrm_mm = mrp_mm
    mrm_mp = mrp_mp
    mr = mp / mm
    
    call kepler_solve(np * (t - t0p), ep, cosf, sinf, f_e, f_M, j)
    
    r = ap * (1.d0 - ep * ep) / (1.d0 + ep * cosf)
    r_cosf = ap * ep * (ep * ep - 1.d0) / (1.d0 + ep * cosf)**2.d0
    r_t0 = r_cosf * sinf * f_M * np
    r_e = -r_cosf * sinf * f_e - ap * (cosf + 2 * ep + cosf * ep * ep) / (1.d0 + ep * cosf)**2.d0
    r_Pp = -(r_cosf * sinf * f_M * (t - t0p) + 2 * r / (3.d0 * np)) * np_Pp
     
    cosfw = cwp * cosf - swp * sinf
    sinfw = swp * cosf + sinf * cwp
    
    ccmss = comegap * cosfw - somegap * sinfw * cip
    cspsc = comegap * sinfw + somegap * cosfw * cip
    scpcs = somegap * cosfw + comegap * sinfw * cip
    ssmcc = somegap * sinfw - comegap * cosfw * cip
    xbc = -r * ccmss
    xbc_m = - r * ap_m * ccmss / ap
    ybc = -r * scpcs
    ybc_m = - r * ap_m * scpcs / ap
    
    dxp(1, :) = xbc_m
    dxm(1, :) = xbc_m
    dxp(2, :) = - r_t0 * ccmss - r * f_M * np * cspsc
    dxm(2, :) = dxp(2, :)
    dxp(3, :) = - r_e * ccmss + r * f_e * cspsc
    dxm(3, :) = dxp(3, :)
    dxp(4, :) = - r_Pp * ccmss + r * f_M * (t - t0p) * np_Pp * cspsc
    dxm(4, :) = dxp(4, :)
    dxp(5, :) = r * scpcs
    dxm(5, :) = dxp(5, :)
    dxp(6, :) = r * cspsc
    dxm(6, :) = dxp(6, :)
    dxp(7, :) = - r * somegap * sinfw * sip
    dxm(7, :) = dxp(7, :)
        
    ! ms
    dyp(1, :) = ybc_m
    dym(1, :) = ybc_m
    ! t0p
    dyp(2, :) = - r_t0 * scpcs - r * f_M * np * ssmcc
    dym(2, :) = dyp(2, :)
    ! ep
    dyp(3, :) = - r_e * scpcs + r * f_e * ssmcc
    dym(3, :) = dyp(3, :)
    ! Pp
    dyp(4, :) = - r_Pp * scpcs + r * f_M * (t - t0p) * np_Pp * ssmcc
    dym(4, :) = dyp(4, :)
    ! Op
    dyp(5, :) = - r * ccmss
    dym(5, :) = dyp(5, :)
    ! wp
    dyp(6, :) = r * ssmcc
    dym(6, :) = dyp(6, :)
    ! ip
    dyp(7, :) = r * comegap * sinfw * sip
    dym(7, :) = dyp(7, :)
    
    call kepler_solve(nm * (t - t0m), em, cosf, sinf, f_e, f_M, j)
    
    r = am * (1.d0 - em * em) / (1.d0 + em * cosf)
    r_cosf = am * em * (em * em - 1.d0) / (1.d0 + em * cosf)**2.d0
    r_t0 = r_cosf * sinf * f_M * nm
    r_e = -r_cosf * sinf * f_e - am * (cosf + 2 * em + cosf * em * em) / (1.d0 + em * cosf)**2.d0
    r_Pm = -(r_cosf * sinf * f_M * (t - t0m) + 2 * r / (3.d0 * nm)) * nm_Pm
     
    cosfw = cwp * cosf - swp * sinf
    sinfw = swp * cosf + sinf * cwp
    
    ccmss = comegam * cosfw - somegam * sinfw * cim
    cspsc = comegam * sinfw + somegam * cosfw * cim
    scpcs = somegam * cosfw + comegam * sinfw * cim
    ssmcc = somegam * sinfw - comegam * cosfw * cim
    x = -r * ccmss
    x_m = - r * am_m * ccmss / am
    y = -r * scpcs
    y_m = - r * am_m * scpcs / am
        
    ! ms, t0p, ep, Pp, Op, wp, ip, mp, t0m, em, Pm, wm, Om, im, mm
    xp = xbc + x * mrp
    dxp(8, :) = xbc_m + x * mrp_mp + x_m * mrp
    dxp(9, :) = (- r_t0 * ccmss - r * f_M * nm * cspsc) * mrp
    dxp(10, :) = (- r_e * ccmss + r * f_e * cspsc) * mrp
    dxp(11, :) = (- r_Pm * ccmss + r * f_M * (t - t0m) * nm_Pm * cspsc) * mrp
    dxp(12, :) = r * cspsc * mrp
    dxp(13, :) = r * scpcs * mrp
    dxp(14, :) = -r * somegam * sinfw * sim * mrp
    dxp(15, :) = xbc_m + x_m * mrp + x * mrp_mm
    
    yp = ybc + y * mrp

    dyp(8, :) = ybc_m + y * mrp_mp + y_m * mrp
    dyp(9, :) = (- r_t0 * scpcs - r * f_M * nm * ssmcc) * mrp
    dyp(10, :) = (- r_e * scpcs + r * f_e * ssmcc) * mrp
    dyp(11, :) = (- r_Pm * scpcs + r * f_M * (t - t0m) * nm_Pm * ssmcc) * mrp
    dyp(12, :) = r * ssmcc * mrp
    dyp(13, :) = - r * ccmss * mrp
    dyp(14, :) = r * comegam * sinfw * sim * mrp
    dyp(15, :) = ybc_m + y_m * mrp + y * mrp_mm

    xm = xbc + x * mrm

    dxm(8, :) = xbc_m + x * mrm_mp + x_m * mrm
    dxm(9, :) = -dxp(9, :) * mr
    dxm(10, :) = -dxp(10, :) * mr
    dxm(11, :) = -dxp(11, :) * mr
    dxm(12, :) = -dxp(12, :) * mr
    dxm(13, :) = -dxp(13, :) * mr
    dxm(14, :) = -dxp(14, :) * mr
    dxm(15, :) = xbc_m + x * mrm_mm + x_m * mrm
    
    ym = ybc + y * mrm

    dym(8, :) = ybc_m + y * mrm_mp + y_m * mrm
    dym(9, :) = -dyp(9, :) * mr
    dym(10, :) = -dyp(10, :) * mr
    dym(11, :) = -dyp(11, :) * mr
    dym(12, :) = -dyp(12, :) * mr
    dym(13, :) = -dyp(13, :) * mr
    dym(14, :) = -dyp(14, :) * mr
    dyp(15, :) = ybc_m + y * mrm_mm + y_m * mrm
    
end

subroutine grad_impacts(t, ms, t0p, ep, Pp, Op, wp, ip, mp, &
    &t0m, em, Pm, wm, Om, im, mm, j, &
    &bp, bpm, theta, dbp, dbpm, dtheta) bind(C, name="grad_impacts")

    integer :: i
    real*8 :: a, b, c, tmp, mu
    real*8, dimension(j) :: bm, bm2, bp2, bpm2
    real*8, dimension(j) :: xp, yp, xm, ym
    real*8, dimension(j) :: bm_xm, bm_ym, bp_xp, bp_yp, bpm_xm, bpm_ym, bpm_xp, bpm_yp
    real*8, dimension(j) :: theta_bp, theta_bpm, theta_bm
    real*8, dimension(15, j) :: dxp, dyp, dxm, dym, dbm
    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: ms, t0p, ep, Pp, Op, wp, ip, mp 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im, mm
    real (c_double), bind(C), dimension(j) :: t
    real (c_double), bind(C), intent(out), dimension(j) :: bp, bpm, theta
    real (c_double), bind(C), intent(out), dimension(15, j) :: dbp, dbpm, dtheta
    
    call grad_coords(t, ms, t0p, ep, Pp, Op, wp, ip, mp, t0m, em, &
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
    
    ! optimization: compute Sin(theta) once and return sth, cth instead of theta directly.
    theta_bm = bm / (bp * bm * Sin(theta))
    theta_bp = ((bpm - bm) * (bpm + bm) - bp2) / (2 * bp2 * bpm * Sin(theta))
    theta_bpm = ((bp - bm) * (bp + bm) - bpm2) / (2 * bpm2 * bp * Sin(theta))
    
    do i=1,15,1
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

end