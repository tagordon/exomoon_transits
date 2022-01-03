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

subroutine solve_kepler_grad(t, n, t0, ecc, a, r, cosf, sinf, &
                           & r_n, r_t0, r_e, r_a, &
                           & cosf_n, cosf_t0, cosf_e, cosf_a, & 
                           & sinf_n, sinf_t0, sinf_e, sinf_a, j)

    integer :: j, i
    real*8 :: n, t0, ecc, ecc2, a, tol, x, err
    real*8, dimension(j) :: t, M, E, sE, cE, inv_tanf, cosf_r, sinf_r
    real*8, dimension(j), intent(out) :: r, cosf, sinf
    real*8, dimension(j), intent(out) :: r_n, r_t0, r_e, r_a 
    real*8, dimension(j), intent(out) :: cosf_n, cosf_t0, cosf_e, cosf_a
    real*8, dimension(j), intent(out) :: sinf_n, sinf_t0, sinf_e, sinf_a
    real*8, dimension(j) :: y, y2, denom
    
    M = n * (t - t0)
    if (ecc .lt. 1.d-5) then
        cosf = Cos(M)
        sinf = Sin(M)
        r = a
        
        r_a = 1.d0
        r_t0 = 0.d0
        r_e = 0.d0
        r_n = 0.d0
        
        cosf_a = 0.d0
        cosf_t0 = sinf * n
        cosf_e = 0.d0
        cosf_n = -sinf * (t - t0)
        
        sinf_a = 0.d0
        sinf_t0 = cosf * n
        sinf_e = 0.d0
        sinf_n = cosf * (t - t0)
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
        sE = Sin(E)
        y = x * sE / (1.d0 + cE)
        y2 = y * y
        denom = 1.d0 / (y2 + 1.d0)
        cosf = (1.d0 - y2) * denom
        sinf = 2 * y * denom
        
        r = a * (1.d0 - ecc * cE)
        r_a = 1.d0 - ecc * cE
        r_t0 = n * a * ecc * sE / (ecc * cE - 1.d0)
        r_e = a * (ecc * sE * sE / (1.d0 - ecc * cE) - cE)
        r_n = (t - t0) * a * ecc * sE / (1.d0 - cE)
        
        ecc2 = ecc * ecc
        x = (1.d0 - ecc2)
        cosf_r = - a * x / (r * r * ecc)
        cosf_a = cosf_r * r_a + x / (r * ecc)
        cosf_t0 = cosf_r * r_t0
        cosf_e = cosf_r * r_e + 1.d0 / ecc2 - a * x / (r * ecc2) - 2 * a / r
        cosf_n = cosf_r * r_n
        
        inv_tanf = cosf / sinf
        sinf_r = - inv_tanf * cosf_r
        sinf_a = - inv_tanf * cosf_a
        sinf_t0 = - inv_tanf * cosf_t0
        sinf_e = - inv_tanf * cosf_e
        sinf_n = - inv_tanf * cosf_n
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

! add derivatives for each cartesian component
subroutine coords_grad(t, ms, t0p, ep, Pp, Op, wp, ip, mp, &
    &t0m, em, Pm, wm, Om, im, mm, j, &
    &xp, yp, zp, xm, ym, zm, dxp, dyp, dzp, dxm, dym, dzm) bind(C, name="coords_grad")

    integer (c_int), bind(C) :: j
    real*8 :: np, nm, ap, am
    real*8 :: mrp, mrm, comegap, somegap, comegam, somegam
    real*8 :: cwp, swp, cip, cwm, swm, cim, sip, sim
    real*8, dimension(j) :: x, y, z, xbc, ybc, zbc
    real*8, dimension(j) :: r, cosf, sinf, cosfw, sinfw
    real (c_double), bind(C) :: ms, t0p, ep, Pp, Op, wp, ip, mp 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im, mm
    real (c_double), bind(C), dimension(j) :: t
    real (c_double), bind(C), dimension(j), intent(out) :: xp, yp, zp, xm, ym, zm
    
    real (c_double), bind(C), dimension(15, j), intent(out) :: dxp, dyp, dzp, dxm, dym, dzm
    real*8, dimension(15, j) :: dx, dy, dz
    real*8, dimension(4, j) :: dcosf, dsinf, dr
    real*8, dimension(j) :: xbc_ms, xbc_t0p, xbc_ep, xbc_Pp, xbc_Op, xbc_wp, xbc_ip, xbc_mp, xbc_mm
    real*8, dimension(j) :: ybc_ms, ybc_t0p, ybc_ep, ybc_Pp, ybc_Op, ybc_wp, ybc_ip, ybc_mp, ybc_mm
    real*8, dimension(j) :: zbc_ms, zbc_t0p, zbc_ep, zbc_Pp, zbc_Op, zbc_wp, zbc_ip, zbc_mp, zbc_mm
    
    real*8, dimension(j) :: xbc_t0m, xbc_em, xbc_Pm, xbc_Om, xbc_wm, xbc_im
    real*8, dimension(j) :: ybc_t0m, ybc_em, ybc_Pm, ybc_Om, ybc_wm, ybc_im
    real*8, dimension(j) :: zbc_t0m, zbc_em, zbc_Pm, zbc_Om, zbc_wm, zbc_im
    
    real*8, dimension(j) :: x_ms, x_t0p, x_ep, x_Pp, x_Op, x_wp, x_ip, x_mp, x_mm
    real*8, dimension(j) :: y_ms, y_t0p, y_ep, y_Pp, y_Op, y_wp, y_ip, y_mp, y_mm
    real*8, dimension(j) :: z_ms, z_t0p, z_ep, z_Pp, z_Op, z_wp, z_ip, z_mp, z_mm
    
    real*8, dimension(j) :: x_t0m, x_em, x_Pm, x_Om, x_wm, x_im
    real*8, dimension(j) :: y_t0m, y_em, y_Pm, y_Om, y_wm, y_im
    real*8, dimension(j) :: z_t0m, z_em, z_Pm, z_Om, z_wm, z_im
    
    real*8 :: np_Pp, nm_Pm
    real*8 :: ap_ms, ap_mp, ap_mm, ap_np, ap_Pp, am_mp, am_mm, am_nm, am_Pm
    real*8 :: comegap_Op, somegap_Op, cwp_wp, swp_wp, cip_ip, sip_ip
    real*8 :: comegam_Om, somegam_Om, cwm_wm, swm_wm, cim_im, sim_im
    real*8 :: mtot, mrp_mm, mrp_mp, mrm_mp, mrm_mm
    real*8, dimension(j) :: r_n, r_t0, r_e, r_a
    real*8, dimension(j) :: cosf_n, cosf_t0, cosf_e, cosf_a
    real*8, dimension(j) :: sinf_n, sinf_t0, sinf_e, sinf_a
    real*8, dimension(j) :: cosfw_wp, cosfw_np, cosfw_Pp, cosfw_t0p, cosfw_ep, cosfw_ap
    real*8, dimension(j) :: cosfw_wm, cosfw_nm, cosfw_Pm, cosfw_t0m, cosfw_em, cosfw_am
    real*8, dimension(j) :: sinfw_wp, sinfw_np, sinfw_Pp, sinfw_t0p, sinfw_ep, sinfw_ap
    real*8, dimension(j) :: sinfw_wm, sinfw_nm, sinfw_Pm, sinfw_t0m, sinfw_em, sinfw_am
    real*8, dimension(j) :: sinfw_ms, sinfw_mp, sinfw_mm, cosfw_ms, cosfw_mp, cosfw_mm
    real*8, dimension(j) :: r_ms, r_mp, r_mm, r_Pp, r_Pm
    
    
    np = 2 * pi / Pp
    nm = 2 * pi / Pm
    np_Pp = -np / Pp
    nm_Pm = -nm / Pm
    
    ap = (G * (ms + mp + mm) / (np ** 2.d0)) ** o3
    am = (G * (mp + mm) / (nm ** 2.d0)) ** o3
    ap_ms = G / (3 * (ap * np) ** 2.d0)
    ap_mp = ap_ms
    ap_mm = ap_ms
    ap_np = - 2 * ap / (3 * np)
    ap_Pp = ap_np * np_Pp
    am_mp = G / (3 * (am * nm) ** 2.d0)
    am_mm = am_mp
    am_nm = -2 * am / (3 * nm)
    am_Pm = am_nm * nm_Pm
    
    comegap = Cos(Op)
    somegap = Sin(Op)
    cwp = Cos(wp)
    swp = Sin(wp)
    cip = Cos(ip)
    sip = Sin(ip)
    comegap_Op = -somegap
    somegap_Op = comegap
    cwp_wp = -swp
    swp_wp = cwp
    cip_ip = -sip
    sip_ip = cip
    
    comegam = Cos(Om)
    somegam = Sin(Om)
    cwm = Cos(wm)
    swm = Sin(wm)
    cim = Cos(im)
    sim = Sin(im)
    comegam_Om = -somegam
    somegam_Om = comegam
    cwm_wm = -swm
    swm_wm = cwm
    cim_im = -sim
    sim_im = cim
    
    mtot = mp + mm
    mrp = -mm / mtot
    mrm = mp / mtot
    mrp_mm = -mrm / mtot
    mrp_mp = -mrp / mtot
    mrm_mm = mrm / mtot
    mrm_mp = mrp / mtot
    
    call solve_kepler_grad(t, np, t0p, ep, ap, r, cosf, sinf, &
                           r_n, r_t0, r_e, r_a, &
                           cosf_n, cosf_t0, cosf_e, cosf_a, & 
                           sinf_n, sinf_t0, sinf_e, sinf_a, j)
     
    cosfw = cwp * cosf - swp * sinf
    cosfw_wp = cwp_wp * cosf - swp_wp * sinf
    cosfw_np = cwp * cosf_n - swp * sinf_n
    !cosfw_Pp = cosfw_np * np_Pp
    cosfw_t0p = cwp * cosf_t0 - swp * sinf_t0
    cosfw_ep = cwp * cosf_e - swp * sinf_e
    cosfw_ap = cwp * cosf_a - swp * sinf_a
    cosfw_ms = cosfw_ap * ap_ms
    cosfw_mp = cosfw_ap * ap_mp
    cosfw_mm = cosfw_ap * ap_mm
    cosfw_Pp = cosfw_ap * ap_Pp + cosfw_np * np_Pp
    
    sinfw = swp * cosf + sinf * cwp
    sinfw_wp = swp_wp * cosf + sinf * cwp_wp
    sinfw_np = swp * cosf_n + sinf_n * cwp
    !sinfw_Pp = sinfw_np * np_Pp
    sinfw_t0p = swp * cosf_t0 + sinf_t0 * cwp
    sinfw_ep = swp * cosf_e + sinf_e * cwp
    sinfw_ap = swp * cosf_a + sinf_a * cwp
    sinfw_ms = sinfw_ap * ap_ms
    sinfw_mp = sinfw_ap * ap_mp
    sinfw_mm = sinfw_ap * ap_mm
    sinfw_Pp = sinfw_ap * ap_Pp + sinfw_np * np_Pp
    
    r_ms = r_a * ap_ms
    r_mp = r_a * ap_mp
    r_mm = r_a * ap_mm
    r_Pp = r_a * ap_Pp + r_n * np_Pp
    
    xbc = -r * (comegap * cosfw - somegap * sinfw * cip)
    xbc_wp = -r * (comegap * cosfw_wp - somegap * sinfw_wp * cip)
    xbc_Pp = - r_Pp * (comegap * cosfw - somegap * sinfw * cip) &
             - r * (comegap * cosfw_Pp - somegap * sinfw_Pp * cip)
    xbc_t0p = - r_t0 * (comegap * cosfw - somegap * sinfw * cip) &
              - r * (comegap * cosfw_t0p - somegap * sinfw_t0p * cip)
    xbc_ep = - r_e * (comegap * cosfw - somegap * sinfw * cip) &
             - r * (comegap * cosfw_ep - somegap * sinfw_ep * cip)
    xbc_ms = - r_ms * (comegap * cosfw - somegap * sinfw * cip) &
             - r * (comegap * cosfw_ms - somegap * sinfw_ms * cip)
    xbc_mp = - r_mp * (comegap * cosfw - somegap * sinfw * cip) &
             - r * (comegap * cosfw_mp - somegap * sinfw_mp * cip)
    xbc_mm = - r_mm * (comegap * cosfw - somegap * sinfw * cip) &
             - r * (comegap * cosfw_mm - somegap * sinfw_mm * cip)
    xbc_Op = - r * (comegap_Op * cosfw - somegap_Op * sinfw * cip)
    xbc_ip = r * somegap * sinfw * cip_ip
    
    ybc = -r * (somegap * cosfw + comegap * sinfw * cip)
    ybc_wp = -r * (somegap * cosfw_wp + comegap * sinfw_wp * cip)
    ybc_Pp = - r_Pp * (somegap * cosfw + comegap * sinfw * cip) &
             - r * (somegap * cosfw_Pp + comegap * sinfw_Pp * cip)
    ybc_t0p = - r_t0 * (somegap * cosfw + comegap * sinfw * cip) &
              - r * (somegap * cosfw_t0p + comegap * sinfw_t0p * cip)
    ybc_ep = - r_e * (somegap * cosfw + comegap * sinfw * cip) &
             - r * (somegap * cosfw_ep + comegap * sinfw_ep * cip)
    ybc_ms = - r_ms * (somegap * cosfw + comegap * sinfw * cip) &
             - r * (somegap * cosfw_ms + comegap * sinfw_ms * cip)
    ybc_mp = - r_mp * (somegap * cosfw + comegap * sinfw * cip) &
             - r * (somegap * cosfw_mp + comegap * sinfw_mp * cip)
    ybc_mm = - r_mm * (somegap * cosfw + comegap * sinfw * cip) &
             - r * (somegap * cosfw_mm + comegap * sinfw_mm * cip)
    ybc_Op = - r * (somegap_Op * cosfw + comegap_Op * sinfw * cip)
    ybc_ip = - r * comegap * sinfw * cip_ip
    
    zbc = r * sinfw * sip
    zbc_wp = r * sinfw_wp * sip
    zbc_Pp = r * sinfw_Pp * sip + r_Pp * sinfw * sip
    zbc_t0p = r * sinfw_t0p * sip + r_t0 * sinfw * sip
    zbc_ep = r * sinfw_ep * sip + r_e * sinfw * sip
    zbc_ms = r * sinfw_ms * sip + r_ms * sinfw * sip
    zbc_mm = r * sinfw_mm * sip + r_mm * sinfw * sip
    zbc_mp = r * sinfw_mp * sip + r_mp * sinfw * sip
    zbc_Op = 0.d0
    zbc_ip = r * sinfw * sip_ip
    
    call solve_kepler_grad(t, nm, t0m, em, am, r, cosf, sinf, &
                           r_n, r_t0, r_e, r_a, &
                           cosf_n, cosf_t0, cosf_e, cosf_a, &
                           sinf_n, sinf_t0, sinf_e, sinf_a, j)
    
    cosfw = cwm * cosf - swm * sinf
    cosfw_wm = cwm_wm * cosf - swm_wm * sinf
    cosfw_nm = cwm * cosf_n - swm * sinf_n
    !cosfw_Pm = cosfw_nm * nm_Pm
    cosfw_t0m = cwm * cosf_t0 - swm * sinf_t0
    cosfw_em = cwm * cosf_e - swm * sinf_e
    cosfw_am = cwm * cosf_a - swm * sinf_a
    cosfw_mp = cosfw_am * am_mp
    cosfw_mm = cosfw_am * am_mm
    cosfw_Pm = cosfw_am * am_Pm + cosfw_nm * nm_Pm
    
    sinfw = swm * cosf + sinf * cwm
    sinfw_wm = swm_wm * cosf + sinf * cwm_wm
    sinfw_nm = swm * cosf_n + sinf_n * cwm
    !sinfw_Pm = sinfw_nm * nm_Pm
    sinfw_t0m = swm * cosf_t0 + sinf_t0 * cwm
    sinfw_em = swm * cosf_e + sinf_e * cwm
    sinfw_am = swm * cosf_a + sinf_a * cwm
    sinfw_mp = sinfw_am * am_mp
    sinfw_mm = sinfw_am * am_mm
    sinfw_Pm = sinfw_am * am_Pm + sinfw_nm * nm_Pm
    
    r_mp = r_a * am_mp
    r_mm = r_a * am_mm
    r_Pm = r_a * am_Pm + r_n * nm_Pm
    
    x = -r * (comegam * cosfw - somegam * sinfw * cim)
    x_wm = -r * (comegam * cosfw_wm - somegam * sinfw_wm * cim)
    x_Pm = - r_Pm * (comegam * cosfw - somegam * sinfw * cim) &
             - r * (comegam * cosfw_Pm - somegam * sinfw_Pm * cim)
    x_t0m = - r_t0 * (comegam * cosfw - somegam * sinfw * cim) &
              - r * (comegam * cosfw_t0m - somegam * sinfw_t0m * cim)
    x_em = - r_e * (comegam * cosfw - somegam * sinfw * cim) &
             - r * (comegam * cosfw_em - somegam * sinfw_em * cim)
    x_mp = - r_mp * (comegam * cosfw - somegam * sinfw * cim) &
             - r * (comegap * cosfw_mp - somegap * sinfw_mp * cip)
    x_mm = - r_mm * (comegam * cosfw - somegam * sinfw * cim) &
             - r * (comegam * cosfw_mm - somegam * sinfw_mm * cim)
    x_Om = - r * (comegam_Om * cosfw - somegam_Om * sinfw * cim)
    x_im = r * somegam * sinfw * cim_im
    
    y = -r * (somegam * cosfw + comegam * sinfw * cim)
    y_wm = -r * (somegam * cosfw_wm + comegam * sinfw_wm * cim)
    y_Pm = - r_Pm * (somegam * cosfw + comegam * sinfw * cim) &
             - r * (somegam * cosfw_Pm + comegam * sinfw_Pm * cim)
    y_t0m = - r_t0 * (somegam * cosfw + comegam * sinfw * cim) &
              - r * (somegam * cosfw_t0m + comegam * sinfw_t0m * cim)
    y_em = - r_e * (somegam * cosfw + comegam * sinfw * cim) &
             - r * (somegam * cosfw_em + comegam * sinfw_em * cim)
    y_mp = - r_mp * (somegam * cosfw + comegam * sinfw * cim) &
             - r * (somegap * cosfw_mp + comegap * sinfw_mp * cip)
    y_mm = - r_mm * (somegam * cosfw + comegam * sinfw * cim) &
             - r * (somegam * cosfw_mm + comegam * sinfw_mm * cim)
    y_Om = - r * (somegam_Om * cosfw + comegam_Om * sinfw * cim)
    y_im = - r * comegam * sinfw * cim_im
    
    z = r * sinfw * Sin(im)
    z_wm = r * sinfw_wm * sim
    z_Pm = r * sinfw_Pm * sim + r_Pm * sinfw * sim
    z_t0m = r * sinfw_t0m * sim + r_t0 * sinfw * sim
    z_em = r * sinfw_em * sim + r_e * sinfw * sim
    z_mm = r * sinfw_mm * sim + r_mm * sinfw * sim
    z_mp = r * sinfw_mp * sim + r_mp * sinfw * sim
    z_Om = 0.d0
    z_im = r * sinfw * sim_im
    
    ! ms, t0p, ep, Pp, Op, wp, ip, mp, t0m, em, Pm, wm, Om, im, mm
    xp = xbc + x * mrp
    dxp(1, :) = xbc_ms
    dxp(2, :) = xbc_t0p
    dxp(3, :) = xbc_ep
    dxp(4, :) = xbc_Pp
    dxp(5, :) = xbc_Op
    dxp(6, :) = xbc_wp
    dxp(7, :) = xbc_ip
    dxp(8, :) = xbc_mp + x * mrp_mp
    dxp(9, :) = x_t0m * mrp
    dxp(10, :) = x_em * mrp
    dxp(11, :) = x_Pm * mrp
    dxp(12, :) = x_wm * mrp
    dxp(13, :) = x_Om * mrp
    dxp(14, :) = x_im * mrp
    dxp(15, :) = x_mm * mrp + x * mrp_mm
    
    yp = ybc + y * mrp
    dyp(1, :) = ybc_ms
    dyp(2, :) = ybc_t0p
    dyp(3, :) = ybc_ep
    dyp(4, :) = ybc_Pp
    dyp(5, :) = ybc_Op
    dyp(6, :) = ybc_wp
    dyp(7, :) = ybc_ip
    dyp(8, :) = ybc_mp + y * mrp_mp
    dyp(9, :) = y_t0m * mrp
    dyp(10, :) = y_em * mrp
    dyp(11, :) = y_Pm * mrp
    dyp(12, :) = y_wm * mrp
    dyp(13, :) = y_Om * mrp
    dyp(14, :) = y_im * mrp
    dyp(15, :) = y_mm * mrp + y * mrp_mm
    
    zp = zbc + z * mrp 
    dzp(1, :) = zbc_ms
    dzp(2, :) = zbc_t0p
    dzp(3, :) = zbc_ep
    dzp(4, :) = zbc_Pp
    dzp(5, :) = zbc_Op
    dzp(6, :) = zbc_wp
    dzp(7, :) = zbc_ip
    dzp(8, :) = zbc_mp + z * mrp_mp
    dzp(9, :) = z_t0m * mrp
    dzp(10, :) = z_em * mrp
    dzp(11, :) = z_Pm * mrp
    dzp(12, :) = z_wm * mrp
    dzp(13, :) = z_Om * mrp
    dzp(14, :) = z_im * mrp
    dzp(15, :) = z_mm * mrp + z * mrp_mm

    xm = xbc + x * mrm
    dxm(1, :) = xbc_ms
    dxm(2, :) = xbc_t0p
    dxm(3, :) = xbc_ep
    dxm(4, :) = xbc_Pp
    dxm(5, :) = xbc_Op
    dxm(6, :) = xbc_wp
    dxm(7, :) = xbc_ip
    dxm(8, :) = xbc_mp + x * mrm_mp
    dxm(9, :) = x_t0m * mrm
    dxm(10, :) = x_em * mrm
    dxm(11, :) = x_Pm * mrm
    dxm(12, :) = x_wm * mrm
    dxm(13, :) = x_Om * mrm
    dxm(14, :) = x_im * mrm
    dxm(15, :) = x_mm * mrm + x * mrm_mm
    
    ym = ybc + y * mrm
    dym(1, :) = ybc_ms
    dym(2, :) = ybc_t0p
    dym(3, :) = ybc_ep
    dym(4, :) = ybc_Pp
    dym(5, :) = ybc_Op
    dym(6, :) = ybc_wp
    dym(7, :) = ybc_ip
    dym(8, :) = ybc_mp + y * mrm_mp
    dym(9, :) = y_t0m * mrm
    dym(10, :) = y_em * mrm
    dym(11, :) = y_Pm * mrm
    dym(12, :) = y_wm * mrm
    dym(13, :) = y_Om * mrm
    dym(14, :) = y_im * mrm
    dyp(15, :) = y_mm * mrm + y * mrm_mm
    
    zm = zbc + z * mrm
    dzm(1, :) = zbc_ms
    dzm(2, :) = zbc_t0p
    dzm(3, :) = zbc_ep
    dzm(4, :) = zbc_Pp
    dzm(5, :) = zbc_Op
    dzm(6, :) = zbc_wp
    dzm(7, :) = zbc_ip
    dzm(8, :) = zbc_mp + z * mrm_mp
    dzm(9, :) = z_t0m * mrm
    dzm(10, :) = z_em * mrm
    dzm(11, :) = z_Pm * mrm
    dzm(12, :) = z_wm * mrm
    dzm(13, :) = z_Om * mrm
    dzm(14, :) = z_im * mrm
    dzm(15, :) = z_mm * mrm + z * mrm_mm
    
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
    
    ! ms, t0p, ep, Pp, Op, wp, ip, mp, t0m, em, Pm, wm, Om, im, mm
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

! compute gradiants for bp, bpm, theta from derivatives of cartesian components 
! returned by coords_grad. 
subroutine input_coords_grad(t, ms, t0p, ep, Pp, Op, wp, ip, mp, &
    &t0m, em, Pm, wm, Om, im, mm, j, &
    &bp, bpm, theta, dbp, dbpm, dtheta) bind(C, name="input_coords_grad")

    integer :: i
    real*8 :: a, b, c, tmp, mu
    real*8, dimension(j) :: bm, bm2, bp2, bpm2
    real*8, dimension(j) :: xp, yp, zp, xm, ym, zm
    real*8, dimension(j) :: bm_xm, bm_ym, bp_xp, bp_yp, bpm_xm, bpm_ym, bpm_xp, bpm_yp
    real*8, dimension(j) :: theta_bp, theta_bpm, theta_bm
    real*8, dimension(15, j) :: dxp, dyp, dzp, dxm, dym, dzm, dbm
    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: ms, t0p, ep, Pp, Op, wp, ip, mp 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im, mm
    real (c_double), bind(C), dimension(j) :: t
    real (c_double), bind(C), intent(out), dimension(j) :: bp, bpm, theta
    real (c_double), bind(C), intent(out), dimension(15, j) :: dbp, dbpm, dtheta
    
    call coords_grad(t, ms, t0p, ep, Pp, Op, wp, ip, mp, t0m, em, &
                Pm, wm, Om, im, mm, j, xp, yp, zp, xm, ym, zm, &
                dxp, dyp, dzp, dxm, dym, dzm)
    
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

    bm = bm * au_rsun
    bpm = bpm * au_rsun
    dbp = dbp * au_rsun
    dbpm = dbpm * au_rsun
    !dtheta = dtheta * au_rsun
            
end

end