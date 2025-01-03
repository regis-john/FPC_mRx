;ASA (27Feb2023): a stripped down version of "fpc_mms.pro",
; specifically looking at the Phan 2018 MRx interval
; specifically looking at the Wilder (2017,2018) MRx interval //25 July 2023
; CDF_TT2000, tdxs_moms, des_yr, des_dmo, des_ddy, des_dhr, des_dmn, des_dsc, des_dmilli, des_dmicros, des_dnanos, /BREAK
; des_yr[a], des_dmo[a], des_ddy[a], des_dhr[a], des_dmn[a], des_dsc[a], des_dmilli[a], des_dmicros[a], des_dnanos[a]

; Find derivative of array (for C'). Endpoints need to be calculated specially
function fpc_mms_2d_deriv, array, delta
  diff_arr = (shift(array, [-1,0]) - shift(array, [1,0])) / $
    (2d * delta)
  diff_arr[0,*] = (array[1,*] - array[0,*]) / delta
  diff_arr[-1,*] = (array[-1,*] - array[-2,*]) / delta

  return, diff_arr
end

; Interpolates but prevents extrapolation
function interp_no_extrap, v, x, xout, _extra = _ex
  int = interpol(v, x, xout, _extra = ex)
  wmin = where( xout lt min(x), nwmin)
  if nwmin gt 0 then int[wmin] = v[0]
  wmax = where( xout gt max(x), nwmax)
  if nwmax gt 0 then int[wmax] = v[-1]
  return, int
end

pro fpc_mrx, setupfile = setupfile, check = check, no_rem_zeros = no_rem_zeros, deltaf = deltaf, $
   mean_fac = mean_fac, ncorr = ncorr, zero_count = zero_count, no_Efilter = no_Efilter
  
  ;..................
;  setupfile="~/MMS_quick/Phan_data/fpc_mrx_Phan2018_setup.sav" ;Phan_data
;  setupfile="~/MMS_quick/Wilder_data/fpc_mrx_Wilder2018_setup_LMN.sav" ;Wilder_data
;  setupfile="/Users/aafshari/Documents/IDL_code/fpc_mrx/fpc_mrx_Wilder2018_setup_LMN.sav" ;Wilder_data
;  setupfile="/Users/aafshari/Documents/IDL_code/fpc_mrx/fpc_mrx_Wilder2018_setup_LMN_NoPlasmaFrame.sav" ;Wilder_data

;Wilder_data, e dist.s shifted to Rx geometry frame and E-field Lorentz transformed to Rx geometry frame
;setupfile="/Users/aafshari/Documents/IDL_code/fpc_mrx/fpc_mrx_Wilder2018_setup_LMN_ionPlasmaFrame.sav"  

;Wilder_data, e dist.s shifted to Rx geometry frame and E-field Lorentz transformed to Rx geometry frame
;velocity coord.s rotated to single LMN frame (vvLMN), and instantaneous FAC (vvfac)
  setupfile="/Users/aafshari/Documents/IDL_code/fpc_mrx/SWilder23_plots_and_codes/fpc_mrx_Wilder2018_23_setup_LMN_and_FAC.sav
  restore, setupfile
  
  ;..................
  ; Need to preserve these so they can be reset after loading the save file
  docheck = keyword_set(check)
  check = docheck
  
  if n_elements(mean_fac) gt 0 then do_mean_fac = mean_fac
  if n_elements(do_mean_fac) gt 0 then mean_fac = do_mean_fac

  remove_zeros = ~keyword_set(no_rem_zeros) ; FLAG: Removes zeros from mean calculation

  if n_elements(zero_count) eq 0 then zero_count = 0 ; Default is to replace 0s with NaN
  
  delta_f = keyword_set(deltaf) ; use delta_f in full f in correlation?
  
  if n_elements(ncorr) eq 0 then ncorr = 47  ; Correlation Interval for S02i ~ 7 seconds
  
  no_Efilter=1
  filter_E = ~keyword_set(no_Efilter) ; use high-pass filtered E_para in correlation?
  ;..................
  
  ;**********************************************************************************
  ; Generalizing variable names to match the species under consideration
  ;**********************************************************************************
  IF KEYWORD_SET(fpc_ions) && ~KEYWORD_SET(fpc_electrons) THEN BEGIN
    vt_mean = vti_mean
    time_slices = dis_time_slices
    dt = 0.150D
    energy_array = dis_energy_data[*,nti:ntf]
  ENDIF

  IF ~KEYWORD_SET(fpc_ions) && KEYWORD_SET(fpc_electrons) THEN BEGIN
    vt_mean = vte_mean
    time_slices = des_time_slices
    dt = 0.03D
    energy_array = des_energy_data[*,nti:ntf]
  ENDIF
  
  ;**********************************************************************************
  ; Selecting the DES time for the exact interval
  ;**********************************************************************************
  destime_INT=des_dist_time_data[nti_des:ntf_des]
;  a=2081
;  CDF_TT2000, des_dist_time_data, des_yr, des_dmo, des_ddy, des_dhr, des_dmn, des_dsc, des_dmilli, des_dmicros, des_dnanos, /BREAK
;  CDF_TT2000, destime_INT, des_yr, des_dmo, des_ddy, des_dhr, des_dmn, des_dsc, des_dmilli, des_dmicros, des_dnanos, /BREAK
;  print, des_yr[a], des_dmo[a], des_ddy[a], des_dhr[a], des_dmn[a], des_dsc[a], des_dmilli[a], des_dmicros[a], des_dnanos[a]
  
  
  ;******************************************************************************
  ; Set up the velocity space bins
  ;******************************************************************************
  ; Find the max velocity the FPI detectors detect in thermal speed units:
  vel_multiplier=DOUBLE(ROUND(MAX(vbin[-1,*]/vt_mean)/2))
  vel_multiplier=6D0
  ; Set velocity bin limits
  ; NOTE(25May2022) ASA: increasing vpar and vprp ranges covers more velocity-space, $
  ;   thus the density is closer to the "true" FPI value
  IF species EQ "i" THEN binfrac = 2.0D-01 ELSE binfrac = 1.0D-01 ;Binwidth as fraction of thermal velocity
  vparmin = -vel_multiplier*vt_mean ; Minimum vparallel (units of v_thermal)
  vparmax =  vel_multiplier*vt_mean ; Maximum vparallel (units of v_thermal)
  
  ;vprp2 is for 3V coordinates
  vprp2min =  -vel_multiplier*vt_mean ; Minimum vperp2 (units of v_thermal)
  vprp2max =  +vel_multiplier*vt_mean ; Maximum vperp2 (units of v_thermal)
  
  ; Set velocity bin size
  dv = DOUBLE(ROUND(binfrac*vt_mean))
  dv_int = FIX(dv) ; Used for plot titles
  ; Set vpar and vprp bin sizes
  dvpar = dv
  dvprp = dv
  
;  ; if we want to set the vperp and vpara binsizes separately
;  dvpar= DOUBLE(ROUND(binfrac*vtipara_mean))
;  dvprp= DOUBLE(ROUND(binfrac*vtiperp_mean))
  
  ; Create velocity bin arrays
  ; vpar: ranges from -vpar_min to +vpar_max
  vpar_bins = createxbins(vparmin,vparmax,dvpar) ; Midpoint (bin center) values
  nvpar = N_ELEMENTS(vpar_bins)
  ; vprp2: for 3V coordinates, ranges from -vprpmin to + vprp_max
  vprp2_bins = createxbins(vprp2min,vprp2max,dvprp) ; Midpoint (bin center) values
  nvprp2 = N_ELEMENTS(vprp2_bins)
  
  ;******************************************************************************
  ; High-pass filter E_parallel and E_perp
  ;******************************************************************************
  fcut_prp_label='0'

  IF (filter_E EQ 1) THEN BEGIN
    n = 5       ; nth order Butterworth Filter
    fcut = 1.0  ; 1.0 Hz cutoff
    fcut_prp = 1.0  ;Hz cutoff for E_perp
    fcut_prp_label='1'

    e_l_para_hipass = BUTTERWORTH_FILT(e_l_fac_dxs_mean[2,*], dt, n, fcut_prp)

    ; Replace the parallel electric field with the filtered
    e_l_fac_dxs_mean[2,*] = e_l_para_hipass

    ;..................
    ; High-pass filter the perpendicular components of the electric field
    ;..................
    eperp1hi=BUTTERWORTH_FILT(REFORM(e_l_fac_dxs_mean[0,*]),dt,n,fcut_prp)
    eperp2hi=BUTTERWORTH_FILT(REFORM(e_l_fac_dxs_mean[1,*]),dt,n,fcut_prp)
    ; put these into e_l_fac_dxs_mean for the FPC
    e_l_fac_dxs_mean[0,*] = eperp1hi ;E_prp1
    e_l_fac_dxs_mean[1,*] = eperp2hi ;E_prp2
  ENDIF
  
  ;******************************************************************************
  ; Compute the Mean VDF over the interval and Compute Delta fe
  ; defualt is delta_f=1, unless explicitly is no_deltaf=1
  ;******************************************************************************
  ; For Mean Interval FAC
  print, 'this is the minmax(ff_v) before removing zeros: ', minmax(ff_v)
  IF (mean_fac) THEN BEGIN ;----------------------------------------------------
    IF (remove_zeros) THEN BEGIN
      ;Save full fe and fe_err (with zeros) in temporary arrays
      tmp_ff_v = ff_v
      tmp_ff_v_err = ff_v_err

      ;NOTE(22April2022) ASA: still reasoning through the process of removing zeros.

      ;Find array values with zero
      izero = WHERE(particle_count_v LE zero_count, nizero)

      ;Set zeros to NaN in fe and fe_err
      if nizero gt 0 then ff_v[izero] = !VALUES.D_NAN
      if nizero gt 0 then ff_v_err[izero] = !VALUES.D_NAN
    ENDIF ; END remove_zeros-------------------------------------

    dff_v = dblarr(nbin, ntot)

    ;Compute mean
    mean_ff_v=MEAN(ff_v, DIMENSION=2, /nan)

    ;Put zeros back to compute delta fe
    IF (remove_zeros) THEN BEGIN
      iinf = WHERE(~FINITE(mean_ff_v), ninf)
      IF ninf GT 0 THEN mean_ff_v[iinf] = 0D0 ;NOTE(01June2022) ASA: this removal of zeroes and putting them back in is contentious
      ;NOTE(25August2022) ASA: contentious only in the tail of the distribution (since the zeroes in the tail are real), but not in the core.
      ff_v_err = TEMPORARY(tmp_ff_v_err)
      ff_v = TEMPORARY(tmp_ff_v)
    ENDIF

    ;Compute delta fe (if delta_f=1)
    IF (delta_f EQ 1) THEN BEGIN  ; Use Delta fe---------------------------
      mean_ff_v_arr = REBIN(mean_ff_v, [nbin,ntot])
      dff_v = ff_v - mean_ff_v_arr
      ff_v = dff_v
    ENDIF ELSE BEGIN  ; Don't use Delta fe---------------------------------
      ;Don't subtract mean fe, just set dfe=fe for correlation calculation
      dff_v=ff_v
    ENDELSE
  ENDIF else begin; END Mean VDF for mean_fac=1-------------------------------------
    ; For Instantaneous FAC
    ;TODO TO BE CODED!
    ;below is temporary code for MRx specific
    dff_v=ff_v
  endelse
  
  ;******************************************************************************
  ; Calculating the velocity-space volume elements (vol_elements) to be 
  ;   multiplied with the phase-space density to yield "velocity-space density 
  ;   distribution"
  ;******************************************************************************
  ;
  ; //NOTE(27Nov2024)ASA: have checked the density calculations using detector 
  ;   energy values uncorrected (i.e. the values directly from the data,
  ;   "energy_array"), the energy values corrected to the "instantaneous" scpot 
  ;   values ("encor"=energy_data[*,nti:ntf]-scpot_dxs_mean_array), and the
  ;   energy values corrected to the interval averaged scpot value
  ;   ("energy_data_corr"=energy_data[*,nti:ntf]-scpot_mean)
  ;   VERDICT: the "best" results (i.e. those closest to the "blessed" ion 
  ;   density values) are given by the uncorrected energy values
  ;   
  ; Calling vspace_vol.pro
  print, ''
  print, ' Running VSPACE_VOL.'
  print, ''
; vol_elements has units of (m/s)^3
  vol_elements=vspace_vol(energy_array, interleave, species)
  print, ''
  print, ' Finished VSPACE_VOL.'
  print, ''
  
;NOTE(24May2022) ASA: the vol_elements_arr below needs to be generalized for other intervals as well,
;  right now it's time-variable is specifically for S00i
;NOTE(31May2022) ASA: below is the update, needs to be checked; also, now need to take in to account interleave
IF (interleave EQ 0) THEN BEGIN
    vol_elements_arr=REBIN(vol_elements, [nbin,ntot]) ;//NOTE(23Jan2023)ASA:is REBIN doing what I think it's doing?
ENDIF ELSE BEGIN
  IF (ntot MOD 2 EQ 0) THEN BEGIN
    vol_elements_arr=REFORM(REBIN(vol_elements, [nbin,interleave+1,ntot/2]),nbin,ntot)
  ENDIF ELSE BEGIN
    vol_elements_arr=REFORM(REBIN(vol_elements, [nbin,interleave+1,(ntot-1)/(interleave+1)]),nbin,ntot-1)
    vol_elements_arr=[[vol_elements_arr],[vol_elements[*,0]]]
  ENDELSE
ENDELSE
print, ' Finished creating VOL_ELEMENTS_ARR.'
print, ''
  
  ;******************************************************************************
  ; Converting the phase space measurements (ff_v) to "velocity-space density
  ;   distribution" measurements (ff_density_nbin_ntot)
  ;******************************************************************************
  ; this is the particle density w.r.t. the integrated volume elements, $
  ;   a.k.a. "velocity-space density distribution"
  ; The final factor of 1D-6 converts from m^-3 to cm^-3
  ff_density_nbin_ntot=vol_elements_arr*ff_v*1D-6 ;dblarr(nbin,ntot)

  
  den1=TOTAL(ff_density_nbin_ntot,1,/nan)
  den1_arr=REBIN(REBIN(den1, [1,ntot]),[3,ntot])
  
  nni_fast_INT_arr=REBIN(REBIN(nni_fast_INT, [1,ntot]),[3,ntot])

  ;******************************************************************************
  ; E-field plotting
  ;******************************************************************************

;btime=[0:(n_elements(b_fac_dxs_mean[0,*])-1)*0.03:0.03]  
;pbl=plot(btime,b_fac_dxs_mean[0,*],'b',LAYOUT=[1,4,1])
;pbm=plot(btime,b_fac_dxs_mean[1,*],/over,'g')
;pbn=plot(btime,b_fac_dxs_mean[2,*],/over,'r')
;pbl.xrange=[0,2.6]
;;pbl.xtitle='Time: 05:03:55.5 + seconds'
;pbl.ytitle='B$_{LMN}$ (nT)'
;pline=plot(pbl.xrange,[0,0],/over)
;  
;  ff_density_3Vbinned_ntot=ff_3V_binned(ff_density_nbin_ntot, vvfac, vt_mean, dvpar, dvprp, vpar_bins, vprp2_bins)
;
; in the following line: 1D6 converts density (den2) from cm^-3 to m^-3; 1D3 converts velocity from km/s to m/s ; units: ampere/m^2 
  j_GSE_INT=!const.e*nni_fast_INT_arr*(vibulk_GSE_fast_INT-vebulk_GSE_INT)*1D6*1D3  
  jdote_GSE_INT=!const.e*nni_fast_INT_arr*(vibulk_GSE_fast_INT-vebulk_GSE_INT)*e_l_GSE_dxs_mean_INT*1D6 ; 1D6 converts density (nni_fast_INT_arr) from cm^-3 to m^-3
;                                                                                   ; 1D3 converts velocity from km/s to m/s 
;                                                                                   ; 1D-3 converts E-field from mV/m to V/m
;                                                                                   ; thus the mutiplier is 1D6
;                                                                                   ; units of energy density transfer rate (J m^-3 s^-1)
jdote_FAC_INT1=jdote_GSE_INT*0D0
j_FAC_INT1=j_GSE_INT*0D0
;NOTE(28Nov2024)ASA: the for-loop below is doing something funky in the rotation 
; so that the magnitude of jdote_FAC_INT1 is not equal to the magnitude of 
; jdote_FAC_INT2. BUT!! the magnitude of jdote_GSE_INT does equal the magnitude 
; of jdote_FAC_INT2. How peculiar!!
;NOTE(03Dec2024)ASA: |jdote_GSE_INT| = |jdote_FAC_INT1| ≠ |jdote_FAC_INT2|

FOR it=0, ntot-1 DO BEGIN
  j_FAC_INT1[*,it]=srot[*,*,it]#j_GSE_INT[*,it]
  jdote_FAC_INT1[*,it]=srot[*,*,it]#jdote_GSE_INT[*,it]
ENDFOR
  j_FAC_INT2=!const.e*nni_fast_INT_arr*(vibulk_FAC_fast_INT-vbulk_FAC_INT)*1D6*1D3
  jdote_FAC_INT2=!const.e*nni_fast_INT_arr*(vibulk_FAC_fast_INT-vbulk_FAC_INT)*e_l_FAC_dxs_mean_INT*1D6

;******************************************************************************
; Compute Instantaneous Correlations C'
; C'_El, C'_Em, and C'_En  in 3V (LMN) space
; C'_Eperp1, C'_Eperp2, and C'_Epara  in 3V (FAC) space
;******************************************************************************
; NOTE: FAC = (eprp1,eprp2,epar) = (0,1,2) in first index
; NOTE: LMN = (l,m,n) = (0,1,2) in first index
; NOTE on UNITS:
;         Converting vpar & vperp from km/s to m/s yield factor 1D3
;         Converting e_l_fac_dxs_mean_INT from mV/m to V/m yields 1D-3
;         TOGETHER, This is a factor of 1

; Vectorizing for speed -- reform vvfac, dff_v, and e_l_fac_dxs_mean_INT to [3, nbin, ntot]
  
; First vvfac
IF (mean_FAC) THEN BEGIN
  if interleave then begin ; Need to do odd/even separately if interleave on
    vvfac_arr = dblarr(3, nbin, ntot)
    neven = n_elements(vvfac_arr[0,0,0:*:2])
    even_ind = lindgen(neven) * 2
    nodd = n_elements(vvfac_arr[0,0,1:*:2])
    odd_ind = (lindgen(nodd) * 2) + 1l
  
    vvfac_arr[*,*,even_ind] = rebin(reform(vvfac[*,*,0], 3, nbin, 1), 3, nbin, neven)
    vvfac_arr[*,*,odd_ind] = rebin(reform(vvfac[*,*,1], 3, nbin, 1), 3, nbin, nodd)
  endif else vvfac_arr = REBIN(REFORM(vvfac,3,nbin,1),3,nbin,ntot) ; in this case, it should be [3, nbin, ntot]
ENDIF ELSE BEGIN
  vvFAC_arr=vvFAC
ENDELSE

; Vectorizing vvLMN from dblarr(3,nbin,interleave+1) to dblarr(3,nbin,ntot]
vvLMN_arr = dblarr(3, nbin, ntot)
neven = n_elements(vvLMN_arr[0,0,0:*:2])
even_ind = lindgen(neven) * 2
nodd = n_elements(vvLMN_arr[0,0,1:*:2])
odd_ind = (lindgen(nodd) * 2) + 1l

vvLMN_arr[*,*,even_ind] = rebin(reform(vvLMN[*,*,0], 3, nbin, 1), 3, nbin, neven)
vvLMN_arr[*,*,odd_ind] = rebin(reform(vvLMN[*,*,1], 3, nbin, 1), 3, nbin, nodd)

; dff_v
dff_v_r = rebin( reform(dff_v, 1, nbin, ntot), 3, nbin, ntot)
  
; full "velocity-space density distribution"
ff_density_3V_nbin_ntot=REBIN(REFORM(ff_density_nbin_ntot,[1,nbin,ntot]),3,nbin,ntot) ;this is in units of cm^-3

elfac_arr = rebin(reform(e_l_fac_dxs_mean_INT, 3, 1, ntot), 3, nbin, ntot)

  ;******************************************************************************
  ;Compute C'E in 3V FAC:
  ;    0:    C'_Eperp1 = q * vprp1 * ff * E_perp1
  ;    1:    C'_Eperp2 = q * vprp2 * ff * E_perp2
  ;    2:    C'_Epara = q * vpara * ff * E_para

  esign = species eq 'e' ? -1d : 1d

  ;vvfac_arr has units km/s
  ;ff_density_3V_nbin_ntot has units cm^-3
  ;elfac_arr has units mV/m
  ;the 1D6 at the end converts from J*s^-1*cm^-3 to J*s^-1*m^-3
  ; thus, cp_e_3v_den_FAC is returned in units of J*s^-1*m^-3 //energy density transfer rate
   cp_e_3v_den_FAC = esign *!Const.e*vvfac_arr*ff_density_3V_nbin_ntot*elfac_arr*1D6
  
  ;FAC coord.s
  ;C'El
  cp_Eperp1_den=REFORM(cp_e_3v_den_FAC[0,*,*])
  ;C'Em
  cp_Eperp2_den=REFORM(cp_e_3v_den_FAC[1,*,*])
  ;C'En
  cp_Epara_den=REFORM(cp_e_3v_den_FAC[2,*,*])
  undefine, cp_e_3v_den_FAC
  ;******************************************************************************
  
  ;******************************************************************************  
   ;Compute C'E in 3V LMN:
   ;    0:    C'_El = q * vprp1 * ff * E_L
   ;    1:    C'_Em = q * vprp2 * ff * E_M
   ;    2:    C'_En = q * vpar * delta_ff * E_N
   ;    
  ;similarly for creating cp_e_3V_den_LMN
;  cp_e_3v_den_LMN = esign *!Const.e*vvLMN_arr*ff_density_3V_nbin_ntot*elLMN_arr*1D6
;  
;  
;  ;LMN coord.s
;  ;C'El
;  cp_El_den=REFORM(cp_e_3v_den_FAC[0,*,*])
;  ;C'Em
;  cp_Em_den=REFORM(cp_e_3v_den_FAC[1,*,*])
;  ;C'En
;  cp_En_den=REFORM(cp_e_3v_den_FAC[2,*,*])
;  undefine, cp_e_3v_den_LMN
  ;******************************************************************************

  ;******************************************************************************
  ; 3V METHODOLOGY
  ;******************************************************************************
  ; First Method: bin results to 3V coordinates (vL,vM,vN) space (at all times)
  ;******************************************************************************
IF(0 EQ 1) THEN BEGIN
  ;try plotting the distribution itself, NOT the "velocity-space density distribution:
  ff_3V_t43=ff_3V_binned_mRx(ff_v[*,43],vvfac,vt_mean, dvpar, dvprp, vpar_bins, vprp2_bins)
  
  ;select the timeslice BEFORE integrating away the third v_chi-axis:
  ; choose ∆v_chi = 2*dv_chi centered about v_chi=0
  ff_3V_vperp1slice_t43=ff_3V_t43[59:60,*,*]
  ff_3V_vperp2slice_t43=ff_3V_t43[*,59:60,*]
  ff_3V_vparaslice_t43=ff_3V_t43[*,*,59:60]

  ;NOW integrate away the third v_chi-axis (which will result in minimal filling in of the central vacuous area in $
  ; the distribution plots:
  ff_vperp1vperp2_vparaslice_t43=TOTAL(ff_3V_vparaslice_t43,3)
  ff_vperp1vpara_vperp2slice_t43=TOTAL(ff_3V_vperp2slice_t43,2)
  ff_vperp2vpara_vperp1slice_t43=TOTAL(ff_3V_vperp1slice_t43,1)
  
  ;plot f(vperp1,vperp2,59:60;43)
  c=alog10(ff_vperp1vperp2_vparaslice_t43)
  minz=-24.
  maxz=-12.
  contourplot_prebinned,vprp2_bins,vprp2_bins,-2.*vt_mean,2.*vt_mean,-2.*vt_mean,2.*vt_mean,dvprp,dvprp,$
    'vperp1 (km/s)','vperp2 (km/s)',c,'ff_e(vperp1,vperp2,∆vpara,t = 43)',minz,maxz,'log(ffe)',0,'SWilder_ffe_vperp1vperp2_vparaslice_t43'

  ;plot f(vperp1,59:60,vpara;43)
  c=alog10(ff_vperp1vpara_vperp2slice_t43)
 minz=-24.
  maxz=-12.
  contourplot_prebinned,vprp2_bins,vprp2_bins,-2.*vt_mean,2.*vt_mean,-2.*vt_mean,2.*vt_mean,dvprp,dvprp,$
    'vpara (km/s)','vperp1 (km/s)',TRANSPOSE(c),'ff_e(vperp1,∆vperp2,vpara,t = 43)',minz,maxz,'log(ffe)',0,'SWilder_ffe_vparavperp1_vperp2slice_t43'

  ;plot f(59:60,vperp2,vpara;43)
  c=alog10(ff_vperp2vpara_vperp1slice_t43)
 minz=-24.
  maxz=-12.
  contourplot_prebinned,vprp2_bins,vprp2_bins,-2.*vt_mean,2.*vt_mean,-2.*vt_mean,2.*vt_mean,dvprp,dvprp,$
    'vpara (km/s)','vperp2 (km/s)',TRANSPOSE(c),'ff_e(∆vperp1,vperp2,vpara,t = 43)',minz,maxz,'log(ffe)',0,'SWilder_ffe_vparavperp2_vperp1slice_t43'

  ;------------------------------------------------------------------------------
  ;below used for APS 2024 presentation of just the distirbutions in LMN coord.s
  ;bin velocity space density distribution into 3V coord.s f_e(vL, vM, vN)
  ff_density_3V_ntot=ff_3V_binned_mRx(ff_density_nbin_ntot,vvfac,vt_mean, dvpar, dvprp, vpar_bins, vprp2_bins)

  ;Integrate away the vperp1 axis for all time:
  ff_density_vperp2vpara_ntot=TOTAL(ff_density_3V_ntot,1)

  ;Integrate away the vperp2 axis for all time:
  ff_density_vperp1vpara_ntot=TOTAL(ff_density_3V_ntot,2)

  ;Integrate away the vpara axis for all time:
  ff_density_vperp1vperp2_ntot=TOTAL(ff_density_3V_ntot,3)

  ;select only the timeslice corresponding to Wilder (2018) Figure 5:
  ff_density_vperp1vpara_t43=ff_density_vperp1vpara_ntot[*,*,43]
  ff_density_vperp2vpara_t43=ff_density_vperp2vpara_ntot[*,*,43]
  ff_density_vperp1vperp2_t43=ff_density_vperp1vperp2_ntot[*,*,43]

  ;select the timeslice BEFORE integrating away the third v_chi-axis:
  ; choose ∆v_chi = 2*dv_chi centered about v_chi=0
  ff_density_3V_vperp1slice_t43=ff_density_3V_ntot[59:60,*,*,43]
  ff_density_3V_vperp2slice_t43=ff_density_3V_ntot[*,59:60,*,43]
  ff_density_3V_vparaslice_t43=ff_density_3V_ntot[*,*,59:60,43]

  ;NOW integrate away the third v_chi-axis (which will result in minimal filling in of the central vacuous area in $
  ; the distribution plots:
  ff_density_vperp1vperp2_vparaslice_t43=TOTAL(ff_density_3V_vparaslice_t43,3)
  ff_density_vperp1vpara_vperp2slice_t43=TOTAL(ff_density_3V_vperp2slice_t43,2)
  ff_density_vperp2vpara_vperp1slice_t43=TOTAL(ff_density_3V_vperp1slice_t43,1)

  ;plot f(vperp1,vperp2,59:60;43)
  c=alog10(ff_density_vperp1vperp2_vparaslice_t43)
  minz=-5.5
  maxz=0.075
  contourplot_prebinned,vprp2_bins,vprp2_bins,-2.*vt_mean,2.*vt_mean,-2.*vt_mean,2.*vt_mean,dvprp,dvprp,$
    'vperp1 (km/s)','vperp2 (km/s)',c,'f_e(vperp1,vperp2,∆vpara,t = 43)',minz,maxz,'log(f(cm!E-3!N))',0,'SWilder_fe_vperp1vperp2_vparaslice_t43'

  ;plot f(vperp1,59:60,vpara;43)
  c=alog10(ff_density_vperp1vpara_vperp2slice_t43)
  minz=-5.5
  maxz=0.075
  contourplot_prebinned,vprp2_bins,vprp2_bins,-2.*vt_mean,2.*vt_mean,-2.*vt_mean,2.*vt_mean,dvprp,dvprp,$
    'vpara (km/s)','vperp1 (km/s)',TRANSPOSE(c),'f_e(vperp1,∆vperp2,vpara,t = 43)',minz,maxz,'log(f(cm!E-3!N))',0,'SWilder_fe_vparavperp1_vperp2slice_t43'

  ;plot f(59:60,vperp2,vpara;43)
  c=alog10(ff_density_vperp2vpara_vperp1slice_t43)
  minz=-5.5
  maxz=0.075
  contourplot_prebinned,vprp2_bins,vprp2_bins,-2.*vt_mean,2.*vt_mean,-2.*vt_mean,2.*vt_mean,dvprp,dvprp,$
    'vpara (km/s)','vperp2 (km/s)',TRANSPOSE(c),'f_e(∆vperp1,vperp2,vpara,t = 43)',minz,maxz,'log(f(cm!E-3!N))',0,'SWilder_fe_vparavperp2_vperp1slice_t43'


  STOP

ENDIF

  ;------------------------------------------------------------------------------
  ;bin velocity space density distribution into 3V coord.s f_e(vperp1, vperp2, vpara)
  ff_density_3V_ntot=ff_3V_binned_mRx(ff_density_nbin_ntot,vvfac,vt_mean, dvpar, dvprp, vpar_bins, vprp2_bins)
  
  ;Integrate away the vperp1 axis for all time:
  ff_density_vperp2vpara_ntot=TOTAL(ff_density_3V_ntot,1)
  
  ;Integrate away the vperp2 axis for all time:
  ff_density_vperp1vpara_ntot=TOTAL(ff_density_3V_ntot,2)
  
  ;Integrate away the vpara axis for all time:
  ff_density_vperp1vperp2_ntot=TOTAL(ff_density_3V_ntot,3)
  
  ;select only the timeslice corresponding to Wilder (2018) Figure 5:
  ff_density_vperp1vpara_t26=ff_density_vperp1vpara_ntot[*,*,26]
  ff_density_vperp2vpara_t26=ff_density_vperp2vpara_ntot[*,*,26]
  ff_density_vperp1vperp2_t26=ff_density_vperp1vperp2_ntot[*,*,26]
  
  ;plot f(vperp1,vperp2)
  c=alog10(ff_density_vperp1vperp2_t26)
  minz=-5.5
  maxz=0.075
;  contourplot_prebinned,vprp2_bins,vprp2_bins,-2.*vt_mean,2.*vt_mean,-2.*vt_mean,2.*vt_mean,dvprp,dvprp,$
;    'vperp1 (km/s)','vperp2 (km/s)',c,'f_e(vperp1,vperp2) t = 26',minz,maxz,'log(f(cm!E-3!N))',0,'SWilder_26_fe_vperp1vperp2_t43'
    
  ;plot f(vpara,vperp1)
  c=alog10(ff_density_vperp1vpara_t26)
  minz=-5.5
  maxz=0.075
;  contourplot_prebinned,vpar_bins,vprp2_bins,-2.*vt_mean,2.*vt_mean,-2.*vt_mean,2.*vt_mean,dvprp,dvprp,$
;    'vpara (km/s)','vperp1 (km/s)',TRANSPOSE(c),'f_e(vpara,vperp1) t = 26',minz,maxz,'log(f(cm!E-3!N))',0,'SWilder_26_fe_vparavperp1_t43'
    
  ;plot f(vpara,vperp2)
  c=alog10(ff_density_vperp2vpara_t26)
  minz=-5.5
  maxz=0.075
;  contourplot_prebinned,vpar_bins,vprp2_bins,-2.*vt_mean,2.*vt_mean,-2.*vt_mean,2.*vt_mean,dvprp,dvprp,$
;    'vpara (km/s)','vperp2 (km/s)',TRANSPOSE(c),'f_e(vpara,vperp2) t = 26',minz,maxz,'log(f(cm!E-3!N))',0,'SWilder_26_fe_vparavperp2_t43'

;--------------
;the following block plots each timeslice of f(v_Chi1, v_Chi2) using contourplot_prebinned
;  minz=-5.5
;  maxz=0.075
  xtit_orig= 'v!DL!N (km/s)'
  ytit_orig= 'v!DM!N (km/s)'
;  
;  FOR i=0, ntot-1 DO BEGIN
;    c=alog10(ff_density_vLvN_ntot[*,*,i])
;    contourplot_prebinned,vpar_bins,vprp2_bins,-1D4,1D4,-1D4,1D4,dvprp,dvprp, $
;      xtit_orig,ytit_orig,c,'f_e(vL,vM) t = '+STRCOMPRESS(i,/REMOVE_ALL),minz,maxz,'log(f(cm!E-3!N))',0,'SWilder_f_e_vLvM_t'+STRCOMPRESS(i,/REMOVE_ALL)
;  ENDFOR
  ;--------------
  ;above used for APS 2024 presentation of just the distirbutions in LMN coord.s
  ;------------------------------------------------------------------------------
  ;--------------
  ;binning C' in 3V
  ;--------------
  ;binning C'_Eperp1
  cp_Eperp1_3Vbinned_ntot=ff_3V_binned_mRx(cp_Eperp1_den, vvfac, vt_mean, dvpar, dvprp, vpar_bins, vprp2_bins)
  ;binning C'_Eperp2
  cp_Eperp2_3Vbinned_ntot=ff_3V_binned_mRx(cp_Eperp2_den, vvfac, vt_mean, dvpar, dvprp, vpar_bins, vprp2_bins)
  ;binning C'_Epara
  cp_Epara_3Vbinned_ntot=ff_3V_binned_mRx(cp_Epara_den, vvfac, vt_mean, dvpar, dvprp, vpar_bins, vprp2_bins)
  
  ; extract times close to the EDR //according to Wilder
  cp_Eperp1_3Vbinned_times=cp_Eperp1_3Vbinned_ntot;[*,*,*,39:43]
  cp_Eperp2_3Vbinned_times=cp_Eperp2_3Vbinned_ntot;[*,*,*,39:43] ; the time index is 44 to correspend to Wilder (2018) JGR Fig.5: 2015-12-09 05:03:56.932592
  cp_Epara_3Vbinned_times=cp_Epara_3Vbinned_ntot;[*,*,*,39:43]
  ; initializing C_{E_chi} arrays
  c_Eperp1_3Vbinned_times=cp_Eperp1_3Vbinned_times*0.0
  c_Eperp2_3Vbinned_times=cp_Eperp2_3Vbinned_times*0.0
  c_Epara_3Vbinned_times=cp_Epara_3Vbinned_times*0.0
  ; initializing the vpar_bins array
  vprp2_bins_arr=REBIN(vprp2_bins,[nvprp2,nvprp2])

  ;--------------
  ; Derivatives of C'_Echi and computing C_Echi
  ;--------------
  FOR intime=0,n_elements(cp_Eperp1_3Vbinned_times[0,0,0,*])-1 DO BEGIN
    FOR invel=0,nvprp2-1 DO BEGIN
      cp_Eperp1_slice=REFORM(cp_Eperp1_3Vbinned_times[*,*,invel,intime])
      cp_Eperp2_slice=REFORM(cp_Eperp2_3Vbinned_times[invel,*,*,intime])
      cp_Epara_slice=TRANSPOSE(REFORM(cp_Epara_3Vbinned_times[invel,*,*,intime])) ;C'_En(vN,vM)
      ; Derivatives of C'_Echi w.r.t. v_chi in 2V coordinates
      diff_cp_Eperp1_2Vbinned_t=fpc_mms_2d_deriv(cp_Eperp1_slice, dvprp * 1D3) ;dC'_El(vL,vM)/dvL
      diff_cp_Eperp2_2Vbinned_t=fpc_mms_2d_deriv(cp_Eperp2_slice, dvprp * 1D3) ;dC'_Em(vM,vN)/dvM
      diff_cp_Epara_2Vbinned_t=TRANSPOSE(fpc_mms_2d_deriv(cp_Epara_slice, dvprp * 1D3)) ;dC'_En(vM,vN)/dvN
      ;..................
      ; Calculation of C_El(vL,vM,vN)
      c_Eperp1_3Vbinned_times[*,*,invel,intime]=(-1D0 * vprp2_bins_arr * 1D3/2 * diff_cp_Eperp1_2Vbinned_t) + (cp_Eperp1_slice / 2D0)
      ; Calculation of C_Em(vL,vM,vN)
      c_Eperp2_3Vbinned_times[invel,*,*,intime]=(-1D0 * vprp2_bins_arr * 1D3/2 * diff_cp_Eperp2_2Vbinned_t) + (cp_Eperp2_slice / 2D0)
      ; Calculation of C_En(vL,vM,vN)
      c_Epara_3Vbinned_times[invel,*,*,intime]=(-1D0 * vprp2_bins_arr * 1D3/2 * diff_cp_Epara_2Vbinned_t) + (cp_Epara_slice / 2D0)
      ;..................
    ENDFOR
  ENDFOR

  ceperp1=total(c_Eperp1_3Vbinned_times,1)
  ceperp1=total(ceperp1,1)
  ceperp1=total(ceperp1,1)
  
  ceperp2=total(c_Eperp2_3Vbinned_times,1)
  ceperp2=total(ceperp2,1)
  ceperp2=total(ceperp2,1)
  
  cepara=total(c_Epara_3Vbinned_times,1)
  cepara=total(cepara,1)
  cepara=total(cepara,1)
  
;  peperp1=plot(ceperp1)
;  peperp2=plot(ceperp2,/over,'r')
;  pepara=plot(cepara,/over,'b')
  
  ;--------------
  ; Manipulations of C_Epara
  ;--------------
    ;Integrate away the vperp1-axis for the selected times:
  c_Epara_vperp2vpara_ntot=TOTAL(c_Epara_3Vbinned_times,1)        ;C_Epara(vperp2,vpara;t)
  
  ;Integrate away the vperp2-axis for the selected times:
  c_Epara_vpara_ntot=TOTAL(c_Epara_vperp2vpara_ntot,1)            ;C_Epara(vpara;t)
  c_Epara_vperp1vpara_ntot=TOTAL(c_Epara_3Vbinned_times,2)        ;C_Epara(vperp1,vpara;t)
  
  ;Integrate away the vpara-axis for the selected times:
  c_Epara_vperp1_ntot=TOTAL(c_Epara_vperp1vpara_ntot,2)           ;C_Epara(vperp1;t)
  c_Epara_vperp2_ntot=TOTAL(c_Epara_vperp2vpara_ntot,2)           ;C_Epara(vperp2;t)

  ;Integrate away the vpara-axis:
  c_Epara_ntot=TOTAL(c_Epara_vpara_ntot,1)                        ;C_Epara(t) = (J.E')_\|

  ;Time average of C_Epara(vpara;t) to get <C_Epara(vpara)>_\tau
  c_Epara_vpara_timeavg=MEAN(c_Epara_vpara_ntot,dim=2)
  
  vt_mean_floor=floor(vt_mean)
;NOTE(25May2022) ASA: when using contourplot_prebinned, the vpar_bins and vprp_bins need to be normalized by vt_mean_floor so that the spacings are even, 
;  else vt_mean has some decimal points that cause contourplot_pribinned to truncate the edge bins in both vpar & vprp.
;  Also, xrange and yrange could be adjusted, but it elongates the plot in one direction making the square bins not look square

  e_index_seconds=indgen(ntot)*0.03
  
  ;below is the timestack plot of just C_Epara(vpara;t); use color table 33 and remove REVERSE(*_orig)
  print,minmax(c_Epara_vpara_ntot)
  d=c_Epara_vpara_ntot
  c=alog10(d)
  check_multiplier=DOUBLE(MIN(ABS(CEIL(ABS(MINMAX(c,/nan))))))
  d=d*10^check_multiplier
  minmaxz=minmax(d)
  minmaxz=ROUND(MEAN(ABS(minmaxz)))
  minz=-FLOAT(minmaxz)
  maxz=FLOAT(minmaxz)
  
  cbar_title='C_Epara(vpara;t) x 10^-'+STRCOMPRESS(FIX(check_multiplier),/REMOVE_ALL)+' (W/m^3)'

    
  contourplot_prebinned,vpar_bins,e_index_seconds,-2.*vt_mean_floor,2.*vt_mean_floor,0.,max(e_index_seconds),dvpar,dt, $
     'vpara (km/s)','t (s)',d,'C_Epara(vpara;t)',minz,maxz,cbar_title,1,'SWilder_23_C_Epara_vpara_timestack'
     
  ;--------------   

  ;--------------
  ; Manipulations of C_Eperp1
  ;--------------
  ;Integrate away the vperp2-axis for the selected times:
  c_Eperp1_vperp1vpara_ntot=TOTAL(c_Eperp1_3Vbinned_times,2)        ;C_Eperp1(vperp1,vpara;t)

  ;Integrate away the vpara-axis for the selected times:
  c_Eperp1_vperp1_ntot=TOTAL(c_Eperp1_vperp1vpara_ntot,2)           ;C_Eperp1(vperp1;t)

  ;Integrate away the vperp1-axis:
  c_Eperp1_ntot=TOTAL(c_Eperp1_vperp1_ntot,1)                       ;C_Eperp1(t) = (J.E')_\bot1

  ;Time average of C_Eperp1(vpara;t) to get <C_Eperp1(vperp1)>_\tau
  c_Eperp1_vperp1_timeavg=MEAN(c_Eperp1_vperp1_ntot,dim=2)

  ;below is the timestack plot of just C_Eperp1(vperp1;t); use color table 33 and remove REVERSE(*_orig)
  print,minmax(c_Eperp1_vperp1_ntot)
  d=c_Eperp1_vperp1_ntot
  c=alog10(d)
;  check_multiplier=DOUBLE(MIN(ABS(CEIL(ABS(MINMAX(c,/nan))))))
  check_multiplier=9. ;to match the multiplier of C_Epara(vpara;t) and put them on the same scale.
  d=d*10^check_multiplier
  minmaxz=minmax(d)
  minmaxz=ROUND(MEAN(ABS(minmaxz)))
  minz=-6. ;both minz and maxz are set to the same values as that of C_Epara(vpara;t) to put them all on the same scale.
  maxz=6.

  cbar_title='C_Eperp1(vperp1;t) x 10^-'+STRCOMPRESS(FIX(check_multiplier),/REMOVE_ALL)+' (W/m^3)'

  contourplot_prebinned,vprp2_bins,e_index_seconds,-2.*vt_mean_floor,2.*vt_mean_floor,0.,max(e_index_seconds),dvprp,dt, $
      'vperp1 (km/s)','t (s)',d,'C_Eperp1(vperp1;t)',minz,maxz,cbar_title,1,'SWilder_23_C_Eperp1_vperp1_timestack'

  ;--------------
 
  ;--------------
  ; Manipulations of C_Eperp2
  ;--------------
  ;Integrate away the vperp1-axis for the selected times:
  c_Eperp2_vperp2vpara_ntot=TOTAL(c_Eperp2_3Vbinned_times,2)        ;C_Eperp2(vperp2,vpara;t)

  ;Integrate away the vpara-axis for the selected times:
  c_Eperp2_vperp2_ntot=TOTAL(c_Eperp2_vperp2vpara_ntot,2)           ;C_Eperp2(vperp2;t)
  
  ;Integrate away the vperp2-axis:
  c_Eperp2_ntot=TOTAL(c_Eperp2_vperp2_ntot,1)                       ;C_Eperp2(t) = (J.E')_\bot2

  ;Time average of C_Epara(vpara;t) to get <C_Epara(vpara)>_\tau
  c_Eperp2_vperp2_timeavg=MEAN(c_Eperp2_vperp2_ntot,dim=2)

  ;below is the timestack plot of just C_Eperp2(vperp2;t); use color table 33 and remove REVERSE(*_orig)
  print,minmax(c_Eperp2_vperp2_ntot)
  d=c_Eperp2_vperp2_ntot
  c=alog10(d)
;  check_multiplier=DOUBLE(MIN(ABS(CEIL(ABS(MINMAX(c,/nan))))))
  check_multiplier=9. ;to match the multiplier of C_Epara(vpara;t) and put them on the same scale.
  d=d*10^check_multiplier
  minmaxz=minmax(d)
  minmaxz=ROUND(MEAN(ABS(minmaxz)))
  minz=-6. ;both minz and maxz are set to the same values as that of C_Epara(vpara;t) to put them all on the same scale.
  maxz=6.


  cbar_title='C_Eperp2(vperp2;t) x 10^-'+STRCOMPRESS(FIX(check_multiplier),/REMOVE_ALL)+' (W/m^3)'

  contourplot_prebinned,vprp2_bins,e_index_seconds,-2.*vt_mean_floor,2.*vt_mean_floor,0.,max(e_index_seconds),dvprp,dt,$
      'vperp2 (km/s)','t (s)',d,'C_Eperp2(vperp2;t)',minz,maxz,cbar_title,1,'SWilder_23_C_Eperp2_vperp2_timestack'
 
 
  print, STRING(7B)
  stop

    ;--------------
  ;the following block plots each timeslice of C_Em(v_Chi1, v_Chi2) using contourplot_prebinned
    minz=-2.0
    maxz=2.0
  xtit_orig= 'v!DL!N (km/s)'
  ytit_orig= 'v!DN!N (km/s)'
  
    FOR i=0, n_elements(c_Em_vLvM_ntot[0,0,*])-1 DO BEGIN
      d=c_Em_vLvM_ntot[*,*,i]
      c=alog10(d)
      check_multiplier=DOUBLE(MIN(ABS(CEIL(ABS(MINMAX(c,/nan))))))
      d=d*10^check_multiplier
      minmaxz=minmax(d)
      minmaxz=ROUND(MEAN(ABS(minmaxz)))
      minz=-2.
      maxz=2.
      print, minz, maxz
      contourplot_prebinned,vprp2_bins,vprp2_bins,-1D4,1D4,-1D4,1D4,dvprp,dvprp, $
        xtit_orig,ytit_orig,d,'C_Em(vL,vM) t = '+STRCOMPRESS(i+40,/REMOVE_ALL),minz,maxz,'C_Em(vL,vM)',1,'SWilder_C_Em_vLvM_t'+STRCOMPRESS(i+40,/REMOVE_ALL)
    ENDFOR
  ;--------------
 stop
  
  ; integrate away the vN axis
  c_Em_vLvM_t=TOTAL(c_Em_3Vbinned_t,3)

  ;below is the plot of just the background phase-space distribution; use color table 33 and remove REVERSE(*_orig)
  print,minmax(c_Em_vLvM_t)
  d=c_Em_vLvM_t
  c=alog10(d)
  check_multiplier=DOUBLE(MIN(ABS(CEIL(ABS(MINMAX(c,/nan))))))
  d=d*10^check_multiplier
  minz=-5
  maxz=5
  contourplot_prebinned,vpar_bins,vprp2_bins,-1D4,1D4,-1D4,1D4,dvpar,dvprp, $
    xtit_orig,ytit_orig,d,'C_Em(vL,vM) t = 14',minz,maxz,'C_Em(vL,vM)',1,'SWilder_C_E_M_vLvM_t14'
  print, STRING(7B)
  stop

  ; integrate away the vM axis
  c_Em_vLvN_t=TOTAL(c_Em_3Vbinned_t,2)

  ;below is the plot of just the background phase-space distribution; use color table 33 and remove REVERSE(*_orig)
  print,minmax(c_Em_vLvN_t)
  d=c_Em_vLvN_t
  c=alog10(d)
  check_multiplier=DOUBLE(MIN(CEIL(ABS(MINMAX(c)))))
  d=d*10^check_multiplier
  minz=-4.
  maxz=4.0
  contourplot_prebinned,vpar_bins,vprp2_bins,-4.*vte_mean,4.*vte_mean,-4.*vte_mean,4.*vte_mean,dvpar,dvprp,'vL (km/s)','vN (km/s)',d,'Wilder(2018)JGR Fig.5',minz,maxz,'C_E_M(vL,vN)',1,'SWilder_C_E_M_vLvN'
  print, STRING(7B)
;  stop

  ; integrate away the vL axis
  c_EM_vMvN_t=TOTAL(c_EM_3Vbinned_t,1)

  ;below is the plot of just the background phase-space distribution; use color table 33 and remove REVERSE(*_orig)
  print,minmax(c_EM_vMvN_t)
  d=c_EM_vMvN_t
  c=alog10(d)
  check_multiplier=DOUBLE(MIN(CEIL(ABS(MINMAX(c)))))
  d=d*10^check_multiplier
  minz=-1.25
  maxz=1.25
  contourplot_prebinned,vprp2_bins,vprp2_bins,-4.*vte_mean,4.*vte_mean,-4.*vte_mean,4.*vte_mean,dvpar,dvprp,'vM (km/s)','vN (km/s)',d,'Wilder(2018)JGR Fig.5',minz,maxz,'C_E_M(vM,vN)',1,'SWilder_C_E_M_vMvN'
  ;--------------
  print, STRING(7B)
  stop
  
  ;--------------
  ;binning C'_El in 3V
  ;--------------
  cp_El_3Vbinned_ntot=ff_3V_binned(cp_El_den, vvfac, vt_mean, dvpar, dvprp, vpar_bins, vprp2_bins)
  ; extract one specific time
  cp_El_3Vbinned_t=cp_El_3Vbinned_ntot[*,*,*,14] ; the time index is 14 to correspend to Wilder (2018) JGR Fig.5: 2015-12-09 05:03:56.932592
  ;..................
  ; initializing C_{E_L} arrays
  c_El_3Vbinned_t=cp_El_3Vbinned_t*0.0
  ; and vpar_bins array
  vpar_bins_arr=REBIN(vpar_bins,[nvpar,nvpar])
  ;--------------
  ; Derivatives of C'_El and computing C_El
  ;--------------
  FOR invN=0,nvprp2-1 DO BEGIN
    cp_El_slice=REFORM(cp_El_3Vbinned_t[*,*,invN])
    ; Derivatives of C' w.r.t. vL in 2V coordinates
    diff_cp_El_2Vbinned_t=fpc_mms_2d_deriv(cp_El_slice, dvpar * 1D3) ;dC'_Eprp1(vL,vM)/dvL
    ;..................
    ; Calculation of C_El(vL,vM,vN)
    c_El_3Vbinned_t[invN,*,*]=(-1D0 * vpar_bins_arr * 1D3/2 * diff_cp_El_2Vbinned_t) + (cp_El_slice / 2D0)
    ;..................
  ENDFOR
  print, STRING(7B)
  
  
  ; integrate away the vN axis
  c_El_vLvM_t=TOTAL(c_El_3Vbinned_t,3)
  
  ;below is the plot of just the background phase-space distribution; use color table 33 and remove REVERSE(*_orig)
  print,minmax(c_El_vLvM_t)
  d=c_El_vLvM_t
  c=alog10(d)
  check_multiplier=DOUBLE(MIN(CEIL(ABS(MINMAX(c)))))
  d=d*10^check_multiplier
  minz=-4.
  maxz=4.0
  contourplot_prebinned,vpar_bins,vprp2_bins,-4.*vte_mean,4.*vte_mean,-4.*vte_mean,4.*vte_mean,dvpar,dvprp,'vL (km/s)','vM (km/s)',d,'Wilder(2018)JGR Fig.5',minz,maxz,'C_E_L(vL,vM)',1,'SWilder_C_E_L_vLvM'
  print, STRING(7B)
  stop
  
  ; integrate away the vM axis
  c_El_vLvN_t=TOTAL(c_El_3Vbinned_t,2)

  ;below is the plot of just the background phase-space distribution; use color table 33 and remove REVERSE(*_orig)
  print,minmax(c_El_vLvN_t)
  d=c_El_vLvN_t
  c=alog10(d)
  check_multiplier=DOUBLE(MIN(CEIL(ABS(MINMAX(c)))))
  d=d*10^check_multiplier
  minz=-5.
  maxz=5.0
  contourplot_prebinned,vpar_bins,vprp2_bins,-4.*vte_mean,4.*vte_mean,-4.*vte_mean,4.*vte_mean,dvpar,dvprp,'vL (km/s)','vN (km/s)',d,'Wilder(2018)JGR Fig.5',minz,maxz,'C_E_L(vL,vN)',1,'SWilder_C_E_L_vLvN'
  print, STRING(7B)
  stop
  
  ; integrate away the vL axis
  c_El_vMvN_t=TOTAL(c_El_3Vbinned_t,1)

  ;below is the plot of just the background phase-space distribution; use color table 33 and remove REVERSE(*_orig)
  print,minmax(c_El_vMvN_t)
  d=c_El_vMvN_t
  c=alog10(d)
  check_multiplier=DOUBLE(MIN(CEIL(ABS(MINMAX(c)))))
  d=d*10^check_multiplier
  minz=-5.
  maxz=5.0
  contourplot_prebinned,vpar_bins,vprp2_bins,-4.*vte_mean,4.*vte_mean,-4.*vte_mean,4.*vte_mean,dvpar,dvprp,'vM (km/s)','vN (km/s)',d,'Wilder(2018)JGR Fig.5',minz,maxz,'C_E_L(vM,vN)',1,'SWilder_C_E_L_vMvN'
  ;--------------
  print, STRING(7B)
  stop

  ;--------------
  ;binning C'_En in 3V
  ;--------------
  cp_En_3Vbinned_ntot=ff_3V_binned(cp_En_den, vvfac, vt_mean, dvpar, dvprp, vpar_bins, vprp2_bins)
  stop
  ;..................
  ; a quick check if ff_density_ntot matches the density from FPI moments
  ff_density_ntot=TOTAL(ff_density_nbin_ntot,1,/nan) ;units of cm^-3
  print, STRING(7B)
  STOP
  STOP

 ;..................
  STOP
END