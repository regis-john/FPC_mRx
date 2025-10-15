;ASA (23Feb2023): a stripped down version of "fpc_mms.pro", 
; specifically for looking at the Phan 2018 MRx interval

; Tests if all elements of an array are "sufficiently close" in value along 2nd dim.
function is_2darr_rep, array
  array0 = array[*,0]
  ; Get dimension sizes
  dimsz = size(array, /dimensions)
  nelem0 = dimsz[0]
  nelem1 = dimsz[1]
  ; Reform, replicate
  testarray = transpose(rebin( reform(array0, 1, nelem0), nelem1, nelem0))
  check = total(abs( array - testarray) le 1e-6) eq n_elements(array) ? 1 : 0
  return, check
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

pro fpc_mrx_setup
;------------------------------------------------------------------------------
;--------------------------
; For Phan (2018) data
;des_sample_specific_start_time = 1670
;mean_fac = 1
;;------------------------------------------------------------------------------
;--------------------------
;For Wilder (2018) data
;2.2 Medium Guide Field Event
;The second example event we show (Figure 4) is the moderate guide field event on 
; 9 December 2015 that was reported by Wilder et al. (2017). In this particular 
; event, the guide field is ~0.5, and the LMN coordinates are 
;   L = [-0.104, -0.59, 0.80] 
;   M = [0.38, -0.77, -0.514] 
;   N = [0.920, 0.250, 0.302]
;****************************
;if running code for Wilder 22, un-comment the block below:
;****************************
;des_sample_specific_start_time = 6416 ;for roughly when the rx event starts
;des_sample_specific_start_time = 6387 ;for when the 2018 paper Fig4 starts
;des_time_slices = 23 ;for roughly when the rx event ends
;des_time_slices = 87
;mean_fac = 0
;LMN = 0 ;we use this since Wilder (2018) JGR gives the LMN coord.s in GSE coord.s in Sec. 2.2
;inst_fac = 1 ;need to run the code in instantaneous FAC for the Wilder interval
;****************************
; From Wilder (2017) PRL
; Time interval 5∶03:56.5 - 5∶03:57.2 UT, between vertical dotted lines in Fig. 4
; let a=6416, then
; des_yr[a], des_dmo[a], des_ddy[a], des_dhr[a], des_dmn[a], des_dsc[a], des_dmilli[a], des_dmicros[a], des_dnanos[a]
; 2015-12-09 05:03:56.482592
; 23 time slices later
; 2015-12-09 05:03:57.202601
; CDF_TT2000, dis_dist_time_data, dis_yr, dis_dmo, dis_ddy, dis_dhr, dis_dmn, dis_dsc, dis_dmilli, dis_dmicros, dis_dnanos, /BREAK
; CDF_TT2000, des_dist_time_data, des_yr, des_dmo, des_ddy, des_dhr, des_dmn, des_dsc, des_dmilli, des_dmicros, des_dnanos, /BREAK
;------------------------------------------------------------------------------
;--------------------------
;For Wilder (2018) data
;2.3 High Guide Field Event
;The final event to be shown is an event with a higher guide field of 1.3, on 
; 4 November 2015. The LMN directions for this event are 
;   L = [0.49, 0.86, 0.13]
;   M = [0.39, 0.35, 0.85]
;   N = [0.77, 0.37, 0.51]
;****************************
;if running code for Wilder 23, un-comment the block below:
;****************************
des_sample_specific_start_time = 2079 ;for when the 2018 paper Fig6 starts
des_time_slices = 60
mean_fac = 0
inst_fac = 1 ;need to run the code in instantaneous FAC for the Wilder interval
;****************************
; From Wilder (2017) PRL
; Time interval 5∶03:56.5 - 5∶03:57.2 UT, between vertical dotted lines in Fig. 4
; let a=6416, then
; des_yr[a], des_dmo[a], des_ddy[a], des_dhr[a], des_dmn[a], des_dsc[a], des_dmilli[a], des_dmicros[a], des_dnanos[a]
; 2015-12-09 05:03:56.482592
; 23 time slices later
; 2015-12-09 05:03:57.202601
; CDF_TT2000, dis_dist_time_data, dis_yr, dis_dmo, dis_ddy, dis_dhr, dis_dmn, dis_dsc, dis_dmilli, dis_dmicros, dis_dnanos, /BREAK
; CDF_TT2000, des_dist_time_data, des_yr, des_dmo, des_ddy, des_dhr, des_dmn, des_dsc, des_dmilli, des_dmicros, des_dnanos, /BREAK
;--------------------------
; Set this before loading so it doesnt get overwritten when sav file loaded
if n_elements(mean_fac) gt 0 then overwrite_mean_fac = mean_fac

; Load variables from save file.
;restore, "~/MMS_quick/Phan_data/inputvars_20161209090304e.sav" ; Phan data
;restore, "~/Documents/IDL_code/fpc_mrx/SWilder22_plots_and_codes/inputvars_20151209050044e.sav" ; Wilder 2.2 data
restore, "~/Documents/IDL_code/fpc_mrx/SWilder23_plots_and_codes/inputvars_20151104043404e.sav" ; Wilder 2.3 data

; Set this now
if n_elements(overwrite_mean_fac) gt 0 then mean_fac = overwrite_mean_fac

;**********************************************************************************
; Index setup
;**********************************************************************************
;--------------------------
; if using Wilder2018_22 data
; Find initial index, nti_des for electrons, nti_dis for ions
; 'tstart_dXs' is the start Epoch value
; If 'des_sample_specific_start_time' is set, nti_des and nti_dis will be set at the 
; closest mutual matching index
; 'des_sample_specific_start_time' will always be chosen w.r.t electrons
IF n_elements(des_sample_specific_start_time) GT 0 THEN BEGIN
  nti_des=des_sample_specific_start_time
  tstart_des=des_dist_time_data[nti_des]
  tstart_dis=fnn3(dis_dist_time_data,tstart_des,index=nti_dis)
  tstart_des=fnn3(des_dist_time_data,tstart_dis,index=nti_des)
  des_sample_specific_start_time=nti_des
  dis_sample_specific_start_time=nti_dis
  print, "des_sample_specific_start_time IS DEFINED"
  CDF_TT2000, tstart_des, des_yr, des_dmo, des_ddy, des_dhr, des_dmn, des_dsc, des_dmilli, des_dmicros, des_dnanos, /BREAK
  help, tstart_des, des_yr, des_dmo, des_ddy, des_dhr, des_dmn, des_dsc, des_dmilli, des_dmicros, des_dnanos
  CDF_TT2000, tstart_dis, dis_yr, dis_dmo, dis_ddy, dis_dhr, dis_dmn, dis_dsc, dis_dmilli, dis_dmicros, dis_dnanos, /BREAK
  help, tstart_dis, dis_yr, dis_dmo, dis_ddy, dis_dhr, dis_dmn, dis_dsc, dis_dmilli, dis_dmicros, dis_dnanos
  PRINT, "*****************************************"
ENDIF ELSE BEGIN
; Find first DIS time which is greater than both any FGM time and any E time
  nti_dis = (where(dis_dist_time_data GE max([min(fgm_time_data),min(e_time_data)]), ndistcheck))[0]
  tstart_dis=dis_dist_time_data[nti_dis]
  tstart_des=fnn3(des_dist_time_data,tstart_dis,index=nti_des)  
  des_sample_specific_start_time = nti_des
  dis_sample_specific_start_time=nti_dis
  print, "des_sample_specific_start_time NOT DEFINED"
  CDF_TT2000, tstart_des, des_yr, des_dmo, des_ddy, des_dhr, des_dmn, des_dsc, des_dmilli, des_dmicros, des_dnanos, /BREAK
  help, tstart_des, des_yr, des_dmo, des_ddy, des_dhr, des_dmn, des_dsc, des_dmilli, des_dmicros, des_dnanos
  CDF_TT2000, tstart_dis, dis_yr, dis_dmo, dis_ddy, dis_dhr, dis_dmn, dis_dsc, dis_dmilli, dis_dmicros, dis_dnanos, /BREAK
  help, tstart_dis, dis_yr, dis_dmo, dis_ddy, dis_dhr, dis_dmn, dis_dsc, dis_dmilli, dis_dmicros, dis_dnanos
  PRINT, "*****************************************"
ENDELSE

; Find final index, ntf_des for electrons, ntf_dis for ions
; 'tend_dXs' is the end Epoch value 
IF KEYWORD_SET(des_time_slices) THEN BEGIN
  ntot_des=des_time_slices
  ntf_des = nti_des + ntot_des - 1
  tend_des = des_dist_time_data[ntf_des]
  tend_dis = fnn3(dis_dist_time_data,tend_des,index=ntf_dis)
  tend_des=fnn3(des_dist_time_data,tend_dis,index=ntf_des)
  print, "des_time_slices IS DEFINED"
  CDF_TT2000, tend_des, des_yr, des_dmo, des_ddy, des_dhr, des_dmn, des_dsc, des_dmilli, des_dmicros, des_dnanos, /BREAK
  help, tend_des, des_yr, des_dmo, des_ddy, des_dhr, des_dmn, des_dsc, des_dmilli, des_dmicros, des_dnanos
  CDF_TT2000, tend_dis, dis_yr, dis_dmo, dis_ddy, dis_dhr, dis_dmn, dis_dsc, dis_dmilli, dis_dmicros, dis_dnanos, /BREAK
  help, tend_dis, dis_yr, dis_dmo, dis_ddy, dis_dhr, dis_dmn, dis_dsc, dis_dmilli, dis_dmicros, dis_dnanos
  PRINT, "*****************************************"
ENDIF ELSE BEGIN
  IF (MIN([MAX(fgm_time_data),MAX(e_time_data)]) GT MIN([MAX(des_dist_time_data),MAX(dis_dist_time_data)])) THEN tend = MIN([MAX(des_dist_time_data),MAX(dis_dist_time_data)]) $
    ELSE tend = MIN([MAX(fgm_time_data),MAX(e_time_data)])
  
  tend_dis=fnn3(dis_dist_time_data,tend,index=ntf_dis)
  tend_dis=dis_dist_time_data[ntf_dis-1]
  tend_dis=fnn3(dis_dist_time_data,tend_dis,index=ntf_dis)
  tend_des=fnn3(des_dist_time_data,tend_dis,index=ntf_des)  
  print, "des_time_slices NOT DEFINED"
  CDF_TT2000, tend_des, des_yr, des_dmo, des_ddy, des_dhr, des_dmn, des_dsc, des_dmilli, des_dmicros, des_dnanos, /BREAK
  help, tend_des, des_yr, des_dmo, des_ddy, des_dhr, des_dmn, des_dsc, des_dmilli, des_dmicros, des_dnanos
  CDF_TT2000, tend_dis, dis_yr, dis_dmo, dis_ddy, dis_dhr, dis_dmn, dis_dsc, dis_dmilli, dis_dmicros, dis_dnanos, /BREAK
  help, tend_dis, dis_yr, dis_dmo, dis_ddy, dis_dhr, dis_dmn, dis_dsc, dis_dmilli, dis_dmicros, dis_dnanos
  PRINT, "*****************************************"
ENDELSE
; Find ntot_des, the total number of data points in the electron interval
; Find ntot_dis, the total number of data points in the ion interval
ntot_des = ntf_des - nti_des + 1
des_time_slices = ntot_des 
ntot_dis = ntf_dis - nti_dis + 1
dis_time_slices = ntot_dis

IF KEYWORD_SET(fpc_ions) && ~KEYWORD_SET(fpc_electrons) THEN BEGIN
  time_slices = dis_time_slices
ENDIF

IF ~KEYWORD_SET(fpc_ions) && KEYWORD_SET(fpc_electrons) THEN BEGIN
  time_slices = des_time_slices
ENDIF
;--------------------------

IF KEYWORD_SET(fpc_ions) && ~KEYWORD_SET(fpc_electrons) THEN BEGIN
  time_slices = dis_time_slices
ENDIF

IF ~KEYWORD_SET(fpc_ions) && KEYWORD_SET(fpc_electrons) THEN BEGIN
  time_slices = des_time_slices
ENDIF

;**********************************************************************************
; Interpolate measurements to match cadences (FGM with EDP cadence)  and  (DIS with DES cadence)
;**********************************************************************************
; Interpolate FGM data to EDP cadence
; the b_fast array is used to Lorentz transform the e*_sc before
; downsampling the latter
; NOTE: bb_fast is the three components of the magnetic field vector
bb_fast=e_sc*0.
FOR i=0,2 DO bb_fast[i,*] = interp_no_extrap(bb_full[i,*], fgm_time_data, e_time_data)

; bbt_fast is the magnitude |B|, interpolated to the cadence of EDP
bbt_fast = interp_no_extrap(bb_full[3,*], fgm_time_data, e_time_data)
bbt_fast=REFORM(bbt_fast)

; Interpolate DIS Vbulk data to DES cadence
vibulk_GSE_fast=vebulk_gse*0.
FOR i=0,2 DO vibulk_GSE_fast[i,*] = interp_no_extrap(vibulk_gse[i,*], dis_moms_time_data, des_moms_time_data)

;the ion density, and temperatures interpolated to match cadence of the DES
nni_fast = interp_no_extrap(nni, dis_moms_time_data, des_moms_time_data)
nni_fast_INT=nni_fast[nti_des:ntf_des]
tipar_fast = interp_no_extrap(tipar, dis_moms_time_data, des_moms_time_data)
tiperp_fast = interp_no_extrap(tiperp, dis_moms_time_data, des_moms_time_data)

;**********************************************************************************
; Compute Bulk Species Flow Speeds (averaged over interval)
;**********************************************************************************
;Mean ion and electron bulk velocities over interval
vebulk_GSE_INT=vebulk_gse[*,nti_des:ntf_des]
vebulk_GSE_mean=MEAN(vebulk_GSE_INT,DIMENSION=2,/DOUBLE)

vibulk_GSE_INT=vibulk_gse[*,nti_dis:ntf_dis]
vibulk_GSE_mean=MEAN(vibulk_GSE_INT,DIMENSION=2,/DOUBLE)
vibulk_GSE_fast_INT=vibulk_GSE_fast[*,nti_des:ntf_des]

vbulk_GSE_INT = keyword_set(fpc_ions) ? vibulk_GSE_INT : vebulk_GSE_INT
vbulk_GSE_mean = keyword_set(fpc_ions) ? vibulk_GSE_mean : vebulk_GSE_mean

;**********************************************************************************
; Species index setup and generalizing variable names to match species
;**********************************************************************************
IF KEYWORD_SET(fpc_ions) && ~KEYWORD_SET(fpc_electrons)THEN BEGIN
  dist_time_data = dis_dist_time_data
  moms_time_data = dis_moms_time_data
  sample_specific_start_time = dis_sample_specific_start_time
  ;indices:
  nti = nti_dis
  ntot = ntot_dis
  ntf = ntf_dis
  ;Epochs:
  tstart = tstart_dis
  tend = tend_dis
  ;Distribution:
  ff = TEMPORARY(ffi)
  ff_err = TEMPORARY(ffi_err)
  ;Energy data:
  energy_data = dis_energy_data
ENDIF
IF ~KEYWORD_SET(fpc_ions) && KEYWORD_SET(fpc_electrons) THEN BEGIN
  dist_time_data = des_dist_time_data
  moms_time_data = des_moms_time_data
  sample_specific_start_time = des_sample_specific_start_time
  ;indices:
  nti = nti_des
  ntot = ntot_des
  ntf = ntf_des
  ;Epochs
  tstart = tstart_des
  tend = tend_des
  ;Distribution:
  ff = TEMPORARY(ffe)
  ff_err = TEMPORARY(ffe_err)
  ;Energy data:
  energy_data = des_energy_data
ENDIF

;**********************************************************************************
; Set start and end times of analysis (Epoch values at nti and ntf)
;**********************************************************************************
;Set time data for DxS over interval
tdxs_moms = moms_time_data[nti:ntf]

;Set time data for DIS over interval //for plotting purposes
tall_dis_dist=dis_dist_time_data[nti_dis:ntf_dis]
tall_dis_moms=dis_moms_time_data[nti_dis:ntf_dis]
IF ARRAY_EQUAL(tall_dis_dist,tall_dis_moms) THEN tall_dis=tall_dis_dist ELSE message, 'dis_dist and dis_moms times not equal'

;Set time data for DES over interval //for plotting purposes
tall_des_dist=des_dist_time_data[nti_des:ntf_des]
tall_des_moms=des_moms_time_data[nti_des:ntf_des]
IF ARRAY_EQUAL(tall_des_dist,tall_des_moms) THEN tall_des=tall_des_dist ELSE message, 'des_dist and des_moms times not equal'

; FGM Data
tstart_mag = fnn3(fgm_time_data, tstart, index=nti_mag)
tend_mag = fnn3(fgm_time_data, tend, index=ntf_mag)
tall_fgm = fgm_time_data[nti_mag:ntf_mag]

; EDP Data
tstart_edp = fnn3(e_time_data, tstart, index = nti_edp)
tend_edp = fnn3(e_time_data, tend, index = ntf_edp)
tall_edp = e_time_data[nti_edp:ntf_edp]

;**********************************************************************************
; Create Array of Values only over interval of analysis
;**********************************************************************************
; Magnetic Fields
; Extract magnetic field from bb_full
bb_INT = bb_full[0:2,nti_mag:ntf_mag]
bbt_INT = REFORM(bb_full[3,nti_mag:ntf_mag])

;**********************************************************************************
; Compute Mean Magnetic Field over Analysis Interval
;**********************************************************************************
bb_mean = MEAN(bb_INT,DIM=2)
bbt_mean = NORM(bb_mean)

;Compute unit vector in direction of bb_mean
bbhat_mean = bb_mean / bbt_mean

;******************************************************************************
; Create srot_LMN variable using LMN matrix given by Wilder (2018) Sec. 2.2
;******************************************************************************
IF (mean_fac) THEN BEGIN
  ;Define single FAC rotation matrix (vbulk_GSE_mean as second direction)
  srot=dblarr(3,3,1) ;NOTE (18FEB2022) ASA: is the last dimension of 1 needed here?
  srot[*,*,0]=create_fac(bbhat_mean,vbulk_GSE_mean)
ENDIF

;Define LMN rotation matrix for SWilder 22
;L = [-0.104, -0.59, 0.80]
;M = [0.38, -0.77, -0.514]
;N = [0.920, 0.250, 0.302]
;srot_LMN=TRANSPOSE([[L],[M],[N]])

;Define LMN rotation matrix for SWilder 23
L = [-0.49, 0.86, 0.13]
M = [-0.39, -0.35, 0.85]
N = [0.77, 0.37, 0.51]
srot_LMN=TRANSPOSE([[L],[M],[N]])
;******************************************************************************
; Rotate e bulk velocity and the ion bulk velocity (both for the specific 
;   interval INT) to LMN coords.    
;******************************************************************************
vebulk_LMN_INT=vebulk_GSE_INT*0D0      ;computed using srot matrix
FOR it=0, N_ELEMENTS(vebulk_GSE_INT[0,*])-1 DO vebulk_LMN_INT[*,it] = srot_LMN[*,*] # vebulk_GSE_INT[*,it]

;rotating ion bulk velocity to LMN coord.s
vibulk_LMN_INT=vibulk_GSE_INT*0D0
FOR it=0, N_ELEMENTS(vibulk_GSE_INT[0,*])-1 DO vibulk_LMN_INT[*,it] = srot_LMN # vibulk_GSE_INT[*,it]

;below is for Swilder 22
;To shift the distributions, as well as Lorentz transform the E-field, both to the reconnection geomtery
; frame (which is the frame where the ion flow is stationary at the B_L reversal), we manually select
; the point where B_L turns from negative to positive, at ~42 in the electron interval, corresponds to
; ~8 in the ion interval
;Rx_frame_vibulk_GSE=vibulk_GSE_INT[*,8]
;Rx_frame_vibulk_LMN=vibulk_LMN_INT[*,8]

;below is for Swilder 23
;To shift the distributions, as well as Lorentz transform the E-field, both to the reconnection geomtery
; frame (which is the frame where the ion flow is stationary at the B_L reversal), we manually select
; the point where B_L turns from negative to positive, at ~26 in the electron interval, corresponds to
; ~8 in the ion interval
Rx_frame_vibulk_GSE=vibulk_GSE_INT[*,5]
Rx_frame_vibulk_LMN=vibulk_LMN_INT[*,5]

;**********************************************************************************
; Compute Lorentz Transform of E-field to the Rx geometry plasma frame
;**********************************************************************************
e_l_GSE=e_sc*0.
e_l_LMN=e_sc*0.
tmp=size(e_sc)
nte=tmp[2]

;NOTE: *1D-3 converts  vxB electric field from V/m to mV/m
FOR i=0,nte-1 DO BEGIN
  e_l_GSE[*,i] = e_sc[*,i] + CROSSP(Rx_frame_vibulk_GSE, bb_fast[*,i]) * 1D-3
  e_l_LMN[*,i] = srot_LMN # e_l_GSE[*,i]
ENDFOR

;**********************************************************************************
; Create Array of Values only over interval of analysis
;**********************************************************************************
;Lorentz Transformed Electric Fields
e_l_GSE_INT = e_l_GSE[*, nti_edp:ntf_edp]
e_l_LMN_INT = e_l_LMN[*, nti_edp:ntf_edp]

; Spacecraft Potential
scpot_INT = scpot_full[nti_edp:ntf_edp]
scpot_mean = MEAN(scpot_INT, /DOUBLE)

;**********************************************************************************
; Correct Energy Bins for Spacecraft Potential
; Uncorrected: energy_data
; Corrected: energy_data_corr
;**********************************************************************************
IF (species EQ 'e') THEN energy_data_corr = energy_data[*,nti:ntf] - scpot_mean $
ELSE energy_data_corr = energy_data[*,nti:ntf] + scpot_mean

;**********************************************************************************
; Interleave mode check
; The interval of interest, which runs from nti:ntf in the full dataset, 
; runs from 0:ntot-1 in the interval dataset.
; The parity is defined relative to the indexing within the interval, 0:ntot-1 
;**********************************************************************************
vbinfactor = (species eq 'i') ? 938.272D6 : 0.51099895D6
interleave_check = is_2darr_rep(energy_data_corr)
if interleave_check eq 1 then begin 
  ; Same energy values across all data/not interleaved
  interleave = 0
  parity = 0
  vbin= dblarr(32,1)
  ; the 1D-3 at the end converts the velocity from m/s to km/s
  vbin[*,parity] = !CONST.C*Sqrt(2.*energy_data_corr[*,0]/vbinfactor)*1D-3
  ;NOTE(05DEC2021) ASA: the above line was previously: vbin[i,parity] = !CONST.C*Sqrt(2.*energy_data_corr[i,0]/vbinfactor)*1D-3  

endif else begin
  ; Different energy values -- check if equal across every other data point
  ;Check that energy bins are interleaved over entire analysis interval 
  evencheck = is_2darr_rep(energy_data_corr[*,0:*:2]) 
  odd_check = is_2darr_rep(energy_data_corr[*,1:*:2]) 
  interleave_check = (2b * (evencheck eq 0)) + (odd_check eq 0)
  
  ; Throw an error if they are not all equal
  case interleave_check of
  1: message, 'ERROR: Energy interleave is NOT on for all time slices (odd)'
  2: message, 'ERROR: Energy interleave is NOT on for all time slices (even)'
  3: message, 'ERROR: Energy interleave is NOT on for all time slices (odd and even)'
     else: break 
  endcase

  interleave = 1
  ; NOTE: parity is relative to interval indexing over 0:ntot-1
  ; the 1D-3 at the end converts the velocity from m/s to km/s
  vbin = !CONST.C*Sqrt(2.*energy_data_corr[*,0:1]/vbinfactor)*1D-3
  ;NOTE(02June2022) ASA: removed this for-loop: for j=0,1 do vbin[*,j] = !CONST.C*Sqrt(2.*energy_data_corr[*,j]/vbinfactor)*1D-3

endelse 

;******************************************************************************
; Compute mean phi bin values
;******************************************************************************
; phi angle locations for the chosen skymap.
phi_array = phi_data[*,nti:ntf]

; phi_mean used to put all of the 3V velocity coordinates at the same location across time
phi_mean=MEAN(phi_data[*,nti:ntf],DIM=2,/DOUBLE)

;******************************************************************************
; Convert from (phi,theta,energy) to (vx,vy,vz)_GSE
;******************************************************************************
; vv has units of km/s
vv=dblarr(3,nbin,interleave+1)
FOR ip=0,interleave DO BEGIN    ;For interleave, executes two loops
   FOR j=0,31 DO BEGIN          ;Energy/velocity loop
      FOR k=0,15 DO BEGIN       ;Polar angle theta loop
         FOR l=0,31 DO BEGIN    ;Azimuthal angle phi loop
            ;Compute linear index over (l,k,j)=(phi,theta,energy)
            i=j*16*32 +k*32 +l
            vv[0,i,ip] = vbin[j,ip]*SIN(theta[k])*COS(phi_mean[l])
            vv[1,i,ip] = vbin[j,ip]*SIN(theta[k])*SIN(phi_mean[l])
            vv[2,i,ip] = vbin[j,ip]*COS(theta[k])
         ENDFOR 
      ENDFOR
   ENDFOR
ENDFOR

;******************************************************************************
; Correct for Look Direction
;******************************************************************************
; Particle velocities are the opposite of the bin look directions
; (vx -> -vx, vy -> -vy, vz -> -vz)
vv=-vv

;******************************************************************************
; Subtract Bulk Flow from Dist Function Velocity Coordinates to shift to the 
;   plasma frame
;******************************************************************************
; vvp has units of km/s
vvp=dblarr(3,nbin,interleave+1)

;we shift to the "reconnection geometry frame, which is where 
;the ion flow is stationary at the B_L reversal. 
Rx_frame_vibulk_GSE_arr=REBIN(Rx_frame_vibulk_GSE,[n_elements(Rx_frame_vibulk_GSE),nbin,interleave+1])
vvp=vv-Rx_frame_vibulk_GSE_arr

;******************************************************************************
; For a rotation to LMN coord.s, rotate velocity coordinates vvp once for each interleave
;******************************************************************************
vvLMN=DBLARR(3,nbin,interleave+1)
vvLMN[*,*,0]=srot_LMN#vvp[*,*,0]
vvLMN[*,*,interleave]=srot_LMN#vvp[*,*,interleave]

;******************************************************************************
; For a single Mean FAC, rotate velocity coordinates vvp once
;******************************************************************************
UNDEFINE, vvfac
; vvfac has units of km/s; vvfac=dblarr(3,nbin) where for dim=0: 0=vprp1, 1=vprp2, 2=vpar
IF (mean_fac) THEN BEGIN
  vvfac=dblarr(3,nbin,interleave+1)
  vvfac[*,*,0]=srot#vvp[*,*,0]

  vvfac[*,*,interleave]=srot#vvp[*,*,interleave]
ENDIF

;**********************************************************************************
;**********************************************************************************
; MAIN LOOP
;**********************************************************************************
;**********************************************************************************
; Define arrays with ntot values over analysis interval (Means at DxS cadence)
b_GSE_dxs_mean_INT=dblarr(3,ntot)
bt_dxs_mean=dblarr(ntot)
b_fac_dxs_mean_INT=dblarr(3,ntot)   ; Computed using rot matrix:
;  can be w.r.t. local B-field (if mean_fac = 0) or ambient B-field (if mean_fac = 1)
b_LMN_dxs_mean_INT=dblarr(3,ntot)   ; array of B-field values in at cadence of DxS; LMN coord.s

e_l_GSE_dxs_mean_INT=dblarr(3,ntot)   ; array of E-field values at cadence of DxS; GSE coord.s
e_l_dxs_min=dblarr(3,ntot)    ; the min values at cadence of DxS; GSE coord.s
e_l_dxs_max=dblarr(3,ntot)    ; the max values at cadence of DxS; GSE coord.s
e_l_FAC_dxs_mean_INT=dblarr(3,ntot)   ; Computed using rot matrix:
                                  ;   can be w.r.t. local B-field (if mean_fac = 0)
                                  ;   or ambient B-field (if mean_fac = 1)
e_sc_dxs_mean=dblarr(3,ntot)  ; array of E-field values in sc-frame at cadence of DxS; GSE coord.s
e_sc_LMN_dxs_mean_INT=dblarr(3,ntot)  ; array of E-field values in sc-frame at cadence of DxS; LMN coord.s
e_l_LMN_dxs_mean_INT=dblarr(3,ntot)  ; array of E-field values in Lorentz transformed to Rx geometry frame,
                                 ;   at cadence of DxS; LMN coord.s

scpot_dxs_mean=dblarr(ntot)   ; array of scpot values at cadence of DxS

vbulk_FAC_INT=dblarr(3,ntot)      ;computed using rot matrix
vibulk_FAC_fast_INT=DBLARR(3,ntot) ;computed using srot matrix

FOR it=0,ntot - 1 DO BEGIN
  ;******************************************************************************
  ; Downsample fields (by averaging) to DxS cadence
  ;******************************************************************************
  ;Set start time and end time of current DES interval
  ttn=dist_time_data[nti+it]
  ttnp1=dist_time_data[nti+it+1]

  ;Find nearest FGM time to start time and end time of current DxS interval
  ;NOTE: the number of FGM points within in_fgm and inp1_fgm should equal to 3 or 4
  ; (for DES) or 19 or 20 (for DIS)
  ttn_fgm = fnn3(fgm_time_data, ttn, index=in_fgm)
  IF (it EQ 0) THEN ttn_fgm0 = ttn_fgm
  ttnp1_fgm = fnn3(fgm_time_data, ttnp1, index=inp1_fgm)
  tdxs_fgm = fgm_time_data[in_fgm:inp1_fgm]

  ;Find nearest EDP time to start time and end time of current DES interval
  ttn_edp = fnn3(e_time_data, ttn, index=in_edp)
  IF (it EQ 0) THEN ttn_edp0 = ttn_edp
  ttnp1_edp = fnn3(e_time_data, ttnp1, index=inp1_edp)

  ;Create Arrays with all measurements within current DxS interval
  b_GSE_dxs = bb_full[0:2,in_fgm:inp1_fgm]
  bt_dxs = bb_full[3,in_fgm:inp1_fgm] ; (10June2022) this is meaningless since we take the mean within DxS cadence (varibale "b_GSE_dxs_mean_INT"),
  ; and we need the magnitude of the mean (variable "bt_dxs_mean")
  e_l_GSE_dxs = e_l_GSE[0:2,in_edp:inp1_edp]
  e_para_dxs = e_sc_para[in_edp:inp1_edp]
  scpot_dxs = scpot_full[in_edp:inp1_edp]
  
  e_sc_dxs = e_sc[0:2,in_edp:inp1_edp]

  ;Compute Means of B- and E- Field over current DES interval
  b_GSE_dxs_mean_INT[0:2,it] = MEAN(b_GSE_dxs, dim = 2)
  bt_dxs_mean[it] = NORM(b_GSE_dxs_mean_INT[*,it]) ;magnitude of the means

  e_l_GSE_dxs_mean_INT[0,it] = MEAN(e_l_GSE_dxs, dim = 2)
  e_l_dxs_min[0,it] = MIN(e_l_GSE_dxs, dim = 2)
  e_l_dxs_max[0,it] = MAX(e_l_GSE_dxs, dim = 2)
  
  e_sc_dxs_mean[0,it] = MEAN(e_sc_dxs, dim = 2)

  scpot_dxs_mean[it] = MEAN(scpot_dxs)

;  ;******************************************************************************
;  ; For a instantaneous FAC, rotate velocity coordinates for each time
;  ;******************************************************************************
  IF (inst_FAC) && (it EQ 0) THEN BEGIN ;this is be executed only once in order to 
                                        ; initialize the vvfac and srot arrays
    
    ;NOTE: With Instantaneous FAC, we don't need two dimensions for interleave,
    ;      since the velocity grid at each timestep will be different
    ;      for dim=0: ; 0=vprp1, 1=vprp2, 2=vpar
    ; the following 2 lines are appropriate if we're doing instantaneous rotations
    vvfac=dblarr(3,nbin,ntot) 
    srot=dblarr(3,3,ntot)
  ENDIF
  IF (~mean_FAC) && (inst_FAC) THEN BEGIN
    ;Compute Rotation Matrix for Instantaneous (DES Averaged) FAC
    srot[*,*,it]=create_fac(b_GSE_dxs_mean_INT[*,it],vbulk_GSE_mean)
    ;Rotate velocity coordinates to FAC for this timeslice
    UNDEFINE, parity
    parity=(it mod 2)
    FOR i=0,nbin-1 DO vvfac[0:2,i,it] = srot[*,*,it] # vvp[0:2,i,parity]
  ENDIF
  
  ;******************************************************************************
  ;Use FAC rotation matrix to rotate electric field
  ; 0=vprp1, 1=vprp2, 2=vpar
  ;******************************************************************************
  srot_mult_ind = mean_fac ? 0 : it
;  srot_mult_ind = 0
  e_l_FAC_dxs_mean_INT[*,it] = srot[*,*,srot_mult_ind] # e_l_GSE_dxs_mean_INT[*,it]
  b_FAC_dxs_mean_INT[*,it] = srot[*,*,srot_mult_ind] # b_GSE_dxs_mean_INT[*,it]
  vbulk_FAC_INT[*,it] = srot[*,*,srot_mult_ind] # vbulk_GSE_INT[*,it]
  vibulk_FAC_fast_INT[*,it] = srot[*,*,srot_mult_ind] # vibulk_GSE_fast_INT[*,it]
  
  b_LMN_dxs_mean_INT[*,it] = srot_LMN[*,*] # b_GSE_dxs_mean_INT[*,it]
  e_sc_LMN_dxs_mean_INT[*,it] = srot_LMN # e_sc_dxs_mean[*,it]
  e_l_LMN_dxs_mean_INT[*,it] = srot_LMN # e_l_GSE_dxs_mean_INT[*,it]
ENDFOR

;**********************************************************************************
;**********************************************************************************
; END MAIN LOOP
;**********************************************************************************
;**********************************************************************************

;******************************************************************************
; Reorder Dist Function from (phi,theta,energy) to single array
;******************************************************************************
; VDF and error in single dimension velocity order
; Can just take the data ranges we need here
; the 1D12 converts from s^3/cm^6 to s^3/m^6
ff_v = reform( ff[*,*,*, nti:ntf], nbin, ntot) * 1d12
ff_v_err = reform( ff_err[*,*,*,nti:ntf], nbin, ntot) * 1d12

; Compute Counts in each velocity bin from Poisson Error
particle_count_v = ROUND( (ff_v / ff_v_err)^2)
; Set elements where ff_v and ff_v_err both equal 0 to 0
ff_ff_err0 = where( ((ff_v eq 0) and (ff_v_err eq 0)) eq 1, nff_ff_err0)
if nff_ff_err0 gt 0 then particle_count_v[ff_ff_err0] = 0

;**********************************************************************************
; Compute plasma parameters and mean plasma parameters (averaged over interval)
;**********************************************************************************
; Mean Ion and electron densities over interval
nni_INT=nni[nti_dis:ntf_dis]
nni_mean=MEAN(nni[nti_dis:ntf_dis],/DOUBLE)
nne_INT=nne[nti_des:ntf_des]
nne_mean=MEAN(nne[nti_des:ntf_des],/DOUBLE)

; Mean Ion and electron parallel and perpendicular temperatures over interval
;ions
tipara_INT=tipar[nti_dis:ntf_dis]
tipara_mean=MEAN(tipara_INT,/DOUBLE)
tiperp_INT=tiperp[nti_dis:ntf_dis]
tiperp_mean=MEAN(tiperp_INT,/DOUBLE)

;electrons
tepara_INT=tepar[nti_des:ntf_des]
tepara_mean=MEAN(tepara_INT,/DOUBLE)
teperp_INT=teperp[nti_des:ntf_des]
teperp_mean=MEAN(teperp_INT,/DOUBLE)

; Mean ion parallel and perpendicular thermal speeds over interval
vtipara_INT=SQRT(2*tipara_INT/938.272D6)*!CONST.c*1D-3
vtipara_mean=MEAN(vtipara_INT,/DOUBLE)
vtiperp_INT=SQRT(2*tiperp_INT/938.272D6)*!CONST.c*1D-3
vtiperp_mean=MEAN(vtiperp_INT,/DOUBLE)

; Mean electron parallel and perpendicular thermal speeds over interval
vtepara_INT=SQRT(2*tepara_INT/938.272D6)*!CONST.c*1D-3
vtepara_mean=MEAN(vtepara_INT,/DOUBLE)
vteperp_INT=SQRT(2*teperp_INT/938.272D6)*!CONST.c*1D-3
vteperp_mean=MEAN(vteperp_INT,/DOUBLE)

;Compute "isotropic" temperatures and thermal velocities
;ions
ti_mean=(tipara_mean + 2*tiperp_mean)/3
Ti_INT=(tipara_INT+2D0*tiperp_INT)/3
vti_INT=SQRT(2*Ti_INT/938.272D6)*!CONST.c*1D-3
vti_mean=SQRT(2*ti_mean/938.272D6)*!CONST.c*1D-3

;electrons
te_mean=(tepara_mean + 2*teperp_mean)/3
Te_INT=(tepara_INT+2D0*teperp_INT)/3
vte_INT=SQRT(2*Te_INT/.511099895D6)*!CONST.c*1D-3 ;returned in km/s
vte_mean=SQRT(2*te_mean/.511099895D6)*!CONST.c*1D-3 ;returned in km/s

; Calculate proton cyclotron frequency (Hz)
;   _INT is at DxS cadence. _mean is w.r.t. the magnitude of the ambient B-field
;   the 1D-9 converts from nT to T
fci_INT=(!const.e*bt_dxs_mean*1D-9)/(2*!const.pi*!const.mp)
fci_mean=(!const.e*bbt_mean*1D-9)/(2*!const.pi*!const.mp)

;calculating the Alfven speed
;   _dxs_arr is at DxS cadence.
;   _mean is w.r.t. the magnitude of the ambient B-field and the mean ion density
; 1D-9 converts bt_dxs_mean from nT to Tesla
; 1D6 converts the density (nni_INT) from cm^-3 to m^-3
; 1D-3 converts the speed from m/s to km/s
; thus, va_dxs_arr and va_mean are returned in units of km/s
va_dxs_arr=bt_dxs_mean*1D-9/(sqrt(!const.mu0*nni_INT*1D6*!const.mp))*1D-3
va_mean=bbt_mean*1D-9/(sqrt(!const.mu0*nni_mean*1D6*!const.mp))*1D-3


;Compute Plasma Beta Values
IF KEYWORD_SET(fpc_ions) && ~KEYWORD_SET(fpc_electrons) THEN BEGIN
  betai_para_INT=1.16D28*(2D0*!const.mu0*nni_INT*!const.k*tipara_INT)/bt_dxs_mean^2D0
  betai_para_mean=1.16D28*(2D0*!const.mu0*nni_mean*!const.k*tipara_mean)/bbt_mean^2D0
  betai_perp_INT=1.16D28*(2D0*!const.mu0*nni_INT*!const.k*tiperp_INT)/bt_dxs_mean^2D0
  betai_perp_mean=1.16D28*(2D0*!const.mu0*nni_mean*!const.k*tiperp_mean)/bbt_mean^2D0
  betai_INT = 1.16D28*(2D0*!const.mu0*nni_INT*!const.k*Ti_INT)/bt_dxs_mean^2D0
  betai_mean = 1.16D28*(2D0*!const.mu0*nni_mean*!const.k*Ti_mean)/bbt_mean^2D0

  betae_para_INT=1.16D28*(2D0*!const.mu0*nne_INT*!const.k*tepara_INT)/bbt_mean^2D0
  betae_para_mean=1.16D28*(2D0*!const.mu0*nne_mean*!const.k*tepara_mean)/bbt_mean^2D0
  betae_perp_INT=1.16D28*(2D0*!const.mu0*nne_INT*!const.k*teperp_INT)/bbt_mean^2D0
  betae_perp_mean=1.16D28*(2D0*!const.mu0*nne_mean*!const.k*teperp_mean)/bbt_mean^2D0
  betae_INT = 1.16D28*(2D0*!const.mu0*nne_INT*!const.k*Te_INT)/bbt_mean^2D0
  betae_mean = 1.16D28*(2D0*!const.mu0*nne_mean*!const.k*Te_mean)/bbt_mean^2D0
ENDIF

IF ~KEYWORD_SET(fpc_ions) && KEYWORD_SET(fpc_electrons) THEN BEGIN
  betai_para_INT=1.16D28*(2D0*!const.mu0*nni_INT*!const.k*tipara_INT)/bbt_mean^2D0
  betai_para_mean=1.16D28*(2D0*!const.mu0*nni_mean*!const.k*tipara_mean)/bbt_mean^2D0
  betai_perp_INT=1.16D28*(2D0*!const.mu0*nni_INT*!const.k*tiperp_INT)/bbt_mean^2D0
  betai_perp_mean=1.16D28*(2D0*!const.mu0*nni_mean*!const.k*tiperp_mean)/bbt_mean^2D0
  betai_INT = 1.16D28*(2D0*!const.mu0*nni_INT*!const.k*Ti_INT)/bbt_mean^2D0
  betai_mean = 1.16D28*(2D0*!const.mu0*nni_mean*!const.k*Ti_mean)/bbt_mean^2D0

  betae_para_INT=1.16D28*(2D0*!const.mu0*nne_INT*!const.k*tepara_INT)/bt_dxs_mean^2D0
  betae_para_mean=1.16D28*(2D0*!const.mu0*nne_mean*!const.k*tepara_mean)/bt_dxs_mean^2D0
  betae_perp_INT=1.16D28*(2D0*!const.mu0*nne_INT*!const.k*teperp_INT)/bt_dxs_mean^2D0
  betae_perp_mean=1.16D28*(2D0*!const.mu0*nne_mean*!const.k*teperp_mean)/bt_dxs_mean^2D0
  betae_INT = 1.16D28*(2D0*!const.mu0*nne_INT*!const.k*Te_INT)/bt_dxs_mean^2D0
  betae_mean = 1.16D28*(2D0*!const.mu0*nne_mean*!const.k*Te_mean)/bt_dxs_mean^2D0
ENDIF

;**********************************************************************************
; Save Variables from Analysis
;**********************************************************************************
STOP
; For Wilder (2018) JGR
SAVE, /variables, filename="fpc_mrx_Wilder2018_23_setup_LMN_and_FAC.sav"
print, STRING(7B)
STOP
END