;NOTE(22April2022) ASA: this is still a work in progress
;NOTE(29OCT2024) ASA: the "_mRx" appendage marks a change from the regular func. "ff_3v_binned.pro"
;  in that the ordering of the coordinates is no longer switched. i.e. we retain
;  the input order of vvfac, which is vL, vM, vN
;
; ff_3Vbinned_ntot returned in units of J/m^3 and dblarr(nvpar,vnvprp2,nvprp2,ntot)

FUNCTION ff_3V_binned_mRx, ff_v, vvfac, vt_mean, dvpar, dvprp, vpar_bins, vprp2_bins
  ;  print,junk
  ; the line below to accounts for if the detector energy bins are in interleave mode or not
  interleave = size(vvfac, /n_dim) eq 3 ? 1 : 0

  nbin = n_elements(ff_v[*,0])
  ntot = n_elements(ff_v[0,*])

  vparmin=vpar_bins[0]-dvpar/2D0
  vprp2min=vprp2_bins[0]-dvprp/2D0

  nvpar = N_ELEMENTS(vpar_bins)
  nvprp2=N_ELEMENTS(vprp2_bins)

  vmap = intarr(3, nbin, interleave+1)   ; 0=vL, 1=vM, 2=vN

  for ip=0,interleave do begin ; For interleave, executes two loops (even/odd)
    vmap[*,*,ip] = find_index_mRx(vvfac[*,*,ip],nvpar,nvprp2,dvpar,dvprp,vparmin,vprp2min)
  endfor

  ff_3Vbinned_ntot=dblarr(nvpar,nvprp2,nvprp2,ntot)

  FOR it=0,ntot-1 DO BEGIN
    FOR iv=0,nbin-1 DO BEGIN
      parity = it MOD (interleave+1)

      ivprp1 = vmap[0,iv,parity]
      ivprp2 = vmap[1,iv,parity]
      ivpar = vmap[2,iv,parity]
      ;Skip this (iv) if index maps to beyond vmin to vmax range
      IF (ivpar EQ -999) OR (ivprp1 EQ -999) OR (ivprp2 eq -999) THEN CONTINUE
      ff_3Vbinned_ntot[ivprp1,ivprp2,ivpar,it] = total( [ff_3Vbinned_ntot[ivprp1,ivprp2,ivpar,it], ff_v[iv,it]], /nan)
    ENDFOR
  ENDFOR

  return, ff_3Vbinned_ntot
END
