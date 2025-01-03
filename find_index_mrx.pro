;NOTE(22April2022) ASA: this is still a work in progress
; Find the indices ivpar,ivprp1,ivprp2 for a given 3V point in FAC system
; This procedure inputs
;     1) vv(3) Vector of 3V velocity point
;     2) nvpar: Number of vpar bins
;     3) nvprp: Number of vprp bins
;     4) vpar_bins: Array with vpar bins
;     5) vprp_bins: Array with vprp bins
; and outputs
;     1) v_index =[ivpar,ivprp1,ivprp2]
; NOTE: IF the value is outside of the V limits, index is -999

FUNCTION find_index_mRx, vv,nvpar,nvprp,dvpar,dvprp,vparmin,vprpmin
  v_index= []

  ;Get non-gyrotropic values of 3V coordinate
  vprp1 = vv[0,*]
  vprp2 = vv[1,*]
  vpar1=vv[2,*]

  ivprp1=FLOOR((vprp1-vprpmin)/dvprp)
  ivprp2 = FLOOR((vprp2-vprpmin)/dvprp)
  ivpar=FLOOR((vpar1-vparmin)/dvpar)

  ivprp1[WHERE(ivprp1 LT 0)]=-999
  ivprp1[WHERE(ivprp1 GT nvprp-1)]=-999

  ivprp2[WHERE(ivprp2 LT 0)]=-999
  ivprp2[WHERE(ivprp2 GT nvprp-1)]=-999

  ivpar[WHERE(ivpar LT 0)]=-999
  ivpar[WHERE(ivpar GT nvpar-1)]=-999
  
  v_index=[ivprp1,ivprp2,ivpar]

  RETURN, v_index
END