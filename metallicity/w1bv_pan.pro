
pro w1bv_pan

;;;;;; using the Pan sample

readcol, 'Pan_2019_tab1_fadded.dat', SNNames, z, dm15, dm15err, Phase, Hostnames,    Ebv_host, Morph, logM, logMperr, logMmerr,   logOH, logOHperr, logOHmerr,   AGNflag, f2535, f2535err, f3025, f3025err, $
format='A, F, F, F, F, A, F, A, F, F, F, F, F, F, A, F, F, F, F'


dm15[where(dm15 eq -99.0)]=!Values.F_NAN
dm15err[where(dm15err eq -99.0)]=!Values.F_NAN

logM[where(logM eq -99.0)]=!Values.F_NAN
logMperr[where(logMperr eq -99.0)]=!Values.F_NAN
logMmerr[where(logMmerr eq -99.0)]=!Values.F_NAN
logMperr[where(dm15err eq -99.0)]=!Values.F_NAN
logOH[where(logOH eq -99.0)]=!Values.F_NAN
logOHperr[where(logOHperr eq -99.0)]=!Values.F_NAN
logOHmerr[where(logOHmerr eq -99.0)]=!Values.F_NAN
f2535[where(f2535 eq -99.0)]=!Values.F_NAN
f2535err[where(f2535err eq -99.0)]=!Values.F_NAN
f3025[where(f3025 eq -99.0)]=!Values.F_NAN
f3025err[where(f3025err eq -99.0)]=!Values.F_NAN

fontsize=16

bpeak_obsmag_array=make_array(n_elements(snnames),6,value=!Values.F_NAN)
bpeak_obsmagerr_array=make_array(n_elements(snnames),6,value=!Values.F_NAN)

bpeak_mag_array=make_array(n_elements(snnames),6,value=!Values.F_NAN)
bpeak_magerr_array=make_array(n_elements(snnames),6,value=!Values.F_NAN)

bpeak_deredmag_array=make_array(n_elements(snnames),6,value=!Values.F_NAN)
bpeak_deredmagerr_array=make_array(n_elements(snnames),6,value=!Values.F_NAN)

colorebv_array=make_array(n_elements(snnames),value=!Values.F_NAN)
mwebv_array=make_array(n_elements(snnames),value=!Values.F_NAN)
hubbletype_array=make_array(n_elements(snnames),/integer, value=!Values.F_NAN)

restore, '../idl/host.sav', verbose=1  
; these are the first order coefficients for the SN1992A spectrum
rlambda_array=[6.20, 8.01,5.43, 4.92,4.16,3.16]


for n=0, n_elements(snnames) -1 do begin

	SNname= snnames[n]
	print, SNname


	panindex=where(host.snname_array eq SNname,count)

	if count ne -1 then bpeak_obsmag_array[n,*]=host.BPEAKAPPMAG_ARRAY[panindex,*]
	if count ne -1 then bpeak_obsmagerr_array[n,*]=host.BPEAKAPPMAGERR_ARRAY[panindex,*]
	if count ne -1 then mwebv_array[n]=host.av_schlafly_array[panindex]/3.1

	if SNname eq 'SN2009Y'  then bpeak_obsmag_array[n,1]=18.81
	if SNname eq 'SN2011fe' then bpeak_obsmag_array[n,3]=!Values.F_NAN

	bpeak_mag_array[n,*]=bpeak_obsmag_array[n,*]-rlambda_array*mwebv_array[n]-rlambda_array*ebv_host[n]
	bpeak_magerr_array[n,*]=sqrt(bpeak_obsmagerr_array[n,*]^2.0)

	colorebv_array[n]=(bpeak_obsmag_array[n,4]-bpeak_obsmag_array[n,5])-(-0.1)

	bpeak_deredmag_array[n,*]=bpeak_obsmag_array[n,*]-rlambda_array*colorebv_array[n]
	bpeak_deredmagerr_array[n,*]= sqrt(bpeak_obsmagerr_array[n,*]^2.0) + sqrt(bpeak_obsmagerr_array[n,4]^2.0+bpeak_obsmagerr_array[n,5]^2.0)*rlambda_array


	if strmatch(Morph[n], '*cE*') eq 1 then hubbletype_array[n] = -6 
	if strmatch(Morph[n], 'E*')   eq 1 then hubbletype_array[n] = -5
	if strmatch(Morph[n], 'E+*')  eq 1 then hubbletype_array[n] = -4

	if strmatch(Morph[n], 'S0-')  eq 1 then hubbletype_array[n] = -3
	if strmatch(Morph[n], 'S0+')  eq 1 then hubbletype_array[n] = -1
	if strmatch(Morph[n], 'S0')   eq 1 then hubbletype_array[n] = -2
	if strmatch(Morph[n], 'S0/a') eq 1 then hubbletype_array[n] = 0
	if strmatch(Morph[n], 'Sa')   eq 1 then hubbletype_array[n] = 1
	if strmatch(Morph[n], 'Sab')  eq 1 then hubbletype_array[n] = 2
	if strmatch(Morph[n], 'Sb')   eq 1 then hubbletype_array[n] = 3
	if strmatch(Morph[n], 'Sbc')  eq 1 then hubbletype_array[n] = 4
	if strmatch(Morph[n], 'Sc')   eq 1 then hubbletype_array[n] = 5
	if strmatch(Morph[n], 'Scd')  eq 1 then hubbletype_array[n] = 6
	if strmatch(Morph[n], 'Sd')   eq 1 then hubbletype_array[n] = 7
	if strmatch(Morph[n], 'Sdm')  eq 1 then hubbletype_array[n] = 8
	if strmatch(Morph[n], 'Sm')   eq 1 then hubbletype_array[n] = 9
	if strmatch(Morph[n], 'I')    eq 1 then hubbletype_array[n] = 10



endfor



print, ' '
for n=0, n_elements(snnames) -1 do print, logOH[n], logOHperr[n], logOHmerr[n], f2535[n], f2535err[n]
print, ' '
for n=0, n_elements(snnames) -1 do print, logOH[n], logOHperr[n], logOHmerr[n], f3025[n], f3025err[n]
print, ' '
for n=0, n_elements(snnames) -1 do print, logOH[n], logOHperr[n], logOHmerr[n], bpeak_mag_array[n,1]-bpeak_mag_array[n,4], sqrt(bpeak_magerr_array[n,1]^2.0+bpeak_magerr_array[n,4]^2.0)
print, ' '
for n=0, n_elements(snnames) -1 do print, logOH[n], logOHperr[n], logOHmerr[n], bpeak_mag_array[n,2]-bpeak_mag_array[n,4], sqrt(bpeak_magerr_array[n,2]^2.0+bpeak_magerr_array[n,4]^2.0)
print, ' '
for n=0, n_elements(snnames) -1 do print, logOH[n], logOHperr[n], logOHmerr[n], bpeak_mag_array[n,3]-bpeak_mag_array[n,4], sqrt(bpeak_magerr_array[n,3]^2.0+bpeak_magerr_array[n,4]^2.0)
print, ' '




for n=0, n_elements(snnames) -1 do print, SNNames[n],  ', ', z[n],  ', ', dm15[n],  ', ', dm15err[n],  ', ', Phase[n], ' ',  ', ', Hostnames[n],     ', ', Ebv_host[n],  ', ', Morph[n],  ', ', hubbletype_array[n],  ', ', logM[n],  ', ', logMperr[n],  ', ', logMmerr[n],    ', ', logOH[n],  ', ', logOHperr[n],  ', ', logOHmerr[n], ', ', AGNflag[n],  ', ', f2535[n],  ', ', f2535err[n],  ', ', f3025[n],  ', ', f3025err[n],  ', ', bpeak_obsmag_array[n,0],  ', ', bpeak_obsmag_array[n,1],  ', ', bpeak_obsmag_array[n,2],  ', ', bpeak_obsmag_array[n,3],  ', ', bpeak_obsmag_array[n,4],  ', ', bpeak_obsmag_array[n,5],  ', ', bpeak_obsmagerr_array[n,0],  ', ', bpeak_obsmagerr_array[n,1],  ', ', bpeak_obsmagerr_array[n,2],  ', ', bpeak_obsmagerr_array[n,3],  ', ', bpeak_obsmagerr_array[n,4],  ', ', bpeak_obsmagerr_array[n,5]






loadct, 33

colortable=intarr(n_elements(snnames))
for i=0,n_elements(snnames)-1 do colortable[i]=floor( ((bpeak_obsmag_array[[i],4]-bpeak_obsmag_array[[i],5])+0.1)/0.4*255)

cgWindow_SetDefs, PS_Decomposed=1



; from http://www.iluvatar.org/~dwijn/idlfigures
!p.font = 1
!p.thick = 2
!x.thick = 2
!y.thick = 2
!z.thick = 2
xsize = 8.8
wall = 0.04
margin=0.16
a = xsize/8.8 - (margin + wall)
b = a * 2d / (1 + sqrt(5))

ysize = (margin + b + wall )*8.8
ticklen = 0.01
xticklen = ticklen/b
yticklen = ticklen/a

nxticks=10

x1 = margin*8.8/xsize
x2 = x1 + a*8.8/xsize
xc = x2 + wall*8.8/xsize
y1 = margin*8.8/ysize
y2 = y1 + b*8.8/ysize
y3 = y2 + wall*8.8/ysize
y4 = y3 + b*8.8/ysize
yc = y4 + wall*8.8/ysize

xdata=[0,1,2,3]
ydata=[2,3,4,5]


restore, '../idl/SN2011fe_redbolmags161.sav'
; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]
febpeak=where(min(feredmags[4,0,*,0]) eq feredmags[4,0,*,0],count)

logOHsplit=8.6
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
plotfilename = 'bpeak_w1vbv_pan.eps'

SET_PLOT, 'PS'

device, filename= plotfilename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize


cgplot, charsize=1, feredmags[4,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[2,*,febpeak,3]-feredmags[5,*,febpeak,3], xrange=[-0.15,0.4], yrange=[0.8,2.2], ystyle=1, xstyle=1, ytitle='(w1-v)!BB!Lpeak', $ 
xtitle=' !A (b-v)!NBpeak   ', $
; double subscripts falling off page
; xtitle='!S!U (b-v) !N B !R!I peak', $
position=[x1,y1,x2,y2], linestyle=0, color=black
;  not getting this to work
;xyouts, 0.02, 0.5, '(b-v) !R!I B peak'

oplot, feredmags[4,*,febpeak,2]-feredmags[5,*,febpeak,2], feredmags[2,*,febpeak,2]-feredmags[5,*,febpeak,2], linestyle=1
oplot, feredmags[4,*,febpeak,1]-feredmags[5,*,febpeak,1], feredmags[2,*,febpeak,1]-feredmags[5,*,febpeak,1], linestyle=2
oplot, feredmags[4,*,febpeak,0]-feredmags[5,*,febpeak,0], feredmags[2,*,febpeak,0]-feredmags[5,*,febpeak,0], linestyle=3

;;;;;;;;;;;;;;;;;
;oploterror, logM, bpeak_mag_array[*,1]-bpeak_mag_array[*,4], sqrt(logMmerr^2.0), sqrt(bpeak_magerr_array[*,1]-bpeak_maerrg_array[*,4]^2.0), symsize=0.3, psym=15

cgoplot, bpeak_mag_array[where(logOH lt logOHsplit),4]-bpeak_mag_array[where(logOH lt logOHsplit),5], bpeak_mag_array[where(logOH lt logOHsplit),2]-bpeak_mag_array[where(logOH lt logOHsplit),5], psym=16, symsize=1, color='red'

cgoplot, bpeak_mag_array[where(logOH gt logOHsplit),4]-bpeak_mag_array[where(logOH gt logOHsplit),5], bpeak_mag_array[where(logOH gt logOHsplit),2]-bpeak_mag_array[where(logOH gt logOHsplit),5], psym=46, symsize=1.2, color='blue'


al_legend, ['logOH_array < 8.6','logOH_array > 8.6'], psym=[16,46], color=['red','blue'], symsize=[1,1.2], $
pos=[0.5,0.45], /norm, charsize=0.8, box=1

device, /close
SET_PLOT, 'X'
;spawn, 'open bpeak_w1vbv_pan.eps'

;;;;;;;;;;;;;;;;;;;;;;;;


for n=0,n_elements(snnames)-1 do print, snnames[n], bpeak_mag_array[n,4]-bpeak_mag_array[n,5], bpeak_mag_array[n,2]-bpeak_mag_array[n,5]


ohsort=sort(logoh)
print, 'Name, Log OH, B-V obs, B-V, W1-B
for n=0,n_elements(snnames)-1 do print, snnames[ohsort[n]], logoh[ohsort[n]], bpeak_obsmag_array[ohsort[n],4]-bpeak_obsmag_array[ohsort[n],5],bpeak_mag_array[ohsort[n],4]-bpeak_mag_array[ohsort[n],5], bpeak_mag_array[ohsort[n],2]-bpeak_mag_array[ohsort[n],4]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;logoh_range=[8.2,9.0]
;xrange=logoh_range
;xvalues=logoh
;err_xhigh=logOHperr
;err_xlow=-logOHmerr
;xtitle='Host log(OH)'

;figurename='uvcolors_logoh.eps'


m2b_range=[1.5,5]
w1b_range=[0,2.5]
ub_range=[-1,1]


;logoh_range=[8.2,9.0]

;xrange=logoh_range
;xvalues=logoh
;err_xhigh=logOHperr
;err_xlow=-logOHmerr
;xtitle='Host log(OH)'

;figurename='uvcolors_logoh.eps'


;readcol, 'Pan_2019_tab1_fadded.dat', SNNames, z, dm15, dm15err, Phase, Hostnames,    Ebv_host, Morph, logM, logMperr, logMmerr,   logOH, logOHperr, logOHmerr,   AGNflag, f2535, f2535err, f3025, f3025err, $
format='A, F, F, F, F, A, F, A, F, F, F, F, F, F, A, F, F, F, F'

;colorplots, figurename='uvcolors_z.eps', xtitle='Redshift', xrange=[0,0.02], xvalues=z, $
; bpeak_mag_array=bpeak_mag_array, bpeak_magerr_array=bpeak_magerr_array, m2b_range=m2b_range, w1b_range=w1b_range, ub_range=ub_range


; colorplots, figurename='uvcolors_phase.eps', xtitle='Phase', xrange=[-6,6], xvalues=phase, $
;bpeak_mag_array=bpeak_mag_array, bpeak_magerr_array=bpeak_magerr_array, m2b_range=m2b_range, w1b_range=w1b_range, ub_range=ub_range



colorplots, xtitle='Pan f2535', figurename='uvcolors_f2535.eps', xrange=[0,0.3], xvalues=f2535, err_xhigh=f2535err, err_xlow=f2535err, bpeak_mag_array=bpeak_mag_array, bpeak_magerr_array=bpeak_magerr_array, m2b_range=m2b_range, w1b_range=w1b_range, ub_range=ub_range


colorplots, xtitle='Pan f3025', figurename='uvcolors_f3025.eps', xrange=[0,1], xvalues=f3025, err_xhigh=f3025err, err_xlow=f3025err, bpeak_mag_array=bpeak_mag_array, bpeak_magerr_array=bpeak_magerr_array, m2b_range=m2b_range, w1b_range=w1b_range, ub_range=ub_range


;colorplots, xtitle='Host EBV', figurename='uvcolors_hostebv.eps', xrange=[0,0.3], xvalues=EBV_host, bpeak_mag_array=bpeak_mag_array, bpeak_magerr_array=bpeak_magerr_array, m2b_range=m2b_range, w1b_range=w1b_range, ub_range=ub_range


;colorplots, xtitle='Total EBV', figurename='uvcolors_totalebv.eps', xrange=[0,0.5], xvalues=EBV_host+mwebv_array, bpeak_mag_array=bpeak_mag_array, bpeak_magerr_array=bpeak_magerr_array, m2b_range=m2b_range, w1b_range=w1b_range, ub_range=ub_range


colorplots, xtitle='B - V', figurename='uvcolors_b-v.eps', xrange=[-0.5,0.6], xvalues=bpeak_mag_array[*,4]-bpeak_mag_array[*,5], bpeak_mag_array=bpeak_mag_array, bpeak_magerr_array=bpeak_magerr_array, m2b_range=m2b_range, w1b_range=w1b_range, ub_range=ub_range




;colorplots, xtitle='DM15(B)', figurename='uvcolors_dm15b.eps', xrange=[0.7,2], xvalues=dm15, err_xhigh=dm15err, err_xlow=dm15err, bpeak_mag_array=bpeak_mag_array, bpeak_magerr_array=bpeak_magerr_array, m2b_range=m2b_range, w1b_range=w1b_range, ub_range=ub_range


colorplots, xtitle='Host log(OH)', figurename='uvcolors_logoh.eps', xrange=[8.2,9], xvalues=logoh, err_xhigh=logOHperr, err_xlow=-logOHmerr, bpeak_mag_array=bpeak_mag_array, bpeak_magerr_array=bpeak_magerr_array, m2b_range=m2b_range, w1b_range=w1b_range, ub_range=ub_range


colorplots, xvalues=hubbletype_array, xrange=[-6,10], xtitle='Hubble Type', figurename='uvcolors_hubbletype.eps', bpeak_mag_array=bpeak_mag_array, bpeak_magerr_array=bpeak_magerr_array, m2b_range=m2b_range, w1b_range=w1b_range, ub_range=ub_range


colorplots, xvalues=logm, xrange=[6.0, 12], err_xhigh=logMperr, err_xlow=-logMmerr, xtitle='Host log(M)', figurename='uvcolors_logm.eps', bpeak_mag_array=bpeak_mag_array, bpeak_magerr_array=bpeak_magerr_array, m2b_range=m2b_range, w1b_range=w1b_range, ub_range=ub_range


;;;;;;;;;  try different ebv measurement


colorplots, xtitle='Host log(OH)', figurename='uvcolors_logoh_bv.eps', xrange=[8.2,9], xvalues=logoh, err_xhigh=logOHperr, err_xlow=-logOHmerr, bpeak_mag_array=bpeak_deredmag_array, bpeak_magerr_array=bpeak_deredmagerr_array, m2b_range=m2b_range, w1b_range=w1b_range, ub_range=ub_range


colorplots, xvalues=hubbletype_array, xrange=[-6,10], xtitle='Hubble Type', figurename='uvcolors_hubbletype_bv.eps', bpeak_mag_array=bpeak_deredmag_array, bpeak_magerr_array=bpeak_deredmagerr_array, m2b_range=m2b_range, w1b_range=w1b_range, ub_range=ub_range


colorplots, xvalues=logm, xrange=[6.0, 12], err_xhigh=logMperr, err_xlow=-logMmerr, xtitle='Host log(M)', figurename='uvcolors_logm_bv.eps', bpeak_mag_array=bpeak_deredmag_array, bpeak_magerr_array=bpeak_deredmagerr_array, m2b_range=m2b_range, w1b_range=w1b_range, ub_range=ub_range


;colorplots, xvalues=bpeak_deredmag_array[*,4]-bpeak_deredmag_array[*,5], xrange=[-0.2,0], err_xhigh=sqrt(bpeak_deredmagerr_array[*,4]^2.0+bpeak_deredmagerr_array[*,5]^2.0), err_xlow=sqrt(bpeak_deredmagerr_array[*,4]^2.0+bpeak_deredmagerr_array[*,5]^2.0), xtitle='B - V', figurename='uvcolors_b-v_bv.eps', bpeak_mag_array=bpeak_deredmag_array, bpeak_magerr_array=bpeak_deredmagerr_array, m2b_range=m2b_range, w1b_range=w1b_range, ub_range=ub_range



plot, ebv_host, bpeak_mag_array[*,4]-bpeak_mag_array[*,5]


print, 'final stop'
stop

end

pro colorplots, m2b_range=m2b_range, w1b_range=w1b_range, ub_range=ub_range, xrange=xrange, xvalues=xvalues, err_xhigh=err_xhigh, err_xlow=err_xlow, xtitle=xtitle, figurename=figurename, bpeak_mag_array=bpeak_mag_array, bpeak_magerr_array=bpeak_magerr_array

nplots=3
; from http://www.iluvatar.org/~dwijn/idlfigures
!p.font = 1
!p.thick = 2
!x.thick = 2
!y.thick = 2
!z.thick = 2
; the default size is given in centimeters
; 8.8 is made to match a journal column width
xsize = 8.8
wall = 0.03
margin=0.14
a = xsize/8.8 - (margin + wall)
b = a * 2d / (1 + sqrt(5))

ysize = (margin + nplots*(b + wall ) )*xsize
ticklen = 0.01
xticklen = ticklen/b
yticklen = ticklen/a

x1 = margin*8.8/xsize
x2 = x1 + a*8.8/xsize
xc = x2 + wall*8.8/xsize
y1 = margin*8.8/ysize
y2 = y1 + b*8.8/ysize




SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=12, bits_per_pixel=8, /color

x=0
cgplot, /noerase, xvalues, bpeak_mag_array[*,3]-bpeak_mag_array[*,4],  err_xhigh=err_xhigh, err_xlow=err_xlow, err_ylow=sqrt(bpeak_magerr_array[*,3]-bpeak_magerr_array[*,4]), err_yhigh=sqrt(bpeak_magerr_array[*,3]-bpeak_magerr_array[*,4]), psym=16, symsize=1, color='violet', $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=xtitle,   ytitle='u - b', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
yrange=ub_range, ystyle=1, xrange=xrange, xstyle=1

;y=5-x
x=1
cgplot, /noerase, xvalues, bpeak_mag_array[*,2]-bpeak_mag_array[*,4],  err_xhigh=err_xhigh, err_xlow=err_xlow, err_ylow=sqrt(bpeak_magerr_array[*,2]-bpeak_magerr_array[*,4]), err_yhigh=sqrt(bpeak_magerr_array[*,2]-bpeak_magerr_array[*,4]), psym=16, symsize=1, color='purple', $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=' ',   ytitle='uvw1 - b', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
yrange=w1b_range, ystyle=1, xrange=xrange, xstyle=1
;,  xtickname=replicate(' ',nxticks+1)

x=2

cgplot, /noerase, xvalues, bpeak_mag_array[*,1]-bpeak_mag_array[*,4],  err_xhigh=err_xhigh, err_xlow=err_xlow, err_ylow=sqrt(bpeak_magerr_array[*,1]-bpeak_magerr_array[*,4]), err_yhigh=sqrt(bpeak_magerr_array[*,1]-bpeak_magerr_array[*,4]), psym=16, symsize=1, color='maroon', $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=' ',   ytitle='uvm2 - b', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
yrange=m2b_range, ystyle=1, xrange=xrange, xstyle=1
;, xtickname=replicate(' ',nxticks+1)

device, /close
SET_PLOT, 'X'
spawn, 'open '+figurename+ ' &'


end


