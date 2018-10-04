pro w1bv

;, dmb15low, dmb15high, bvhigh

dmb15low=0.5
dmb15high=2.5
bvhigh=2

fontsize=16

restore, 'host.sav'  

index11fe=where(host.snname_array eq 'SN2011fe')
print, host.dm15_array[index11fe,4] 

restore, 'SN2011fe_redbolmags161.sav'
fedm=29.04

; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]



;normal=where(host.sntype2_array eq 'Ia' and host.dm15_array[*,4] gt dmb15low and host.dm15_array[*,4] lt dmb15high and (host.bpeakappmag_array[*,4]-host.bpeakappmag_array[*,5]) lt bvhigh)
 normal=where(host.sntype2_array eq 'Ia' and (host.bpeakappmag_array[*,4]-host.bpeakappmag_array[*,5]) lt bvhigh)
 
ia=where(host.sntype2_array eq 'Ia')

readcol,"snialist.txt",swiftsnia,quality,format='A,A',/silent
;normal=swiftsnia

bpeak_array=make_array(6,n_elements(normal),value=!Values.F_NAN)
bpeak_err_array=make_array(6,n_elements(normal),value=!Values.F_NAN)
firstepoch_array=make_array(n_elements(normal),value=!Values.F_NAN)
hstspectra_array=make_array(n_elements(normal),value=' ')
swiftspectra_array=make_array(n_elements(normal),value=' ')
nuvb_array=make_array(n_elements(normal),value=' ')
nuvr_array=make_array(n_elements(normal),value=' ')
milne_array=make_array(n_elements(normal),value=' ')

stritcolor_array=make_array(n_elements(normal),value=' ')

readcol,"swiftspectralist.txt",swiftspectra,format='A',/silent
readcol,"hstspectralist.txt",hstspectra,format='A',/silent

;readcol, 'Stritzinger_2018_redblue.txt', SN, Host, Redshift, EBV_MW, EBV_host, t_first, t_rise, DeltamB15, pm, dmerr, M_B, pm2, mberr, SpectralType, Color, References, format='(A, A, F, F, A, F, F, F, F, F, A, F, F, A, F, A, A, A, A)'

readcol, 'Stritzinger_2018_redblue.txt', strit_SNname, strit_Host, strit_Redshift, strit_EBV_MW, strit_EBV_host, pm3, strit_ebvhosterr, strit_t_first, strit_t_rise, strit_DeltamB15, pm, dmerr, strit_M_B, pm2, strit_mberr, strit_Spectralcode, strit_SpectralType, strit_Color, References, format='(A, A, A, A, A,A,A,A,A,A, A, A, A, A, A, A, A, A, A)', comment='#'

milnenuvb=['SN2011ia','SN2011fe','SN2011by','SN2008hv','SN2008Q']
milnenuvr=['SN2015F','SN2013gy','SN2013gs','SN2013ex','SN2013cs','SN2012hr', 'SN2011im','SN2008ec','SN2007co','SN2007af','SN2005df','SN2005cf']
;;; adding in those too red
milnenuvr=['SN2015F','SN2013gy','SN2013gs','SN2013ex','SN2013cs','SN2012hr', 'SN2011im','SN2008ec','SN2007co','SN2007af','SN2005df','SN2005cf', 'SN2013eu','SN2010kg','SN2010gp','SN2010ev']

milnemuvb=['SN2012ht', 'SN2006ej','SN2006dm','SN2010gn']



for n=0, n_elements(host.snname_array[normal]) -1 do begin

	SNname= host.snname_array[normal[n]]
	print, SNname

	Bpeaktime=host.bpeakmjd_array[normal[n]]


	bpeak_array[*,n]=host.bpeakappmag_array[normal[n],*]
	bpeak_err_array[*,n]=host.bpeakappmagerr_array[normal[n],*]
	;for n=0,n_elements(normal)-1 do print, host.snname_array[normal[n]], w1vbluest_array[*,n], host.dm15_array[normal[n],4], firstepoch_array[n]

	swiftindex=where(swiftspectra eq SNname)
	hstindex=where(hstspectra eq SNname)
	if swiftindex[0] ne -1 then swiftspectra_array[n]='swift'
	if hstindex[0] ne -1 then hstspectra_array[n]='hst'
	stritindex=where(strit_snname eq SNname)	
	if stritindex[0] ne -1 then stritcolor_array[n]=strit_color[stritindex]

	muvbindex  =where(milnemuvb eq SNname)
	nuvbindex  =where(milnenuvb eq SNname)
	nuvrindex  =where(milnenuvr eq SNname)
	if nuvbindex[0] ne -1 then milne_array[n]='nuvb'
	if muvbindex[0] ne -1 then milne_array[n]='muvb'
	if nuvrindex[0] ne -1 then milne_array[n]='nuvr'
endfor


hst=where(hstspectra_array eq 'hst')
swift=where(swiftspectra_array eq 'swift')
stritred=where(stritcolor_array eq 'red')
stritblue=where(stritcolor_array eq 'blue')


muvb=where(milne_array eq 'muvb')
nuvb=where(milne_array eq 'nuvb')
nuvr=where(milne_array eq 'nuvr')

; plot, host.dm15_array[normal[n],4],w1vbluest_array[2,n]-w1vbluest_array[5,n]

;for n=0, n_elements(host.snname_array[normal]) -1 do print, host.snname_array[normal[n]], w1vbluest_array[2,n]-w1vbluest_array[5,n]
; save, 

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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


plotfilename = 'bpeak_w1vbv_spectra.eps'

SET_PLOT, 'PS'

device, filename= plotfilename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

restore, 'SN2011fe_redbolmags161.sav'
; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]
febpeak=where(min(feredmags[4,0,*,0]) eq feredmags[4,0,*,0],count)

; print, ((feredmags[2,10,febpeak,0]-feredmags[5,10,febpeak,0])-(feredmags[2,0,febpeak,0]-feredmags[5,0,febpeak,0]))/((feredmags[4,10,febpeak,0]-feredmags[5,10,febpeak,0])-(feredmags[4,0,febpeak,0]-feredmags[5,0,febpeak,0]))

;print, ((feredmags[2,10,febpeak,1]-feredmags[5,10,febpeak,1])-(feredmags[2,0,febpeak,1]-feredmags[5,0,febpeak,1]))/((feredmags[4,10,febpeak,1]-feredmags[5,10,febpeak,1])-(feredmags[4,0,febpeak,1]-feredmags[5,0,febpeak,1]))
;print, ((feredmags[2,10,febpeak,2]-feredmags[5,10,febpeak,2])-(feredmags[2,0,febpeak,2]-feredmags[5,0,febpeak,2]))/((feredmags[4,10,febpeak,2]-feredmags[5,10,febpeak,2])-(feredmags[4,0,febpeak,2]-feredmags[5,0,febpeak,2]))
;print, ((feredmags[2,10,febpeak,3]-feredmags[5,10,febpeak,3])-(feredmags[2,0,febpeak,3]-feredmags[5,0,febpeak,3]))/((feredmags[4,10,febpeak,3]-feredmags[5,10,febpeak,3])-(feredmags[4,0,febpeak,3]-feredmags[5,0,febpeak,3]))

cgplot, charsize=1, feredmags[4,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[2,*,febpeak,3]-feredmags[5,*,febpeak,3], xrange=[-0.15,0.3], yrange=[0.8,2.2], ystyle=1, xstyle=1, ytitle='(w1-v)!BB!Lpeak', $ 
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
oploterror, bpeak_array[4,*]-bpeak_array[5,*], bpeak_array[2,*]-bpeak_array[5,*], sqrt(bpeak_err_array[4,*]^2.0+bpeak_err_array[5,*]^2.0), sqrt(bpeak_err_array[2,*]^2.0+bpeak_err_array[5,*]^2.0), symsize=0.3, psym=15

if hst[0] ne -1 then cgoplot, bpeak_array[4,hst]-bpeak_array[5,hst], bpeak_array[2,hst]-bpeak_array[5,hst], psym=16, symsize=1, color='red'

if swift[0] ne -1 then cgoplot, bpeak_array[4,swift]-bpeak_array[5,swift], bpeak_array[2,swift]-bpeak_array[5,swift], psym=46, symsize=1.2, color='blue'

xyouts, bpeak_array[4,hst]-bpeak_array[5,hst] - 0.06, bpeak_array[2,hst]-bpeak_array[5,hst]+0.02, host.snname_array[normal[hst]], charsize=0.5
xyouts, bpeak_array[4,swift]-bpeak_array[5,swift] - 0.06, bpeak_array[2,swift]-bpeak_array[5,swift]+0.02, host.snname_array[normal[swift]], charsize=0.5

al_legend, ['HST','Swift'], psym=[16,46], color=['red','blue'], symsize=[1,1.2], $
pos=[0.7,0.45], /norm, charsize=0.8, box=1

device, /close
SET_PLOT, 'X'
spawn, 'open bpeak_w1vbv_spectra.eps'

;;;;;;;;;;;;;;;;;;;;;;;;

plotfilename = 'bpeak_m2w1vw1v_spectra.eps'

SET_PLOT, 'PS'

device, filename= plotfilename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

restore, 'SN2011fe_redbolmags161.sav'
; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]

cgplot, charsize=1, feredmags[1,*,febpeak,3]-feredmags[2,*,febpeak,3], feredmags[2,*,febpeak,3]-feredmags[5,*,febpeak,3], xrange=[1,3.5], yrange=[0.8,2.2], $ 
xtitle=' !A (m2-w1)!NBpeak   ', $
ystyle=1, xstyle=1, ytitle='(w1-v)!BB!Lpeak', position=[x1,y1,x2,y2], linestyle=0, color='black'
oplot, feredmags[1,*,febpeak,2]-feredmags[2,*,febpeak,2], feredmags[2,*,febpeak,2]-feredmags[5,*,febpeak,2], linestyle=1
oplot, feredmags[1,*,febpeak,1]-feredmags[2,*,febpeak,1], feredmags[2,*,febpeak,1]-feredmags[5,*,febpeak,1], linestyle=2
oplot, feredmags[1,*,febpeak,0]-feredmags[2,*,febpeak,0], feredmags[2,*,febpeak,0]-feredmags[5,*,febpeak,0], linestyle=3

;;;;;;;;;;;;;;;;;
oploterror, bpeak_array[1,*]-bpeak_array[2,*], bpeak_array[2,*]-bpeak_array[5,*], sqrt(bpeak_err_array[1,*]^2.0+bpeak_err_array[2,*]^2.0), sqrt(bpeak_err_array[2,*]^2.0+bpeak_err_array[5,*]^2.0), symsize=0.3, psym=15

if hst[0] ne -1 then cgoplot, bpeak_array[1,hst]-bpeak_array[2,hst], bpeak_array[2,hst]-bpeak_array[5,hst], psym=16, symsize=1, color='red'

if swift[0] ne -1 then cgoplot, bpeak_array[1,swift]-bpeak_array[2,swift], bpeak_array[2,swift]-bpeak_array[5,swift], psym=46, symsize=1.2, color='blue'

xyouts, bpeak_array[1,hst]-bpeak_array[2,hst] - 0.06, bpeak_array[2,hst]-bpeak_array[5,hst]+0.02, host.snname_array[normal[hst]], charsize=0.5
xyouts, bpeak_array[1,swift]-bpeak_array[2,swift] - 0.06, bpeak_array[2,swift]-bpeak_array[5,swift]+0.02, host.snname_array[normal[swift]], charsize=0.5

al_legend, ['HST','Swift'], psym=[16,46], color=['red','blue'], symsize=[1,1.2], $
pos=[0.7,0.45], /norm, charsize=0.8, box=1

device, /close
SET_PLOT, 'X'
spawn, 'open bpeak_m2w1vw1v_spectra.eps'

;;;;;;;;;;;;;;;;;;;;;;;;

plotfilename = 'bpeak_w1vbv_strit.eps'

SET_PLOT, 'PS'

device, filename= plotfilename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

restore, 'SN2011fe_redbolmags161.sav'
; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]
febpeak=where(min(feredmags[4,0,*,0]) eq feredmags[4,0,*,0],count)


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
oploterror, bpeak_array[4,*]-bpeak_array[5,*], bpeak_array[2,*]-bpeak_array[5,*], sqrt(bpeak_err_array[4,*]^2.0+bpeak_err_array[5,*]^2.0), sqrt(bpeak_err_array[2,*]^2.0+bpeak_err_array[5,*]^2.0), symsize=0.3, psym=15

cgoplot, bpeak_array[4,stritred]-bpeak_array[5,stritred], bpeak_array[2,stritred]-bpeak_array[5,stritred], psym=16, symsize=1, color='red'

cgoplot, bpeak_array[4,stritblue]-bpeak_array[5,stritblue], bpeak_array[2,stritblue]-bpeak_array[5,stritblue], psym=46, symsize=1.2, color='blue'

;xyouts, bpeak_array[4,hst]-bpeak_array[5,hst] - 0.06, bpeak_array[2,hst]-bpeak_array[5,hst]+0.02, host.snname_array[normal[hst]], charsize=0.5
;xyouts, bpeak_array[4,swift]-bpeak_array[5,swift] - 0.06, bpeak_array[2,swift]-bpeak_array[5,swift]+0.02, host.snname_array[normal[swift]], charsize=0.5

al_legend, ['Strit+18 red','Strit+18 blue'], psym=[16,46], color=['red','blue'], symsize=[1,1.2], $
pos=[0.5,0.45], /norm, charsize=0.8, box=1

device, /close
SET_PLOT, 'X'
spawn, 'open bpeak_w1vbv_strit.eps'

;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

SET_PLOT, 'PS'

device, filename='w1v_bv_bpeak_nuvbr.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

xrange=[-0.2,0.2]
nxticks=4

cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='b-v',   ytitle='(w1-v)!BB!Lpeak', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[0.5,2.5], ystyle=1, xrange=[-0.2,0.3], xstyle=1, $
xticks=5, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, color='black'


;;;;;;;;;;   at b band maximum light


cgoplot, color='black', feredmags[4,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[2,*,febpeak,3]-feredmags[5,*,febpeak,3], linestyle=0
cgoplot, color='black', feredmags[4,*,febpeak,2]-feredmags[5,*,febpeak,2], feredmags[2,*,febpeak,2]-feredmags[5,*,febpeak,2], linestyle=1
cgoplot, color='black', feredmags[4,*,febpeak,1]-feredmags[5,*,febpeak,1], feredmags[2,*,febpeak,1]-feredmags[5,*,febpeak,1], linestyle=2
cgoplot, color='black', feredmags[4,*,febpeak,0]-feredmags[5,*,febpeak,0], feredmags[2,*,febpeak,0]-feredmags[5,*,febpeak,0], linestyle=3




w1vsne=where(finite(bpeak_array[2,*]) eq 1 and finite(bpeak_array[4,*]) eq 1 and finite(bpeak_array[5,*]) eq 1 )

nuvrw1vsne=where(finite(bpeak_array[2,*]) eq 1 and finite(bpeak_array[4,*]) eq 1 and finite(bpeak_array[5,*]) eq 1  and milne_array eq 'nuvr')

nuvbw1vsne=where(finite(bpeak_array[2,*]) eq 1 and finite(bpeak_array[4,*]) eq 1 and finite(bpeak_array[5,*]) eq 1  and milne_array eq 'nuvb')

muvbw1vsne=where(finite(bpeak_array[2,*]) eq 1 and finite(bpeak_array[4,*]) eq 1 and finite(bpeak_array[5,*]) eq 1  and milne_array eq 'muvb')

oploterror, bpeak_array[4,w1vsne]-bpeak_array[5,w1vsne], bpeak_array[2,w1vsne]-bpeak_array[5,w1vsne], sqrt(bpeak_err_array[4,w1vsne]^2.0+bpeak_err_array[5,w1vsne]^2.0), sqrt(bpeak_err_array[2,w1vsne]^2.0+bpeak_err_array[5,w1vsne]^2.0), psym=15, symsize=0.5, color='black'

oploterror, bpeak_array[4,nuvrw1vsne]-bpeak_array[5,nuvrw1vsne], bpeak_array[2,nuvrw1vsne]-bpeak_array[5,nuvrw1vsne], sqrt(bpeak_err_array[4,nuvrw1vsne]^2.0+bpeak_err_array[5,nuvrw1vsne]^2.0), sqrt(bpeak_err_array[2,nuvrw1vsne]^2.0+bpeak_err_array[5,nuvrw1vsne]^2.0), psym=16, color='red'

oploterror, bpeak_array[4,nuvbw1vsne]-bpeak_array[5,nuvbw1vsne], bpeak_array[2,nuvbw1vsne]-bpeak_array[5,nuvbw1vsne], sqrt(bpeak_err_array[4,nuvbw1vsne]^2.0+bpeak_err_array[5,nuvbw1vsne]^2.0), sqrt(bpeak_err_array[2,nuvbw1vsne]^2.0+bpeak_err_array[5,nuvbw1vsne]^2.0), psym=15, color='royal blue'


oploterror, bpeak_array[4,muvbw1vsne]-bpeak_array[5,muvbw1vsne], bpeak_array[2,muvbw1vsne]-bpeak_array[5,muvbw1vsne], sqrt(bpeak_err_array[4,muvbw1vsne]^2.0+bpeak_err_array[5,muvbw1vsne]^2.0), sqrt(bpeak_err_array[2,muvbw1vsne]^2.0+bpeak_err_array[5,muvbw1vsne]^2.0), psym=46, color='dark green', symsize=1.2

al_legend, ['NUV-red','MUV-blue','NUV-blue'], psym=[16,46,15], color=['red', 'dark green','royal blue'], position=[-0.2,2.5], box=0, charsize=0.9

;xyouts, bpeak_array[4,w1vsne]-bpeak_array[5,w1vsne] + 0.01, bpeak_array[2,w1vsne]-bpeak_array[5,w1vsne] - 0.02, host.snname_array[normal[w1vsne]], charsize=0.5

device, /close
SET_PLOT, 'X'


spawn, 'open w1v_bv_bpeak_nuvbr.eps'

;;;;;;;;;

SET_PLOT, 'PS'

device, filename='w1m2_bv_bpeak_nuvbr.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

xrange=[-0.2,0.2]
nxticks=4

cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='b-v',   ytitle='(m2-w1)!BB!Lpeak', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[0.5,4.5], ystyle=1, xrange=[-0.2,0.3], xstyle=1, $
xticks=5, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, color='black'


;;;;;;;;;;   at b band maximum light


cgoplot, color='black', feredmags[4,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[1,*,febpeak,3]-feredmags[2,*,febpeak,3], linestyle=0
cgoplot, color='black', feredmags[4,*,febpeak,2]-feredmags[5,*,febpeak,2], feredmags[1,*,febpeak,2]-feredmags[2,*,febpeak,2], linestyle=1
cgoplot, color='black', feredmags[4,*,febpeak,1]-feredmags[5,*,febpeak,1], feredmags[1,*,febpeak,1]-feredmags[2,*,febpeak,1], linestyle=2
cgoplot, color='black', feredmags[4,*,febpeak,0]-feredmags[5,*,febpeak,0], feredmags[1,*,febpeak,0]-feredmags[2,*,febpeak,0], linestyle=3




w1bvsne=where(finite(bpeak_array[1,*]) eq 1 and finite(bpeak_array[2,*]) eq 1 and finite(bpeak_array[4,*]) eq 1 and finite(bpeak_array[5,*]) eq 1 )

nuvrw1vsne=where(finite(bpeak_array[1,*]) eq 1 and finite(bpeak_array[2,*]) eq 1 and finite(bpeak_array[4,*]) eq 1 and finite(bpeak_array[5,*]) eq 1  and milne_array eq 'nuvr')

nuvbw1vsne=where(finite(bpeak_array[1,*]) eq 1 and finite(bpeak_array[2,*]) eq 1 and finite(bpeak_array[4,*]) eq 1 and finite(bpeak_array[5,*]) eq 1  and milne_array eq 'nuvb')

muvbw1vsne=where(finite(bpeak_array[1,*]) eq 1 and finite(bpeak_array[2,*]) eq 1 and finite(bpeak_array[4,*]) eq 1 and finite(bpeak_array[5,*]) eq 1  and milne_array eq 'muvb')

oploterror, bpeak_array[4,w1vsne]-bpeak_array[5,w1vsne], bpeak_array[1,w1vsne]-bpeak_array[2,w1vsne], sqrt(bpeak_err_array[4,w1vsne]^2.0+bpeak_err_array[5,w1vsne]^2.0), sqrt(bpeak_err_array[2,w1vsne]^2.0+bpeak_err_array[1,w1vsne]^2.0), psym=15, symsize=0.5, color='black'

oploterror, bpeak_array[4,nuvrw1vsne]-bpeak_array[5,nuvrw1vsne], bpeak_array[1,nuvrw1vsne]-bpeak_array[2,nuvrw1vsne], sqrt(bpeak_err_array[4,nuvrw1vsne]^2.0+bpeak_err_array[5,nuvrw1vsne]^2.0), sqrt(bpeak_err_array[2,nuvrw1vsne]^2.0+bpeak_err_array[1,nuvrw1vsne]^2.0), psym=16, color='red'

oploterror, bpeak_array[4,nuvbw1vsne]-bpeak_array[5,nuvbw1vsne], bpeak_array[1,nuvbw1vsne]-bpeak_array[2,nuvbw1vsne], sqrt(bpeak_err_array[4,nuvbw1vsne]^2.0+bpeak_err_array[5,nuvbw1vsne]^2.0), sqrt(bpeak_err_array[2,nuvbw1vsne]^2.0+bpeak_err_array[5,nuvbw1vsne]^2.0), psym=15, color='royal blue'


oploterror, bpeak_array[4,muvbw1vsne]-bpeak_array[5,muvbw1vsne], bpeak_array[1,muvbw1vsne]-bpeak_array[2,muvbw1vsne], sqrt(bpeak_err_array[4,muvbw1vsne]^2.0+bpeak_err_array[5,muvbw1vsne]^2.0), sqrt(bpeak_err_array[2,muvbw1vsne]^2.0+bpeak_err_array[1,muvbw1vsne]^2.0), psym=46, color='dark green', symsize=1.2

al_legend, ['NUV-red','MUV-blue','NUV-blue'], psym=[16,46,15], color=['red', 'dark green','royal blue'], position=[-0.2,4.5], box=0, charsize=0.9

;xyouts, bpeak_array[4,w1vsne]-bpeak_array[5,w1vsne] + 0.01, bpeak_array[2,w1vsne]-bpeak_array[5,w1vsne] - 0.02, host.snname_array[normal[w1vsne]], charsize=0.5

device, /close
SET_PLOT, 'X'


spawn, 'open w1m2_bv_bpeak_nuvbr.eps'



;;;;;;;;;


xsize = 8.8
wall = 0.04
margin=0.16
a = xsize/8.8 - (margin + wall)
b = a * 2d / (1 + sqrt(5))

ysize = (margin + b + wall + b + wall + b + wall)*8.8
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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

SET_PLOT, 'PS'

device, filename='w1v_bpeak_colors_abs_color.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

xrange=[-0.2,0.2]
nxticks=4

loadct, 33

colortable=intarr(n_elements(normal))
for i=0,n_elements(normal)-1 do colortable[i]=floor( ((bpeak_array[4,i]-bpeak_array[5,i])+0.1)/0.4*255)
colortable[where(colortable lt 0.0)]=0
colortable[where(colortable gt 255)]=255
                                               
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2+y2-y1,x2,y2+y2+y2-y1-y1], $
xrange=[-0.1,0.3], yrange=[3.0,0.5], ytitle='(w1-v)!BB!Lpeak',  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
ystyle=1, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color='black'


cgcolorbar, ncolors=256 , POSITION=[x1,y2+y2+y2-y1-y1, x2, 0.95], range=[-0.1,0.3], /top, charsize=1

print, 'test 3'



cgoplot, color='black', feredmags[4,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[2,*,febpeak,3]-feredmags[5,*,febpeak,3], linestyle=0
cgoplot, color='black', feredmags[4,*,febpeak,2]-feredmags[5,*,febpeak,2], feredmags[2,*,febpeak,2]-feredmags[5,*,febpeak,2], linestyle=1
cgoplot, color='black', feredmags[4,*,febpeak,1]-feredmags[5,*,febpeak,1], feredmags[2,*,febpeak,1]-feredmags[5,*,febpeak,1], linestyle=2
cgoplot, color='black', feredmags[4,*,febpeak,0]-feredmags[5,*,febpeak,0], feredmags[2,*,febpeak,0]-feredmags[5,*,febpeak,0], linestyle=3

; could maybe replace with SN2017erp
;
;restore, 'SN2005cf_redbolmags161.sav'
; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]

;oplot, cfredmags[4,*,0,3]-cfredmags[5,*,0,3], cfredmags[2,*,0,3]-cfredmags[5,*,0,3], linestyle=0
;oplot, cfredmags[4,*,0,2]-cfredmags[5,*,0,2], cfredmags[2,*,0,2]-cfredmags[5,*,0,2], linestyle=1
;oplot, cfredmags[4,*,0,1]-cfredmags[5,*,0,1], cfredmags[2,*,0,1]-cfredmags[5,*,0,1], linestyle=2
;oplot, cfredmags[4,*,0,0]-cfredmags[5,*,0,0], cfredmags[2,*,0,0]-cfredmags[5,*,0,0], linestyle=3



for i=0,n_elements(normal)-1 do if host.snname_array[normal[i]] ne 'SN2008Q' and  host.snname_array[normal[i]] ne 'SN2008hv'  and  host.snname_array[normal[i]] ne  'SN2011by'  and  host.snname_array[normal[i]] ne  'SN2011fe'  and  host.snname_array[normal[i]] ne  'SN2011ia' then cgplots, bpeak_array[4,i]-bpeak_array[5,i], bpeak_array[2,i]-bpeak_array[5,i], psym=16, color=colortable[i], symsize=1.0


index11fe=where(host.snname_array[normal] eq 'SN2011fe')
index11ia=where(host.snname_array[normal] eq 'SN2011ia')
index11by=where(host.snname_array[normal] eq 'SN2011by')
index08Q=where(host.snname_array[normal] eq 'SN2008Q')
index08hv=where(host.snname_array[normal] eq 'SN2008hv')

cgplots, bpeak_array[4,index11ia]-bpeak_array[5,index11ia], bpeak_array[2,index11ia]-bpeak_array[5,index11ia], psym=34, color=colortable[index11ia], symsize=1.0
cgplots, bpeak_array[4,index11fe]-bpeak_array[5,index11fe], bpeak_array[2,index11fe]-bpeak_array[5,index11fe], psym=46, color=colortable[index11fe], symsize=1.5
cgplots, bpeak_array[4,index11by]-bpeak_array[5,index11by], bpeak_array[2,index11by]-bpeak_array[5,index11by], psym=17, color=colortable[index11by], symsize=1.0
cgplots, bpeak_array[4,index08hv]-bpeak_array[5,index08hv], bpeak_array[2,index08hv]-bpeak_array[5,index08hv], psym=15, color=colortable[index08hv], symsize=1.0
cgplots, bpeak_array[4,index08Q]-bpeak_array[5,index08Q], bpeak_array[2,index08Q]-bpeak_array[5,index08Q], psym=18, color=colortable[index08Q], symsize=1.0

al_legend, ['SN2008Q','SN2008hv','SN2011by','SN2011fe','SN2011ia'], color=[colortable[index08Q],colortable[index08hv],colortable[index11by],colortable[index11fe],colortable[index11ia]], psym=[18,15,17,46,34], position=[-0.1, 1.8], charsize=0.6, box=0


 oploterror, bpeak_array[4,*]-bpeak_array[5,*], bpeak_array[2,*]-bpeak_array[5,*], sqrt(bpeak_err_array[4,*]^2.0+bpeak_err_array[5,*]^2.0), sqrt(bpeak_err_array[2,*]^2.0+bpeak_err_array[5,*]^2.0), psym=16, color='black', symsize=0.1


;;;;;;;;;;;;;;;;;;;;;;
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2,x2,y2+y2-y1], $
 ytitle='w1!BB!Lpeak !N- $\mu$', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[-15.0,-19], ystyle=1, xrange=[-0.1,0.3], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=4, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color='black'


cgoplot, color='black', feredmags[4,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[2,*,febpeak,3]-fedm, linestyle=0
cgoplot, color='black', feredmags[4,*,febpeak,2]-feredmags[5,*,febpeak,2], feredmags[2,*,febpeak,2]-fedm, linestyle=1
cgoplot, color='black', feredmags[4,*,febpeak,1]-feredmags[5,*,febpeak,1], feredmags[2,*,febpeak,1]-fedm, linestyle=2
cgoplot, color='black', feredmags[4,*,febpeak,0]-feredmags[5,*,febpeak,0], feredmags[2,*,febpeak,0]-fedm, linestyle=3


al_legend, ['R!LV!N=3.1','R!LV!N=1.7','SMC', 'CSLMC'], linestyle=[0,1,2,3], position=[0.1,-18.9], charsize=0.8, box=0, background='white'

;oploterror, bpeak_array[4,w1vsne]-bpeak_array[5,w1vsne], bpeak_array[2,w1vsne]-host.dm_best_array[normal[w1vsne]], sqrt(bpeak_err_array[4,w1vsne]^2.0+bpeak_err_array[5,w1vsne]^2.0), sqrt(host.dm_best_err_array[normal[w1vsne]]^2.0+bpeak_err_array[5,w1vsne]^2.0), psym=15

for i=0,n_elements(normal)-1 do if host.snname_array[normal[i]] ne 'SN2008Q' and  host.snname_array[normal[i]] ne 'SN2008hv'  and  host.snname_array[normal[i]] ne  'SN2011by'  and  host.snname_array[normal[i]] ne  'SN2011fe'  and  host.snname_array[normal[i]] ne  'SN2011ia' then  cgplots, bpeak_array[4,i]-bpeak_array[5,i], bpeak_array[2,i]-host.dm_best_array[normal[i]], psym=16, color=colortable[i], symsize=1.0

cgplots, bpeak_array[4,index11ia]-bpeak_array[5,index11ia], bpeak_array[2,index11ia]-host.dm_best_array[normal[index11ia]], psym=34, color=colortable[index11ia], symsize=1.0
cgplots, bpeak_array[4,index11fe]-bpeak_array[5,index11fe], bpeak_array[2,index11fe]-host.dm_best_array[normal[index11fe]], psym=46, color=colortable[index11fe], symsize=1.5
cgplots, bpeak_array[4,index11by]-bpeak_array[5,index11by], bpeak_array[2,index11by]-host.dm_best_array[normal[index11by]], psym=17, color=colortable[index11by], symsize=1.0
cgplots, bpeak_array[4,index08hv]-bpeak_array[5,index08hv], bpeak_array[2,index08hv]-host.dm_best_array[normal[index08hv]], psym=15, color=colortable[index08hv], symsize=1.0
cgplots, bpeak_array[4,index08Q]-bpeak_array[5,index08Q], bpeak_array[2,index08Q]-host.dm_best_array[normal[index08Q]], psym=18, color=colortable[index08Q], symsize=1.0

oploterror, bpeak_array[4,*]-bpeak_array[5,*], bpeak_array[2,*]-host.dm_best_array[normal[*]], sqrt(bpeak_err_array[4,*]^2.0+bpeak_err_array[5,*]^2.0), sqrt(host.dm_best_err_array[normal[*]]^2.0+bpeak_err_array[2,*]^2.0), psym=16, color='black', symsize=0.1



cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='!A(b-v)!NBpeak       ',   ytitle='v!BB!Lpeak !N- $\mu$', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[-17.0,-20], ystyle=1, xrange=[-0.1,0.3], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=3, ytickv=ytickvalues, color='black'



cgoplot, color='black', feredmags[4,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[5,*,febpeak,3]-fedm, linestyle=0
cgoplot, color='black', feredmags[4,*,febpeak,2]-feredmags[5,*,febpeak,2], feredmags[5,*,febpeak,2]-fedm, linestyle=1
cgoplot, color='black', feredmags[4,*,febpeak,1]-feredmags[5,*,febpeak,1], feredmags[5,*,febpeak,1]-fedm, linestyle=2
cgoplot, color='black', feredmags[4,*,febpeak,0]-feredmags[5,*,febpeak,0], feredmags[5,*,febpeak,0]-fedm, linestyle=3


;oploterror, bpeak_array[4,w1vsne]-bpeak_array[5,w1vsne], bpeak_array[5,w1vsne]-host.dm_best_array[normal[w1vsne]], sqrt(bpeak_err_array[4,w1vsne]^2.0+bpeak_err_array[5,w1vsne]^2.0), sqrt(host.dm_best_err_array[normal[w1vsne]]^2.0+bpeak_err_array[5,w1vsne]^2.0), psym=16

cgplots, bpeak_array[4,index11ia]-bpeak_array[5,index11ia], bpeak_array[5,index11ia]-host.dm_best_array[normal[index11ia]], psym=34, color=colortable[index11ia], symsize=1.0
cgplots, bpeak_array[4,index11fe]-bpeak_array[5,index11fe], bpeak_array[5,index11fe]-host.dm_best_array[normal[index11fe]], psym=46, color=colortable[index11fe], symsize=1.5
cgplots, bpeak_array[4,index11by]-bpeak_array[5,index11by], bpeak_array[5,index11by]-host.dm_best_array[normal[index11by]], psym=17, color=colortable[index11by], symsize=1.0
cgplots, bpeak_array[4,index08hv]-bpeak_array[5,index08hv], bpeak_array[5,index08hv]-host.dm_best_array[normal[index08hv]], psym=15, color=colortable[index08hv], symsize=1.0
cgplots, bpeak_array[4,index08Q]-bpeak_array[5,index08Q], bpeak_array[5,index08Q]-host.dm_best_array[normal[index08Q]], psym=18, color=colortable[index08Q], symsize=1.0

for i=0,n_elements(normal)-1 do if host.snname_array[normal[i]] ne 'SN2008Q' and  host.snname_array[normal[i]] ne 'SN2008hv'  and  host.snname_array[normal[i]] ne  'SN2011by'  and  host.snname_array[normal[i]] ne  'SN2011fe'  and  host.snname_array[normal[i]] ne  'SN2011ia' then  cgplots, bpeak_array[4,i]-bpeak_array[5,i], bpeak_array[5,i]-host.dm_best_array[normal[i]], psym=16, color=colortable[i], symsize=1.0



oploterror, bpeak_array[4,*]-bpeak_array[5,*], bpeak_array[5,*]-host.dm_best_array[normal[*]], sqrt(bpeak_err_array[4,*]^2.0+bpeak_err_array[5,*]^2.0), sqrt(host.dm_best_err_array[normal[*]]^2.0+bpeak_err_array[5,*]^2.0), psym=16, color='black', symsize=0.1


device, /close
SET_PLOT, 'X'


spawn, 'open w1v_bpeak_colors_abs_color.eps'


;for n=0,n_elements(normal)-1 do print, host.snname_array[normal[n]], bpeak_array[5,n]-host.dm_best_array[normal[n]], sqrt(host.dm_best_err_array[normal[n]]^2.0+bpeak_err_array[5,n]^2.0)

;;;;;;;;;;;;;;;;;;;;;
SET_PLOT, 'PS'

device, filename='w1v_bpeak_colors_abs_dmb15_nocolorbar.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, /color

xrange=[-0.2,0.2]
nxticks=4
                                                  
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2+y2-y1,x2,y2+y2+y2-y1-y1], $
xrange=[1,1.4], yrange=[3.0,0.5], ytitle='(w1 - v)!BB!Lpeak',  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
ystyle=1, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color='black'


for i=0,n_elements(normal)-1 do if host.snname_array[normal[i]] ne 'SN2008Q' and  host.snname_array[normal[i]] ne 'SN2008hv'  and  host.snname_array[normal[i]] ne  'SN2011by'  and  host.snname_array[normal[i]] ne  'SN2011fe'  and  host.snname_array[normal[i]] ne  'SN2011ia' then  cgplots, host.dm15_array[normal[i],4], bpeak_array[2,i]-bpeak_array[5,i], psym=16, color=colortable[i], symsize=1.0


cgplots, host.dm15_array[normal[index11ia],4], bpeak_array[2,index11ia]-bpeak_array[5,index11ia], psym=34, color=colortable[index11ia], symsize=1.0
cgplots, host.dm15_array[normal[index11fe],4], bpeak_array[2,index11fe]-bpeak_array[5,index11fe], psym=46, color=colortable[index11fe], symsize=1.5
cgplots, host.dm15_array[normal[index11by],4], bpeak_array[2,index11by]-bpeak_array[5,index11by], psym=17, color=colortable[index11by], symsize=1.0
cgplots, host.dm15_array[normal[index08hv],4], bpeak_array[2,index08hv]-bpeak_array[5,index08hv], psym=15, color=colortable[index08hv], symsize=1.0
cgplots, host.dm15_array[normal[index08Q],4],  bpeak_array[2,index08Q] -bpeak_array[5,index08Q],  psym=18, color=colortable[index08Q], symsize=1.0

oploterror, host.dm15_array[normal,4], bpeak_array[2,*]-bpeak_array[5,*], host.dm15err_array[normal,4], sqrt(bpeak_err_array[2,*]^2.0+bpeak_err_array[5,*]^2.0), psym=16, symsize=0.1, color='black'

;;;;;;;;;;;;;;;;;;;;;;
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2,x2,y2+y2-y1], $
 ytitle='w1!BB!Lpeak !N- $\mu$', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[-15.0,-19], ystyle=1, xrange=[1,1.4], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=4, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color='black'

for i=0,n_elements(normal)-1 do if host.snname_array[normal[i]] ne 'SN2008Q' and  host.snname_array[normal[i]] ne 'SN2008hv'  and  host.snname_array[normal[i]] ne  'SN2011by'  and  host.snname_array[normal[i]] ne  'SN2011fe'  and  host.snname_array[normal[i]] ne  'SN2011ia' then  cgplots, host.dm15_array[normal[i],4], bpeak_array[2,i]-host.dm_best_array[normal[i]], psym=16, color=colortable[i], symsize=1.0

cgplots, host.dm15_array[normal[index11ia],4], bpeak_array[2,index11ia]-host.dm_best_array[normal[index11ia]], psym=34, color=colortable[index11ia], symsize=1.0
cgplots, host.dm15_array[normal[index11fe],4], bpeak_array[2,index11fe]-host.dm_best_array[normal[index11fe]], psym=46, color=colortable[index11fe], symsize=1.5
cgplots, host.dm15_array[normal[index11by],4], bpeak_array[2,index11by]-host.dm_best_array[normal[index11by]], psym=17, color=colortable[index11by], symsize=1.0
cgplots, host.dm15_array[normal[index08hv],4], bpeak_array[2,index08hv]-host.dm_best_array[normal[index08hv]], psym=15, color=colortable[index08hv], symsize=1.0
cgplots, host.dm15_array[normal[index08Q],4], bpeak_array[2,index08Q]-host.dm_best_array[normal[index08Q]], psym=18, color=colortable[index08Q], symsize=1.0

oploterror, host.dm15_array[normal,4], bpeak_array[2,*]-host.dm_best_array[normal[*]], host.dm15err_array[normal,4], sqrt(host.dm_best_err_array[normal[*]]^2.0+bpeak_err_array[5,*]^2.0), psym=16, symsize=0.1, color='black'

cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='!A$\Delta$ M!N15!A(B)',   ytitle='v!BB!Lpeak !N- $\mu$',  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[-17.0,-20], ystyle=1, xrange=[1,1.4], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=3, ytickv=ytickvalues, color='black'

for i=0,n_elements(normal)-1 do if host.snname_array[normal[i]] ne 'SN2008Q' and  host.snname_array[normal[i]] ne 'SN2008hv'  and  host.snname_array[normal[i]] ne  'SN2011by'  and  host.snname_array[normal[i]] ne  'SN2011fe'  and  host.snname_array[normal[i]] ne  'SN2011ia' then  cgplots, host.dm15_array[normal[i],4], bpeak_array[5,i]-host.dm_best_array[normal[i]], psym=16, color=colortable[i], symsize=1.0

cgplots, host.dm15_array[normal[index11ia],4], bpeak_array[5,index11ia]-host.dm_best_array[normal[index11ia]], psym=34, color=colortable[index11ia], symsize=1.0
cgplots, host.dm15_array[normal[index11fe],4], bpeak_array[5,index11fe]-host.dm_best_array[normal[index11fe]], psym=46, color=colortable[index11fe], symsize=1.5
cgplots, host.dm15_array[normal[index11by],4], bpeak_array[5,index11by]-host.dm_best_array[normal[index11by]], psym=17, color=colortable[index11by], symsize=1.0
cgplots, host.dm15_array[normal[index08hv],4], bpeak_array[5,index08hv]-host.dm_best_array[normal[index08hv]], psym=15, color=colortable[index08hv], symsize=1.0
cgplots, host.dm15_array[normal[index08Q],4],  bpeak_array[5,index08Q] -host.dm_best_array[normal[index08Q]],  psym=18, color=colortable[index08Q], symsize=1.0

oploterror, host.dm15_array[normal,4],  bpeak_array[5,*]-host.dm_best_array[normal[*]], host.dm15err_array[normal,4],  sqrt(host.dm_best_err_array[normal[*]]^2.0+bpeak_err_array[5,*]^2.0), psym=16, symsize=0.1, color='black'

device, /close
SET_PLOT, 'X'

spawn, 'open w1v_bpeak_colors_abs_dmb15_nocolorbar.eps'




SET_PLOT, 'PS'

device, filename='w1v_bpeak_colors_abs_dmb15.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, /color

nxticks=4
                                                  
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2+y2-y1,x2,y2+y2+y2-y1-y1], $
xrange=[1,1.4], yrange=[3.0,0.5], ytitle='(w1-v)!BB!Lpeak',$
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
ystyle=1, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color=black


for i=0,n_elements(normal)-1 do cgplots, host.dm15_array[normal[i],4], bpeak_array[2,i]-bpeak_array[5,i], psym=16, color=colortable[i], symsize=1.0


oploterror, host.dm15_array[normal,4], bpeak_array[2,*]-bpeak_array[5,*], host.dm15err_array[normal,4], sqrt(bpeak_err_array[2,*]^2.0+bpeak_err_array[5,*]^2.0), psym=16, symsize=0.1

;cgLOADCT, 33, NCOLORS=100
cgcolorbar, ncolors=256 , POSITION=[x1,y2+y2+y2-y1-y1, x2, 0.95], range=[-0.1,0.3], /top, charsize=1

;;;;;;;;;;;;;;;;;;;;;;
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2,x2,y2+y2-y1], $
 ytitle='w1!BB!Lpeak !N- $\mu$', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[-16.0,-19], ystyle=1, xrange=[1,1.4], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=3, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color=black


for i=0,n_elements(normal)-1 do cgplots, host.dm15_array[normal[i],4], bpeak_array[2,i]-host.dm_best_array[normal[i]], psym=16, color=colortable[i], symsize=1.0



oploterror, host.dm15_array[normal,4], bpeak_array[2,*]-host.dm_best_array[normal[*]], host.dm15err_array[normal,4], sqrt(host.dm_best_err_array[normal[*]]^2.0+bpeak_err_array[5,*]^2.0), psym=16, symsize=0.1


cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='!A$\Delta$ M!N15!A(B)',   ytitle='v!BB!Lpeak !N- $\mu$', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[-17.0,-20], ystyle=1, xrange=[1,1.4], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=3, ytickv=ytickvalues, color=black

for i=0,n_elements(normal)-1 do cgplots, host.dm15_array[normal[i],4], bpeak_array[5,i]-host.dm_best_array[normal[i]], psym=16, color=colortable[i], symsize=1.0


oploterror, host.dm15_array[normal,4],  bpeak_array[5,*]-host.dm_best_array[normal[*]], host.dm15err_array[normal,4],  sqrt(host.dm_best_err_array[normal[*]]^2.0+bpeak_err_array[5,*]^2.0), psym=16, symsize=0.1


device, /close
SET_PLOT, 'X'


spawn, 'open w1v_bpeak_colors_abs_dmb15.eps'



print, 'final stop'
stop

end
