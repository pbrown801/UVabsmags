restore, 'host.sav', verbose=1

for n=0,n_elements(host.snname_array)-1 do if finite(host.appmag_array[n,0]) eq 1 then  print, host.snname_array[n], host.dm15_array[n,4], host.dm15err_array[n,4],                       host.appmag_array[n,0]-host.dm_best_array[0], host.appmagerr_array[n,0], host.appmag_array[n,1]-host.dm_best_array[0], host.appmagerr_array[n,1], host.appmag_array[n,2]-host.dm_best_array[0], host.appmagerr_array[n,2], host.appmag_array[n,3]-host.dm_best_array[0], host.appmagerr_array[n,3], host.appmag_array[n,4]-host.dm_best_array[0], host.appmagerr_array[n,4], host.appmag_array[n,5]-host.dm_best_array[0], host.appmagerr_array[n,5], host.bpeakappmag_array[n,0]-host.dm_best_array[0], host.bpeakappmagerr_array[n,0], host.bpeakappmag_array[n,1]-host.dm_best_array[0], host.bpeakappmagerr_array[n,1], host.bpeakappmag_array[n,2]-host.dm_best_array[0], host.bpeakappmagerr_array[n,2], host.bpeakappmag_array[n,3]-host.dm_best_array[0], host.bpeakappmagerr_array[n,3], host.bpeakappmag_array[n,4]-host.dm_best_array[0], host.bpeakappmagerr_array[n,4], host.bpeakappmag_array[n,5]-host.dm_best_array[0], host.bpeakappmagerr_array[n,5]
