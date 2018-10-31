;A function for phase chi^2 in the presence of wrapping
function phase_chi2_gpi, p
common phasestuff, fitmat, ph_mn, ph_err
if (size(ph_mn))[0] eq 2 then dim3=(size(ph_mn))[2] else dim3=1
chi2=0
for wav=0,dim3-1 do begin
   chi2+=modsq(1.-exp(complex(0,[ph_mn[*,wav],0] - reform(p#fitmat))))/[ph_err[*,wav],0.01]^2
endfor
return,total(chi2)
;return, total(modsq(1.-exp(complex(0,ph_mn - reform(p#fitmat))))/ph_err^2)
end
