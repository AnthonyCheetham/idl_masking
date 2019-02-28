function copy_oitable, oitab
;
; Make a copy of oi-table, making sure that the associated heap
; variables are also copied, so that even if the original table was
; erased, the copy will remain valid.
;
; 2012-10-12  M. Kishimoto
;

newoitab = oitab

for k=0, n_elements(oitab)-1 do begin
  for i=0, n_tags(oitab)-1 do begin
    if size(oitab[k].(i),/tname) eq 'POINTER' then begin
      ;
      ; this will make a copy of the heap variable *oitab[k].(i)
      ;
      newoitab[k].(i) = ptr_new(*oitab[k].(i))
    endif
  endfor; loop over each structure member
endfor; loop over each table

return, newoitab
end

