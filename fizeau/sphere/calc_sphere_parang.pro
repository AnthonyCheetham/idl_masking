function calc_sphere_parang,head
;; IDL function to calculate the parallactic angle from a SPHERE observation
;; This works directly from the RA/Dec as opposed to the parang header keywords,
;; since there were issues with the header values until mid 2016.
;; Two main issues:
;;     1. Proper motion was not included in the RA/Dec used to calculate the parang
;;     2. The computer controlling the derotator was not synched to the correct time,
;;         so it was calculating the wrong value

;; Conversion factors
d2r = !dpi/180.
r2d = 180./!dpi

;; First, get the correct RA and Dec
actual_ra = double(sxpar_conica(head,'4DROT2RA',nchar=8))
actual_dec = double(sxpar_conica(head,'4DROT2DEC',nchar=9))

;; These values were in weird units: HHMMSS.ssss
actual_ra_hr = floor(actual_ra/10000.)
actual_ra_min = floor(actual_ra/100. - actual_ra_hr*100.)
actual_ra_sec = (actual_ra - actual_ra_min*100. - actual_ra_hr*10000.)

actual_ra_deg = (actual_ra_hr + actual_ra_min/60. + actual_ra_sec/60./60.) * 360./24.

;; the sign makes this complicated, so remove it now and add it back at the end
sgn = sign(actual_dec)
actual_dec *= sgn

actual_dec_deg = floor(actual_dec/10000.)
actual_dec_min = floor(actual_dec/100. - actual_dec_deg*100.)
actual_dec_sec = (actual_dec - actual_dec_min*100. - actual_dec_deg*10000.)

actual_dec_deg = (actual_dec_deg + actual_dec_min/60. + actual_dec_sec/60./60.)*sgn
actual_dec *= sgn

;; This is copied from GRAPHIC (the Geneva imaging pipeline), which was borrowed from some code by Arthur Vigan

geolat_rad=sxpar_conica(head,'ESOTELGEOLAT',nchar=12)*d2r

ha_deg=(sxpar_conica(head,'LST',nchar=3)*15./3600.)-actual_ra_deg

;; VLT TCS formula
f1 = cos(geolat_rad) * sin(d2r*ha_deg)
f2 = sin(geolat_rad) * cos(d2r*actual_dec_deg) - cos(geolat_rad) * sin(d2r*actual_dec_deg) * cos(d2r*ha_deg)

parang_deg = -r2d*atan(-f1,f2)

;; Now, check the date, and if prior to 12 July 2016, correct for a separate issue
alt = sxpar_conica(head, 'TELALT',nchar=6)
derot_begin = sxpar_conica(head,'4DROT2BEGIN',nchar=11)
derot_err = atan( tan( (alt-2.*derot_begin)*!pi/180.) )*180./!pi
date_mjd = sxpar_conica(head,'MJD-OBS')

if (date_mjd le 57582.) then begin
    parang_deg += derot_err ;; this number is the approx MJD for 12 July 2016.
    print,'  Correcting for derotator problems due to data being taken before July 2016'
endif
return,parang_deg
end