Color table gallery
http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
Font gallery
http://www.ncl.ucar.edu/Document/Graphics/font_tables.shtml

Local machine examples:
/usr/local/ncarg/lib/ncarg/nclex/xyplot

strings use double quotes
! access named dimensions, u!0="time"
new() makes new variable
delete() deletes old variable
list_files()
list_filevars()
if(x.gt.y) then a else b end if
NB: NCL has no "else if" construction
| follows dimension name to indicate named subscripting 

String variables can be used to reference attributes and coordinates
by enclosing the variable reference within dollar signs '$'
http://www.ncl.ucar.edu/Document/Manuals/Ref_Manual/NclVariables.shtml

Named subscripting does dimension permutation
Change T(lon,lat,lev) to Tnew(lev,lat,lon):
Tnew = T(lev|:,lat|:,lon|:)

conform() does broadcasting
onedtond() does re-dimensioning
dim_avg_Wrap() averages rightmost dimension and retains metadata
dim_sum() 

Assignment:
->	file variable assignment
!	dimension name assignment
&	coordinate variable assignment
@	attribute assignment
Otherwise normal-value-to-variable assignment occurs. 



