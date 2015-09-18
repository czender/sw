************************************************************************
                    THIS IS UNSUPPORTED SOFTWARE
    but documented bug reports and enhancements are very welcome.
************************************************************************

IBP 3.7 -- An interactive graphical analyzer for GCM data in netCDF
format 

Purpose:
o IBP runs under IDL and creates graphs from your netCDF GCM
datafiles.  
o IBP is purely widget-driven and is intended for color workstations. 
o IBP displays color or black and white graphs on X-window, GIF, and
Postscript output devices; the hardcopy is generally publication
quality. 
o IBP has a few built-in data processing features including: zonal
averages of single and multi-level fields, area-averages of lat-lon
regions, field differencing, linear regression analysis, and
geographic mapping of points in a histogram frequency bin or
scatterplot rectangle.  
o IBP is fairly flexible as far as axis labelling, contour-levels, etc.

IBP can display 6 major graph types: 
1. color fill or color image cylindrical projection contour maps
2. color image polar azimuthal projection contour maps
3. color fill lat-lev zonal average contours
4. B&W single level zonal averages (up to 2 fields overplotted) from a
lat-lon region 
5. B&W histograms from a lat-lon region
6. B&W scatterplots from a lat-lon region

CAVEATS: 1. IBP does some automatic scaling of fields with common
names---it assumes the input field is in CCM2 output dimensions (SI,
except for cloud water paths), and scales the field to user-friendly
units (e.g., mm/day) automatically. IBP also has many built-in
defaults that attempt to pick standard axes ranges and plot labels
based on file and field names. These defaults are intended to be
user-configurable (see ibp_sct.pro file), but in practice you may
just want to work around them if they get in the way. Pre-defined
scaling and defaults are very convenient for producing quick and dirty
graphs but often cause confusion when user's attempt to pre-scale the
data themselves using the CCM processor, for example.

2. IBP is not documented, but is fairly intuitive to use once the
basic button sequences have been learned. Until you master this
ordering, expect IBP to crash.

3. Although the text in all the labels is widget-configurable, text
size and location is fixed.

4. IBP does not support regional datasets, the input should be global.
The results of inputting regional datasets are unpredictable.  IBP
works at any truncation type and resolution except for the area
averaging (which depends on the gaussian weights) which has only been
tested for T42 data.  IBP assumes the existence of dimensions named
'lat' and 'lon' (and 'lev' and 'time' if you're working with > 2
dimensional data).

5. IBP is written in IDL, an interpreted language, and thus adding (or
removing) features (like new widgets, graph types, defaults) is fairly
painless if you know IDL.

6. An inexhaustive list of some nice NCAR graphics features missing
from IBP: 
 o high/lo labels
 o cross-hatching
 o batch processing/macros
 o user-defined regional color fills
 o graph windowing (combining multiple graphs onto one page)
 o user-selected font size
 o user-selected X window size
 o NCAR graphics is standalone; IBP requires proprietary IDL software.

7. No one knows what IBP stands for, but possibilities include
   Itty Bitty Processor 
   Illustrating Bad Programming
   IDL Bug Producer
   Instead of Beautiful Plots
   I Better Procrastinate

NCAR/CGD local user notes:
/data2/zender/pub/src/ibp* contains all the files needed to run IBP.
The files in this distribution directory will generally be more bug
free than the files in the development directory, /home/zender/idl,
but you are welcome to snoop around there for any recent
changes. Files in this distribution directory will be updated after
each RCS revision.

Converting history tapes to netCDF format:

ccm2nc is Brian Eaton's history tape to netCDF translator.
the syntax is ccm2nc hist_tape netCDF_file, e.g.,
ccm2nc /CCM2/T42/%data%/SEP1 /usr/tmp/joeschmoe/SEP1.nc
AIX and UNICOS versions of ccm2nc are located in
/crestone/u1/eaton/bin/AIX/ccm2nc and
alpine.ucar.edu:/home/alpine1/eaton/bin/ccm2nc, respectively

cond, maintained by John Truesdale, performs netCDF translation with
more batch-oriented capabilities. the Cray executable is located in
/crestone/u1/jet/craybin/condnew < cond.$$ and a sample cond
invocation from a UNIX shell script follows: 
#
cat > cond.$$ << 'END1'
;Get initial data history tape into netCDF format
E$HIST
     FORMAT  = 'netcdf'
     DELINP  = .true.
     MSDIRI  = '/CCM2/T42/%data%'
     FILESI  = 'SEP1'
     RCPDIR  = 'ra.cgd.ucar.edu:/data2/zender/ccm'
     FILESO  = 'SEP1.nc'
     DAYFRST = -1.0
     NTHIST  = 50
     RETPD   = 365
     WRPSWD  = ' ' $
'END1'
#
# Execute program
#
/crestone/u1/jet/craybin/condnew < cond.$$
#

INSTALLATION:

Files:
ibp-3.7.tar.gz			the entire source tree:
./ibp-3.7/README.ibp		this file
./ibp-3.7/ibp.pro		top level procedure to load all other procedures
./ibp-3.7/ibp_srt.txt		important startup file for IDL (see below)
./ibp-3.7/ibp_widgets.pro	ibp main procedure and event handler.
./ibp-3.7/ibp_Xdefaults.txt	example IBP X resources for ~/.Xdefaults
./ibp-3.7/ibp_sct.pro		user editable defaults 
./ibp-3.7/ibp_commons.com	ibp common blocks
./ibp-3.7/ibp_figures.pro	most data manipulation utilities
./ibp-3.7/ibp_analyze.pro	mostly obsolete
./ibp-3.7/ibp_anim.pro		animation capabilities
./ibp-3.7/ibp_average.pro	gaussian weight averager
./ibp-3.7/ibp_colors.com	common block for colors		
./ibp-3.7/ibp_colors.pro	color generation and color bar procedures
./ibp-3.7/ibp_custom.pro	customized graphics procedures
./ibp-3.7/ibp_functions.pro	utility and printer procedures
./ibp-3.7/ibp_hist.pro		histogram procedures
./ibp-3.7/ibp_lat_lev.pro	latitude-level procedures
./ibp-3.7/ibp_lon_lat.pro	procedures for various map-styles
./ibp-3.7/ibp_map_set.pro	IDL 3.1 map_set() procedures
./ibp-3.7/ibp_scat.pro		scatterplot procedures
./ibp-3.7/ibp_time_lat.pro	time latitude (Hovmuller) procedures
./ibp-3.7/ibp_zon_avg.pro	zonal average (line figure) procedures
./ibp-3.7/SEP1.nc		default data file for IBP
./ibp-3.7/*.tbl			example user-defined color tables

You can either install IBP from the compressed tarfile,
ibp-3.7.tar.gz, or by copying all the source files to your own working
directory, $MY_IBP_DIR.  The compressed tarfile ibp-3.7.tar.gz must be
uncompressed (with GNU gunzip) and untarred (tar -xvf ibp-3.7.tar).
The files will be untarred into a subdirectory called ibp-3.7. If you
do not wish this to be your working directory, then move all the files
to your working directory, called $MY_IBP_DIR.  In the following
examples, $MY_IBP_DIR= /data2/zender/pub/src/ibp-3.7 (remember to use
the full pathname).  Set an environment variable in your .cshrc:

	setenv IDL_STARTUP $MY_IBP_DIR/ibp_startup.txt

This file, ibp_startup.txt, contains a command to increase the default
amount of memory IDL uses to hold code. This is necessary because the
IBP package is a lot of code.  You must use IDL version 3.6 or
higher. You can use packages to make sure you get version 3.6 on your
path, or just do it in your .cshrc with

	setenv IDL_DIR /opt/idl3.6
	setenv IDL_PATH +/opt/idl3.6/lib

If you encounter errors such as "unknown osf keysym" (typically under
SunOS with twm window manager), then you need to set an additional
environment variable to tell IDL where to find the keysym database,

	setenv XKEYSYMDB /usr/local/lib/X11/XKeysymDB

Other IDL installation hints are documented in the IDL distribution
subdirectory notes, e.g., /opt/idl3.6/notes/sun.doc

There are some example X-resources for IBP which you might choose to
install in your ~/.Xdefaults file located in the ibp_Xdefaults.txt
file.  These resources add color highlighting (and, eventually,
keyboard accelators) to the IBP menu-buttons. If you make a nicer set
of these, please let me know so I can include them in the next
distribution.

RUNNING IBP:

Your are now ready to start IBP. The recommended way to run IBP is to
change directories to $MY_IBP_DIR and then start IDL by typing "idl":

/data2/zender/pub/src/ibp-3.6: idl

IDL will run in the shell window you execute it from, this is where
all the error messages and diagnostics will appear.  When you get to
the IDL prompt, type "ibp":

IDL> ibp

After about a minute the IBP window will come up.  At this point you
shouldn't need the IDL shell window anymore, so you can iconize it. If
IBP crashes or you want to quit the IDL session, just type "exit" in
the IDL window:

IDL> exit
/data2/zender/pub/src/ibp-3.6:

If IBP crashes and you know why it crashed, it's much easier to
restart IBP by typing 'xm' at the IDL> prompt. 'xm' will restart the
event manager and you can usually pick up exactly where you were
before the crash. This is much quicker method of crash recovery than
exiting IDL and then restarting IDL and IBP.

IBP must have a file to load when you start it, the filename should
end in .nc to be viewed by the default configuration of the file
selector widget.  The default file is currently SEP1.nc, the canonical
CCM2 initial dataset supplied with the IDP distribution.  You can
change this default at startup by giving IBP a filename with the file
keyword:

ibp,file='/data2/zender/ccm/388_jul.nc'

The single quotes are required, the pathname is relative to the
current directory, so '388_jul.nc' could also work.  To change the
default start file permanently, edit the file ibp_widgets.pro replace
the string SEP1.nc with whatever file you want as default. This is
recommended since it allows you to specify a default path to your data
that is separate from the directory with all the IBP source files.

When you are working in production mode, it is helpful to have a version
of ghostview running in the background on the file $MY_IDP_DIR/idl.ps.
You can preview all your graphics with this ghostview, which is much
faster than spawning a new ghostview from within IBP.




