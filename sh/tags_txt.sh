# Change paper names:
# old=CaZ091;new=CaZ09;mv ${DATA}/ppr/ppr_${old}.pdf ${DATA}/ppr/ppr_${new}.pdf;mv ${DATA}/ppr/ppr_${old}_csz.pdf ${DATA}/ppr/ppr_${new}_csz.pdf;mv ${DATA}/ppr/ppr_${old}.bib ${DATA}/ppr/ppr_${new}.bib;
# Dust only: mv /var/www/html/ppr/ppr_${old}.pdf /var/www/html/ppr/ppr_${new}.pdf;mv /var/www/html/ppr/ppr_${old}_csz.pdf /var/www/html/ppr/ppr_${new}_csz.pdf

cd ~/sh
# etags `find ~/. -name Makefile -print`
etags -f TAGS ~/tex/*.bib ~/tex/*.tex ~/tex/*.sty ~/crr/*.tex ~/crr/*.sty ~/job/*.tex ~/job/*.html ~/www/*.html

