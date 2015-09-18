# Purpose: Makefile demonstrates problem with environment variable expansions
# Results differ on IRIX gmake 3.79.1 vs. GNU/Linux 3.79.1
# IRIX version prints value of $HOME environment variable and ignores Makefile modification
# Linux version prints $HOME environment variable as modified in Makefile.tst
# Usage: make -f Makefile.tst
HOME := ${HOME}_plus_foo_set_in_Makefile
tst:
	@printf "HOME = ${HOME}\n"
