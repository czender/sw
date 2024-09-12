// $Id$

// Purpose: Demonstrate portable internationalization i18n with gettext()

/* Copyright (C) 1997--present Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* Notes: i18n is maintained separately from ccc for convenience
   internationalization i18n stands for "i followed by 18 letters followed by n"
   localization l10n stands for "l followed by 10 letters followed by n"
   ccc is important to always have running on any architecture
   But non-Linux/Solaris architectures often have no gettext() support
   Also, i18n requires special Makefile support (msgfmt) for message catalogs and libraries (libintl) */

/* Usage:
   Extract strings into PO (portable object) file:
   xgettext --help
   xgettext --default-domain=i18n --join-existing i18n.cc
   Generate template PO (portable object) file:
   xpot --help
   xpot --output-file=${HOME}/c++/i18n.po i18n.cc
   Compile PO (portable object) file to MO (message object or machine object) file:
   msgfmt --help
   msgfmt --output-file=${HOME}/share/locale/es/LC_MESSAGES/i18n.mo --statistics ${HOME}/c++/i18n.po
   Compile code:
   g++ -I. -I/usr/local/include -O -Wall -c i18n.cc -o ${HOME}/obj/LINUX/i18n.o
   g++ -o /home/zender/bin/LINUX/i18n /home/zender/obj/LINUX/i18n.o
   fxm: 20010114 Program does not compile with -g
   Change language:
   The LANG token seems to automagically change all the LC_* tokens correspondingly
   export LANG='es'
   export LANG='en_US'
   export LOCALE='POSIX' # Locale
   export LOCALE='C' # Locale

   i18n Localization individual tokens
   LC_CTYPE # Regular expressions and case-modification functions
   LC_NUMERIC # Formatting functions printf(), sprintf(), write()
   LC_TIME # POSIX date formatting function strftime()
   NLS stands for "Native Language Support"

   gettext() installation is its own best example:
   /usr/share/gettext/intl

   Newsgroups:
   comp.software.international (no activity)
   comp.std.internat (no activity)

   i18n links:
   ftp://sunsite.unc.edu/pub/Linux/utils/nls/catalogs/Incoming/locale-tutorial-0.8.txt.gz (good introduction to C locale using catgets(), but no gettext())
   ftp://ftp.ora.com/pub/examples/nutshell/ujip/doc/i18n-books.txt (Way, way, way out of date and useless)
   http://www.uni-ulm.de/~s_smasch/Locale/ (appears to be very useful, specific, need to read thoroughly) */

// Standard C++ headers
#include <iostream> // Standard C++ I/O streams cout, cin
#include <string> // Standard C++ string class

// Internationalization i18n
// libintl.h header is required in every file with gettext() functions
#if ( !defined ALPHA ) && ( !defined MACOS ) && ( !defined SGI6 ) && ( !defined SGI64 ) && ( !defined SGIMP64 )
# include <libintl.h> // Internationalization i18n
#endif // !OS has libintl.h
// Linux and SGI libintl.h load locale.h themselves, but Solaris libintl.h does not so do it manually
#include <locale.h> // Locale setlocale()
#ifndef _LIBINTL_H
// Define stub for gettext to allow compiling when libintl.h is not available
# define gettext(foo) foo
#endif // !_LIBINTL_H

// Standard C headers
#include <cstdlib> // strtod, strtol, malloc, getopt, getenv

int main(int argc,char **argv)
{
  /* RedHat GNU/Linux defaults to en_US for all LC_* tokens
     Instead of LC_ALL, it is possible to set LC_* tokens individually */

  std::string nvr_HOME((std::getenv("HOME")) ? std::getenv("HOME") : ""); // [sng] Environment variable HOME
  std::string lcl_dir(nvr_HOME+"/share/locale"); // [sng] LOCALEDIR is e.g., /usr/share/locale
#ifdef _LIBINTL_H
  // Next three lines invoke i18n, are same for every program
  setlocale(LC_ALL,""); // LC_ALL sets all localization tokens to same value
  // gettext() manual recommends supplying PACKAGE and LOCALEDIR in Makefile or config.h
  // bindtextdomain(PACKAGE,LOCALEDIR); // LOCALEDIR is e.g., /usr/share/locale
  bindtextdomain("i18n",lcl_dir.c_str()); // LOCALEDIR is e.g., /usr/share/locale
  // MO files should be in $LOCALEDIR/es/LC_MESSAGES
  textdomain("i18n"); // PACKAGE is name of program
#endif // not _LIBINTL_H
  std::cout << "Testing internationalization..." << std::endl;
  std::cout << "export LANG=\'es\';i18n;export LANG=\'en_US\';export LC_ALL=\'C\'" << std::endl;
  std::cout << "setenv LANG \'es\';i18n;setenv LANG \'en_US\';setenv LC_ALL \'C\'" << std::endl;
  std::cout << gettext("The tree is green") << std::endl;
  std::cout << gettext("I thought we were friends") << std::endl;

  if(argc == argc+1){;} // CEWU Compiler Error Warning Usage
  if(argv != (char **)73){;} // CEWU Compiler Error Warning Usage
} // !main()
  
