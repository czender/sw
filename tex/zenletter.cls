% $Id$

% Purpose: Personal letter class for LaTeX2e

% Based on pletter.cls from KoD95 p. 308

% Usage: 
% /bin/cp $HOME/tex/zenletter.cls $TEXMFDIR/tex/latex2e/csz

% NB: To add a new class or style file to the TeX tree, install in csz
% directory and then update $TEXMFDIR/ls-R by running `texhash' as root. 

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{zenletter}
\RequirePackage{ifthen}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{letter}}
\ProcessOptions
\LoadClass[a4paper]{letter}

\setlength{\textheight}{8.25in} \setlength{\textwidth}{6.5in}
\setlength{\oddsidemargin}{0pt} \setlength{\evensidemargin}{0pt}
\setlength{\topmargin}{-20pt} \setlength{\headheight}{36pt}
\setlength{\headsep}{35pt} \setlength{\footskip}{25pt} 

% \DeclareFixedFont{command}{code}{family}{series}{shape}{size}
\DeclareFixedFont{\xiirn}{OT1}{cmr}{m}{n}{12}
\DeclareFixedFont{\xiissn}{OT1}{cmss}{m}{n}{12}
\DeclareFixedFont{\viiisssl}{OT1}{cmss}{m}{sl}{8}
\DeclareFixedFont{\xiisssl}{OT1}{cmss}{m}{sl}{12}
\DeclareFixedFont{\xssn}{OT1}{cmss}{m}{n}{10}
\DeclareFixedFont{\viitin}{OT1}{cmti}{m}{n}{7}
\newfont{\poetfont}{cmdunh10 scaled \magstep3}

% Sum of \leftfield and \rightfield is \textwidth
\newlength{\wholefield} \setlength{\wholefield}{\textwidth} 
\newlength{\leftfield} \setlength{\leftfield}{4.75in} 
\newlength{\rightfield} \setlength{\rightfield}{1.75in}

\newsavebox{\household} 
\sbox{\household}
{\parbox[b]{\leftfield}{\poetfont Charlie, Robynn, Olivia, \& Ruby Zender}}  

\newcommand{\Stroke}[1]{\raisebox{0.7\height}%
        {\underline{\textit{\tiny #1}}}}

\newsavebox{\householdaddress}
\sbox{\householdaddress}
{\parbox[b]{\rightfield}{\viiisssl\baselineskip10pt
%Tel.: (303) 417-9834\\
%717 9th St.\\
%%717 9\Stroke{th} St.\\
%%717 9$^{\hbox{th}}$ St.\\
%%717 9$^{\underline{\text{th}}}$ St.\\
%%717 9\textsuperscript{\underline{\itshape th}} St.\\
%Boulder CO~~80302}}
%Tel.: (303) 417-9834\\
%2004A Los Trancos Dr.\\
3 Whistler Ct.\\
%%717 9\Stroke{th} St.\\
%%717 9$^{\hbox{th}}$ St.\\
%%717 9$^{\underline{\text{th}}}$ St.\\
%%717 9\textsuperscript{\underline{\itshape th}} St.\\
Irvine CA~~92617}}

% This header places the short address next to the banner header
\newsavebox{\letterhead} 
\sbox{\letterhead}
{\parbox{\textwidth}{\usebox{\household}%
{\usebox{\householdaddress}}\\
\rule{\textwidth}{1pt}}} 

\renewcommand{\ps@firstpage}
{\setlength{\headheight}{41pt}\setlength{\headsep}{25pt}%
\renewcommand{\@oddhead}{\usebox{\letterhead}}%
\renewcommand{\@oddfoot}{}%
\renewcommand{\@evenhead}{}\renewcommand{\@evenfoot}{}}

\renewcommand{\ps@headings}
{\setlength{\headheight}{12pt}%
\renewcommand{\@oddhead}
{\parbox{\textwidth}{
\slshape \headtoname{} \toname\hfill\@date\hfill 
\pagename{} \thepage\\
\rule[3pt]{\textwidth}{1pt}}}
\renewcommand{\@oddfoot}{}
\renewcommand{\@evenhead}{\@oddhead}
\renewcommand{\@evenfoot}{\@oddfoot}
} % end \renewcommand (NB: this is missing from KoD95 p. 307)

\providecommand{\@evenhead}{}\providecommand{\@oddhead}{}
\providecommand{\@evenfoot}{}\providecommand{\@oddfoot}{}

\pagestyle{headings}

\renewcommand{\opening}[1]{\thispagestyle{firstpage}%
\hfill\parbox[t]{\rightfield}{\@date}

#1 \par\nobreak
} % end \renewcommand




