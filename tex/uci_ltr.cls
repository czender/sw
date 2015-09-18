% $Id$

% Purpose: Letter class for LaTeX2e with UCI letterhead

% Based on mpletter.cls from KoD95 p. 300

% Usage: 
% /bin/cp ${HOME}/tex/uci_ltr.cls ${TEXMFDIR}/tex/latex2e/csz

% NB: To add a new class or style file to the TeX tree, install in csz
% directory and then update ${TEXMFDIR}/ls-R by running `texhash' as root 

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{uci_ltr}
\RequirePackage{ifthen}
\newboolean{@spanish}
\setboolean{@spanish}{false}
\DeclareOption{spanish}{\setboolean{@spanish}{true}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{letter}}
\ProcessOptions
\LoadClass[a4paper]{letter}

\newcommand{\englishnames}{%
\newcommand{\yrefname}{\textsl{Your Ref.}}
\newcommand{\ymailname}{\textsl{Your Letter from}}
\newcommand{\myrefname}{\textsl{Our Ref.}}
\newcommand{\subjectname}{\textsl{Subject}}
\newcommand{\telephonename}{Telephone}
\newcommand{\stationname}{Mesa Lab}
\newcommand{\USAname}{USA}
\newcommand{\areacode}{949}
\newcommand{\postcode}{92612-3100}
} % end \newcommand

\newcommand{\spanishnames}{%
\newcommand{\yrefname}{\textsl{Su Ref.}}
\newcommand{\ymailname}{\textsl{Su Letro de}}
\newcommand{\myrefname}{\textsl{Nuestro Ref.}}
\newcommand{\subjectname}{\textsl{Subjecto}}
\newcommand{\telephonename}{Telephono}
\newcommand{\stationname}{Mesa Lab}
\newcommand{\USAname}{Estados Unidos}
\newcommand{\areacode}{949}
\newcommand{\postcode}{92612-3100}
} % end \newcommand

\ifthenelse{\boolean{@spanish}}
{\RequirePackage{spanish}\spanishnames}{\englishnames}

\newcommand{\@yref}{}
\newcommand{\@myref}{}
\newcommand{\@fax}{}
\newcommand{\@email}{}
\newcommand{\@ymail}{}
\newcommand{\@subject}{}

\newcommand{\yref}[1]{\renewcommand{\@yref}{\yrefname: #1}}
\newcommand{\ymail}[1]{\renewcommand{\@ymail}{\ymailname: #1}}
\newcommand{\myref}[1]{\renewcommand{\@myref}{\myrefname: #1}}
\newcommand{\subject}[1]{\renewcommand{\@subject}{\subjectname: #1}}
\newcommand{\fax}[1]{\renewcommand{\@fax}{#1}}
\newcommand{\email}[1]{\renewcommand{\@email}{#1}}

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
\newfont{\poetfont}{cmdunh10 scaled \magstep4}

% Sum of \leftfield and \rightfield is \textwidth
\newlength{\wholefield} \setlength{\wholefield}{\textwidth} 
\newlength{\leftfield} \setlength{\leftfield}{4.75in} 
\newlength{\rightfield} \setlength{\rightfield}{1.75in}

\newsavebox{\FIRM} \newsavebox{\firmaddress}
\newsavebox{\firm} \newsavebox{\firmreturn}

\sbox{\FIRM}
{\parbox[t]{\wholefield}{
\begin{center}
{\xiirn UNIVERSITY\, OF\, CALIFORNIA\, AT\, IRVINE}\\
{\xiissn School of Physical Sciences}\\
{\xssn Departments of Earth System Science and Computer Science}\\
{\viiisssl ESS Department $\bullet$ 3228 Croul Hall $\bullet$ Irvine, California 92697-3100
$\bullet$ USA}
\end{center}
}}

\sbox{\firm}
{\parbox[t]{\leftfield}{\xiissn UNIVERSITY\, OF\, CALIFORNIA\, AT\, IRVINE}}  

\sbox{\firmreturn}
{\viiisssl\underline{ESS Department, 3228 Croul Hall, Irvine CA~~92697-3100, \USAname}}

\sbox{\firmaddress}
{\parbox[t]{\rightfield}{\viiisssl\baselineskip10pt
Department of Earth System Science\\3228 Croul Hall\\Irvine, CA~~92697-3100\\\USAname}}

\newsavebox{\firmhead} \newsavebox{\firmfoot}

% This header has only the NCAR banner
%\sbox{\firmhead}
%{\parbox{\textwidth}{\usebox{\FIRM}\\[3pt] \rule{\textwidth}{1pt}}} 

% This header places the short NCAR address next to the NCAR banner
\sbox{\firmhead}
{\parbox{\textwidth}{\usebox{\FIRM}
%\raisebox{3pt}{\usebox{\firmaddress}}
\\[3pt] \rule{\textwidth}{1pt}}} 

\sbox{\firmfoot}
{\parbox{\textwidth}{\rule{\textwidth}{0.6pt}\\
\viitin\setlength{\baselineskip}{8pt}%
\begin{center}
The National Center for Atmospheric Research is operated by the
University Corporation \\
for Atmospheric Research under sponsorship of the National Science 
Foundation.
%\vskip.2em
%An Equal Opportunity/Affirmative Action Employer \\
\end{center}
}% end \parbox
}% end \sbox

\renewcommand{\ps@firstpage}
% headsep sets the amount of vertical space between the NCAR
% letterhead and the main body of the page, which, in this case, is
% the to/from address block.
%{\setlength{\headheight}{50pt}\setlength{\headsep}{-35pt}%
%{\setlength{\headheight}{70pt}\setlength{\headsep}{25pt}%
{\setlength{\headheight}{25pt}\setlength{\headsep}{50pt}%
\renewcommand{\@oddhead}{\usebox{\firmhead}}%
%\renewcommand{\@oddfoot}{\raisebox{-20pt}[0pt]{\usebox{\firmfoot}}}%
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
%\raisebox{0in}{%
\parbox[t]{\textwidth}{
%\usebox{\firmreturn}\\
\parbox[t][3.5cm][c]{\leftfield}{\toname\\\toaddress}%
%\raisebox{-72pt}{%
\parbox[t]{\rightfield}{\fromname
\ifthenelse{\equal{\telephonenum}{}}
{}{\\ Tel.: (\areacode)\,\telephonenum}
\ifthenelse{\equal{\@fax}{}}
{}{\\ FAX: (\areacode)\,\@fax}
\ifthenelse{\equal{\@email}{}}
{}{\\{\xssn E-mail: \@email}}
\\[5mm] \@date
}% end \parbox
%}% end \raisebox
}% end \parbox
%}% end \raisebox
\par
\rule{\textwidth}{0.6pt}
\makebox[\leftfield][l]
{\ifthenelse{\equal{\@yref}{}}
{\@ymail}{\@yref\hfill\@ymail\hfill}
}% end \makebox
\@myref\par
\ifthenelse{\equal{\@subject}{}}
{}{\@subject\par}
%\vspace{\parskip} #1 \par\nobreak
#1 \par\nobreak
} % end \renewcommand




