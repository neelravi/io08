\ProvidesClass{TREX}
\LoadClass[a4paper,twoside,final,12pt]{article}

\usepackage[paper=a4paper,
top=1.18in,bottom=1in,right=0.75in,left=0.75in]{geometry}% http://ctan.org/pkg/geometry
\usepackage[utf8]{inputenc}
\usepackage{lastpage}
\usepackage{graphicx}
\usepackage{fancyhdr}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{titlesec}
\usepackage{colortbl}

\renewcommand{\familydefault}{\sfdefault}

\usepackage[printonlyused,nohyperlinks]{acronym}


\usepackage[pdftex, pdfusetitle, hyperfigures, hyperindex,
bookmarksnumbered,colorlinks, bookmarks, breaklinks,
linktocpage, citecolor=blue, urlcolor=blue, linkcolor=black]{hyperref}

\usepackage[titles]{tocloft}
\renewcommand{\contentsname}{Table of Contents}
\renewcommand{\cftsecfont}{\normalfont}% titles in bold
\renewcommand{\cftsecpagefont}{\normalfont}% page numbers in bold
\renewcommand{\cftdotsep}{1}
\renewcommand{\cftsecleader}{\cftdotfill{\cftsecdotsep}}% dot leaders in bold

\renewcommand{\cftsecdotsep}{\cftdotsep}



\setlength{\headheight}{1cm}
\setlength{\footskip}{2cm}
\addtolength{\textheight}{-\footskip}


%%--Set the section number depth.
\setcounter{secnumdepth}{3}
\setcounter{tocdepth}{4}
\setcounter{page}{2}

\bibliographystyle{IEEEtran}
\sloppy

%%-----------------------------------------------
%%-- Float related environments and commands.
%%-----------------------------------------------
%%--Modify behavior of the figure and table environment.
\def\fps@figure{hbt}
\def\fnum@figure{\textbf{\figurename~\thefigure}}
\def\fps@table{htb}
\def\fnum@table{\textbf{\tablename~\thetable}}
\def\mathfigure{\thesisfiguretrue\figure}
\def\endmathfigure{\endfigure\thesisfigurefalse}

% --------------------------------------------
% Headers
% --------------------------------------------
\pagestyle{fancyplain}
\renewcommand{\sectionmark}[1]{\markright{\thesection.\ #1}}


\usepackage{xcolor}
\definecolor{myblue}{HTML}{2F5496}
\titleformat{\section}
{\color{myblue}\normalfont\Large}
{\color{myblue}\thesection}{1em}{}


\makeatletter
\makeindex
\makeglossary
\setcounter{tocdepth}{4}
\setcounter{secnumdepth}{1}

% Table of contents with section*

\usepackage{xparse}

\let\oldsection\section
\makeatletter
\newcounter{@secnumdepth}
\RenewDocumentCommand{\section}{s o m}{%
  \IfBooleanTF{#1}
    {\setcounter{@secnumdepth}{\value{secnumdepth}}% Store secnumdepth
     \setcounter{secnumdepth}{0}% Print only up to \chapter numbers
     \oldsection{#3}% \section*
     \setcounter{secnumdepth}{\value{@secnumdepth}}}% Restore secnumdepth
    {\IfValueTF{#2}% \section
       {\oldsection[#2]{#3}}% \section[.]{..}
       {\oldsection{#3}}}% \section{..}
}
\makeatother

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Floating figures
\renewcommand\floatpagefraction{.9}
\renewcommand\topfraction{.9}
\renewcommand\bottomfraction{.9}
\renewcommand\textfraction{.1}
\setcounter{totalnumber}{50}
\setcounter{topnumber}{50}
\setcounter{bottomnumber}{50}

\endinput


