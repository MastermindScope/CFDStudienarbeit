@echo off
Rem  htcontext.bat (2020-10-25-17:53), generated from tex4ht-mkht.tex
Rem  Copyright 2009-2020 TeX Users Group
Rem  Copyright 2003-2009 Eitan M. Gurari
Rem 
Rem  This work may be distributed and/or modified under the
Rem  conditions of the LaTeX Project Public License, either
Rem  version 1.3 of this license or (at your option) any
Rem  later version. The latest version of this license is in
Rem    http://www.latex-project.org/lppl.txt
Rem  and version 1.3 or later is part of all distributions
Rem  of LaTeX version 2003/12/01 or later.
Rem 
Rem  This work has the LPPL maintenance status "maintained".
Rem 
Rem  The Current Maintainer of this work
Rem  is the TeX4ht Project <https://tug.org/tex4ht>.
Rem 
Rem  If you modify this file, changing the
Rem  version identification be appreciated.


        call texexec --arg="ht-1=%2" --use=tex4ht --dvi --nobackend %5 %1 
        tex4ht %1 -i/tex4ht/ht-fonts/%3 -ewin32/tex4ht.env
        t4ht %1 %4 -ewin32/tex4ht.env 



