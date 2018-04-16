(TeX-add-style-hook "bemlib"
 (lambda ()
    (TeX-add-symbols
     "J")
    (TeX-run-style-hooks
     "pdftex"
     "hyperref"
     "ps2pdf"
     "times"
     "doxygen"
     "alltt"
     "float"
     "graphicx"
     "fancyhdr"
     "makeidx"
     "a4wide"
     "latex2e"
     "bk10"
     "book"
     "a4paper")))

