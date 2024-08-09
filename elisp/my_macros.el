; $Id$

; Purpose: Custom macros always loaded at Emacs startup time

; Usage: 
; (require 'my_macros) ; Custom macros

; These macros are created as follows:
; 1. Record the macro
; 2. Name the macro using M-x name-last-kbd-macro
; 3. Insert the macro in this file using M-x insert-kbd-macro 
; 4. Edit the keyboard macro edit-kbd-macro C-x C-k. 
;    This pops-up a buffer with a nicely commented macro format

; This procedure converts fixed-format *.F files to free-format *.F90 files
; Do this procedure in order

; Replace "     $" by "      " and place " &" at end of previous line
; query-replace-regexp "[\n]     \$" "&\n     " nil
(fset 'dlr2f90
   [escape ?x ?s ?e ?a ?r ?c ?h ?  ?f ?  ?  ?r ?  return ?^ ?  ?  ?  ?  ?  ?\\ ?$ return ?\C-b ?\C-d ?\C-p ?\C-e ?  ?& ?\C-n ?\C-a])

; Replace "! ... &" by "& ! ..."
; (query-replace-regexp "! \(.*\) &" "& ! \1" nil)

; Replace floating point number with floating point number f "1.23" by "1.23f"
; (query-replace-regexp "\([+-]?[0-9]*\.[0-9]*[LlDd]?\)" "\1f" nil)

; Replace "c    Comment" by "! Comment"
; query-replace-regexp "^c     " "! " nil
(fset 'cmt2f90
   [escape ?x ?s ?e ?a ?r ?c ?h ?  ?f ?  ?  ?r ?  return ?^ ?c ?  ?  ?  ?  ?  return ?\C-a ?\C-d ?\C-d ?\C-d ?\C-d ?\C-d ?! ?\C-n ?\C-a])

; Replace static_cast<prc_cmp>(foo) with PRC_CMP(foo)

; Replace "dbg.com" by "use dbg_mdl"
(fset 'dbg2f90
   [?\C-s ?d ?b ?g ?. ?c ?o ?m ?\C-a ?\C-k ?\C-k escape ?x ?s ?e ?a ?r ?c ?h ?  ?b ?a ?  return ?i ?m ?p ?l ?i ?c ?i ?t ?  ?n ?o ?n ?e return ?\C-a return ?\C-p tab ?u ?s ?e ?  ?d ?b ?g ?_ ?m ?d ?l ?  ?! ?  ?[ ?m ?d ?l ?] ?  ?D ?e ?b ?u ?g ?g ?i ?n ?g ?  ?c ?o ?n ?s ?t ?a ?n ?t ?s ?, ?  ?p ?r ?g ?_ ?n ?m ?, ?  ?d ?b ?g ?_ ?l ?v ?l ?\C-a])

; Replace "CCM.com" by "use ccm_grd"
(fset 'ccm2f90
   [?\C-s ?C ?C ?M ?. ?c ?o ?m ?\C-a ?\C-k ?\C-k escape ?x ?s ?e ?a ?r ?c ?h ?  ?b ?a ?  return ?i ?m ?p ?l ?i ?c ?i ?t ?  ?n ?o ?n ?e return ?\C-a return ?\C-p tab ?u ?s ?e ?  ?c ?c ?m ?_ ?g ?r ?d ?  ?! ?  ?[ ?m ?d ?l ?] ?  ?C ?C ?M ?  ?g ?r ?i ?d ?  ?p ?a ?r ?a ?m ?e ?t ?e ?r ?s ?\C-a])

; Replace "xtr.com" by "use xtr_mdl"
(fset 'xtr2f90
   [?\C-s ?x ?t ?r ?. ?c ?o ?m ?\C-a ?\C-k ?\C-k escape ?x ?s ?e ?a ?r ?c ?h ?  ?b ?a ?  return ?i ?m ?p ?l ?i ?c ?i ?t ?  ?n ?o ?n ?e return ?\C-a return ?\C-p tab ?u ?s ?e ?  ?x ?t ?r ?_ ?m ?d ?l ?, ?o ?n ?l ?y ?: ?x ?t ?r ?_ ?i ?n ?i ?  ?! ?  ?[ ?m ?d ?l ?] ?  ?E ?x ?t ?r ?a ?p ?o ?l ?a ?t ?i ?o ?n ?  ?c ?o ?n ?s ?t ?a ?n ?t ?s ?\C-a])

; Place lisp commands in buffer, execute using M-x eval-buffer
; Replace "foo :: bar" by "foo::bar"
; (query-replace-regexp "\(\w*)\)\W*::" "\1::" nil)
; Replace "character foo*80" by "character(80)::foo"
; (query-replace-regexp "character \(.*\)\*\([0-9]*\) " "character(\2)::\1 " nil)
; Replace "!foo" by "! foo"
; (query-replace-regexp "!\([^ ]\)" "! \1" nil)
; Replace "  &" by " &"
; (query-replace-regexp " *\&" " &" nil)
; Replace "end do     !" by "end do !"
; (query-replace-regexp "end \(do\|if\) *!" "end \1 !" nil)
; Replace "foo   !" by "foo !" as long as "!" is not followed by "="
; (query-replace-regexp "\([^ !\t\n]\) + +! " "\1 ! " nil)
; Replace non-advancing I/O usage
; (tags-query-replace ",$)'" ")',advance=\"no\"" nil)
; Replace logical comparisons
; (tags-query-replace "\\.gt\\." " > " nil)
; (tags-query-replace "\\.ge\\." " >= " nil)
; (tags-query-replace "\\.lt\\." " < " nil)
; (tags-query-replace "\\.le\\." " >= " nil)
; (tags-query-replace "\\.eq\\." " == " nil)
; (tags-query-replace "\\.ne\\." " /= " nil)

; Convert automatic arrays to dynamic arrays
(fset 'f772f90
   [?\C-a ?\C-k ?\C-k ?\C-y ?\C-y ?\C-y ?\C-p ?\C-p ?\C-p escape ?f ?\C-b ?\C-f ?, ?d ?i ?m ?e ?n ?s ?i ?o ?n ?( ?: ?) ?, ?a ?l ?l ?o ?c ?a ?t ?a ?b ?l ?e ?: ?: ?\C-d ?\C-s ?( ?\C-b ?\C-  ?\C-s ?) ?\C-b ?\C-f ?\C-w ?\C-n ?\C-a escape ?f escape ?b ?a ?l ?l ?o ?c ?a ?t ?e ?( ?\C-  ?\C-s ?  ?\C-f ?\C-b ?\C-w ?\C-s ?) ?\C-b ?\C-f ?, ?s ?t ?a ?t ?= ?r ?c ?d ?) ?\C-n ?\C-a escape ?f escape ?b ?\C-  ?\C-s ?  ?\C-b ?\C-f ?\C-w ?i ?f ?  ?( ?a ?l ?l ?o ?c ?a ?t ?e ?d ?( ?\C-s ?( ?\C-b ?) ?) ?  ?d ?e ?a ?l ?l ?o ?c ?a ?t ?e ?( ?, ?s ?t ?a ?t ?= ?r ?c ?d ?) ?\C-  ?\C-s ?) ?\C-b ?\C-f ?\C-w ?\C-a ?\C-s ?( ?\C-s ?\C-b ?\C-f ?\C-  ?\C-s ?) ?\C-b ?\C-w ?\C-y ?\C-s ?c ?a ?t ?e ?( ?, ?\C-b ?\C-y ?\C-n ?\C-a])

; Next three commands attempt to define ctl-s to search forward

;(global-set-key "\C-s" 'isearch-forward)
;(define-key global-map "\C-s" 'isearch-forward)
(fset 'ctls
   "xgloba s k is for ")
(global-set-key "\C-s"  'isearch-forward)
(global-set-key "\C-r"  'isearch-backward)
(global-set-key "\M-q"     'query-replace)
(global-set-key "\M-Q"     'tags-query-replace)

; Turn on Auto Fill mode automatically in Text mode and related modes
(add-hook 'text-mode-hook 'turn-on-auto-fill)

; Set the esc-q key to override the c-mode esc-q (which is fill paragraph)
(add-hook 'c-mode-hook
          (function
           (lambda ()
             (define-key esc-map "q" 'query-replace))))

; Define ctl-z map
(defvar ctl-z-map (make-keymap)
  "Default keymap for C-z commands.
The normal global definition of the character C-z indirects to this keymap.")
(global-set-key "\C-z" ctl-z-map)

; Insert frequently used pathnames at cursor in the command line
; Typing, e.g., C-zccr (C-z C-c C-c C-r) in the command line editor
; will insert /home/zender/crm at cursor
(fset 'fsh
   "/fs/cgd/home0/zender/")
(define-key ctl-z-map "\C-f\C-s\C-h" 'fsh)
(fset 'fsd
   "/fs/cgd/data0/zender/")
(define-key ctl-z-map "\C-f\C-s\C-d" 'fsd)
(fset 'sds
   "/fs/cgd/home0/zender/dst/dst")
(define-key ctl-z-map "\C-s\C-d\C-s" 'sds)
(fset 'sdm
   "/fs/cgd/home0/zender/dmr/dmr")
(define-key ctl-z-map "\C-s\C-d\C-m" 'sdm)
(fset 'ccr
   "/home/zender/crm/")
(define-key ctl-z-map "\C-c\C-c\C-r" 'ccr)
(fset 'cte
   "/home/zender/tex/")
(define-key ctl-z-map "\C-c\C-t\C-e" 'cte)
(fset 'cmi
   "/home/zender/mie/")
(define-key ctl-z-map "\C-c\C-m\C-i" 'cmi)
(fset 'cds
   "/home/zender/dst/")
(define-key ctl-z-map "\C-c\C-d\C-s" 'cds)
(fset 'cdm
   "/home/zender/dmr/")
(define-key ctl-z-map "\C-c\C-d\C-m" 'cdm)
(fset 'dtm
   "/data/zender/tmp/")
(define-key ctl-z-map "\C-d\C-t\C-m" 'dtm)
(fset 'tze
   "/data/zender/_aux0_/")
(define-key ctl-z-map "\C-t\C-z\C-e" 'tze)

; Insert anonymous mail header lines
(fset 'anon
   "anon@anon.penet.fiX-Anon-To: X-Anon-Password: ")
(define-key ctl-z-map "a" 'anon)
(setq c-mode-hook
      (function
       (lambda ()
         (define-key c-mode-map "\C-ca" 'anon)
         )))

; Insert a C comment at cursor
(fset 'c_comment
   "/*  */\M-3\C-b")
(define-key ctl-z-map "c" 'c_comment)
(setq c-mode-hook
      (function
       (lambda ()
         (define-key c-mode-map "\C-cc"    "/*  */\M-3\C-b")
         )))

; Insert an Asciidoc keyword emphasis at cursor
(fset 'c_keyword
   "**  **\M-3\C-b")
(define-key ctl-z-map "k" 'c_keyword)
(setq c-mode-hook
      (function
       (lambda ()
         (define-key c-mode-map "\C-cc"    "**  **\M-3\C-b")
         )))

; Insert an editing begin comment at cursor 
(fset 'f_edit_begin
   "\C-ac++csz\C-a")
(define-key ctl-z-map "5" 'f_edit_begin)
(setq c-mode-hook
      (function
       (lambda ()
         (define-key c-mode-map "\C-c5"    "\C-ac++csz\C-a")
         )))

; Insert an editing end comment at cursor 
(fset 'f_edit_end
   "\C-ac--csz\C-a")
(define-key ctl-z-map "6" 'f_edit_end)
(setq c-mode-hook
      (function
       (lambda ()
         (define-key c-mode-map "\C-c6"    "\C-ac--csz\C-a")
         )))

; Change $     foo, ! Comment into integer foo    ! Comment
(fset 'f_int
   [?\C-s ?$ ?\C-b ?\C-d ?  ?i ?n ?t ?e ?g ?e ?r ?\C-d ?\C-d ?\C-d ?\C-d ?\C-s ?, ?\C-b ?\C-d tab ?\C-n ?\C-a])

; Insert a starting fortran ruler at cursor
(fset 'f_ruler_start
   "\C-ac++csz789012345678901234567890123456789012345678901234567890123456789012\C-a")
(define-key ctl-z-map "7" 'f_ruler_start)
(setq c-mode-hook
      (function
       (lambda ()
         (define-key c-mode-map "\C-c7"    "\C-ac++csz789012345678901234567890123456789012345678901234567890123456789012\C-a")
         )))
; Insert a \verb'' at cursor
(defalias 'latex_verb
  (read-kbd-macro "\\ verb'' C-b"))
(define-key ctl-z-map "v" 'latex_verb)
; Insert a \textcolor{blue}{} at cursor
(defalias 'latex_blue
  (read-kbd-macro "\\ textcolor{blue}{} C-b"))
(define-key ctl-z-map "b" 'latex_blue)
; Insert a \csznote{} % end csznote at cursor
(defalias 'latex_note
  (read-kbd-macro "\\ csznote{ } % end csznote"))
(define-key ctl-z-map "n" 'latex_note)
; Insert a \textcolor{red}{} at cursor
(defalias 'latex_red
  (read-kbd-macro "\\ textcolor{red}{} C-b"))
(define-key ctl-z-map "r" 'latex_red)
; Insert a \textcolor{green}{} at cursor
(defalias 'latex_grn
  (read-kbd-macro "\\ textcolor{green}{} C-b"))
(define-key ctl-z-map "g" 'latex_grn)
; Insert a LaTeX equation (\ref{eqn:}) at cursor
(defalias 'latex_eqn
  (read-kbd-macro "(\\ ref{eqn:}) C-b C-b"))
(define-key ctl-z-map "e" 'latex_eqn)
; Insert a LaTeX figure (Figure~\ref{fgr:}) at cursor
(defalias 'latex_fgr
  (read-kbd-macro "Figure~\\ ref{fgr:} C-b"))
(define-key ctl-z-map "f" 'latex_fgr)
; Insert an ending fortran ruler at cursor
(fset 'f_ruler_end
   "\C-ac--csz789012345678901234567890123456789012345678901234567890123456789012\C-a")
(define-key ctl-z-map "8" 'f_ruler_end)
(setq c-mode-hook
      (function
       (lambda ()
         (define-key c-mode-map "\C-c8"    "\C-ac--csz789012345678901234567890123456789012345678901234567890123456789012\C-a")
         )))
; Insert an if(dbg_lvl == ){;}// end dbg
(fset 'c_dbg
   "	if(dbg_lvl == ){	;	} // end if dbg")
(define-key ctl-z-map "d" 'c_dbg)
(setq c-mode-hook
      (function
       (lambda ()
         (define-key c-mode-map "\C-cd" 'c_dbg)
         )))
; Insert an if(){;}else{;}// end else
; (fset 'c_insert_else
;    "	if(){	;	}else{	;	} // end if-else */	")
; (define-key ctl-z-map "e"  'c_insert_else)
; (setq c-mode-hook
;       (function
;        (lambda ()
;          (define-key c-mode-map "\C-ce"  'c_insert_else)
;          )))
; Insert a foo loop
; (fset 'c_foo_loop
;    "	for(idx=0;idx<foo_nbr;idx++){	;	} // end loop over foo	")
; (define-key ctl-z-map "f" 'c_foo_loop)
; (setq c-mode-hook
;       (function
;        (lambda ()
;          (define-key c-mode-map "\C-cf" 'c_foo_loop)
;          )))
; Insert an if(){;}// end if
(fset 'c_insert_if
   "	if(){	;} // end if		")
(define-key ctl-z-map "i" 'c_insert_if)
(setq c-mode-hook
      (function
       (lambda ()
         (define-key c-mode-map "\C-ci" 'c_insert_if)
         )))
; C-comment out from point until the end of the current line
(fset 'c_comment_out
   "/**/")
(define-key ctl-z-map "o" 'c_comment_out)
(setq c-mode-hook
      (function
       (lambda ()
         (define-key c-mode-map "\C-co" 'c_comment_out)
         )))
; Insert a print statement
(fset 'c_print
   "	(void)fprintf(stderr,\"\\n\");\\n")
(define-key ctl-z-map "p" 'c_print)
(setq c-mode-hook
      (function
       (lambda ()
         (define-key c-mode-map "\C-cp" 'c_print)
         )))
; Insert a size loop
(fset 'c_sz_loop
   "	for(sz_idx=0;sz_idx<sz_nbr;sz_idx++){	;	} // end loop over sz	")
(define-key ctl-z-map "s" 'c_sz_loop)
(setq c-mode-hook
      (function
       (lambda ()
         (define-key c-mode-map "\C-cs" 'c_sz_loop)
         )))
; Tab then move down a line, for getting correct indents in C mode
(fset 'c_tab
   [tab 14 1])
(define-key ctl-z-map "t" 'c_tab)
(setq c-mode-hook
      (function
       (lambda ()
         (define-key c-mode-map "\C-ct" 'c_tab)
         )))
; Undo the next pair of C comment indicators
(fset 'c_comment_undo
   "/**/")
(define-key ctl-z-map "u" 'c_comment_undo)
(setq c-mode-hook
      (function
       (lambda ()
         (define-key c-mode-map "\C-cu" 'c_comment_undo)
         )))

; Lowercase a file the hard way
(fset 'lowercase_file
   "rAa<rBb<rCc<rDd<rEe<rFf<rGg<rHh<rIi<rJjj<rKk<rLl<rMm<rNn<rOo<rPp<rQqrRr<rSs<rTt<rUu<rVv<rWw<rXx<rYy<rZz<")

; Lowercase an HTML file
(fset 'lowercase_html_file
   "^[r<DD><dd>[^")

(provide 'my_macros)

; end of my_macros.el


