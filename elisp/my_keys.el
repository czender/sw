; $Id$

; Purpose: Emacs key mapping

; Usage: 
; (require 'my_keys) ; Custom key maps

(global-unset-key "\e\e")

(define-key esc-map "." 'end-of-buffer)
(define-key esc-map "," 'beginning-of-buffer)
(define-key esc-map "og" 'goto-line)
(define-key esc-map "ol" 'what-line)

; Add these two lines to avoid ctrl-s use with vt100
(define-key esc-map "g" 'goto-line)
(define-key esc-map "q" 'query-replace)
(define-key esc-map "r" 'replace-string)
(define-key esc-map "p" 'fill-paragraph)
(define-key esc-map "j" 'fill-paragraph)
(define-key esc-map "s" 'isearch-forward)
(define-key ctl-x-map "" 'save-buffer)

; Keypad definition
(define-key esc-map "OQ" 'describe-key-briefly)             ; PF2
(define-key esc-map "OR" 'isearch-forward)                  ; PF3
(define-key esc-map "OS" 'kill-line)                        ; PF4
(define-key esc-map "Op" 'next-line)                        ; 0
(define-key esc-map "Oq" 'forward-word)                     ; 1
(define-key esc-map "Or" 'end-of-line)                      ; 2
(define-key esc-map "Os" 'forward-char)                     ; 3
(define-key esc-map "Ot" 'end-of-buffer)                    ; 4
(define-key esc-map "Ou" 'beginning-of-buffer)              ; 5
(define-key esc-map "Ov" 'kill-region)                      ; 6
(define-key esc-map "Ow" 'execute-extended-command)         ; 7
(define-key esc-map "Ox" 'scroll-up)                        ; 8
(define-key esc-map "Oy" 'query-replace)                    ; 9
(define-key esc-map "Om" 'kill-word)                        ; _
(define-key esc-map "Ol" 'delete-char)                      ; ,
(define-key esc-map "On" 'set-mark-command)                 ; .

; VT100 arrow keys
(define-key esc-map "OB" 'next-line)                  ;DOWN arrow key
(define-key esc-map "OA" 'previous-line)              ;UP arrow key
(define-key esc-map "OD" 'backward-char)              ;LEFT arrow key
(define-key esc-map "OC" 'forward-char)               ;RIGHT arrow key

; SUN arrow keys
; NB: M-[ is defined as nil (it used to be backwards-paragraph)
(define-key esc-map "[" 'nil)
(define-key esc-map "[B" 'next-line)                  ;DOWN arrow key
(define-key esc-map "[A" 'previous-line)              ;UP arrow key
(define-key esc-map "[D" 'backward-char)              ;LEFT arrow key
(define-key esc-map "[C" 'forward-char)               ;RIGHT arrow key
(define-key esc-map "[2" 'execute-extended-command)

(provide 'my_keys)

; end of my_keys.el
