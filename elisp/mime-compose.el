;;; --------------------------------------------------------------------------
;;; File: --- mime-compose.el ---
;;; Author: Marc Andreessen (marca@ncsa.uiuc.edu)
;;; Additional code: Keith Waclena (k-waclena@uchicago.edu).
;;; Ronald Florence <ron@mlfarm.com>, 21 Aug 1993:  patches for emacs-19 from
;;;   wmperry@mango.ucs.indiana.edu, changed anon-ftp inclusions to
;;;   query transfer-type, added richtext parsing code from
;;;   gnus-mime.el.
;;; ron@auda.mlfarm.com, 31 Oct 1994: added application types
;;; ron@auda.mlfarm.com, 25 Jul 1996: fixed on-the-fly audio recording
;;;   for Suns.
;;; Ralf Fassel <ralf@natlab.research.philips.com>, 25 Jul 1996:
;;;    patched code for included raw-binaries.
;;;
;;; Copyright (C) National Center for Supercomputing Applications, 1992.
;;;
;;; This program is free software; you can redistribute it and/or modify
;;; it under the terms of the GNU General Public License as published by
;;; the Free Software Foundation; either version 1, or (at your option)
;;; any later version.
;;;
;;; This program is distributed in the hope that it will be useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;;; GNU General Public License for more details.
;;;
;;; You should have received a copy of the GNU General Public License
;;; along with your copy of Emacs; if not, write to the Free Software
;;; Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
;;;
;;; -------------------------------- CONTENTS --------------------------------
;;;
;;; mime-compose: Utility routines for composing MIME-compliant mail.
;;; $Revision$
;;; $Date$
;;;
;;; Canonical list of features:
;;;   Automatic MIME header construction.
;;;   Include GIF/JPEG image.
;;;   Include audio file.
;;;   Include PostScript file.
;;;   Include raw binary/nonbinary file.
;;;   Include xwd window dump.
;;;   Include reference to anonymous/regular FTP.
;;;   Include audio snippet recorded on the fly.
;;;   Convert region to MIME richtext.
;;;   Convert region to any ISO 8859 charset.
;;;   Optional conversion of plaintext bodyparts to quoted-printable
;;;     with arbitrary charset when messages are sent.
;;;   Deemphasizing/highlighting of MIME headers.
;;;   Completion on content type and charset.
;;;   Automatic encoding in base64 and quoted-printable formats.
;;;   Selective display hides raw data.
;;;   Works with mail-mode and mh.
;;;
;;; ------------------------------ INSTRUCTIONS ------------------------------
;;;
;;; Use the normal Emacs mail composer (C-x m).
;;;
;;; (Or, use with Emacs mh-e by loading this file *after* loading mh-e.
;;; Try putting (require 'mime-compose) in mh-letter-mode-hook.)
;;;
;;; Do nothing special to prepare a message to have MIME elements
;;; included in it.
;;;
;;; The basic commands to add MIME elements (images, audio, etc.) to a
;;; message are as follows:
;;;
;;; mail-mode  (mh-e)       function                 what happens
;;; ~~~~~~~~~  (~~~~~~~~~)  ~~~~~~~~                 ~~~~~~~~~~~~
;;; C-c g      (C-c C-m g)  mime-include-gif         Add a GIF file.
;;; C-c j      (C-c C-m j)  mime-include-jpeg        Add a JPEG file.
;;; C-c a      (C-c C-m a)  mime-include-audio       Add an audio file.
;;; C-c p      (C-c C-m p)  mime-include-postscript  Add a PostScript file.
;;;
;;; (Note that mime-compose assumes you have the 'mmencode' program
;;; installed on your system.  See 'WHAT MIME IS' below for more
;;; information on mmencode and the metamail distribution.)
;;;
;;; Some mime-compose commands create data themselves; these follow:
;;;
;;; C-c x      (C-c C-m x)
;;;   mime-include-xwd-dump
;;;   Add the result of an X-window dump.  The program named in
;;;   mime-xwd-command will be run, and the resulting dump will be
;;;   inserted into the message.
;;; C-c s      (C-c C-m s)
;;;   mime-include-audio-snippet
;;;   Add an audio snippet, recorded on the fly.  CURRENTLY THIS WORKS
;;;   ONLY FOR SILICON GRAPHICS INDIGO AND 4D/35's.  Recording begins
;;;   immediately; press 'y' to end recording or 'n' to abort the
;;;   whole process.  The resulting audio file will be converted to
;;;   standard mulaw format and incorporated into the message.
;;;
;;; If you have a raw binary file and MIME or mime-compose doesn't
;;; have built-in support for its format (e.g. an Emacs Lisp
;;; byte-compiled file), you can use:
;;;
;;; C-c r      (C-c C-m r)
;;;   mime-include-raw-binary
;;;   Add a raw binary file.  You will be prompted for both the
;;;   filename and the content type of the file; if you do not give a
;;;   content type, the default (application/octet-stream) will be
;;;   used, and the recipient will be able to have his/her MIME mail
;;;   handler extract the raw binary file from the message.
;;;
;;; Similarly, to include nonbinary (text) files using
;;; quoted-printable encoding, use:
;;;
;;; C-c n       (C-c C-m n)
;;;   mime-include-raw-nonbinary 
;;;   Add a raw nonbinary (text) file.  You will be prompted for both
;;;   the filename and the content type of the file (which defaults to
;;;   text/plain).  With prefix arg, you will also be prompted for the
;;;   character set (default is US-ASCII).
;;;
;;; In addition to including files and generating inclusions on the
;;; fly, you can also point to external elements: files that will not
;;; be included in the document, but can be accessed by the recipient
;;; in some other way (most commonly, via FTP).  The following
;;; commands handle this:
;;;
;;; C-c e       (C-c C-m e)
;;;   mime-include-external-anonftp
;;;   Point to an external file (assumed to be accessable via
;;;   anonymous FTP).  You will be prompted for the name of the FTP
;;;   site, the remote directory name, and remote filename, the remote
;;;   file's content type, and a description of the remote file.
;;; C-c f       (C-c C-m f)
;;;   mime-include-external-ftp
;;;   This is the same as 'C-c e', except that the file will be
;;;   accessed via regular FTP rather than anonymous FTP -- a username
;;;   and password will have to be provided by the recipient to gain
;;;   access to the file.
;;;
;;; Note that whenever you are prompted for a content type, Emacs'
;;; completion feature is active: press TAB for a list of valid types.
;;; You can also enter a type not in the completion list.
;;;
;;; If you type in text that belongs in a character set other than the
;;; default (US-ASCII), you can use the following function to encode
;;; the text and generate appropriate MIME headers:
;;;
;;; C-c C-r i   (C-c C-m C-r i)
;;;   mime-region-to-charset 
;;;   Encode region in an alternate character set.  (MIME only
;;;   sanctions the use of ISO charsets; thus, the command key for
;;;   this function is 'i'.)  You will be prompted for a character set
;;;   (minibuffer completion is provided).
;;;
;;; MIME also defines a 'richtext' format; you can encode the current
;;; region as richtext with:
;;;
;;; C-c C-r r   (C-c C-m C-r r)
;;;   mime-region-to-richtext
;;;   Encode region as richtext.  With prefix arg, you will be
;;;   prompted for a character set, else the default (US-ASCII) is
;;;   used.
;;;
;;; If you regularly use 8-bit characters in your messages, you will
;;; probably want all of your plaintext bodyparts automatically
;;; encoded in quoted-printable and labeled as belonging to the
;;; character set that you're using when a message is sent.  To have
;;; this happen, set this variable:
;;;
;;; mime-encode-plaintext-on-send  (default NIL)
;;;   If T, all text/plain bodyparts in the message will be encoded in
;;;   quoted-printable and labeled with charset mime-default-charset
;;;   (by default, US-ASCII) when a message is sent.  If NIL,
;;;   text/plain bodyparts will not be touched.
;;;
;;; ---------------------------- ADDITIONAL NOTES ----------------------------
;;;
;;; mime-compose uses Emacs' selective-display feature: only the first
;;; line of any encoded data file will be displayed, followed by
;;; ellipses (indicating that some data is not being shown).  See the
;;; variable 'mime-use-selective-display' below.
;;;
;;; If you are running Lucid Emacs, the mail-mode popup menu (attached
;;; to the third mouse button) will include mime-compose entries.
;;;
;;; If you are running Lucid Emacs or Epoch, highlighting will be used
;;; to deemphasize the various MIME headers (but emphasize the various
;;; MIME content types).  You can turn this feature off; see the
;;; variable 'mime-use-highlighting'.
;;;
;;; After your message has been `mimified' (by including a MIME
;;; element), it is best not to put trailing text outside the final
;;; boundary at the end of the file -- such text will not be
;;; considered to be part of the message by MIME-compliant mail
;;; readers (although it will still be sent).
;;;
;;; As you compose a complex MIME message, you may notice useless
;;; bodyparts accumulating: extra text/plain bodyparts, in particular,
;;; containing no text.  These bodyparts will be stripped from the
;;; message before the message is sent, so you (and I) won't look like
;;; a moron to the recipient.
;;;
;;; A command that usually isn't necessary, but is provided in case
;;; you wish to send a plaintext message with the various MIME headers
;;; and boundaries, is:
;;;
;;; C-c m     (C-c C-m m)    mime-mimify-message   Mimify a message.
;;;
;;; MIME messages can contain elements and structures not yet
;;; supported by mime-compose.  If you have ideas or code for support
;;; that should be provided by mime-compose, please send them to the
;;; author.
;;;
;;; ------------------------ WHAT MIME-COMPOSE IS NOT ------------------------
;;;
;;; mime-compose is not a MIME message handler.  It will not interpret
;;; MIME messages, display images, or anything similar.
;;;
;;; mime-compose is not intelligent enough (yet) to construct complex
;;; MIME messages (with nested boundaries, parallel message elements,
;;; and so on).
;;;
;;; mime-compose will not enforce correctness (MIME compliance) on
;;; your messages.  mime-compose generates MIME-compliant message
;;; elements, but will sit quietly if you alter them or add your own
;;; incorrect elements.
;;;
;;; In particular, note that the MIME specification demands a blank
;;; line following the Content declarations for a bodypart.
;;; mime-compose will give you that blank line, but will not demand
;;; that you leave it blank; if you don't, your message will not be
;;; happy.
;;;
;;; ------------------------------ WHAT MIME IS ------------------------------
;;;
;;; MIME defines a format for email messages containing non-plaintext
;;; elements (images, audio, etc.).  MIME is detailed in Internet RFC
;;; 1341, by N. Borenstein and N. Freed.  You can FTP this RFC from
;;; many archive sites, including uxc.cso.uiuc.edu.
;;;
;;; Few mail readers handle MIME messages, yet.  However, most popular
;;; mail readers can be easily patched to feed MIME messages to a
;;; program called 'metamail', which can handle MIME messages.  You
;;; can FTP metamail from thumper.bellcore.com in /pub/nsb as
;;; mm.tar.Z.  Since mime-compose requires the existence of the
;;; program 'mmencode' (from the metamail distribution) to insert
;;; binary and nonbinary files into messages, it is a Good Idea to
;;; have metamail installed on your system.
;;;
;;; --------------------------------------------------------------------------
;;; LCD Archive Entry:
;;; mime-compose|Marc Andreessen|marca@ncsa.uiuc.edu|
;;; MIME-compliant message generation utilities.|
;;; 1992-11-21|1.47|~/misc/mime-compose.el.Z|
;;; --------------------------------------------------------------------------

(provide 'mime-compose)

(defvar mime-running-mh-e (boundp 'mh-letter-mode-map)
  "Non-nil if running under mh-e.")

(if (not mime-running-mh-e)
    (require 'sendmail))

;;; ---------------------- User-customizable variables -----------------------

(defvar mime-use-selective-display t
  "*Flag for using selective-display to hide bodies of MIME enclosures.
If non-NIL, selective-display will be used; if NIL, it will not be used.")

(defvar mime-default-charset "us-ascii"
  "*Default character set for MIME messages elements.  According to the
MIME specification, this can be either US-ASCII or iso-8859-x, where x
must be between 1 and 9 inclusive.")

(defvar mime-encode-plaintext-on-send nil
  "*Non-NIL if plaintext bodyparts should be encoded in quoted-printable
and labeled with mime-default-charset when a message is sent; NIL
otherwise.")

(defvar mime-use-highlighting t
  "*Flag to use highlighting for MIME headers and content types in
Epoch or Lucid Emacs; if non-NIL, highlighting will be used.")

(defvar mime-deemphasize-color "grey80"
  "*Color for de-highlighting MIME headers in Epoch or Lucid Emacs.")

(defvar mime-emphasize-color "yellow"
  "*Color for highlighting MIME content types in Epoch or Lucid Emacs.")

(defvar mime-name-included-files t
  "*If non-NIL, use name attribute for included files.")

(defvar mime-use-waiting-messages t
  "*If non-NIL, enable waiting messages feature.")

(defvar mime-primary-boundary "mysteryboxofun"
  "*Word used as the primary MIME boundary.")

(defvar mime-xwd-command "xwd -frame"
  "*Command used to do a window dump under the X Window System.")

(defvar mime-encode-base64-command "mmencode"
  "*Command used to encode data in base64 format.")

(defvar mime-encode-qp-command "mmencode -q"
  "*Command used to encode data in quoted-printable format.")

(defvar mime-babbling-description "talking"
  "*Adjective(s) applying to audio snippets.")

;;; ---------------------------- Other variables -----------------------------

(defvar mime-valid-include-types
  '(("image/gif" 1)
    ("image/jpeg" 2)
    ("image/x-xbm" 3)
    ("image/x-xwd" 4)
    ("application/postscript" 5)
    ("application/andrew-inset" 6)
    ("application/octet-stream" 7)
    ("application/gnuplot" 8)
    ("application/xfig" 9)
    ("application/oleo" 10)
    ("application/pdf" 11)
    ("application/html" 12)
    ("text/richtext" 13)
    ("text/plain" 14)
    ("audio/basic" 15)
    ("video/mpeg" 16)
    ("message/rfc822" 17))
  "A list of valid content types for minibuffer completion.")

(defvar mime-valid-charsets
  '(("us-ascii" 1)
    ("iso-8859-1" 2)
    ("iso-8859-2" 3)
    ("iso-8859-3" 4)
    ("iso-8859-4" 5)
    ("iso-8859-5" 6)
    ("iso-8859-6" 7)
    ("iso-8859-7" 8)
    ("iso-8859-8" 9)
    ("iso-8859-9" 10))
  "A list of valid charset names for minibuffer completion.")

(defvar mime-using-silicon-graphics (or (eq system-type 'silicon-graphics-unix)
					(eq system-type 'irix))
  "Flag to indicate use of Silicon Graphics platform.  If T, Emacs is being
run on a Silicon Graphics workstation; else it is not.")

(defvar mime-running-fsf19 
  (and (or (string= "19" (substring emacs-version 0 2))
           (string= "20" (substring emacs-version 0 2)))
       (not (string-match "Lucid" emacs-version))
       (eq window-system 'x))
  "Non-nil if running GNU Emacs 19 or 20 in xwindows")
 
(defvar mime-running-lemacs (string-match "Lucid" emacs-version)
  "Non-nil if running Lucid Emacs.")

(defvar mime-running-epoch (boundp 'epoch::version)
  "Non-nil if running Epoch.")

(if (and mime-running-epoch mime-use-highlighting)
    (progn
      (defvar mime-deemphasize-style (make-style))
      (set-style-foreground mime-deemphasize-style mime-deemphasize-color)
      (defvar mime-emphasize-style (make-style))
      (set-style-foreground mime-emphasize-style mime-emphasize-color)))

(if (and mime-running-lemacs mime-use-highlighting)
    (progn
      (defvar mime-deemphasize-style (make-face 'mime-deemphasize-face))
      (set-face-foreground mime-deemphasize-style mime-deemphasize-color)
      (defvar mime-emphasize-style (make-face 'mime-emphasize-face))
      (set-face-foreground mime-emphasize-style mime-emphasize-color)))

(if (and mime-running-fsf19 mime-use-highlighting)
    (progn
      (defvar mime-deemphasize-style 'mime-deemphasize-face)
      (make-face 'mime-deemphasize-face)
      (set-face-foreground mime-deemphasize-style mime-deemphasize-color)
      (defvar mime-emphasize-style 'mime-emphasize-face)
      (make-face 'mime-emphasize-face)
      (set-face-foreground mime-emphasize-style mime-emphasize-color)))

(defvar mime-audio-file "/tmp/.fooblatz"
  "Filename to store audio snippets recorded on the fly.")

(defvar mime-audio-tmp-file "/tmp/.fooblatz.aiff"
  "Filename to store audio snippets recorded on the fly.")

(defconst mime-waiting-message-lines
  '("Mail mime-compose bug reports to marca@ncsa.uiuc.edu and pray for help."
    "For the daring: ftp.ncsa.uiuc.edu:/outgoing/marca/mime-compose.el"
    "Feature requests?  Fervent wishes?  Unfulfilled desires?  Write code!"
    "mime-compose.el: the Kitchen Sink(tm) of mail composers."
    "Q: How many Elisp hackers does it take to change a light bulb?"
    "A: None -- we glow in the dark."
    ".gnol oot yaw rof scamE gnisu neeb ev'uoy ,siht daer nac uoy fI"
    "Macs?  We don' need no steenkin Macs!  We got MIME!"
    "All hail MIME.  All hail MIME.  Yay.  Yay.  Woo.  Woo.")
  "List of stupid strings to display while waiting for more to do.")

;;; --------------------------- Utility functions ----------------------------

(defun mime-primary-boundary ()
  "Return the current primary boundary.  Note that in the current version
of mime-compose.el, there is no support for secondary boundaries (for
parallel or alternate bodyparts, etc.).  In the future, there may be."
  (let ((date (current-time-string)))
    (concat (if (not (string= (substring date 8 9) " "))
		(substring date 8 9))
	    (substring date 9 10) "." (substring date 4 7) 
	    "." (substring date -4 nil) mime-primary-boundary)))

(defun mime-hide-region (from to hideflag)
  "Hides or shows lines from FROM to TO, according to HIDEFLAG:
If T, region is hidden, else if NIL, region is shown."
  (let ((old (if hideflag ?\n ?\^M))
        (new (if hideflag ?\^M ?\n))
        (modp (buffer-modified-p)))
    (unwind-protect (progn
                      (subst-char-in-region from to old new t))
      (set-buffer-modified-p modp))))

(defun mime-maybe-hide-region (start end)
  "Hide the current region if mime-use-selective-display is T."
  (if mime-use-selective-display
      (mime-hide-region start end t)))

(defun mime-add-description (description)
  "Add a description to the current MIME message element."
  (interactive "sDescription: ")
  (save-excursion
    (if (re-search-backward (concat "--" (mime-primary-boundary))
                            (point-min) t)
        (progn
          (next-line 2)
          (insert "Content-Description: " description "\n")))))

(defun mime-display-waiting-messages ()
  "Display cute messages until input arrives.  Shamelessly stolen
>from VM, the Kitchen Sink(tm) of mail readers."
  (if mime-use-waiting-messages
      (progn
        (if (sit-for 2)
            (let ((lines mime-waiting-message-lines))
              (message
               "mime-compose.el $Revision$, by marca@ncsa.uiuc.edu")
              (while (and (sit-for 4) lines)
                (message (car lines))
                (setq lines (cdr lines)))))
        (message "")
        (if (not (input-pending-p))
            (progn
              (sit-for 2)
              (if (not (input-pending-p))
                  (mime-display-waiting-messages)))))))

;;; ------------------------------ Highlighting ------------------------------

(if mime-use-highlighting
    (progn
      (if mime-running-fsf19
	  (defun mime-add-zone (start end style)
	    (let ((ovl (make-overlay start end)))
	      (overlay-put ovl 'face style))))
      (if mime-running-lemacs
          (defun mime-add-zone (start end style)
            "Add a Lucid Emacs extent from START to END with STYLE."
            (let ((extent (make-extent start end)))
              (set-extent-face extent style)
              (set-extent-data extent 'mime-compose))))
      (if mime-running-epoch
          (defun mime-add-zone (start end style)
            "Add an Epoch zone from START to END with STYLE."
            (let ((zone (add-zone start end style)))
              (epoch::set-zone-data zone 'mime-compose))))))

(defun mime-maybe-highlight-region (start end)
  "Maybe highlight a region of text.  Region is from START to END."
  (if (and (or mime-running-epoch mime-running-lemacs mime-running-fsf19)
           mime-use-highlighting)
      (progn
        (mime-add-zone start end mime-deemphasize-style)
        (save-excursion
          (goto-char start)
          (if (re-search-forward "Content-Type: " end t)
              (let ((s (match-end 0)))
                (re-search-forward "[;\n]")
                (mime-add-zone 
                 s (- (match-end 0) 1) mime-emphasize-style)))))))

;;; -------------------------- mime-mimify-message ---------------------------

(defun mime-mimify-message ()
  "Add MIME headers to a message.  Add an initial informational message
for mail readers that don't process MIME messages automatically.  Add
an initial area for plaintext.  Add a closing boundary at the end of
the message.

This function is safe to call more than once."
  (interactive)
  (let ((mail-header-separator (if (eq major-mode 'mh-letter-mode)
                                   "\n\n\\|^-+$"
                                 mail-header-separator)))
    (or
     (save-excursion
       (goto-char (point-min))
       (re-search-forward "^Mime-Version: "
                          (save-excursion
                            (goto-char (point-min))
                            (re-search-forward mail-header-separator)
                            (point))
                          t))
     (let ((mime-virgin-message (save-excursion
                                  (next-line -1)
                                  (looking-at mail-header-separator))))
       (if mime-virgin-message
           (insert "\n"))
       (save-excursion
         (save-excursion
           (goto-char (point-min))
           (re-search-forward mail-header-separator)
           (beginning-of-line)
           (insert "Mime-Version: 1.0\n")
           (insert "Content-Type: multipart/mixed; boundary=" (mime-primary-boundary) "\n")
           (mime-maybe-highlight-region (save-excursion (next-line -3) (point))
                                        (- (point) 1))
           (next-line 1)
           (let ((start (point)) end)
             (insert
              "> If you are reading this, your mail reader may not support MIME.\n")
             (insert
              "> Some parts of this message will be readable as plain text.\n")
             (setq end (point))
             (mime-maybe-hide-region start (- end 1)))
           (insert "\n")
           (goto-char (point-max))
           (insert "\n")
           (insert "\n")
           (insert "--" (mime-primary-boundary) "--\n")
           (mime-maybe-highlight-region (save-excursion (next-line -1) (point))
                                        (- (point) 1)))
         (save-excursion
           (goto-char (point-min))
           (re-search-forward mail-header-separator)
           (beginning-of-line)
           ;; THIS HAS TO MATCH the number of lines of text included
           ;; as a message ``header'' above.
           (if mime-use-selective-display
               (next-line 3)
             (next-line 4))
           (insert "--" (mime-primary-boundary) "\n")
           (insert "Content-Type: text/plain\n")
           (mime-maybe-highlight-region
            (save-excursion (next-line -2) (point))
            (- (point) 1))
           (insert "\n"))
         (if mime-virgin-message
             (backward-delete-char 1))))))
  (if (interactive-p)
      (mime-display-waiting-messages)))

(defun mime-open-text-bodypart ()
  "At current point, just open up a new plaintext bodypart."
  (interactive)
  (mime-mimify-message)
  (push-mark)
  (let ((start (point)) end)
    (insert "--" (mime-primary-boundary) "\n")
    (insert "Content-Type: text/plain")
    (setq end (point))
    (insert "\n\n")
    (mime-maybe-highlight-region start end))
  (mime-display-waiting-messages))

;;; ---------------------------- file inclusions -----------------------------

(defun mime-include-file (filename content-type binary &optional charset)
  "Include a file named by FILENAME and with MIME content type
CONTENT-TYPE.  If third argument BINARY is T, then the file is binary;
else it's text.  Optional fourth arg CHARSET names character set for
data.  Data will be encoded in base64 or quoted-printable format as
appropriate."
  (mime-mimify-message)
  (push-mark)
  (insert "--" (mime-primary-boundary) "\n")
  (insert "Content-Type: " content-type)
  (if charset
      (insert "; charset=" charset))
  (if (and mime-name-included-files (not (string= filename mime-audio-file)))
      (insert "; name=\"" (file-name-nondirectory filename) "\""))
  (insert "\n")
  (if (not (string= filename mime-audio-file))
      (insert "Content-Description: " filename "\n"))
  (if binary
      (insert "Content-Transfer-Encoding: base64\n")
    (insert "Content-Transfer-Encoding: quoted-printable\n"))
  (mime-maybe-highlight-region 
   (save-excursion (re-search-backward 
                    (concat "--" (mime-primary-boundary))) (point))
   (- (point) 1))
  (insert "\n")
  (let ((start (point)) (seldisp selective-display))
    (setq selective-display nil)
    (if binary
	(shell-command (concat mime-encode-base64-command " < " filename) t)
      (shell-command (concat mime-encode-qp-command  " < " filename) t))
    (goto-char (mark t))
    (setq selective-display seldisp)
    (next-line 1)
    (mime-maybe-hide-region start (1- (point)))
    (insert "\n")
    (insert "--" (mime-primary-boundary) "\n")
    (insert "Content-Type: text/plain\n")
    (mime-maybe-highlight-region 
     (save-excursion (re-search-backward 
                      (concat "--" (mime-primary-boundary))) (point))
     (- (point) 1))
    (insert "\n\n")
    (next-line -1)))

(defun mime-include-binary-file (filename content-type)
  "Include a binary file named by FILENAME at point in a MIME message.
CONTENT-TYPE names MIME content type of file.  Data will be encoded in
base64 format."
  (mime-include-file filename content-type t))

(defun mime-include-nonbinary-file (filename content-type &optional charset)
  "Include a nonbinary file named by FILENAME at point in a MIME
message.  CONTENT-TYPE names MIME content type of file; optional third
arg CHARSET names MIME character set.  Data will be encoded in
quoted-printable format."
  (mime-include-file filename content-type nil charset))

;;; -------------------------- external references ---------------------------

;; ron@auda.mlfarm.com, 12 Aug 1993
;; added prompt for transfer mode

(defun mime-include-external (site directory name content-type description 
                                   mode access-type)
  "Include an external pointer in a MIME message.  Args are SITE,
DIRECTORY, NAME, CONTENT-TYPE, DESCRIPTION, MODE, and ACCESS-TYPE; these are
all strings."
  (mime-mimify-message)
  (push-mark)
  (insert "--" (mime-primary-boundary) "\n")
  (insert "Content-Type: message/external-body;\n")
  (insert "\taccess-type=\"" access-type "\";\n")
  (insert "\tsite=\"" site "\";\n")
  (insert "\tdirectory=\"" directory "\";\n")
  (if (not (string= mode ""))
      (insert "\tmode=\"" mode "\";\n"))
  (insert "\tname=\"" name "\"\n")
  (insert "Content-Description: " description "\n")
  (insert "\n")
  (insert "Content-Type: " content-type "\n")
  (mime-maybe-highlight-region 
   (save-excursion (re-search-backward 
                    (concat "--" (mime-primary-boundary))) (point))
   (- (point) 1))
  (insert "\n")
  (insert "\n")
  (insert "--" (mime-primary-boundary) "\n")
  (insert "Content-Type: text/plain\n")
  (mime-maybe-highlight-region 
   (save-excursion (re-search-backward 
                    (concat "--" (mime-primary-boundary))) (point))
   (- (point) 1))
  (insert "\n"))

(defun mime-include-external-anonftp (site directory name description mode)
  "Include an external pointer (anonymous FTP) in a MIME message.  Args
are SITE, DIRECTORY, MODE, NAME, and DESCRIPTION; these are all strings, and
if interactive, will be prompted for."
  (interactive 
   "sFTP site: \nsRemote directory name: \nsRemote filename: \nsDescription: \nsTransfer mode: ")
  (let ((content-type 
         (completing-read "Content type: " mime-valid-include-types
                          nil nil nil)))
    ;; Unadvertised default.
    (if (string= content-type "")
        (setq content-type "application/octet-stream"))
    (mime-include-external site directory name content-type 
                           description mode "anon-ftp"))
  (mime-display-waiting-messages))

(defun mime-include-external-ftp (site directory name description mode)
  "Include an external pointer (regular FTP) in a MIME message.  Args
are SITE, DIRECTORY, NAME, and DESCRIPTION; these are all strings, and
if interactive, will be prompted for."
  (interactive 
   "sFTP site: \nsRemote directory name: \nsRemote filename: \nsDescription: \nsTransfer mode: ")
  (let ((content-type 
         (completing-read "Content type: " mime-valid-include-types
                          nil nil nil)))
    ;; Unadvertised default.
    (if (string= content-type "")
        (setq content-type "application/octet-stream"))
    (mime-include-external site directory name content-type 
                           description mode "ftp"))
  (mime-display-waiting-messages))

;;; ------------------------------ window dumps ------------------------------

(defun mime-include-xwd-dump ()
  "Run program named by 'mime-xwd-command' and include the results in
a MIME message."
  (interactive)
  (mime-mimify-message)
  (push-mark)
  (insert "--" (mime-primary-boundary) "\n")
  (insert "Content-Type: image/x-xwd\n")
  (insert "Content-Description: Window dump from " (system-name) "\n")
  (insert "Content-Transfer-Encoding: base64\n")
  (mime-maybe-highlight-region 
   (save-excursion (re-search-backward 
                    (concat "--" (mime-primary-boundary))) (point))
   (- (point) 1))
  (insert "\n")
  (let ((start (point)) end (seldisp selective-display))
    (next-line 1)
    (save-excursion
      (next-line -1)
      (message "When crosshair cursor appears, click on window...")
      (sit-for 0)
      (call-process "/bin/sh" nil t nil "-c" mime-xwd-command)
      (message "")
      (sit-for 0))
    (setq end (point))
    (setq selective-display nil)
    (shell-command-on-region start end mime-encode-base64-command t)
    (setq selective-display seldisp)
    (setq end (point))
    (mime-maybe-hide-region start (- end 1))
    (insert "\n")
    (insert "--" (mime-primary-boundary) "\n")
    (insert "Content-Type: text/plain\n")
    (mime-maybe-highlight-region 
     (save-excursion (re-search-backward 
                      (concat "--" (mime-primary-boundary))) (point))
     (- (point) 1))
    (insert "\n\n")
    (next-line -1))
  (mime-display-waiting-messages))

;;; ----------------------------- audio snippets -----------------------------

(defun mime-sgi-grab-audio-snippet ()
  "Grab an audio snippet into file named in 'mime-audio-file'.
This routine works on SGI Indigo's and 4D/35's."
  (let (audio-process done-flag)
    (setq audio-process 
          (start-process "snippet" "snippet" 
                         "/usr/sbin/recordaiff" "-n" "1" "-s" "8" "-r" "8000"
                         mime-audio-tmp-file))
    ;; Quick hack to make Emacs sit until recording is done.
    (setq done-flag
          (y-or-n-p "Press y when done recording (n to abort): "))
    (interrupt-process "snippet")
    ;; Wait until recordaiff has written data to disk.
    (while (eq (process-status "snippet") 'run)
      (message "Waiting...")
      (sleep-for 1))
    (message "Done waiting.")
    ;; Kill off recordaiff and our buffer.
    (delete-process "snippet")
    (kill-buffer "snippet")
    ;; Remove the old mulaw file and do the conversion.
    (call-process "/bin/rm" nil nil nil "-f" mime-audio-file)
    (if done-flag
        (call-process "/usr/sbin/sfconvert" nil nil nil mime-audio-tmp-file
                      mime-audio-file "-o" "mulaw"))
    (call-process "/bin/rm" nil nil nil "-f" mime-audio-tmp-file)
    ;; Return done flag.  If nil, mime-include-audio-snippet should
    ;; clean up.
    done-flag))

(defun mime-sun-grab-audio-snippet ()
  "Grab an audio snippet into file named in 'mime-audio-file'.
This is the Sun version.  I don't know if it works.  I don't have a
SPARCstation to test on at the moment."
  (let (audio-process done-flag)
    (setq audio-process
          (start-process "snippet" "snippet"
;                         "/bin/sh" "-c" "/bin/cat" "<" "/dev/audio"
;                         ">" mime-audio-file))
			 "/usr/demo/SOUND/record" mime-audio-file))

    ;; Quick hack to make Emacs sit until recording is done.
    (setq done-flag
          (y-or-n-p "Press y when done recording (n to abort): "))
    (interrupt-process "snippet")
    ;; Wait until the record process is done.
    (while (eq (process-status "snippet") 'run)
      (message "Waiting...")
      (sleep-for 1))
    (message "Done waiting.")
    ;; Kill off the record process and our buffer.
    (delete-process "snippet")
    (kill-buffer "snippet")
    ;; Return done flag.  If nil, mime-include-audio-snippet should
    ;; clean up.
    done-flag))

(defun mime-include-audio-snippet ()
  "Record a snippet of audio in a MIME message.  This should work on
both Silicon Graphics and Sun platforms.  Code contributions for other
platforms are welcome."
  (interactive)
  (let ((mime-grab-audio-snippet
         (if mime-using-silicon-graphics
             'mime-sgi-grab-audio-snippet
           'mime-sun-grab-audio-snippet)))
    (if (eq (funcall mime-grab-audio-snippet) t)
        (progn
          (mime-include-binary-file mime-audio-file "audio/basic")
          (save-excursion
            (next-line -4)
            (mime-add-description 
             (concat (user-full-name) " " 
                     mime-babbling-description ".")))
	  (delete-file mime-audio-file))))
  (mime-display-waiting-messages))

;;; ------------------------- Basic include commands -------------------------

(defun mime-include-gif (filename)
  "Include a GIF file named by FILENAME."
  (interactive "fGIF image filename: ")
  (mime-include-binary-file filename "image/gif")
  (mime-display-waiting-messages))

(defun mime-include-jpeg (filename)
  "Include a JPEG file named by FILENAME."
  (interactive "fJPEG image filename: ")
  (mime-include-binary-file filename "image/jpeg")
  (mime-display-waiting-messages))

(defun mime-include-acrobat (filename)
  "Include an Acrobat file named by FILENAME."
  (interactive "fAcrobat filename: ")
  (mime-include-binary-file filename "application/pdf")
  (mime-display-waiting-messages))

(defun mime-include-audio (filename)
  "Include an audio file named by FILENAME.  Note that to match the
MIME specification for audio/basic, this should be an 8-bit mulaw file."
  (interactive "fAudio filename: ")
  (mime-include-binary-file filename "audio/basic")
  (mime-display-waiting-messages))

(defun mime-include-postscript (filename)
  "Include a PostScript file named by FILENAME."
  (interactive "fPostScript filename: ")
  (mime-include-nonbinary-file filename "application/postscript")
  (mime-display-waiting-messages))

(defun mime-include-raw-binary (filename)
  "Include a raw binary file named by FILENAME."
  (interactive "fRaw binary filename: ")
  (let ((content-type 
         (completing-read "Content type (RET for default): " 
                          mime-valid-include-types
                          nil nil nil)))
    (if (string= content-type "")
        (setq content-type "application/octet-stream"))
    (mime-include-binary-file filename content-type))
  (mime-display-waiting-messages))

(defun mime-include-raw-nonbinary (filename &optional prefix-arg)
  "Include a raw nonbinary file named by FILENAME.  With prefix arg,
prompt for character set."
  (interactive "fRaw nonbinary filename: \nP")
  (let ((charset
         (if prefix-arg
             (completing-read "Character set: " mime-valid-charsets
                              nil nil nil)
           mime-default-charset))
        (content-type 
         (completing-read "Content type (RET for default): " 
                          mime-valid-include-types
                          nil nil nil)))
    (if (string= content-type "")
        (setq content-type "text/plain"))
    (if (string= charset "")
        (setq charset "asdfasdfdfsdafs"))
    (mime-include-nonbinary-file filename content-type charset))
  (mime-display-waiting-messages))

;;; ---------------------------- Region commands -----------------------------

(defun mime-encode-region (start end content-type charset)
  "Encode a region specified by START and END.  CONTENT-TYPE and
CHARSET name the content type and character set of the data in the
region."
  ;; Start by encoding the region in quoted-printable.  This will
  ;; move end, but not start.
  (goto-char end)
  (let ((seldisp selective-display))
    (setq selective-display nil)
    (shell-command-on-region start end mime-encode-qp-command t)
    (setq selective-display seldisp))
  ;; Now pick up the new end.
  (setq end (point))
  ;; Pop up to start and insert the header; this will also change
  ;; end, but with save-excursion we'll end up at the new end.
  (save-excursion
    (goto-char start)
    (push-mark)
    (insert "--" (mime-primary-boundary) "\n")
    (insert "Content-Type: " content-type "; charset=" charset "\n")
    (insert "Content-Transfer-Encoding: quoted-printable\n")
    (mime-maybe-highlight-region 
     (save-excursion (re-search-backward 
                      (concat "--" (mime-primary-boundary))) (point))
     (- (point) 1))
    (insert "\n"))
  ;; Pick up the new end again.
  (setq end (point))
  ;; Insert the trailing boundary and the new text/plain header.
  (insert "\n")
  (insert "--" (mime-primary-boundary) "\n")
  (insert "Content-Type: text/plain\n")
  (mime-maybe-highlight-region 
   (save-excursion (re-search-backward 
                    (concat "--" (mime-primary-boundary))) (point))
   (- (point) 1))
  (insert "\n")
  ;; Last but not least, add MIME headers if necessary.
  (save-excursion
    (mime-mimify-message)))

(defvar rich-substitutions
  '(
    ("< "    "<lt>")
    ("\\b\\*\\B" "</bold>")
    ("\\b~\\B" "</italic>")
    ("\\B~\\b" "<italic>")
    ("\\B\\*\\b" "<bold>")
    ("\\B_\\b" "<underline>")
    ("\\b_\\B" "</underline>")
    ("\n\n\n" "<np>")
    ("\n\n" "<nl>")
    ("\n"   " ")
    )
  "A table of REGEXP to translate text to MIME's text/richtext format.")

(defun mime-region-to-richtext (start end &optional prefix-arg)
  "Convert the current region to MIME richtext.  MIME headers are
added if necessary; a MIME boundary is added at the start of the
region to indicate richtext; the conversion (see below) is done; a new
boundary is added for more text.

With prefix arg, prompt for character set; else use value of
mime-default-charset.

Text parsing: *foo* becomes <bold>foo</bold>, ~foo~ becomes
<italic>foo</italic>, and _foo_ becomes <underline>foo</underline>.
The font-directives do not work after punctuation, like *foo!*.  Blank
lines are converted to <nl>."
  (interactive "r\nP")
  (let ((charset
         (if (not prefix-arg)
             mime-default-charset
           (completing-read "Character set: " mime-valid-charsets
                            nil nil nil)))
	(subs rich-substitutions)
        pat rep)
    (save-excursion
      (save-restriction
	(narrow-to-region start end)
	(while subs
	  (setq pat (car (car subs)))
	  (setq rep (car (cdr (car subs))))
	  (setq subs (cdr subs))
	  (goto-char start)
	  (while (re-search-forward pat (point-max) t)
	    (replace-match rep)))
	(setq end (point-max))))
    ;; Unadvertised default.
    (if (string= charset "")
        (setq charset mime-default-charset))
    (mime-encode-region start end "text/richtext" 
                        charset))
  (mime-display-waiting-messages))


(defun mime-region-to-charset (start end)
  "Convert the current region to plaintext in a non-default character
set.  You are prompted for a character set, and the text in the region
is encoded in quoted-printable format and identified as being in that
character set."
  (interactive "r")
  (let ((charset
         (completing-read "Character set: " mime-valid-charsets
                          nil nil nil)))
    ;; Unadvertised default.
    (if (string= charset "")
        (setq charset mime-default-charset))
    (mime-encode-region start end "text/plain" charset))
  (mime-display-waiting-messages))

;;; -------------------------------- Keymaps ---------------------------------

;;; Add functions to MH letter mode.
(if mime-running-mh-e
    ;; Running mh-e.
    (if (or (not (boundp 'mh-letter-mode-mime-map)) 
            (not mh-letter-mode-mime-map))
        (progn
          (setq mh-letter-mode-mime-map (make-sparse-keymap))
          (define-key mh-letter-mode-map "\C-c\C-m" mh-letter-mode-mime-map)
          (define-key mh-letter-mode-mime-map "m" 'mime-mimify-message)
          (define-key mh-letter-mode-mime-map "g" 'mime-include-gif)
          (define-key mh-letter-mode-mime-map "j" 'mime-include-jpeg)
          (define-key mh-letter-mode-mime-map "a" 'mime-include-audio)
          (define-key mh-letter-mode-mime-map "p" 'mime-include-postscript)
          (define-key mh-letter-mode-mime-map "r" 'mime-include-raw-binary)
          (define-key mh-letter-mode-mime-map "n" 'mime-include-raw-nonbinary)
          (define-key mh-letter-mode-mime-map "x" 'mime-include-xwd-dump)
          (define-key mh-letter-mode-mime-map "e" 
            'mime-include-external-anonftp)
          (define-key mh-letter-mode-mime-map "f" 
            'mime-include-external-ftp)
          (define-key mh-letter-mode-mime-map "s"
            'mime-include-audio-snippet)
          (define-key mh-letter-mode-mime-map "\C-r" 'mime-region-map)))
  ;; Not running mh-e.
  (progn
    (define-key mail-mode-map "\C-cm" 'mime-mimify-message)
    (define-key mail-mode-map "\C-cg" 'mime-include-gif)
    (define-key mail-mode-map "\C-cj" 'mime-include-jpeg)
    (define-key mail-mode-map "\C-ca" 'mime-include-audio)
    (define-key mail-mode-map "\C-cA" 'mime-include-acrobat)
    (define-key mail-mode-map "\C-cp" 'mime-include-postscript)
    (define-key mail-mode-map "\C-cr" 'mime-include-raw-binary)
    (define-key mail-mode-map "\C-cn" 'mime-include-raw-nonbinary)
    (define-key mail-mode-map "\C-cx" 'mime-include-xwd-dump)
    (define-key mail-mode-map "\C-ce" 'mime-include-external-anonftp)
    (define-key mail-mode-map "\C-cf" 'mime-include-external-ftp)
    (define-key mail-mode-map "\C-cs" 'mime-include-audio-snippet)
    
    ;; Functions that operate on regions.
    (defvar mime-region-map (make-sparse-keymap))
    (define-key mail-mode-map "\C-c\C-r" mime-region-map)
    (define-key mime-region-map "r" 'mime-region-to-richtext)
    (define-key mime-region-map "i" 'mime-region-to-charset)))

  
;;; -------------------------------- Menubar ---------------------------------

;;; Change the menubar for emacs19 - must be done after loading sendmail
(defun mime-emacs19-menubar ()
  (define-key mail-mode-map [menu-bar mime]
	(cons "Mime" (make-sparse-keymap "Mime")))
  (define-key mail-mode-map [menu-bar mime mimeftp]
	'("Include External FTP" . mime-include-external-ftp))
  (define-key mail-mode-map [menu-bar mime mimeanon]
	'("Include External AnonFTP" . mime-include-external-anonftp))
  (define-key mail-mode-map [menu-bar mime mimenorm]
	'("Include Raw Nonbinary File" . mime-include-raw-nonbinary))
  (define-key mail-mode-map [menu-bar mime mimebin]
	'("Include Raw Binary File" . mime-include-raw-binary))
  (define-key mail-mode-map [menu-bar mime mimeaudio]
	'("Include Audio Snippet" . mime-include-audio-snippet))
  (define-key mail-mode-map [menu-bar mime mimeaud]
	'("Include Audio File" . mime-include-audio))
  (define-key mail-mode-map [menu-bar mime mimexwd]
	'("Include XWD Dump" . mime-include-xwd-dump))
  (define-key mail-mode-map [menu-bar mime mimepos]
	'("Include Postscript" . mime-include-postscript))
  (define-key mail-mode-map [menu-bar mime mimeacro]
	'("Include Acrobat" . mime-include-acrobat))
  (define-key mail-mode-map [menu-bar mime mimejpg]
	'("Include JPEG File" . mime-include-jpeg))
  (define-key mail-mode-map [menu-bar mime mimegif]
	'("Include GIF File" . mime-include-gif)))

;; All we do at the moment is replace the popup menu defined in
;; Lucid Emacs 19.3's sendmail.el.
(and 
 mime-running-lemacs 
 (setq 	mail-mode-menu
       '("Mail Mode"
	 "Sending Mail:"
	 "----"
	 ["Send and Exit"		mail-send-and-exit		t]
	 ["Send Mail"			mail-send			t]
	 ["Sent Via"			mail-sent-via			t]
	 "----"
	 "Go to Field:"
	 "----"
	 ["To:"				mail-to				t]
	 ["Subject:"			mail-subject			t]
	 ["CC:"				mail-cc				t]
	 ["BCC:"			mail-bcc			t]
	 ["Text"			mail-text			t]
	 "----"
	 "Miscellaneous Commands:"
	 "----"
	 ["Yank Original"		mail-yank-original		t]
	 ["Fill Yanked Message"		mail-fill-yanked-message	t]
	 ["Insert Signature"		mail-signature			t]
	 "----"
	 "MIME Inclusions:"
	 "----"
	 ["Include GIF File"		mime-include-gif		t]
	 ["Include JPEG File"		mime-include-jpeg		t]
	 ["Include Audio File"		mime-include-audio		t]
	 ["Include PostScript File"	mime-include-postscript		t]
	 ["Include Acrobat File"	mime-include-acrobat		t]
	 ["Include XWD Dump"		mime-include-xwd-dump		t]
	 ["Include Audio Snippet"	mime-include-audio-snippet	t]
	 ["Include Raw Binary File"	mime-include-raw-binary		t]
	 ["Include Raw Nonbinary File"	mime-include-raw-nonbinary	t]
	 ["Include External AnonFTP"	mime-include-external-anonftp	t]
	 ["Include External FTP"	mime-include-external-ftp	t]
	 "----"
	 ["Abort" kill-buffer t]
	 )))

;;; ----------------------------- New mail-send ------------------------------

;; If we're not running Lemacs, pop in a new mail-send routine.
(if (not mime-running-lemacs)
    (defun mail-send ()
      "Send the message in the current buffer.
If  mail-interactive  is non-nil, wait for success indication
or error messages, and inform user.
Otherwise any failure is reported in a message back to
the user from the mailer."
      (interactive)
      (message "Sending...")
      (run-hooks 'mail-send-hook)
      (funcall send-mail-function)
      (set-buffer-modified-p nil)
      (delete-auto-save-file-if-necessary)
      (message "Sending...done")))

;;; --------------------------------- Hooks ----------------------------------

;; Author: Daniel LaLiberte (liberte@cs.uiuc.edu).
(defun mime-postpend-unique-hook (hook-var hook-function)
  "Postpend HOOK-VAR with HOOK-FUNCTION, if it is not already an element.
hook-var's value may be a single function or a list of functions."
  (if (boundp hook-var)
      (let ((value (symbol-value hook-var)))
        (if (and (listp value) (not (eq (car value) 'lambda)))
            (and (not (memq hook-function value))
                 (set hook-var (append value (list hook-function))))
          (and (not (eq hook-function value))
               (set hook-var (append value (list hook-function))))))
    (set hook-var (list hook-function))))

(defun mime-unfrob-selective-display ()
  "Turn off selective display throughout this buffer."
  (if mime-use-selective-display
      (progn
        (message "Unfrobbing selective-display...")
        (mime-hide-region (point-min) (point-max) nil))))

(defun mime-strip-useless-bodyparts ()
  "Strip useless (empty) bodyparts out of a message."
  (save-excursion
    (goto-char (point-min))
    (while (re-search-forward
            (concat "^--" (mime-primary-boundary)
                    "\nContent-Type: text.*[\n]*--" (mime-primary-boundary))
            (point-max) t)
					;added second argument to
					;replace-match to prevent
					;upcasing of initial letter
					;ron@auda.mlfarm.com, 12 Aug 1993
      (replace-match (concat "--" (mime-primary-boundary)) t)
      ;; Go all the way back up to start over.
      (goto-char (point-min)))))

(defun mime-encode-region-qp (start end)
  "Encode a region specified by START and END in quoted-printable
format.  Return the new endpoint.  Do not use save-excursion."
  ;; Start by encoding the region in quoted-printable.  This will
  ;; move end, but not start.
  (goto-char end)
  (let ((seldisp selective-display))
    (setq selective-display nil)
    (shell-command-on-region start end mime-encode-qp-command t)
    (setq selective-display seldisp)))

(defun mime-encode-plaintext ()
  "Encode all plaintext bodyparts in the message in quoted-printable
and set the charset to mime-default-charset."
  (save-excursion
    (goto-char (point-min))
    ;; We're looking for text/plain bodyparts with no extra fields.
    (while (re-search-forward
            (concat "^--" (mime-primary-boundary)
                    "\nContent-Type: text/plain\n") (point-max) t)
      (let* ((head (match-beginning 0))
             (start (match-end 0))
             ;; Assume there's a closing boundary; go find it.
             (end (save-excursion (re-search-forward
                                   (concat "^--" (mime-primary-boundary)))
                                  (- (match-beginning 0) 1))))
        ;; Maybe there's already a Content-Transfer-Encoding.  If so,
        ;; never mind.
        (or (re-search-forward "^Content-Transfer-Encoding: " end t)
            (let ((new-end (save-excursion
                             (mime-encode-region-qp start end))))
              (save-excursion
                (goto-char head)
                (next-line 1)
                (end-of-line)
                (let ((s (point)))
                  (insert "; charset=" mime-default-charset "\n")
                  (insert "Content-Transfer-Encoding: quoted-printable")
                  (mime-maybe-highlight-region s (point))))))))))

(defun mime-send-hook-function ()
  "Function to be called from mail-send-hook.  Unfrob selective
display if active, strip out empty (useless) bodyparts, and optionally
encode plaintext bodyparts in quoted-printable with a given charset."
  (mime-unfrob-selective-display)
  (mime-strip-useless-bodyparts)
  (and mime-encode-plaintext-on-send
       (mime-encode-plaintext)))

;; Before the message is sent, remove the selective display crap.
(if mime-running-mh-e
    (mime-postpend-unique-hook 'mh-before-send-letter-hook
                               'mime-send-hook-function)
  (mime-postpend-unique-hook 'mail-send-hook 'mime-send-hook-function))

; ron@auda.mlfarm.com, 11 Aug 1993 protests, also below
;(defun mime-setup-hook-function ()
;  (if mime-use-selective-display
;      (setq selective-display t)))

;; During mail setup, activate selective-display if necessary.  We use
;; mail-mode-hook rather than mail-setup-hook because if a message is
;; being composed and C-x m gets hit again, mail-mode will be
;; reentered, causing selective-display to revert to nil and possibly
;; screwing up the display bigtime unless mail-mode-hook knows what to
;; do.
; (if mime-running-mh-e
;    (mime-postpend-unique-hook 'mh-letter-mode-hook
;                               'mime-setup-hook-function)
;  (mime-postpend-unique-hook 'mail-mode-hook 'mime-setup-hook-function))


(if mime-running-fsf19
	(mime-emacs19-menubar)) ;; take care of the menubar in emacs19

;; Kai.Grossjohann@cs.uni-dortmund.de: In order to set up mime-compose
;; with Gnus in Emacs 20, you need the following lines.
;
; (add-hook 'message-send-hook 'mime-send-hook-function)
; (define-key message-mode-map "\C-cm" 'mime-mimify-message)
; (define-key message-mode-map "\C-cg" 'mime-include-gif)
; (define-key message-mode-map "\C-cj" 'mime-include-jpeg)
; (define-key message-mode-map "\C-ca" 'mime-include-audio)
; (define-key message-mode-map "\C-cA" 'mime-include-acrobat)
; (define-key message-mode-map "\C-cp" 'mime-include-postscript)
; (define-key message-mode-map "\C-cr" 'mime-include-raw-binary)
; (define-key message-mode-map "\C-cn" 'mime-include-raw-nonbinary)
; (define-key message-mode-map "\C-cx" 'mime-include-xwd-dump)
; (define-key message-mode-map "\C-ce" 'mime-include-external-anonftp)
; (define-key message-mode-map "\C-cf" 'mime-include-external-ftp)
; (define-key message-mode-map "\C-cs" 'mime-include-audio-snippet)
;
; ;; Functions that operate on regions.
; (defvar mime-region-map (make-sparse-keymap))
; (define-key message-mode-map "\C-c\C-r" mime-region-map)
; (define-key mime-region-map "r" 'mime-region-to-richtext)
; (define-key mime-region-map "i" 'mime-region-to-charset)
;
; (define-key message-mode-map [menu-bar mime]
;       (cons "Mime" (make-sparse-keymap "Mime")))
; (define-key message-mode-map [menu-bar mime mimeftp]
;       '("Include External FTP" . mime-include-external-ftp))
; (define-key message-mode-map [menu-bar mime mimeanon]
;       '("Include External AnonFTP" . mime-include-external-anonftp))
; (define-key message-mode-map [menu-bar mime mimenorm]
;       '("Include Raw Nonbinary File" . mime-include-raw-nonbinary))
; (define-key message-mode-map [menu-bar mime mimebin]
;       '("Include Raw Binary File" . mime-include-raw-binary))
; (define-key message-mode-map [menu-bar mime mimeaudio]
;       '("Include Audio Snippet" . mime-include-audio-snippet))
; (define-key message-mode-map [menu-bar mime mimeaud]
;       '("Include Audio File" . mime-include-audio))
; (define-key message-mode-map [menu-bar mime mimexwd]
;       '("Include XWD Dump" . mime-include-xwd-dump))
; (define-key message-mode-map [menu-bar mime mimepos]
;       '("Include Postscript" . mime-include-postscript))
; (define-key message-mode-map [menu-bar mime mimeacro]
;       '("Include Acrobat" . mime-include-acrobat))
; (define-key message-mode-map [menu-bar mime mimejpg]
;       '("Include JPEG File" . mime-include-jpeg))
; (define-key message-mode-map [menu-bar mime mimegif]
;       '("Include GIF File" . mime-include-gif)))
