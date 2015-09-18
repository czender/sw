# .cshrc

# User specific aliases and functions

alias rm 'rm -i'
alias cp 'cp -i'
alias mv 'mv -i'

alias m  'less'
alias h  'history'
alias csrc 'source ~/.tcshrc'
alias cd 'cd \!*; set prompt=${cwd}" ROOT"#" "'
alias dir 'ls -lga'
setenv PATH "/usr/sbin:/sbin:/bin:/usr/bin:/usr/local/bin:$PATH"
