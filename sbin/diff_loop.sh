#! /bin/sh

if [ -t 0 ]; then
   echo "Usage: $0 [gdiff options...] < starlab_diff_files
Produce side-by-side diffs for given list of file pairs
as produced by starlab_compare.  Default gdiff options:
  --width=160 --left-column
" >&2
   exit 1
fi

if [ $# = 0 ]; then
  set -- --width=160 --left-column
fi

while read left right; do
  echo "

==========
"
  echo "****  `ls -l $left`"
  echo "****  `ls -l $right`"
  echo "****  diff -y $* $left $right"
  diff -y $* $left $right
done
