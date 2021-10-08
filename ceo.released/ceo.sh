#!/bin/sh
dir=`dirname "$0"`

msa=$dir/example/msa1.txt
msaout=$dir/example/msa1.flt.txt
msacluout=$dir/example/msa1.flt.clu.txt
msahtmlout=$dir/example/msa1.clu.html

ceo cmd=clu fs=0.7 fi=0.98 msa=$msa msaout=$msaout msacluout=$msacluout msahtmlout=$msahtmlout cfrom=0.65 cto=0.95 csteps=6
