#!/bin/bash

# https://gist.github.com/aaronwolen/d7db57b83211eacd4f61
# Aaron Wolen 2016
# 
# Partially convert Rnw presentations to Rmd syntax
# - [x]: Code chunks
# - [x]: Section headers
# - [x]: Slide headers
# - [x]: Presenter notes
# - [x]: Lists
# - [ ]: Inlinde code (sort of, not really)

input=cmst.Rnw
output=cmst.Rmd

cat "$input" \
	| perl -pe 's/^\<{2}(.*)\>{2}=/```{r $1}/' \
	| perl -pe 's/@/```/' \
	| perl -pe 's/\\section{(.*)}/# $1/' \
	| perl -pe 's/(\\begin{frame})(\[fragile\])?{(.*)}/## $3/g' \
	| perl -pe 's/\\end{frame}//' \
	| perl -pe 's/\\note{/<div class="notes">/' \
	| perl -pe 's/^}/<\/div>/' \
	| perl -pe 's/\\(begin|end){itemize}//' \
	| perl -pe 's/\s*\\item\s*/* /' \
| perl -pe 's/\\src{(.*)}/`$1`/g' > "$output"