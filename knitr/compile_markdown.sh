#!/bin/sh

if [ ! -f $1 ]; then
    echo "File $1 not found"
    exit 1
fi

name=$(basename $1 .md)
path=$(dirname $1)

# create pdf version
pandoc --smart -f markdown-markdown_in_html_blocks   -o $path/$name.pdf $1
# create HTML version
pandoc --toc  --mathjax  -c github-pandoc.css -f markdown -t html -o $path/$name.html $1

# convert links to HTML
sed -i "" 's/\.md/.html/g' $path/$name.html
